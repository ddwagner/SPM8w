function vxdata = spm8w_roigen(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: April, 2013 - DDW
% ==============================================================================
% spm8w_roigen([xyz],radius or mask,['boldfile.nii'],write,fullmonty,noprog)
%
% Swiss army knife for ROI generation and value extraction. Supports single 
% voxels (radius = 1), spheres (radius in voxels) and image masks (must be in 
% same space as the source image). Volumes from which to extract data can be 
% either 3D (e.g., beta images, anatomical MRIs) or 4D (e.g., pre-processed bold
% data). Finally, roigen can be told to write the ROI file to disk in the
% same space as the original data (i.e., writefile = 1) either for display 
% purposes or for bug checking roigen.
%
% Finally, the fullmonty option (i.e., fullmonty = 1) will ouput the
% following structure:
%   data.Source -the source file from which values were extracted
%   data.Def    -the ROI definition (Coordinate, sphere, radius or img mask)
%   data.XYZmm  -the MNI coordinates of every voxel in the ROI
%   data.volHDR -the volume header information from the source file.    
%   data.Values -the voxel values or voxel timeseries of every voxel per input
%               file in the ROI (this could be used for mvpa analyses). 
%   data.Avg    -the average voxel value or voxel timeseries within the ROI. 
%
% NB: If submitting more than one file, file must be in same space (i.e.,
% don't mix and match MRI and fMRI files in one go!). 
%
% Examples: 
% Single Voxel:
%   vxdata = spm8w_roigen([12,12,12],1,{'swuabold1.nii','swuabold2.nii'})
% 6mm sphere writing sphere as img to disk:
%   vxdata = spm8w_roigen([-30,45,10],6,'beta_0001.img', 1)
% Apriori ROI:
%   vxdata = spm8w_roigen([],'L_AMYG.nii','con_00001.img')
% Fullmonty:
%   data = spm8w_roigen([],'L_AMYG.nii',{'beta_0001.img','beta_0002.img'},0,1)
% 
% ==============================================================================
% CHANGE LOG:
% -1st version, spm8w_roigen integrates and extends aspectes of spm8w_roitool 
%  and spm8w_rawvx so that all ROI based value extractions can eb done with
%  spm8w_roigen. DDW April 2013. 
% OLD CHANGE LOG for make_sphere_mask: 
% -Voxels are now rounded to nearest location gridspace and user is warned.
%  This helps catch errors but also allows for users to input coordinates
%  that are not on the native gridspace (e.g. peaks from other studies).
%  DDW June 2008
% -Removes some flipping fixes that were uncessary and causing crashing on 
%  dartmouth data (our philips data isn't flipped... I think... I hope).
%  DDW June 2008
% -Modified and added as part of spm8w_getdata.m (since we never use it
%  outside of this context. Makes for cleaner files.). Also added some
%  info to the volume description field -DDW Jan 2010
% =======1=========2=========3=========4=========5=========6=========7=========8

%---Input checks
switch (nargin)
  case nargin < 2
    error('Please specify a coordinate: [x,y,z] and a radius (1 or higher)...');
  case 2
    roi       = varargin{1};
    radius    = varargin{2};
    vfiles    = [];
    writefile = 0;
    fullmonty = 0;
    noprog    = 0;
  case 3
    roi       = varargin{1};
    radius    = varargin{2};
    vfiles    = varargin{3};
    writefile = 0;
    fullmonty = 0;
    noprog    = 0;
  case 4
    roi       = varargin{1};
    radius    = varargin{2};
    vfiles    = varargin{3};
    writefile = varargin{4};
    fullmonty = 0;
    noprog    = 0;
  case 5
    roi       = varargin{1};
    radius    = varargin{2};
    vfiles    = varargin{3};
    writefile = varargin{4};
    fullmonty = varargin{5}; 
    noprog    = 0;
  case 6
    roi       = varargin{1};
    radius    = varargin{2};
    vfiles    = varargin{3};
    writefile = varargin{4};
    fullmonty = varargin{5}; 
    noprog    = varargin{6}; 
  otherwise 
    error('Too many parameters.'); 
end
%---Set the mask to either user specified anatomical mask (XYZ as char)
%---Or to bigmask. If want to change bigmask to some other default, here is
%---the spot to edit.
if ischar(radius)
    mask = radius; %replace bigmask with the anatomy mask
else
    mask = which('bigmask.nii');
end

%---Load header info
if isempty(vfiles)
    vfiles = spm_select(Inf,'image','Please select an image or volume file');
end
%--Check if cell(user input) or char array(spm_select output). 
if iscell(vfiles); 
    vfiles = char(vfiles); 
else
    %remove ',1' ot force load full volume
    for vfiles_i = 1:size(vfiles,1)
        vfiles(vfiles_i,:) = strrep(vfiles(vfiles_i,:),',1','  ');
    end     
end

%---Process volumes
for vol_i = 1:size(vfiles,1)
    %tic %timeit
    %--Load one volumne header for img space details (more is slow!)
    vol_hdr = spm_vol([deblank(vfiles(vol_i,:)),',1']);
    %sneak the dimensions from nifti extension since we're only loading one
    %vol of the 4d file for speed the 4th dimensions won't be in the regular 
    %spm_vol structure. 
    N = vol_hdr.private.dat.dim;  %-number of images
    if length(N) < 4
        N = 1;    %If there's no 4th dimension then there's only 1 volume
    else
        N = N(4); %Number of volumes = size of 4th dimension
    end
    fname = spm_str_manip(vol_hdr.fname,'t'); 
    if ~noprog  %check if no progress information is disabled.
        fprintf('Extracting data from file: %s...', fname);
    end
    %--Generate ROIs (only on the first date)
    if vol_i == 1
        %--Get space details (borrowed from spm_searchlight.m)
        %--------------------------------------------------------------------------               
        M            = vol_hdr.mat;                        %-voxels to mm matrix
        iM           = inv(M);                             %-mm to voxels matrix
        DIM          = vol_hdr.dim;                        %-image dimensions
        [x,y,z]      = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
        XYZ          = [x(:)';y(:)';z(:)']; clear x y z    %-voxel coordinates {vx}
        XYZmm        = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))]; clear XYZ %-voxel coordinates {mm}
        %--Check if XYZ is coordinates or mask
        if ischar(roi)
            %Load mask, find where mask = 1 and apply to Y
            mask_hd  = spm_vol(mask);
            mask_vol = spm_read_vols(mask_hd);
            mask_fname = spm_str_manip(mask_hd.fname,'t');
            mask_vx  = length(find(mask_vol));
        else
            %-Check if XYZ is in voxel step size.
            loc=round(roi/abs(M(1,1)))*abs(M(1,1)); 
            if (roi(1)==loc(1))==0 || (roi(2)==roi(2))==0 || (roi(3)==loc(3))==0
                fprintf('Warning coordinates have been rounded to nearest voxel in image space\n');
                fprintf('coordinates %s have been rounded to %s \n',num2str(roi),num2str(loc));
                roi = loc;
            end  
            %-Check that XYZ is on the voxel grid in case value is out of range
            %(i.e. -300 is in voxel step size but off grid). 
            if isempty(find(XYZmm(1,:)==loc(1),1)) || isempty(find(XYZmm(2,:)==loc(2),1)) || isempty(find(XYZmm(3,:)==loc(3),1)),
                error('the specified location is not on the grid!');
            end
            %-Get voxel indices for sphere (formula from spm_ROI.m)  
            Q           = ones(1,size(XYZmm,2));
            sphere_idxs = find(sum((XYZmm - roi'*Q).^2) <= radius^2);
            clear Q %free up ram
            %sphere      = XYZmm(:,sphere_idxs); %for debugging
            mask_vol    = zeros(DIM);
            mask_vol(sphere_idxs) = 1;
            mask_fname  = sprintf('sphereRadius%dmm_%d_%d_%d.nii',radius,roi);
            mask_vx     = length(find(mask_vol));
            if writefile == 1
                 %get info from source image
                 vwrite = vol_hdr;
                 vwrite.fname   = mask_fname;
                 vwrite.descrip =sprintf('%dmm ROI mask at coordinates %d %d %d',radius,roi);
                 if ~noprog  %check if no progress information is disabled.
                     fprintf('Writing out ROI file: %s(%d voxels)...',  vwrite.fname,...
                         mask_vx);
                 end
                 spm_write_vol(vwrite,mask_vol);                  
            end        
        end
    end %if vol_i = 1
    %--Extract value(s)
    %-old code (about 5x slower than spm_sample_vol)
%     Y = spm_read_vols(spm_vol(vfiles(vol_i,:)));  %So far everything seems to go in memory fine 
%     if size(Y,4) > 1
%         if ~noprog  %check if no progress information is disabled.
%             fprintf('using %s (%d voxels)...', ...
%                 spm_str_manip(mask_fname,'r'), mask_vx);
%         end
%         %replicate the mask in 4D to apply to Y
%         %mask_vol4D = repmat(mask_vol,[1,1,1,length(vvols{vol_i})]);
%         values = Y(repmat(mask_vol,[1,1,1,N])==1);
%         %clear mask_vol4D 
%         %rehape to voxelXtimeseries array
%         values = reshape(values, mask_vx,N)';       
%     else
%         if ~noprog  %check if no progress information is disabled.
%             fprintf('using %s (%d voxels)...', ...
%                 spm_str_manip(mask_fname,'r'), mask_vx);
%         end
%         values = Y(mask_vol==1)';            
%     end     
%     %cleanup
%     clear Y
    %-New code using spm_sample_vol
    %-Create X Y Z matrices.
    if ~noprog  %check if no progress information is disabled.
        fprintf('using %s (%d voxels)...', ...
                 spm_str_manip(mask_fname,'r'), mask_vx);
    end
    %Convert mm to vx coordinates
    mask_XYZ = XYZmm(:,mask_vol(:,:,:)==1);
    vox_XYZ  = iM * [mask_XYZ; ones(1,size(mask_XYZ,2))];  
    %load the vol info
    vvols = spm_vol(vfiles(vol_i,:));
    %create empty data mat
    values = zeros(N,mask_vx);
    for n_vol = 1:N
       values(n_vol, :) = spm_sample_vol(vvols(n_vol),...
           vox_XYZ(1,:),vox_XYZ(2,:),vox_XYZ(3,:),0);
    end
    %Generate fullmonty output (TODO make sure XYZmm and Values match)
    data(vol_i).Source = deblank(vfiles(vol_i,:));
    data(vol_i).Def    = spm_str_manip(mask_fname,'r'); 
    data(vol_i).XYZmm  = XYZmm(:,mask_vol(:,:,:)==1);
    data(vol_i).volHDR = vol_hdr; %save the volhdr for later writing in same space
    data(vol_i).Values = values;
    data(vol_i).Avg    = nanmean(values,2); %average timeseries in a region
    if ~noprog  %check if no progress information is disabled.
        fprintf('done...');
    end
    if ~noprog  %check if no progress information is disabled.
        %toc
    end
end

%Return output (either fullmonty or just the result)
vxdata = [];
if fullmonty
    vxdata = data;
else
    for vol_i = 1:size(vfiles,1)
        vxdata = [vxdata; data(vol_i).Avg];
    end
end

