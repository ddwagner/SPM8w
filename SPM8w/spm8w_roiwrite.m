function spm8w_roiwrite(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: May, 2013 - DDW
% ==============================================================================
% spm8w_roiwrite(ROI, fname)
%
% Spm8w_roiwrite takes a matlab structure defining the coordinates, image
% space and values for an ROI. This structure was ideally formed by 
% spm8w_roigen and contains, at a minimum, the fields: XYZmm, Values and
% volHDR. XYZmm is a 3x#voxels array, Values is a nX#voxels array and
% volHDR contains the output of spm_vol on the original files (thus
% containing the data defining the space in which the ROI data will be
% written into). Optionally, the user can provide a filename string specifying
% the path and name of the new file.
%
% roiwrite can be used for writing an ROI file (values = ones) to visualize
% a given ROI but can also be used to write a volumn containing voxel
% values calculated by any given function (i.e. in Matlab, or even the
% output of R). The only stipulation is that the data are stored in an
% nX#voxels array where n = 1 for a 3D map (i.e., 1X33 voxels) or n>1 for a
% 4D map. The columns in Values must correspond to the coordinates of the
% columns in XYZmm. Finally, the information in volHDR will be used to
% create a valid hdr into which to write the data.
%
% NOTE: At present roiwrite does not write out 4D files but a series of 3D
% files (if the values array is bigger than 1X#voxels). In a rush now... 
% When/if we need 4D I'll add it... 
%
% Examples: 
%   spm8w_roiwrite(data, './name of some directory/h8tjazz.nii');
% ==============================================================================
% CHANGE LOG:
% =======1=========2=========3=========4=========5=========6=========7=========8

%---Input checks
switch (nargin)
  case nargin < 1
    error('Please specify provide an ROI structure (see spm8w_roigen.m)...');
  case 1
    roi    = varargin{1}.XYZmm;
    values = varargin{1}.Values;
    volhdr = varargin{1}.volHDR;
    fname  = sprintf('./ROI_%dVoxels.nii',size(roi,2));
  case 2
    roi    = varargin{1}.XYZmm;
    values = varargin{1}.Values;
    volhdr = varargin{1}.volHDR;
    [fpath, fname, fext] = fileparts(varargin{2});
    if isempty(fext)
        fext = '.nii';
    end
    fname = [fpath, fname, fext];
  otherwise 
    error('Too many parameters.'); 
end

%---Get space details (borrowed from spm_searchlight.m)
%--------------------------------------------------------------------------               
M            = volhdr.mat;                         %-voxels to mm matrix
iM           = inv(M);                             %-mm to voxels matrix
DIM          = volhdr.dim;                        %-image dimensions
[x,y,z]      = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
XYZ          = [x(:)';y(:)';z(:)']; clear x y z    %-voxel coordinates {vx}
XYZmm        = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))]; clear XYZ %-voxel coordinates {mm}
 
%---Assign value to appropriate array 
values_vol = zeros(DIM);
%--Go through voxels and assign values to approriate places. 
for vx_i = 1:size(roi,2)
    idx_x = find(ismember(XYZmm(1,:),roi(1,vx_i)));
    idx_y = find(ismember(XYZmm(2,:),roi(2,vx_i)));
    idx_z = find(ismember(XYZmm(3,:),roi(3,vx_i)));
    val_idx = intersect(intersect(idx_y,idx_x),idx_z); 
    values_vol(val_idx) = values(1,vx_i);  
end

%---Prep for writing.
volhdr         = rmfield(volhdr, 'private');
volhdr.fname   = fname;
volhdr.dt      = [16 0]; %change to flaoting point if not there already. float32 
                         %should be enough precision
volhdr.descrip = sprintf('%d voxel(s) written out by spm8w_roiwrite.m', size(roi,2));
%volhdr.pinfo   = ?; %not sure if we need to edit pinfo yet.

spm_write_vol(volhdr,values_vol);  
