function spm8w_voitool(varargin) 
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: May, 2010 - DDW
% ==============================================================================
% function spm8w_voitool('mode','sub_list','voi_para_file')
%
% spm8w_voitool performes a number of functions related to timeseries analysis 
% (VOI based timeseries exctraction, PPI regressor creation and PPI plotting).
% Function is defined in 'type'. 
%
% Type: 'VOI'
% extract the first eigenvariate (or mean) from a volume of interest and saves 
% the resulting timeseries as a text file in a user specified onsets directory.
% VOI can be extracted from group map or on a subject by subject basis. 
% However, note that once a threshold is used to define a VOI, voitool will
% proceed to use the entirety of the voxels within the VOI radius, irrespective
% of whether they met the voi defining threshold or not. This is to avoid 
% uneven number of voxels within VOIs across subjects as normally SPM will 
% only calculate eigenvariates on the above threshold voxels within a VOI 
% which can lead to VOIs with only one voxel contributing to the timeseries.
% Note that voxels at edge of bigmask will have ROIs with less voxels than 
% they should (since the voi extends beyond bigmask).
%
% Type: 'PPI'
% Creates PPI regressors (PPI, Y and P) based on a predefined VOI and 
% additional specs (see VOI_studyname.m). PPI, Y and P are saved as text files
% in a user specified onsets directory.
%
% Type: 'PPI_PLOT'
% Creates a PPI plot per subject and saves as gif or png file. In addition it
% saves the necessary data as an excel compatible text file in the user VOI
% directory (in the future I'll move everything into a studyroot/VOI directory).
% This plot is the standard way fo plotting PPI data and consists of the PPI
% regressor for each condition and region. 
%
% spm8w_voitool takes a list of subjects (one per line), analysis type (i.e. 
% VOI or PPI or PPI_PLOT) and VOI_name.m file as input, this file contains all
% the parameters defining VOI extraction, PPI regressor creation or PPI 
% plotting parameters.
% ==============================================================================
% CHANGE LOG:
% -VOI matches spm8 output from eigenvariate button!!! -DDW May /10
% -Added PPI support -DDW May/10
% -Improved code and added PPI plotting capability. -DDW June /10
% -Rewrote sections for more sensible batching based on VOI_xxx.m -DDW June/10
% -Major re-write. Fixed loc_max (now properly searched for nearest local
%  maxima rather than nearest suprathreshold voxel). -DDW Nov/12
% -Added support for image based ROI (.nii) files. -DDW Nov/12
% -Added routine to optimize local maxima (taking top of all loc max within
%  a pre-defined distance). -DDW Nov/12
% =======1=========2=========3=========4=========5=========6=========7=========8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the current working dir
cwd = pwd;
% Load user defaults, otherwise load vanilla defaults
if exist(which('spm8w_userdefaults.m'), 'file')
    [def_dir,def_file] = fileparts(which('spm8w_userdefaults.m'));
else
    [def_dir,def_file] = fileparts(which('spm8w_defaults.m'));
end
% Goto dir, eval defaults, check for errors, come back to cwd
cd(def_dir); 
try
    eval(def_file); 
catch
    defmsg = sprintf('The defaults file %s.m in %s is not evaluable, please check your syntax\n',def_file,def_dir); 
    error_reporter(cwd, defmsg, 1);
    error('Exiting...');
end
cd(cwd); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch (nargin)
    case 0 
        error('Please specify subjects and analysis type, e.g. spm8w_voitool(''PPI'',''09jan09aa'')');
    case 1
        voi_type    = varargin{1};
        subjects    = spm8w_getsub;
        v           = spm8w_getp('VOI'); 
    case 2
        voi_type    = varargin{1};
        subjects    = varargin{2};
        v           = spm8w_getp('VOI'); 
    case 3  
        voi_type    = varargin{1};
        subjects    = varargin{2};
        v           = spm8w_getp('VOI',[],varargin{3});
    otherwise
        error('You''ve specified too many inputs. Either 1 or 2 or 3 only please.');
end

%%%Check cell type and fix if necessary (we want cells to be consistent)
if ~iscell(subjects)
    subjects = cellstr(subjects);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Super switch! VOI or PPI or PPI_PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch(voi_type)
    case('VOI')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Adjust VOI specs for all_subjects and 
        %%% generate a new v.voi_specs with 1 subj 
        %%% per row. Rather than loop through subjects
        %%% we loop through the new v.voi_specs. This
        %%% allows us to mix and match 'all_subjects'
        %%% and subj specific coordinates.May2010 DDW
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        indx = 1; %Sets the index for adjusting size of v.voi_specs
        for i = 1:size(v.voi_specs,1)
            if (strmatch(v.voi_specs{i,1},'all_subjects'))
                tmp_voi_specs(indx:indx+length(subjects)-1,1) = subjects;
                tmp_voi_specs(indx:indx+length(subjects)-1,2) = v.voi_specs(i,2);
                tmp_voi_specs(indx:indx+length(subjects)-1,3) = v.voi_specs(i,3);
                tmp_voi_specs(indx:indx+length(subjects)-1,4) = v.voi_specs(i,4);
                indx = indx + length(subjects);
            else
                tmp_voi_specs(indx,1) = v.voi_specs(i,1);
                tmp_voi_specs(indx,2) = v.voi_specs(i,2);
                tmp_voi_specs(indx,3) = v.voi_specs(i,3);
                tmp_voi_specs(indx,4) = v.voi_specs(i,4);
                indx = indx + 1;
            end    
        end
        %voi_specs wants her kit back.
        v.voi_specs=tmp_voi_specs;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Check directories
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %%% Check all Results directories
        for i = 1:size(v.voi_specs,1)
                if ~isdir(deblank(eval(v.glm_dir)))
                    error(['The GLM directory at %s does not exist '...
                           'please check your paths.'], eval(v.glm_dir));
                end
        end
        %%% Check Onsets directory
        if ~isdir(v.ons_dir)
            fprintf('The directory %s does not exist, creating...\n', v.ons_dir);
            [s,w] = system(spm8w_osbabel(sprintf('mkdir "%s"',v.ons_dir)));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Loops through v.voi_specs creating VOIs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        spm('Defaults','fmri')
        start_time = datestr(now);     
        fprintf('\n==========Creating VOIs at %s\n', start_time); 
        %%% Init the diagnostics variables
        v.voi_diag_d(size(v.voi_specs,1),1:3)= {'','',''}; 
        %%% Drive by each subject's GLM dir and deliver a voi.
        for i = 1:size(v.voi_specs,1)
            fprintf('===Working on subject %s and region %s...\n',v.voi_specs{i,1},v.voi_specs{i,2});
            %%% Set the diag
            v.voi_diag_d{i,1} = v.voi_specs{i,2};   %VOIname
            v.voi_diag_d{i,2} = v.voi_specs{i,1};   %Currentsubject
            %%% If defining VOI based on GLM then need below. When we do
            %%% just pure VOI extraction (timeseries from preprocessed data 
            %%% maybe with minimal adjustment (motion reg) then we can put
            %%% an IF statemetn here. -Nov 2012
            %%% Goto glm dir
            cd(eval(v.glm_dir));
            if ~isdir('VOI')
                [s,w] = system(spm8w_osbabel(sprintf('mkdir "%s"','VOI')));                
            end
            % Find the contrast index for the named contrast
            load('SPM.mat');
            % Check for contrast if only 1 then error
            if length(SPM.xCon) <= 1
                error('GLM for subject %s appears to be estimated but no contrasts have been calculated...', v.voi_specs{i,1});
            end
            for ii = 1:length(SPM.xCon)
                if strcmpi(v.con_name, SPM.xCon(ii).name)
                    con_num = ii;
                    break;
                end
            end
            %Define xSPM structure prior to loading contrast data
            xSPM = struct('swd', eval(v.glm_dir),...    %location of GLM
                          'title', [],...               %title will be filled in
                          'Ic', con_num,...             %con number found above
                          'Im', [],...                  %Implict masking 
                          'u', v.con_pval,...           %p value
                          'k', v.con_ext,...            %extent threshold
                          'thresDesc', 'none');         %correction type
            [SPM, xSPM] = spm_getSPM(xSPM);             %Performs contrast and puts above threshold voxels into xSPM
            %%% Quick check if coords are part of contrast (for
            %%% diagnostics). This find vector in a matrix (but transposes
            %%% col to row vector to use ismember appropriately).
            % find(ismember(xSPM.XYZmm',str2num(v.voi_specs{i,3}),'rows')==1)
            %%% Check for suprathreshold voxels
            if isempty(xSPM.XYZ)
                fprintf('No suprathreshold voxels for subjects:%s... Skipping...', v.voi_specs{i,1}); 
                v.voi_diag_d(i,3)={'No Suprathreshold Voxels'};
                return 
            end
            %%% DEFINE ROI (image, sphere, search locmax)
            xY_xyz = str2num(v.voi_specs{i,3})';  %init what will become the final xyz if adjustemnets are made for locmax
            %%% Check if v.voi_specs 4 is a number for sphere using cheap hack.     
            if isempty(str2num(v.voi_specs{i,4}))  
                fprintf('Performing VOI extraction using image mask...\n');
                maskfile = spm8w_osbabel(fullfile(v.roi_img,v.voi_specs{i,4}));
                xY = struct('name', v.voi_specs{i,2},...          %Name for VOI file
                            'Ic', 1,...                           %1 is effects of interest
                            'Sess',1,...                          %We always concat sessions
                            'def','mask',...                      %Sphere or mask
                            'spec', maskfile,...                  %sphere size or filename
                            'str', sprintf('ROI mask: %s',v.voi_specs{i,4}));      
                fprintf('Recalculating contrast to include all voxels in image MASK...\n');
                clear xSPM
                xSPM = struct('swd', eval(v.glm_dir),...    %location of GLM
                          'title', [],...               %title will be filled in
                          'Ic', con_num,...             %con number found above
                          'Im', [],...                  %Implict masking 
                          'u', 1,...                    %p value
                          'k', 0,...                    %extent threshold
                          'thresDesc', 'none');         %correction type
                [SPM, xSPM] = spm_getSPM(xSPM);     
                v.voi_diag_d(i,3) = {'No Search'};
            else 
                %%% Search for local max of distance is over v.con_maxdist, drop warning
                if v.con_locmax
                    fprintf('Searching for nearest local maxima to %s...\n',v.voi_specs{i,3});
                    % Below is code to find nearest loc max. Might be worth
                    % wraping this in its own matlab function for future use. 
                    XYZfull = SPM.xVol.XYZ; %All voxel indeces from stat map volume
    %                 ResMS   = SPM.VResMS;   %Volume info of ResMS image
    %                 %This next command feeds the con or resms vol info to spm_get_data
    %                 %with the XYZ coordiantes of the full volume 
    %                 cBetas  = spm_get_data(SPM.xCon(con_num).Vcon, XYZfull);
    %                 resids  = spm_get_data(ResMS, XYZfull);
    %                 %The following is some scaling variable?  
    %                 VcB = SPM.xCon(con_num).c'*SPM.xX.Bcov*SPM.xCon(con_num).c;
    %                 Z   = cBetas./sqrt(resids*VcB);
                    %%% So the Z above is IDENTICAL to using spm_get_data on
                    %%% the Tstat map so let's just do that since it's
                    %%% computed.
                    %T = spm_get_data(SPM.xCon(con_num).Vspm, XYZfull);
                    T = spm_get_data(SPM.xCon(con_num).Vspm, xSPM.XYZ);                
                    %%% The t-stat map is unthresholded, so we should
                    %%% threshold. We know the significnat voxel coordinates from xSPM.XYZ
                    %T3 = T(ismember(XYZfull',xSPM.XYZ','rows')');
                    %%% Find nearest locmax
                    target     = str2num(v.voi_specs{i,3}); %Coordinates to start search from
                    [~,~,XYZsig,~] = spm_max(T, xSPM.XYZ);  %Find local maximas (replace T with Z)
                    %Transform XYZsig to XYZsigmm for use with spm_XYZreg (which
                    %requires vox in mm). 
                    M    = SPM.xCon(con_num).Vspm.mat;
                    iM   = inv(M); %inverted transform in case we ever need it later
                    XYZsig_mm = M(1:3,:)*[XYZsig; ones(1,size(XYZsig,2))]; %From voxel to mm space
                    %Just take nearest
                    [xyz,~,locmax_dist] = spm_XYZreg('NearestXYZ',target,XYZsig_mm);               
                    %Nearest locmac and its distance from original xyz
                    if locmax_dist > v.con_maxdist
                        fprintf(['Distance between local maxima %d,%d,%d & coord %s '...
                                'is larger than %dmm...\n'],xyz(1),xyz(2),xyz(3),v.voi_specs{i,3},v.con_maxdist);
                    else
                        fprintf('Distance from %s & local maxima %d,%d,%d is %4.1fmm...\n',...
                                v.voi_specs{i,3}, xyz(1),xyz(2),xyz(3), locmax_dist);
                    end
                    xY_xyz = xyz;   %Assign new xyz to xY structure
                    v.voi_diag_d(i,3) = {num2str(locmax_dist)};
                    %%%IMPORTANT FIX - July 2010
                    %%%Because we're using a contrast to define the center of
                    %%%the sphere it can happen that not all voxels surrounding
                    %%%the peak found by spm_XYZreg above are actually supra-
                    %%%threshold and so they will not be in the xSPM structure
                    %%%Therefore we recalc xSPM using a p=1 so that the VOI 
                    %%%actually contains all voxels defined by the radius 
                    %%%v.voi_specs{i,4}. Otherwise you can get VOIs with only 
                    %%%1 voxel in them. -DDW
                    fprintf('Recalculating contrast to include all voxels in VOI...\n');
                    clear xSPM
                    xSPM = struct('swd', eval(v.glm_dir),...    %location of GLM
                              'title', [],...               %title will be filled in
                              'Ic', con_num,...             %con number found above
                              'Im', [],...                  %Implict masking 
                              'u', 1,...                    %p value
                              'k', 0,...                    %extent threshold
                              'thresDesc', 'none');         %correction type
                    [SPM, xSPM] = spm_getSPM(xSPM);  
                 else
                    v.voi_diag_d(i,3) = {'No Search'};
                end                       
                %%% Previously I did it all by hand, see prior versions of spm
                %%% for that code (or bottom of this file). in case we ever
                %%% want to revert. xY is the voi structure, xY.XYZm is the
                %%% voxel (mm) that encompass the ROI. Q is the voxel indices 
                % Make xY structure for sphere using xyz figured out above
                xY = struct('xyz', xY_xyz,...                     %Coordinates determine above (either original or shifted)
                            'name', v.voi_specs{i,2},...          %Name for VOI file
                            'Ic', 1,...                           %1 is effects of interest
                            'Sess',1,...                          %We always concat sessions
                            'def','sphere',...                    %Sphere or mask
                            'spec', str2num(v.voi_specs{i,4}),... %sphere size or filename
                            'str', sprintf('%0.1fmm sphere',str2num(v.voi_specs{i,4})));    
            end
            %%% Define ROI based on xY from above.                                      
            [xY,xY.XYZmm,Q] = spm_ROI(xY, xSPM.XYZmm);
            %%% EXTRAC THE DATUMS
            %%% Could use spm_regions.m but that forces us to use
            %%% eignvariates instead of mean timeseries
            y = spm_get_data(SPM.xY.VY,xSPM.XYZ(:,Q)); %Get data using indices from spm_ROI.
            %Since we don't filter, this does nothing. In the future we may
            %want to enable this, but in our case SPM.xX.K is inf. 
            y = spm_filter(SPM.xX.K,SPM.xX.W*y);  %No harm in leaving it on.
            %%%Adjust for null space of effects of interest contrast.
            beta = spm_get_data(SPM.Vbeta,xSPM.XYZ(:,Q));
            y = y - spm_FcUtil('Y0',SPM.xCon(xY.Ic),SPM.xX.xKXs,beta);
            %%%the following lines might be only for session stuff
            xY.X0 = SPM.xX.xKXs.X(:,[SPM.xX.iB SPM.xX.iG]);
            j = SPM.Sess(xY.Sess).row;
            y = y(j,:);
            xY.X0 = xY.X0(j,:);
            xY.X0 = [xY.X0 SPM.xX.K(xY.Sess).X0];
            xY.X0 = xY.X0(:,any(xY.X0));
            %%% So far we have the same xY.X0 that we started with???
            %%% Compute first eigenvariate (vtmp otehrwise overwrite v struct)
            [m n] = size(y);
            if m > n
                [vtmp s vtmp] = svd(y'*y);
                s       = diag(s);
                vtmp    = vtmp(:,1);
                u       = y*vtmp/sqrt(s(1));
            else
                [u s u] = svd(y*y');
                s       = diag(s);
                vtmp    = u(:,1);
                u       = y'*u/sqrt(s(1));
            end
            d    = sign(sum(vtmp));
            u    = u*d;
            vtmp = vtmp*d;
            Y    = u*sqrt(s(1)/n);
            %%% Add final data to structure
            %%% Corrcoeff the xY.y_mean with xY.u to 
            %%% see that the mean timeseries and first eigenvariate are similar.
            xY.y      = y;
            xY.m      = mean(y,2); %mean timeseries along 2nd dimension
            xY.u      = Y;         %eigenvariate
            xY.v      = vtmp;      %eigenimage
            xY.s      = s;         %eigenvalues
            % Print correltion between mean and eigen
            r = corrcoef(xY.m,xY.u);
            fprintf('Correlation (eigenvariate & mean timeseries): %4.3f\n', r(1,2));
            % Save VOI.mat
            cd('VOI');
            save([xY.name, '.mat'],'xY');  
            fprintf('VOI saved to %s... Size of VOI is: %d voxels...\n', [xY.name,'.mat'], size(xY.y,2));
            % Save text file
            cd(v.ons_dir);
            if ~v.eigenvar
                tmp = xY.m;
                save([v.voi_specs{i,1}, '_',xY.name,'.txt'],'tmp','-ASCII');
            else
                tmp = xY.u;
                save([v.voi_specs{i,1}, '_',xY.name,'.txt'],'tmp','-ASCII');
            end
            fprintf('Timeseries saved in: %s\n\n', v.ons_dir);
            % ET PHONE HOME
            cd(v.root);
            % CLEANUP
            clear Q SPM Y beta con_num d indx ii j m n r s tmp_voi_specs tmp u vtmp tmp xSPM xY y
        end
        fprintf('Completed VOI extraction...\n');
        if v.con_pval < 1
            fprintf('VOI Diagnostics\n===============\n')
            if v.con_locmax == 1; tmp = 'Yes';else tmp = 'No';end;
            fprintf('Con p-value:%4.4f\tCon extent:%d\nCon Search:%s\t\tCon Distance:%d\n',...
                    v.con_pval,v.con_ext,tmp,v.con_maxdist);
            %%%Sortation magic parsing the results of the VOI extraction
            %%%sorting into keep and reject bins. Not very elegant but works
            %%%june 2010 -ddw
            v.voi_diag_d = sortrows(v.voi_diag_d,[1,2]);
            cur_region   = v.voi_diag_d{1,1};
            cur_region_i = 1;
            cur_keep_i   = 1;
            cur_reject_i = 1;
            diags.keep = {};
            diags.reject = {};
            for i = 1:size(v.voi_diag_d,1)
                if ~isempty(strmatch(v.voi_diag_d{i,1},cur_region))  %while working on same region
                    if str2num(v.voi_diag_d{i,3}) <= v.con_maxdist
                        diags(cur_region_i).keep(cur_keep_i,:) = v.voi_diag_d(i,:);
                        cur_keep_i = cur_keep_i + 1;
                    else
                        diags(cur_region_i).reject(cur_reject_i,:) = v.voi_diag_d(i,:);
                        cur_reject_i = cur_reject_i + 1;          
                    end                    
                else
                    cur_region = v.voi_diag_d{i,1};       %set region to next region
                    cur_region_i = cur_region_i + 1;      %add to region counter
                    cur_keep_i = 1; cur_reject_i = 1;     %reset counters
                    if str2num(v.voi_diag_d{i,3}) <= v.con_maxdist
                        diags(cur_region_i).keep(cur_keep_i,:) = v.voi_diag_d(i,:);
                        cur_keep_i = cur_keep_i + 1;
                    else
                        diags(cur_region_i).reject(cur_reject_i,:) = v.voi_diag_d(i,:);
                        cur_reject_i = cur_reject_i + 1;          
                    end 
                end
            end
            %%%equate the field of diag
            for i = 1:length(diags)
                if size(diags(i).keep,1) > size(diags(i).reject,1)
                    diags(i).reject(size(diags(i).keep,1),:) = {'','',''};
                elseif size(diags(i).keep,1) < size(diags(i).reject,1)
                    diags(i).keep(size(diags(i).reject,1),:) = {'','',''};                
                end    
            end            
            %%%
            for i = 1:length(diags)
                %Print the VOI
                if isempty(diags(i).keep{1,1}); 
                    fprintf('\nVOI: %s\n',diags(i).reject{1,1}); 
                else
                    fprintf('\nVOI: %s\n',diags(i).keep{1,1});
                end
                if v.con_locmax == 1
                    fprintf('Subjects within %dmm\tSubject beyond %dmm\t\n',v.con_maxdist,v.con_maxdist);
                else
                    fprintf('Subjects with suprathreshold voxels\tSubjects without suprathreshold voxels\n');
                end
                %print data
                for ii = 1:size(diags(i).keep,1)
                    fprintf('%s(%smm)\t\t\t%s(%sfmm)\n',diags(i).keep{ii,2},diags(i).keep{ii,3},diags(i).reject{ii,2},diags(i).reject{ii,3});
                end
            end 
            %%%Cleanup
            clear cur_keep_i cur_region cur_region_i cur_reject_i diags distance i ii tmp xyz
        end

       
    case('PPI')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Check directories
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Check all Results directories
        for i = 1:length(subjects)
                if ~isdir(deblank(eval(v.ppi_dir)))
                    error(['The GLM directory at %s does not exist '...
                           'please check your paths.'], eval(v.ppi_dir));
                end
        end
        %%% Check Onsets directory
        if ~isdir(v.ons_dir)
            fprintf('The directory %s does not exist, creating...\n', v.ons_dir);
            [s,w]=unix(['mkdir ', v.ons_dir]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Loops through subjects creating PPIs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        start_time = datestr(now);   
        fprintf('\n==========Creating PPI regressors at %s\n\n', start_time); 
        %%% Drive by each subject's GLM dir and create PPI regressors.
        for i = 1:length(subjects)
            fprintf('======Working on subject %s...\n',subjects{i});
            for ii = 1:size(v.ppi_specs,1)
                fprintf('===Creating PPI for region: %s\n',v.ppi_specs{ii,1});
                cd(eval(v.ppi_dir));
                if ~isdir('VOI')
                    error(['The VOI directory at %s does not exist '...
                           'please check your paths.'], eval(v.ppi_dir));
                end
                %Create the PPI con array
                uU = zeros(length(v.ppi_specs{ii,3}),3);
                for iii = 1:length(v.ppi_specs{ii,3})
                    uU(iii,:) = [iii,1,v.ppi_specs{ii,3}(iii)];
                end
                %Prepare for PPI blastoff
                load('SPM.mat');
                %Run the PPI
                PPI = spm_peb_ppi(SPM,'ppi',[eval(v.ppi_dir),'/VOI/',v.ppi_specs{ii,1},'.mat'],uU,v.ppi_specs{ii,2},v.ppi_graphics);
                %Close the interactive window
                G = spm_figure('FindWin','Interactive'); close(G);
                if(v.ppi_graphics)
                    %print figure to pdf if ppi_graphics is at 1
                    F   = spm_figure('FindWin',1); 
                    %Get figure 1 (the PPI figure doesn't have the graphics tag as it should, bad spm programmers, baaaaad)
                    %%%Check for previous model files.
                    if exist(['../../',subjects{i},'_ppi.pdf'],'file')
                    eval('print -depsc2 -r300 -painters -noui ../../tmp_ppi.ps');
                    [s,w]     = unix('ps2pdf ../../tmp_ppi.ps ../../tmp_ppi.pdf');
                    mergepdfs = ['gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=../../tmp_ppi_merge.pdf -dBATCH ../../',subjects{i},'_ppi.pdf ../../tmp_ppi.pdf'];   
                    [s,w]=unix(mergepdfs);
                    mergename = ['mv ../../tmp_ppi_merge.pdf ../../',subjects{i},'_ppi.pdf'];
                    [s,w]     = unix(mergename);
                    [s,w]     = unix(['rm ../../tmp_ppi*.*']);
                    else
                    eval('print -depsc2 -r300 -painters -noui ../../tmp_ppi.ps');
                    [s,w]     = unix(['ps2pdf ../../tmp_ppi.ps ../../',subjects{i},'_ppi.pdf']);
                    [s,w]     = unix(['rm ../../tmp_ppi.ps']);
                    end
                    %%% Close figures
                    close(F);
                end
                %Move PPI.mat to PPI directory
                if ~isdir('PPI')
                   [s,w]=unix('mkdir PPI'); 
                end
                [s,w]=unix(['mv PPI_',v.ppi_specs{ii,2},'.mat ./PPI']); 
                %Print correltion between PPI, Y and P
                r = corrcoef(PPI.ppi,PPI.Y);
                fprintf('Correlation PPI & Y: %4.3f\n',abs(r(1,2)));
                r = corrcoef(PPI.ppi,PPI.P);
                fprintf('Correlation PPI & P: %4.3f\n',abs(r(1,2)));
                r = corrcoef(PPI.Y,PPI.P);
                fprintf('Correlation Y and P: %4.3f\n',abs(r(1,2)));
                %Export PPI regressors to txt files in ONSETS directory          
                cd(v.ons_dir);
                tmp = PPI.ppi; save([subjects{i}, '_',v.ppi_specs{ii,2},'_PPI.txt'],'tmp','-ASCII');
                tmp = PPI.P; save([subjects{i}, '_',v.ppi_specs{ii,2},'_P.txt'],'tmp','-ASCII');
                tmp = PPI.Y; save([subjects{i}, '_',v.ppi_specs{ii,2},'_Y.txt'],'tmp','-ASCII');
                fprintf('PPI, P & Y files saved in: %s\n\n', v.ons_dir);
                % ET PHONE HOME
                cd(v.root); 
                % Cleanup
                clear SPM mergename mergepdfs iii s w xY uU SPM PPI F G tmp r
            end
        end

    case('PPI_PLOT')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Check directories
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Check all Results directories
        for i = 1:length(subjects)
                if ~isdir(deblank(eval(v.ppi_dir)))
                    error(['The GLM directory at %s does not exist '...
                           'please check your paths.'], eval(v.ppi_dir));
                end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Loops through subjects creating reduced
        %%% PPIs for plotting.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        start_time = datestr(now);   
        fprintf('\n==========Creating PPI plots at %s\n\n', start_time); 
        %%% Drive by each subject's GLM dir and create PPI plot
        for i = 1:length(subjects)
            fprintf('======Working on subject %s...\n',subjects{i});
            for ii = 1:size(v.plot_specs,1)
                fprintf('===Creating PPI plot for regions %s & %s\n',v.plot_specs{ii,1}, v.plot_specs{ii,2});
                cd(eval(v.ppi_dir));
                if ~isdir('VOI')
                    error(['The VOI directory at %s does not exist '...
                           'please check your paths.'], eval(v.ppi_dir));
                end               
                %Create the PPI con array
                uU1 = zeros(length(v.plot_specs{ii,3}),3);
                uU2 = zeros(length(v.plot_specs{ii,4}),3);
                for iii = 1:length(v.plot_specs{ii,3})      %For contrast 1
                    uU1(iii,:) = [iii,1,v.plot_specs{ii,3}(iii)];
                end
                for iii = 1:length(v.plot_specs{ii,4})      %For contrast 2
                    uU2(iii,:) = [iii,1,v.plot_specs{ii,4}(iii)];
                end    
                %Prepare for a PPI Plotting party
                load('SPM.mat');
                %Run the 4 reduced PPIs
                % warning on
                PPIR1C1 = spm_peb_ppi(SPM,'ppi',[eval(v.ppi_dir),'/VOI/',v.plot_specs{ii,1},'.mat'],uU1,['plot1',v.plot_specs{ii,1}],0);
                PPIR1C2 = spm_peb_ppi(SPM,'ppi',[eval(v.ppi_dir),'/VOI/',v.plot_specs{ii,1},'.mat'],uU2,['plot2',v.plot_specs{ii,1}],0);
                PPIR2C1 = spm_peb_ppi(SPM,'ppi',[eval(v.ppi_dir),'/VOI/',v.plot_specs{ii,2},'.mat'],uU1,['plot1',v.plot_specs{ii,2}],0);
                PPIR2C2 = spm_peb_ppi(SPM,'ppi',[eval(v.ppi_dir),'/VOI/',v.plot_specs{ii,2},'.mat'],uU2,['plot2',v.plot_specs{ii,2}],0);
                %Close the interactive window
                G = spm_figure('FindWin','Interactive'); close(G);
                %Plot it!
                figure; 
                plot(PPIR1C1.ppi, PPIR2C1.ppi,'k.');
                hold on
                plot(PPIR1C2.ppi, PPIR2C2.ppi,'r.');
                %Add fit lines
                plot(PPIR1C1.ppi, polyval(polyfit(PPIR1C1.ppi, PPIR2C1.ppi,1),PPIR1C1.ppi),'k-','linewidth',1.5)
                plot(PPIR1C2.ppi, polyval(polyfit(PPIR1C2.ppi, PPIR2C2.ppi,1),PPIR1C2.ppi),'r-','linewidth',1.5)
                %Tweak figure
                a = legend([v.plot_specs{ii,5},': ',v.plot_specs{ii,1}, ' vs. ',v.plot_specs{ii,2}],...
                    [v.plot_specs{ii,6},': ',v.plot_specs{ii,1}, ' vs. ',v.plot_specs{ii,2}]);
                set(a,'Fontsize',10,'location','NorthEast','fontname','Times');
                xlabel([v.plot_specs{ii,1},' activity'],'fontsize', 14,'fontname','Times');
                ylabel([v.plot_specs{ii,2},' activity'],'fontsize', 14,'fontname','Times');
                title(['PPI plot for subject',subjects{i}],'fontsize',16,'fontname','Times');      
                %Save figure as image
                if isempty(strmatch(v.plot_export,'none'))
                    name = [subjects{i},'_',v.plot_specs{ii,5},'vs.',v.plot_specs{ii,6},'_',v.plot_specs{ii,1},'vs.',v.plot_specs{ii,2},'.txt'];
                    saveas(gcf,name,v.plot_export); 
                    [s,w]=unix(['mv ',subjects{i},'_',v.plot_specs{ii,5},'*.*','./PPI']); 
                    
                end
                %Export PPI regressors to txt files for offline plotting        
                tmp = PPIR1C1.ppi; save([subjects{i}, '_PLOT_',v.plot_specs{ii,1},'_',v.plot_specs{ii,5},'.txt'],'tmp','-ASCII');
                tmp = PPIR1C2.ppi; save([subjects{i}, '_PLOT_',v.plot_specs{ii,1},'_',v.plot_specs{ii,6},'.txt'],'tmp','-ASCII');
                tmp = PPIR2C1.ppi; save([subjects{i}, '_PLOT_',v.plot_specs{ii,2},'_',v.plot_specs{ii,5},'.txt'],'tmp','-ASCII');
                tmp = PPIR2C2.ppi; save([subjects{i}, '_PLOT_',v.plot_specs{ii,2},'_',v.plot_specs{ii,6},'.txt'],'tmp','-ASCII');
                [s,w]=unix(['mv ',subjects{i},'_PLOT_*.txt ./PPI']); 
                [s,w]=unix(['rm PPI_plot*.mat']); 
                fprintf('PPI plot text files saved in %s PPI directory\n\n', subjects{i});
                % ET PHONE HOME
                cd(v.root); 
                % Cleanup
                clear G PPIR1C1 PPIR1C2 PPIR2C1 PPIR2C2 SPM a s w uU1 uU2 tmp
            end
        end
end