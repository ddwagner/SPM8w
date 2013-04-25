function spm8w_dartel(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: September, 2012 - DDW
% ==============================================================================
% spm8w_dartel('mode','sub_list','parameter_file')
% 
% spm8w_dartel is used to perform study-wide dartel normalization. The first 
% step involves generating a study-wide template based on segmented anatomical
% MRIs (run spm8w_seg8.m first). After this, invidual anatomical and functional
% files can be warped to template space (study-wise intersubject alignment) 
% and then the template space can be registered to MNI space. 
%
% spm8w_dartel requies a list of subjects, a parameters file and the mode 
% specification (either 'template' or 'mprage2mni' or 'epi2mni'). Before 
% normalizing to MNI, a template must be created.Templates will be stored in 
% the SUBJECTS/DARTEL directory. The tempate needs to be regenerated everytime 
% a subject is added or removed.
%
% TODO: Unfinished attempt to incorporate VBM8 based normalization into the
% pipeline. Paramters file is mostly setup. But VBM8 segmentation and
% flowfield calculation was done using the SPM gui + vbm8 toolbox. Results
% are stored in anatomy_vbm8 (to segregate from dartel based norm). Now we
% need to normalize EPI using this flowfield (which requires hacking the
% vbm8 toolbox since it's not setup for 4D files) and then saving these
% files to the FUNCTIONAL dir with the VBM8 instead of DARTEL token. 
%
% Some tricks for calling dartel w/o spm gui were borrowed from aamod_dartel_*
% by Rhodri Cusack and Jonathan peelle.
%
% Todo: Figure out how to get flowfield without generating a template for
% those special cases where we don't want to include a particular subject
% in template creation (due to bad or missing anatomical).
% ==============================================================================
% CHANGE LOG:
% -First version -DDW Sept/12
% -major edits to accomodate VBM8 style seg, template and norm - DDW Sept/12
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
        error('Please specify the type of dartel process, e.g. spm8w_dartel(''template'')');
    case 1
        dart_type   = varargin{1};
        subjects    = spm8w_getsub;
        para_file   = spm_select(1,'^P_.*\.m$',bbc.pfile,[],[cwd,'/SCRIPTS']);
    case 2
        dart_type   = varargin{1};
        subjects    = varargin{2};
        para_file   = spm_select(1,'^P_.*\.m$',bbc.pfile,[],[cwd,'/SCRIPTS']);
    case 3  
        dart_type   = varargin{1};
        subjects    = varargin{2};
        para_file   = varargin{3};
    otherwise
        error('You''ve specified too many inputs. Either 1 or 2 or 3 only please.');
end

%%%Check cell type and fix if necessary (we want cells to be consistent)
if ~iscell(subjects)
    subjects = cellstr(subjects);
end

switch(dart_type)
    case('template')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Set Defaults and Check files
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% GET A PARA_FILE FOR INITIAL CHECKS
        p = spm8w_getp('P',subjects{1},para_file); %get a para_file for first checks.
        %%% SETUP PARAMS BASED ON DARTEL OR VBM8             
        if p.coreg2epi; mprage = [p.mprage,'_coreg']; else mprage = p.mprage; end
        if strcmp(p.whichseg,'DARTEL'); dpfx = 'rc'; dimp = [dpfx,'1',mprage,'.nii']; end 
        if strcmp(p.whichseg,'VBM8'); dpfx = 'rp'; dimp = [dpfx,'1',mprage,'_affine.nii']; end       
        if strcmp(p.whichtemplate,'DARTEL'); prntstr = 'DARTEL study template generation'; end
        if strcmp(p.whichtemplate,'VBM8'); prntstr = 'DARTEL Registration using the VBM8 template'; end 
        fprintf('=============Beginning %s for %d subjects=============\n', prntstr,length(subjects));
        spm('Defaults','fmri')
        %%% Convert subjects variable to rc1 paths and make sure files exist. 
        fprintf('Checking subjects for segmentation files...');
        for i = 1:length(subjects)
            p = spm8w_getp('P',subjects{i},para_file);
            if strcmp(p.whichseg,'DARTEL'); anatdir = p.anat; end 
            if strcmp(p.whichseg,'VBM8'); anatdir = p.vbm8; end        
            subj_rc{i} = fullfile(anatdir,dimp);           
            if ~exist(subj_rc{i},'file')
                error('Subject %s has not been segmented or is missing dartel ready segmentation files...',subjects{i});
            end
        end
        fprintf('Done...\n');
        %%%Delete any prior dartel files if doing DARTEL style template generation
        if exist([p.dartdir,'/Template_0.nii'],'file') && strcmp(p.whichtemplate,'DARTEL')
            fprintf('Prior dartel templates found... deleting files and proceeding...\n');  
            try
                [s,w] = system(spm8w_osbabel(['rm ',sprintf('"%s/"%s',p.dartdir,'Template_*.nii')]));
            catch
            end
        end
        %%%Add current date to p structure
        start = datestr(now);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Get file locations and setup JOB structure
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        addpath(fullfile(spm('dir'),'toolbox','DARTEL')); %Put DARTEL in path
        %Create a var containing all images for each tissue class
        all_subj_rc{1} = subj_rc;   %Always use grey matter
        for i = 2:p.tissues
            for ii = 1:length(subj_rc)
                subj_rc{ii} = strrep(subj_rc{ii},[dpfx,num2str(i-1)],[dpfx,num2str(i)]);
            end
            all_subj_rc{i} = subj_rc;
        end
        if strcmp(p.whichtemplate,'DARTEL')
            %Setup job structure using dartel defaults (from SPM8 tbx_cfg_dartel)
            param = struct('its',{3,3,3,3,3,3},...
                           'rparam',{[4 2 1e-6],[2 1 1e-6],[1 0.5 1e-6],[0.5 0.25 1e-6],...
                                     [0.25 0.125 1e-6],[0.25 0.12 1e-6]},...
                           'K',{0,0,1,2,4,6},...
                           'slam',{16,8,4,2,1,0.5});
            settings = struct('template', 'Template',...
                              'rform',0,...
                              'param',param,...
                              'optim',struct('lmreg',0.01,'cyc',3,'its',3));
            job = struct('images',{all_subj_rc},'settings',settings);
            fprintf('%s, this may take quite awhile...\n',prntstr);
            spm8w_dartel_template(job);
            fprintf('done\n');
            try
                [s,w] = system(spm8w_osbabel(sprintf('mkdir "%s"',p.dartdir)));
            end
            p = spm8w_getp('P',subjects{1},para_file);
            [s,w] = system(spm8w_osbabel(sprintf('mv "%s"/%s "%s"',p.anat,'Template_*.nii',p.dartdir)));     
        elseif strcmp(p.whichtemplate,'VBM8')
            %Setup job structure for DARTEL existing template
            %get path to VBM8 template files 
            for i = 1:6
                t{i} = which(['Template_',num2str(i),'_IXI550_MNI152.nii']);
            end
            param = struct('its',{3,3,3,3,3,3},...
                           'rparam',{[4 2 1e-6],[2 1 1e-6],[1 0.5 1e-6],[0.5 0.25 1e-6],...
                                     [0.25 0.125 1e-6],[0.25 0.12 1e-6]},...
                           'K',{0,0,1,2,4,6},...
                           'template',{t(1),t(2),t(3),t(4),t(5),t(6)});
            settings = struct('rform',0,...
                              'param',param,...
                              'optim',struct('lmreg',0.01,'cyc',3,'its',3));
            job = struct('images',{all_subj_rc},'settings',settings);
            fprintf('%s, this may take quite awhile...\n',prntstr);      
                spm8w_dartel_warp(job);    
            fprintf('done\n');
        end 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Calculate Time and Save Log
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        stop = datestr(now); 
        p.time_elapsed = etime(datevec(stop),datevec(start)); %time elapsed in seconds
        [hours, minutes, seconds] = spm8w_timecalc(p.time_elapsed);
        fprintf('==============================================================\n');
        fprintf('%s for %d subjects completed at %s\nand took %d hours, %d minutes and %d seconds...\n',prntstr, length(subjects),stop,hours,minutes,seconds);
        if strcmp(p.whichtemplate,'DARTEL'); fprintf('DARTEL Template files are saved to %s\n',p.dartdir); end
    
    case{'mprage2mni','epi2mni'}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Set Defaults and Check files
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        addpath(fullfile(spm('dir'),'toolbox','DARTEL')); %Put DARTEL in path
        %Check if this is single subject or if list was provided
        if length(subjects) == 1; tmpstr = subjects{1};
        else tmpstr = [num2str(length(subjects)), ' subjects'];
        end    
        %%% GET A PARA_FILE FOR INITIAL CHECKS
        p = spm8w_getp('P',subjects{1},para_file); %get a para_file for first checks.
        %%% SETUP PARAMS BASED ON DARTEL OR VBM8  
        if p.coreg2epi; mprage = [p.mprage,'_coreg']; else mprage = p.mprage; end
        if strcmp(p.whichseg,'DARTEL'); dpfx = 'c'; dimp = [dpfx,'1',mprage,'.nii']; end 
        if strcmp(p.whichseg,'VBM8'); dpfx = 'p'; dimp = [dpfx,'1',mprage,'.nii']; end       
        if strcmp(p.whichtemplate,'DARTEL'); normtok = 'd'; end
        if strcmp(p.whichtemplate,'VBM8'); normtok = 'v'; end 
        if strcmp(p.whichseg,'DARTEL') && strcmp(p.whichtemplate,'DARTEL')
            flowfile = ['u_rc1',mprage,'_Template.nii'];
        elseif strcmp(p.whichseg,'DARTEL') && strcmp(p.whichtemplate,'VBM8')
            flowfile = ['u_rc1',mprage,'_affine.nii'];
        elseif strcmp(p.whichseg,'VBM8') && strcmp(p.whichtemplate,'DARTEL')
            flowfile = ['u_rp1',mprage,'_affine_Template.nii'];
        elseif strcmp(p.whichseg,'VBM8') && strcmp(p.whichtemplate,'VBM8')
            flowfile = ['u_rp1',mprage,'_affine.nii'];
        end 
        %%% mprage2mni filechecks
        if strcmp(dart_type,'mprage2mni')
            fprintf('=============Beginning DARTEL normalization of ANATOMY to MNI for %s=============\n', tmpstr);       
            spm('Defaults','fmri')
            %%% Determine if all appropriate files exist.
            fprintf('Checking files...\n');
            for i = 1:length(subjects)
                p = spm8w_getp('P',subjects{i},para_file);   
                if strcmp(p.whichtemplate,'DARTEL'); anatdir = p.anat; end 
                if strcmp(p.whichtemplate,'VBM8'); anatdir = p.vbm8; end        
                % populate file location vars
                subj_c1{i}   = fullfile(anatdir,dimp);
                subj_c2{i}   = fullfile(anatdir,strrep(dimp,['1',mprage],['2',mprage]));
                subj_anat{i} = fullfile(anatdir,[mprage,'.nii']);
                dartflow{i}  = fullfile(anatdir,flowfile);               
                % check if files exist.
                if ~exist(subj_c1{i},'file') && p.coreg2epi
                    error('Subject %s coregistered anatomical files are missing, did you run segmentation?...',subjects{i});
                elseif ~exist(subj_c1{i},'file') && ~p.coreg2epi
                    error('Subject %s segmented tissues (i.e., c1) are missing, did you run segmentation?...',subjects{i});
                end
                if ~exist(dartflow{i},'file')
                    error('Subject %s flow field is missing, did you generate a template with this subject?',subjects{i});
                end             
                %%%Delete any leftover or prior dartel files                 
                if exist(fullfile(anatdir,[normtok,mprage,'.nii']),'file')
                    fprintf('Prior dartel normalization found... deleting files and proceeding...');  
                    try
                        [s,w] = system(spm8w_osbabel(['rm ',fullfile(anatdir,[normtok,mprage,'.nii'])]));
                        [s,w] = system(spm8w_osbabel(['rm ',fullfile(anatdir,['sm',normtok,dpfx,'1',mprage,'.nii'])]));
                        [s,w] = system(spm8w_osbabel(['rm ',fullfile(anatdir,['sm',normtok,dpfx,'2',mprage,'.nii'])]));                          
                    catch
                    end
                end
            end
        fprintf('Done...\n') 
        %%%Start time
        start = datestr(now);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Get file locations and setup JOB structure
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Setup job structure. If Dartel we need to align template to MNI
        %If VBM8 the template is already aligned (i.e. IXI_MN152).
        if strcmp(p.whichtemplate,'DARTEL'); job.template{1} = fullfile(p.dartdir,'Template_6.nii'); end
        if strcmp(p.whichtemplate,'VBM8'); job.template{1} = which('Template_6_IXI550_MNI152.nii'); end     
        %Setup images      
        fprintf('Normalizing %s.nii with DARTEL:\n',mprage);
        job.normtok  = normtok;
        job.bb       = nan(2,3);
        job.vox      = ones(1,3);
        job.fwhm     = zeros(1,3);
        job.preserve = 0; %Preserve concentration for anat and epi
        for i = 1:length(subjects)
            job.data.subj(i).flowfield(1) = dartflow(i);
            job.data.subj(i).images       = subj_anat(i);
        end
        spm8w_dartel_norm_fun(job);          
        fprintf('Normalizing segmented tissues with DARTEL:');
        job.vox      = ones(1,3);
        job.fwhm     = ones(1,3) * p.dartelsmooth;
        job.preserve = 1;   %Modulate for tissue classes
        for i = 1:length(subjects)
            job.data.subj(i).flowfield(1) = dartflow(i);
            job.data.subj(i).images       = [subj_c1(i);subj_c2(i)];
        end
        spm8w_dartel_norm_fun(job); 
        %%%Rename warped files to DARTEL (don't need this anymore, added it
        %%%spm8w_dartel_norm_fun.
%             for i = 1:length(subj_c1)
%                 [s,w] = system(spm8w_osbabel(sprintf('mv "%s/"%s "%s/"%s',...
%                         spm_str_manip(subj_c1{i},'h'),['smw',spm_str_manip(subj_c1{i},'t')],...
%                         spm_str_manip(subj_c1{i},'h'),['smDARTEL',spm_str_manip(subj_c1{i},'t')]))); 
%                 [s,w] = system(spm8w_osbabel(sprintf('mv "%s/"%s "%s/"%s',...
%                         spm_str_manip(subj_c2{i},'h'),['smw',spm_str_manip(subj_c2{i},'t')],...
%                         spm_str_manip(subj_c2{i},'h'),['smDARTEL',spm_str_manip(subj_c2{i},'t')])));                        
%             end  
        elseif strcmp(dart_type,'epi2mni')
            fprintf('=============Beginning DARTEL normalization of EPI to MNI for %s=============\n', tmpstr);       
            spm('Defaults','fmri')
            if strcmp(p.whichnorm,'DARTEL'); epinormtok = 'd'; end
            if strcmp(p.whichnorm,'VBM8'); epinormtok = 'v'; end 
            %%% Determine if all appropriate files exist.
            fprintf('Checking files...\n');
            for i = 1:length(subjects)
                p = spm8w_getp('P',subjects{i},para_file);   
                if strcmp(p.whichnorm,'DARTEL'); anatdir = p.anat; end 
                if strcmp(p.whichnorm,'VBM8'); anatdir = p.vbm8; end
                % Ensure coreg
                if ~p.coreg2epi
                    error(['Your parameter file indicates that the anatomical was not coregistered',...
                        ' to the EPI (p.coreg2epi = 0), template is therefore invalid'])
                end                   
                % Determine boldtok    
                boldtok = []; 
                if p.despike;   boldtok = ['k',boldtok]; end
                if p.slicetime; boldtok = ['a',boldtok]; end
                if p.unwarp;    boldtok = ['u',boldtok]; end
                if ~p.unwarp;   boldtok = ['r',boldtok]; end
                boldtok = [boldtok, p.bold];
                % populate file location vars
                dartflow{i}  = fullfile(anatdir,flowfile); 
                for ii = 1:p.nses
                    subj_bold{ii}{i} = fullfile(p.func,[boldtok,num2str(ii),'.nii']);
                end
                % check if files exist.
                if ~exist(subj_bold{1}{i},'file') 
                     error('Subject %s preprocessed bold files are missing, did you preprocess this subject?',subjects{i});              
                end
                if ~exist(dartflow{i},'file')
                    error('Subject %s flow field is missing, did you generate a template with this subject?',subjects{i});
                end             
                %%%Delete any leftover or prior dartel files                 
                if exist(fullfile(p.func,['s',epinormtok,boldtok,'1.nii']),'file')
                    fprintf('Prior dartel normalization found... deleting files and proceeding...');  
                    try
                        [s,w] = system(spm8w_osbabel(['rm ',fullfile(p.func,['s',epinormtok,boldtok,'*.*'])]));  
                    catch
                    end
                end        
            end
            fprintf('Done...\n') 
            %%%Start time
            start = datestr(now);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Get file locations and setup JOB structure
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            %Setup job structure. If Dartel we need to align template to MNI
            %If VBM8 the template is already aligned (i.e. IXI_MN152).
            if strcmp(p.whichtemplate,'DARTEL'); job.template{1} = fullfile(p.dartdir,'Template_6.nii'); end
            if strcmp(p.whichtemplate,'VBM8'); job.template{1} = which('Template_6_IXI550_MNI152.nii'); end     
            %Setup images      
            fprintf('Normalizing %s runs with DARTEL:\n',boldtok);
            job.normtok  = epinormtok;
            fprintf('Normalizing functional volumes with DARTEL:');
            job.vox      = p.voxsize;
            job.bb       = p.boundbox;
            job.fwhm     = ones(1,3) * p.dartelsmooth;
            job.preserve = 0; %Preserve concentration for anat and epi
            for i = 1:length(subjects)
                job.data.subj(i).flowfield(1) = dartflow(i);
                job.data.subj(i).images       = subj_bold{1}(i);
                for ii = 2:p.nses
                    job.data.subj(i).images   = [job.data.subj(i).images; subj_bold{ii}(i)];
                end                
            end
            spm8w_dartel_norm_fun(job);         
        end             
        %Should we check toe delete the 2mni.mat in case of prior norm?
        %looks like it works except that the flow fields have 5 instead of
        %4 dimensions. Run again with SPMgui to double check.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Calculate Time and Save Log
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        stop = datestr(now); 
        p.time_elapsed = etime(datevec(stop),datevec(start)); %time elapsed in seconds
        [hours, minutes, seconds] = spm8w_timecalc(p.time_elapsed);
        fprintf('==============================================================\n');
        fprintf('DARTEL normalization for %s completed at %s\nand took %d hours, %d minutes and %d seconds...\n',tmpstr,stop,hours,minutes,seconds);
end
cd(p.root)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=============ADDITIONAL FUNCTIONS============
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%REPORT LAST ERROR MESSAGE AND DUMP USER BACK TO STUDYROOT
function error_reporter(rootdir, bbc_error, sparelines)
    cd(rootdir);
    fprintf('\n-----------------------------ERROR!-----------------------------\n%s\n',bbc_error);
    fprintf('----------------------------------------------------------------\n');
    err = lasterror;
    fprintf('\nLast error message: %s\n',err.message);
    for i = 1:length(err.stack)-sparelines
        fprintf('\nError occurred in function %s, at line %d of file %s\n',err.stack(i).name,err.stack(i).line,err.stack(i).file);
    end
return
