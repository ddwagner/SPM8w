function spm8w_nifticonverter(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: March, 2013 - DDW
% Created: March, 2013 - DDW
% ==============================================================================
% spm8w_nifticonverter(parameter_file, subjects);
%
% spm8w_nifticonverter converts raw philips format PAR/REC into NIFTI files
% and organizes them into a subid/NIFTI directory ready for use with SPM8w.
% At present, spm8w_nifticonverter only works with raw PAR/REC files and thus
% is suited for Dartmouth data, however if people at other labs can tip me 
% off to their centre's raw data format, we can add functionality to convert 
% that data to NIFTI or simply to unzip RAW data and organize it into the
% appropriate folders. 
%
% spm8w_nifticonverter requires that your subject data exists inside a 
% raw_subID.tar.gz file which contains a directory of same name as subID and 
% inside that directory a series of PAR/REC data files. The location of these
% files is set in the Parameter file, but usually corresponds to:
% STUDYNAME/SUBJECTS/RAW. See example dataset for organization. 
%
% Raw data needs to be stored in a tar.gz file which can be made in unix/linux
% using the rawtar command in the BIN directory provided with SPM8w. Or by
% using matlab's tar and gzip commands (as implemented in spm8w_rawtar.m).  
%
% NOTES: Philips from 2011 stored data differently than 2012 namely, order of 
% slices volumes in the header are different. Thus neither dcm2nii nor Michel's 
% command line tool work. R2AGUI appears to do the job, but it is a little 
% unwiedly to use. spm8w_nifticonverter facilitates the process by organizing 
% the r2agui files into the appropriate directories (for SPM8w), 
% renaming them and deleting temporary files.
%
% For the moment this is designed to work with Dartmouth data, 
% although it can be extended to work with other data conversion
% schemes (nibabel, etc.). The primary purpose of spm8w_nifticonverter
% is to take the output of other converters and rename/arrange files
% so that they're ready to be processed with spm8w.
%
% To use, simply run spm8w_nifticonverter. Optional arguments are to provide
% the path to your parameter_file as well as providing the subjID. 
% ==============================================================================
% CHANGE LOG:
% -Batched up the old spm8w_r2agui.
% =======1=========2=========3=========4=========5=========6=========7=========8

%---Setup
%--Get the current working dir
cwd = pwd;
%--Load user defaults, otherwise load vanilla defaults
if exist(which('spm8w_userdefaults.m'), 'file')
    [def_dir,def_file] = fileparts(which('spm8w_userdefaults.m'));
else
    [def_dir,def_file] = fileparts(which('spm8w_defaults.m'));
end
%--Goto dir, eval defaults, check for errors, come back to cwd
cd(def_dir); 
try
    eval(def_file); 
catch
    defmsg = sprintf('The defaults file %s.m in %s is not evaluable, please check your syntax\n',def_file,def_dir); 
    error_reporter(cwd, defmsg, 1);
    error('Exiting...');
end
cd(cwd); 

%--Input checks
switch (nargin)
  case 0 
    p        = spm8w_getp('P');
    subjects = spm_select(Inf,'^raw_.*\.gz$',bbc.subjects,[],p.rawdir);
    cd(p.root) 
  case 1
    p        = spm8w_getp('P',[],varargin{1});
    subjects = spm_select(Inf,'^raw_.*\.gz$',bbc.subjects,[],p.rawdir);
    cd(p.root)
  case 2
    p        = spm8w_getp('P',[],varargin{1});
    subjects = varargin{2};
    cd(p.root)  
  otherwise
    error('You''ve specified too many inputs. Either 0, 1 or 2 only please.');
end

%---Check for epty subjects in case of abort.
if isempty(subjects) %Check for empty subjects in case user aborted
    error(bbc.error_subjects);
end

%---Check that rawdir exists
if ~exist(p.rawdir,'dir')
    error(['I can''t seem to get to the the raw data dir. Are you sure ' ...
        'your\nraw data dir:%s is correct?\n'], p.rawdir)
else
    cd(p.rawdir);
end

%---Clean up subjects list and cast into cell array
%--Remove dir and file extensions and raw_ prefix
subjects_tmp = spm_str_manip(subjects,'trr');
clear subjects
subjects(1:length(subjects_tmp(:,1)))={zeros}; %Preallocate subjects
for i = 1:size(subjects_tmp,1)
    [toktmp, subjects_tok] = strtok(subjects_tmp(i,:),'_');
    subjects{i}=deblank(subjects_tok(2:end));
end

%---R2AGUI OPTIONS
options.subaan = 1;
options.usealtfolder = 0;
options.altfolder = '~/';
options.prefix= [];
options.angulation= 1;
options.rescale= 1;
options.usefullprefix= 0;
options.outputformat= 1;
options.dim= 4;
options.dti_revertb0= 0;

%---RUN R2AGUI
for i = 1:length(subjects)
    cd(p.rawdir)
    subjpath        = fullfile(p.rawdir,subjects{i});
    options.pathpar = [subjpath,filesep];
    %--Extract the raw files (duplicating rawextract)
    fprintf('Extracting raw data from archive, this may take awhile...\n');
    untar(['raw_',subjects{i},'.tar.gz']);
    %--Find the raw PAR files for this subject
    files = dir(fullfile(subjpath,'*.PAR'));
    count = 1;
    for ii = 1:length(files)
        %--Prune scout scans
        if files(ii).bytes > 10000 %scouts tend to be 8078 bytes
            filelist{count} = files(ii).name;
            count = count + 1;
        end
    end
    fprintf('Converting...\n');
    convert_r2a(filelist,options);
    fprintf('Done.\n');
    %--Make various subject dirs
    %--First reload p file with subjID
    p = spm8w_getp('P',subjects{i},p.para_file);
    if exist(fullfile(p.subdir,subjects{i}),'dir')
        fprintf(['Previous Subject directory exists... spm8w will delete the directory\n',...
          'and convert from raw.\n'])
        if (~p.overridedel)
            dodelete = input('Are you sure you want to delete prior preprocessing directory? Y/N [Y]:','s');
            if isempty(dodelete)
                dodelete = 'Y';
            end
            if strcmp(lower(dodelete),'n')
                error('Cancelling nifticonverter...')
            end
        end
    %--Delete prior subject directory
        fprintf('Deleting prior subject directory...\n')
        rmdir(fullfile(p.subdir,subjects{i}),'s')
    end
    %--Make dirs, move files and rename them
    mkdir(fullfile(p.subdir,subjects{i}))
    mkdir(p.nifti)
    for i_flist = 1:size(filelist,2)
        [tmp,file_noext] = fileparts(filelist{i_flist});
        moveme = fullfile(p.rawdir,subjects{i},file_noext);
        movefile([moveme,filesep,'*.nii'],p.nifti);
    end    
    %--Anatomical moving and renaming
    mrifile = dir(fullfile(p.nifti,'*T1TFE*.nii'));
    if ~isempty(mrifile)
        for i_mri = 1:size(mrifile,1)
            if i_mri == 2
                mriname = 'mprage_halfandhalf.nii';
                fprintf('WARNING: Found a 2nd T1TFE file, assuming it is mprage half-and-half...\n');
            else
                mriname = 'mprage.nii';
            end
            fprintf('Renaming %s to %s and compressing...\n',mrifile(i_mri).name,mriname); 
            movefile(fullfile(p.nifti,mrifile(i_mri).name),fullfile(p.nifti,mriname));
            gzip(fullfile(p.nifti,mriname));
            delete(fullfile(p.nifti,mriname));
        end
    end
    %--DTI moving and renaming 
    dtifiles = dir(fullfile(p.nifti,'*DwiSE*.nii'));
    if ~isempty(dtifiles)
        fprintf('Renaming and compressing dti files...\n');
        dtirun = 1;
        lastgradnum = 0;
        for x = 1:length(dtifiles)
            %-Figure out the appropriate DTI name (2 runs plus the g
            %-numbers)
            gradnum = dtifiles(x).name;
            gradnum = gradnum(strfind(gradnum,'-g')+2:end-4);
            gradnum = str2num(gradnum);
            if gradnum < lastgradnum
                dtirun = dtirun + 1;
            end
            lastgradnum = gradnum;
            dtiname = sprintf('dti%d_g%02d.nii',dtirun,gradnum);
            movefile(fullfile(p.nifti,dtifiles(x).name),fullfile(p.nifti,dtiname));
            gzip(fullfile(p.nifti,dtiname));
            delete(fullfile(p.nifti,dtiname));
        end
    end
    %--BOLD moving and renaming
    boldfiles = dir(fullfile(p.nifti,'*FEEPI*.nii'));
    restrun = 1;
    boldrun = 1;
    %--jazzchk is a hacky way of finding out if we're dealing with the 
    %--example dataset. Unfotunately, I used resting state data in the 
    %--example dataset, so nifticonverter will convert it to rest.nii
    %--unless we stop it. This will ruin anyone who actually uses 
    %--H8TJAZZ in their study name. Hopefully users won't!
    jazzchk = sum(strfind(p.para_file,'H8TJAZZ')) > 1; 
    for x = 1:length(boldfiles)
        %-check for ress (rest filesize = 55296352)
        %this should always be the case with rest
        if boldfiles(x).bytes == 55296352 && jazzchk == 0 
            boldname = sprintf('rest%d.nii',restrun);
            restrun = restrun + 1;
        else
            boldname = sprintf('bold%d.nii',boldrun);
            boldrun = boldrun + 1;
        end
        fprintf('Renaming %s to %s and compressing...\n',boldfiles(x).name,boldname); 
        movefile(fullfile(p.nifti,boldfiles(x).name),fullfile(p.nifti,boldname));
        gzip(fullfile(p.nifti,boldname));
        delete(fullfile(p.nifti,boldname));
    end
    %--Delete the r2agui dirs
    rmdir(subjpath,'s');
    clear filelist newfiles files
end

%---Drop Info
fprintf('Conversion to nifti on %d subjects completed...\n', length(subjects));                
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


