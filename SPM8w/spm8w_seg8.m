function spm8w_seg8(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: August, 2012 - DDW
% ==============================================================================
% spm8w_seg8('sub_id','parameter_file')
% 
% spm8w_seg8 will segment an MRI using SPM8's New Segment toolbox First argument
% is subject ID, second is optional parameter file (full path required). If 
% second argument is unspecified, a window will prompt for the file. 
%
% spm8w_seg8 will begin by coregistering anatomical to EPI if p.coreg2epi is
% set in the P file. This will estimate the registration but not reslice the 
% anatomical. After this, the anatomical will be segmented. Upon completion,
% segmented images can be found in the subject's anatomical MRI directory 
% (specified in the parameter file). 
%
% In addition a log file (subID_log_seg8.mat) will exist in subject's top-level
% dir, containing a structure, p, with relevant information regarding the 
% variables used for segmentation.
%
% Some tricks for calling new segment were borrowed from aamod_segment8 
% by Rhodri Cusack and Jonathan peelle.
%
% NOTE: For now using the spm8 TPM.nii. But could conceivably allow user to 
% set it to the VBM8 tissue probaility map. Albeit, that might not work.
% TODO: Figure out how to allow coregister to rest vs. coregister to
% functional. Also need a generalized coregister script to coreg functional
% to rest. 
% ==============================================================================
% CHANGE LOG:
% -First version -DDW August/12
% =======1=========2=========3=========4=========5=========6=========7=========8

%---Setup
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

%---Input checks
switch (nargin)
  case 0 
    error('Please specify a subject id, e.g. spm8w_preprocess(''01jan09ab'')');
  case 1
    p.subj = varargin{1};
    p      = spm8w_getp('P',p.subj);
    cd(p.root)
  case 2
    p.subj    = varargin{1};
    para_file = varargin{2};
    p         = spm8w_getp('P',p.subj,para_file);
    cd(p.root)
  otherwise
    error('You''ve specified too many inputs. Either 1 or 2 only please.');
end
  
try
  cd(fullfile(p.subdir,p.subj));
catch
    subjerror = sprintf(['I can''t seem to get to the subject''s dir.\nAre you sure ' ...
        'your subject dir (%s)\nand subject (%s) names are correct?'], p.subdir, p.subj);
    error_reporter(cwd, subjerror, 1);
    error('Exiting...');
end

%---Set Defaults
spm('Defaults','fmri')
%--spm('Defaults','fmri').m makes defaults global, but for some reason this 
%--doesn't stick by the time we get to unwarping unless we make defaults
%--global from within spm8w_preprocess.m -DDW Apr/10
global defaults
%--String to save log information
saver=['save ',p.subj,'_log_seg8 p;']; 
%--Bold and mprage token
bold   = p.bold;
mprage = p.mprage;
%--Add current date to p structure
p.start = datestr(now);

%---Cleanup previous preprocessing and checkout a fresh copy
if exist(p.anat, 'dir')
    try
        delete(sprintf('%s/%s',p.anat,'mprage.nii')) %in case coreg put new transforms
    end
    if exist([p.anat,'/c1',mprage,'.nii'], 'file') | exist([p.anat,'/',mprage,'_coreg.nii'], 'file')
      fprintf(['Previous segmentati and/or coregistration found... spm8w will delete the prior\n',...
              'segmentation/coregistration files and start over.\n'])
      if (~p.overridedel)
        dodelete = input('Are you sure you want to delete prior segmentation/coregistration? Y/N [Y]:','s');
        if isempty(dodelete)
          dodelete = 'Y';
        end
        if strcmp(lower(dodelete),'n')
          error('Cancelling segmentation...')
        end
      end
      %--Delete prior prepro
      fprintf('Deleting prior segmentation and/or coregistration...\n')
      try
        delete(sprintf('%s/%s',p.anat,'c*.nii'))
        delete(sprintf('%s/%s',p.anat,['mean',bold,'1.nii']))
        delete(sprintf('%s/%s',p.anat,'rc*.nii'))
        delete(sprintf('%s/%s',p.anat,'y_*.nii'))
        delete(sprintf('%s/%s',p.anat,'*_seg8.mat'))
        delete(sprintf('%s/%s',p.anat,'*_coreg.nii'))
      catch
      end
    end
else
    mkdir(p.anat)
end

%---Checkout fresh files from NIFTI dir
fprintf('Checking out fresh copy of anatomical MRI from:\n%s...',p.nifti)
gunzip([p.nifti,'/mprage*.gz'],p.anat);
fprintf('Done...\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Coregister Anatomical to Mean EPI image for Run 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To use DARTEL on EPI images it's important that Anatomicals and EPI
%are in alignment prior to segmentation and calculation of flow fields.
%Because we should limit reslicing on the EPI data we'll coregister the 
%anatomical to EPI first.
if(p.coreg2epi)
    fprintf('=============Coregistering Anatomical to Mean EPI for %s=============\n', p.subj);
    %Assume session 1 hasn't been realigned (because we don't know at what
    %stage user decided to do segmentation). Realign, generate mean EPI. 
    fprintf('Realigning to generate mean EPI from run 1...\n');
    copyfile(fullfile(p.func,[bold,'1.nii']),p.anat);
    %Read in functional files
    P{1} = spm_select('FPList','ANATOMY',['^',bold,'1.*\.nii']);
    realign_rtm = 0;        %Just need 1-pass realignment
    realign_flags = struct('rtm',realign_rtm, 'wrap', p.wrap_r, 'interp', p.interp_r); %-DDW Nov/11
    spm8w_realign(P, realign_flags);
    fprintf('Writing mean EPI...\n');  
    flags = struct('which',[0 1]); %Sets output to no reslice (0) output mean (1)
    spm_reslice(P,flags);    
    fprintf('Deleting intermediate files...\n'); 
    delete(sprintf('%s/%s*.*',p.anat,bold))
    delete(sprintf('%s/rp_%s*.*',p.anat,bold))
    %Coregister Anatomy to run1 mean bold.  
    fprintf('Coregistering Anatomical to Mean EPI...\n'); 
    %Make a copy of the MPRAGE
    copyfile(fullfile(p.anat,[mprage,'.nii']),fullfile(p.anat,[mprage,'_coreg.nii']))
    T = spm_vol(spm_select('FPList','ANATOMY',['^','mean',bold,'1.*\.nii']));
    S = spm_vol(spm_select('FPList','ANATOMY',['^',mprage,'_coreg.*\.nii']));   
    coreg_flags          = spm_get_defaults('coreg.estimate');
    coreg_flags.params   = [0 0 0  0 0 0];
    coreg_flags.graphics = spm('CmdLine');
    x = spm_coreg(T,S,coreg_flags);   
    M   = spm_matrix(x);           %convert coreg trans to matrix
    MMS = spm_get_space(S.fname);  %Get current image space
    spm_get_space(S.fname, M\MMS); %Set coreg matrix to sform
    %nifti sform now has transformation stored.
    mprage = [mprage, '_coreg'];
    fprintf('Coregistration transformation saved to header of %s.\n', [mprage,'.nii'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get file locations and setup JOB structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('=============Beginning Segmentation on %s=============\n', p.subj);
addpath(fullfile(spm('dir'),'toolbox','Seg'));    %Put the new seg toolbox files in path
addpath(fullfile(spm('dir'),'toolbox','DARTEL')); %Put DARTEL in path since newseg needs optimNn
tpm      = which('TPM.nii');                      %Get the TPM location
if size(tpm,1) ~= 1                               %Make sure only one TPM location
    error('More than one tissue probability map (TPM.nii) in path!')
end
anat_img = fullfile(p.anat,[mprage,'.nii']);  %Get anat image location
if ~exist(anat_img,'file')                      %Make sure file exists
    error('Unable to find %s.nii in %s',mprage,p.anat)
end

%User editable settings
job.channel(1).biasreg = p.biasreg;   %recommends increasing reg if good MRI intensity across volume
job.warp.mrf           = p.mrf;       %JA recommends an MRF of 2 for cleaning up segmentations
job.warp.reg           = p.warpreg;   %JA suggests decreasing warpreg for T1s by about a factor of 10 (but check QA)
job.warp.samp          = p.sampling;  %JA suggests increasing sampling for slightly better seg but at a big cost to cpu time

%Setup job using segmentation defaults (from SPM8)
job.channel(1).vols{1}  = anat_img;
job.channel(1).biasfwhm = 60;          %If intensity bias is not smooth then decrease, if very smooth then increase
job.channel(1).write    = [0,0];
job.channel(1).tpm      = tpm;
%Setup tissues
ngaus   = [2 2 2 3 4 2];                         %Num gaus for that tissue
%nval = [Native, DARTEL imported] I believe that native is the segmentation in
%native space whereas Dartel is aligned to the TPM and lower resolution but
%note that we bump the resolution back up to one using job.warp.vox. Cause we're crazy like that.
nval    = {[1 1],[1 1],[1 0],[0 0],[1 0],[0 0]}; 
%nwarped = Write out images for VBM [modulated, unmodulated] 
%Note that JA's tutorial says not to export them for VBM as we do this at a
%later stage (Normalise to MNI preserve amount)
nwarped = {[0 0],[0 0],[0 0],[0 0],[0 0],[0 0]}; 
for i = 1:numel(ngaus)
    tissue(i).tpm    = [tpm ',' num2str(i)];
    tissue(i).ngaus  = ngaus(i);
    tissue(i).native = nval{i};
    tissue(i).warped = nwarped{i};
end
job.tissue = tissue;
%Setup warp
job.warp.affreg         = 'mni'; %Easter european brains
job.warp.vox            = 1;     %This does nothing!
job.warp.write          = [0 1]; %Write forward transform image, if ever need inverse as well, then [1 1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run New Segment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
segout = spm8w_preproc_run(job);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Delete intermediate and unncessary files (ARE THERE ANY FOR SEG WE
%%% DON'T NEED? I won't know until I figure out how to analyse this
%%% stuff!)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% switch(p.cleanup)
%     case 0
%         fprintf('No cleanup, all preprocessed files will remain in %s/FUNCTIONAL\n\n',p.subj);
%     case 1
%         fprintf('Cleaning up.... deleting all but bold, swuabold/swubold and wuabold/wubold\n\n');
%         boldtok = bold(1:strfind(bold,p.bold)-1);
%         for i = 3:length(boldtok)
%             [s,w] = system(spm8w_osbabel(sprintf('rm FUNCTIONAL/%sbold*',boldtok(i:end))));
%         end     
%     case 2
%          fprintf('Cleaning up.... deleting all but bold and swuabold/swubold\n\n');
%          boldtok = bold(1:strfind(bold,p.bold)-1);
%          for i = 2:length(boldtok)
%             [s,w] = system(spm8w_osbabel(sprintf('rm FUNCTIONAL/%sbold*',boldtok(i:end))));
%          end  
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate Time and Save Log
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.stop = datestr(now); 
p.time_elapsed = etime(datevec(p.stop),datevec(p.start)); %time elapsed in seconds
eval(saver);
[hours, minutes, seconds] = spm8w_timecalc(p.time_elapsed);
fprintf('==============================================================\n');
fprintf('Subject %s segmentation finished at %s\nand took %d hours, %d minutes and %d seconds...\n',p.subj,p.stop,hours,minutes,seconds);
fprintf('A log file of all the session variables was saved to %s_log_seg8.mat...\n',p.subj);
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
