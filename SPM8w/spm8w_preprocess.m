function spm8w_preprocess(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: March, 2006
% ==============================================================================
% spm8w_preprocess('sub_id','parameter_file')
%
% spm8w_preprocess will perform image preprocessing for anatomical and
% functional imaging data. The first argument is a subjectID. The second 
% argument is the full path the a parameter files (P_studyname.m) and is 
% optional. If the second argument is unspecified, matlab will allow you to
% select the file.
%
% Preprocessing output will be located in the subject's FUNCTIONAL/REST
% and ANATOMICAL diretories depending on the preprocessing steps chosen.
%
% A log file (subID_log.mat) will exist in the subject's top-level
% dir, containing a structure (p) with relevant information regarding when and
% what happened during preprocessing.
%
% A PDF file (subID_prepro.pdf) will also be output to the subject's top-level 
% directory and contains matlab figure outputs of individual preprocessing
% steps.
% ==============================================================================
% CHANGE LOG:
% -Made some adjustments to figures to ensure realignment params get printed
%  when slice_timing is turned off. - DDW Jan/08
% -Added a fix for when runs are greater than 6 to allow for additional
%  figures to get added to the spm8.ps - DDW Nov/09
% -Converted spm2_preprocess to spm8 and made further adjustments to naming
%  conventions, shuffle check etc. - DDW Dec/09
% -Small adjustment to bounding box to fix spm8 origin problem - DDW Jan/10
% -Added sliceorder spec to the P file, rather than hardcoded here - DDW Jan/10
% -Made changes throughout this and spm8 functions to remove
%  progress bars and only spawn figures when necessary for printing - DDW Feb/10
% -Housecleaning - DDW Apr/10
% -Realign defaults updated to spm norms (rtm is 1 instead of 0). - DDW Nov/11
% -Fixed issue where mean image not written if no unwarping - DDW Nov/11
% -Pushed some defaults to the P file (p.realign_rtm, p.wrap, p.interp,
%  p.voxsize, p.boundbox) - DDW Nov/11
% -Added 3dDespike and some diagnostic figures using code from slices_analyse. 
%  DDW March/12
% -Cleaned up code for opneing, making and closing figures. -DDW March/12
% -Windows compatible -DDW March/12
% -Fixed prior data cleanup (wasn't working on unix) -DDW August/12
% -Number of additions for resting state preprocessing and fixes for
%  flexibility -DDW August/12
% -Added code for "checking out" unpreprocessed nifti files form subject's NIFTI
%  dir. This is so that we always have a virgin set since spm8 will ocassionally 
%  re-write the headers, even for raw bold1.nii files. -DDW February/13 
% -Fixed a pernicious bug whereby calling spm_write_sn with a argout was 
% causing spm to zero the scaling info (v.pinfo = [1 0] on line 239 of 
% spm_write_sn.m. The fix was to not call it with an ouput argument which we 
% were doing due to legacy pre-nifti code.- DDW March/13
% -Added scale factor adjustments to snr and slice noise calculation by
% pulling the scale factor from the private.dat.scl_slope variable in the
% nifti header.
% -Added support for different number of TRs per session. -DDW March/13
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
%--printstr hack, since it's missing from SPM8, allows us to 
%--print our own figs and gets around the SPM8 "feature"
%--whereby spm_print brings up a gui window and awaits user input. 
defaults.printstr = ['print -dpsc2 -painters -append -noui tmp_prepro_',p.logsuffix,'.ps']; 
%--Spawn and hide the spm graphics window
F = spm_figure('Create','Graphics','SPM8w Preprocessing','off');
%--String to save log information
saver=['save ',p.subj,'_log_',p.logsuffix,' p;']; 
%--Bold token
bold = p.bold;
%--Add current date to p structure
p.start = datestr(now);
%--Save diary info
diary_name = [p.subj,'_diary_',p.logsuffix,'.txt']; %string for diary name
diary (diary_name);

%---Print info to cmd window
fprintf('=============Beginning Preprocessing on %s=============\n', p.subj);

%---Cleanup previous preprocessing and checkout a fresh copy
if exist(p.func,'dir')
  fprintf(['Previous FUNCTIONAL directory exists... spm8w will delete the directory\n',...
          'and check out fresh copies from your nifti dir.\n'])
  fprintf('Any outlier txt files you had will be preserved...\n');
  if (~p.overridedel)
    dodelete = input('Are you sure you want to delete prior preprocessing directory? Y/N [Y]:','s');
    if isempty(dodelete)
      dodelete = 'Y';
    end
    if strcmp(lower(dodelete),'n')
      error('Cancelling preprocessing...')
    end
  end
  %--Delete prior prepro
  fprintf('Deleting prior functional preprocessing files...\n')
  try delete([p.subj,'*_',p.logsuffix,'.*']); end
  try delete(['*_',p.logsuffix]); end
  try movefile([p.func,'/outliers*.txt'],'./'); end
  try rmdir(p.func,'s'); end
  try mkdir(p.func); end
  try movefile('./outliers*.txt',p.func); end
end

%---Checkout fresh files from NIFTI dir
fprintf('Checking out fresh copies of bold files from:\n%s...',p.nifti)
gunzip([p.nifti,'/bold*.gz'],p.func)
fprintf('Done...\n')

%---Shuffle Check
if(p.shuffle)
  start_time = datestr(now);
  fprintf('==========Performing shuffle check on %s at %s\n', p.subj,start_time); 
  for ses = 1:p.nses
    p.shuff{ses} = spm8w_shuffle_check(p.func,[bold,num2str(ses)],p.nTE);
  end
  %--Determine and print time it took to complete.
  time_elapsed              = etime(datevec(datestr(now)),datevec(start_time));
  [hours, minutes, seconds] = spm8w_timecalc(time_elapsed);
  fprintf('Shuffle check finished. Time to complete: %d hours, %d minutes and %d seconds ...\n\n',hours,minutes,seconds);
end

%---3dDespike using AFNI 
if(p.despike)
  %--Check OS, as AFNI and FSL not available on windows
  if(~ispc)
      start_time = datestr(now);
      fprintf('==========Performing 3dDespike on %s at %s\n', p.subj,start_time); 
      for ses=1:p.nses
         fprintf('Despiking run %d...', ses);
         %--remove previous files or afni panic
         [s,w]          = system(sprintf('rm %s',fullfile(p.func,['d',bold,num2str(ses),'.nii'])));
         [tmpa,doutput] = system(sprintf('3dDespike -nomask -prefix %s %s',...
                          fullfile(p.func,['k',bold,num2str(ses),'.nii']),...
                          fullfile(p.func,[bold,num2str(ses),'.nii'])));
         %--parse output
         tmploc = strfind(doutput(strfind(doutput, '++ FINAL:'):end),'edits [');
         sedit = [doutput(strfind(doutput, '++ FINAL:')+tmploc(1)+6:strfind(doutput, '++ FINAL:')+tmploc(1)+8),'%'];
         bedit = [doutput(strfind(doutput, '++ FINAL:')+tmploc(2)+6:strfind(doutput, '++ FINAL:')+tmploc(2)+8),'%'];
         %--output results
         fprintf('finished\nOut of all timepoints %s edited and %s were heavily edited.\n',sedit,bedit);
         %--Adjust for AFNI screwing the headers by copying bold hdr
         %--(this may not be necessary with non dartmouth nifti but can't hurt)
         system(sprintf('fslchfiletype NIFTI_PAIR %s',fullfile(p.func,[bold,num2str(ses)])));
         system(sprintf('fslchfiletype NIFTI_PAIR %s',fullfile(p.func,['d',bold,num2str(ses)])));
         system(sprintf('cp %s.hdr %s.hdr',fullfile(p.func,[bold,num2str(ses)]),fullfile(p.func,['k',bold,num2str(ses)])));
         system(sprintf('fslchfiletype NIFTI %s',fullfile(p.func,[bold,num2str(ses)])));
         system(sprintf('fslchfiletype NIFTI %s',fullfile(p.func,['k',bold,num2str(ses)]))); 
      end
      %--Set bold token to dbold
      bold = ['k',bold];
      %--Determine and print time it took to complete.
      eval(spm8w_osbabel(['!touch despiked_',p.logsuffix]));
      time_elapsed              = etime(datevec(datestr(now)),datevec(start_time));
      [hours, minutes, seconds] = spm8w_timecalc(time_elapsed);
      fprintf('3dDespike finished. Time to complete: %d hours, %d minutes and %d seconds ...\n\n',hours,minutes,seconds);
  else
      fprintf('==========Your platform (Win) does not support AFNI and FSL, skipping 3dDespike\n'); 
  end
end

%---Slice Time Correct EPI images
if(p.slicetime)
  start_time = datestr(now);   
  fprintf('==========Performing slice time correction on %s at %s\n', p.subj,start_time); 
  %--read in functional files
  P = cell(1,p.nses);
  for ses = 1:p.nses
    P{ses} = spm_select('FPList',p.func,['^',bold,num2str(ses),'.*\.nii']);
  end  
  %--Slice timing correction
  fprintf('Slice acquisition order for this scan is: ');
  fprintf('%d ',p.sliceorder);
  fprintf('\nReference slice for slice timing is: %d', p.refslice);
  TA = p.TR-p.TR/p.nTE;    
  spm8w_slice_timing(P,p.sliceorder,p.refslice,[TA/(p.nTE-1) p.TR-TA]);
  %--Set bold token to abold
  bold = ['a',bold];
  %--Determine and print time it took to complete.
  p.slicetimed = datestr(now);
  time_elapsed = etime(datevec(p.slicetimed),datevec(start_time));
  [hours, minutes, seconds] = spm8w_timecalc(time_elapsed);
  fprintf('Slice Time Correction finished. Time to complete: %d hours, %d minutes and %d seconds ...\n\n',hours,minutes,seconds);
  eval(saver);
  eval(spm8w_osbabel(['!touch slicetime_',p.logsuffix]));
end

%---Realign EPI images
if(p.realign)
  start_time = datestr(now);   
  fprintf('==========Performing realignment on %s at %s\n', p.subj,start_time);
  %--Read in functional files
  P = cell(1,p.nses);
  for ses = 1:p.nses
    P{ses} = spm_select('FPList',p.func,['^',bold,num2str(ses),'.*\.nii']);
  end
  realign_flags = struct('rtm',p.realign_rtm, 'wrap', p.wrap_r, 'interp', p.interp_r); %-DDW Nov/11
  spm8w_realign(P, realign_flags);
  %--No need to write our realigned data as we'll use the transformation matrix at the next 
  %--reslice (i.e. unwarping or normalization). Might need to be careful when it comes to Dartel
  %--It should be able to pick up on the .mat files, but cleanup could
  %--harm the realignment if mat files are deleted prior to Dartel!!!!
  %--Print figure
  eval(defaults.printstr);
  %--Create mean image if no unwarping is selected (otherwise unwarp does it). 
  if(~p.unwarp) && (p.normalize)  %DDW - Sep/12
    flags = struct('which',[0 1]); %Sets output to no reslice (0) output mean (1)
    spm_reslice(P,flags);
  elseif (~p.unwarp) && (~p.normalize)
    flags = struct('which',[1 1]); %If not unwarping and not normalizing, better reslice here!
    fprintf('Writing out realigned images...\n');
    spm_reslice(P,flags); 
  end
  %--Determine and print time it took to complete.
  p.realigned  = datestr(now);
  time_elapsed = etime(datevec(p.realigned),datevec(start_time));
  [hours, minutes, seconds] = spm8w_timecalc(time_elapsed);
  fprintf('\nRealignment finished. Time to complete: %d hours, %d minutes and %d seconds ...\n\n',hours,minutes,seconds);
  eval(saver);
  %--Make shufflecheck figure
  if (p.shuffle)
      F=spm_figure('GetWin','Graphics'); set(F,'visible','off'); spm_figure('Clear',F);
      ys = (1:p.nTE)*0.25+4;  %yticks
      %-Layout depends on number of sessions (big for 3, small for >3)
      subplot_ses = 1;
      for ses=1:p.nses
        ra = load([p.func,'/rp_',bold,num2str(ses),'.txt']); 
        p.ra{ses} = ra;
        if p.nses <= 2
            subplot(2,1,subplot_ses); %changed ses to subplot_ses
        elseif p.nses == 3
            subplot(3,1,subplot_ses); %changed ses to subplot_ses
        else 
            subplot(3,2,subplot_ses); %changed ses to subplot_ses
        end
        imagesc(1:p.nTR(ses),ys,p.shuff{ses});
        set(gca,'YDir','Normal');
        hold on
        colormap jet
        dra   = diff(ra(:,1:3));
        [x,y] = find(abs(dra)>1);
        plot(ra-repmat(ra(1,:),p.nTR(ses),1));
        for j=1:length(x)
            plot([x(j),x(j)],[-5,20],'k--');
        end
        axis([0,p.nTR(ses),-4,max(ys)]);
        title([p.subj,' ',bold,' ',num2str(ses)]);
        %FIX to allow for additional pages if more than 6 runs - DDW Nov 2009
        if ses == 6 || ses == 12 || ses == 18 || ses == 24 || ses == 30
            eval(defaults.printstr)
            spm_figure('clear',F)  
            subplot_ses = 1;  %reset the subplot
        else
            subplot_ses = subplot_ses + 1;
        end 
      end
  end
  %--Set bold token to rbold if realigned images were written
  if (~p.unwarp) && (~p.normalize)
    bold = ['r',bold];
  end
  %--write figure to ps file
  eval(defaults.printstr);
  eval(spm8w_osbabel(['!touch realigned_',p.logsuffix]));
 end

%---Unwarp EPI images
if(p.unwarp)
  start_time = datestr(now);   
  fprintf('==========Performing unwarping on %s at %s\n', p.subj,start_time);
  %--read in functional files
  P = cell(1,p.nses);
  for ses = 1:p.nses
    P{ses} = spm_select('FPList',p.func,['^',bold,num2str(ses),'.*\.nii']);
  end
  %--Previously we set some defaults for uwe_flags, but they just
  %--duplicated the spm defaults, so I've removed them now that the
  %--the latest spm_uw_estimate grabs things from spm_defaults.m 
  %--DDW Nov/11            
  %--Estimate unwarping parameters
  for ses=1:p.nses
    fprintf('Estimating unwarping parameters for run %d\n',ses);
    Ptmp                = spm_vol(P{ses});
    uwe_flags.M         = Ptmp(1).mat;
    ds                  = spm8w_uw_estimate(P{ses});   %% this version has no graphics
    ads(ses)            = ds;
    [path,name]         = fileparts(P{ses});
    pefile              = fullfile(path,[name '_uw.mat']);
    save(pefile,'ds');
    fprintf('Estimation for run %d complete\n',ses);
  end
  %--Write unwarped images (more efficient if spm_uw_apply returned VO)
  disp(['Writing out unwarped images at ',datestr(now),' ...']);
  uwr_flags = struct(  'wrap',        p.wrap_r,...     %Default is [0 0 0]
                       'interp',      p.interp_r);     %Default is 4  
  spm_uw_apply(ads,uwr_flags);
  %--Determine and print time it took to complete.
  p.unwarped   = datestr(now);
  time_elapsed = etime(datevec(p.unwarped),datevec(start_time));
  [hours, minutes, seconds] = spm8w_timecalc(time_elapsed);
  fprintf('Unwarping finished. Time to complete: %d hours, %d minutes and %d seconds ...\n\n',hours,minutes,seconds);
  eval(saver);
  eval(spm8w_osbabel(['!touch unwarped_',p.logsuffix]));
  bold =['u',bold];
end

%---Mode 1000 Normalization
if isfield(p,'gms1k') && p.gms1k == 1 || isfield(p,'gms1k') && p.gms1k == 2
  start_time = datestr(now);   
  fprintf('==========Normalizing data to a mode of 1,000 on %s at %s\n', p.subj,start_time);
  %--read in functional files
  P = cell(1,p.nses);
  for ses = 1:p.nses
    P{ses} = spm_select('FPList',p.func,['^',bold,num2str(ses),'.*\.nii']);
  end
  %%% Calculate Mean or Media
  %%% Notes for myself: spm_global calculates the global mean over a 3D
  %%% image. It's in mex/c code unfortunately so can't view m source.
  %%% Reverse engineering it: For a given image, spm_global gives 180.04
  %%% It first calculates the mean (including zeros) then calcualtes
  %%% it again by taking only elements > mean/8. We can replicate this
  %%% function using mean(z1(z1>mean2(z1)/8)) or mean(z1(z1>mean(z1(:))/8) 
  %%% which gives 180.04. So that's what the global does, the mode would 
  %%% work the same way. (by the way this works for 3d or 4d). But be
  %%% careful when claculating variance. var of z1(:) can't use that trick.
  
%Ptmp = P{1};
%img = spm_read_vols(spm_vol(Ptmp));
%   
%   
%   z1   = img(:,:,:,1)
%   z2   = img(:,:,:,150)
%   
%   
%   spm_global(z1) = 180.0403
%   mean of z1(z1~=0) = 109.7267
%   mean(z1(z1>109.7267/8)) = 211.8839 %too big
%   mean2(z1) = 51.6278 %includes zeros
%   mean(z1(z1>51.6278/8)) = 180.0403 %so this is what spm_global does!
   %mean 2 will mean entire matrix
 
  %--Determine and print time it took to complete.
  p.unwarped   = datestr(now);
  time_elapsed = etime(datevec(p.unwarped),datevec(start_time));
  [hours, minutes, seconds] = spm8w_timecalc(time_elapsed);
  fprintf('Mode 1,000 normalisation finished. Time to complete: %d hours, %d minutes and %d seconds ...\n\n',hours,minutes,seconds);
  eval(saver);
  eval(spm8w_osbabel(['!touch mode1k_',p.logsuffix]));
  bold =['m',bold];
end

%---Normalise EPI images
if(p.normalize)
    start_time = datestr(now);   
    fprintf('==========Performing normalisation on %s at %s\n', p.subj,start_time); 
    %--Reload unwarped images
    P = cell(1,p.nses);
    for ses = 1:p.nses
        P{ses} = spm_select('FPList',p.func,['^',bold,num2str(ses),'.*\.nii']);
    end
    V = spm_vol(P);
    V = cat(1,V{:});
    %--Memory map the mean image 
    meanf   = [spm_str_manip(P{1}(1,:),'h') '/mean',...
              spm_str_manip(P{1}(1,:),'t')];
    Vm      = spm_vol(meanf);
    %--Normalization
    defaults.normalise.write.vox    = p.voxsize;
    defaults.normalise.write.bb     = p.boundbox;
    defaults.normalise.write.interp = p.interp_w;
    defaults.normalise.write.wrap   = p.wrap_w;
    matname  = [spm_str_manip(Vm.fname,'sd') '_sn.mat'];
    VG       = fullfile(spm('Dir'),'templates','EPI.nii');
    disp('Determining normalisation parameters ...');
    %don't need to add defaults.normalise.estimate unless we decide to
    %change it
    params   = spm8w_normalise(VG,Vm,matname,'','',defaults.normalise.estimate); 
    %--Print and close
    F=spm_figure('GetWin','Graphics'); set(F,'visible','off');       
    eval(defaults.printstr);
    msk      = spm_write_sn(V,params,defaults.normalise.write,'mask');
    spm_write_sn(Vm,params,defaults.normalise.write,msk);	% write nmean
    disp(['Writing normalised images at ',datestr(now),' ...']);
    diagmsg = 'normalised:';
    %--Write normalised
    %--Previosuly we wrote the normalized vol one file at a time (a
    %throwback to img/hdr days?) but this causes an issue whereby spm8
    %will zero the scaling info if called with an output argument, 
    %the fix is simply to write it in one shot. 
    %OLD CODE:
    %for i = 1:length(V)
    %    V1 = spm_write_sn(V(i),params,defaults.normalise.write,msk);
    %    spm_write_vol(V1,V1.dat);
    %    fprintf([repmat('\b',1,(length(diagmsg) + 1 +length(num2str(i)))),'%s %d'], diagmsg, i); %tweaked the msg so that it doesn't take up so much screen space
    %end;
    %NEW CODE 2013:
    spm_write_sn(V,params,defaults.normalise.write,msk);
    %--Determine and print time it took to complete.
    p.normalised = datestr(now);
    time_elapsed = etime(datevec(p.normalised),datevec(start_time));
    [hours, minutes, seconds] = spm8w_timecalc(time_elapsed);
    fprintf('\nNormalization finished. Time to complete: %d hours, %d minutes and %d seconds ...\n\n',hours,minutes,seconds);
    eval(saver);
    eval(spm8w_osbabel(['!touch normalised_',p.logsuffix]));
    bold =['w',bold];
end

%---Smoothing
if(p.smoothing)
    start_time = datestr(now);   
    fprintf('==========Smoothing images on %s at %s\n', p.subj,start_time); 
    %--Reload normalized images
    fprintf('Loading files...\n');
    P = cell(1,p.nses);
    for ses = 1:p.nses
        P{ses} = spm_select('FPList',p.func,['^',bold,num2str(ses),'.*\.nii']);
    end
    V = spm_vol(P);
    fprintf('Smoothing session: ');
    for i=1:length(V)
        name_tmp = V{i}.fname;
        Q        =['s',spm_str_manip(name_tmp,'t')]; % append 's' to smoothed filenames.
        fprintf([repmat('\b',1,1),'%1.0f'],i);  %DDW
        spm_smooth([spm_str_manip(name_tmp,'h'),'/',spm_str_manip(name_tmp,'t')],[spm_str_manip(name_tmp,'h'),'/',Q],p.smooth_kernel)    
    end
    %--Determine and print time it took to complete.
    p.smoothed   = datestr(now);
    time_elapsed = etime(datevec(p.smoothed),datevec(start_time));
    [hours, minutes, seconds] = spm8w_timecalc(time_elapsed);
    fprintf('\nSmoothing finished. Time to complete: %d hours, %d minutes and %d seconds ...\n\n',hours,minutes,seconds);
    eval(saver);
    eval(spm8w_osbabel(['!touch smoothed_',p.logsuffix])); 
    bold =['s',bold];
end;

%---SNR 
if(p.snr)
    %--Turn off warning to avoid divide by zero warning
    warning off
    %--Print start info
    start_time = datestr(now);   
    fprintf('==========Calculating SNR on %s at %s', p.subj,start_time); 
    %--Reload smoothed images
    P = cell(1,p.nses);
    for ses = 1:p.nses
        P{ses} = spm_select('FPList',p.func,['^',bold,num2str(ses),'.*\.nii']);
    end
    V = spm_vol(P);
    for ses = 1:p.nses
        files  = V{ses};
        scalf  = files(1).private.dat.scl_slope;
        if scalf > 0
            fprintf('\nData will be scaled by %.2f prior to SNR calcualtion...\n', scalf);
        end
        data   = spm_read_vols(files(1));
        avg    = zeros(size(data));
        sd_tmp = zeros(size(data));
        n      = size(files,1);
        fprintf('\n');
        for i = 1:n
             diagmsg = 'Calculating average signal on volume:';
             fprintf([repmat('\b',1,(length(diagmsg) + 1 +length(num2str(i)))),'%s %d'], diagmsg, i);             
             data    = spm_read_vols(files(i))*scalf; %scale data prior to calc
             avg     = avg+data/n;
        end
        fprintf('\n');
        for i = 1:n
             diagmsg = 'Calculating standard deviation on volume:';
             fprintf([repmat('\b',1,(length(diagmsg) + 1 +length(num2str(i)))),'%s %d'], diagmsg, i);             
             data    = spm_read_vols(files(i))/scalf; %scale data prior to calc
             sd_tmp  = sd_tmp+(avg-data).^2;  
        end
        fprintf('\n');
        sd = sqrt(sd_tmp/(n-1));
        snr = avg./sd;
        %-Clean up snr varaible
        snr(isnan(snr)) = 0;
        snr(snr>5000)   = 0; 
        %-Make sure SNR directories exists
        try
            mkdir(fullfile(p.root,['/SNR_',upper(p.logsuffix)]));
            mkdir(fullfile(p.root,['/SNR_',upper(p.logsuffix),'/AVG']));
            mkdir(fullfile(p.root,['/SNR_',upper(p.logsuffix),'/SD']));
            mkdir(fullfile(p.root,['/SNR_',upper(p.logsuffix),'/SNR']));
        end
        %-Output volumes
        fprintf('Writing out volumes...');
        vol_out       = files(1);
        vol_out.pinfo = [1 0]'; %set scaling to 1 (we don't want the original data's scaling)
        vol_out.fname = fullfile(p.root,['SNR_',upper(p.logsuffix),'/AVG'],[p.subj,'_avg_run',num2str(ses),'.nii']);
        spm_write_vol(vol_out, avg);
        vol_out.fname = fullfile(p.root,['SNR_',upper(p.logsuffix),'/SD'],[p.subj,'_sd_run',num2str(ses),'.nii']);
        spm_write_vol(vol_out, sd);
        vol_out.fname = fullfile(p.root,['SNR_',upper(p.logsuffix),'/SNR'],[p.subj,'_snr_run',num2str(ses),'.nii']);
        spm_write_vol(vol_out, snr);
    end
    %--Make figures
    F=spm_figure('GetWin','Graphics'); set(F,'visible','off'); spm_figure('Clear',F);
    currentplot = 1;
    for ses=1:p.nses
        %-load files
        files_avg   = spm_vol(fullfile(p.root,['SNR_',upper(p.logsuffix),'/AVG'],[p.subj,'_avg_run',num2str(ses),'.nii']));
        files_sd    = spm_vol(fullfile(p.root,['SNR_',upper(p.logsuffix),'/SD'],[p.subj,'_sd_run',num2str(ses),'.nii']));
        files_snr   = spm_vol(fullfile(p.root,['SNR_',upper(p.logsuffix),'/SNR'],[p.subj,'_snr_run',num2str(ses),'.nii']));
        data_avg    = spm_read_vols(files_avg);
        data_sd     = spm_read_vols(files_sd);
        data_snr    = spm_read_vols(files_snr);
        % Scale data for display purposes /  different scanners 
	scaleavg = max(data_avg(:))/900;
	data_avg = data_avg/scaleavg;
	data_sd = data_sd/scaleavg;       
        %-Slices
        slice{1} = squeeze(data_avg(:,:,20));   stitles{1} = 'Average 1';   sthresh{1} = [10,900];
        slice{2} = squeeze(data_avg(:,:,24));   stitles{2} = 'Average 2';   sthresh{2} = [10,900];
        slice{3} = squeeze(data_sd(:,:,20));    stitles{3} = 'Stdev 1';     sthresh{3} = [2,20];
        slice{4} = squeeze(data_sd(:,:,24));    stitles{4} = 'Stdev 2';     sthresh{4} = [2,20];
        slice{5} = squeeze(data_snr(:,:,20));   stitles{5} = 'SNR 1';       sthresh{5} = [10,350];
        slice{6} = squeeze(data_snr(:,:,24));   stitles{6} = 'SNR 1';       sthresh{6} = [10,350];
        slice{7} = squeeze(data_snr(26,:,:));   stitles{7} = 'SNR Sagital'; sthresh{7} = [10,350];
        slice{8} = squeeze(data_snr(:,32,:));   stitles{8} = 'SNR Coronal'; sthresh{8} = [10,350];
        %-Plot defaults
        s_xpos = [0.05,0.26,0.54,0.75,0.05,0.26,0.54,0.75];
        s_tpos = [27,27,27,27,27,27,33,27];
        colormap hot
        %-Set Currentplot defaults (1 = TOP)
        if currentplot == 1 %Top plot
            t_ypos = 0.75;
            splot = 1;
            s_ypos = [0.75,0.75,0.75,0.75,0.55,0.55,0.55,0.55];
            
        else                %Bottom plot
            t_ypos = 0.27;
            splot = 9;
            s_ypos = [0.27,0.27,0.27,0.27,0.07,0.07,0.07,0.07];    
        end
        %-Plots - 8 per halfpage. So subplot is 4,4, and 1 to 16.
        for iplot = 1:8
            subplot(4,4,splot)
                set(gca,'position',[s_xpos(iplot),s_ypos(iplot),0.2,0.2])
                imagesc(flipud(slice{iplot}'),sthresh{iplot})
                axis equal 
                axis off
                title(stitles{iplot},'fontweight','bold','position',[s_tpos(iplot),0.5])
                splot = splot + 1;
        end
        %Titles
        titlestr = ['SNR Subject: ',p.subj, ' Run: ',num2str(ses)];
        titleax = axes('Position',[0.1 t_ypos 0.8 0.2],'Parent',F,'Visible','off');
        set(get(titleax,'Title'),'String',titlestr,'FontSize',16,'FontWeight','Bold','Visible','on');
        %-Print Check
        if ses == p.nses
           eval(defaults.printstr);
           break;
        end            
        if currentplot == 1
            currentplot = 2;
        else
           eval(defaults.printstr);  
           spm_figure('clear',F);  
           currentplot = 1;
       end     
    end  
    %--Turn warning back on
    warning on
    %--Determine and print time it took to complete.
    p.snr        = datestr(now);
    time_elapsed = etime(datevec(p.snr),datevec(start_time));
    [hours, minutes, seconds] = spm8w_timecalc(time_elapsed);
    fprintf('\nSNR calculation finished. Time to complete: %d hours, %d minutes and %d seconds ...\n\n',hours,minutes,seconds);
    eval(saver);
    eval(spm8w_osbabel(['!touch snr_',p.logsuffix]));     
end;

%---Slice Noise Analysis 
%---borrowed from slices_analyse by Antonia Hamilton
if(p.slices)
      start_time = datestr(now);
      fprintf('==========Performing slice noise check on %s at %s\n', p.subj,start_time); 
      fprintf('Loading unsmoothed normalised volumes...\n');
      %--Check for unsmoothed data (in case user decided not to smooth)
      if strcmp(bold(1),'s')
        boldtok = bold(2:end);
      else
        boldtok = bold;
      end
      %--Load Images  
      P = cell(1,p.nses);
      for ses=1:p.nses
        P{ses} = spm_select('FPList',p.func,['^',boldtok,num2str(ses),'.*\.nii']);
        if p.slicetime == 1; rabold = 'abold'; else rabold = 'bold'; end
        ra = load([p.func,'/rp_',rabold,num2str(ses),'.txt']); 
        p.ra{ses} = ra;
      end
      %--Concat runs or not? For now I'm concating.  
      V      = spm_vol(P);
      V      = cat(1,V{:}); %Remove this and put the rest of code in for loop above for non concat runs
      ra     = cat(1,p.ra{:});
      nslice = V(1).dim(3);
      nscan  = length(V);
      scalf  = V(1).private.dat.scl_slope;
      if scalf > 1
          fprintf('Data will be scaled by %.2f prior to slice noise calcualtion...\n', scalf);
      end
      %--So much for using p.mask, we need bigmask_3x3x3 here instead of 1x1x1
      %--Ergo voila a dumb hack which will break if non bigmask people don't
      %--have an additional mask in 3x3x3
      mask_name = which(p.mask); 
      [maskpath, maskfile, maskext] = fileparts(mask_name);
      mask_name = fullfile(maskpath,[maskfile,'_3x3x3',maskext]);
      M         = spm_vol(mask_name);
      M.dat     = spm_read_vols(M);
      fprintf('Calculating slice wise noise'); 
      %--work out good neighbours
      neigh              = repmat(-5:2:5,nscan,1)+repmat(1:nscan,6,1)';
      neigh(neigh<1)     = NaN;
      neigh(neigh>nscan) = NaN;
      for kk=1:nscan      
        %-calculate neighbours to pull from for noise calc
        df      = abs(ra-repmat(ra(kk,:),nscan,1));
        trans   = sum(df(:,1:3),2);
        rot     = sum(df(:,4:6),2);
        good    = find(trans<1 & (rot*180/pi)<1);
        overlap = intersect(good,neigh(kk,:)); 
        if(length(overlap)<3)      %% really bad scans
          neigh(neigh==kk) = NaN;     
          neigh(kk,:)      = NaN;
        elseif(length(overlap<6))  %% moderate bad scans
          missing          = setdiff(neigh(kk,:),good);
          for i=1:length(missing)
            mind           = find(neigh(kk,:)==missing(i));
            neigh(kk,mind) = NaN;
          end
        end
      end
      %--Calculate slice noise
      for kk=1:nscan  %% for every input file
          gv(kk) = spm_global(V(kk));
          [dat1,loc] = spm_read_vols(V(kk),1);  %% read with zero masking   
          dat1 = dat1/scalf;                    %% Apply scaling
          nn = neigh(kk,:);
        for jj=1:length(nn)    %% read some neighbours
            if(isfinite(nn(jj)))  %% for good neighbours
                jind = nn(jj);
                [dat2,loc] = spm_read_vols(V(jind),1);  %% read with zero masking
                dat2 = dat2/scalf;                      %% Apply scaling
                for i=1:nslice
                    slice1 = squeeze(dat1(:,:,i));
                    slice2 = squeeze(dat2(:,:,i));
                    msk = squeeze(M.dat(:,:,i));
                    df = (slice1-slice2).*(msk>0.5);   %mask the data      
                    scan_noise(jj,i) = nanmean(abs(df(msk>0.5)));
                end
            else  %% for bad neighbours don't calc slicenoise
                scan_noise(jj,1:nslice) = NaN;
            end
        end
        %-Average the mean of the difference between current slice
        %-and all good slice neighbors
        noise(:,kk) = nanmean(scan_noise)'; 
        fprintf('.');
      end
      fprintf('\n');
      %--Make slices figure (now 200% more better with spm figurage)
      %--DDW march/12
      F=spm_figure('GetWin','Graphics'); set(F,'visible','on'); spm_figure('Clear',F);
      th  = 15;  %default was 20 from slices_defaults.m
      wth = 25;  %default was 30 lowering defaults since new coil -DDW Feb/12 
      
      colormap('default');
      subplot(3,1,1)
          imagesc(noise,[0,80])
          ttl = ['Subject: ',p.subj, ' File: ',boldtok];
          title(ttl);
      subplot(3,1,2)
          [n,b] = hist(noise(:),[0:80]);
          bar(b,n);
          set(gca,'xlim',[0,80])
          title('distribution of slice noise')
          hold on
          plot([th,th],[0,max(n)],'r-')
          h=text(th,max(n),[num2str(100*sum(noise(:)>th)./prod(size(noise)),3),...
                          '% of slices are over ',num2str(th)]);
          set(h,'Color',[1 0 0])
          plot([wth,wth],[0,max(n)/2],'g-')
          h=text(wth,max(n)/2,[num2str(100*sum(noise(:)>wth)./prod(size(noise)),3),...
                          '% of slices are over ',num2str(wth)]);
          set(h,'Color',[0 1 0])
      subplot(6,1,5)
          plot(ra(:,1:3))
          title('translation (mm)')
          set(gca,'XTick',[])
      subplot(6,1,6)
          plot(ra(:,4:6)*180/pi)
          title('rotation (deg)')
      %--Title
      titleax = axes('Position',[0.12 0.75 0.8 0.2],'Parent',F,'Visible','off');
      set(get(titleax,'Title'),'String','Slice Noise Analysis','FontSize',16,'FontWeight','Bold','Visible','on');
      %--Print
      eval(defaults.printstr);  
      %--Determine and print time it took to complete.
      time_elapsed              = etime(datevec(datestr(now)),datevec(start_time));
      [hours, minutes, seconds] = spm8w_timecalc(time_elapsed);
      fprintf('Slice noise check finished. Time to complete: %d hours, %d minutes and %d seconds ...\n\n',hours,minutes,seconds);
end

%---Delete intermediate and unncessary files
switch(p.cleanup)
    case 0
        fprintf('No cleanup, all preprocessed files will remain in %s\n\n',p.func);
    case 1
        fprintf('Cleaning up.... deleting all but bold, swuabold/swubold and uabold/ubold/abold\n\n');
        %--now need to be crafty to delete only appropriate files
        if p.unwarp == 1
            boldtok = 'ward'; %normalize, slicetime, realign and despike
        elseif p.realign == 1
            boldtok = 'wad';  %normalize, slicetime and despike
        else
            boldtok = 'wd';   %normalize and despike
        end
        for k = 1:p.nses
            for i = 1:length(boldtok)
                try
                    [s,w] = system(spm8w_osbabel(sprintf('rm "%s/"%s',p.func,[boldtok(i),'*',p.bold,num2str(k),'.nii'])));
                    [s,w] = system(spm8w_osbabel(sprintf('rm "%s/"%s',p.func,[boldtok(i),'*',p.bold,num2str(k),'*.mat'])));
                catch
                end
            end
        end
    case 2
         fprintf('Cleaning up.... deleting all but bold and swuabold/swubold\n\n');
         boldtok = bold(1:strfind(bold,p.bold)-1);
         for i = 2:length(boldtok)
            try
            [s,w] = system(spm8w_osbabel(sprintf('rm "%s/"%s*',p.func,[boldtok(i:end),p.bold])));
            catch
            end  
         end
end

%---Close hidden figure
spm_figure('Close',F); 

%---Convert multipage ps file to pdf 
%Switching from unix only ps2pdf to matlab ps2pdf.m which
%should be more cross-platform. -DDW Nov/11
ps2pdf('psfile',['tmp_prepro_',p.logsuffix,'.ps'],'pdffile',[p.subj,'_prepro_',p.logsuffix,'.pdf'],...
       'gspapersize','a4','deletepsfile',1);
%[s,w]=unix(['ps2pdf tmp_prepro.ps ',p.subj,'_prepro.pdf']);
%--delete flak.
[s,w] = system(spm8w_osbabel('rm tmp_*.*'));

%---Calculate Time and Save Log
p.stop = datestr(now); 
p.time_elapsed = etime(datevec(p.stop),datevec(p.start)); %time elapsed in seconds
eval(saver);
[hours, minutes, seconds] = spm8w_timecalc(p.time_elapsed);
fprintf('==============================================================\n');
fprintf('Subject %s preprocessing finished at %s\nand took %d hours, %d minutes and %d seconds...\n',p.subj,p.stop,hours,minutes,seconds);
fprintf('A log file of all the session variables was saved to %s_log_%s.mat...\n',p.subj,p.logsuffix);
fprintf('A txt file of the command window output was saved to %s_diary_%s.txt...\n',p.subj,p.logsuffix);
fprintf('A pdf file of the preprocessing figures was saved to %s_prepro_%s.pdf...\n',p.subj,p.logsuffix);
diary off;
cd(p.root)

%---ADDITIONAL FUNCTIONS

%---Report last error message and dump user back to studyroot
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
