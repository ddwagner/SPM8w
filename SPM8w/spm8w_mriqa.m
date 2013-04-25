function spm8w_mriqa(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: March, 2013 - DDW
% ==============================================================================
% spm8w_mriqa('subjects', 'parameter_file')
%
% spm8w_mriqa is a simple tool to generate a PDF of a single slice of an
% an antomical MRI in every plane. Its main purpose is to spot check
% MRI quality prior to segmentation (e.g. for DARTEL). 
%
% spm8w_mriqa takes a parameters file and a list of subjectIDs. If none are 
% specified, it will ask you for them. 
% A PDF file (studyname_mriqa.pdf) will be output to the study's root dir.
% ==============================================================================
% CHANGE LOG:
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
    subjects = spm8w_getsub;
    p        = spm8w_getp('P');
  case 1
    subjects = varargin{1};
    p        = spm8w_getp('P');
  case 2
    subjects = varargin{1};
    p        = spm8w_getp('ROI','',varargin{2});
  otherwise
    error('You''ve specified too many inputs. Either 0, 1 or 2 only please.');
end

%Goto root. 
cd(p.root);

%---Set Defaults
defaults.printstr = ['print -dpsc2 -painters -append -noui tmp_mriqa.ps']; 
qa_start = datestr(now);

%---Delete prior files
if exist(fullfile(p.root,'tmp_mriqa.ps'),'file')
    delete(fullfile(p.root,'tmp_mriqa.ps'));
end
if exist(fullfile(p.root,'MRI_QA.pdf'),'file')
    delete(fullfile(p.root,'MRI_QA.pdf'));
end

%---Print info to cmd window
fprintf('=============Beginning MRI QA at %s=============\n', qa_start);

%---Loop through subjects to get paths to anatomical MRIs
mri_list = {};    %List of paths to anatomical files
i_mri_list = 1;   %Counter for list
for i = 1:length(subjects)
    %--Reload p with subject.
    p = spm8w_getp('P', subjects{i}, p.para_file);
    %--Exist dir and mprage?
    if exist(p.anat,'dir') && exist(fullfile(p.anat,'mprage.nii'),'file')
        mri_list{i_mri_list} = fullfile(p.anat,'mprage.nii');
        i_mri_list = i_mri_list + 1;
    elseif exist(fullfile(p.nifti,'mprage.nii.gz'),'file')
        %---If not anatomical, check for and checkout a fresh one from NIFTI dir
        fprintf('Anatomical file not found, checking out fresh copy from:\n%s...',p.nifti)
        mkdir(p.anat)
        gunzip([p.nifti,'/mprage.nii.gz'],p.anat)
        fprintf('Done...\n')
        mri_list{i_mri_list} = fullfile(p.anat,'mprage.nii');
        i_mri_list = i_mri_list + 1;
    else
        fprintf('Subjects: %s, Cannot find anatomical MRI in subjects anatomical or nifti directories...\n', subjects{i});
        fprintf('Continuing with next subject...\n');
    end
end

%---Generate figure. 
fprintf('=============Generating QA figure=============\n');
V = spm_vol(mri_list);
%--Spawn and hide the spm graphics window
F = spm_figure('Create','Graphics','SPM8w MRI QA','off');
currentplot = 1;
for i_files = 1:length(mri_list)
    fprintf('.');
    mri_data = spm_read_vols(V{i_files});
    mritop   = max(max(max(mri_data)))/2; %Trial and error says this is a nice threhsold
    %-Slices
    slice{1} = squeeze(mri_data(128,:,:)); stitles{1} = 'Sagital'; sthresh{1} = [10,mritop];
    slice{2} = squeeze(mri_data(:,128,:)); stitles{2} = 'Coronal'; sthresh{2} = [10,mritop];
    slice{3} = squeeze(mri_data(:,:,80));  stitles{3} = 'Axial Bottom'; sthresh{3} = [10,mritop];
    slice{4} = squeeze(mri_data(:,:,110)); stitles{4} = 'Axial Top'; sthresh{4} = [10,mritop];
    %-Plot defaults
    s_xpos = [0.05,0.26,0.54,0.75];
    s_tpos = [27,27,27,27];
    colormap gray
    %-Set Currentplot defaults (1 = TOP, 4 = BOTTOM)
    switch(currentplot)
        case 1
             t_ypos = 0.75; %title y-position
             splot = 1; %subplot
             s_ypos = [0.75,0.75,0.75,0.75]; 
        case 2
             t_ypos = 0.51; %title y-position
             splot = 5; %subplot
             s_ypos = [0.51,0.51,0.51,0.51]; 
        case 3
             t_ypos = 0.27;
             splot = 9;
             s_ypos = [0.27,0.27,0.27,0.27];    
        case 4
             t_ypos = 0.02;
             splot  = 13;
             s_ypos = [0.02,0.02,0.02,0.02];    
    end
    %-Plots - 8 per halfpage. So subplot is 4,4, and 1 to 16.
    for iplot = 1:4
        subplot(4,4,splot)
            set(gca,'position',[s_xpos(iplot),s_ypos(iplot),0.2,0.2])
            imagesc(flipud(slice{iplot}'),sthresh{iplot})
            axis equal 
            axis off
            title(stitles{iplot},'fontweight','bold','position',[s_tpos(iplot),0.5])
            splot = splot + 1;
    end
    %-Titles
    %get subname from mri_list (and not subjects since subjects discards
    %people without MRIs so the match from subjects to mri_list isn't
    %exact.
    subname = fileparts(mri_list{i_files});
    subname = subname(length(p.subdir)+1:end);
    %cleanup and be careful of filesep
    if strcmp(subname(1),filesep)
        subname = subname(2:end);
    end
    subname = strtok(subname,filesep);  
    titlestr = sprintf('MRI Subject: %s', subname);
    titleax = axes('Position',[0.1 t_ypos 0.8 0.2],'Parent',F,'Visible','off');
    set(get(titleax,'Title'),'String',titlestr,'FontSize',16,'FontWeight','Bold','Visible','on');
    %Nex tplot
    if currentplot == 4
        eval(defaults.printstr);  
        spm_figure('clear',F);  
        currentplot = 1;
        fprintf('Page complete...\n');
    else
        currentplot = currentplot + 1;
    end
    %-Print Check
    if i_files == length(mri_list)
        eval(defaults.printstr);
        fprintf('\n');
        break;
    end
end     
%---Close hidden figure
spm_figure('Close',F); 

%---Convert multipage ps file to pdf 
%Switching from unix only ps2pdf to matlab ps2pdf.m which
%should be more cross-platform. -DDW Nov/11
ps2pdf('psfile',['tmp_mriqa.ps'],'pdffile',['MRI_QA.pdf'],...
       'gspapersize','a4','deletepsfile',1);
%--delete flak.
if exist(fullfile(p.root,'tmp_mriqa.ps'),'file')
    delete(fullfile(p.root,'tmp_mriqa.ps'));
end

%---Calculate Time and Save Log
qa_stop = datestr(now); 
p.time_elapsed = etime(datevec(qa_stop),datevec(qa_start)); %time elapsed in seconds
[hours, minutes, seconds] = spm8w_timecalc(p.time_elapsed);
fprintf('==============================================================\n');
fprintf('MRI QA finished at %s\nand took %d hours, %d minutes and %d seconds...\n',qa_stop,hours,minutes,seconds);
fprintf('A pdf file of the QA figures was saved to MRIQA.pdf...\n');
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
