function spm8w
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: August, 2010
% ==============================================================================
% spm8w 
%
% spm8w will initialize the SPM8w "GUI" and is the primary means of selecting
% different preprocessing and analysis steps. spm8w should be run from a study's 
% top-level directory and for most steps it will require a parameters file 
% (i.e., P, ROI, VOI). 
%
% An example of how a study directory should be organized and of various
% parameter files is available in:
% EXAMPLE_DATA/2013_H8TJAZZ_SPM8w_r5236
% 
% OPTIONS
% ---------------------------
%   *Preprocess:  Preprocess bold fMRI data according to 'P' file. 
%
%   *Compute GLMs: Comput GLMs using specifications in 'P' file. 
%
%   *Compute Contrasts: Compute linear contrasts of parameter estimates 
%                       using contrasts in 'CON' file. 
%                       NB: If your 'P' file already specifies the 'CON' file 
%                       name for that analysis then spm8w will not ask you 
%                       for the con file. 
%
%                       Once contrasts are complete spm8w will copy the 
%                       con_* files to the RFX folder (named based 
%                       on the con names).
%
%   *Random Effects: Computes a simple one-sample t-test random effects 
%                    analysis on con_* files from 1st level analyses 
%                    and stored in the RFX directory.
% 
%                    In its current form, RFX analysis proceeds without 
%                    global calculation, masking or grand mean scaling. 
%
%   *VOI Extraction: spm8w will ask you for a VOI_study.m (see VOI_H8TJAZZ.m)
%                    file with VOI extraction specifications. This file has 
%                    subject names, coordinates and sphere size for VOI 
%                    timeseries extraction. Importantly this can only be run 
%                    on GLM estimates and contrast computed data. Moreover, 
%                    this requires that you estimated your data with spm8w_3684 
%                    (april 2010) or later as earlier versions had improper 
%                    "effects of interest" contrasts which are required for 
%                    adjusting the timeseries for nuissance variables.
%
%   *PPI Reg Maker:  spm8w will ask you for a VOI_study.m file which
%                    contains specifications for PPI regressor creation.
%                    The regressors will be saved to a user specified 
%                    directory (usually the onsets directory) as 
%                    text files (Y,P and PPI) which can be used in a 
%                    regular GLM as user regressors. 
%                              
%   *PPI Plotter: spm8w will ask you for a VOI_study.m file which
%                 contains the specifics for plotting PPI graphs for each
%                 subject. PPI graphs can be saved in a number of file
%                 formats, in addition the generated regressors used for
%                 plotting are saved as simple txt files in the VOI
%                 directory which can be used for graphing in SPSS or
%                 excel.
%
%   *ROI Analysis: spm8w will ask you for a ROI_study.m 
%                  (see ROI_H8TJAZZ.m) file with ROI analysis 
%                  specifications. This file has coordinates and/or 
%                  image names for ROI analysis. spm8w will then ouput
%                  a tab delimited file with extracted parameter estimates 
%                  and a file of basic statistics as setup in ROI_study.m
%
%                  Currently ROI Analysis can do basic descriptions, 
%                  correlations and one and two sample t-tests. See 
%                  spm8w_roitool.m for more help. 
%
%   *Conjunction Analysis: Simple script to run a conjunction analysis
%                          (logical 'AND' conjunction analysis, see
%                          Tom Nichols paper on conjunctions). This
%                          script is originally from Joe Moran and
%                          calculates the minimum t-value from a number
%                          of spm_t files. Modified to play nice with
%                          spm8w.
%
%   *Monte Carlo Simulations: spm8w will ask you for a P_study.m 
%                             which has specifications for running the 
%                             spm8w wrapper for AFNI's AlphaSim. Masks 
%                             other than bigmask can be specified but 
%                             must be in proper image space and in nifti
%                             format. spm8w no longer supports masks in
%                             afni BRIK format.
% 
%                             Output will go to command line and append to a 
%                             text file in your study root directory.
%
%   *Design Search: spm8w will ask you for a D_study.m file which has
%                   has specifications for the study design and the number
%                   of iterations to perform in order to determine the most
%                   optimal randomization scheme. Output of Design Search
%                   results will go to commmand window. PDFs of estimated
%                   models will be saved in the DESIGN/NAME/MODELS 
%                   directory. Design matrices will be saved in the 
%                   DESIGN/NAME/MATRICES directory.
% 
%   *Design Check: spm8w will ask for a design matrix file and will estimate 
%                  a GLM and calculate design efficiency based on
%                  colliniarity between conditions. Output of Design Search 
%                  results will go to commmand window. No PDFs will be saved.                         
%
%   *Design Build: spm8w will ask you for a D_study.m file which has
%                  has specifications for the study design and for a
%                  design matrix file specifying condition randomization.
%                  Output will go to the directory defined in D_study.m and
%                  will consist of a numbered list of stimuli files ordered
%                  according to the design matrix and split by run.
%
%   *Onset Maker: spm8w will ask you to specify a csv file containing 
%                 onset information (can be stored anywhere but spm8w will 
%                 start by looking in your onsets folder.
%
%                 Output will be a series of txt files with onsets for 
%                 conditions, parametrics and subject specific duration 
%                 files (e.g. RT).
%
%   *PAR/REC Converter: spm8w will ask you for your P_study.m file 
%                       from which it will get the study root dir.
%                       It will pass this on to spm8w_nifticonverter.m
%                       which will ask for which subjects you want to 
%                       convert and will convert the analyse then nifti
%                       and store the results in the appropriate SUBJECTS
%                       dir in your root study dir. 
%                       The idea here is that spm8w_nifticonverter.m
%                       is a wrapper for other converters and can
%                       change over time to accomadate different formats
%                       (DICOM etc.) and different converters (nibabel). 
%
% ==============================================================================
% CHANGE LOG:
% -First version SPM2 scripts - Joe Moran 02/17/06.
% -Added Random Effects - Joe Moran 03/04/06.
% -Added additional batching options and additional ONSET specs - DDW Jan/08
% -spm2_compute and spm2_compute_mixed have been merged (mixed designs are
%  now specified in the params file) - DDW Jan/08
% -Fixed a problem with having more than 100 contrasts (e.g. in mixed
%  designs) -DDW April/08
% -Upgraded the scripts to SPM8! - DDW Dec/09
% -Added ROI analysis tool - DDW Jan/10
% -Put modules into functions to allow for more flexible batching - DDW Jan/10
% -Con file can now be specd in P for the ultra lazy - DDW Feb/10
% -Added AlphaSim wrapper, conjunctions and onset maker to the options - DDW Feb/10
% -Moved redundant code to functions. Re-wrote error handling. Housecleaning - DDW Apr/10
% -Added an editable defaults file to live in users matlab dir for changing
%  stock questions/error messages and file filters - DDW Apr/10
% -RFX analysis will now recursively search through top level RFX directory
% (only one level) - DDW Apr/10
% -Added spm8w_voitool (VOI/PPI) support to the menu - DDW May/10
% -Added support for PPI plotting and DCM (just the button, in prep for
% future inclusion). - DDW June/10
% -Added button for par/rec conversion straight to nifti with proper
% folder structure. - DDW June/10
% -Changed the spm_defaults to spm('Defaults','fmri') - DDW Aug/10
% -Added buttons for Design Search, Check and Build (i.e.
% spm8w_designtool.m) - DDW Mar/11
% -Made windows compatible (using spm8w_osbabel.m and putting paths in
% quotes, which Unix doesn't mind but is necessary for windows paths with
% spaces - DDW Mar/12
% -Added buttons for Resting state preprocessing, Segmentation, Template 
% and mprage and epi normalization - DDW Sep/12
% -Switched to using spm8w_getp for Params. - DDW Sep/12
% ------------------------------------------------------------------------------
% NOTES: 
% There is a different number of planes at GLM (spm2=54, spm8=53).
% but I can find no indication why or that it matters.
% File checks: 
% bold1 = 80x80x30(spm2 and spm8)
% SWUABOLD = 53x65x54 (spm2) 53x65x53 (spm8)
% T2 template = 91x109x91 (spm2 and spm8)
% check bigmask orienation and number of planes -OK!
% 02jun06dr: Models (including motion regressors) are identical
%            Stat maps are identical, occasionnaly small
%            differences, often favoring spm8.
% If/When we move to matlab R2008 use onCleanup for more graceful exiting
% on ctr-c.
% =======1=========2=========3=========4=========5=========6=========7=========8

%---SETUP
%--Clean any previous crashes
close all; clear all; 
%--Set the globals
global cwd bbc
%--Set current working dir
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

%---INIT GUI
initgui()

%---FUNCTIONS
%--Function to start the gui
function initgui()
    close all
    %Generate figure
    spm8wfig = figure('Name','SPM8w_r5236', ...
                      'MenuBar','none', ...
                      'Toolbar','none', ...
                      'NumberTitle','off', ...
                      'Resize', 'off', ...
                      'Position',[50,400,552,616]); %X,Y then Width, Height
    %Change background
    set(spm8wfig, 'Color', ([202,198,191] ./255));
    %Set image
    imgfile    = which('spm8w_title.jpg');
    titleimage = axes('Units', 'pixels','position', [2,536,550,107]); %X,Y then Width, Height
    titleimg   = imread(imgfile);
    image(titleimg); 
    axis off 
    axis image
    %Set Panels
    ypos = 500; %Sets starting Y position of panels. All other position is relative
    p0 = uipanel('BackgroundColor',([47,48,48] ./255),'Units','pixels','Position',[52, 1,450,568]);
    p1 = uipanel('BackgroundColor',([82,30,27] ./255),'Units','pixels','Position',[60, ypos,435,62]);
    p2 = uipanel('BackgroundColor',([92,34,26] ./255),'Units','pixels','Position',[60, (ypos - 66),435,62]);
    p3 = uipanel('BackgroundColor',([102,38,25] ./255),'Units','pixels','Position',[60, (ypos - 104),435,34]);
    p4 = uipanel('BackgroundColor',([112,42,24] ./255),'Units','pixels','Position',[60, (ypos - 170),435,62]);
    p5 = uipanel('BackgroundColor',([122,46,23] ./255),'Units','pixels','Position',[60, (ypos - 208),435,34]);
    p6 = uipanel('BackgroundColor',([122,46,23] ./255),'Units','pixels','Position',[60, (ypos - 274),435,62]);
    p7 = uipanel('BackgroundColor',([112,42,24] ./255),'Units','pixels','Position',[60, (ypos - 312),435,34]);   
    p8 = uipanel('BackgroundColor',([132,50,22] ./255),'Units','pixels','Position',[60, (ypos - 350),435,34]);
    p9 = uipanel('BackgroundColor',([142,54,21] ./255),'Units','pixels','Position',[60, (ypos - 388),435,34]);
    p10 = uipanel('BackgroundColor',([152,58,20] ./255),'Units','pixels','Position',[60, (ypos - 454),435,62]);
    p11 = uipanel('BackgroundColor',([162,64,19] ./255),'Units','pixels','Position',[60, (ypos - 492),435,34]);   
    
    process={'Preprocess',...                         %1
             'Preprocess/GLM',...                     %2
             'Preprocess/GLM/Contrasts',...           %3
             'Estimate GLM',...                       %4
             'Estimate GLM/Contrasts',...             %5
             'Compute Contrasts',...                  %6
             'Random Effects',...                     %7
             'VOI Extraction',...                     %8
             'PPI Reg Maker',...                      %9
             'PPI Plotter',...                        %10
             'ROI Analysis',...                       %11
             'Conjunction Analysis',...               %12
             'Segment Anatomy',...                    %13
             'Template & Flow Fields',...             %14
             'Normalize Anatomy',...                  %15
             'Normalize Functional',...               %16
             'REST Preprocess',...                    %17
             'ART Outlier Detection',...              %18
             'Monte Carlo Simulations (AlphaSim)',... %19
             'Design Search',...                      %20
             'Design Check',...                       %21
             'Design Build',...                       %22
             'Onset Maker',...                        %23
             'PAR/REC Converter'};                    %24
     %Set Buttons
     %Panel 01-Preprocessing
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{1},'Units', 'pixels','Position',[44,31,160,26],'Parent',p1,'Callback',{@choice, 1});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{2},'Units', 'pixels','Position',[229,31,160,26],'Parent',p1,'Callback',{@choice, 2});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{3},'Units', 'pixels','Position',[84,3,260,26],'Parent',p1,'Callback',{@choice, 3});
     %Panel 02-GLM
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{4},'Units', 'pixels','Position',[44,31,160,26],'Parent',p2, 'Callback',{@choice, 4});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{5},'Units', 'pixels','Position',[229,31,160,26],'Parent',p2, 'Callback',{@choice, 5});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{6},'Units', 'pixels','Position',[84,3,260,26],'Parent',p2, 'Callback',{@choice, 6});
     %Panel 03-RFX
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{7},'Units', 'pixels','Position',[84,3,260,26],'Parent',p3,'Callback',{@choice, 7});
     %Panel 04-PPI
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{8},'Units', 'pixels','Position',[44,31,160,26],'Parent',p4, 'Callback',{@choice, 8});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{9},'Units', 'pixels','Position',[229,31,160,26],'Parent',p4, 'Callback',{@choice, 9});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{10},'Units', 'pixels','Position',[84,3,260,26],'Parent',p4, 'Callback',{@choice, 10});      
     %Panel 05-ROI and Conjunction
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{11},'Units', 'pixels','Position',[44,3,160,26],'Parent',p5, 'Callback',{@choice, 11});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{12},'Units', 'pixels','Position',[229,3,160,26],'Parent',p5,'Callback',{@choice, 12});
     %Panel 06-Segment and DARTEL
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{13},'Units', 'pixels','Position',[44,31,160,26],'Parent',p6, 'Callback',{@choice, 13});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{14},'Units', 'pixels','Position',[229,31,160,26],'Parent',p6, 'Callback',{@choice, 14});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{15},'Units', 'pixels','Position',[44,3,160,26],'Parent',p6, 'Callback',{@choice, 15});      
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{16},'Units', 'pixels','Position',[229,3,160,26],'Parent',p6, 'Callback',{@choice, 16});           
     %Panel 07-REST Preprocessing
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{17},'Units', 'pixels','Position',[84,3,260,26],'Parent',p7,'Callback',{@choice, 17});
     %Panel 08-Artifact Detection
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{18},'Units', 'pixels','Position',[84,3,260,26],'Parent',p8,'Callback',{@choice, 18});
     %Panel 09-AlphaSim
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{19},'Units', 'pixels','Position',[84,3,260,26],'Parent',p9,'Callback',{@choice, 19});
     %Panel 10-Design Search, Check & Build
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{20},'Units', 'pixels','Position',[44,31,160,26],'Parent',p10,'Callback',{@choice, 20});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{21},'Units', 'pixels','Position',[229,31,160,26],'Parent',p10,'Callback',{@choice, 21});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{22},'Units', 'pixels','Position',[84,3,260,26],'Parent',p10,'Callback',{@choice, 22});
     %Panel 11-Utilities
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{23},'Units', 'pixels','Position',[44,3,160,26],'Parent',p11, 'Callback',{@choice, 23});
     uicontrol('Style','PushButton','HorizontalAlignment','left','String',process{24},'Units', 'pixels','Position',[229,3,160,26],'Parent',p11, 'Callback',{@choice, 24});
return

%--The Choice function starts the chosen analysis
function choice (h, eventdata, choice)
    %-Load neccessary files according to chosen operation
    close(gcf);
    global bbc cwd

    %-Load subjects for choices 1,2,3,4,5,6,8,9,10,13,14,15,16,18
    if choice == 1 || choice == 2 || choice == 3 || choice == 4 || choice == 5 || choice == 6 || choice == 8 || choice == 9 || choice == 10 || choice == 13 || choice == 14 || choice == 15 || choice == 16 || choice == 17 || choice == 18
        subjects = spm8w_getsub;
    end

    %-Load Parameters File for choices 1,2,3,4,5,6,13,14,15,16,18,19,24
    if choice == 1 || choice == 2 || choice == 3 || choice == 4 || choice == 5 || choice == 6 || choice == 13 || choice == 14 || choice == 15 || choice == 16 || choice == 18 || choice == 19 || choice == 24
        p = spm8w_getp('P');
    end

    %-Load Contrasts File
    if choice == 3 || choice == 5 || choice == 6
        try
            contrast_file = p.confile;
        catch
            contrast_file = spm_select(1,'^CON_.*\.m',bbc.confile,[],[cwd,'/SCRIPTS']);
            %Check for indecisive user
            if isempty(contrast_file)
                error(bbc.error_confile);
            end
        end
        cd(p.root); 
    end

    %-Load VOI/PPI Parameters File
    if choice == 8 || choice == 9 || choice == 10
        v = spm8w_getp('VOI');
    end

    %-Load ROI Parameters File
    if choice == 11
        r = spm8w_getp('ROI');
    end
    
    %-Load REST Parameters File
    if choice == 17
        r = spm8w_getp('R');
    end  
    
    %-Load Design Parameters File
    if choice == 20 || choice == 22
         d = spm8w_getp('D');
    end
    
    %-Perform chosen operations
    switch choice
        case 1 
            % PREPROCESSING     
            vindaloops(subjects, p.para_file, [], p.root, bbc.error_prepro, 'spm8w_preprocess', 'Preprocessing');

        case 2    
            % PREPROCESSING & GLM COMPUTATION    
            vindaloops(subjects, p.para_file, [], p.root, bbc.error_prepro, 'spm8w_preprocess', 'Preprocessing');
            fprintf('==============================================================\n');
            vindaloops(subjects, p.para_file, [], p.root, bbc.error_compute, 'spm8w_compute', 'GLM');

        case 3
            % PREPROCESSING & GLM COMPUTATION & CONTRASTS     
            vindaloops(subjects, p.para_file, [], p.root, bbc.error_prepro, 'spm8w_preprocess', 'Preprocessing');
            fprintf('==============================================================\n');
            vindaloops(subjects, p.para_file, [], p.root, bbc.error_compute, 'spm8w_compute', 'GLM');
            fprintf('==============================================================\n');
            vindaloops(subjects, p.para_file, contrast_file, p.root, bbc.error_contrasts, 'spm8w_contrasts', 'Contrasts');

        case 4 
            % GLM COMPUTATION
            vindaloops(subjects, p.para_file, [], p.root, bbc.error_compute, 'spm8w_compute', 'GLM');

        case 5
            % GLM COMPUTATION & CONTRASTS
            vindaloops(subjects, p.para_file, [], p.root, bbc.error_compute, 'spm8w_compute', 'GLM');
            fprintf('==============================================================\n');
            vindaloops(subjects, p.para_file, contrast_file, p.root, bbc.error_contrasts, 'spm8w_contrasts', 'Contrasts');

        case 6 
            % CONTRASTS
            vindaloops(subjects, p.para_file, contrast_file, p.root, bbc.error_contrasts, 'spm8w_contrasts', 'Contrasts');

        case 7 
            % RFX ANALYSES
            spm8w_rfx();

        case 8
            % VOI EXTRACTION
            try
                spm8w_voitool('VOI', subjects, v.para_file);
            catch
                error_reporter(cwd, bbc.error_voitool, 1);
                error('Exiting...');
            end 
            cd(v.root);
        
        case 9
            % PPI MAKER
            try
                spm8w_voitool('PPI', subjects, v.para_file);
            catch
                error_reporter(cwd, bbc.error_voitool, 1);
                error('Exiting...');
            end 
            cd(v.root);
            
        case 10
            % PPI PLOTTER
            try
                spm8w_voitool('PPI_PLOT', subjects,  v.para_file);
            catch
                error_reporter(cwd, bbc.error_voitool, 1);
                error('Exiting...');
            end 
            cd(v.root);
         
        case 11
            % ROI ANALYSIS
            try
                spm8w_roitool(r.para_file);
            catch
                error_reporter(cwd, bbc.error_roitool, 1);
                error('Exiting...');
            end 
            cd(r.root);

        case 12
            % CONJUNCTION ANALYSIS
            try
                spm8w_conjunction(); 
            catch
                error_reporter(cwd, bbc.error_conj, 1);
                error('Exiting...');
            end 
            cd(cwd);
            
        case 13 
            % SEGMENTATION   
            vindaloops(subjects, p.para_file, [], p.root, bbc.error_prepro, 'spm8w_seg8', 'Segmentation');

        case 14 
            % DARTEL TEMPLATE    
            try
                spm8w_dartel('template', subjects, p.para_file);
            catch
                error_reporter(cwd, bbc.error_voitool, 1);
                error('Exiting...');
            end 
            cd(p.root)
            
        case 15 
            % NORMALIZE ANATOMY     
            vindaloops(subjects, p.para_file, 'mprage2mni', p.root, bbc.error_prepro, 'spm8w_dartel', 'DARTEL Normalizing');            
%             Originally we gave the full list, now we let vindaloops split
%             it up. Same thing. Maybe switch back some day.
%             try
%                 spm8w_dartel('mprage2mni', subjects, p.para_file);
%             catch
%                 error_reporter(cwd, bbc.error_voitool, 1);
%                 error('Exiting...');
%             end 
%             cd(p.root)
%         
        case 16
            % NORMALIZE EPI   
            vindaloops(subjects, p.para_file, 'epi2mni', p.root, bbc.error_prepro, 'spm8w_dartel', 'DARTEL Normalizing');            
      
        case 17 
            % REST PREPROCESSING   
            vindaloops(subjects, r.para_file, [], p.root, bbc.error_prepro, 'spm8w_preprocess', 'REST Preprocessing');

        case 18
            % ART Outlier Detection
            vindaloops(subjects, p.para_file, [], p.root, bbc.error_art, 'spm8w_art', 'ART Outlier Detection');

        case 19
            % ALPHASIM
            try
                spm8w_alphasim(p.para_file);
            catch
                error_reporter(cwd, bbc.error_alphasim, 1);
                error('Exiting...');
            end
            cd(p.root);

         case 20
            % DESIGN SEARCH
            try
                spm8w_designtool('search', design_para_file);
            catch
                error_reporter(cwd, bbc.error_designtool, 1);
                error('Exiting...');
            end
            cd(d.root);
        
        
         case 21
            % DESIGN CHECK
            try
                spm8w_designtool('check');
            catch
                error_reporter(cwd, bbc.error_designtool, 1);
                error('Exiting...');
            end
            cd(cwd);
            
        case 22
            % DESIGN BUILD
            try
                spm8w_designtool('build', design_para_file);
            catch
                error_reporter(cwd, bbc.error_designtool, 1);
                error('Exiting...');
            end
            cd(d.root);
        
        case 23
            % ONSETMAKER
            try
                spm8w_onsmaker();
            catch
                error_reporter(cwd, bbc.error_onsetmaker, 1);
                error('Exiting...');
            end 
            cd(cwd);
            
        case 24
            % NIFTI CONVERTER
            try
                spm8w_nifticonverter(p.para_file);
            catch
                error_reporter(cwd, bbc.error_nifticonv, 1);
                error('Exiting...');
            end 
            cd(cwd);
    end % case statement
    %-Clear all variables to reset subject queue.    
    clear all;
return

%--Loops through subjects running requested command... all hail vindaloops. 
function vindaloops(subjects, para_file, add_arg, rootdir, bbc_error, analysis_func, analysis_name)
    fprintf('Beginning %s on %d subjects\n',analysis_name, length(subjects));
    for z=1:length(subjects);
        x=char(subjects{z});
        try
            if strcmp(analysis_func, 'contrasts')
                evalthis = sprintf('%s (x,para_file,add_arg);',analysis_func);
            elseif strcmp(analysis_func, 'spm8w_dartel')
                evalthis = sprintf('%s (add_arg,x,para_file);',analysis_func);
            else
                evalthis = sprintf('%s (x,para_file);',analysis_func);
            end
            eval(evalthis);
            fprintf('Completed %s on %s (%d of %d subjects)\n',analysis_name, x, z, length(subjects));  
        catch
            error_reporter(rootdir, bbc_error, 2);
            error('Exiting...'); 
        end     
    end
    clear z x
    close all
    cd(rootdir);
    fprintf('%s calculation on %d subjects finished\n',analysis_name, length(subjects));
return
  
%--Report last error message and dump user back to studyroot
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