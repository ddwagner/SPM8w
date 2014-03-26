function spm8w_compute(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: March, 2006
% ==============================================================================
% spm8w_compute(sub_id, parameter_file)
%
% spm8w_compute is an SPM8 adaptation of our spm2_compute scripts for 
% specifying and estimating models.
%
% First argument is subject id, second can be full path to the params file, 
% or else a window will prompt you for it.
% ==============================================================================
% CHANGE LOG:
% -Made some changes to allow for more advanced design specification - DDW Jun/06
% -Added edits by GW for making mixed designs. - DDW Jan/08
% -Added option to have more than 1 parametric/condition, up to 3. 
%  These should be orthogonal (e.g. factors from a PCA). - DDW Dec/09
% -Dragged spm2_compute into the spm8 world. DDW - Dec/09
% -Set to manually specify SPM.xCon since SPM8 no longer specifies this
%  during estimation. - DDW Dec/09
% -Added tweaks so that duration and para can be subject specific but other
% regressors are global. - DDW Jan/10
% -Added option to demean the regressors of interest to make collinearity
% checks easier to interpret. - DDW Feb/10
% -p.cond renamed to p.event_cond to match with more flexible design
% specifications (e.g. regular block designs now possible) in the P file 
% and in spm8w_mk_regressors.m - DDW Feb/10
% -Fixed effects of interest so that it's only user specified regressors 
% (i.e. conditions but not nuissance). This probably won't be correct if the
% user uses the new p.reg_cond and expects those to be nuissance. - DDW Apr/10
% -Housecleaning: re-wrote make_regressors and folded it in as make_nuissance.
% Block and state-item blocks are now inserted at the design spec and not as a 
% hack of make_regressors.
% -Added support for custom regressors. Decided to constrain this to txt
% files only rather than have a bunch of hacks to detect PPI.mat or VOI.mat
% or custom txt files. Use spm8w_voitool to make the appropriate txt
% files. - DDW May/10
% -Changed the spm_defaults to spm('Defaults','fmri') - DDW Aug/10
% -Added ability to model runs seperately (add a p.include_run = #run) to P
% file. - DDW Aug/10
% -TODO: Add ability to compute multiple models based on different onsets
% glm_dir p.event_cond/block_cond and reg_cond must be cell and contain 
% multiple cells. Additionally all the prefix options can contain 
% multiple cells. 
% -PDFs of models are now based on model name (rather than just printing
% them to a single file). -DDW Nov/11
% -Changed behavior to delete prior GLMs if user attempts to overwrite. GLM
% calcs are short and most of the time this is the desired behavior. Not
% doing so will kill batch. -DDW Nov/11
% -Added preliminary code for orthogonality checking which relies on some
% tricky color coding to get positive nad negative correlations to show up
% different. -DDW March/11
% -Mask is now explicitly set in p file (p.mask) - DDW March/12
% -Fixed bug where files don't get spec if not standard prepro -DDW March/12
% -Windows compatible -DDW March/12
% -Now supports addition of outliers(again, only took 2 years!) -DDW March/12
% -Supports inclusion of polynomials as nuissance regressors 
% (p.polyorder) -DDW April/12
% -Fixed bug where design correlations crash on regressor only designs 
%  DDW April/12
% -Added check to ensure onsets file exists -DDW April/12
% -Tweaked outliers to only run if outlier txt files exist & p.outliers=1 
% DDW April/12
% -Fixed bug with outliers code (oops, wrong exist var) -DDW June/12
% -Added DARTEL & VBM8 support -DDW Sept/12
% -Added support for different number of TRs per run. 
% -Added p.duration for cases where evetns/blocks are constant dur. -DDW
% May/13
% -Chnaged the way p.include_run behaves (you know keep p.nses to the full
% runs, see manual). Sep/13
% =======1=========2=========3=========4=========5=========6=========7=========8

%--Input checks
%--Set some vars, clear SPM and make P global
global p; 
GLM_start = datestr(now);
clear SPM;
switch (nargin)
  case 0 
    error('Please specify a subject id, e.g. spm8w_compute(''01jan09ab'')');
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

%---Make p.orth_fix global (for spm_fMRI_design.m)
global disable_orth
disable_orth = p.disable_orth;

%---Adjust p.nses, p.nTR and p.include_run for the included sessions
%---This should go in the private part of the P file but I don't want to
%---make people change until next full release
if ~strcmp(p.include_run,'all')
    p.glm_nTR = p.nTR(p.include_run);
    p.nses = length(p.include_run);
else
    p.glm_nTR = p.nTR;
    p.include_run = 1:p.nses;
end

%---Set Defaults
%%%---------------------------------------------------
%%% User-defined parameters for this analysis
%%%---------------------------------------------------
SPM.nscan          = sum(p.glm_nTR);        % number of scans for each of nsess sessions
SPM.xY.RT          = p.TR;                  % experiment TR in seconds
SPM.xGX.iGXcalc    = 'None';                % global normalization: OPTIONS:'Scaling'|'None'
SPM.xX.K.HParam    = p.hpf;                 % high-pass filter cutoff (secs) [Inf = no filtering] -------DDW CHANGE
SPM.xVi.form       = p.autocorr;            % intrinsic autocorrelations: OPTIONS: 'none'|'AR(1) + w' ------DDW CHANGE
%%%---------------------------------------------------
%%% basis functions and timing parameters
%%%---------------------------------------------------

SPM.xBF.name       = p.hrf;         % OPTIONS:'hrf'
                                    %         'hrf (with time derivative)'
                                    %         'hrf (with time and dispersion derivatives)'
                                    %         'Fourier set'
                                    %         'Fourier set (Hanning)'
                                    %         'Gamma functions'
                                    %         'Finite Impulse Response'
                          
SPM.xBF.T          = 16;            % number of time bins per scan
SPM.xBF.T0         = 1;             % reference time bin 
SPM.xBF.UNITS      = p.time;        % OPTIONS: 'scans'|'secs' for onsets ----------------DDW CHANGE
SPM.xBF.Volterra   = 1;             % OPTIONS: 1|2 = order of convolution; 1 = no Volterra

%--Check for presence of FIR variables
if isfield(p,'hrfwindow') && isfield(p,'hrfbasis')
    SPM.xBF.length     = p.hrfwindow;            % Length of HRF window
    SPM.xBF.order      = p.hrfbasis;             % Order of Basis Set
end

%---Setup directories
if ~exist(p.res, 'dir')
    mkdir(p.res)
    cd(p.res);
else
    cd(p.res) ;
end
if ~exist(p.glm, 'dir')
    mkdir(p.glm)
    cd(p.glm);
else
    fprintf('Prior GLM directory exists at: %s\nDeleting the directory and proceeding...\n', p.glm);
    rmdir(p.glm, 's')
    mkdir(p.glm)
    cd(p.glm);
end

%---Determine boldtok and rptok and add to p structure
p.boldtok = []; p.rptok = [];
if p.despike;   p.boldtok = ['k',p.boldtok]; p.rptok = ['k',p.rptok]; end
if p.slicetime; p.boldtok = ['a',p.boldtok]; p.rptok = ['a',p.rptok]; end
if p.unwarp;    p.boldtok = ['u',p.boldtok]; end
if p.normalize && strcmp(p.normtype,'OLD'); p.boldtok = ['w',p.boldtok]; end
if p.normalize && strcmp(p.normtype,'DARTEL'); p.boldtok = ['d',p.boldtok]; end
if p.normalize && strcmp(p.normtype,'VBM8'); p.boldtok = ['v',p.boldtok]; end
if p.smoothing; p.boldtok = ['s',p.boldtok]; end
p.boldtok = [p.boldtok, p.bold];
p.rptok   = [p.rptok, p.bold];

%---Get GLM name from p.glm for display purposes
[glm_dir, glm_name] = fileparts(p.glm);

%---Setup design matrix 
fprintf('==========Beginning Design and GLM estimation for %s on %s\n', glm_name, p.subj); 

%---Init Trial and User Regressors 
SPM.Sess.C.C    = [];
SPM.Sess.C.name = {};
SPM.Sess.U      = [];   %May 2010-DDW

%---EVENT-RELATED DESIGNS (LOAD ONSETS, PAR, AND DUR FOR EVENT CONDITIONS)
if ~isempty(p.event_cond)
    for x = 1:length(p.event_cond)
        %--PRINT INFO 
        onsets = fullfile(p.onsets,[p.onsets_prefix,p.event_cond{x},p.onsets_ext]);
        if ~exist(onsets,'file')
            error('Error, Cannot find the file: %s... Are you sure it exists?\n',spm_str_manip(onsets,'t'));     
        end
        fprintf('===Loading event onset file: %s...\n',spm_str_manip(onsets,'t'));     
        %--LOAD ONSETS
        onsfile = spm_load(onsets);
        %--Make column vector if row vector
        if size(onsfile,1) < size(onsfile,2); onsfile=onsfile'; end
        %--Durations
        durat = fullfile(p.onsets,[p.duration_prefix,p.event_cond{x},'_dur',p.onsets_ext]);
        %--Check for existence of _dur files.
        if exist(durat, 'file')
            %-PRINT INFO 
            fprintf('\nLoading event duration file: %s... \n',spm_str_manip(durat,'t'));     
            %-LOAD DURATIONS
            durations = spm_load(durat);
            %-Make column vector if row vector
            if size(durations,1) < size(durations,2); durations=durations'; end   
        else
            durations = p.duration;
        end
        %--If dur in seconds, divide by TR.
        if (p.durtime)
            if ~exist(durat,'file')
                fprintf(['Durations (%d) in seconds, dividing by '...
                        'TR of%4.1f\n'],durations,p.TR);
            else
                fprintf(['Durations in seconds, dividing by '...
                        'TR of%4.1f\n'],p.TR);  
            end
            durations = durations/p.TR;
        else           
            if ~exist(durat,'file')
                fprintf('Duration: %dTRs...\n',durations);
            end
        end
        %--Assign names and onsets to SPM structure
        SPM.Sess.U(x).name = {strrep(p.event_cond{x},'_','-')};   % string in cell
        SPM.Sess.U(x).ons  = onsfile;                             % onsets in scans        
        SPM.Sess.U(x).dur = durations;                            % durations in scans
        %--CHECK FOR PARAM FILES    
        param = fullfile(p.onsets,[p.para_prefix,p.event_cond{x},'_par', p.onsets_ext]);
        par1  = fullfile(p.onsets,[p.para_prefix,p.event_cond{x},'_par1',p.onsets_ext]);    
        par2  = fullfile(p.onsets,[p.para_prefix,p.event_cond{x},'_par2',p.onsets_ext]);   
        par3  = fullfile(p.onsets,[p.para_prefix,p.event_cond{x},'_par3',p.onsets_ext]);  
        %--PARAMETRIC MODULATION
        %--Check for existence of _par files. Modified to allow for more than
        %--one parametric (currently a max of 3). -Dec 2009 DDW
        if exist(param, 'file')
            fprintf('===Loading event par file: %s...\n',spm_str_manip(param,'t'));     
            parfile = spm_load(param);
            %-Make column vector if row vector
            if size(parfile,1)<size(parfile,2); parfile=parfile'; end
            SPM.Sess.U(x).P.name = 'other';         % 'none' | 'time' | 'other'
            SPM.Sess.U(x).P.P    = parfile;         % vector same length as onsets
            SPM.Sess.U(x).P.h    = 1;               % order of polynomial expansion
        elseif exist(par1,'file')
            for i = 1:3
                if exist(eval(['par',num2str(i)]),'file')
                    fprintf('===Loading event par file: %s...\n',spm_str_manip(eval(['par',num2str(i)]),'t'));
                    parfile = spm_load(eval(['par',num2str(i)]));  % Get the first par file. Dec2009 - DDW
                    %-Make column vector if row vector
                    if size(parfile,1) < size(parfile,2); parfile=parfile'; end
                    %-Modified the SPM.Sess.U(x).P(i) to make it a structure to
                    %-accomidate multiple regressors. -Dec 2009 DDW
                    SPM.Sess.U(x).P(i).name = ['other',num2str(i)];  % 'none' | 'time' | 'other'  
                    SPM.Sess.U(x).P(i).P    = parfile;               % vector same length as onsets
                    SPM.Sess.U(x).P(i).h    = 1;                     % order of polynomial expansion
                end 
            end
        else
            SPM.Sess.U(x).P.name        = 'none';                    % none for events without mod   
        end
    end %for x=1:length(p_event_cond)
end % condition x loop    
      
%---BLOCK DESIGNS
%---CHECK BLOCK DESIGN TYPE (BLOCK OR STATE ITEM)
if ~isempty(p.block_cond) %If there are no events, then pure block design
    for x=1:length(p.block_cond)
        %--LOAD ONSETS AND DURATIONS
        onsets = fullfile(p.onsets,[p.onsets_prefix,p.block_cond{x},p.onsets_ext]);
        if ~exist(onsets,'file')
            error('Error, Cannot find the file: %s... Are you sure it exists?\n',spm_str_manip(onsets,'t'));     
        end
        fprintf('===Loading block onset file: %s...',spm_str_manip(onsets,'t'));          
        onsfile = spm_load(onsets);
        durat = fullfile(p.onsets,[p.duration_prefix,p.block_cond{x},'_dur',p.onsets_ext]);
        %--Check for existence of _dur files.
        if exist(durat, 'file')
            %-PRINT INFO 
            fprintf('\nLoading block duration file: %s... \n',spm_str_manip(durat,'t'));     
            %-LOAD DURATIONS
            durations = spm_load(durat);
            %-Make column vector if row vector
            if size(durations,1) < size(durations,2); durations=durations'; end          
        else
            durations = p.duration;
        end
        %--If dur in seconds, divide by TR.
        if (p.durtime)
            if ~exist(durat,'file')
                fprintf(['Durations (%d) in seconds, dividing by '...
                        'TR of%4.1f\n'],durations,p.TR);
            else
                fprintf(['Durations in seconds, dividing by '...
                        'TR of%4.1f\n'],p.TR);  
            end
            durations = durations/p.TR;
        else           
            if ~exist(durat,'file')
                fprintf('Duration: %dTRs...\n',durations);
            end
        end
        %--Make column vector if row vector
        if size(onsfile,1) < size(onsfile,2); onsfile=onsfile'; end
        %--Determine design type
        if ~isempty(p.block_cond) && isempty(p.event_cond)      %BLOCK DESIGN
            %-Assign names and onsets to SPM structure
            SPM.Sess.U(x).name   = {strrep(p.block_cond{x},'_','-')}; % string in cell
            SPM.Sess.U(x).ons    = onsfile;               % onsets in scans
            SPM.Sess.U(x).P.name = 'none';                % no paras for block designs  
            SPM.Sess.U(x).dur    = durations; 
        elseif ~isempty(p.block_cond) && ~isempty(p.event_cond) %MIXED DESIGN
            fprintf('Events and blocks are specified, treating %s as user regressor (no HRF convolution)...\n',spm_str_manip(onsets,'t'));     
            onsfile    = onsfile + 1; %%%Fix for block indexing problem - march 2010 - ddw
            %block(x).name = {p.block_cond{x}};  % string in cell
            %-Specify Condition Block Variables
            state_data = zeros((sum(p.glm_nTR)),1);
            %-Add Condition Block Onsets to Condition Block Variable
            for i = 1:length(onsfile);
                state_data(onsfile(i):(onsfile(i)+durations(i)),1) = 1;
            end
            %-Assign names and onsets to SPM structure
            SPM.Sess.C.C(:,end+1)  = state_data;
            SPM.Sess.C.name{end+1} = strrep(p.block_cond{x},'_','-');
            %-Removed convolving blocks with HRF based on new info from GW&LHS -DDW
            %-Truncate State Variable so that its the Same Length as the Scan Length
            %-block(x).state=state_tmp(1:(p.nses*p.nTR)); %not necessary
            %-anymore
        end
    end %for x=1:length(p_block_cond)
end % ~isempty

%---USER REGRESSOR DESIGNS
%---Making this .txt only, use spm8w_voitool in
%---conjunction with a VOI_XXXX.m file in your 
%---SCRIPTS directory, this way there's nothing special
%---about PPI or Timecourse or other regression analysis
%---It's just an onset file with row = # of TRs.
if ~isempty(p.reg_cond) && isempty(p.event_cond) && isempty(p.block_cond) %then pure reg design
    for x=1:length(p.reg_cond)
        %--PRINT INFO
        onsets = fullfile(p.onsets,[p.onsets_prefix,p.reg_cond{x},p.onsets_ext]);
        if ~exist(onsets,'file')
            error('Error, Cannot find the file: %s... Are you sure it exists?\n',spm_str_manip(onsets,'t'));     
        end
        fprintf('====Loading regressor file: %s ...\n',spm_str_manip(onsets,'t'));     
        %--LOAD REGRESSOR
        onsfile = spm_load(onsets);
        %--Make column vector if row vector
        if size(onsfile,1) < size(onsfile,2); onsfile=onsfile'; end
        %--Check size of regressor and drop error if mismatch
        if length(onsfile) ~= sum(p.glm_nTR)
           error(['Length of regressor %s is not the same as size of design\n' ...
                  'Size of regressor = %d size of design = %d\n'],spm_str_manip(onsets,'t'),length(onsfile),sum(p.glm_nTR));
        end
        %--Assign names and regressors to SPM structure  
        SPM.Sess.C.C(:,end+1)  = onsfile;
        SPM.Sess.C.name{end+1} = strrep(p.reg_cond{x},'_','-');   
    end %for x=1:length(p_reg_cond)
elseif ~isempty(p.reg_cond)
    error(['For now event/block and regressor designs are disabled\n' ...
           'since I can''t think of a use, email me if you need this -ddw']);   
end % condition x loop     
   
%---BUILD AND ADD NUISSANCE REGRESSORS
[regress, regname] = make_nuissance(p);  %now calls make_nuissance - DDW Apr/10
for i=1:size(regress{1},2)
    SPM.Sess.C.C(:,end+1)  = regress{1}(:,i);
    SPM.Sess.C.name{end+1} = regname{i};
end

%---Setup estimation 
%--make sure defaults are present in workspace
spm('Defaults','fmri');

%--Get filenames for functional scans
clear scans;
if strcmp(p.include_run, 'all')  %Convert 'all' to a list of runs
    p.include_run = 1:p.nses;
end
fprintf('Loading images of type: %s\n', [p.boldtok,'.nii']);
for srun_i = 1:p.nses
    fprintf('Loading images for run %d (nTR = %d)...\n',p.include_run(srun_i),p.glm_nTR(srun_i));
    scans{srun_i} = spm_select('FPList',p.func,['^',p.boldtok,num2str(p.include_run(srun_i)),'\.nii']);
end

%--specify data: matrix of filenames and TR
SPM.xY.P = char(scans);  %SPM expects a char array not a cell array
%--Configure and print design matrix.
%--This is where the SPM.xX structure gets filled in that we need for
%--demeaning 2010 DDW
SPM   = spm_fmri_spm_ui(SPM);
G     = spm_figure('FindWin','Interactive'); close(G);
F     = spm_figure('FindWin','Graphics'); set(F,'visible','off');

%%% Generate orothogonality figure -March/12 DDW
FS    = spm('FontSizes');
%%% Freeze Colormaps using freezeColors.m v2.3
freezeColors       
%%% Add correlation matrix of predictors to SPM figure
%Returns the upper triangular correlation matrix. 
if isempty(SPM.Sess.U)    %Check if regressor only design (i.e. PPI)
    cormat_tmp = corrcoef(SPM.xX.X(:,1:length(p.reg_cond)));        %Correlation matrix between user regressors
    cormat_tmp(logical(tril(ones(length(p.reg_cond)),-1))) = 0; %Replace lower diag with zeros (which will be white)
else
    cormat_tmp = corrcoef(SPM.xX.X(:,1:numel(SPM.Sess.U)));        %Correlation matrix between user regressors
    cormat_tmp(logical(tril(ones(numel(SPM.Sess.U)),-1))) = 0; %Replace lower diag with zeros (which will be white)
end
nPar       = length(cormat_tmp);                          %Get number of regressors
PTick      = spm_DesRep('ScanTick',nPar,32);
hDesO      = axes('Position',[.50 .1 .3 .14]);
hDesOIm    = imagesc(cormat_tmp);
caxis([-1 1])                           %set caxis limites to 1 -1
load('cormap.mat');                     %cormat that ships with spm8w for correlations.
colormap(cormap)
freezeColors                            %freeze and reset colors back to gray.
colormap(gray)
%Clean axes
set(hDesO,'Box','off','TickDir','out',...
        'XaxisLocation','top','XTick',PTick,'XTickLabel','',...
        'YaxisLocation','right','YTick',PTick,'YTickLabel','',...
        'YDir','reverse')
tmp = [1,1]'*[[0:nPar]+0.5];
line('Xdata',tmp(1:end-1)','Ydata',tmp(2:end)')
axes('Position',[.81 .1 0.01 .139],'Visible','off',...
     'DefaultTextFontSize',FS(8),'DefaultTextInterpreter','TeX',...
     'YDir','reverse','YLim',[0,nPar]+0.5)
for i=PTick   
    if isempty(SPM.Sess.U)  %fix for PPI/User reg designs -DDW April/12 
        text(0,i,SPM.Sess.C.name(i),'HorizontalAlignment','left')   %add condition names
    else
        text(0,i,SPM.Sess.U(i).name,'HorizontalAlignment','left')   %add condition names    
    end
end
%Set title
str = 'Design Correlations...';   
hAx = axes('Position',[0.38,0.05,0.94,0.22],'Visible','off');
set(hAx,'Units','points');
AxPos = get(hAx,'Position');
set(hAx,'YLim',[0,AxPos(4)])
dy = FS(9); y0 = floor(AxPos(4)) -dy; y = y0;
text(0.3,y,str,...
    'HorizontalAlignment','Center',...
    'FontWeight','Bold','FontSize',FS(11))   

%--Print figure, overwriting previous models of same name. Latest matlab
%--seems to bungle figure printing so this has gotten hackier...
%--The idea is we have to set the figure's paperunits to correspond to
%--the size of the figure in inches, plus some flak to position things 
%--appropriately. This makes a figure that's bigger than 8.5x11 but 
%--c'est la vie! -DDW March/13
set(F,'Units','inches');
set(F,'PaperUnits','inches');
pos = get(F,'Position');
set(F,'PaperSize',[pos(3) pos(4)]);
set(F,'PaperPositionMode','manual');
set(F,'PaperPosition',[0 0 pos(3) pos(4)]);
printstr = ['print -dpdf -zbuffer -noui ../../',p.subj,'_',lower(glm_name),'.pdf'];
eval(printstr);
spm_figure('close',F); drawnow;

%--Insert the bigmask into the design matrix. 
%--This assumes you have only one copy of bigmask.nii in your matlab path
%--If you have more than one copy then what the hell are you doing?
%--Swapped bigmask.nii to p.mask which users can alter in the P file.
mask_name = which(p.mask);
mask_vol  = spm_vol(mask_name);
SPM.xM.VM = mask_vol;
SPM.xM.T  = [];
SPM.xM.TH = ones(size(SPM.xM.TH))*(-Inf);
SPM.xM.I  = 0;
SPM.xM.xs = struct('Masking', sprintf('explicit masking only - using %s',p.mask));  %DDW 03/12
%--Print mask msg
fprintf('\n==========The SPM.mat file has been modfied. Estimating the model...');
fprintf('\nThe model will include all voxels in the brain as defined by the mask at:\n');
fprintf('%s\n', mask_name);
%--Manually demean all regressors in order for SPM8 orthogonality checks to
%--be equivalent to spm99! Everything but constant and conditions are already
%--mean centered, reg_interest is calculated by counting regressor names in
%--SPM.Sess and amount of regressors whose P has a name other than none.
%--Still not sure if valid to do this for FIR models or mixed designs.
%--Feb 2010 DDW
if(p.demean)
    reg_interest = size(SPM.Sess.U,2);
    for i = 1:size(SPM.Sess.U,2)
        for ii = 1:size(SPM.Sess.U(i).P,2)
            if strfind(SPM.Sess.U(i).P(ii).name,'other')
                reg_interest = reg_interest + 1;
            end
        end
    end
    for i = 1:reg_interest
        %SPM.xX.X(:,i)=SPM.xX.X(:,i)-mean(SPM.xX.X(:,i));
        SPM.xX.X(:,i)=spm_detrend(SPM.xX.X(:,i)); %same as above but uses spm function
    end
end
%--Put the mask and the demeaned design regressors into the SPM structure
save SPM SPM;

%%% Removed call to specmask as all we need is this small segment. If ever 
%%% we need to make a new bigmask we'll have to go back and reexamine
%%% specmask.m
%SPM = specmask(SPM, mask_name);
%%% Specify new mask image
%%%
%%%     - This option could be used if you don't like the default mask SPM
%%%       creates - you can supply a mask file here or ask it to create one
%%%       from this subject's anatomy.  See specmask.m for more details.
%%% If you have a mask file already, use this.  mask_name is a full-path
%%% filename to a mask image specifying which voxels to include in estimation.
%%% If you don't have a mask file already, but would like to create one based
%%% one your subject's anatomy to explicitly include the entire brain rather
%%% than simply using SPM's default mask, use this.  anat_name is a full-path
%%% filename to an anatomical image (normalized if you've normalized your 
%%% functionals); template is the appropriate SPM template (for segmentation).
%%%
%%% SPM = specmask(SPM, anat_name, fullfile(spm('Dir'),'templates', 'T1.img'));

%---Estimate Design
SPM = spm_spm(SPM);
F   = spm_figure('FindWin','Graphics');    %Get the graphics window
G   = spm_figure('FindWin','Interactive'); %Get the progress bar window
close(F);
close(G);drawnow;
%---Calculate collinearity and taskXmotion interactions
%---This is for the future. -DDW
%%% Load a fresh copy of the SPM structure.
%SPM = load(SPM)
%%%Convolved regressors are in SPM.xX.X
%%%FIgure out how many are task
%length(p.event_cond);
%corr(1:length(p.event_cond),1:length(p.event_cond)) = zeros;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Specify SPM.xCon variable
%%% SPM8 no longer specifies the SPM.xCon var after 
%%% estimation so we have to do it manually. -Dec 2009 DDW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%% Set contrasts in SPM struture
%%% Create 1st contrast for 'effects of interest
%%% Fixed so that effects of interest includes only 
%%% conditions of interest. This will make it easier to
%%% adjust timecourses with voi button and do PPI among
%%% others. 
Fcname = 'effects of interest';
iX0    = 1:SPM.xX.iB;   %Set iX0 to span all regressors
%--Size of SPM.Sess.U doesn't work because of user regressors and FIR models
%--so have to use the 'r' name trick from the con file. - DDW Apr/10
%--len = sum(cell2mat(strfind(SPM.Sess.C.name,'r')))+1; %%% Count number of nuissance 'r'egs and add +1 for constant
%iX0(1:size(SPM.Sess.U,2)) = [];            %Remove all but nuissance (this may not work for all designs)      
len = sum(strcmp(SPM.Sess.C.name,'r'))+1;   %Count number of nuissance 'r'egs and add +1 for constant
iX0(1:end-len) = [];
xCon = spm_FcUtil('Set',Fcname,'F','iX0',iX0,SPM.xX.xKXs);

%--Append contrasts for fMRI - specified by SPM.Sess(s).Fc(i)
%--This is nice but we delete these at the Contrasts level...
%--This gives us free F tests for all conditions but 
%--Would make the cons for FIR designs overly long so
%--Commenting out since we never used it anyway... -DDW Apr/10
% if isfield(SPM,'Sess')
%     for s = 1:length(SPM.Sess)
% 	for i = 1:length(SPM.Sess(s).Fc)
% 	    iX0           = 1:SPM.xX.iB;
% 	    iX            = SPM.Sess(s).col(SPM.Sess(s).Fc(i).i);
% 	    iX0(iX)       = [];
% 	    Fcname        = sprintf('%s',SPM.Sess(s).Fc(i).name); %Removed the 'sess(1):' prefix put it back if
%                                                               %we batch in non concatenated sessions                                                             
% 	    xCon(end + 1) = spm_FcUtil('Set',Fcname,'F','iX0',iX0,SPM.xX.xKXs);
% 	end
%     end
% end

SPM.xCon = xCon;
save SPM SPM;

%---Touch it and return to root!
eval(spm8w_osbabel(sprintf('!touch "%s"',fullfile(p.subdir,p.subj,['estimated_',p.logsuffix]))));
cd(p.root);

%---Spit Time
GLM_stop                  = datestr(now);
GLM_time_elapsed          = etime(datevec(GLM_stop),datevec(GLM_start)); %time elapsed in seconds
[hours, minutes, seconds] = spm8w_timecalc(GLM_time_elapsed);
fprintf('Subject %s GLM computation finished at %s\n and took %d hours, %d minutes and %d seconds...\n',p.subj,GLM_stop,hours,minutes,seconds);

function[regress, regname] = make_nuissance(p)
% Function to make nuissance regressors 
% 
% This function will build regressors based on outlier scans, linear
% trends, sessions and realignment parameters.
% Basically same as Joe Moran's mk_regressors.m (7/19/04).
% 
% To prepare for future inclusion of user regressor batching (timecourse,
% PPI designs) I've removed my state-item hacks from mk_regressors.m and 
% folded everything in spm8w_compute.m It was getting a little ridiculous
% with all kinds of daft checks to figure out if we were running block,
% mixed or and er study. - DDW Apr/10
%   
% MODIFICATIONS:
% -Added some of JM's updates (e.g. outlier regressors) - DDW - Jan/08
% -Added a few lines to grab rp_bold if slicetime was not performed - DDW Jan/08
% -Dragged spm2_mk_regressors into the spm8 world. - DDW Dec/09
% -Folded in spm2_state_regressors so no longer need to eval external file - DDW Dec/09
% -Length of block regressors is now dicated by a blockname_dur.txt file.
%  Allowing for more flexible sate/item and block designs. - DDW Feb/10
% -Added p.block_convolve so that we can do both block (convolve) and mixed
% designs (do not convolve) - DDW Feb/10
% -Added fix for error in block onsets. Indexing started at 1 when should
% start at 0 like spm. So added +1 to all block onses. - DDW Mar/10
% -Removed mixed design code from mk_regressors. Also retired
% p.block_convolve and just use a check to see if state-item or regular
% block. - DDW Apr/10
% -Cleaned code (append rather than guess size of final regressor matrix)
%  and added back support for outliers. - DDW Mar/12
% -Removed linear trends and converted to polynomials (linear trend = 1st
% order poly; Quadratic = 2nd order; Cubic = 3rd order; Constant = 0).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify default variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_tps  = sum(p.glm_nTR);
regressors = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make outlier regressors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if p.outliers
    count = 1;
    for i = 1:p.nses  %Make a variable per outlier, adjusting for run lenghths
        if exist(fullfile(p.func,['outliers_run',num2str(p.include_run(i)),'.txt']),'file')    
            outliers{i} = load(fullfile(p.func,['outliers_run',num2str(p.include_run(i)),'.txt']));
            fprintf('===Loading outliers file for run: %d...\n',p.include_run(i));     
            for j = 1:length(outliers{i})
                str=['out',num2str(count),' = zeros(',num2str(total_tps),',1);'];            eval(str);
                str=['out',num2str(count),'(',num2str(outliers{i}(j)+(i-1)*p.glm_nTR(i)),') = 1;']; eval(str);
                count = count+1;
            end
        end
    end
    if exist('out1','var')
        out_regs          = zeros(total_tps,length([outliers{:}])); %one col per outlier
        for i = 1:size(out_regs,2);
            str           = ['out',num2str(i)];
            out_regs(:,i) = eval(str);
        end  
        regressors = [regressors, out_regs];     
    end
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
% Build polynomial regressors and session means
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
polys = zeros((total_tps),(p.polyord+1)*p.nses-1); % init polynomial regs
count = 1;
for i = 1:p.nses
    fprintf('===Generating polynomial regressors for run: %d...\n',p.include_run(i));     
    for ii = 1:p.polyord
        polys(i*p.glm_nTR(i)-p.glm_nTR(i)+1:i*p.glm_nTR(i),count) = linspace(-1,1,p.glm_nTR(i)).^ii';
        count = count + 1;
    end
end
%Now do 0th order poly (i.e. session mean)
for i = 1:p.nses-1
    ii = 0;
    polys(i*p.glm_nTR(i)-p.glm_nTR(i)+1:i*p.glm_nTR(i),count) = linspace(-1,1,p.glm_nTR(i)).^ii';
    count = count + 1;
end
regressors = [regressors, polys];      % add linear trends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make realignment regressors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if p.move   
    realign = [];
    for i = 1:p.nses
        fprintf('===Loading realignment parameters file for run: %d...\n',p.include_run(i));     
        realign{i} = load(fullfile(p.func,['rp_',p.rptok,num2str(p.include_run(i)),'.txt'])); 
    end
    realign  = cat(1,realign{:});   
    %%% Insert realignment regressors into regressors variable
    %%% fill last n cols of regressor matrix, leaving first 
    %%% (2*runs)-1 free for trends and means 
    regressors = [regressors, realign];      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make regressor names variable, one name per regressor
% initialise regressors variable with rows = total_TPs,
% columns = # of total regressors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:size(regressors,2);
    reg_names(j,:) = 'r';
end

regress{1}=regressors;                               % convert regressor variables to
regname=cellstr(reg_names);                          % cells, 
return
