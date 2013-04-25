function varargout = spm8w_do_rfx(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: March, 2006 - JM
% ==============================================================================
% function spm8w_do_rfx('rfx_dir',D_structure)
% 
% spm8w_do_rfx will compute a randome ffects analysis on the con  images in 
% argument 1 using the stats defined in the stats_structure variable. If you 
% leave argument 2 blank it will run a one-sample t-test, no grand mean 
% scaling, no masking, no global calculation. If you want to change these 
% settings, alter the stats_structure in spm8w.m If there is a demand, we can 
% move this structure to the Params file at the cost of forcing the user to 
% specify the Params file everytime they run an RFX.
%
% spm_spm_ui for details on how this structure can be built and specified.
%
% The stuff following is taken from spm_spm_ui from SPM8 (see actual
% spm_spm_ui.m file for full details).
% ==============================================================================
% CHANGE LOG:
% -Joe Moran 03/05/06
% -Pushed the design structure into the spm8w file rather than its own
%  seperate .mat file in the path -DDW Dec/09
% -Converted to spm8 -DDW Dec/10
% -Housecleaning - DDW Apr/10
%
% FROM SPM_SPM_UI.M
% Variables saved in the SPM stucture
%
% xY.VY         - nScan x 1 struct array of memory mapped images
%                 (see spm_vol for definition of the map structure)
% xX            - structure describing design matrix
% xX.D          - design definition structure
%                 (See definition in main body of spm_spm_ui.m)
% xX.I          - nScan x 4 matrix of factor level indicators
%                 I(n,i) is the level of factor i corresponding to image n
% xX.sF         - 1x4 cellstr containing the names of the four factors
%                 xX.sF{i} is the name of factor i
% xX.X          - design matrix
% xX.xVi        - correlation constraints for non-spericity correction
% xX.iH         - vector of H partition (condition effects) indices,
%                 identifying columns of X correspoding to H
% xX.iC         - vector of C partition (covariates of interest) indices
% xX.iB         - vector of B partition (block effects) indices
% xX.iG         - vector of G partition (nuisance variables) indices
% xX.name     - p x 1 cellstr of effect names corresponding to columns
%                 of the design matrix
% 
% xC            - structure array of covariate details
% xC(i).rc      - raw (as entered) i-th covariate
% xC(i).rcname  - name of this covariate (string)
% xC(i).c       - covariate as appears in design matrix (after any scaling,
%                 centering of interactions)
% xC(i).cname   - cellstr containing names for effects corresponding to
%                 columns of xC(i).c
% xC(i).iCC     - covariate centering option
% xC(i).iCFI    - covariate by factor interaction option
% xC(i).type    - covariate type: 1=interest, 2=nuisance, 3=global
% xC(i).cols    - columns of design matrix corresponding to xC(i).c
% xC(i).descrip - cellstr containing a description of the covariate
% 
% xGX           - structure describing global options and values
% xGX.iGXcalc   - global calculation option used
% xGX.sGXcalc   - string describing global calculation used
% xGX.rg        - raw globals (before scaling and such like)
% xGX.iGMsca    - grand mean scaling option
% xGX.sGMsca    - string describing grand mean scaling
% xGX.GM        - value for grand mean (/proportional) scaling
% xGX.gSF       - global scaling factor (applied to xGX.rg)
% xGX.iGC       - global covariate centering option
% xGX.sGC       - string describing global covariate centering option
% xGX.gc        - center for global covariate
% xGX.iGloNorm  - Global normalisation option
% xGX.sGloNorm  - string describing global normalisation option
% 
% xM            - structure describing masking options
% xM.T          - Threshold masking value (-Inf=>None,
%                 real=>absolute, complex=>proportional (i.e. times global) )
% xM.TH         - nScan x 1 vector of analysis thresholds, one per image
% xM.I          - Implicit masking (0=>none, 1=>implicit zero/NaN mask)
% xM.VM         - struct array of explicit mask images
%                 (empty if no explicit masks)
% xM.xs         - structure describing masking options
%                 (format is same as for xsDes described below)
% 
% xsDes         - structure of strings describing the design:
%                 Fieldnames are essentially topic strings (use "_"'s for
%                 spaces), and the field values should be strings or cellstr's
%                 of information regarding that topic. spm_DesRep.m
%                 uses this structure to produce a printed description
%                 of the design, displaying the fieldnames (with "_"'s 
%                 converted to spaces) in bold as topics, with
%                 the corresponding text to the right
% 
% SPMid         - String identifying SPM and program versions
%=======================================================================
% - FORMAT specifications for programers
%=======================================================================
%( This is a multi function function, the first argument is an action  )
%( string, specifying the particular action function to take.          )
%
% FORMAT spm_spm_ui('CFG',D)
% Configure design
% D       - design definition structure - see format definition below
% (If D is a struct array, then the user is asked to choose from the
% design definitions in the array. If D is not specified as an
% argument, then user is asked to choose from the standard internal
% definitions)
%
% FORMAT [P,I] = spm_spm_ui('Files&Indices',DsF,Dn,DbaTime)
% PET/SPECT file & factor level input
% DsF     - 1x4 cellstr of factor names (ie D.sF)
% Dn      - 1x4 vector indicating the number of levels (ie D.n)
% DbaTime - ask for F3 images in time order, with F2 levels input by user?
% P       - nScan x 1 cellsrt of image filenames
% I       - nScan x 4 matrix of factor level indices
%
% FORMAT D = spm_spm_ui('DesDefs_Stats')
% Design definitions for simple statistics
% D      - struct array of design definitions (see definition below)
%
% FORMAT D = spm_spm_ui('DesDefs_PET')
% Design definitions for PET/SPECT models
% D      - struct array of design definitions (see definition below)
%
% FORMAT D = spm_spm_ui('DesDefs_PET96')
% Design definitions for SPM96 PET/SPECT models
% D      - struct array of design definitions (see definition below)

%=======================================================================
% Design definitions specification for programers & power users
%=======================================================================
% Within spm_spm_ui.m, a definition structure, D, determines the
% various options, defaults and specifications appropriate for a
% particular design. Usually one uses one of the pre-specified
% definitions chosen from the menu, which are specified in the function
% actions at the end of the program (spm_spm_ui('DesDefs_Stats'),
% spm_spm_ui('DesDefs_PET'), spm_spm_ui('DesDefs_PET96')). For
% customised use of spm_spm_ui.m, the design definition structure is
% shown by the following example:
%
% D = struct(...
%       'DesName','The Full Monty...',...
%       'n',[Inf Inf Inf Inf],  'sF',{{'repl','cond','subj','group'}},...
%       'Hform',                'I(:,[4,2]),''-'',{''stud'',''cond''}',...
%       'Bform',                'I(:,[4,3]),''-'',{''stud'',''subj''}',...
%       'nC',[Inf,Inf],'iCC',{{[1:8],[1:8]}},'iCFI',{{[1:7],[1:7]}},...
%       'iGXcalc',[1,2,3],'iGMsca',[1:7],'GM',50,...
%       'iGloNorm',[1:9],'iGC',[1:11],...
%       'M_',struct('T',[-Inf,0,0.8*sqrt(-1)],'I',Inf,'X',Inf),...
%       'b',struct('aTime',1));
%
% ( Note that the struct command expands cell arrays to give multiple    )
% ( records, so if you want a cell array as a field value, you have to   )
% ( embed it within another cell, hence the double "{{"'s.               )
%
%                           ----------------
%
% Design structure fields and option definitions
% ==============================================
%
% D.Desname  - a string naming the design
%
% In general, spm_spm_ui.m accomodates four factors. Usually these are
% 'group', 'subject', 'condition' & 'replication', but to allow for a
% flexible interface these are dynamically named for different designs,
% and are referred to as Factor4, Factor3, Factor2, and Factor1
% respectively. The first part of the D definition dictates the names
% and number of factor levels (i.e. number of subjects etc.) relevant
% for this design, and also how the H (condition) and B (block)
% partitions of the design matrix should be constructed.
%
% D.n        - a 1x4 vector, indicating the number of levels. D.n(i)
%              for i in [1:4] is the number of levels for factor i.
%              Specify D.n(i) as 1 to ignore this factor level,
%              otherwise the number of levels can be pre-specified as a
%              given number, or given as Inf to allow the user to
%              choose the number of levels.
%
% D.sF       - a 1x4 cellstr containing the names of the four
%              factors. D.sF{i} is the name of factor i.
%
% D.b.aTime  - a binary indicator specifying whether images within F3
%              level (subject) are selected in time order. For time
%              order (D.b.aTime=1), F2 levels are indicated by a user
%              input "condition" string (input handled by spm_input's
%              'c' type). When (D.b.aTime=0), images for each F3 are
%              selected by F2 (condition). The latter was the mode of
%              SPM95 and SPM96. (SPM94 and SPMclassic didn't do
%              replications of conditions.)
%
% Once the user has entered the images and indicated the factor levels,
% a nScan x 4 matrix, I, of indicator variables is constructed
% specifying for each scan the relevant level of each of the four
% factors. I(n,i) is the level of factor i corresponding to image n.
% This I matrix of factor indicators is then used to construct the H
% and B forms of the design matrix according to the prescripton in the
% design definition D:
%
% D.Hform    - a string specifying the form of the H partition of the
%              design matrix. The string is evaluated as an argument
%              string for spm_DesMtx, which builds design matrix
%              partitions from indicator vectors.
%             (eval(['[H,Hnames] = spm_DesMtx(',D.Hform,');']))
%
% D.BForm   - a string specifying the form of the G partition.
%
% ( Note that a constant H partition is dropped if the B partition can   )
% ( model the constant effect.                                           )
%
% The next part of the design definition defines covariate options.
% Covariates are split into covariates (of interest) and nuisance
% variables. The covariates of interest and nuisance variables are put
% into the C & G partitions of the design matrox (the final design
% matrix is [H,C,B,G], where global nuisance covariates are appended to
% G). In SPM94/5/6 the design matrix was partitioned into effects of
% interest [H,C] and effects of no interest [B,G], with an F-test for
% no effects of interest and adjusted data (for effects of no interest)
% following from these partitions. SPM99 is more freestyle, with
% adjustments and F-tests specified by contrasts. However, the concept
% of effects of interest and of no interest has been maintained for
% continuity, and spm_spm_ui.m computes an F-contrast to test for "no
% effects of interest".
%
% D.nC       - a 1x2 vector: D.nC(1) is the number of covariates,
%              D.nC(2) the number of nuisance variables. Specify zero
%              to skip covariate entry, the actual number of
%              covariates, or Inf to let the user specify the number of
%              covariates. As with earlier versions, blocks of design
%              matrix can be entered. However, these are now treated as
%              a single covariate entity, so the number of
%              covariates.nuisance variables is now the number of items
%              you are prompted for, regardless of their dimension. (In
%              SPM95-6 this number was the number of covariate vectors
%              that could be entered.)
%
% D.iCC      - a 1x2 cell array containing two vectors indicating the
%              allowable covariate centering options for this design.
%              These options are defined in the body of spm_spm_ui.m,
%              in variables sCC & CFIforms. Use negative indices to
%              indicate the default, if any - the largest negative
%              wins.
%
% D.iCFI     - a 1x2 cell array containing two vectors indicating the
%              allowable covariate by factor interactions for this
%              design. Interactions are only offered with a factor if
%              it has multiple levels. The options are defined in the
%              body of spm_spm_ui.m, in variables sCFI & CFIforms. Use
%              negative indicies to indicate a default.
%
% The next part defines global options:
%
% D.iGXcalc  - a vector of possible global calculation options for
%              this design, as listed in the body of spm_spm_ui.m in
%              variable sGXcalc. (If other global options are chosen,
%              then the "omit" option is not offered.) Again, negative
%              values indicate a default.
%
% D.iGloNorm - a vector of possible global normalisation options for
%              this design, as described in the body of spm_spm_ui.m in
%              variable sGloNorm.
%
% D.iGMsca   - a vector of possible grand mean scaling options, as
%              described in the body of spm_spm_ui.m in variable
%              sGMsca. (Note that grand mean scaling is redundent when
%              using proportional scaling global flow normalisation.)
%
% D.iGC      - a vector of possible global covariate centering
%              options, corresponding to the descriptions in variable
%              iCC given in the body of spm_spm_ui.m. This is only
%              relevant for AnCova type global normalisation, and even
%              then only if you're actually interested in constraining
%              the values of the parameters in some useful way.
%              Usually, one chooses option 10, "as implied by AnCova".
%
% The next component specifies masking options:
%
% D.M_.T     - a vector defining the analysis threshold: Specify
%              "-Inf" as an element to offer "None" as an option. If a
%              real element is found, then absolute thresholding is
%              offered, with the first real value proffered as default
%              threshold. If an imaginary element is found, then
%              proportional thresholding if offered (i.e. the threshold
%              is a proportion of the image global), with the (abs of)
%              the first imaginary element proffered as default.
%
% D.M_.I     - Implicit masking? 0-no, 1-yes, Inf-ask. (This is
%              irrelevant for image types with a representation of NaN,
%              since NaN is then the mask value, and NaN's are always
%              masked.)
%
% D.M.X      - Explicit masking? 0-no, 1-yes, Inf-ask.
% 
%                           ----------------
%
% To use a customised design structure D, type spm_spm_ui('cfg',D) in the
% Matlab command window.
%
% The easiest way to generate a customised design definition structure
% is to tweak one of the pre-defined definitions. The following code
% will prompt you to select one of the pre-defined designs, and return
% the design definition structure for you to work on:
%
% D = spm_spm_ui(char(spm_input('Select design class...','+1','m',...
%       {'Basic stats','Standard PET designs','SPM96 PET designs'},...
%       {'DesDefs_Stats','DesDefs_PET','DesDefs_PET96'},2)));
% D = D(spm_input('Select design type...','+1','m',{D.DesName}'))
% =======1=========2=========3=========4=========5=========6=========7=========8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The following is a hacked form of SPM8's spm_spm_ui.m following
%%% the structure of Joe's original do_rfx but adding anything unique
%%% to spm8 - Dec 2009 DDW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-Condition arguments
%-----------------------------------------------------------------------
Action = 'CFG';

switch lower(Action)
case 'cfg'
%=======================================================================
% - C O N F I G U R E   D E S I G N
%=======================================================================
% spm_spm_ui('CFG',D)
switch (nargin)
  case nargin < 2  
    error('Not enough input arguments');
  case 2
    rfx_path = varargin{1};
    D        = varargin{2};
  otherwise
    error('You''ve specified too many inputs. 2 only please.');
end
%%% Move to rfx dir.
cd(rfx_path);

%-Option definitions
    %-------------------------------------------------------------------
    %-Generic factor names
    sF = {'sF1','sF2','sF3','sF4'};
    
    %-Covariate by factor interaction options
    sCFI = {'<none>';...                            %-1
            'with sF1';'with sF2';'with sF3';'with sF4';...         %-2:5
            'with sF2 (within sF4)';'with sF3 (within sF4)'};       %-6,7
    
    %-DesMtx argument components for covariate by factor interaction options
    % (Used for CFI's Covariate Centering (CC), GMscale & Global normalisation)
    CFIforms = {    '[]',       'C',    '{}';...            %-1
            'I(:,1)',       'FxC',  '{D.sF{1}}';...         %-2
            'I(:,2)',       'FxC',  '{D.sF{2}}';...         %-3
            'I(:,3)',       'FxC',  '{D.sF{3}}';...         %-4
            'I(:,4)',       'FxC',  '{D.sF{4}}';...         %-5
            'I(:,[4,2])',   'FxC',  '{D.sF{4},D.sF{2}}';... %-6
            'I(:,[4,3])',   'FxC',  '{D.sF{4},D.sF{3}}' };  %-7
    
    %-Centre (mean correction) options for covariates & globals            (CC)
    % (options 9-12 are for centering of global when using AnCova GloNorm) (GC)
    sCC = {     'around overall mean';...               %-1
            'around sF1 means';...                  %-2
            'around sF2 means';...                  %-3
            'around sF3 means';...                  %-4
            'around sF4 means';...                  %-5
            'around sF2 (within sF4) means';...         %-6
            'around sF3 (within sF4) means';...         %-7
            '<no centering>';...                    %-8
            'around user specified value';...           %-9
            '(as implied by AnCova)';...                %-10
            'GM';...                        %-11
            '(redundant: not doing AnCova)'}';          %-12
    %-DesMtx I forms for covariate centering options
    CCforms = {'ones(nScan,1)',CFIforms{2:end,1},''}';
    
    
    %-Global normalization options (options 1-7 match CFIforms)       (GloNorm)
    sGloNorm = {    'AnCova';...                        %-1
            'AnCova by sF1';...                 %-2
            'AnCova by sF2';...                 %-3
            'AnCova by sF3';...                 %-4
            'AnCova by sF4';...                 %-5
            'AnCova by sF2 (within sF4)';...            %-6
            'AnCova by sF3 (within sF4)';...            %-7
            'proportional scaling';...              %-8
            '<no global normalisation>'};               %-9
    
    %-Grand mean scaling options                                        (GMsca)
    sGMsca = {  'scaling of overall grand mean';...         %-1
            'scaling of sF1 grand means';...            %-2
            'scaling of sF2 grand means';...            %-3
            'scaling of sF3 grand means';...            %-4
            'scaling of sF4 grand means';...            %-5
            'scaling of sF2 (within sF4) grand means';...       %-6
            'scaling of sF3 (within sF4) grand means';...       %-7
            '(implicit in PropSca global normalisation)';...    %-8
            '<no grand Mean scaling>'   };          %-9
    %-NB: Grand mean scaling by subject is redundent for proportional scaling
    
    
    %-Global calculation options                                       (GXcalc)
    sGXcalc  = {    'omit';...                      %-1
            'user specified';...                    %-2
            'mean voxel value (within per image fullmean/8 mask)'}; %-3

%-Set factor names for this design
%-----------------------------------------------------------------------
sCC      = sf_estrrep(sCC,[sF',D.sF']);
sCFI     = sf_estrrep(sCFI,[sF',D.sF']);
sGloNorm = sf_estrrep(sGloNorm,[sF',D.sF']);
sGMsca   = sf_estrrep(sGMsca,[sF',D.sF']);

%-Get filenames & factor indicies
%-----------------------------------------------------------------------
%%% \d\d\D\D\D\d\d = 11JAN11*.img because otherwise spm will suck up stray con
%%% files leftover from prior RFX. Very very very evil behavior.
%%% Fixed again again to account for abnormal naming.
%%% Two ways of doing this ['^con_.*_.*\.img'] will only grab cons of form
%%% con_XXXXXXX_XXXX.img. Which should cover everything. Can also use the
%%% bbc.filter here to specify the formula, but so far I don't think one
%%% needs to. Only subject cons should have the form con_XXX_XXX.img.
%%% Thus the con_0002.img should never get included.
P     = spm_select('FPList',rfx_path,['^con_.*_.*\.img']); 
len   = size(P,1);
P     = cellstr(P);

%-Additional design parameters
%-----------------------------------------------------------------------
I     = [[1:len]' ones(len,3)]; % setup factor indices for DesMtx
nScan = size(I,1);				% number of observation
bL       = any(diff(I,1),1); 	%-Multiple factor levels?
	    % NB: bL(2) might be thrown by user specified f1 levels
	    %     (D.b.aTime & D.n(2)>1) - assumme user is consistent?
bFI      = [bL(1),bL(2:3)&~bL(4),bL(4),bL([2,3])&bL(4)];
	    %-Allowable interactions for covariates
	    %-Only offer interactions with multi-level factors, and
	    % don't offer by F2|F3 if bL(4)!


%-Build Condition (H) and Block (B) partitions
%===================================================================
H=[];Hnames=[];
B=[];Bnames=[];
eval(['[H,Hnames] = spm_DesMtx(',D.Hform,');'])
if rank(H)==nScan, error('unestimable condition effects'), end
eval(['[B,Bnames] = spm_DesMtx(',D.Bform,');'])
if rank(B)==nScan, error('unestimable block effects'), end

%-Drop a constant H partition if B partition can model constant
if size(H,2)>0 & all(H(:)==1) & (rank([H B])==rank(B))
    H = []; Hnames = {};
    warning('Dropping redundant constant H partition')
end

%-Covariate partition(s): interest (C) & nuisance (G) excluding global
%=======================================================================
nC = D.nC;			%-Default #covariates
C  = {[],[]}; Cnames = {{},{}};	%-Covariate DesMtx partitions & names
xC = [];			%-Struct array to hold raw covariates


dcname = {'CovInt','NusCov'};	%-Default root names for covariates
dstr   = {'covariate','nuisance variable'};

%GUIpos = spm_input('!NextPos');  %%%NO GFX

nc     = [0,0];
for i  = 1:2			% 1:covariates of interest, 2:nuisance variables

    %if isinf(nC(i)), nC(i)=spm_input(['# ',dstr{i},'s'],GUIpos,'w1'); end
    %%% NOGFX

    while nc(i) < nC(i)

	%-Create prompt, get covariate, get covariate name
            %-----------------------------------------------------------
            if nC(i)==1
                str=dstr{i};
            else
                str=sprintf('%s %d',dstr{i},nc(i)+1);
            end
            c      = spm_input(str,GUIpos,'r',[],[nScan,Inf]);
            if any(isnan(c(:))), break, end     %-NaN is dummy value to exit
            nc(i)  = nc(i)+1;                   %-#Covariates (so far)
            if nC(i)>1, tstr = sprintf('%s^{%d}',dcname{i},nc(i));
            else,       tstr = dcname{i}; end
            cname  = spm_input([str,' name?'],'+1','s',tstr);
            rc     = c;                         %-Save covariate value
            rcname = cname;                     %-Save covariate name	%-Save covariate name

        %-Interaction option? (if single covariate vector entered)?
            %-----------------------------------------------------------
            if size(c,2) == 1
                %-User choice of interaction options, default is negative
                %-Only offer interactions for appropriate factor combinations
                if length(D.iCFI{i})>1
                    iCFI = intersect(abs(D.iCFI{i}),find([1,bFI]));
                    dCFI = max([1,intersect(iCFI,-D.iCFI{i}(D.iCFI{i}<0))]);
                    iCFI = spm_input([str,': interaction?'],'+1','m',...
                        sCFI(iCFI),iCFI,find(iCFI==dCFI));
                else
                    iCFI = abs(D.iCFI{i});     %-AutoSelect default option
                end
            else
                iCFI = 1;
            end

        % Always offer "no centering" as default for design matrix blocks
            %-----------------------------------------------------------
            DiCC = D.iCC{i};
            if size(c,2)>1, DiCC = union(DiCC,-8); end
            if length(DiCC)>1
                %-User has a choice of centering options
                %-Only offer factor specific for appropriate factor combinations
                iCC = intersect(abs(DiCC),find([1,bFI,1]) );
                %-Default is max -ve option in D, overridden by iCFI if CFI
                if iCFI == 1, dCC = -DiCC(DiCC<0); else, dCC = iCFI; end
                dCC = max([1,intersect(iCC,dCC)]);
                iCC = spm_input([str,': centre?'],'+1','m',...
                    sCC(iCC),iCC,find(iCC==dCC));
            else
                iCC = abs(DiCC);    %-AutoSelect default option
            end
            %-Centre within factor levels as appropriate
            if any(iCC == [1:7]), c = c - spm_meanby(c,eval(CCforms{iCC})); end
            
            %-Do any interaction (only for single covariate vectors)
            %-----------------------------------------------------------
            if iCFI > 1             %-(NB:iCFI=1 if size(c,2)>1)
                tI        = [eval(CFIforms{iCFI,1}),c];
                tConst    = CFIforms{iCFI,2};
                tFnames   = [eval(CFIforms{iCFI,3}),{cname}];
                [c,cname] = spm_DesMtx(tI,tConst,tFnames);
            elseif size(c,2)>1          %-Design matrix block
                [null,cname] = spm_DesMtx(c,'X',cname);
            else
                cname = {cname};
            end

        %-Store raw covariate details in xC struct for reference
            %-Pack c into appropriate DesMtx partition
            %-----------------------------------------------------------
            %-Construct description string for covariate
            str = {sprintf('%s: %s',str,rcname)};
            if size(rc,2)>1, str = {sprintf('%s (block of %d covariates)',...
                        str{:},size(rc,2))}; end
            if iCC < 8, str=[str;{['used centered ',sCC{iCC}]}]; end
            if iCFI> 1, str=[str;{['fitted as interaction ',sCFI{iCFI}]}]; end
            
            tmp       = struct( 'rc',rc,    'rcname',rcname,...
                'c',c,      'cname',{cname},...
                'iCC',iCC,  'iCFI',iCFI,...
                'type',i,...
                'cols',[1:size(c,2)] + ...
                size([H,C{1}],2) +  ...
                size([B,C{2}],2)*(i-1),...
                'descrip',{str}             );
            if isempty(xC), xC = tmp; else, xC = [xC,tmp]; end
            C{i}      = [C{i},c];
            Cnames{i} = [Cnames{i}; cname];
            
        end % (while)
        
    end % (for)
    
clear c tI tConst tFnames
%spm_input('!SetNextPos',GUIpos); %%% NO GFX

%-Unpack into C & G design matrix sub-partitions
G = C{2}; Gnames = Cnames{2};
C = C{1}; Cnames = Cnames{1};


 %-Options...
    %===================================================================
    %-Global normalization options                             (GloNorm)
    %-------------------------------------------------------------------
    if length(D.iGloNorm)>1
        %-User choice of global normalisation options, default is negative
        %-Only offer factor specific for appropriate factor combinations
        iGloNorm = intersect(abs(D.iGloNorm),find([1,bFI,1,1]));
        dGloNorm = max([0,intersect(iGloNorm,-D.iGloNorm(D.iGloNorm<0))]);
        iGloNorm = spm_input('GloNorm: Select global normalisation','+1','m',...
            sGloNorm(iGloNorm),iGloNorm,find(iGloNorm==dGloNorm));
    else
        iGloNorm = abs(D.iGloNorm);
    end
    
    
    %-Grand mean scaling options                                 (GMsca)
    %-------------------------------------------------------------------
    if iGloNorm==8
        iGMsca=8;   %-grand mean scaling implicit in PropSca GloNorm
    elseif length(D.iGMsca)==1
        iGMsca = abs(D.iGMsca);
    else
        %-User choice of grand mean scaling options
        %-Only offer factor specific for appropriate factor combinations
        iGMsca = intersect(abs(D.iGMsca),find([1,bFI,0,1]));
        %-Default is max -ve option in D, overridden by iGloNorm if AnCova
        if iGloNorm==9, dGMsca=-D.iGMsca(D.iGMsca<0); else, dGMsca=iGloNorm; end
        dGMsca = max([0,intersect(iGMsca,dGMsca)]);
        iGMsca = spm_input('GMsca: grand mean scaling','+1','m',...
            sGMsca(iGMsca),iGMsca,find(iGMsca==dGMsca));
    end
    
    
    %-Value for PropSca / GMsca                                     (GM)
    %-------------------------------------------------------------------
    if iGMsca == 9                      %-Not scaling (GMsca or PropSca)
        GM = 0;                         %-Set GM to zero when not scaling
    else                                %-Ask user value of GM
        if iGloNorm==8
            str = 'PropSca global mean to';
        else
            str = [strrep(sGMsca{iGMsca},'scaling of','scale'),' to'];
        end
        GM = spm_input(str,'+1','r',D.GM,1);
        %-If GM is zero then don't GMsca! or PropSca GloNorm
        if GM==0, iGMsca=9; if iGloNorm==8, iGloNorm=9; end, end
    end
    
    %-Sort out description strings for GloNorm and GMsca
    %-------------------------------------------------------------------
    sGloNorm = sGloNorm{iGloNorm};
    sGMsca   = sGMsca{iGMsca};
    if iGloNorm==8
        sGloNorm = sprintf('%s to %-4g',sGloNorm,GM);
    elseif iGMsca<8
        sGMsca   = sprintf('%s to %-4g',sGMsca,GM);
    end
    

%-Global centering (for AnCova GloNorm)                             (GC)
%-----------------------------------------------------------------------
%-Specify the centering option for the global covariate for AnCova
%-Basically, if 'GMsca'ling then should centre to GM (iGC=11). Otherwise,
% should centre in similar fashion to AnCova (i.e. by the same factor(s)),
% such that models are seperable (iGC=10). This is particularly important
% for subject specific condition effects if then passed on to a second-level
% model. (See also spm_adjmean_ui.m) SPM96 (& earlier) used to just centre
% GX around its (overall) mean (iGC=1).

%-This code allows more general options to be specified (but is a bit complex)
%-Setting D.iGC=[-10,-11] gives the standard choices above

%-If not doing AnCova then GC is irrelevant
if ~any(iGloNorm == [1:7])
	iGC = 12;
	gc  = [];
else
	%-Annotate options 10 & 11 with specific details
	%---------------------------------------------------------------
	%-Tag '(as implied by AnCova)' with actual AnCova situation
	sCC{10} = [sCC{iGloNorm},' (<= ',sGloNorm,')'];
	%-Tag 'GM' case with actual GM & GMsca case
	sCC{11} = sprintf('around GM=%g (i.e. %s after grand mean scaling)',...
		GM,strrep(sCC{iGMsca},'around ',''));

	%-Constuct vector of allowable iGC
	%---------------------------------------------------------------
	%-Weed out redundent factor combinations from pre-set allowable options
	iGC = intersect(abs(D.iGC),find([1,bFI,1,1,1,1]));
	%-Omit 'GM' option if didn't GMsca (iGMsca~=8 'cos doing AnCova)
	if any(iGMsca==[8,9]), iGC = setdiff(iGC,11); end
	%-Omit 'GM' option if same as '(as implied by AnCova)'
	if iGloNorm==iGMsca, iGC = setdiff(iGC,11); end

	%-If there's a choice, set defaults (if any), & get answer
	%---------------------------------------------------------------
	if length(iGC)>1
		dGC = max([0,intersect(iGC,-D.iGC(D.iGC<0))]);
		str = 'Centre global covariate';
		if iGMsca<8, str = [str,' (after grand mean scaling)']; end
		iGC = spm_input(str,'+1','m',sCC(iGC),iGC,find(iGC==dGC));
	elseif isempty(iGC)
		error('Configuration error: empty iGC')
	end

	%-If 'user specified' then get value
	%---------------------------------------------------------------
	if iGC==9
		gc     = spm_input('Centre globals around','+0','r',D.GM,1);
		sCC{9} = sprintf('%s of %g',sCC{iGC},gc);
	else
		gc  = 0;
	end
end


%-Thresholds & masks defining voxels to analyse                   (MASK)
%=======================================================================
%GUIpos = spm_input('!NextPos'); %%% NO GFX

%-Analysis threshold mask
%-----------------------------------------------------------------------
%-Work out available options:
% -Inf=>None, real=>absolute, complex=>proportional, (i.e. times global)
M_T = D.M_.T; if isempty(M_T), M_T = [-Inf, 100, 0.8*sqrt(-1)]; end
M_T = {	'none',		M_T(min(find(isinf(M_T))));...
	'absolute',	M_T(min(find(isfinite(M_T)&(M_T==real(M_T)))));...
	'relative',	M_T(min(find(isfinite(M_T)&(M_T~=real(M_T)))))	};

%-Work out available options
%-If there's a choice between proportional and absolute then ask
%-----------------------------------------------------------------------
q = ~[isempty(M_T{1,2}), isempty(M_T{2,2}), isempty(M_T{3,2})];

if all(q(2:3))
	%tmp = spm_input('Threshold masking',GUIpos,'b',M_T(q,1),find(q));
	%%%%NOGFX
	q(setdiff([1:3],tmp))=0;
end

%-Get mask value - note that at most one of q(2:3) is true
%-----------------------------------------------------------------------
if ~any(q)				%-Oops - nothing specified!
	M_T = -Inf;
elseif all(q==[1,0,0])			%-no threshold masking
	M_T = -Inf;
else					%-get mask value
	if q(1),	args = {'br1','None',-Inf,abs(M_T{1+find(q(2:3)),2})};
	else,		args = {'r',abs(M_T{1+find(q(2:3)),2})}; end
	if q(2)
		M_T = spm_input('threshold',GUIpos,args{:});
	elseif q(3)
		M_T = spm_input('threshold (relative to global)',GUIpos,...
								args{:});
		if isfinite(M_T) & isreal(M_T), M_T=M_T*sqrt(-1); end
	else
		error('Shouldn''t get here!')
	end
end

%-Make a description string
%-----------------------------------------------------------------------
if isinf(M_T)
	xsM.Analysis_threshold = 'None (-Inf)';
elseif isreal(M_T)
	xsM.Analysis_threshold = sprintf('images thresholded at %6g',M_T);
else
	xsM.Analysis_threshold = sprintf(['images thresholded at %6g ',...
		'times global'],imag(M_T));
end



%-Implicit masking: Ignore zero voxels in low data-types?
%-----------------------------------------------------------------------
% (Implicit mask is NaN in higher data-types.)
type = getfield(spm_vol(P{1,1}),'dt')*[1,0]';
    if ~spm_type(type,'nanrep')
        switch D.M_.I
        case Inf,    M_I = spm_input('Implicit mask (ignore zero''s)?',...
                '+1','y/n',[1,0],1);        %-Ask
        case {0,1}, M_I = D.M_.I;           %-Pre-specified
        otherwise,  error('unrecognised D.M_.I type')
        end
        
        if M_I, xsM.Implicit_masking = 'Yes: zero''s treated as missing';
        else,   xsm.Implicit_masking = 'No'; end
    else
        M_I = 1;
        xsM.Implicit_masking = 'Yes: NaN''s treated as missing';
    end
   
    
%-Explicit mask images (map them later...)
    %-------------------------------------------------------------------
    switch(D.M_.X)
    case Inf,   M_X = spm_input('explicitly mask images?','+1','y/n',[1,0],2);
    case {0,1}, M_X = D.M_.X;
    otherwise,  error('unrecognised D.M_.X type')
    end
    if M_X, M_P = spm_select(Inf,'image','select mask images'); else, M_P = {}; end
    
%-Global calculation                                            (GXcalc)
%=======================================================================
iGXcalc = abs(D.iGXcalc);
%-Only offer "omit" option if not doing any GloNorm, GMsca or PropTHRESH
if ~(iGloNorm==9 & iGMsca==9 & (isinf(M_T)|isreal(M_T)))
	iGXcalc = intersect(iGXcalc,[2:size(sGXcalc,1)]);
end
if isempty(iGXcalc)
	error('no GXcalc options')
elseif length(iGXcalc)>1
	%-User choice of global calculation options, default is negative
	dGXcalc = max([1,intersect(iGXcalc,-D.iGXcalc(D.iGXcalc<0))]);
	iGXcalc = spm_input('Global calculation','+1','m',...
	    	sGXcalc(iGXcalc),iGXcalc,find(iGXcalc==dGXcalc));
else
	iGXcalc = abs(D.iGXcalc);
end

if iGXcalc==2				%-Get user specified globals
	g = spm_input('globals','+0','r',[],[nScan,1]);
end
sGXcalc = sGXcalc{iGXcalc};


% Non-sphericity correction (set xVi.var and .dep)
    %===================================================================
    xVi.I   = I;
    nL      = max(I);                                % number of levels
    mL      = find(nL > 1);                          % multilevel factors
    xVi.var = sparse(1,4);                           % unequal variances
    xVi.dep = sparse(1,4);                           % dependencies

    if length(mL) > 1

        % repeated measures design
        %---------------------------------------------------------------
        if spm_input('non-sphericity correction?','+1','y/n',[1,0],0)

            % make menu strings
            %-----------------------------------------------------------
            for i = 1:4
                mstr{i} = sprintf('%s (%i levels)',D.sF{i},nL(i));
            end
            mstr = mstr(mL);
            
            % are errors identical
            %-----------------------------------------------------------
            if spm_input('are errors identical','+1','y/n',[0,1],0)
                str    = 'unequal variances are between';
                [i j]  = min(nL(mL));
                i      = spm_input(str,'+0','m',mstr,[],j);

                % set in xVi and eliminate from dependency option
                %-------------------------------------------------------
                xVi.var(mL(i)) = 1;
                mL(i)          = [];
                mstr(i)        = [];
            end
            
            % are errors independent
            %-----------------------------------------------------------
            if spm_input('are errors independent','+1','y/n',[0,1],0)
                str        = ' dependencies are within';
                [i j]      = max(nL(mL));
                i          = spm_input(str,'+0','m',mstr,[],j);

                % set in xVi
                %-------------------------------------------------------
                xVi.dep(mL(i)) = 1;
            end

        end
    end
    
    %-Place covariance components Q{:} in xVi.Vi
    %-------------------------------------------------------------------
    xVi     = spm_non_sphericity(xVi);


%=======================================================================
% - C O N F I G U R E   D E S I G N
%=======================================================================
%spm('FigName','Stats: configuring',Finter,CmdLine);
%spm('Pointer','Watch');


%-Images & image info: Map Y image files and check consistency of
% dimensions and orientation / voxel size
%=======================================================================
fprintf('%-40s: ','Mapping files')                                   %-#
VY    = spm_vol(char(P));


%-Check compatability of images (Bombs for single image)
%-------------------------------------------------------------------
spm_check_orientations(VY);
fprintf('%30s\n','...done')                                      %-#


%-Global values, scaling and global normalisation
%===================================================================
%-Compute global values
%-------------------------------------------------------------------
switch iGXcalc, case 1
    %-Don't compute => no GMsca (iGMsca==9) or GloNorm (iGloNorm==9)
    g = [];
case 2
    %-User specified globals
case 3
    %-Compute as mean voxel value (within per image fullmean/8 mask)
    g     = zeros(nScan,1 );
    fprintf('%-40s: %30s','Calculating globals',' ')             %-#
    for i = 1:nScan
        str = sprintf('%3d/%-3d',i,nScan);
        fprintf('%s%30s',repmat(sprintf('\b'),1,30),str)%-#
        g(i) = spm_global(VY(i));
    end
    fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')       %-#
otherwise
    error('illegal iGXcalc')
end
rg = g;

fprintf('%-40s: ','Design configuration')                        %-#

%-Scaling: compute global scaling factors gSF required to implement proportional
% scaling global normalisation (PropSca) or grand mean scaling (GMsca),
% as specified by iGMsca (& iGloNorm)
%-----------------------------------------------------------------------
switch iGMsca, case 8
	%-Proportional scaling global normalisation
	if iGloNorm~=8, error('iGloNorm-iGMsca(8) mismatch for PropSca'), end
	gSF    = GM./g;
	g      = GM*ones(nScan,1);
case {1,2,3,4,5,6,7}
	%-Grand mean scaling according to iGMsca
	gSF    = GM./spm_meanby(g,eval(CCforms{iGMsca}));
	g      = g.*gSF;
case 9
	%-No grand mean scaling
	gSF    = ones(nScan,1);
otherwise
	error('illegal iGMsca')
end

%-AnCova: Construct global nuisance covariates partition (if AnCova)
%-------------------------------------------------------------------
if any(iGloNorm == [1:7])

    %-Centre global covariate as requested
    %---------------------------------------------------------------
    switch iGC, case {1,2,3,4,5,6,7}    %-Standard sCC options
        gc = spm_meanby(g,eval(CCforms{iGC}));
    case 8                  %-No centering
        gc = 0;
    case 9                  %-User specified centre
        %-gc set above
    case 10                 %-As implied by AnCova option
        gc = spm_meanby(g,eval(CCforms{iGloNorm}));
    case 11                 %-Around GM
        gc = GM;
    otherwise               %-unknown iGC
        error('unexpected iGC value')
    end


    %-AnCova - add scaled centred global to DesMtx `G' partition
    %---------------------------------------------------------------
    rcname     = 'global'; 
    tI         = [eval(CFIforms{iGloNorm,1}),g - gc];
    tConst     = CFIforms{iGloNorm,2};
    tFnames    = [eval(CFIforms{iGloNorm,3}),{rcname}];
    [f,gnames]  = spm_DesMtx(tI,tConst,tFnames);
    clear tI tConst tFnames

    %-Save GX info in xC struct for reference
    %---------------------------------------------------------------
    str     = {sprintf('%s: %s',dstr{2},rcname)};
    if any(iGMsca==[1:7]), str=[str;{['(after ',sGMsca,')']}]; end
    if iGC ~= 8, str=[str;{['used centered ',sCC{iGC}]}]; end
    if iGloNorm > 1
        str=[str;{['fitted as interaction ',sCFI{iGloNorm}]}]; 
    end
    tmp  = struct(  'rc',rg.*gSF,       'rcname',rcname,...
        'c',f,          'cname' ,{gnames},...
        'iCC',iGC,      'iCFI'  ,iGloNorm,...
        'type',         3,...
        'cols',[1:size(f,2)] + size([H C B G],2),...
        'descrip',      {str}       );

    G = [G,f]; Gnames = [Gnames; gnames];
    if isempty(xC), xC = tmp; else, xC = [xC,tmp]; end


elseif iGloNorm==8 | iGXcalc>1

    %-Globals calculated, but not AnCova: Make a note of globals
    %---------------------------------------------------------------
    if iGloNorm==8
        str = { 'global values: (used for proportional scaling)';...
                '("raw" unscaled globals shown)'};
    elseif isfinite(M_T) & ~isreal(M_T)
        str = { 'global values: (used to compute analysis threshold)'};
    else
        str = { 'global values: (computed but not used)'};
    end

    rcname ='global';
    tmp     = struct(   'rc',rg,    'rcname',rcname,...
        'c',{[]},   'cname' ,{{}},...
        'iCC',0,    'iCFI'  ,0,...
        'type',     3,...
        'cols',     {[]},...
        'descrip',  {str}           );

    if isempty(xC), xC = tmp; else, xC = [xC,tmp]; end
end


%-Save info on global calculation in xGX structure
%-------------------------------------------------------------------
xGX = struct(...
    'iGXcalc',iGXcalc,  'sGXcalc',sGXcalc,  'rg',rg,...
    'iGMsca',iGMsca,    'sGMsca',sGMsca,    'GM',GM,'gSF',gSF,...
    'iGC',  iGC,        'sGC',  sCC{iGC},   'gc',   gc,...
    'iGloNorm',iGloNorm,    'sGloNorm',sGloNorm);



%-Construct masking information structure and compute actual analysis
% threshold using scaled globals (rg.*gSF)
%-------------------------------------------------------------------
if isreal(M_T), M_TH =      M_T  * ones(nScan,1);   %-NB: -Inf is real
else,       M_TH = imag(M_T) * (rg.*gSF); end

if ~isempty(M_P)
    VM  = spm_vol(char(M_P));
    xsM.Explicit_masking = [{'Yes: mask images :'};{VM.fname}'];
else
    VM  = [];
    xsM.Explicit_masking = 'No';
end
xM     = struct('T',M_T, 'TH',M_TH, 'I',M_I, 'VM',{VM}, 'xs',xsM);


%-Construct full design matrix (X), parameter names and structure (xX)
%===================================================================
X      = [H C B G];
tmp    = cumsum([size(H,2), size(C,2), size(B,2), size(G,2)]);
xX     = struct(    'X',        X,...
    'iH',       [1:size(H,2)],...
    'iC',       [1:size(C,2)] + tmp(1),...
    'iB',       [1:size(B,2)] + tmp(2),...
    'iG',       [1:size(G,2)] + tmp(3),...
    'name',     {[Hnames; Cnames; Bnames; Gnames]},...
    'I',        I,...
    'sF',       {D.sF});


%-Design description (an nx2 cellstr) - for saving and display
%===================================================================
tmp = { sprintf('%d condition, +%d covariate, +%d block, +%d nuisance',...
        size(H,2),size(C,2),size(B,2),size(G,2));...
        sprintf('%d total, having %d degrees of freedom',...
        size(X,2),rank(X));...
        sprintf('leaving %d degrees of freedom from %d images',...
        size(X,1)-rank(X),size(X,1))                };
xsDes = struct( 'Design',           {D.DesName},...
    'Global_calculation',       {sGXcalc},...
    'Grand_mean_scaling',       {sGMsca},...
    'Global_normalisation',     {sGloNorm},...
    'Parameters',           {tmp}           );


fprintf('%30s\n','...done')                                      %-#



%-Assemble SPM structure
%===================================================================
SPM.xY.P    = P;            % filenames
SPM.xY.VY   = VY;           % mapped data
SPM.nscan   = size(xX.X,1); % scan number
SPM.xX      = xX;           % design structure
SPM.xC      = xC;           % covariate structure
SPM.xGX     = xGX;          % global structure
SPM.xVi     = xVi;          % non-sphericity structure
SPM.xM      = xM;           % mask structure
SPM.xsDes   = xsDes;        % description
%SPM.SPMid   = SPMid;        % version

%-Save SPM.mat and set output argument
%-------------------------------------------------------------------
fprintf('%-40s: ','Saving SPM configuration')                    %-#

if spm_matlab_version_chk('7') >=0
    save('SPM', 'SPM', '-V6');
else
    save('SPM', 'SPM');
end;
fprintf('%30s\n','...SPM.mat saved')                             %-#
varargout = {SPM};

%-Display Design report
%===================================================================
%fprintf('%-40s: ','Design reporting')                            %-#
%fname     = cat(1,{SPM.xY.VY.fname}');
%spm_DesRep('DesMtx',SPM.xX,fname,SPM.xsDes)
%fprintf('%30s\n','...done')     


%-End: Cleanup GUI
%===================================================================
%spm_clf(Finter)
%spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                 %-#
%spm('FigName','Stats: configured',Finter,CmdLine);
%spm('Pointer','Arrow')
fprintf('\n\n')
 
%NOW ESTIMATE!!!
spm_spm(SPM);
varargout = {SPM};

otherwise
%=======================================================================
% - U N K N O W N   A C T I O N
%=======================================================================
warning(['Illegal Action string: ',Action])


%=======================================================================
% - E N D
%=======================================================================
end

%=======================================================================
%- S U B - F U N C T I O N S
%=======================================================================

function str = sf_estrrep(str,srstr)
%=======================================================================
for i = 1:size(srstr,1)
	str = strrep(str,srstr{i,1},srstr{i,2});
end
