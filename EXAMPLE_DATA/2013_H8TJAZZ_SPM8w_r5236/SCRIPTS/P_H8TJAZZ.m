% ==============================================================================
% SPM8w r5236
% Parameters File
% Last update: February 2013 - DDW
% =======1=========2=========3=========4=========5=========6=========7=========8

%---DIRECTORY AND FILE SPECIFICATIONS
study_dir = '/flash/KELLERTON/DDW/2013_H8TJAZZ_SPM8w_r5236';
raw_dir   = 'SUBJECTS/RAW';
glm_dir   = 'GLM_H8TJAZZ';                                                      
rfx_dir   = 'RFX_H8TJAZZ';
confile   = 'CON_H8TJAZZ.m'; %Name of the confile for this study <optional>

%---ONSETS SPECIFICATIONS
onsets_dirname    = 'H8TJAZZ';
p.onsets_prefix   = []; 	 %if subject specific onsets: [p.subj,'_']       
p.duration_prefix = []; 	 %if subject specific durations: [p.subj,'_']
p.para_prefix     = []; 	 %if subject specific parametrics: [p.subj,'_']
p.onsets_ext      = '.txt';  %if no extension leave blank (ie: [];)    
    
%---CONDITIONS - Name of the conditions will be used to find the onset/reg files
p.event_cond = {'human';'animal';'vegetable';'mineral'};
p.block_cond = {};
p.reg_cond   = {};
    
%---SCANNING PARAMETERS
p.nses = 2;   %Number of runs                                                       
p.nTR  = 120; %Number of TRs per run                           
p.TR   = 2.5; %Length of TR in seconds                           
p.nTE  = 36;  %Number of slices                           
   
%---PREPROCESSING AND QA ROUTINES - 1=YES 0=NO
p.shuffle       = 1; %Calculate shuffle check 
p.despike       = 0; %Despike (AFNI's 3dDespike must be in path)
p.slicetime     = 1; %Slicetime correction                                     
p.realign       = 1; %Motion correction     
p.unwarp        = 1; %Unwarping (correct field inhomogeneties)                                
p.normalize     = 1; %Normalize                                     
p.smoothing     = 1; %Smooth 
p.smooth_kernel = 8; %Size of smoothing kernel in FWHM                                 
p.snr           = 1; %Calculate SNR
p.slices        = 1; %Calculate slice noise
 
%---MODEL DESIGN AND ESTIMATION PARAMETERS - 1=YES 0=NO
p.normtype      = 'OLD'; %Which norm routine to input to model (OLD|DARTEL|VBM8)
p.include_run   = 'all'; %Specify run to model: 'all' or run# (e.g. 2)
p.polyord       = 1;     %Order of polynomials: 0:Const|1:Linear|2:Quad|3:Cubic
p.move          = 1;     %Include motion regressors?
p.outliers      = 1;     %Include outlier scans as nuissance? (use spm8w_art.m)
p.demean        = 0;     %Demean condition regressors ala SPM99?
p.durtime       = 0;     %Event duration in seconds? (will be divided by p.TR)
p.disable_orth  = 0;     %Disable w/i trial orth (1) or enable it (0, SPM default)

%---SEGMENTATION, DARTEL AND VBM8 PARAMETERS - 1=YES 0=NO
p.coreg2epi     = 1;        %Coregister Anat to EPI data (mean of session 1)?
p.biasreg       = 0.0001;   %Intensity bias regularisation
							%(0 | 0.00001 | *0.0001 | 0.001 | 0.01 | 0.1 | 1 | 10) 
p.mrf           = 2;        %Segmentation cleanup via MRF (0 | 1 | *2 | 3 | 4)
p.warpreg       = 4;        %Warp regularisation (spm default:4|can use 0.4 for T1)                        
p.sampling      = 3;        %Sampling distance (1|2|*3) < 3mm will slow computation 
p.tissues       = 2;        %Tissue classes for template (1:GM|*2:GM+WM|3:GM+WM+CSF)
p.dartelsmooth  = 6;        %Smooting kernel for DARTEL to MNI (GM,WM,EPI only)
p.whichseg      = 'DARTEL'; %Segmentation type? (DARTEL | VBM8)
p.whichtemplate = 'DARTEL'; %Template style? VBM8 uses IXI550 (DARTEL | VBM8)
p.whichnorm     = 'DARTEL'; %Target for EPI normalization? (DARTEL | VBM8)
    
%---POST-PREPROCESSING CLEANUP
%--0 = keep all files (waste of space)
%--1 = delete all but bold, swuabold/swubold & uabold/ubold/abold (req for DARTEL) 
%--2 = delete all but bold & swuabold/swubold (if not using DARTEL)
p.cleanup = 1;    
  
%---MONTE CARLO SIMULATIONS (ALPHASIM)
%--Optional specifications for using the spm8w wrapper 
%--for Afni's AlphaSim (spm8w_alphasim.m).
p.as_mask  = 'bigmask_3x3x3_nocsf'; %Mask must be in path and Nifti
p.as_thr   = '0.001';               %Threshold for simulations (0.001|0.005|etc.)
p.as_iter  = '10000';               %Num  of iterations (1,000|10,000|etc.) 
p.as_fast  = 0;                     %1=YES 0=NO; Speeds up simulations x2
p.as_fname = 'as_H8TJAZZ.txt';      %Filename for alphasim output.

%---HRF
%--OPTIONS:'hrf'
%          'hrf (with time derivative)'
%          'hrf (with time and dispersion derivatives)'
%          'Fourier set'
%          'Fourier set (Hanning)'
%          'Gamma functions'
%          'Finite Impulse Response'
p.hrf           = 'hrf';  %Set to Finite Impulse Response for MIXED designs.
%p.hrfwindow    = 20;     %For FIR: Length of HRF window 
%p.hrfbasis     = 8;      %For FIR: Number of bins per window

%---ADVANCED OPTIONS   
p.hpf         = Inf;     %HPF inf=no cutoff|otherwise cutoff in secs i.e. 128)                               
p.autocorr    = 'none';  %Autocorrelation correction (none | 'AR(1) + w')                                   
p.time        = 'scans'; %Onsets specified in 'scans' or 'secs' 
p.refslice    = 1;       %Reference slice for slice time correction
p.realign_rtm = 1;       %Realign to 1st only (0) or 2pass, 1st then sess mean (1)
p.overridedel = 1;       %Override check when deleting previous preprocessing (0 | 1)

% ==============================================================================
% DO NOT EDIT BEYOND THIS POINT OR ELSE UNIVERSE RESETS
% ==============================================================================

%---DEFAULTS
subjects_dir    = 'SUBJECTS';
anatomy_dir     = 'ANATOMY';   
anatvbm8_dir    = 'ANATOMY_VBM8';
functional_dir  = 'FUNCTIONAL';   
onsets_dir      = 'ONSETS';
results_dir     = 'RESULTS';                     
rfxroot         = 'RFX';
dartel_template = 'DARTEL_TEMPLATE';
scripts         = 'SCRIPTS';
p.bold          = 'bold';  
p.mprage        = 'mprage';
p.logsuffix     = 'bold';

%---ADJUST p.nTR TO NUMBER OF TRs PER SESSION
if size(p.nTR,2) == 1 && p.nses > 1
    p.nTR = repmat(p.nTR,p.nses,1)';
end

%---SLICE ORDER 
sliceorder=[];
%--Formula for interleaved on Philips Achieva 3T -DDW July 2009
for i = 1:round(sqrt(p.nTE))
    sliceorder = [sliceorder i:round(sqrt(p.nTE)):p.nTE];
end
%--Formula for regular interleaved sequences
%sliceorder=[1:2:p.nTE 2:2:p.nTE];
p.sliceorder = sliceorder;
    
%---WRITE PARAMETERS
%--Interpolation set to 7 bsplines gives a less blury normalize EPI.
%--After smoothing that extra spatial resolution is lost.
%--However, for studies with low smoothing, it might be worth bumping up
%--the interp_r and interp_w.  
p.wrap_r     = [0 0 0];     %SPM8 default [0 0 0] | spm2 scripts was [0 1 0]
p.wrap_w     = [0 0 0];     %SPM8 default [0 0 0] | spm2 scripts was [0 1 0]
p.interp_r   = 7;           %SPM8 default 4       | spm2 scripts was 7
p.interp_w   = 7;           %SPM8 default 1       | spm2 scripts was 7
p.voxsize    = [3 3 3];     %SPM8 default [2 2 2] | don't resample voxel sizes
p.boundbox   = [-78 -117 -72;...
                 78 76 86]; %SPM8 default is [-78 -112 -50; 78 76 85] 
							%but that's for 2x2x2 (I think!)

%---GLM MASK
%--Main mask for Dartmouth is bigmask.nii (1x1x1 or 3x3x3)
p.mask = 'bigmask.nii';

%---PATHS
%--p.root: The root dir for your study
%--p.rawdir: The dir where raw data (subid.tar.gz) are stored.
%--p.subdir: The dir where subjects are stored   
%--p.anat: The dir name containing Anatomical scans 
%--p.vbm8: The dir name containing Anatomical scans (VBM8 output)   
%--p.func: The dir name containing Functional scans  
%--p.res: The dir where the analysis results will be written     
%--p.glm: The dir where the SPM.mat will be written to 
%--p.rfxroot: The dir where all RFX will be stored
%--p.rfx: The dir for this specific rfx
%--p.dartdir: The dir where DARTEL templates are stored
%--p.onsets: The dir where the onset files are stored.
%--p.confile: Name of the confile to use
p.root    = study_dir;
p.rawdir  = fullfile(p.root, raw_dir);
p.subdir  = fullfile(p.root, subjects_dir);
p.nifti   = fullfile(p.subdir, p.subj, 'NIFTI');
p.anat    = fullfile(p.subdir, p.subj, anatomy_dir);
p.vbm8    = fullfile(p.subdir, p.subj, anatvbm8_dir); 
p.func    = fullfile(p.subdir, p.subj, functional_dir);     
p.res     = fullfile(p.subdir, p.subj, results_dir);                       
p.glm     = fullfile(p.res, glm_dir);                 
p.rfxroot = fullfile(p.root, rfxroot);                             
p.rfx     = fullfile(p.rfxroot, rfx_dir);                          
p.dartdir = fullfile(p.root, subjects_dir, dartel_template);       
p.onsets  = fullfile(p.root, onsets_dir, onsets_dirname);
try
	p.confile = fullfile(p.root, scripts, confile);     
catch
end