function spm8w_alphasim(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: January, 2010 - DDW
% ==============================================================================
% spm8w_alphasim();
% 
% spm8w_alphasim is a matlab wrapper for AFNI's alphasim allowing for the use 
% of alphasim with SPM8 generated data. spm8w_alphasim requires specifications
% from your parameters file (see P_H8TJAZZ.m for example). Currently 
%
% spm8w_alphasim has been tested with a 3x3x3mm version of bigmask.nii 
% and masks of the amygdala (in nifti format). It appears that the native
% support for nifti is good and there's no need to convert to brik format.
% Volume information is derived from mask while smoothness estimates require 
% you load the SPM.mat for the analysis you wish to run the monte carlos for.
% ==============================================================================
% CHANGE LOG:
% -Allowed specification of mask. -DDW Feb/10
% -Paramaters now come for P file. -DDW Feb/10
% -Now loads nifti masks (no longer supports brik). -DDW July/10
% -Fixed for latest alphasim, some of the output parsing might need
% tweaking depending on your version of AFNI. -DDW March/12
% -Now gets mean residuals across all subject from 1st level GLM (more
% valid than 2nd level which tends to be overly smooth). -DDW Sept/12
%
% DETAILS: Smoothness is in voxels not mm. Alphasim wants it in mm. So need
% to multiply voxel x voxel size to mm. Do not use the RESMS.img.You need the 
% actual residual images (which spm deletes after GLL estimation). However, 
% SPM saves the fwhm into the spm structure (see spm_spm.m).
% =======1=========2=========3=========4=========5=========6=========7=========8


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear SPM;                               %Just in case
p.subj = [];                             %Or else eval will crash.

switch (nargin)
  case 0 
    p = spm8w_getp('P');
  case 1
    try
        p = spm8w_getp('P',[],varargin{1});
    catch
        error(['Please be sure to provide the full path to your parameters file, or leave ' ...
              'this argument blank to select via a window.']);
    end      
    otherwise
        error('You''ve specified too many inputs. Only specify ROI_para.m file.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set Defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask = which([p.as_mask,'.nii']);                     %Hack to get path to mask (e.g. bigmask_3x3x3+orig)
pthr = p.as_thr;                                      %Threshold for simulation (usually 0.001 or 0.005)
iter = p.as_iter;                                     %Numer of iterations
fwhm = spm8w_getavgsmooth(p.para_file);               %Multiply by voxel size to get fwhm in mm.
diary_name = p.as_fname;                              %Sets the output filename to p.as_fname
if p.as_fast                                          %Check if user wants SPEEEEEEEED!
    fast = '-fast';
else
    fast = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run AlphaSim!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Build command
AlphaSim = ['AlphaSim -mask ',mask,' -iter ',num2str(iter),' -pthr ',num2str(pthr)...
    ' -fwhmx ',num2str(fwhm(1)),' -fwhmy ',num2str(fwhm(2)),' -fwhmz ',num2str(fwhm(3)),' -quiet ',fast];
tmpdate=datestr(now);
%Running AlphaSim
fprintf('==========Running AFNI''s AlphaSim...\n');
fprintf('Using mask at: %s\n',mask);
fprintf('Threshold is: %s\n',num2str(pthr))
fprintf('FWHM (XYZ) of residuals is: %smm %smm %smm\n',num2str(fwhm(1)),num2str(fwhm(2)),num2str(fwhm(3)));
fprintf('Running %s simulations, this may take awhile...\n',num2str(iter));
fprintf('==============================================================\n');
%Save output to diary file (easy but not elegant).
%%% Setup basic stats name
%%% Delete previous file if exists
try
    [s,w]=unix(['rm ',diary_name]);   
end
diary (diary_name);
system(AlphaSim);
diary off;
%%%Parse the diaryfile
%Super cludge! might break with alphasim changes?
[cluster,pvalue]=textread(diary_name,'%n%*n%*n%*n%*n%n','delimiter','\t','headerlines',35);
index = find(pvalue<0.050001);  %find the pvalue a bit below 0.05 (allows for 0.05 exactly)
sigclus = cluster(index(1));        %clustersize at that location.
diary (diary_name);
fprintf('\n==============================================================\n');
fprintf('Output of AlphaSim saved to file %s\n',[p.root,diary_name]);
fprintf('Whole brain correction at %s requires clusters of %d voxels or greater\n',num2str(pthr),sigclus);
diary off;

        


