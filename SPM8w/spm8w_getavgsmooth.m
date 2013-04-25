function mean_fwhm = spm8w_getavgsmooth(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: September, 2012 - DDW
% ==============================================================================
% spm8w_getavgsmooth(para_file,subjects)
% 
% spm8w_getavgsmooth extracts the average smoothness from the residuals for
% a given contrast
% ==============================================================================
% CHANGE LOG:
% =======1=========2=========3=========4=========5=========6=========7=========8

%---Input checks
switch (nargin)
  case 0 
    p         = spm8w_getp('P');
    subj      = spm8w_getsub;
  case 1
    p         = spm8w_getp('P',[],varargin{1});
    subj      = spm8w_getsub;
  case 2
    p         = spm8w_getp('P',[],varargin{1});
    subj      = varargin{3};    
  otherwise
    error('You''ve specified too many inputs. Specify only a subjid.');
end

%---Get estimated FWHM of residuals from SPM.mat file
fprintf('Loading FWHM for %d subjects...\n',length(subj));
fwhm_list = [];
for i = 1:length(subj)
    p       = spm8w_getp('P', subj{i},p.para_file);
    spm_mat = load(fullfile(p.glm,'SPM.mat'));
    fwhm_list(i,:) = spm_mat.SPM.xVol.FWHM * abs(spm_mat.SPM.xVol.M(1,1));
end
mean_fwhm = mean(fwhm_list);
min_fwhm  = min(fwhm_list);
max_fwhm  = max(fwhm_list);
fprintf('The mean smoothness across all subjects is [%2.2f,%2.2f,%2.2f]\nwith a min of [%2.2f,%2.2f,%2.2f] and a max of [%2.2f,%2.2f,%2.2f]\n',...
        mean_fwhm,min_fwhm,max_fwhm);
    
    
    
    
    
    
    