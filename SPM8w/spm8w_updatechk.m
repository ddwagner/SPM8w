function spm8w_updatechk(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: February 2013 - DDW
% ==============================================================================
% spm8w_updatechk()
% 
% Simple script to print out the revision number of various core SPM
% files used by spm8w in order to faciliate figuring out which need
% to be combed over for assimiliation in SPM8w.
% ==============================================================================
% CHANGE LOG:
% =======1=========2=========3=========4=========5=========6=========7=========8

%---Set defaults
spm8_path  = fileparts(which('spm.m'));
spm8w_path = fileparts(which('spm8w.m'));
addpath([spm8_path,'/toolbox/Seg/']);
addpath([spm8_path,'/toolbox/DARTEL/']);

edit_list = {'realign','slice_timing','uw_estimate','normalise',...
    'normalise_disp','preproc_run','maff8','preproc8','preproc_write8',...
    'dartel_template','dartel_norm_fun','dartel_warp'};
%---Print version match
[v,r] = spm('Ver');
fprintf('SPM8 revision: r%s\n', r);
for i = 1:length(edit_list)
    file8  = sprintf('spm_%s.m',edit_list{i}); 
    file8w = sprintf('spm8w_%s.m',edit_list{i}); 
    r8 = spm('Ver',file8);
    r8w = spm('Ver',file8w);
    fprintf('Revision: r%s | SPM8 file: %s\n', r8, file8);
    fprintf('Revision: r%s | SPM8 file: %s\n', r8w, file8w);
end