% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: August, 2010
% ==============================================================================
% spm8path.m
%
% Simple script to add paths for SPM8/SPM8w/SPM8 tools and inform the user of
% their version and file location. Can be copied and modified to setup
% different revisions of SPM8/SPM8w depending on the study currently being 
% worked on.
%
% To run on matlab start: matlab -nosplash -nodesktop -r spm8path
% =======1=========2=========3=========4=========5=========6=========7=========8
fprintf('=== Adding paths for SPM8/SPM8w/SPM8 tools ... \n');

%SPM8 and SPM8w path vars (edit these)
spm8_path  = '/afs/dbic.dartmouth.edu/usr/kelley/kelley/kelley-tools/SPM8/SPM8_5236';
spm8w_path = '/afs/dbic.dartmouth.edu/usr/kelley/kelley/kelley-tools/SPM8/SPM8w/SPM8w';

%Optional tools paths (edit these)
tools_path   = '/afs/dbic.dartmouth.edu/usr/kelley/kelley/kelley-tools/SPM8/TOOLS/';
r2agui_path  = [tools_path, 'R2AGUI_2.7'];
xjview_path  = [tools_path, 'XJVIEW_8.11'];
conn_path    = [tools_path, 'CONN_12'];
mni2tal_path = [tools_path, 'MNI2TAL'];
snpm8b_path  = [tools_path, 'SNPM8b'];
glmflex_path = [tools_path, 'GLM_Flex'];
glmflex_uts  = [tools_path, 'GLM_Flex_Utilities'];

%Add all paths
addpath(spm8_path, spm8w_path) 
addpath(r2agui_path, xjview_path, conn_path, mni2tal_path, snpm8b_path, ...
        glmflex_path, glmflex_uts)

%Check SPM8 revision
[r,v] = spm('Ver');
spm8loc = fileparts(which('spm'));

%Check SPM8w location
spm8wloc = fileparts(which('spm8w'));

%Print path checks
fprintf(['=== You are using %s revision %s located in ... \n' ...
         '=== %s \n'], r, v,spm8loc);
fprintf(['=== You are using SPM8w located in ... \n' ...
         '=== %s \n'], spm8wloc);
fprintf(['=== You are using fMRI tools located in ... \n' ...
         '=== %s \n'], tools_path);

clear r v *_path glmflex_uts spm8loc spm8wloc