function data = spm8w_shuffle_check(epi_dir, file_stub, nslice)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: June, 2002 - Petr Janata
% ==============================================================================
% data = spm8w_shuffle_check(epi_dir, file_stub, nslice);
%
% Checks for shuffling of slices in EPI images.
%
% You must specify:
%  1) the EPI directory, e.g. '/data2/modfmri/28jan02PJ/epi/run1'
%  2) the root of the name of the EPI files, e.g. 'r28jan02PJ'
%  3) the number of slices, e.g. 27
%
% A color image will be displayed which shows the sum activity in each slice
% through time.  The color of each slice primarily reflects the number of
% voxels in that slice with signal.  Thus, shuffled or duplicated slices appear
% as a vertical line of shuffled colors.
%
% Specification of nslice is optional (spm8w_shuffle_check can figure it
% out on its own).
%
% v.1.0 06/20/02 Petr Janata
% ==============================================================================
% CHANGE LOG:
% -Overhauled code to level up shuffle check to SPM8 compatbility -DDW Dec/09
% =======1=========2=========3=========4=========5=========6=========7=========8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch (nargin)
  case nargin < 2 
    error('Too few paramters. You must specify the EPI directory and root filename');
  case nargin > 3
    error('Too many paramters. You must specify the EPI directory and root filename');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grab Files and determine dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P = spm_select('FPList',epi_dir,['^',file_stub,'.*\.nii']); %Load the bold nifti file.

%Figure out number of volumes
fprintf('Loading file...\n');
V = spm_vol(P);
nvol = length(V);

%Determine slices or check user input slices for accuracy
if (nargin < 3)
    nslice = V(1).dim(3);
elseif (nslice ~= V(1).dim(3))
    error (['User input slices (e.g. ', num2str(nslice), ') do not match slices in file (#slices = ',num2str(vol_tmp.dat.dim(3)),')']);
end

%Tell user what's going on
fprintf('Analyzing %d volumes beginning with %s in directory %s\n', nvol, file_stub, epi_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start processing the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = zeros(nslice, nvol);
for ivol = 1:nvol
    %fprintf(' .') Uses up too much space - Dec 2009 DDW
    for islice = 1:nslice    
        tmp = spm_slice_vol(V(ivol),spm_matrix([0 0 islice]),V(ivol).dim(1:2),0);
        data(islice,ivol) = sum(tmp(:));
    end
end

%Drop a new line
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct a figure (NB: in the Kellerton pipeline we 
% build our shuffle figure later so this section is disabled.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%imagesc(data);
%xlabel('Volume')
%ylabel('Slice')
%title(sprintf('%s/%s*.nii',epi_dir, file_stub))




