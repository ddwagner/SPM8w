function vxdata = spm8w_rawvx(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: April, 2012 - DDW
% ==============================================================================
% spm8w_rawvx([xyz],'boldfile.nii')
%
% Simple function to extract raw timeseries from a single MNI coordinate
% XYZ. Requires at least an MNI coordinate. If no bold file is specified
% it will allow you to select one. 
%
% Example: 
% vx = spm8w_rawvox([21,-6,21],'swuabold1.nii');
% vx = spm8w_rawvox([21,-6,21]);
% ==============================================================================
% CHANGE LOG:
% 
% =======1=========2=========3=========4=========5=========6=========7=========8

%---Input checks
switch (nargin)
  case nargin < 1
    error('Please specify a coordinate: [x,y,z]');
  case 1
    XYZ     = varargin{1};
    vfile   = [];
  case 2
    XYZ     = varargin{1};
    vfile   = varargin{2};
  otherwise 
    error('Too many paramters.'); 
end

%---Load volumes
if isempty(vfile)
    vfile = spm_select(Inf,'image','Please select a normalized bold file');
    [a,b,c] = fileparts(vfile(1:end-2));
    vols = spm_vol(spm_select('FPList',a,['^',b,'.*\',c]));
else
    [a,b,c] = fileparts(vfile);
    if isempty(a) 
        a = './';
    end
    vols = spm_vol(spm_select('FPList',a,['^',b,'.*\',c]));
end

%---Convert XYZ to vx index and extact timeseries
vox    = inv(vols(1).mat) * [XYZ 1]';
vxdata = zeros(length(vols),1);
fprintf('Extracting timeseries at voxel %d,%d,%d: ',XYZ);
for i = 1:length(vxdata)
    vxdata(i) = spm_sample_vol(vols(i),vox(1),vox(2),vox(3),1);
    fprintf('.');
end
fprintf('\n');
