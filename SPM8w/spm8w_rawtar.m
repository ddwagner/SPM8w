function spm8w_nifticonverter(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: March, 2013 - DDW
% Created: March, 2013 - DDW
% ==============================================================================
% spm8w_rawtar(subjects);
%
% spm8w_rawtar is a simple script designed to tar/gzip directories containing
% raw data. It's purpose is to provide a simple and/or windows alternative to
% unix based rawtar command provided with SPM8w. 
%
% To use, simply provide a cell array of directory names or if called directly
% spm8w_rawtar will prompt you to select subjects. 
%
% spm8w_rawtar should be run from the directory containing raw data directories.
% Afterwards copy the resulting raw_subID.tar.gz files to a directory within
% your study path to use with spm8w_nifticonverter.m
% ==============================================================================
% CHANGE LOG:
% -Batched up the old spm8w_r2agui.
% =======1=========2=========3=========4=========5=========6=========7=========8

%---Setup
%--Get the current working dir
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
%--Goto root
cd(cwd); 

%--Input checks
switch (nargin)
  case 0 
    subjects = spm_select(Inf, 'dir',bbc.subjects,[],[cwd],bbc.filter);
  case 1
    subjects = varargin{1};
  otherwise
    error('You''ve specified too many inputs. Either 0 or 1 only please.');
end

%---Check for epty subjects in case of abort.
if isempty(subjects) %Check for empty subjects in case user aborted
    error('No subjects have been provided. Exiting...');
end

%---Clean up subjects list and cast into cell array
%--Remove dir and file extensions and raw_ prefix
subjects_tmp = spm_str_manip(subjects,'htr');
clear subjects
subjects(1:length(subjects_tmp(:,1)))={zeros}; %Preallocate subjects
for i = 1:size(subjects_tmp,1)
     subjects(i)={deblank(subjects_tmp(i,:))};
end

for i = 1:length(subjects)
    fprintf('Compressing subject: %s\n', subjects{i});
    tar(['raw_',subjects{i}], subjects{i});
    gzip(['raw_',subjects{i},'.tar']);
    delete(['raw_',subjects{i},'.tar']);
end