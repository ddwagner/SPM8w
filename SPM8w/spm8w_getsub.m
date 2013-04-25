function subjects = spm8w_getsub(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: September, 2012 - DDW
% ==============================================================================
% spm8w_getsubjects()
%
% spm8w_getsubjects gets subjects with a smile.
% ==============================================================================
% CHANGE LOG:
% - 
% =======1=========2=========3=========4=========5=========6=========7=========8

%%%GET DEFAULTS
% Set the current working dir
cwd = pwd;
% Load user defaults, otherwise load vanilla defaults
if exist(which('spm8w_userdefaults.m'), 'file')
    [def_dir,def_file] = fileparts(which('spm8w_userdefaults.m'));
else
    [def_dir,def_file] = fileparts(which('spm8w_defaults.m'));
end
% Goto dir, eval defaults, check for errors, come back to cwd
cd(def_dir); 
try
    eval(def_file); 
catch
    defmsg = sprintf('The defaults file %s.m in %s is not evaluable, please check your syntax\n',def_file,def_dir); 
    error_reporter(cwd, defmsg, 1);
    error('Exiting...');
end
%%% Goto root
cd(cwd); 

%%%GET SUBJECTS
%%% Get subjects using a regexp to search only for directories
%%% that match the pattern in defaults or userdefaults - April 2010 DDW
%%% Unlike SPM2, SPM8's spm_str_manip has problems with directories so have to 
%%% use h to remove trailing slash before using TR - Dec 2009 DDW 
subjects_tmp = spm_select(Inf, 'dir',bbc.subjects,[],[cwd,'/SUBJECTS'],bbc.filter);
subjects_tmp = spm_str_manip(subjects_tmp, 'htr');  
if isempty(subjects_tmp)
    error(bbc.error_subjects);
end
subjects(1:length(subjects_tmp(:,1)))={zeros}; %Preallocate subjects
for i = 1:length(subjects_tmp(:,1))
    subjects(i)={deblank(subjects_tmp(i,:))};
end

%%%REPORT LAST ERROR MESSAGE AND DUMP USER BACK TO STUDYROOT
function error_reporter(rootdir, bbc_error, sparelines)
    cd(rootdir);
    fprintf('\n-----------------------------ERROR!-----------------------------\n%s\n',bbc_error);
    fprintf('----------------------------------------------------------------\n');
    err = lasterror;
    fprintf('\nLast error message: %s\n',err.message);
    for i = 1:length(err.stack)-sparelines
        fprintf('\nError occurred in function %s, at line %d of file %s\n',err.stack(i).name,err.stack(i).line,err.stack(i).file);
    end
return