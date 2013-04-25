function params = spm8w_getp(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: May, 2012 - DDW
% ==============================================================================
% spm8w_getp(type, subjID, para_file)
%
% spm8w_getp puts the parameters structure into the workspace. Optionally takes
% subjectID as well. 
%
% Type is either 'P','D','R','ROI','VOI'. Default is 'P'.
%
% example: p = spm8w_getp('P','01feb01aa')
% ==============================================================================
% CHANGE LOG:
% - Removed this function out of spm8w.m so it can be used by other scripts. 
% =======1=========2=========3=========4=========5=========6=========7=========8

%---Input checks
switch (nargin)
  case 0 
    type      = 'P';
    p.subj    = [];    
  case 1
    type      = varargin{1};
    p.subj    = [];
  case 2
    type      = varargin{1};
    p.subj    = varargin{2};   
  case 3
    type      = varargin{1};
    p.subj    = varargin{2};   
    para_file = varargin{3};
  otherwise
    error('You''ve specified too many inputs.');
end

%---GET DEFAULTS
%--Set the current working dir
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

%---SET MESSAGE
switch type
    case 'P'
        msg       = bbc.pfile;
        errmsg    = bbc.error_pfile;
        err_eval  = bbc.error_pfile_eval;
        paramstok = 'p';
    case 'D'
        msg       = bbc.dfile;
        errmsg    = bbc.error_dfile;
        err_eval  = bbc.error_dfile_eval;
        paramstok = 'd';
    case 'R'
        msg       = bbc.pfile;
        errmsg    = bbc.error_pfile;
        err_eval  = bbc.error_pfile_eval;
        paramstok = 'p';
    case 'ROI'
        msg       = bbc.roifile;
        errmsg    = bbc.error_roifile;
        err_eval  = bbc.error_roifile_eval;
        paramstok = 'r';
    case 'VOI'
        msg       = bbc.voifile;
        errmsg    = bbc.error_voifile;
        err_eval  = bbc.error_voifile_eval;
        paramstok = 'v';
end

%---GET PARAMS
%--Get and eval params.
if ~exist('para_file','var')
    msg = 'Please select your Parameters file.';
    para_file = spm_select(1,['^',type,'_.*\.m$'],msg,[],[cwd,'/SCRIPTS']);
end
%--Check for indecisive user
if isempty(para_file)
    error(errmsg);
end
%--Try to eval Pfile, come back to cwd if fail
cd(spm_str_manip(para_file,'h'));
try
    eval(spm_str_manip(para_file,'tr')); 
catch
    error_reporter(cwd, err_eval, 1);
    error('Exiting...');
end
eval(['cd(',paramstok,'.root)']); 
eval(sprintf('%s.para_file = ''%s'';',paramstok, para_file));  %Assign para-file location
eval(sprintf('params = %s;',paramstok))

%---Report last error message and dump user back to studyroot
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