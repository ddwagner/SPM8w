function spm8w_art(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: January, 2010 - DDW
% ==============================================================================
% spm8_art('sub_id','parameter_file')
%
% spm8w_art is an SPM8 adaptation of our older art_batch.m (ddw) which used
% art.m (written by Sue Whitefield, Shay Mozes, Paul Mazaika, and Jeff Cooper 
% and modified by Joe Moran) which formerly worked with SPM2w but I've been 
% slow to bring it into SPM8w.
%
% spm8w_art uses an earlier version of art_detect in order to detect outliers 
% due to motion or global signal. The volumes of these outliers are saved to 
% a textfile in the subID/FUNCTIONAL directory which then gets fed into 
% GLM computation in order to downweight those volumes.
%
% First argument is subject id, second can be full path to
% the params file, or else a window will prompt you for it.
% ==============================================================================
% CHANGE LOG:
% -SPM2w version - DDW May/07
% -SPM8w version - DDW Mar/12
% -Prints output to PDF - DDW Mar/12
% =======1=========2=========3=========4=========5=========6=========7=========8

%---Setup
%--Get the current working dir
cwd = pwd;
% Load user defaults, otherwise load vanilla defaults
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
cd(cwd); 

%---Input checks
switch (nargin)
  case 0 
    error('Please specify a subject id, e.g. spm8w_art(''01jan09ab'')');
  case 1
    p.subj = varargin{1};
    p      = spm8w_getp('P',p.subj);
    cd(p.root)
  case 2
    p.subj    = varargin{1};
    para_file = varargin{2};
    p         = spm8w_getp('P',p.subj,para_file);
    cd(p.root)
  otherwise
    error('You''ve specified too many inputs. Either 1 or 2 only please.');
end
  
try
  cd(fullfile(p.subdir,p.subj));
catch
    subjerror = sprintf(['I can''t seem to get to the subject''s dir.\nAre you sure ' ...
        'your subject dir (%s)\nand subject (%s) names are correct?'], p.subdir, p.subj);
    error_reporter(cwd, subjerror, 1);
    error('Exiting...');
end

%---Set Defaults
spm('Defaults','fmri')
%warning off MATLAB:divideByZero;
defaults.printstr = 'print -dpsc2 -painters -append -noui tmp_art.ps'; 
%--Delete any leftover PS files
try
    [s,w] = system(spm8w_osbabel('rm tmp_*.*'));
catch
end

%---Determine boldtok and rptok and add to p structure
p.boldtok = []; p.rptok = [];
if p.despike;   p.boldtok = ['d',p.boldtok]; p.rptok = ['d',p.rptok]; end
if p.slicetime; p.boldtok = ['a',p.boldtok]; p.rptok = ['a',p.rptok]; end
if p.unwarp;    p.boldtok = ['u',p.boldtok]; end
if p.normalize; p.boldtok = ['w',p.boldtok]; end
if p.smoothing; p.boldtok = ['s',p.boldtok]; end
p.boldtok = [p.boldtok, p.bold];
p.rptok   = [p.rptok, p.bold];

%---Run ART
for ii = 1:p.nses
      fprintf('=============Beginning Outlier Detection on %s Run %d=============\n', p.subj, ii);
      art(p,ii)
      %-Add Title
      F = gcf;
      titleax = axes('Position',[0.12 0.73 0.8 0.2],'Parent',F,'Visible','off');
      set(get(titleax,'Title'),'String',sprintf('Art Outlier Detection | Subject: %s Run: %d',p.subj,ii),'FontSize',16,'FontWeight','Bold','Visible','on');
      %-Print figures
      eval(defaults.printstr);
      pause;
      close(F);
      %-Move outlier files to Functional dir
      [s,w] = system(spm8w_osbabel(sprintf('mv "%s" "%s"',[fullfile(p.subdir,p.subj),'/outliers_run',num2str(ii),'.txt'],p.func)));    
end
 
%---Convert multipage ps file to pdf 
ps2pdf('psfile','tmp_art.ps','pdffile',[p.subj,'_art.pdf'],...
       'gspapersize','a4','deletepsfile',1);
%---delete flak.
delete('tmp_*.*')

%---Final Thoughts
%--Touch it and return to root!
eval(spm8w_osbabel(sprintf('!touch "%s"',fullfile(p.subdir,p.subj,['artdetect_',p.logsuffix]))));
cd(p.root);
fprintf('==============================================================\n');
fprintf('A txt file of timeseries outliers was saved to outlier_run#.txt...\n');
fprintf('A pdf file of Outlier Detection plots was saved to %s_art.pdf...\n',p.subj);


%---ADDITIONAL FUNCTIONS

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

