function spm8w_contrasts(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: March, 2006
% ==============================================================================
% spm8w_contrasts('sub_id','parameter_file', 'con_file')
%
% spm8w_contrasts calculates contrasts between weighted parameter estimates
% defined in the CON parameters file. 
%
% First argument is subject id, second can be full path to the params file, 
% or else a window will prompt you for it. The third argument is optional if and
% only if the contrasts file is defined in the paramters file (i.e., p.confile)
% ==============================================================================
% CHANGE LOG:
% -Seperated from the additional functions section of spm8w.m so that it can
% be called seperately by other functions.
% =======1=========2=========3=========4=========5=========6=========7=========8

%---Input checks
%--Set some vars, clear SPM and make P global
switch (nargin)
  case 0 
    error('Please specify a subject id, e.g. spm8w_compute(''01jan09ab'')');
  case 1
    p.subj = varargin{1};
    p      = spm8w_getp('P',p.subj);
  case 2
    p.subj    = varargin{1};
    para_file = varargin{2};
    p         = spm8w_getp('P',p.subj,para_file);
  case 3
    p.subj    = varargin{1};
    para_file = varargin{2};
    p         = spm8w_getp('P',p.subj,para_file);    
    p.confile = varargin{3};
  otherwise
    error('You''ve specified too many inputs. Either 1 or 2 or 3 only please.');
end

%---Check that RFX directories exist.
%-Check for DIRs
cd(p.root)
if ~isdir(fullfile(p.root,'RFX'))
    mkdir(fullfile(p.root,'RFX'));
    mkdir(p.rfx);
elseif ~isdir([p.rfx])
    mkdir(p.rfx);
end

%---Run contrasts
%--Get GLM name from p.glm for display purposes
[glm_dir, glm_name] = fileparts(p.glm);
%--Set start time
con_start = datestr(now);
%--Goto GLM dir.
cd(p.glm);
%--Print info
fprintf('==========Running contrasts for %s on %s ...\n',glm_name, p.subj);
copyfile(p.confile, p.glm);
%--Chop off path and extension and run the contrast function. 
eval(spm_str_manip(p.confile,'tr')); 
%--Copy con_*img files to RFX
%-Find where to start since at least con_1 is an F contrast
con_num     = find([SPM.xCon.STAT]=='T');
con_num     = con_num(1);
for y = 1:length(con_name)
    %-make RFX dir for each contrast
    con_dir = fullfile(p.rfx,con_name{y});
    if ~exist(con_dir,'dir')
        mkdir(con_dir);
    end
    %-Make sure we get the file names right
    %-Modified April 2009/DDW to allow for 
    %-more than 100 Contrasts (e.g. FIR models) 
    if con_num < 10
        suffix(y,:) = ['00',int2str(con_num)];
    elseif con_num < 100
        suffix(y,:) = ['0',int2str(con_num)];
    else    
        suffix(y,:) = int2str(con_num);
    end
    %-Copy the con*img and hdr files to the
    %-relevant RFX directory
    fprintf('\nCopying contrast %s to %s',num2str(con_num),con_dir);
    copyfile(['con_0',suffix(y,:),'.hdr'],fullfile(con_dir,['con_',p.subj,'_0',suffix(y,:),'.hdr']));
    copyfile(['con_0',suffix(y,:),'.img'],fullfile(con_dir,['con_',p.subj,'_0',suffix(y,:),'.img']));
    %-Update contrast number
    con_num=con_num+1;   
end % y contrast loop
fprintf('\n');

%---Touch it and return to root!
eval(spm8w_osbabel(sprintf('!touch "%s"',fullfile(p.subdir,p.subj,['contrasts_',p.logsuffix]))));
cd(p.root);

%---Spit Time
con_stop                  = datestr(now);
con_time_elapsed          = etime(datevec(con_stop),datevec(con_start)); %time elapsed in seconds
[hours, minutes, seconds] = spm8w_timecalc(con_time_elapsed);
fprintf('Subject %s contrasts computation finished at %s\nand took %d hours, %d minutes and %d seconds...\n',p.subj,con_stop,hours,minutes,seconds);
