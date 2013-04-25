function spm8w_rfx(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: March, 2006
% ==============================================================================
% spm8w_rfx(rfxdir)
%
% spm8w_rfx will perform a random effects analysis (one-sample t-test) across
% a set of contrast maps from a 1st level analysis (using spm8w_contrasts). 
%
% If spm8w_rfx is supplied with a directory, it will automagically determine if
% this is a top level RFX direcotry (containing multiple contrast directories)
% or an individual contrast directory. In the frormer case it will traverse 
% through all the subdirectories, in the latter it will perform an rfx analysis
% on only that directory.
% ==============================================================================
% CHANGE LOG:
% -Seperated from the additional functions section of spm8w.m so that it can
% be called seperately by other functions.
% =======1=========2=========3=========4=========5=========6=========7=========8

%---Input checks
%--Set some vars, clear SPM and make P global
cwd = pwd;
switch (nargin)
  case 0 
    rfx_dirs = spm_select(Inf, 'dir','Please select directories for RFX analysis',[],[cwd,'/RFX'],'^\w.*');   
  case 1
    rfx_dirs = varargin{1};
    %fix trailing slash needed later
    if rfx_dirs(end) ~= '/' && rfx_dirs(end) ~= '\'
        rfx_dirs(end+1) = '/';
    end    
    %-Check user input exists
    if ~exist(rfx_dirs, 'dir')
        error('The specified directory %s does not exist!', rfx_dir);
    end
  otherwise
    error('You''ve specified too many inputs. Specify only the path to a single directory.');
end

%---RFX analysis setup
D = struct('DesName','One sample t-test', ...
           'n',[Inf 1 1 1], 'sF',{{'obs','','',''}}, ...
           'Hform', 'I(:,2),''-'',''mean''', ...
           'Bform', '[]', ...
           'nC',[0,0],'iCC',{{[8],[8]}},'iCFI',{{[1],[1]}}, ...
           'iGXcalc',[-1],'iGMsca',[-9],'GM',[], ...
           'iGloNorm',[9],'iGC',[12], ...
           'M_',struct('T',[-Inf],'I',Inf,'X',0), ...
           'b',struct('aTime',0));
con_val=1; % Set up t contrast value

%---RFX analysis
%--Check if this is top level RFX directory, if so find contents and make a 
%--new list. - DDW Apr/10
rfx_dirs_tmp = {};
for i = 1:length(rfx_dirs(:,1))
    dirlist      = dir(deblank(rfx_dirs(i,:)));        
    dirlist(1:2) = [];                          % Clear first two entries (i.e. '.' and '..')   
    for j = 1:length(dirlist)
        if dirlist(j).isdir == 1                % If it's a dir, add it to the list
            rfx_dirs_tmp{end+1} = [deblank(rfx_dirs(i,:)),dirlist(j).name,'/'];
        elseif strfind(dirlist(j).name,'.hdr')
            rfx_dirs_tmp{end+1} = deblank(rfx_dirs(i,:));
            break                               % Get me out of this damn loop y'all.
        end
    end  
end 
rfx_dirs = rfx_dirs_tmp;                        % Transfer the magic back to rfx_dirs. 

start_time = datestr(now);
for j=1:length(rfx_dirs)
    fprintf('==========Performing RFX on %s at %s\n', spm_str_manip(rfx_dirs{j},'ht'),datestr(now)); 
    cd(rfx_dirs{j});  
    %-Delete previous RFX files 
    deletefiles = {'spm';'mask';'SPM';'ess';'ResMS';'RPV';'beta'};
    for i = 1:length(deletefiles)
        try
            delete([deletefiles{i},'*.*']);
        catch
        end
    end            
    %-Compute the model
    spm8w_do_rfx('./',D);
    %-Make a t contrast for this model
    fprintf('==========Running t contrast on group model in %s \n',spm_str_manip(rfx_dirs{j},'ht'));
    con_name = spm_str_manip(pwd,'t');
    load SPM;
    %-Manually set xCon since SPM8 no longer does this and 
    %-we get a structure error if its empty.
    Fcname          = 'effects of interest';
    iX0             = [SPM.xX.iB SPM.xX.iG];
    SPM.xCon        = spm_FcUtil('Set',Fcname,'F','iX0',iX0,SPM.xX.xKXs);
    SPM.xCon(end+1) = spm_FcUtil('Set',con_name,'T','c',con_val,SPM.xX.xKXs);
    spm_contrasts(SPM);
    clear SPM;
    fprintf('\n\n');
    cd(cwd);
end
time_elapsed              = etime(datevec(datestr(now)),datevec(start_time));
[hours, minutes, seconds] = spm8w_timecalc(time_elapsed);
fprintf('\nRFX finished. Time to complete: %d hours, %d minutes and %d seconds ...\n\n',hours,minutes,seconds);