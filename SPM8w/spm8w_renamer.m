function spm8w_renamer(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: May, 2012 - DDW
% ==============================================================================
% spm8w_renamer(list,type,prefix)
% 
% spm8w_reanmer is a simpl utility that would probably be better implemented
% in bash! Anyhow, sicne we're constantly renaming things for clarity 
% here's a tool that takes a tab delimited list of dir or filenames and 
% prefixes and renames the associated dir/files such that the prefix comes
% first. Simple simple. 
% 
% Type is either 'dir' or 'files'
% If 'dir' spm8w_renamer will rename dirs as they appear in the list adding
% the prefix in column two to the original name.
% If 'files' spm8w_renamer will use the 1st column as a wildcard to find files 
% beginning with that name. Each of those files will in turn have the prefix 
% in column 2 added to it's name. Errr... maybe i should make an example.
%
% If preifx is 0, then renamer simply appends the prefix in column 2 to the 
% dir/file name. This is the default behavior, such as when turning dir/files
% called 03feb12aa to control_03feb12aa (i.e. specifying group membership via 
% a prefix).
% 
% If prefix is set to 1, then the "prefix" now becomes the new name. So if the 
% prefix was control_aa, then the dir/files with 03feb12aa will become 
% control_aa and the 03feb12aa name will be lost forever! As with prefix = 0. 
% If the type is a dir, then the dir is simply renamed to whatever is in 
% column 2 (i.e., no append). For type = files, and prefix = 1, renamer will 
% instead remove the original wildcard and append the new prefix. I.e. 
% 03feb12aa_condition1.txt becomes aa_condition1.txt if the prefix in second 
% column as aa... Phew this is getting complicated. Rewrite this help section 
% sometime, yo.
%
% the text file should look like: 
% 03feb12aa ListA
% 03feb12bb ListB
% etc.
%
% if using a variable as input the var should be a cell, for instance:
% a = {'03feb12aa','ListA';'03feb12bb','ListB'};
%
% To call it: spm8w_renamer(a,'files')
% ==============================================================================
% CHANGE LOG:
% -Added support for prefix (i.e. append prefix, or rename to prefix).
% =======1=========2=========3=========4=========5=========6=========7=========8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch (nargin)
    case 2       
    listfile = varargin{1};
    type     = varargin{2};
    prefix    = 0;
    case 3   
    listfile = varargin{1};
    type     = varargin{2};
    prefix   = varargin{3};   
    otherwise
        error('Please specify a list variable or the filename of a list (i.e. list.txt) and a type of renaming');
end

%Check if file or var and populate the list./
if ~iscell(listfile)
    if exist(listfile,'file')
        list = textscan(fopen(listfile), '%s %s');  
        list = [list{1},list{2}]; %reformat the textscan output
    else
        error('The file %s does not appear to exist.',listfile);
    end
else
    list = listfile;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Renamer loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(type,'dir')
    for i = 1:size(list,1)
        if exist(list{i,1},'dir')
            if prefix == 0
                newname = [list{i,2},'_',list{i,1}];
            elseif prefix == 1
                newname = [list{i,2}];
            else
                error('Invalid prefix setting, please use 0 or 1');
            end
            fprintf('Renaming %s to %s...\n',list{i,1},newname);   
            [s,w] = system(spm8w_osbabel(sprintf('mv %s %s',list{i,1},newname)));    
        else
           fprintf('Directory %s does not exist, skipping...\n',list{i,1}); 
        end
    end
    fprintf('done\n');
elseif strcmp(type,'files')
    for i = 1:size(list,1)
        files = dir([list{i,1},'*']);
        for ii = 1:length(files)
            if prefix == 0
                newname = [list{i,2},'_',files(ii).name];
            elseif prefix == 1
                 newname = [list{i,2},files(ii).name(length(list{i,1})+1:end)];
            else
                error('Invalid prefix setting, please use 0 or 1');
            end
            fprintf('Renaming %s to %s...\n',files(ii).name,newname);   
            [s,w] = system(spm8w_osbabel(sprintf('mv %s %s',files(ii).name,newname)));       
        end
        clear files
    end  
    fprintf('done\n');  
else
    error('Type(%s) is incorrect. Only dir or files is allowed',type);
end








