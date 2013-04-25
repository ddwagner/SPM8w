function spm8w_conjunction(varargin) 
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: February 2009 - JM/DDW
% ==============================================================================
% function spm8w_conjunction('output.img', NumConditions)
% Got this from Joe Moran. A script to create a map of the minimum 
% t-statistic across conditions (resulting in a logical AND conjunction)
%
% spm8w_conjunction('ThisIsMyFilenameYo',4)
% If no input arguments are given, spm8w_conjunction will ask for output 
% filename and numnber of conditions.
% ==============================================================================
% CHANGE LOG:
% - Added support for the conjunction of more than 2 conditions - DDW Feb/09
% - SPM8erized it! -DDW Jan/10
% - Changed the spm_defaults to spm('Defaults','fmri') - DDW Aug/10
% =======1=========2=========3=========4=========5=========6=========7=========8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check inputs
if nargin == 0
    fname   = spm_input('Name of output file?','+1','s');
    con_num = spm_input('Number of t-files?','+1','i');
elseif nargin == 1
    fname   = varargin{1};
    con_num = spm_input('Number of t-files?','+1','i');
elseif nargin == 2
    fname   = varargin{1};
    con_num = varargin{2};
else   
    msg1    = 'You''ve specified too many inputs. Simply enter and output filename and number of conjunctions to run';
    error(msg1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load t-stat files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%turn on spm_defaults to solve flipping messages
spm('Defaults','fmri')
%get pwd and try moving to rfx
cwd = pwd;
try 
    cd('RFX')
end

%Loop through condition number and get all t-files
vt = spm_vol(spm_select(con_num,'image',['Please select all ',num2str(con_num),' spm_t files'],[],[],'^spmT.*'));
%for i = 1:con_num
%   vt{i}=spm_vol(spm_select(1,'image',['Select all',num2str(i),'spm_t files for', num2str(i)]));
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop through t-stat files and create minimums
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make love to the command window
fprintf('==========Creating conjunction map from all %d conditions\n',con_num);
%Make an array of just the t_values.
for i = 1:con_num
   t{i} = spm_read_vols(vt(i));
end

%Calculate the minimum t values for all conditions
for i = 1:(con_num-1)
   if i==1
       Yout = max(0,min(t{1},t{2})) + min(0,max(t{1},t{2}));
   else
       Yout = max(0,min(Yout,t{i+1})) + min(0,max(Yout,t{i+1}));
   end
end

%Command windows came back for some more.
fprintf('Done... Writing out conjunction volume...\n');

%Go home
cd(cwd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Write out conjunction volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setup writing
Vout = vt(1);

%check for .img extension in case of twit user
[tmp,tmp,ext] = fileparts(fname);
if isempty(ext)
    fname = [fname,'.img'];
end
Vout.fname = fname;
spm_write_vol(Vout,Yout);
fprintf('Conjunction of all %d t-stat files saved to %s\n',con_num,fname);
close all;