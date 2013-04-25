function spm8w_onsmaker(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: February, 2010 - DDW
% ==============================================================================
% FORMAT spm8w_onsmaker('rtfile.csv)
%
% Takes a csv formatted file (see example onsets_h8tjazz.xls) with subjects, 
% conditions, parametrics, TRs and RTs and outputs condition onsets, par files, 
% and subject_condition_dur.txt files. Will also converts any NR to mean of 
% condition. To use you must format file accorting to onsets_h8tjazz.csv
% ==============================================================================
% CHANGE LOG:
% -Modified to make all onsets, not just duration files Feb 2010 - DDW
% -Switching to csv files because xlsread is giving too many 
%  errors. June 2010 - DDW
% =======1=========2=========3=========4=========5=========6=========7=========8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cwd = pwd;

switch (nargin)
  case 0 
    try 
        cd('ONSETS')
    end
    onsmatrix = spm_select(1,'^.*\.csv$','Please select your onsets csv file',[],[]);
    cd(['/',spm_str_manip(onsmatrix,'h')]);
  case 1
    onsmatrix = varargin{1};
  otherwise
    error('You''ve specified too many inputs. Only specify csv file.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('Loading file %s\n',onsmatrix);
% warning off
% [num,txt]=xlsread(onsmatrix);
% warning on
%Load csv file and assign to a num,txt format.
fid = fopen(onsmatrix, 'r');
    header = fgetl(fid);
    %parse header
    i = 1;
    while ~isempty(header)
        [t, header] = strtok(header,',');
        txt(i) = {t};
        i = i + 1;
    end
    %parse for txt
    col2 = textscan(fid, '%*s%s%*[^\n]','Delimiter',',');
    %pase for num
    frewind(fid);
    fgetl(fid); %skips header
    numtmp = textscan(fid, ['%n%*s',repmat('%n',1,length(txt)-2)],'Delimiter',',','TreatAsEmpty',{'NA','na','NR','nr'});
fclose(fid);
clear ans fid header i onsmatrix t 
%Reconstruct txt and num as it was with xlsread
%Put col2 into column 2 of txt var
txt(2:length(col2{1})+1,2) = col2{1}(1:end);
%Fix num so that it includes NANs in 2nd column
num      = numtmp{:,1};
num(:,2) = NaN(length(num),1);
for i = 2:length(numtmp)
    num(:,i+1) = numtmp{:,i};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign names and condition sizes to structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Put subject names in an array
fprintf('Submitting onset data to a rather invasive search');
%Assuming onsfile is formatted correctly, subjects name should start at column 6
onsdata.subjects = txt(1,6:length(txt(1,:)));
%Get parametric names and put into structure
%onsdata.paranames = txt(1,3:5); NOT NECESSARY
%Put conditions and condition length into an array (should be col2)
tmp = txt(2:length(txt),2);
condname(1) = tmp(1);
condlength(1) = 1;
current=1;
for i = 1:length(tmp)
    if i > 1
        if strcmp(condname{current},tmp{i})
            condlength(current)= condlength(current) + 1;
        else
            current=current+1;
            condname(current) = tmp(i);
            condlength(current) = 1;
        end
    end
    fprintf('.'); 
end
%dump to structure
onsdata.cond=condname;
onsdata.condsize=condlength;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign TRs and parametrics to appropriate conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currentTR=1;
for i = 1:length(onsdata.cond)
    onsdata.tr{i}=num(currentTR:currentTR+onsdata.condsize(i)-1,1);
    onsdata.para1{i}=num(currentTR:currentTR+onsdata.condsize(i)-1,3);
    onsdata.para2{i}=num(currentTR:currentTR+onsdata.condsize(i)-1,4);
    onsdata.para3{i}=num(currentTR:currentTR+onsdata.condsize(i)-1,5);
    currentTR=currentTR+onsdata.condsize(i);
    fprintf('.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correct RTs for NANs and assign to structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nReplacing missing RTs with within condition mean');
%Make a data structure
currentTR=1;
for i = 1:length(onsdata.cond)
    onsdata.rtdata{1,i}=num(currentTR:currentTR+onsdata.condsize(i)-1,6:length(onsdata.subjects)+5);    
    %calc mean without nans
    tmpmean=nanmean(onsdata.rtdata{i});
    %find out where is nan
    notnan=isnan(onsdata.rtdata{i});
    %replace with mean per subject
    tmpdata=onsdata.rtdata{i};
    for ii = 1:length(onsdata.subjects)
        tmpdata(notnan(:,ii),ii)=tmpmean(ii);
        fprintf('.')
    end
    onsdata.rtdata{i}= tmpdata; 
    currentTR=currentTR+onsdata.condsize(i);
    fprintf('.')
end
fprintf('\nNo responses destroyed and replaced with within condition mean\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spit out the onsets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Writing out files\n');
%Condition onsets
for i = 1:length(onsdata.cond)
    fname=[onsdata.cond{i},'.txt'];
    tmpdata=onsdata.tr{i};
    save(fname,'tmpdata','-ASCII');
    fprintf('%s written\n',fname);
end
%Parametric onsets MAX 3 para files per condition is hardcoded, can change if needed
for i = 1:3
    for ii = 1:length(onsdata.cond)
        %check first row for NAN, if nan, assume all rows are nan and skip
        if ~isnan(eval(['onsdata.para',num2str(i),'{ii}(1)']))
            fname=[onsdata.cond{ii},'_par',num2str(i),'.txt'];
            tmpdata=eval(['onsdata.para',num2str(i),'{ii}']);
            save(fname,'tmpdata','-ASCII');   
            fprintf('%s written\n',fname);
        end
    end
end
%Duration onsets.
for i = 1:length(onsdata.subjects)
    for ii = 1:length(onsdata.cond)
        fname=[onsdata.subjects{i},'_',onsdata.cond{ii},'_dur.txt'];
        tmpdata=onsdata.rtdata{ii};
        tmpdata=tmpdata(:,i);
        save(fname,'tmpdata','-ASCII');
        fprintf('%s written\n',fname);
    end
end
cd(cwd);



