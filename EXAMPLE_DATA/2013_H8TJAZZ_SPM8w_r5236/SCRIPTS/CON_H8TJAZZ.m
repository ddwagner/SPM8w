% ==============================================================================
% SPM8w r5236
% Contrasts File
% Last update: February 2013 - DDW
% =======1=========2=========3=========4=========5=========6=========7=========8

%---LOAD SPM STRUCTURE AND FIND NUMBER OF NUISSANCE REGRESSORS
load SPM;
SPM.xCon=SPM.xCon(1);                     %Reset contrasts
len = sum(strcmp(SPM.Sess.C.name,'r'))+1; %Number of nuissance regressors

%---SPECIFY CONTRASTS
%--Only specify contrasts of interests, use zeros for nuisance regressors.
%--Use len variable to set zeroes in large designs (i.e., zeros(1,len))
%--Contrast weights should sum to zero (except all vs baseline) 
%--Contrast weights should never exceed one 
%--(i.e., [1 -1/2 -1/2] and not [2 -1 -1])
% ==============================================================================
% Contrast order: [Human, Animal, Vegetable, Mineral]

con_vals{1}         = [1 1 1 1 zeros(1,len)];
con_name{1}         = 'allVSbas';

con_vals{2}         = [1 0 0 0 zeros(1,len)];
con_name{2}         = 'humVSbas';

con_vals{3}         = [0 1 0 0 zeros(1,len)];
con_name{3}         = 'aniVSbas';

con_vals{4}         = [0 0 1 0 zeros(1,len)];
con_name{4}         = 'vegVSbas';

con_vals{5}         = [0 0 0 1 zeros(1,len)];
con_name{5}         = 'minVSbas';

con_vals{6}         = [1 -1/3 -1/3 -1/3 zeros(1,len)];
con_name{6}         = 'humVSavm';

con_vals{7}         = [1 -1 0 0 zeros(1,len)];
con_name{7}         = 'humVSani';

con_vals{8}         = [1 0 -1 0 zeros(1,len)];
con_name{8}         = 'humVSveg';

con_vals{9}         = [1 0 0 -1 zeros(1,len)];
con_name{9}         = 'humVSmin';

% ==============================================================================
% DO NOT EDIT BEYOND THIS POINT OR THERE WILL BE NO CHRISTMAS
% ==============================================================================

%---ADD CONTRASTS TO SPM STRUCTURE
for j = 1:length(con_vals)
    SPM.xCon(end + 1) = spm_FcUtil('Set', con_name{j}, 'T', 'c', ...
    (con_vals{j})', SPM.xX.xKXs);
end

%---EVALUATE CONTRASTS
spm_contrasts(SPM);

