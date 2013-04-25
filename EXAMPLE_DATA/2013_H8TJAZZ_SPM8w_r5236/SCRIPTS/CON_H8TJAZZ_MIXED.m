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
% For mixed block/er designs you can use some matlab tricks to make contrasts
% Conditions:     [-1to8---9to16---17to24---25to36-----37--------38----]
% Contrast order: [-human, animal, oncue,   offcue, implicit, explicit-]
effects = 34; %Number of columns of interest (all events + states)

con_name{1}         = 'implicitState';
con_vals{1}         = [zeros(1,32) 1 0 zeros(1,len)];

con_name{2}         = 'explicitState';
con_vals{2}         = [zeros(1,32) 0 1 zeros(1,len)];
 
con_name{3}         = 'allStates';
con_vals{3}         = [zeros(1,32) 1 1 zeros(1,len)];
 
con_name{4}         = 'humanT1';
con_vals{4}         = [1 0 0 0 0 0 0 0 ... 
                      zeros(1,effects-8+len)];

con_name{5}         = 'humanT2';
con_vals{5}         = [0 1 0 0 0 0 0 0 ... 
                      zeros(1,effects-8+len)];
                                     
con_name{6}         = 'humanT3';
con_vals{6}         = [0 0 1 0 0 0 0 0 ... 
                      zeros(1,effects-8+len)];                  
                  
con_name{7}         = 'humanT4';
con_vals{7}         = [0 0 0 1 0 0 0 0 ... 
                      zeros(1,effects-8+len)];
  
con_name{8}         = 'humanT5';
con_vals{8}         = [0 0 0 0 1 0 0 0 ... 
                      zeros(1,effects-8+len)];                
     
con_name{9}         = 'humanT6';
con_vals{9}         = [0 0 0 0 0 1 0 0 ... 
                      zeros(1,effects-8+len)];             
                                  
con_name{10}         = 'humanT7';
con_vals{10}         = [0 0 0 0 0 0 1 0 ... 
                      zeros(1,effects-8+len)];                                   

con_name{11}         = 'humanT8';
con_vals{11}         = [0 0 0 0 0 0 0 1 ... 
                      zeros(1,effects-8+len)];

con_name{12}         = 'animalT1';
con_vals{12}         = [zeros(1,8) ...
                        1 0 0 0 0 0 0 0 ... 
                      zeros(1,effects-16+len)];

con_name{13}         = 'animalT2';
con_vals{13}         = [zeros(1,8) ...
                        0 1 0 0 0 0 0 0 ... 
                      zeros(1,effects-16+len)];
                                     
con_name{14}         = 'animalT3';
con_vals{14}         = [zeros(1,8) ...
                        0 0 1 0 0 0 0 0 ... 
                      zeros(1,effects-16+len)];                  
                  
con_name{15}         = 'animalT4';
con_vals{15}         = [zeros(1,8) ...
                        0 0 0 1 0 0 0 0 ... 
                      zeros(1,effects-16+len)];
  
con_name{16}         = 'animalT5';
con_vals{16}         = [zeros(1,8) ...
                        0 0 0 0 1 0 0 0 ... 
                      zeros(1,effects-16+len)];                
     
con_name{17}         = 'animalT6';
con_vals{17}         = [zeros(1,8) ...
                        0 0 0 0 0 1 0 0 ... 
                      zeros(1,effects-16+len)];             
                                  
con_name{18}         = 'animalT7';
con_vals{18}         = [zeros(1,8) ...
                        0 0 0 0 0 0 1 0 ... 
                      zeros(1,effects-16+len)];                                   

con_name{19}         = 'animalT8';
con_vals{19}         = [zeros(1,8) ...
                        0 0 0 0 0 0 0 1 ... 
                      zeros(1,effects-16+len)];

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

