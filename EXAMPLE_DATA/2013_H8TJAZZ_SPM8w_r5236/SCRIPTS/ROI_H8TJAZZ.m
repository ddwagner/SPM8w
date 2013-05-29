% ==============================================================================
% SPM8w r5236
% ROI Parameters File
% Last update: March, 2013 - DDW
% =======1=========2=========3=========4=========5=========6=========7=========8

%---DIRECTORY AND FILE SPECIFICATIONS
study_dir     = '/flash/KELLERTON/DDW/2013_H8TJAZZ_SPM8w_r5236';                     
roi_dir       = 'H8TJAZZ';           %dir where ROI results will be written                           
rfx_dir       = 'RFX_H8TJAZZ';       %dir containing con files for ROI 
roi_img_dir   = 'ROI_images';        %dir containing anatomical or functional masks <optional>
roi_spec_file = 'h8tjazz_spec.xlsx'; %File containing the roi specifications
roi_var_file  = 'h8tjazz.xlsx';      %Name of excel file with variables for ttest and correl
data_name     = 'H8TJAZZ.txt';       %Name of data file
    
%---ROI SPECFICIATIONS
%--Conditions must be the same name as condition directories in RFX dir
r.conditions = { 
		'humVSbas'
		'aniVSbas'
		'vegVSbas'
		'minVSbas'
				};

r.roi_specs = { 
        %If using a spec file this variable will be ignored.
		%To use mask write filename ('AAL_LAMYG.nii') instead of sphere size
		%Do not use '/' in REGION name or you'll irk the matlab                                  
		%REGION------------------COORDINATES----------SIZE|FILE
		'L.DMPFC.BA10'           '-3,63,18'              '8'
		'R.PRECUNEUS.BA31'       '3,-60,24'              '8'
		'L.AMYGDALA'             ''                'AAL_LAMYG.nii'
              };
    
r.roi_stats = {                      
		%Performs one or two sample t-tests, correlations and descriptives
        %'all_conditions is a reserved word. Otherwise type the
        %condition name or a formula.
        %valid stats: descriptives, t-test1,t-test2,correl,correl2
  	    %STATISTIC-------CONDITION/FORMULA
    	'descriptives'   'all_conditions' 
	    't-test1'        'humVSbas-(aniVSbas+vegVSbas+minVSbas)/3'
	    't-test1'        'humVSbas'
	    'correl'         'vegVSbas-minVSbas'
        'correl2'        'humVSbas'
   		't-test2'        'humVSbas-aniVSbas'
    		  };
                             
% ==============================================================================
% DO NOT EDIT BEYOND THIS POINT OR YOUR FEW REMAINING FRIENDS WILL ABANDON YOU 
% ==============================================================================
%---DIRS
roi_var_dir  = 'VARIABLES'; %dir containing vars for correlations and t-tests
roi_spec_dir = 'SPEC';      %dir containing the roi excel spec files
roi_rfx_dir  = 'RFX';       %dir where RFX results are stored, relative to root

%---PATHS
%r.root: The root dir for your study
%r.roi: The dir where ROI data will be written
%r.roi_img: The dir containing optional ROI mask images
%r.var_dir: The dir containing the vars files for corr/ttest2
%r.rfx_dir: The dir containing condition directories for parameter extraction
r.root       = study_dir;                           
r.roi        = fullfile(r.root,'ROI',roi_dir);     
r.roi_img    = fullfile(r.root,'ROI',roi_img_dir); 
r.rfx_dir    = fullfile(r.root,roi_rfx_dir,rfx_dir);           
r.spec_file  = fullfile(r.root,'ROI',roi_spec_dir,roi_spec_file); 
r.var_file   = fullfile(r.root,'ROI',roi_var_dir,roi_var_file); 

%--ASSIGN DATA NAME TO R STRUCTURE
r.data_name = data_name;         
                                
%--LOCATION OF STANDARD IMAGE WITH VOL SPECS (FOR ROI GENERATION)
r.standard_space = which('standard_space.nii');