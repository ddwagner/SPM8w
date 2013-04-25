% ==============================================================================
% SPM8w r5236
% VOI Parameters File
% Last update: February 2013 - DDW
% =======1=========2=========3=========4=========5=========6=========7=========8

%---DIRECTORY AND FILE SPECIFICATIONS
study_dir   = '/flash/KELLERTON/DDW/2013_H8TJAZZ_SPM8w_r5236';                     
ons_dir     = 'PPI';         %dir where timeseries/PPI regressors will be written                 
glm_dir     = 'GLM_H8TJAZZ'; %GLM dir from which the VOI/PPI is defined
roi_img_dir = 'ROI_images';  %dir containing anatomical or functional masks <optional>

%---VOI SPECFICIATIONS
v.con_name    = 'humVSani';     %VOI defining contrast
v.con_pval    = 1;              %p value for contrast (1, 0.05, 0.001 etc.)
v.con_ext     = 0;              %extent trehshold
v.con_locmax  = 0;              %Search for nearest local maxima 0=NO|1=YES
                                %if there is none, then no VOI will be
                                %created for that subject.
v.con_maxdist = 15;             %Maximum distance from source coordinates (in mm)                          
v.eigenvar    = 1;              %output eigenvariate otherwise mean 1=YES|0=MEAN                             
v.voi_specs = {                 
		%VOI definitions. 'all_subjects' is special word
		%To use mask write filename ('AAL_LAMYG.nii') instead of sphere size
		%Do not use '/' in VOI name or you'll irk the matlab        
   		%SUBJECT------------VOINAME------COORDINATES------SIZE|FILE-----
	   'all_subjects'      'L.DMPFC'      '-3,63,18'         '8'
   			  };

%---PPI SPECFICIATIONS
v.ppi_graphics = 0;   %Will display and save the PPI deconvolution figure to a PDF
v.ppi_specs    = {
		%VOINAME   = Name of VOI to use
        %PPINAME   = Name to give PPI ('_PPI' is auto appended)
        %PSYVECTOR = Psy defining contrast
	    %VOINAME-----------PPINAME-----------PSYVECTOR-----
 	    %                              	   Hum Ani Veg Min
        'L.DMPFC'       'DMPFC_HvsA'         [1 -1 0 0]           
     			 };
        
%---PPI PLOT SPECFICATIONS
v.plot_export  = 'PNG';  %PDF|PS|PNG|EPS|JPG|TIFF or none for no export
v.plot_specs   = {
		%REGION1   = Name of VOI file (seed)
        %REGION2   = Name of VOI file (target)
        %CONTRAST1 = Pos con from PSYVECTOR 
        %CONTRAST2 = Neg con from PSYVECTOR (absolute value)
 	    %REGION1----------REGION2-----CONTRAST1--CONTRAST2----CON1NAME-------CON2NAME-----
   		'DMPFC'           'FFA'       [1 0 0 0]  [0 -1 0 0] '  Human'        'Animal'
    			 };
       
% ==============================================================================
% DO NOT EDIT BEYOND THIS POINT OR THEY WILL CANCEL TELEVISION INCLUDING HBO
% ==============================================================================
%---PATHS
%v.root: The root dir for your study
%v.ons_dir: The dir where VOI data will be written 
%v.roi_img: The dir containing optional ROI mask images
%v.glm_dir: GLM dir from which the VOI/PPI is defined
%v.ppi_dir: Where VOI/PPI regressors will be written/calculated
v.root      = study_dir;                           
v.ons_dir   = fullfile(v.root,'ONSETS',ons_dir);  
v.roi_img   = fullfile(v.root,'ROI',roi_img_dir);
%%-The glm_dir has to be a string since we need to eval it against
%-each row of v.voi_specs, which we don't do until spm8w_voitool.m 
v.glm_dir   = ['fullfile(''',v.root,''',''SUBJECTS'',v.voi_specs{i,1},''RESULTS',filesep,glm_dir,''')'];
v.ppi_dir   = ['fullfile(''',v.root,''',''SUBJECTS'',subjects{i},''RESULTS',filesep,glm_dir,''')'];