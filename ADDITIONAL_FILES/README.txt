SPM8w
=====
The files in this directory should be added to your locally installed copy of SPM8.

canonical.7z
------------
The archive canonical.7z contains high quality brains in MNI space that may be used
for visualization. These files need to be placed inside the SPM8/canonical directory.

spm8
----
Dead simple bash script to start matlab in terminal mod and run the spm8 path script.
Only useful for Linux and possibly Mac.

spm8path.m
----------
Example script of how to setup your matlab paths for SPM8 and SPM8w. Although you can
manually add SPM8/SPM8w and its associated directories to your Matlab path via the gui,
this method is preferable as it allows multiple versions of SPM8 and multiple branches
of SPM8w to coexist (simply change the paths in this script to point to other versions).

spm_fMRI_design.m
-----------------
This is a modified version of the spm_fMRI_design.m file that comes with SPM8. You
need to overwrite the spm8 file with this one. The only modification is to add a 
catch to allow disabling of within-condition orthogonolization of regressors. 
Whether you do so or not depends on the variable disable_orth in your Parameters 
file. The SPM default is to enable within-condition orthogonalization, however
under certain circumnstances this may be undesirable. 
