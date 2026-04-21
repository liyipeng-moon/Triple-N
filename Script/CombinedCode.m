% modify this for your computer
root_dir = 'C:\Users\moonl\Desktop\NNN';
H5_dir = 'C:\Users\moonl\Desktop\NNN\Data\Raw\H5FILES';
prep_dir = 'C:\Users\moonl\Desktop\NNN\Data\processed';
code_dir = 'C:\Users\moonl\Desktop\Triple-N-main\Script';
% extract others/FMRI.zip and others/ModelFeature.zip and paste here!
Model_Dir = 'C:\Users\moonl\Desktop\NNN\Data\others\ModelFeature';
FMRI_DIR = 'C:\Users\moonl\Desktop\NNN\Data\FMRI';
cd(code_dir);
save('DIRS.mat',"root_dir","H5_dir","prep_dir","code_dir","Model_Dir","FMRI_DIR")
addpath(genpath(pwd))
CheckEnv
%% From raw to processed
% if you've download 'processed' folder from scienceDB, you can directly
% skip this part.
analysis_S1_NSD_NC
%% Figure1, EF2 related
plot_F0_BasicInfo
plot_F1_b
plot_F1_e
plot_F1_fgh_F2_f
%% EF2 related
plot_EF2_abe
% EF2cd is produced by plot_EF2_cd.ipynb
plot_EF2fghijkl
%% Figure2, EF3 related
plot_F2_EF3
%% Figure3, EF45 related - clustering psths
analysis_S3_psth 
% Note: S3_psth must be run beforew other Figure 3 analysis, 
% as the following code depends on the clustering results.
plot_F3_image_wise % this is figure 3 lower part
plot_F3S_pref % this extended figure 5
%% Figure4, EF6 related
% please uncompress the human fMRI data first.
% it't located in /others/FMRI.zip
analysis_and_plot_F4a_and_corrmap
% surface plot if reproduced by plot_F4_PY_plot_surface_data.ipynb
plot_F4ghij
plot_EF6
% To plot human surface, refer to plot_F4_PY_plot_surface_data.ipynb
%% Figure5, EF78 related
warning('Start unit-wise and voxel-wise encoding analysis! This could take hours...')
analysis_S5_encoding_forloop
analysis_S5_encoding_NSD
% please revise line 110 as 'decode_tin = -10:X:360 if memory if not enough
% in our test, X can be 1 ms with 48G RAM, and 5ms for 24 G RAM
analysis_S5_decoding
% note that above code (fMRI encoding, ephys encoding, and decoding) need 
% to run first before we plot figure 5 related data.
plot_F5_EF78
%% EF9 related
plot_EF9_EVC