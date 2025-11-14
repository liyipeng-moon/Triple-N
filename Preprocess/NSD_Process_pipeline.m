clear
root_dir = 'C:\Users\moonl\Desktop\NNN';
raw_data_dir = fullfile(root_dir,'Data\Raw\SesFolder');
cd(root_dir)
addpath(genpath('C:\Users\admin\AppData\Roaming\MathWorks\MATLAB Add-Ons\Apps\NIMHMonkeyLogic22'))
addpath(genpath(root_dir))
%%
interested_path = {};
raw_data_dir = 'F:\NSD_FULL_RAW';
cd(raw_data_dir)
all_dir = dir('2*');
for dir_idx = 1:length(all_dir)
    interested_path{end+1} = fullfile(raw_data_dir,all_dir(dir_idx).name);
end
%%
for path_now = [1:90]
    Load_Data_function(interested_path{path_now});
    PostProcess_function_raw(interested_path{path_now});
    PostProcess_function(interested_path{path_now});
    PostProcess_function_LFP(interested_path{path_now});
end