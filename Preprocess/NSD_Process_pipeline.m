clear
root_dir = 'C:\Users\moonl\Desktop\Triple-N-main\Preprocess';
cd(root_dir);addpath(genpath(root_dir))

addpath(genpath('C:\Users\admin\AppData\Roaming\MathWorks\MATLAB Add-Ons\Apps\NIMHMonkeyLogic22'))
% install MonkeyLogic at https://monkeylogic.nimh.nih.gov/download.html

% this is where you extract raw session folders from 3 links.
raw_data_dir = fullfile(root_dir,'Data\Raw\SesFolder'); 
cd(raw_data_dir);
interested_path = {}; 
all_dir = dir('2*');
for dir_idx = 1:length(all_dir)
    interested_path{end+1} = fullfile(raw_data_dir,all_dir(dir_idx).name);
end
%%
for path_now = 1:90
    Load_Data_function(interested_path{path_now});
    s(interested_path{path_now});
    PostProcess_function(interested_path{path_now});
    PostProcess_function_LFP(interested_path{path_now});
end