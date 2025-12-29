function [LVR_m,LVR_h,TimeLag] = gen_species_difference(v_idx, l_idx)
root_dir = 'C:\Users\moonl\Desktop\NNN';
[proc_dir,~] = gen_dirs(root_dir);
manual_data = readtable("exclude_area.xls");
%% Analysis for alexnet_fc6 and mpnet
load('fMRI_result.mat')
Interested_Voxel = find(voxel_info_1~=0);
for vertex = 1:length(Interested_Voxel)
    fMRI_performance(vertex,1) = max(r(Interested_Voxel(vertex),v_idx));
    fMRI_performance(vertex,2) = r(Interested_Voxel(vertex),l_idx);
end
%%
AVG_data = [];
TimeCourseData = [];
TimeLag = [];
for ses = 1:102
    if(strcmp(manual_data.Area{ses},'IT'))
        EncodingData = dir(fullfile(proc_dir,sprintf("s5_encode_ses_%03d.mat",ses)));
        load(fullfile(proc_dir,EncodingData.name));
        Visual_Performance = max(pred_r_array(:,v_idx),[],2);
        LLM_Performance = pred_r_array(:,l_idx);
        AVG_data = [AVG_data; [Visual_Performance,LLM_Performance]];
        [~,v_time]=max(squeeze(max(pred_r2_array_t(:,v_idx,:),[],2)),[],2);
        [~,l_time]=max(squeeze(max(pred_r2_array_t(:,l_idx,:),[],2)),[],2);
        diff_t = l_time-v_time;
        good_unit = find(min(Visual_Performance,LLM_Performance)>0.05);
        TimeLag = [TimeLag, diff_t(good_unit)'];
    end
end
%%

x_here = AVG_data(:,1);
y_here = AVG_data(:,2);
[LVR_m,~] = demingRegression(x_here, y_here);

x_here = fMRI_performance(:,1);
y_here = fMRI_performance(:,2);
[LVR_h,~] = demingRegression(x_here, y_here);

TimeLag = mean(TimeLag);