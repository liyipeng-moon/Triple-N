% %
% Description:
% This script processes session-wise  data and extracts neural responses for manually defined brain areas. 
% For each session:
%   1. Load preprocessed neural data (Processed_sesXX*)
%   2. Load corresponding raw/meta data
%   3. Select units based on:
%        - reliability threshold (default: 0.4)
%        - spatial position constraints (manual annotation)
%   4. Compute:
%        - raw response (area_rsp)
%        - z-scored response across images (area_norm_rsp)
%        - population average response (area_avg_rsp)
%   5. Rank images based on average response (top 1000 stimuli)
%   6. Visualize top-4 preferred images for each area
%   7. Save:
%        - visualization (.png)
%        - data (.mat)
%%
clear;clc;close all
load DIRS.mat
cd(root_dir);
addpath(genpath(pwd));

reliability_thres = 0.4;
load img_pool.mat
face_idx = 1001:1024; body_idx = 1000+[26:31,43:48,50:61]; obj_idx = setdiff(1025:1072, body_idx);
mkdir area_wise_data
%% Gen probe-wise data
manual_data = readtable("exclude_area.xls");
interested_ses = 1:90;
for ses_now = 1:90
    proc1_file_name = dir(fullfile(prep_dir,sprintf('Processed_ses%02d*', interested_ses(ses_now))));
    proc1_file_name = proc1_file_name.name;
    pro1_data = load(fullfile(prep_dir,proc1_file_name));
    filename_here = dir(fullfile(H5_dir,sprintf('ses%02d*h5',interested_ses(ses_now))));
    filename_here = filename_here.name;
    metaname_here = dir(fullfile(H5_dir,sprintf('ses%02d*mat',interested_ses(ses_now))));
    metaname_here = metaname_here.name;
    meta_data = load(fullfile(H5_dir,metaname_here));
    area = find(manual_data.SesIdx==interested_ses(ses_now));
    for area_idx = 1:length(area)
        x1 = manual_data.y1(area(area_idx));
        x2 = manual_data.y2(area(area_idx));
        interested_unit = find(pro1_data.reliability_best>0.4 & pro1_data.pos > x1 & pro1_data.pos < x2);
        area_rsp = pro1_data.response_best(interested_unit,:);
        area_norm_rsp = zscore(area_rsp,0,2);
        area_avg_rsp = mean(area_norm_rsp);
        [a,b] = sort(area_avg_rsp(1:1000),'descend');

        rr = pro1_data.reliability_best(interested_unit);
        figure; img=[];
        for i = 1:4
            img =[img, imresize(img_pool{b(i)},[224,224])];
        end
        imshow(img)
        file_name = sprintf('%s_ses%02d_Idx%03d.mat', manual_data.AREALABEL{area(area_idx)},ses_now,area(area_idx));
        sgtitle(file_name(1:end-4),Interpreter='none')
        saveas(gcf,fullfile("area_wise_data",sprintf('%s.png',file_name(1:end-4))))
        
        area_rsp=single(area_rsp);area_norm_rsp=single(area_norm_rsp);area_avg_rsp=single(area_avg_rsp);rr=single(rr);
        save(fullfile('area_wise_data/', file_name),'area_rsp','area_norm_rsp','area_avg_rsp','rr');
    end
end