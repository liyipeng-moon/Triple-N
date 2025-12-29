clear;clc
root_dir = 'C:\Users\moonl\Desktop\NNN';
cd(root_dir)
addpath(genpath(pwd));
[proc_dir,raw_dir] = gen_dirs(root_dir);
reliability_thres = 0.4;
all_neuron_data = [];
all_neuron_data_t1 = [];
all_neuron_data_t2 = [];
all_neuron_data_t3 = [];
all_neuron_data_t4 = [];
all_neuron_data_t5 = [];
info_data = load(dir(fullfile(raw_dir,'*ses01*info*')).name);

IT_session = [1:70, 88];
for ses_idx = IT_session
    
    proc1_file_name = dir(fullfile(proc_dir,sprintf('Processed_ses%02d*', ses_idx)));
    proc1_file_name = proc1_file_name.name;
    pro1_data = load(fullfile(proc_dir,proc1_file_name));
    rr_here = pro1_data.reliability_best;
    rsp_here = pro1_data.response_best(rr_here>0.4, 1:1000);
    rsp_here = zscore(rsp_here,0,2);
    all_neuron_data = [all_neuron_data; rsp_here];


    filename_here = dir(fullfile(raw_dir,sprintf('ses%02d*h5',ses_idx)));
    filename_here = filename_here.name;
    PSTHData = h5read(fullfile(raw_dir,filename_here), '/response_matrix_img');

    good_units = find(rr_here>0.4);
    ses_data_t3 = [];ses_data_t4 = [];ses_data_t5 = [];
    for uu = 1:length(good_units)
        uu_here = good_units(uu);
        t1 = pro1_data.best_r_time1(uu_here);
        t2 = pro1_data.best_r_time2(uu_here);
        r3 = squeeze(sum(PSTHData(uu_here,1:1000,info_data.global_params.pre_onset+(t1:t2)),3));
        ses_data_t3 = [ses_data_t3;zscore(r3)];
    end
    all_neuron_data_t3 = [all_neuron_data_t3;ses_data_t3];
    fprintf('%d \n', ses_idx)
end

ROI_info=load('ROI_info.mat');
subject_pool = {1,2,5,7};
for SS = 1:4
    interested_subject=subject_pool{SS};
    hemi_here = 'lh';
    fMRI_data = load(fullfile(root_dir,"NNN_Data/FMRI",sprintf('S%d_%s_Rsp.mat',interested_subject,hemi_here)));
    fMRI_data = fMRI_data.mean_brain_data;
    fMRI_data = double(fMRI_data)./300;
    lh_data = fMRI_data(find(getfield(ROI_info, sprintf('S%d_%s_General',interested_subject,hemi_here))),:);
    hemi_here = 'rh';
    fMRI_data = load(fullfile(root_dir,'NNN_Data/FMRI/',sprintf('S%d_%s_Rsp.mat',interested_subject,hemi_here)));
    fMRI_data = fMRI_data.mean_brain_data;
    fMRI_data = double(fMRI_data)./300;
    rh_data = fMRI_data(find(getfield(ROI_info, sprintf('S%d_%s_General',interested_subject,hemi_here))),:);
    all_data_subj{SS}=zscore([lh_data;rh_data],0,2);
end
all_brain_data = [all_data_subj{1}', all_data_subj{2}',all_data_subj{3}',all_data_subj{4}'];

load alexnet_resp.mat
LLM = load('LLM_allv2.mat');
llm_ebd = zeros([1000, 768]);
loc = 0;
for img = 1:1000
    loc = find(LLM.image_id==LLM.image_id(loc+1));
    llm_ebd(img,:) = mean(LLM.embeddings(loc,:),1);
    loc = loc(end);
end
max_repetition = 50000;
num_series = [2:5:99,100:20:1000];
fprintf('Load Data Finished!!!')
%% Decoding for different space
y_pool = {fc6_4096', llm_ebd};
ytitle_pool = {'Visual','Language'};
x_pool = {all_neuron_data_t3',all_brain_data};
xtitle_pool = {'Macaque','Human'};
ldg_pool = {};
acc_save = {};
for brain_here = 1:length(x_pool)
    [coeff,neuron_score,latent,tsquared,explained,mu] = pca(x_pool{brain_here},NumComponents=500);
    neuron_score = neuron_score';
    for space_here = 1:length(y_pool)
        [coeff,FC_score,latent,tsquared,explained,mu] = pca(y_pool{space_here},NumComponents=100);
        [acc_here,r_here] = fn_decoding_acc(neuron_score, FC_score,max_repetition,num_series,1:1000);
        ldg_pool{end+1} = sprintf('%s %s', xtitle_pool{brain_here}, ytitle_pool{space_here});
        acc_save{brain_here,space_here}=acc_here(num_series);
        r_save{brain_here, space_here} = r_here;
    end
end
ClusInfo = load('ClusInfo.mat'); 
clus_size = max(ClusInfo.Cluster_idx);
space_here = 1;[coeff,FC_score,latent,tsquared,explained,mu] = pca(y_pool{space_here},NumComponents=100);
acc_pool = [];
x_pool = {all_neuron_data_t3',all_brain_data}; xtitle_pool = {'Macaque','Human'};
for brain_here = 1:length(x_pool)
    [coeff,neuron_score,latent,tsquared,explained,mu] = pca(x_pool{brain_here},NumComponents=500);
    neuron_score = neuron_score';
    for cc = 1:clus_size
        [acc_here,r_here] = fn_decoding_acc(neuron_score, FC_score,max_repetition,12,find(ClusInfo.Cluster_idx==cc));
        acc_pool(brain_here,cc) = acc_here(12);
    end
    acc_12(brain_here) = fn_decoding_acc_cross(neuron_score, FC_score,max_repetition,ClusInfo.Cluster_idx);
end
save(fullfile(proc_dir,'S6.mat'))
%% Decoding Across Time

clear;clc
root_dir = 'C:\Users\moonl\Desktop\NNN';
[proc_dir,raw_dir] = gen_dirs(root_dir);
cd(root_dir)
addpath(genpath(pwd));
reliability_thres = 0.4;


IT_session = [1:70, 88];

decode_tin = -10:1:360;
allt_binned_data = zeros([30000, length(decode_tin),1000],'single');
unit_idx = 1;
info_data = load(dir(fullfile(raw_dir,'*ses01*info*')).name);
for ses_idx = IT_session
    proc1_file_name = dir(fullfile(proc_dir,sprintf('Processed_ses%02d*', ses_idx)));
    proc1_file_name = proc1_file_name.name;
    pro1_data = load(fullfile(proc_dir,proc1_file_name));
    rr_here = pro1_data.reliability_best;
    rsp_here = pro1_data.response_best(rr_here>0.4, 1:1000);
    rsp_here = zscore(rsp_here,0,2);
    filename_here = dir(fullfile(raw_dir,sprintf('ses%02d*h5',ses_idx)));
    filename_here = filename_here.name;
    PSTHData = h5read(fullfile(raw_dir,filename_here), '/response_matrix_img');
    good_units = find(rr_here>0.4);
    ses_bin_data = [];
    for uu = 1:length(good_units)
        uu_here = good_units(uu);
        tt_data_unit = zeros([length(decode_tin),1000]);
        for tpoint = 1:length(decode_tin)
            tt_here = decode_tin(tpoint);
            tt_data_unit(tpoint,:) = zscore(PSTHData(uu_here,1:1000,info_data.global_params.pre_onset+tt_here));
        end
        allt_binned_data(unit_idx,:,:) = single(tt_data_unit);
        unit_idx = unit_idx+1;
    end
    fprintf('%d %d\n', ses_idx,unit_idx)
    clear PSTHData
end

allt_binned_data(unit_idx:end,:,:)=[];
load alexnet_resp.mat
LLM = load('LLM_allv2.mat');
llm_ebd = zeros([1000, 768]);
loc = 0;
for img = 1:1000
    loc = find(LLM.image_id==LLM.image_id(loc+1));
    llm_ebd(img,:) = mean(LLM.embeddings(loc,:),1);
    loc = loc(end);
end
max_repetition = 50000;
num_series = [10,50,100,200];
fprintf('Load Data Finished!!!')
y_pool = {fc6_4096', llm_ebd};
ytitle_pool = {'Visual','Language'};
for tt = 1:size(allt_binned_data,2)
    x_pool{tt} = squeeze(allt_binned_data(:,tt,:))';
end
clear allt_binned_data
%% Decoding Across Time

tic
acc_here = {};
r_here = {};
start_parfor
for space_here = 1:length(y_pool)
    tic
    [coeff,FC_score,latent,tsquared,explained,mu] = pca(y_pool{space_here},NumComponents=100);
    ppm = ParforProgMon('Par...', length(x_pool) , 1, 1200, 200);
    for brain_here = 1:length(x_pool)
        [coeff,neuron_score,latent,tsquared,explained,mu] = pca(x_pool{brain_here},NumComponents=500);
        neuron_score = neuron_score';
        [acc_now, r_now] = fn_decoding_acc(neuron_score, FC_score,max_repetition,num_series,1:1000);
        acc_here{brain_here,space_here} = acc_now;
        r_here{brain_here,space_here} =r_now;
        ppm.increment();
    end
    space_here
end
toc
save(fullfile(proc_dir,'S6_time_decoding.mat'), "acc_here","r_here","decode_tin","num_series")

return