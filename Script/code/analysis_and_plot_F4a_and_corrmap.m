clear;clc
root_dir = 'C:\Users\moonl\Desktop\NNN';
cd(root_dir)
addpath(genpath(pwd));
[proc_dir,raw_dir] = gen_dirs(root_dir);
reliability_thres = 0.4;
manual_data = readtable("exclude_area.xls");
stats_dir = fullfile(root_dir,"Figs/stats/");
mkdir Figs\
mkdir Figs\F4\
UnitType = {};
neuron_area = [];
bigROI = [];
bigROI_areawise=[];
interested_time_point = -49:5:350;
ephys_time_course = [];
meta_example = load("ses10_240726_M5_3_info.mat");
for area_now = 1:size(manual_data,1)
    ses_idx = manual_data.SesIdx(area_now);
    proc1_file_name = dir(fullfile(proc_dir,sprintf('Processed_ses%02d*', ses_idx)));proc1_file_name = proc1_file_name.name;
    pro1_data = load(fullfile(proc_dir,proc1_file_name));
    filename_here = dir(fullfile(raw_dir,sprintf('ses%02d*h5',ses_idx)));filename_here = filename_here.name;
    PSTHData = h5read(fullfile(raw_dir,filename_here), '/response_matrix_img');
    x1 = manual_data.y1(area_now);x2 = manual_data.y2(area_now);pos_now = pro1_data.pos;

    unit_here = find(pos_now>x1 & pos_now<x2 & pro1_data.reliability_best>reliability_thres);

    neuron_area = [neuron_area, area_now*ones([1, length(unit_here)])];
    bigROI = [bigROI, repmat(manual_data.Area{area_now}(1),[1,length(unit_here)])];
    bigROI_areawise(area_now) = manual_data.Area{area_now}(1);
    rsp_here{area_now} = pro1_data.response_best(unit_here,1:1000);
    data_here = PSTHData(unit_here,1:1000,meta_example.global_params.pre_onset+interested_time_point);
    ephys_time_course = [ephys_time_course; data_here];

    fprintf('Loaded Data From %d \n', area_now)
    clear PSTHData
end
%% Load Human NSD
ROI_info=load('ROI_info.mat');
subject_pool = {1,2,5,7};
fMRI_dataset_path = 'C:\Users\moonl\Desktop\NNN\NNN_Data\FMRI';
fMRI_data_subject_pool ={};
all_data_fmri_general=[];
all_data_fmri_EVC=[];
for SS = 1:4
    interested_subject=subject_pool{SS};
    hemi_here = 'lh';
    fMRI_data = load(fullfile(fMRI_dataset_path,sprintf('S%d_%s_Rsp.mat',interested_subject,hemi_here)));
    fMRI_data = fMRI_data.mean_brain_data; fMRI_data = double(fMRI_data)./300;
    lh_data = fMRI_data;
    hemi_here = 'rh';
    fMRI_data = load(fullfile(fMRI_dataset_path,sprintf('S%d_%s_Rsp.mat',interested_subject,hemi_here)));
    fMRI_data = fMRI_data.mean_brain_data; fMRI_data = double(fMRI_data)./300;
    rh_data = fMRI_data;

    LOC_data=[];
    % LH general
    ROI_map = find(getfield(ROI_info, sprintf('S%d_lh_General',interested_subject)));
    all_data_fmri_general = [all_data_fmri_general;lh_data(ROI_map,:)];
    LOC_data = [LOC_data; lh_data(ROI_map,:)];
    % RH general
    ROI_map = find(getfield(ROI_info, sprintf('S%d_rh_General',interested_subject)));
    all_data_fmri_general = [all_data_fmri_general;rh_data(ROI_map,:)];
    LOC_data = [LOC_data; rh_data(ROI_map,:)];

    EVC_data = [];
    %LH EVC
    ROI_map = find(getfield(ROI_info, sprintf('S%d_lh_EVC',interested_subject)));
    all_data_fmri_EVC = [all_data_fmri_EVC;lh_data(ROI_map,:)];
    EVC_data = [EVC_data; lh_data(ROI_map,:)];
    %RH EVC
    ROI_map = find(getfield(ROI_info, sprintf('S%d_rh_EVC',interested_subject)));
    all_data_fmri_EVC = [all_data_fmri_EVC;rh_data(ROI_map,:)];
    EVC_data = [EVC_data; rh_data(ROI_map,:)];

    % 1 is VTC, 2 is EVC
    all_data_subj{SS,1}=LOC_data;
    all_data_subj{SS,2}=EVC_data;
end

%% RSA part
for human_brain = 1:2
    for ss = 1:4
        all_data_fmri =  zscore(all_data_subj{ss,human_brain},0,2);
        rdm_human{human_brain,ss} = 1-squareform(pdist(all_data_fmri',"correlation"));
        vec_human{human_brain,ss} = gen_rdm_vec(rdm_human{human_brain,ss});
    end
    sim_within_human = [];
    for s1 = 1:4
        for s2 = 1:4
            if(s1~=s2)
                sim_within_human = [sim_within_human, corr(vec_human{human_brain,s1},vec_human{human_brain,s2},'type','Spearman')];
            end
        end
    end
    SIM_cross_human{human_brain} = sim_within_human;
end

for macaque_brain = 1:2
    switch macaque_brain
        case 1
            interested_area = 'I'; % IT cortex
        case 2
            interested_area = 'E'; % EVC
    end
    time_course_data = ephys_time_course(bigROI==interested_area,:,:);
    rsa_vec = [];
    for tt = 1:length(interested_time_point)
        data_here = squeeze(time_course_data(:, :, tt));
        data_here = zscore(data_here, 0, 2);
        rdm_this_time =  1-squareform(pdist(data_here',"correlation"));
        rsa_vec(tt,:) = gen_rdm_vec(rdm_this_time);
        fprintf('finished %d s \n', interested_time_point(tt))
    end
    vec_monkey_time{macaque_brain} = rsa_vec;
    all_data_ephys = [];
    for area_here = 1:max(neuron_area)
        if(bigROI_areawise(area_here)==interested_area)
            all_data_ephys = [all_data_ephys; rsp_here{area_here}];
        end
    end
    all_data_ephys = zscore(all_data_ephys, 0, 2);
    rdm_monkey{macaque_brain} =  1-squareform(pdist(all_data_ephys',"correlation"));
    vec_monkey{macaque_brain} = gen_rdm_vec(rdm_monkey{macaque_brain});
end

for human_brain = 1:2
    for monkey_brain = 1:2
        for ss = 1:4
            cross_sim(human_brain,monkey_brain,ss) = corr(vec_monkey{monkey_brain},vec_human{human_brain,ss},'type','Spearman');
            for tt =  1:length(interested_time_point)
                cross_sim_t(human_brain,monkey_brain,ss,tt) = corr(vec_monkey_time{monkey_brain}(tt,:)',vec_human{human_brain,ss},'type','Spearman');
                tt
            end
        end
    end
end
%%
ldg = {};
data_to_show = squeeze(cross_sim(1,1,:));
figure;set(gcf,'Position',[550 450 550 280])
cc = colormap_matplotlib('Set1');

subplot(1,3,1); hold on
mean_data = mean(SIM_cross_human{1});
std_data = std(SIM_cross_human{1});
n = length(SIM_cross_human{1});
ci_lower = mean_data - (std_data / sqrt(n));
ci_upper = mean_data + (std_data / sqrt(n));
patch([0.1,3,3,0.1],[ci_lower,ci_lower,ci_upper,ci_upper],[0.8,0.8,0.8],'EdgeColor',[1,1,1]);
bar(mean(data_to_show),'FaceAlpha',0.3,'FaceColor',cc(1,:),'EdgeColor','k','LineWidth',1)
mean_data = mean(data_to_show);
std_data = std(data_to_show);
errorbar(1,mean_data,std_data / sqrt(4),'LineWidth',2,'Color','k')
scatter([0.8,0.8,1.2,1.2],data_to_show,'filled','MarkerEdgeColor','k');
ylabel('Cross-species similarity (r)')
xlabel('Average')
xlim([0,2])
ylim([-0.02,0.6])
xticks([])
set_font

subplot(1,3,[2,3]);
hold on
data_to_show = squeeze(cross_sim_t(1,1,:,:));
shadedErrorBar(interested_time_point,mean(data_to_show),std(data_to_show)/sqrt(4),'patchSaturation',0.3,'lineprops',{'Color',cc(1,:)})
plot(interested_time_point,mean(data_to_show),'LineWidth',1,'Color',cc(1,:))
ylim([-0.02,0.6])
xlabel('Time from onset(ms)')
t1 = interested_time_point(1);
t2 = interested_time_point(end);
patch([t1,t2,t2,t1],[ci_lower,ci_lower,ci_upper,ci_upper],[0.8,0.8,0.8],'EdgeColor',[1,1,1]);
xlim(interested_time_point([1,end]))
set_font
figsave(fullfile(root_dir,'Figs/F4/'),sprintf('Fig4B'))
%% Generate example rdm
all_data_fmri_general = zscore(all_data_fmri_general,0,2);
[coeff,score,latent,tsquared,explained,mu] = pca(all_data_fmri_general',NumComponents=2);
[~,order_here] = sort(score(:,1));
rdm1 = 1-squareform(pdist(all_data_fmri_general(:,order_here(1:10:1000))',"correlation"));
all_data_ephys = [];
for area_here = 1:max(neuron_area)
    if(bigROI_areawise(area_here)=='I')
        all_data_ephys = [all_data_ephys; rsp_here{area_here}];
    end
end
all_data_ephys = zscore(all_data_ephys,0,2);
rdm2 = 1-squareform(pdist(all_data_ephys(:,order_here(1:10:1000))',"correlation"));

cm=colormap_matplotlib('bwr');
figure;set(gcf,'Position',[20 250 1650 520])
subplot(1,2,1)
rdm1=tril(rdm1,-1);rdm1(rdm1==0)=max(rdm1(:));
imagesc(rdm1)
axis equal;colorbar;colormap(gca,cm(255:-1:128,:))
clim([-0.3,0.3]);set_font
subplot(1,2,2)
rdm2=tril(rdm2,-1);rdm2(rdm2==0)=max(rdm2(:));
imagesc(rdm2)
axis equal;colorbar;colormap(gca,cm(1:128,:))
clim([-0.3,0.3]);set_font
figsave(fullfile(root_dir,'Figs','F4'),'F4A')

%% Generate correlation map..
SS = 5;
hemi_here = 'lh';
fMRI_data = load(fullfile(fMRI_dataset_path,sprintf('S%d_%s_Rsp.mat',SS,hemi_here)));
fMRI_data = fMRI_data.mean_brain_data;
fMRI_data = double(fMRI_data)./300;
lh_data = fMRI_data;
hemi_here = 'rh';
fMRI_data = load(fullfile(fMRI_dataset_path,sprintf('S%d_%s_Rsp.mat',SS,hemi_here)));
fMRI_data = fMRI_data.mean_brain_data;
fMRI_data = double(fMRI_data)./300;
rh_data = fMRI_data;
area_array = neuron_area;ev_result = {};
rdm_lh = zeros([max(neuron_area),length(ROI_info.S5_lh_General)]);
rdm_rh = zeros([max(neuron_area),length(ROI_info.S5_rh_General)]);
for area_idx = 1:max(neuron_area)
    response_area = rsp_here{area_idx};
    response_area = zscore(response_area,0,2);
    mean_rsp = mean(response_area);
    sim_to_population = corr(mean_rsp', response_area',type="Spearman");
    minium_data = quantile(sim_to_population, 0.1);
    neuron_set{1}=find(sim_to_population<minium_data);
    max_data = quantile(sim_to_population, 0.333);
    neuron_set{2}=find(sim_to_population>max_data);
    neuron_set{3} = setdiff(1:length(sim_to_population), [neuron_set{1},neuron_set{2}]);
    for cc = 2
        mean_rsp_here = response_area(neuron_set{cc},:);
        mean_rsp_here = mean(mean_rsp_here);
        rdm_lh(area_idx,:)= corr(lh_data', mean_rsp_here','type','Pearson');
        rdm_rh(area_idx,:)=  corr(rh_data', mean_rsp_here','type','Pearson');
    end
    area_idx
end
save(fullfile(proc_dir,sprintf('RDM_S%d_splits.mat', 5)), 'rdm_lh', 'rdm_rh');