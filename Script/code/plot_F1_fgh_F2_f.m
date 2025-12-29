clear
clc
close all
root_dir = 'C:\Users\moonl\Desktop\NNN';
[proc_dir,raw_dir] = gen_dirs(root_dir);
addpath(genpath(pwd));

reliability_thres = 0.4;
manual_data = readtable("exclude_area.xls");

fanof = [];
reliability = [];
Sparseness = [];
UnitType = [];
session_array = [];
area_array = [];
rsp_all = [];
snr_array = [];
for ses_now = [1:90]
    proc1_file_name = dir(fullfile(proc_dir,sprintf('Processed_ses%02d*', ses_now)));
    proc1_file_name = proc1_file_name.name;
    pro1_data = load(fullfile(proc_dir,proc1_file_name));
    unit_size = length(pro1_data.reliability_basic);
    session_array = [session_array, ses_now*ones([1, unit_size])];

    ses_area_array = zeros([1, unit_size]);
    all_area = find(manual_data.SesIdx==ses_now);

    % Compute Fano Factor
    filename_here = dir(fullfile(raw_dir,sprintf('ses%02d*h5',ses_now)));
    filename_here = filename_here.name;
    metaname_here = dir(fullfile(raw_dir,sprintf('ses%02d*mat',ses_now)));
    metaname_here = metaname_here.name;
    meta_data = load(fullfile(raw_dir,metaname_here));
    RasterData = h5read(fullfile(raw_dir,filename_here), '/raster_matrix_img');
    t1 = median(pro1_data.best_r_time1(pro1_data.reliability_best>0.4));
    t2 = median(pro1_data.best_r_time2(pro1_data.reliability_best>0.4));
    t_array = meta_data.global_params.pre_onset + [t1:t2];
    RasterSum = squeeze(sum(RasterData(:,:,t_array),3));
    Fanofactor_across_img = zeros([size(RasterSum,1),1000]);
    for img = 1:1000
        trial_here = find(meta_data.img_idx==img);
        Rsp_here = RasterSum(:,trial_here);
        Var_across_trial = var(Rsp_here,[],2);
        Mean_across_trial = mean(Rsp_here,2);
        Fanofactor_across_img(:,img) = Var_across_trial./Mean_across_trial;
    end

    %Add...
    for area_idx = all_area'
        if(strcmpi(manual_data.AREALABEL{area_idx},'Unknown'))
            continue;
        end
        x1 = manual_data.y1(area_idx);
        x2 = manual_data.y2(area_idx);
        unit_here = pro1_data.pos>x1 & pro1_data.pos<x2;
        ses_area_array(unit_here)=area_idx;
    end

    area_array = [area_array, ses_area_array];
    reliability = [reliability, pro1_data.reliability_best];
    rsp_here = pro1_data.response_best(:,1:1000);
    rsp_all = [rsp_all;rsp_here];

    fanof = [fanof, mean(Fanofactor_across_img','omitnan')];
    mean_rsp = mean(rsp_here,2).^2;
    mean_rsp_sq = mean(rsp_here.^2,2);
    Sparseness = [Sparseness, 1 - mean_rsp'./mean_rsp_sq'];

    UnitType = [UnitType, pro1_data.UnitType];
    snr_array = [snr_array, pro1_data.snrmax];

    ses_now
end
%%
interested_area = {'MB','AB','MF','AF','MO','AO','LPP','PITP','CLC','AMC','Unknown'};
interested_name = {'MiddleBody','AnteriorBody','MiddleFace','AnteriorFace','MiddleObject','AnteriorObject','Scene1','Scene2','MiddleColor','AnteriorColor','Unknown'};

interested_area = fliplr(interested_area);
interested_name = fliplr(interested_name);
reliability_here = [];
fanof_here = [];
Sparseness_vec = [];
sim_to_population = [];
area_label = [];
snr_here = [];
area_location_pool = {};
for aa = 1:length(interested_area)
    area_location = [];
    if(aa==1) % unknown
        area_location = [area_location, find(area_array'==0)];
    end
    for search_area = 1:height(manual_data)
        try
            if(strcmp(manual_data.AREALABEL{search_area}(1:length(interested_area{aa})),interested_area{aa}))
                area_location = [area_location; find(area_array'==search_area)];
            end
        end
    end
    area_location_pool{aa}=area_location;
    reliability_here = [reliability_here, reliability(area_location)];
    area_label = [area_label, aa*ones([1, length(area_location)])];
    LOC = intersect(area_location, find(reliability>0.4));
    fanof_here = [fanof_here, fanof(LOC)];
    Sparseness_vec = [Sparseness_vec, Sparseness(LOC)];
    snr_here = [snr_here, snr_array(LOC)];
end

cm_here = colormap_matplotlib('Set3',length(interested_area));
figure
set(gcf,'Position',[500 600 1000 220])
subplot(1,5,1); hold on
boxplot(reliability_here,area_label,'Symbol','','Whisker',0.5,'orientation','horizontal');
yticks(1:length(interested_area));
yticklabels(interested_name)
set_font

subplot(1,5,2); hold on
boxplot(reliability_here,area_label,'Symbol','','Whisker',0.5,'orientation','horizontal');
xlim([-0.1,1])
xticks([0:0.2:1])
pretty_box(cm_here)
ylim(ylim+0.25)
yticks([]);
title('Reliability (r)')
set_font

subplot(1,5,3); hold on
boxplot(fanof_here,area_label(reliability_here>0.4),'Symbol','','Whisker',0.5,'orientation','horizontal');
yticks([]);
ylim(ylim+0.25)
xlim([0.8,1.6])
pretty_box(cm_here)
title('Fano factor')
set_font

subplot(1,5,4); hold on
boxplot(Sparseness_vec,area_label(reliability_here>0.4),'Symbol','','Whisker',0.5,'orientation','horizontal');
yticks([]);
ylim(ylim+0.25)
xlim([0,1])
pretty_box(cm_here)
title('Sparseness')
set_font

subplot(1,5,5); hold on
boxplot(snr_here,area_label(reliability_here>0.4),'Symbol','','Whisker',0.5,'orientation','horizontal');
pretty_box(cm_here)
yticks([]);
xlim([3,50])
ylim(ylim+0.25)
title('SNRmax')
set_font

figsave(fullfile(root_dir,'Figs\F1'),sprintf('F1low'))
%%
figure
interested_area = {'MB1','MB2','MB3','AB1','AB3','MF1','MF3','AF1','AF3','MO1s1','MO1s2','MO2','MO5','AO2','AO5','LPP4','PITP3','PITP4','CLC3','AMC3'};
interested_name = {'MBody-M1','MBody-M2','MBody-M3','ABody-M1','ABody-M3','MFace-M1','MFace-M3','AFace-M1','AFace-M3','MObjectSite1-M1','MObjectSite2-M1','MObject-M2','MObject-M5','AObject-M2','AObject-M5','Place1-M4','Place2-M3','Place2-M4','CentralColor-M3','AnteriorColor-M3'};
area_wise_rsp = [];
for aa = 1:length(interested_area)
    area_location = [];
    for search_area = 1:height(manual_data)
        try
            if(strcmp(manual_data.AREALABEL{search_area}(1:length(interested_area{aa})),interested_area{aa}))
                area_location = [area_location; find(area_array'==search_area)];
            end
        end
    end
    area_location_pool{aa}=area_location;
    area_rsp = zscore(rsp_all(area_location,:),0,2);
    area_wise_rsp(aa,:)=mean(zscore(rsp_all(area_location,:),0,2));
    within_sim = corr(mean(area_rsp(1:2:end,:))',mean(area_rsp(2:2:end,:))');
    within_area(aa)=within_sim;
end


BIG_RDM = 1-squareform(pdist(area_wise_rsp,'correlation'));
for aa = 1:length(interested_area)
    BIG_RDM(aa,aa)=within_area(aa);
end
imagesc(BIG_RDM)
cc = colorbar;
colormap(give_me_orange_bao)
clim([-0.6, 0.6])
xticks(1:length(interested_area))
yticks(1:length(interested_area))
xticklabels(interested_name)
yticklabels(interested_name)
axis equal
fig_here = gca;
fig_here.XAxis.TickDirection="none";
fig_here.YAxis.TickDirection="none";
fig_here.XAxis.FontSize=8;
fig_here.YAxis.FontSize=8;
xlim([0.5,length(interested_area)+0.5])
ylim([0.5,length(interested_area)+0.5])
cc.Ticks=[-0.6,0,0.6];
cc.Label.String='Correlation';
set_font
box on
xtickangle(90)
figsave(fullfile(root_dir,'Figs\F2'),sprintf('F2F'))

return
%% Do tSNE
rsp_pool = zscore(rsp_all(reliability>0.4,:),0,2);
rng(1009)
t_emb = tsne(rsp_pool');
t_emb = t_emb-min(t_emb(:));
t_emb = floor(t_emb*250+100);
load img_pool.mat
big_img = 255*uint8(ones([8000,8000,3]));
for iii = 1:1000
    t1 = t_emb(iii,1);
    t2 = t_emb(iii,2);
    big_img(t1:t1+226, t2:t2+226,:) = img_pool{iii};
end
figure;imshow(big_img)
imwrite(big_img,fullfile(root_dir,"Figs/F2S",'F2tsne.png'))