close all
clear;clc
root_dir = 'C:\Users\moonl\Desktop\NNN';
cd(root_dir)
addpath(genpath(pwd));
[proc_dir,raw_dir] = gen_dirs(root_dir);
reliability_thres = 0.4;
manual_data = readtable("exclude_area.xls");
stats_dir = fullfile(root_dir,"Figs/stats/");
mkdir Figs\F3\
%%
reliability = {};
UnitType = {};
psth_ses = {};
psth_pool = {};
interested_time_point = 1:300;
interested_time_point_rdm = 1:5:350;
meta_example = load("ses10_240726_M5_3_info.mat"');
%%
for area_now = 1:102
    if(strcmp(manual_data.Area{area_now},'EVC'))
        continue
    end
    ses_idx = manual_data.SesIdx(area_now);
    proc1_file_name = dir(fullfile(proc_dir,sprintf('Processed_ses%02d*', ses_idx)));
    proc1_file_name = proc1_file_name.name;
    pro1_data = load(fullfile(proc_dir,proc1_file_name));

    filename_here = dir(fullfile(raw_dir,sprintf('ses%02d*h5',manual_data.SesIdx(area_now))));
    PSTH_HERE = h5read(fullfile(raw_dir,filename_here.name), '/response_matrix_img');

    x1 = manual_data.y1(area_now);
    x2 = manual_data.y2(area_now);
    pos_now = pro1_data.pos;
    unit_here = find(pos_now>x1 & pos_now<x2);

    reliability{area_now} = pro1_data.reliability_best(unit_here);
    rsp_here{area_now} = pro1_data.response_best(unit_here,1:1000);
    UnitType{area_now} = pro1_data.UnitType(unit_here);
    psth_ses{area_now} = pro1_data.mean_psth(unit_here, meta_example.global_params.pre_onset+interested_time_point);


    pos_ses{area_now} = pro1_data.pos(unit_here);
    psth_pool{area_now}=PSTH_HERE(unit_here, 1:1000, meta_example.global_params.pre_onset+interested_time_point_rdm);

    disp(area_now)
end
%%
rng(1009)
interested_area = {'MB','AB','MF','AF','MO','AO','LPP','PITP','CLC','AMC','Unknown'};
interested_name = {'MiddleBody','AnteriorBody','MiddleFace','AnteriorFace','MiddleObject','AnteriorObject','Scene1','Scene2','MiddleColor','AnteriorColor','Unknown'};
all_mean_psth = [];
reliability_here = [];
UnitType_vec = [];
area_idx = [];
brain_area = [];
pos_all = [];
psth_all = [];
psth_even = [];
psth_odd = [];
for aa = 1:length(interested_area)
    for search_area = 1:height(manual_data)
        try
            if(strcmp(manual_data.AREALABEL{search_area}(1:length(interested_area{aa})),interested_area{aa}))
                reliability_here = [reliability_here,reliability{search_area}];
                UnitType_vec = [UnitType_vec, UnitType{search_area}];
                all_mean_psth = [all_mean_psth; psth_ses{search_area}];
                area_idx = [area_idx, aa*ones([1, length(reliability{search_area})])];
                brain_area = [brain_area, search_area*ones([1, length(reliability{search_area})])];
                pos_all = [pos_all, pos_ses{search_area}];
                psth_all = [psth_all; psth_pool{search_area}];
                psth_even = [psth_even; psth_pool{search_area}(:,1:2:1000,:)];
                psth_odd = [psth_odd; psth_pool{search_area}(:,2:2:1000,:)];
            end
        end
    end
end
Interested_units = find(reliability_here>0.4);
un_nornalize_path = all_mean_psth(Interested_units,:);
all_mean_psth = zscore(all_mean_psth(Interested_units,:),0,2);
pos_all = pos_all(Interested_units);
area_all = brain_area(Interested_units);
psth_all = psth_all(Interested_units,:,:);

%% CV
keyboard
psth_u1 = squeeze(mean(psth_even(Interested_units,:,:),2));
psth_u2 = squeeze(mean(psth_odd(Interested_units,:,:),2));
psth_um = [psth_u2+psth_u1]/2;
psth_1 = zscore(psth_u1,0,2);
psth_2 = zscore(psth_u2,0,2);

RDM1 = squareform(pdist(psth_1,'correlation'));
RDM2 = squareform(pdist(psth_2,'correlation'));

[Y1,~] = cmdscale(double(RDM1),2);
[Y2,~] = cmdscale(double(RDM2),2);
rng(1009)
idx = kmeans(psth_1, 3);



idx(idx==1)=4;
idx(idx==3)=1;
idx(idx==2)=3;
idx(idx==4)=2;

%%
figure; 
set(gcf,'Position',[100 1 200 1000])
cm_here = colormap_matplotlib('Set1',9);

cm_here = cm_here([2,3,1],:);

for k = 1:3
    subplot(4,1,1); hold on
    plot(interested_time_point_rdm,mean((psth_1(idx==k,:))),Color=cm_here(k,:),LineWidth=2);
    
    subplot(4,1,2); hold on
    plot(interested_time_point_rdm,mean((psth_2(idx==k,:))),Color=cm_here(k,:),LineWidth=2);
    
    subplot(4,1,4); hold on
    scatter(Y1(idx==k,1),Y1(idx==k,2),2,MarkerFaceColor=cm_here(k,:),MarkerEdgeAlpha=0)
    
    subplot(4,1,3); hold on
    plot(interested_time_point_rdm,mean((psth_um(idx==k,:))),Color=cm_here(k,:),LineWidth=2);
end

subplot(4,1,1); hold on
xlim([0,300]);title('Resoonse to odd stimuli (train)');xlabel('Time (ms)');ylabel('Norm. Rsp. (a.u.)')

subplot(4,1,2); hold on
xlim([0,300]);title('Resoonse to even stimuli (test)');xlabel('Time (ms)');ylabel('Norm. Rsp. (a.u.)')

subplot(4,1,4); hold on
xticks([]);yticks([]);xlabel('MDS1');ylabel('MDS2');title('PSTH MDS to odd stimuli')
       
subplot(4,1,3); hold on
xlim([0,300]);title('Raw Firing Rate');xlabel('Time (ms)');ylabel('Un-Norm. Rsp. (Hz)')

figsave(fullfile(root_dir,"Figs/F3/"),'CrossValidate')
%%
rng(1009)
best_k = 3;
idx = kmeans(all_mean_psth, best_k);
idx(idx==1)=4;
idx(idx==3)=1;
idx(idx==4)=3;
%%

%%
rdm_tt = 1:size(psth_all,3);
rdm_vec_cc = [];
for cc = 1:best_k
    this_cc = find(idx==cc);
    pop_rsp = psth_all(this_cc,:,:);
    for t1 = 1:length(rdm_tt)
        t1_here = rdm_tt(t1);
        rsp1 = zscore(squeeze(pop_rsp(:,:,t1_here)),0,2);
        rdm1 = squareform(pdist(rsp1','correlation'));
        rdm_vec_cc(t1,cc,:) = rdm1(find(tril(rdm1)));
        t1
    end
end

figure
set(gcf,'Position',[750 200 850 850])
for x1 = 1:3
    for y1 = 1:3
        subplot(4,3,(x1-1)*3+y1); hold on
        rdm = [];
        for t1 = 1:length(rdm_tt)
            for t2 = 1:length(rdm_tt)
                rdm(t1,t2) = corr(squeeze(rdm_vec_cc(t1,x1,:)),squeeze(rdm_vec_cc(t2,y1,:)));
            end
            t1
        end
        tt_plot = interested_time_point_rdm(rdm_tt);
        imagesc(tt_plot,tt_plot,rdm)
        clim([0,0.8])
        contour(tt_plot,tt_plot,rdm,[0.2,0.4],'LineColor',[1,1,1]*0.5,'LineStyle','-')
        ylabel(sprintf('Cluster %d time', x1))
        xlabel(sprintf('Cluster %d time', y1))
        xticks([50:100:350])
        xticks([50:100:350])
        axis equal
        drawnow
    end
end
subplot(4,3,[10,11,12])
cb=colorbar;
cb.Location='northoutside';
colormap(flipud(colormap_matplotlib('spectral')))
clim([0,0.8])
cb.Label.String='Correlation (r)';
axis off

figsave(fullfile(root_dir,"Figs/F3/"),sprintf('F3S_Cross'))

%%
rdm_tt = 1:70;
for cc = 1:best_k
    this_cc = find(idx==cc);
    pop_rsp = psth_all(this_cc,:,:);
    rdm_vec = [];
    for t1 = 1:length(rdm_tt)
        t1_here = rdm_tt(t1);
        rsp1 = zscore(squeeze(pop_rsp(:,:,t1_here)),0,2);
        rdm1 = squareform(pdist(rsp1','correlation'));
        rdm_vec(t1,:) = rdm1(find(tril(rdm1)));
        t1
    end
    rdm_matrix{cc} = squareform(pdist(rdm_vec,'spearman'));
end
%%
figure(987)
set(gcf,'Position',[500 200 1250 425])
subplot(2,5,[1,6]); hold on
[a,b]=sort(idx);
yy_loc = find(diff(a)~=0);
yline(yy_loc,LineWidth=1,LineStyle="--")
imagesc(all_mean_psth(b,:));
ylim([0, length(idx)])
% colorbar('northoutside')
colormap(give_me_orange_bao)
clim([-3,3])
xlabel('Time (ms)')
yticks([]);
reliability_cluster = {};
cm_here = colormap_matplotlib('Set1',9);
set_font
ylb = [];
title_pool = {};
for cc_idx = 1:best_k
    subplot(2,5,2); hold on
    data_here = all_mean_psth(idx==cc_idx,:);
    data_here = squeeze(mean(psth_all(cc_idx==idx,:,:),2));
    data_here = zscore(data_here,0,2);
    m = mean(data_here);
    e = 1.96*std(data_here)./sqrt(size(data_here,1));
    shadedErrorBar(interested_time_point_rdm, m,e,'Lineprops',{'Color', cm_here(cc_idx,:),'LineWidth',1},'Patchsaturation',1)
    reliability_cluster{cc_idx} = reliability_here(Interested_units(cc_idx==idx));
    ylb = [ylb, sprintf('Cluster%d    ',cc_idx)];
    xlabel('Time (ms)')
    add_onset_area(20,150)
    xticks(50:100:350);
    xlim([0,350])
    if(cc_idx==1)
        ylabel('Norm. firing rate(a.u.)')
    end
    title_pool{cc_idx} = sprintf('Cluster%d, n = %d',cc_idx,size(data_here,1));
end
xlim([0,300])
legend({'Cluster1','','Cluster2','','Cluster3'},'Box','off',Location='best')
set_font

for cc = 1:3
    subplot(2,5,6+cc);hold on
    ylabel('Time (ms)')
    imagesc(interested_time_point_rdm,interested_time_point_rdm,1-rdm_matrix{cc});
    contour(interested_time_point_rdm,interested_time_point_rdm,1-rdm_matrix{cc},[0.3,0.5],'LineColor',[1,1,1]*0.4,'LineStyle','-')
    set(gca,'Colormap',flipud(colormap_matplotlib('spectral')))
    clim([0,1])
    axis equal
    xlim([0,350])
    ylim([0,350])
    xticks(50:100:350);
    yticks(50:100:350);
    xlabel('Time (ms)')
    box off
    title(title_pool{cc})
    set_font
end

subplot(2,5,[1,6])
ylabel(ylb)
ratio_here = [];
for aa = 1:max(area_idx)
    area_loc_here = find(area_idx(Interested_units)==aa);
    for cc = 1:best_k
        ratio_here(aa,cc) = sum(idx(area_loc_here)==cc)./length(area_loc_here);
    end
end
set_font

subplot(2,5,3)
b=barh(flipud(ratio_here),'stacked','FaceColor','flat',EdgeAlpha=0);
for area_size = 1:length(interested_area)
    for cluster_here = 1:best_k
        b(cluster_here).CData(area_size,:) = cm_here(cluster_here,:);
    end
end
xticks([0,1])
yticks(1:length(interested_area))
yticklabels(fliplr(interested_name))
xlabel('Ratio')
set_font

subplot(2,5,4)
violinplot(reliability_here(Interested_units),idx, ...
    'ViolinColor',cm_here,'MarkerSize',2,'ShowWhiskers',true,'ShowData',true,'ViolinAlpha',0.7);
xlim([0.5, best_k+0.5])
ylabel('Reliability')
xticks([1,2,3])
ylim([0.4,1])
xticklabels({'C1','C2','C3'})
box off
set_font
figsave(fullfile(root_dir,'Figs/F3/'),sprintf('F3A'))
%%

MI_area = [];
MI_perm = [];
for area_here = setdiff(1:102,[18,19,21,79])
    if(strcmp(manual_data.Area{area_here},'EVC'))
        continue
    end

    figure
    set(gcf,'Position',[700 750 900 240])
    subplot(1,3,1);hold on
    depth_here = [manual_data.y1(area_here),manual_data.y2(area_here)];

    location_here = find(area_all==area_here);
    pos_area = pos_all(location_here);
    ed = min(pos_area):150:max(pos_area);
    idx_area = idx(location_here);
    for cc = 1:best_k
        cc_Loc = find(cc==idx_area);
        histogram(pos_area(cc_Loc),ed,FaceColor=cm_here(cc,:),EdgeAlpha=0, FaceAlpha=0.5,Normalization="count",Orientation="horizontal")
    end

    yl=depth_here;
    yl = [min(pos_area), max(pos_area)];
    ylim(yl)
    ylabel('Depth (um)')
    xlabel('# Units')
    legend('Cluster1','Cluster2','Cluster3',Location='best',box='off')

    filename_here = dir(fullfile(raw_dir,sprintf('ses%02d*h5',manual_data.SesIdx(area_here))));
    LFPData = h5read(fullfile(raw_dir,filename_here.name), '/LFP_Data');
    meta_here = dir(fullfile(raw_dir,sprintf('ses%02d*info.mat',manual_data.SesIdx(area_here))));
    meta_data = load(meta_here.name);


    subplot(1,3,2); hold on
    lfp_data_points = floor((meta_data.global_params.pre_onset+interested_time_point)/2);
    lfp_data_points = unique(lfp_data_points);
    lfp_all_data = squeeze(mean(LFPData(1:1000,:,lfp_data_points),1))';
    imagesc(0:2:300,meta_data.LFP_META.depth_vals,lfp_all_data')
    ylim(yl);
    xlabel('Time (ms)')
    set(gca,"YLim",yl);
    cc_bar=colorbar;
    set(gca,"Colormap",colormap_matplotlib('plasma'))
    title('Local Field Potential (uv)')
    
    MI_area(area_here)= mutualinfo(idx_area,pos_area');
    pp_MI = [];
    for pp = 1:100
        pp_MI(pp) = mutualinfo(idx_area,pos_area(randperm(length(pos_area)))');
    end
    MI_perm(area_here)= median(pp_MI);

    sgtitle(sprintf('session %d, area = %s, MI=%.02f',manual_data.SesIdx(area_here),manual_data.AREALABEL{area_here}, MI_area(area_here)));
    set_font
    figsave(fullfile(root_dir,'Figs/F3/F3S/'),sprintf('F3S_%d',area_here))
end
%%
cc_here = colormap_matplotlib('Oranges',30);
cc_here = cc_here(12,:);
figure(987)

subplot(2,5,[5,10]); hold on

MI_area(MI_area==0)=[];
MI_perm(MI_perm==0)=[];

boxplot(MI_area,'positions', 1, 'Widths', 0.3, 'Colors', [1,1,1],'Whisker',0,'Symbol','')
boxplot(MI_perm,'positions', 2, 'Widths', 0.3, 'Colors', [1,1,1],'Whisker',0,'Symbol','')
h = findobj(gca, 'Tag', 'Box');
for i = 1:length(h)
    patch(get(h(i), 'XData'), get(h(i), 'YData'), cc_here, 'FaceAlpha',0.6);
end
scatter(ones([1,length(MI_area)]),MI_area,12,"filled",MarkerFaceColor=[0.5,0.5,0.5],MarkerFaceAlpha=0.4);
scatter(2*ones([1,length(MI_area)]),MI_perm,12,"filled",MarkerFaceColor=[0.5,0.5,0.5],MarkerFaceAlpha=0.4);
for aa = 1:length(MI_area)
    plot([1,2],[MI_area(aa),MI_perm(aa)],'Color',[0.5,0.5,0.5,0.25]);
end
xlim([0.5, 3.5])
xticks([1,2])
xticklabels({'Observed','Permutation'})
ylabel('MI(cluster, position)')
[p,h,stats] = ranksum(MI_area,MI_perm,tail="both");
formatted_p = sprintf('p-value: %e', p);
title(formatted_p,Interpreter="none")
ylim([0,0.5])
box off
set_font

set(gcf,'Position',[200 250 1200 375])

figsave(fullfile(root_dir,'Figs/F3'),sprintf('F3AS3'))

%% ClusterSize
figure
set(gcf,'Position',[500 600 300 300])
all_methods = {'DaviesBouldin'};
kk_list = 1:8;
for method_here = 1:length(all_methods)
    subplot(1, length(all_methods), method_here)
    hold on
    eval_data = evalclusters(all_mean_psth,"kmeans",all_methods{method_here},"KList",kk_list);
    plot(kk_list, eval_data.CriterionValues,'k','LineWidth',1)
    scatter(kk_list, eval_data.CriterionValues,'k')
    scatter(eval_data.OptimalK, eval_data.CriterionValues(eval_data.OptimalK),'r','filled')
    xlim([1.5, kk_list(end)-0.5])
    xticks([2:kk_list(end)]);
    xlabel('Cluster number')
    title(all_methods{method_here})
    drawnow
end
set_font
figsave(fullfile(root_dir,"Figs/F3"),sprintf('F3AS1'))

%% For TrialWise
all_cluster = zeros([1,length(reliability_here)]);
all_cluster(Interested_units)=idx;
for brain_area_idx = 1:max(brain_area)
    brain_loc = find(brain_area==brain_area_idx);
    clus_save{brain_area_idx}=all_cluster(brain_loc);
end
save(fullfile(proc_dir,'Clus_savee.mat'),'all_cluster','clus_save')