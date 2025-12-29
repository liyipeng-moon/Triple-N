clear;clc;close all
root_dir = 'C:\Users\moonl\Desktop\NNN';
cd(root_dir)
addpath(genpath(pwd));
[proc_dir,raw_dir] = gen_dirs(root_dir);
manual_data = readtable("exclude_area.xls");
stats_dir = fullfile(root_dir,"Figs/stats/");
fid = fopen(fullfile(stats_dir, 'F5.txt'),'w');
[s,s_name] = load_embedding;
%% Analysis for alexnet_fc6 and mpnet
Alex_FC_idx = 14;
Mpnet_idx = 6;
load('fMRI_result.mat')
Interested_Voxel = find(voxel_info_1~=0);
for vertex = 1:length(Interested_Voxel)
    fMRI_performance(vertex,1) = r(Interested_Voxel(vertex),Alex_FC_idx);
    fMRI_performance(vertex,2) = r(Interested_Voxel(vertex),Mpnet_idx);
end
%%
AVG_data = [];
TimeCourseData = [];
for ses = 1:102
    if(strcmp(manual_data.Area{ses},'IT'))
        EncodingData = dir(fullfile(proc_dir,sprintf("s5_encode_ses_%03d.mat",ses)));
        load(fullfile(proc_dir,EncodingData.name));
        Visual_Performance = pred_r_array(:,Alex_FC_idx);
        LLM_Performance = pred_r_array(:,Mpnet_idx);
        AVG_data = [AVG_data; [Visual_Performance,LLM_Performance]];

        R2 = min(pred_r2_array(:,[Mpnet_idx,Alex_FC_idx])')';
        TimeData = pred_r2_array_t(R2>0.05, [Mpnet_idx,Alex_FC_idx],:);
        TimeCourseData = [TimeCourseData; TimeData];

    end
    ses
end
%%

figure(1);
set(gcf,'Position',[30 200 1825 600])
subplot(2,5,6); pbaspect([1,1,1])
hold on
x_here = AVG_data(:,1);
y_here = AVG_data(:,2);
[beta,beta_CI] = demingRegression(x_here, y_here);

bin_scatter(x_here,y_here,[0:0.01:1],[0:0.01:1]);
plot([0,1], [0,beta],'Color','k','LineWidth',1)
xlabel('Visual model (r)');ylabel('Language model (r)')
title(sprintf('NNN performance\nLVR: %.02f',beta))
set_font

subplot(2,5,7); pbaspect([1,1,1])
hold on
x_here = fMRI_performance(:,1);
y_here = fMRI_performance(:,2);
[beta,beta_CI] = demingRegression(x_here, y_here);
bin_scatter(x_here,y_here,0:0.01:1,0:0.01:1);
plot([0,1], [0,beta],'k','LineWidth',1)
xlabel('Visual model (r)');ylabel('Language model (r)')
title(sprintf('NSD performance\nLVR: %.02f',beta))
set_font

clear a b
for mm = [1,2]
    data_here = squeeze(TimeCourseData(:,mm,:));
    [a(mm,:),b(mm,:)] = max(data_here');
end

figure(52)
subplot(3,4,11);
pbaspect([1,1,1])
hold on
bin_scatter(b(1,:),b(2,:),50:2:200,50:2:200);
cm_here = colormap_matplotlib('bwr');
cm_here = flipud(cm_here(floor([2:257]/2),:));
set(gca,'Colormap',cm_here)
xlabel('Visual peak latency (ms)')
ylabel('Language peak latency (ms)')
xticks(50:50:200);yticks(50:50:200)
set_font

figure(1)
%%
% ROI wise
interested_area = {'MB','AB','MF','AF','MO','AO','LPP','PITP','CLC','AMC'};
interested_area_name = {'Middle body','Anterior body','Middle face','Anterior face','Middle object','Anterior object','Scene1','Scene2','Middle Color','Anterior Color','Human vertex'};
figure(2)
set(gcf,'Position',[0 300 1100 330])
betas_bootstrap=[];
for aa = 1:length(interested_area)
    AVG_data = [];
    TimeCourseData = [];
    for search_area = 1:height(manual_data)
        try
            if(strcmp(manual_data.AREALABEL{search_area}(1:length(interested_area{aa})),interested_area{aa}))
                EncodingData = dir(fullfile(proc_dir,sprintf("s5_encode_ses_%03d.mat",search_area)));
                load(fullfile(proc_dir,EncodingData.name));
                Visual_Performance = pred_r_array(:,Alex_FC_idx);
                LLM_Performance = pred_r_array(:,Mpnet_idx);
                AVG_data = [AVG_data; [Visual_Performance,LLM_Performance]];

                R2 = min(pred_r2_array(:,[Mpnet_idx,Alex_FC_idx])')';
                TimeData = pred_r2_array_t(R2>0.05,[Mpnet_idx,Alex_FC_idx],:);
                TimeCourseData = [TimeCourseData; TimeData];

            end
        end
    end
    x_here = AVG_data(:,1);y_here = AVG_data(:,2);
    [beta_area_wise(aa),beta_CI_area(aa,:),betas_bootstrap(aa,:)] = demingRegression(x_here,y_here);
    ll = mod(1+aa,2)*5 + round((aa/2));
    subplot(2,5,ll); hold on
    for mm = [1,2]
        ccall = colormap_matplotlib('tab20',20);
        cc = ccall(mm,:);
        data_here = squeeze(TimeCourseData(:,mm,:));
        m = median(data_here,'omitnan');
        e = std(data_here,'omitnan')./sqrt(size(TimeCourseData,1));
        shadedErrorBar(1:350, m,e,'lineprops',{'Color',cc,'Linewidth',2});
    end
    yl = ylim;
    ylim([-0.05, yl(2)])
    xlabel('Time (ms)')
    xlim([0,300])
    ylabel('Accuracy (r)')
    title(interested_area_name{aa})
    add_onset_area(10,150)
    set_font
    [vmax1,vloc1]=max(squeeze(TimeCourseData(:,2,:))');
    [vmax2,vloc2]=max(squeeze(TimeCourseData(:,1,:))');
    delt=vloc2-vloc1;
    delt_save{aa} = delt;
end

legend({'Language','Visual'},'Box','off',Location='best')
drawnow
figsave(fullfile(root_dir,'Figs/F5'),'EncodingSupp')

figure(1)
x_here = fMRI_performance(:,1);
y_here = fMRI_performance(:,2);
[beta_area_wise(aa+1),beta_CI_area(aa+1,:)] = demingRegression(x_here, y_here);
subplot(2,5,8);
pbaspect([1,1,1])
hold on

CI = (beta_CI_area(:,2)-beta_CI_area(:,1))/2;

xticks_here = [1,2,1,2,1,2,1,2,1,2,3];
cc = colormap_matplotlib('tab20',20);
for bb_data = 1:length(xticks_here)
    scatter(xticks_here(bb_data), beta_area_wise(bb_data),'MarkerFaceColor',cc(4+bb_data,:),MarkerEdgeAlpha=0,Marker='square')
end
for bb_data = 1:length(xticks_here)
    errorbar(xticks_here(bb_data), beta_area_wise(bb_data), CI(bb_data),'CapSize',0,'LineWidth',2,Color=cc(4+bb_data,:));
end
for aa = 1:5
    m_idx = aa*2-1;
    a_idx = aa*2;
    plot([1,2], beta_area_wise([m_idx,a_idx]),LineWidth=1,Color=[0.5,0.5,0.5])
end
legend({'Body','','Face','','Object','','Scene','','Color','','Human Vertex'},Location="best")
xticks([1,2,3]); xlim([0.5,6.5])
xticklabels({'Middle IT', 'Anterior IT', 'Human VTC'})
ylim([0.55,1.05])
yticks([0.5:0.1:1]);
xtickangle(90)
ylabel('Slope')
colorbar

for aa = 1:5
    m_idx = aa*2-1;
    a_idx = aa*2;
    delta = betas_bootstrap(a_idx,:)-betas_bootstrap(m_idx,:);
    ci = prctile(delta, [0.5 99.5]);
    fprintf(fid,'%s, 99 CI = %.02f - %.02f\n', interested_area_name{aa*2-1}, ci(1), ci(2));

    beta_to_test(1,aa)=beta_area_wise(m_idx);
    beta_to_test(2,aa)=beta_area_wise(a_idx);

end

[h,p,ci,stats] = ttest(beta_to_test(1,:),beta_to_test(2,:));
fprintf(fid, 'MacaqueMA compare: t(%d)=%.02f, p=%.07f\n', stats.df, stats.tstat, p);
[h,p,ci,stats] = ttest(beta_to_test(:), beta_area_wise(end));
fprintf(fid, 'MacaqueHuman compare: t(%d)=%.02f, p=%.07f\n', stats.df, stats.tstat, p);

figure(52)
subplot(3,4,10);
pbaspect([1,1,1])
hold on
xticks_here = [1:10];
cc = colormap_matplotlib('tab20',20);
data_now = [];
array_now = [];
interested_combine = {[1,2],[3,4],[5,6],[7,8],[9,10]};
for combine_category = 1:length(interested_combine)
    data_diff = [];
    ll = interested_combine{combine_category};
    for bb_data = 1:length(ll)
        data_diff = delt_save{ll(bb_data)};
        data_now = [data_now,data_diff];
        array_now = [array_now,combine_category*ones([1,length(data_diff)])];
    end
end
boxplot(data_now,array_now,'Symbol','',MedianStyle='line',Whisker=0.5);
ylim([-20,50])
cc = colormap_matplotlib('tab20',20);
pretty_box(cc(13:-2:1,:))
lab = {'Body','Face','Object','Scene','Color'};
for cc  = 1:5
    test_data = data_now(array_now==cc);
    [h,p,ci,stats] = ttest(test_data);
    fprintf(fid,'Time Diff: %s, t(%d)=%.02f, p=%.1e \n', lab{cc}, stats.df, stats.tstat, p)
    scatter(cc, mean(test_data),'MarkerEdgeColor','k',MarkerFaceColor=[1,1,1])
end
yline(0,LineStyle="--",LineWidth=1)
ylabel('Time lag (ms)')
xticklabels(lab)
xtickangle(90)
set_font
%% Example unit
for aa = [1:2]
    figure(1);
    subplot(2,5,aa+3);
    pbaspect([1,1,1])
    hold on
    TimeCourseData = [];
    for search_area = 1:height(manual_data)
        try
            if(strcmp(manual_data.AREALABEL{search_area}(1:length(interested_area{aa})),interested_area{aa}))
                EncodingData = dir(fullfile(proc_dir,sprintf("s5_encode_ses_%03d.mat",search_area)));
                load(fullfile(proc_dir,EncodingData.name));
                R2 = min(pred_r2_array(:,[Mpnet_idx,Alex_FC_idx])')';
                TimeData = pred_r2_array_t(R2>0.05, [Mpnet_idx,Alex_FC_idx],:);
                TimeCourseData = [TimeCourseData; TimeData];
            end
        end
    end
    kk = squeeze(TimeCourseData(:,2,:));
    [kk, ~] =  max(kk');
    [~,order_unit] = sort(kk,'descend');
    LLM_units = order_unit(1);
    vv = []; xx = [];cc = colormap_matplotlib('tab20',20);
    for mm = [2,1]
        dd = conv(squeeze(TimeCourseData(LLM_units,mm,:))',ones([1,1]),'same')./1;
        [vv(mm),xx(mm)] = max(dd);
        plot(dd,'Color',cc(2*mm-1,:),'LineWidth',2)
    end
    for mm = [2,1]
        scatter(xx(mm),vv(mm)+0.05,'filled',Marker='v',MarkerEdgeAlpha=0,MarkerFaceColor=cc(2*mm-1,:));
    end
    xlim([0,250]);xlabel('Time (ms)');ylabel('Performance (r)')
    title(sprintf('Example %s unit\ntime lag = %d ms', interested_area_name{aa},diff(xx)))
    legend({'Visual','Language'},'Location','northeast',Box='off')
    colorbar
end

%%
load(fullfile(proc_dir,'S6.mat'))
t_result = load(fullfile(proc_dir,"S6_time_decoding.mat"));
fid = fopen(fullfile(root_dir,'Figs/stats','F5Decode.txt'),'w');

plotx = 2;ploty = 5;
cc_here = colormap_matplotlib('tab20',20);
for brain_here = 1:size(acc_save,1)
    for space_here = 1:size(acc_save,2)
        figure(1);
        cm_now = cc_here(2*(brain_here-1)+space_here,:);
        subplot(plotx,ploty,9);pbaspect([1,1,1]);hold on
        plot(num_series, acc_save{brain_here,space_here}, 'LineWidth', 2,'Color',cm_now)

        figure(52)
        subplot(3,4,space_here);hold on
        plot(1:size(FC_score,2), r_save{brain_here,space_here}, 'LineWidth', 2,'Color',cm_now)
        scatter(1:size(FC_score,2), r_save{brain_here,space_here}, 10, 'MarkerFaceColor',cm_now,'MarkerEdgeAlpha',0)
    end
end

figure(1)
subplot(2,5,9);hold on
plot(num_series, 100./num_series,'LineWidth',2,LineStyle=':',Color=[0.5 0.5 0.5])
ldg_pool{2} = 'Macaque Language';
ldg_pool{4} = 'Human Language';
legend(ldg_pool,'Box','off',Location='best')
ylim([0,100]);xlim([0,1000]);xticks(0:200:1000);xlabel('Number of images');ylabel('Decoding Accuracy (%)');xtickangle(0);set_font
set_font
colorbar

figure(52)
subplot(3,4,1);hold on
title('Decoding visual feature')
legend({'Macaque IT','','Human IT',''},Location="northeast",Box="off")
xlabel('# Component');ylabel('Accuracy (r)');xlim([0,50]);ylim([0,1]);set_font
set_font

subplot(3,4,2);hold on
title('Decoding language feature')
legend({'Macaque IT','','Human IT',''},Location="northeast",Box="off")
xlabel('# Component');ylabel('Accuracy (r)');xlim([0,50]);ylim([0,1]);set_font
set_font

subplot(3,4,3); hold on
max_pc = 20;
scatter(r_save{1,1}(1:max_pc),r_save{2,1}(1:max_pc),12,'filled');
scatter(r_save{1,2}(1:max_pc),r_save{2,2}(1:max_pc),12,'filled');
minmin = 0.1;
plot([minmin,0.9], [minmin,0.9],'LineStyle',':','Color','k')
legend({'Visual','Language'},'box','off',Location='best',Orientation='horizontal')
xlim([minmin,1]);ylim([minmin,1])
ylabel('Human performance')
xlabel('Macaque performance')
set_font

figure(52)
subplot(3,4,9);pbaspect([1,1,1]); hold on
mm = mean(acc_pool');ee = std(acc_pool')./sqrt(clus_size);
bar(3.5,mm(1),'FaceAlpha',0,EdgeColor=cc_here(1,:),LineWidth=2);
bar(4.5,mm(2),'FaceAlpha',0,EdgeColor=cc_here(3,:),LineWidth=2);
for cc = 1:clus_size
    plot([3.5,4.5], acc_pool(:,cc),'LineWidth',0.8,'Color',[0.6,0.6,0.6])
end
bar(1,acc_12(1),'EdgeAlpha',0,FaceColor=cc_here(1,:));
bar(2,acc_12(2),'EdgeAlpha',0,FaceColor=cc_here(3,:));
errorbar([3.5,4.5], mm,ee,'LineStyle','none',LineWidth=2,Color=[0,0,0])
xticks([1,2,3.5,4.5]);
xticklabels([xtitle_pool,xtitle_pool])
ylim([-5, 115]);yticks(0:20:100);ylabel('Decoding Accuracy (%)')
xlim([0.4,5])
plot([0.5,5],[100/12,100/12],'LineStyle',':','Color','k')
xtickangle(45)
set_font

[h,p,ci,stats] = ttest(acc_pool(1,:), acc_pool(2,:));
text(4,98,'***','HorizontalAlignment','center')
plot([3.5,4.5],[97,97],'LineWidth',1,'Color','k')
script = sprintf('cross-species: log-pval = %.1e, t(%d)=%.04f \n', p, stats.df,stats.tstat);
fprintf(fid, script)

[h,p,ci,stats] = ttest(acc_pool(1,:), acc_12(1));
script = sprintf('Macaque Compare: pval = %.1e, t(%d)=%.04f \n', p, stats.df,stats.tstat);
fprintf(fid, script)
text([1+3.5]/2,106,'*','HorizontalAlignment','center')
plot([1,3.5],[105,105],'LineWidth',1,'Color','k')

[h,p,ci,stats] = ttest(acc_pool(2,:), acc_12(2));
script = sprintf('Human Compare: pval = %.1e, t(%d)=%.04f \n', p, stats.df,stats.tstat);
fprintf(fid, script)
text([2+4.5]/2,114,'***','HorizontalAlignment','center')
plot([2,4.5],[113,113],'LineWidth',1,'Color','k')

set_font

acc_plot = [];
rr_plot = [];
for tt = 1:size(t_result.acc_here,1)
    for ss = 1:size(t_result.acc_here,2)
        acc_plot(tt,ss,:)=t_result.acc_here{tt,ss}(t_result.num_series);
        rr_plot(tt,ss,:)=t_result.r_here{tt,ss};
    end
end
decode_tin = t_result.decode_tin;
tnum_series = t_result.num_series;

for nn = 1:length(tnum_series)
    bb = 100/tnum_series(nn);
    ccbin = 5;
    figure(52);subplot(3,4,nn+4); hold on
    data_here = squeeze(acc_plot(1:end,1,nn));
    data_here = conv(data_here,ones([1,ccbin]),'same')/ccbin;
    plot(decode_tin,data_here,'LineWidth',2,Color='k',LineStyle='-')
    [a,b1] = max(data_here);
    scatter(decode_tin(b1),1.05*a,'filled','Marker','v',MarkerEdgeAlpha=0,MarkerFaceColor=[0,0,0])
    data_here = squeeze(acc_plot(1:end,2,nn));
    data_here = conv(data_here,ones([1,ccbin]),'same')/ccbin;
    plot(decode_tin,data_here,'LineWidth',2,Color=[0.7,0.7,0.7],LineStyle='-')
    [a,b2] = max(data_here);
    scatter(decode_tin(b2),1.05*a,'filled','Marker','v',MarkerEdgeAlpha=0,MarkerFaceColor=[0,0,0])
    xlim([0,300])
    xlabel('Time (ms)')
    ylabel('Performance (%)')
    yline(bb,'LineStyle','--','Color',[0.5,0.5,0.5])
    yl = ylim;
    ylim([-5, yl(2)])
    add_onset_area(20,150)
    if(nn==4)
        legend({'Visual','','Language','','',''},'Location','best',Box='off')
    end
    set_font
    title(sprintf('Number of images: %d \n Peak time lag = %d ms', tnum_series(nn),b2-b1));

    if(nn==1)
        figure(1)
        subplot(plotx,ploty,10); pbaspect([1,1,1]);hold on
        data_here = squeeze(acc_plot(1:end,1,nn));
        data_here = conv(data_here,ones([1,ccbin]),'same')/ccbin;
        plot(decode_tin,data_here,'LineWidth',2,Color='k',LineStyle='-')
        [a,b1] = max(data_here);
        scatter(decode_tin(b1),1.05*a,'filled','Marker','v',MarkerEdgeAlpha=0,MarkerFaceColor=[0,0,0])
        data_here = squeeze(acc_plot(1:end,2,nn));
        data_here = conv(data_here,ones([1,ccbin]),'same')/ccbin;
        plot(decode_tin,data_here,'LineWidth',2,Color=[0.7,0.7,0.7],LineStyle='-')
        [a,b2] = max(data_here);
        scatter(decode_tin(b2),1.05*a,'filled','Marker','v',MarkerEdgeAlpha=0,MarkerFaceColor=[0,0,0])
        xlim([0,300])
        xlabel('Time (ms)')
        ylabel('Performance (%)')
        yline(bb,'LineStyle','--','Color',[0.5,0.5,0.5])
        yl = ylim;
        ylim([0, yl(2)])
        add_onset_area(20,150)
        if(nn==4)
            legend({'Visual','','Language','','',''},'Location','best',Box='off')
        end
        title(sprintf('Number of images: %d \n Peak time lag = %d mss', tnum_series(nn),b2-b1));
        colorbar
        set_font
    end
end
drawnow

%%
figure(1)
figsave(fullfile(root_dir,'Figs/F5'),'F5')
figure(2)
figsave(fullfile(root_dir,'Figs/F5'),'F5S2')
figure(52)
set(gcf,'Position',[500 60 1150 850])
figsave(fullfile(root_dir,'Figs/F5'),'F5S1')

return
%% Supp
for ll = 1:8
    Language_model_series{ll} = s_name{ll};
end

Vision_model_series{1} = 'Resnet50-DINO';
Vision_model_idx{1} = 17:34;

Vision_model_series{2} = 'Alexnet-IN1K';
Vision_model_idx{2} = 9:16;

Vision_model_series{3} = 'Resnet50-IN1K';
Vision_model_idx{3} = 35:52;

Vision_model_series{4} = 'ViT-B-16-IN1K';
Vision_model_idx{4} = 53:65;

Vision_model_series{5} = 'InceptionV3-IN1K';
Vision_model_idx{5} = 66:83;
tic
for v = 1:5
    parfor language_idx = 1:8
        [LVR_m(v,language_idx),LVR_h(v,language_idx),TimeLag(v,language_idx)] = gen_species_difference(Vision_model_idx{v},language_idx);
    end
end
toc
%%
close all; figure(999); set(gcf,'Position',[1 41 1050 210])

colors  = colormap_matplotlib('Set2',8); markers = {'o','s','^','d','v','>','<','p'};
subplot(1,4,1); hold on
for v = 1:5
    for l = 1:8
        scatter(LVR_h(v,l),LVR_m(v,l),'filled','MarkerFaceColor',colors(v,:),'Marker',markers{l})
    end
end
xlim([0.6,1]);ylim([0.6,1]);xticks([0.6:0.2:1]);yticks([0.6:0.2:1]);xlabel('Human LVR');ylabel('Macaque LVR');set_font

plot([0.88,1],[0.64,0.64],'Color','k')
plot([0.88,1],[0.75,0.75],'Color','k')
plot([0.88,0.88],[0.64,0.75],'Color','k')
plot([1,1],[0.64,0.75],'Color','k')


subplot(1,4,2); hold on
for v = 1:5
    for l = 1:8
        scatter(LVR_h(v,l),LVR_m(v,l),'filled','MarkerFaceColor',colors(v,:),'Marker',markers{l})
    end
end
title(sprintf('corr=%.02f',corr(LVR_m(:),LVR_h(:))))
xlim([0.89,1]);ylim([0.64,0.75]);yticks([0.65:0.05:0.75]);xticks([0.9:0.05:1]);xlabel('Human LVR');ylabel('Macaque LVR');set_font

subplot(1,4,3); hold on
for v = 1:5
    scatter(0,0,'filled','MarkerFaceColor',colors(v,:),'Marker',markers{1})
end
legend(Vision_model_series); axis off; set_font

subplot(1,4,4); hold on
for l = 1:8
    scatter(0,0,'filled','MarkerFaceColor','k','Marker',markers{l})
end
legend(Language_model_series); axis off; set_font


figsave(fullfile(root_dir,"Figs/F5/"),'compareModel')