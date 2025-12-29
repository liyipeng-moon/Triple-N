clear;clc;
close all
root_dir = 'C:\Users\moonl\Desktop\NNN';
cd(root_dir);
addpath(genpath(pwd));

[proc_dir,raw_dir] = gen_dirs(root_dir);
mkdir Figs\
mkdir Figs\F2\

interested_ses = 1:90;

reliability_thres = 0.4;
color_pool = my_colormap('FOB');
load img_pool.mat
time_to_plot = 0:300;
load Clusinfo.mat
face_idx = 1001:1024;
body_idx = 1000+[26:31,43:48,50:61];
obj_idx = setdiff(1025:1072, body_idx);
loc_idx = [face_idx, body_idx, obj_idx];
noise_corr = []; noise_corr_old_ver = [];
%% Gen PC1 for visualization
rng(1009)
clus_mean = [];
for ses_now = 1:length(interested_ses)
    proc1_file_name = dir(fullfile(proc_dir,sprintf('Processed_ses%02d*', interested_ses(ses_now))));
    proc1_file_name = proc1_file_name.name;
    pro1_data = load(fullfile(proc_dir,proc1_file_name));
    unit_plot_here = find(pro1_data.reliability_best>reliability_thres);
    rsp_here = pro1_data.response_best(unit_plot_here,1:1000);
    rsp_here = zscore(rsp_here,0,2);
    area_mean = [];
    for cc = 1:max(Cluster_idx)
        area_mean(cc,:) = mean(rsp_here(:,Cluster_idx==cc),2);
    end
    clus_mean = [clus_mean,area_mean];
end
[coeff,score,latent,tsquared,explained,mu] = pca(clus_mean,NumComponents=1);
[~,clus_order] = sort(score);
big_img = [];
all_color = colormap_matplotlib('jet',max(Cluster_idx));
for cc = 1:max(Cluster_idx)
    imgs_here = find(Cluster_idx==clus_order(cc));
    mm = [];
    for ee = 1:5
        mm = [mm, img_pool{imgs_here(ee)}];
    end
    big_img = [big_img;add_edge(mm,floor(255*all_color(cc,:)),10)];
end
figure;
set(gcf,'Position',[400 400 1475 600])
subplot(1,2,1)
imshow(big_img)
subplot(1,2,2)
hold on
for cc = 1:max(Cluster_idx)
    clus_loc = find(Cluster_idx==clus_order(cc));
    scatter(LLM_tsne(clus_loc,1),LLM_tsne(clus_loc,2),15,MarkerFaceColor=all_color(cc,:),MarkerEdgeAlpha=0)
end

figsave(fullfile(root_dir, 'Figs\F2'),sprintf('F2S_exampleClus'))
close all
%%
nc_array = [];
sc_array = [];
%%
close all

manual_data = readtable("exclude_area.xls");
interested_ses = 1:90;
for ses_now = 1:90
    figure
    set(gcf,'Position',[300 90 380 780])
    proc1_file_name = dir(fullfile(proc_dir,sprintf('Processed_ses%02d*', interested_ses(ses_now))));
    proc1_file_name = proc1_file_name.name;
    pro1_data = load(fullfile(proc_dir,proc1_file_name));

    filename_here = dir(fullfile(raw_dir,sprintf('ses%02d*h5',interested_ses(ses_now))));
    filename_here = filename_here.name;
    metaname_here = dir(fullfile(raw_dir,sprintf('ses%02d*mat',interested_ses(ses_now))));
    metaname_here = metaname_here.name;
    meta_data = load(fullfile(raw_dir,metaname_here));

    RasterData = h5read(fullfile(raw_dir,filename_here), '/raster_matrix_img');
    PSTHData = h5read(fullfile(raw_dir,filename_here), '/response_matrix_img');
    LFPData = h5read(fullfile(raw_dir,filename_here), '/LFP_Data');

    time_to_find = meta_data.global_params.pre_onset+time_to_plot;
    subplot(4,2,[1,2]); hold on
    area = find(manual_data.SesIdx==interested_ses(ses_now));
    for area_idx = 1:length(area)
        x1 = manual_data.y1(area(area_idx));
        x2 = manual_data.y2(area(area_idx));
        patch([x1,x2,x2,x1],[reliability_thres,reliability_thres,1,1],[0,0,0],'FaceAlpha',0.15,'EdgeAlpha',0)
    end
    scatter(pro1_data.pos, pro1_data.reliability_best,5,'filled','MarkerEdgeAlpha', 0 , MarkerFaceAlpha=0.3,MarkerFaceColor='k')
    SI_now = [pro1_data.F_SI;pro1_data.B_SI;pro1_data.O_SI];
    [a,b]=max(SI_now);
    for category_now = 1:3
        all_idx = find(b==category_now & a>=0.2 & pro1_data.reliability_best>reliability_thres);
        if(~isempty(all_idx))
            scatter(pro1_data.pos(all_idx), pro1_data.reliability_best(all_idx),5,'filled','MarkerFaceColor',color_pool{category_now})
        end
    end
    ylabel('Reliability')
    xlim([min(pro1_data.pos), max(pro1_data.pos)])
    set_font

    subplot(4,2,[5,6]);
    x1 = manual_data.y1(area(1));x2 = manual_data.y2(area(1));
    unit_plot_here = find(pro1_data.reliability_best>reliability_thres & pro1_data.pos>x1 & pro1_data.pos<x2);
    all_data = pro1_data.response_best(unit_plot_here,1:1072);
    all_data = zscore(all_data, 0 ,2);
    mean_rsp = [];
    for cc = 1:max(Cluster_idx)
        data_here = all_data(:,find(Cluster_idx==clus_order(cc)));
        mean_rsp(cc,:) = mean(data_here');
    end
    imagesc(mean_rsp)
    clim([-1,1])
    colormap(give_me_orange_bao)
    box off
    xlabel('# Unit')
    ylabel('# Cluster')
    yticks([])
    set_font

    subplot(4,2,8)
    rr = pro1_data.reliability_best(unit_plot_here);
    pop_all_data = mean(all_data(:,1:1000));
    all_unit_type = pro1_data.UnitType(unit_plot_here);
    [~,a] = sort(pop_all_data,'descend');
    [~,b] = sort(pop_all_data,'ascend');
    big_img = [];
    for e = 1:3
        big_img = [big_img,[img_pool{a(e)};img_pool{a(e+4)}; 255*ones([40,227,3]); img_pool{b(e)};img_pool{b(e+4)}]];
    end
    imshow(big_img)
    imwrite(big_img,fullfile(root_dir,'Figs/F2',sprintf('%d.png',ses_now)))
    ylabel('Least       Most')
    set_font

    categoty_to_show = manual_data.Categoty(area(1));
    switch categoty_to_show{1}
        case 'F'
            SI_Data = pro1_data.F_SI;idx_here = face_idx;
        case 'B'
            SI_Data = pro1_data.B_SI;nm='Body';
        case 'O'
            SI_Data = pro1_data.O_SI;nm='Object';
    end

    subplot(4,2,7); hold on
    title(sprintf('%.1f%% d-prime>0.2', 100*sum(SI_Data(unit_plot_here)>0.2)./length(SI_Data(unit_plot_here))))
    histogram(SI_Data(unit_plot_here),20,'EdgeAlpha',0,'FaceColor',[0.5,0.5,0.5])
    xlabel(sprintf('%s selectivity', nm))
    ylabel('# Units')
    xlim([-max(abs(SI_Data)), max(abs(SI_Data))])
    xline(0.2,'LineStyle','--','LineWidth',1)
    filename_here(filename_here=='_')='-';
    sgtitle(filename_here(1:end-3))
    set_font


    position_lim = floor([min(pro1_data.pos),max(pro1_data.pos)]/10)*10;
    position_array = position_lim(1):25:position_lim(2);
    position_to_plot = position_array+200/2;
    m_array = [];
    for p_idx = 1:length(position_array)
        try
            p1 = position_array(p_idx);
            p2 = p1+200;
            nn_here = find(pro1_data.reliability_best>0.4 & pro1_data.pos>p1 & pro1_data.pos<p2);
            data_here = pro1_data.response_best(nn_here, :);
            pop_data = sum(data_here);
            zdata_here = zscore(data_here,0,2);
            m_array(p_idx,4) = CalcSI(pop_data(face_idx),pop_data(setdiff(1001:1072,face_idx)));
            m_array(p_idx,5) = CalcSI(pop_data(body_idx),pop_data(setdiff(1001:1072,body_idx)));
            m_array(p_idx,6) = CalcSI(pop_data(obj_idx),pop_data(setdiff(1001:1072,obj_idx)));
        catch
            m_array(p_idx,4) = 0;
            m_array(p_idx,5) = 0;
            m_array(p_idx,6) = 0;
        end
    end
    set_font
    subplot(4,2,[3,4]); hold on
    plot(position_to_plot,m_array(:,4),LineWidth=2,Color=color_pool{1})
    plot(position_to_plot,m_array(:,5),LineWidth=2,Color=color_pool{2})
    plot(position_to_plot,m_array(:,6),LineWidth=2,Color=color_pool{3})
    xlim([min(pro1_data.pos), max(pro1_data.pos)])
    ylim([-max(abs(m_array(:))),max(abs(m_array(:)))]*1.1)
    yline(0)
    ylabel('Population selectivity')
    xlabel('Depth (um)')
    set_font

    figsave(fullfile(root_dir,'Figs/F2'),sprintf('F2B_%02d', interested_ses(ses_now)))
    try
        unit_here = find(pro1_data.reliability_best>0.4 & pro1_data.UnitType==1);
        Selected_Raster = RasterData(unit_here,:,:);
        t1 = median(pro1_data.best_r_time1(unit_here));t2 = median(pro1_data.best_r_time2(unit_here));
        t_array = meta_data.global_params.pre_onset + [t1:t2];
        RasterSum = squeeze(sum(Selected_Raster(:,:,t_array),3));
        TrialIdx = meta_data.meta_data.trial_valid_idx(meta_data.meta_data.trial_valid_idx~=0);
        MTX = Organize_Trial_Data(RasterSum,TrialIdx);
        results = performgsn(MTX,struct('wantshrinkage',1));
        
        noise_rdm = results.cNb; signal_rdm = results.cSb;


        residual_response = []; signal_response=[]; trial_mean = squeeze(mean(MTX(:, :, :),3));
        for uu = 1:size(MTX,1)
            rsp_this_unit = squeeze(MTX(uu,:,:));
            mean_this_unit =  trial_mean(uu,:);
            residual=[];signal=[];
            for trial_index = 1:size(rsp_this_unit,2)
                res = rsp_this_unit(:,trial_index)-mean_this_unit';
                residual = [residual; res];
                signal = [signal; rsp_this_unit(:,trial_index)-res];
            end
            residual_response(uu,:) = residual;
            signal_response(uu,:) = signal;
        end

        noise_rdm_my = 1-squareform(pdist(residual_response,'correlation'));
        signal_rdm_my = 1-squareform(pdist(signal_response,'correlation'));

        empty_rdm = ones(size(noise_rdm));empty_rdm = tril(empty_rdm,-1);        

        nsc_here = corr(noise_rdm(empty_rdm==1),signal_rdm(empty_rdm==1),type='Pearson');
        noise_corr = [noise_corr; nsc_here];

        nsc_old_here = corr(noise_rdm_my(empty_rdm==1),signal_rdm_my(empty_rdm==1),'type','Pearson');
        noise_corr_old_ver = [noise_corr_old_ver; nsc_old_here];

        empty_rdm = eye(size(empty_rdm));
        noise_rdm(empty_rdm==1)=0;signal_rdm(empty_rdm==1)=0; % do not show horizontal line for better visualization
        figure(ses_now+200)
        set(gcf,'Position',[750 750 1000 200])
        subplot(1,4,1)
        imagesc(noise_rdm);colorbar
        cl = max(abs(noise_rdm(:)));
        clim([-cl,cl]);axis square;
        xticks([]);yticks([])
        title('Noise Correlation')
        subplot(1,4,2)
        imagesc(signal_rdm);colorbar
        cl = max(abs(signal_rdm(:)));
        clim([-cl,cl]);axis square;
        xticks([]);yticks([])
        colormap(give_me_orange_bao);
        title('Signal Correlation')
        sgtitle(sprintf('Session %d, corr = %.02f',ses_now,nsc_here))
    end
end


for ses_now = interested_ses
    figure(ses_now+200)
    subplot(1,4,3); hold on
    [bb,counts]=GoodHist(noise_corr,0:0.1:1,[0.25,0.25,0.25]);
    [a,b] = max(counts);
    scatter(median(noise_corr),a+3,'filled','LineWidth',2,MarkerEdgeAlpha=0,MarkerFaceColor=[0.25,0.25,0.25],Marker='v')
    xlabel('Correlation')
    ylabel('Session Counts')
    set_font

    subplot(1,4,4); hold on
    scatter(noise_corr,noise_corr_old_ver,10,MarkerEdgeColor='k');
    xlim([0,1]);ylim([0,1]);plot([0,1],[0,1],Color='k');xlabel('GSN'); ylabel('Original Method')
    set_font

    figsave('Figs\F2',sprintf('F2N_%02d',ses_now))
end