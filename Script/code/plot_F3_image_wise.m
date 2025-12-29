clearvars -except latency_save
interested_area={'AF3'};
cd C:\Users\moonl\Desktop\NNN
addpath(genpath(pwd))
manual_data = readtable("exclude_area.xls");
interested_time_point = 1:400;
H5_dir = 'C:\Users\moonl\Desktop\NNN\NNN_Data\Raw\H5FILES';
meta_example = load(fullfile(H5_dir,"ses10_240726_M5_3_info.mat"));
Proc_dir = 'C:\Users\moonl\Desktop\NNN\NNN_Data\Processed';
load img_pool.mat
load Clus_savee.mat
for aa = 1:length(interested_area)
    PSTH_CombinedData = [];
    ses_number = [];
    avg_rsp = [];
    for search_area = 1:height(manual_data)
        if(strcmp(manual_data.AREALABEL{search_area},interested_area{aa}))
            ThisSes_idx = manual_data.SesIdx(search_area);
            proc1_file_name = dir(fullfile(Proc_dir,sprintf('Processed_ses%02d*', ThisSes_idx)));
            proc1_file_name = proc1_file_name.name;
            pro1_data = load(fullfile(Proc_dir,proc1_file_name));
            cd(H5_dir)
            filename_here = dir(sprintf('ses%02d*h5',ThisSes_idx));
            filename_here = filename_here.name;
            metaname_here = dir(sprintf('ses%02d*mat',ThisSes_idx));
            metaname_here = metaname_here.name;
            meta_data = load(metaname_here);
            PSTHData = h5read(filename_here, '/response_matrix_img');
            cd ..

            x1=manual_data.y1(search_area);x2=manual_data.y2(search_area);
            good_neuron_idx = find(pro1_data.pos>x1 & pro1_data.pos<x2 & pro1_data.reliability_best>0.4);
            idx = clus_save{search_area}(find(clus_save{search_area}));
            interested_unit = 3;
            sprintf('Area %d, have %d neuron, Cluster2 n=%d ',search_area,length(good_neuron_idx),sum(idx==interested_unit))
            good_neuron_idx = good_neuron_idx(idx==interested_unit);
            data_here = PSTHData(good_neuron_idx,:, meta_example.global_params.pre_onset+interested_time_point);
            PSTH_CombinedData = [PSTH_CombinedData; data_here];
            ses_number = [ses_number, ones([1, length(good_neuron_idx)])*search_area];
            avg_rsp = [avg_rsp; pro1_data.response_best(pro1_data.reliability_best>0.4 & pro1_data.F_SI>0.2,1:1000)];
        end

    end
    
    for uu = 1:size(PSTH_CombinedData,1)
        dd = PSTH_CombinedData(uu,:,:);
        PSTH_CombinedData(uu,:,:) = dd./max(dd(:));
    end
    avg_rsp = zscore(avg_rsp,0,2);
    analysis_image = find(mean(avg_rsp)>0.7);
    latency = {};pop_rsp_mean={};
    all_this_ses = unique(ses_number);
    for ss = 1:length(all_this_ses)
        pop_rsp = squeeze(mean(PSTH_CombinedData(all_this_ses(ss)==ses_number, analysis_image, :),1));
        [~, latency{ss}] = max(pop_rsp(:,20:end)');
        latency{ss} = latency{ss}+20;
        latency_save{ss,interested_unit}=latency{ss};
        pop_rsp_mean{ss} = mean(avg_rsp(all_this_ses(ss)==ses_number,analysis_image));
    end

    ss = 1;
    width = find(pop_rsp_mean{ss}>2.3);

    [x,order]=sort(latency{ss}(width),'descend');
    all_img_here = length(width);
    cm_here = colormap_matplotlib('plasma');
    cm_here = flipud(cm_here(floor(linspace(40,200,all_img_here)),:));
    img = [];
    figure;
    set(gcf,'Position',[500 600 1000 450])
    subplot(2,3,1); hold on
    illustration_points = [length(order):-3:1,1];
    scatter(latency{ss},pop_rsp_mean{ss},8,MarkerFaceAlpha=0,MarkerEdgeColor=[0,0,0])

    interested_images = [47    68    98   116   147];
    for ii = illustration_points
        scatter(latency{ss}(width(order(ii))),pop_rsp_mean{ss}(width(order(ii))),20,'filled','MarkerFaceColor',cm_here(ii,:),MarkerEdgeAlpha=0);
    end
    xlabel('Peak Latency (ms)')
    ylabel('Firing Rate (a.u.)')

    subplot(2,3,2)
    hold on
    for ii = illustration_points
        m = mean(squeeze(PSTH_CombinedData(all_this_ses(ss)==ses_number, analysis_image(width(order(ii))), :)));
        e = std(squeeze(PSTH_CombinedData(all_this_ses(ss)==ses_number, analysis_image(width(order(ii))), :)))./sqrt(sum(all_this_ses(ss)==ses_number));
        shadedErrorBar(interested_time_point,m,e,'lineprops',{'Color',cm_here(ii,:)},'PatchSaturation',0.5);
        img = [img, add_edge(img_pool{analysis_image(width(order(ii)))}, floor(cm_here(ii,:).*255),20)];
    end
    xlim([50,250])
    xlabel('Time (ms)')
    ylabel('Norm. firing rate')
    title('Cluster 2')

    subplot(2,3,4)
    hold on
    scatter(latency_save{1,2},latency_save{2,2},8,'k');
    tmp_data = [latency_save{1,2},latency_save{2,2}];
    xl_here = [min(tmp_data)-5,max(tmp_data)+5];
    axis equal;xlim(xl_here);ylim(xl_here)
    plot(xl_here,xl_here,Color=[0.5,0.5,0.5],LineWidth=1,LineStyle="--")
    text(180,120,sprintf('corr %.02f', corr(latency_save{1,2}',latency_save{2,2}',type="Spearman")))
    xlabel('Latency Clus2 Ses1 (ms)')
    ylabel('Latency Clus2 Ses2 (ms)')

    subplot(2,3,5);  hold on
    scatter(latency_save{1,1},latency_save{1,2},8,'k');
    tmp_data = [latency_save{2,1},latency_save{2,2}];
    xl_here = [min(tmp_data)-5,max(tmp_data)+5];
    axis equal;xlim(xl_here);ylim(xl_here)
    plot(xl_here,xl_here,Color=[0.5,0.5,0.5],LineWidth=1,LineStyle="--")
    text(180,120,sprintf('corr %.02f', corr(latency_save{1,1}',latency_save{1,2}',type="Spearman")))
    xlabel('Latency Clus1 Ses1 (ms)')
    ylabel('Latency Clus2 Ses1 (ms)')

    load alexnet_resp.mat
    rsp_here = fc6_4096(:,analysis_image);
    [coeff,pc_score,latent] = pca(rsp_here',NumComponents=5);
    for uu = 1:length(analysis_image)
        train_set = setdiff(1:length(analysis_image), uu);
        x = lscov([pc_score(train_set,:), ones([size(pc_score,1)-1,1])],latency{2}(train_set)');
        predicted_latency(uu) = [pc_score(uu,:),1]*x;
    end

%     scatter(latency{1},predicted_latency,8,'k');
%     axis equal;xlim(xl_here);ylim(xl_here)
%     plot(xl_here,xl_here,Color=[0.5,0.5,0.5],LineWidth=1,LineStyle="--")
%     text(180,120,sprintf('corr %.02f', corr(latency{1}',predicted_latency',type="Spearman")))
%     xlabel('Latency session1 (ms)')
%     ylabel('Latency predicted by alexnet (ms)')

    pc_score ={};
    rsp_here = cv5_data(:,analysis_image);
    [coeff,pc_score{1},latent] = pca(rsp_here',NumComponents=5);
    rsp_here = fc6_4096(:,analysis_image);
    [coeff,pc_score{2},latent] = pca(rsp_here',NumComponents=5);
    rsp_here = fc7_4096(:,analysis_image);
    [coeff,pc_score{3},latent] = pca(rsp_here',NumComponents=5);
    rsp_here = fc8_4096(:,analysis_image);
    [coeff,pc_score{4},latent] = pca(rsp_here',NumComponents=5);
    rsp_here = softmax_1000(:,analysis_image);
    [coeff,pc_score{5},latent] = pca(rsp_here',NumComponents=5);

    data_here = [];
    for cc = 1:5
        for uu = 1:length(analysis_image)
            train_set = setdiff(1:length(analysis_image), uu);
            x = lscov([pc_score{cc}(train_set,:), ones([size(pc_score{cc},1)-1,1])],latency{2}(train_set)');
            predicted_latency(uu) = [pc_score{cc}(uu,:),1]*x;
        end
        data_here(cc) = corr(latency{1}',predicted_latency',type="Spearman");
    end
    subplot(2,3,6)
    bar(1:5, data_here,'EdgeAlpha',0)
    xticklabels({'conv5','fc6','fc7','fc8','output'})
    set(gca,"XTickLabelRotation",45)
    xlabel('layer')
    ylabel('Predict Accuracy (r)')
    xlim([-1,7])
    set_font
    subplot(2,3,3)
    imshow(img)
    figsave('Figs\F3A',sprintf('F3E'))
end
% sub

%%
data_here = [latency_save{1,1};latency_save{1,2};latency_save{2,1};latency_save{2,2}];
figure
% scatter(latency_save{,})
imagesc(1-squareform(pdist(data_here,'correlation')))