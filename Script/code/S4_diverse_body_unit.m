%% Second Part of Fig.4
%% Play with ROI
clear
manual_data = readtable("exclude_area.xls");
root_dir = 'C:\Users\moonl\Desktop\NNN';
[proc_dir,raw_dir] = gen_dirs(root_dir);
load ROI_data.mat
load img_pool.mat

order = randperm(12);
face_idx = 1000+[1:5,9:15,6,7,8,16:24];
face_idx = face_idx([order,order+12]);
body_idx = 1000+[26:31,43:48,50:61];
body_idx = body_idx([order,order+12]);
obj_idx = setdiff(1025:1072, body_idx);
obj_idx = obj_idx([order,order+12]);
loc_idx = [face_idx, body_idx, obj_idx];
HHs = ['L','R'];

h_here = 1; % corresponding to Left
ROI_Here = [1,2]; % corresponding to EBA and FBA
BIG_x = [];
ROI_LGD = {};
for hh = h_here
    for ROI_idx = ROI_Here
        roi_data = [];
        for ss = [1,2,5,7]
            brain_data = ROI_data{ss,ROI_idx,hh};
            brain_data = zscore(brain_data,0,2);
            roi_data = [roi_data; mean(ROI_data{ss,ROI_idx,hh})];
        end
        roi_data = mean(roi_data,'omitnan');
        BIG_x = [BIG_x;roi_data];
        ROI_LGD{end+1} = sprintf('%s%s',HHs(hh),all_interested_roi{ROI_idx});
    end
end
%
interested_area =  {'MB3','AB3'};
figure
set(gcf,'Position',[250 350 1050 360])
for aa = 1:length(interested_area)
    plot_idx = aa;
    all_rsp = [];
    avg_rsp = [];
    F_si = [];
    pos = [];
    ses_id = [];
    fob_rsp = [];
    for search_area = 1:height(manual_data)
        if(strcmp(manual_data.AREALABEL{search_area},interested_area{aa}))
            ThisSes_idx = manual_data.SesIdx(search_area);
            proc1_file_name = dir(fullfile(proc_dir,sprintf('Processed_ses%02d*', ThisSes_idx)));
            proc1_file_name = proc1_file_name.name;
            pro1_data = load(fullfile(proc_dir,proc1_file_name));
            filename_here = dir(fullfile(raw_dir,sprintf('ses%02d*h5',ThisSes_idx)));
            filename_here = filename_here.name;
            metaname_here = dir(fullfile(raw_dir,sprintf('ses%02d*mat',ThisSes_idx)));
            metaname_here = metaname_here.name;
            meta_data = load(metaname_here);
            x1=manual_data.y1(search_area);x2=manual_data.y2(search_area);
            switch interested_area{aa}(2)
                case 'F'
                    SI = pro1_data.F_SI;
                case 'B'
                    SI = pro1_data.B_SI;
            end
            good_neuron_idx = find(pro1_data.pos>x1 & pro1_data.pos<x2 & pro1_data.reliability_best>0.4 & SI>0.2);

            avg_rsp = [avg_rsp; pro1_data.response_best(good_neuron_idx,1:1000)];
            F_si = [F_si, pro1_data.F_SI(good_neuron_idx)];
            pos = [pos, pro1_data.pos(good_neuron_idx)];
            ses_id = [ses_id, search_area*ones([1,length(good_neuron_idx)])];
            fob_rsp = [fob_rsp;zscore(pro1_data.response_best(good_neuron_idx,loc_idx),0,2)];
        end
    end

    avg_rsp = zscore(avg_rsp,0,2);

    pr1 = [];pr2 = []; pr = []; pr_p = [];
    for uu = 1:size(avg_rsp,1)
        % first 500 stim
        data_combined = [avg_rsp(uu,1:500)',BIG_x(:,1:500)'];
        [partial_corr_matrsix, p_values] = partialcorr(data_combined,'Type','Pearson');
        pr1(uu,:)=partial_corr_matrsix(1,2:end);
        % second 500 stim
        data_combined = [avg_rsp(uu,501:1000)',BIG_x(:,501:1000)'];
        [partial_corr_matrsix, p_values] = partialcorr(data_combined,'Type','Pearson');
        pr2(uu,:)=partial_corr_matrsix(1,2:end);
        % 1000 stim

        data_combined = [avg_rsp(uu,:)',BIG_x(:,:)'];
        [partial_corr_matrsix, p_values] = partialcorr(data_combined,'Type','Pearson');
        pr(uu,:) = partial_corr_matrsix(1,2:end);
        pr_p(uu,:) = p_values(1,2:end);
        uu
    end
    [a,b]=sort(pr1(:,1)-pr1(:,2));
    [c,d]=sort(pr2(:,1)-pr2(:,2));
    subplot(2,5,1+5*(plot_idx-1)); hold on
    cm = colormap_matplotlib('Set1',8);
    scatter(1:length(a), pr2(b,1),5,'filled',MarkerEdgeAlpha=0,MarkerFaceColor=cm(1,:),MarkerFaceAlpha=0.3)
    scatter(1:length(a), pr2(b,2),5,'filled',MarkerEdgeAlpha=0,MarkerFaceColor=cm(2,:),MarkerFaceAlpha=0.3)
    
    [r1,p]=corr([1:length(b)]', pr2(b,1),'type','Spearman');
    [r2,p]=corr([1:length(b)]', pr2(b,2),'type','Spearman');
    title(sprintf('corr=%.02f, corr=%.02f',r1,r2))
    
    edges = 1:length(a);
    for ee = 1:length(edges)
        data1 = pr2(b(max(edges(ee)-25,1):min(edges(ee)+25,length(a))) ,1);
        bin_med1(ee) = mean(data1);
        se_1(ee) = std(data1)./sqrt(length(data1));
    
        data2 = pr2(b(max(edges(ee)-25,1):min(edges(ee)+25,length(a))) ,2);
        bin_med2(ee) = mean(data2);
        se_2(ee) = std(data2)./sqrt(length(data2));
    end
    errorbar(edges(1:25:end), bin_med1(1:25:end), se_1(1:25:end),'LineWidth',1.5, 'CapSize',0,Color=cm(1,:))
    errorbar(edges(1:25:end), bin_med2(1:25:end), se_2(1:25:end),'LineWidth',1.5, 'CapSize',0,Color=cm(2,:))
    title(interested_area{aa})
    ylim([-0.2, 0.32])
    xlim([-1,length(a)+1])
    xticks([])
    xlabel('Relative rank (cross-validated)')
    ylabel('PCC (r)')
    clear bin_med1 bin_med2 se_1 se_2


    subplot(2,5,2+5*(plot_idx-1)); hold on
    scatter(pr(:,1),pr(:,2),5,'filled',MarkerFaceColor=[0.5,0.5,0.5],MarkerEdgeAlpha=0);
    
    eba_neuron = find(pr(:,1)>0 & pr_p(:,1)<0.01 & pr_p(:,2)>0.01);
    scatter(pr(eba_neuron,1),pr(eba_neuron,2),5,'filled',MarkerEdgeAlpha=0,MarkerFaceColor=cm(3,:));
    
    fba_neuron = find(pr(:,2)>0 & pr_p(:,2)<0.01 & pr_p(:,1)>0.01);
    scatter(pr(fba_neuron,1),pr(fba_neuron,2),5,'filled',MarkerEdgeAlpha=0,MarkerFaceColor=cm(4,:));

    xlim([-0.15, 0.25]);ylim([-0.15, 0.25])
    colorbar off
    xline(0,'LineStyle','--','LineWidth',1);yline(0,'LineStyle','--','LineWidth',1)
    xlabel(sprintf('PCC %s',ROI_LGD{1}));ylabel(sprintf('PCC %s',ROI_LGD{2}))
    legend({'',sprintf('EBA %.0f%%', 100*length(eba_neuron)/size(pr,1)),sprintf('FBA %.0f%%', 100*length(fba_neuron)/size(pr,1)) },'Box','off')
    
    subplot(2,5,3+5*(plot_idx-1)); hold on
    body_img = find(mean(avg_rsp)>0.5);
    corr(mean(avg_rsp(eba_neuron,:))',mean(avg_rsp(fba_neuron,:))')
    corr(mean(avg_rsp(eba_neuron,body_img))',mean(avg_rsp(fba_neuron,body_img))')
    eba_neuron_rsp{plot_idx} = avg_rsp(eba_neuron,:);
    fba_neuron_rsp{plot_idx} = avg_rsp(fba_neuron,:);
    
    a_here = mean(avg_rsp(eba_neuron,:));
    b_here = mean(avg_rsp(fba_neuron,:));
    scatter(a_here,b_here,6,'filled',MarkerEdgeAlpha=0,MarkerFaceColor=[0.5,0.5,0.5],MarkerFaceAlpha=0.5)
    p = polyfit(a_here, b_here, 1);
    y_fit = polyval(p, a_here);
    plot(a_here, y_fit, 'k','LineWidth',1);
    text(0.8,-0.5,sprintf('corr=%.02f', corr(a_here',b_here','type','Spearman')),Color=[0,0,0],FontWeight="bold")
    
    a_here = mean(avg_rsp(eba_neuron,body_img));
    b_here = mean(avg_rsp(fba_neuron,body_img));
    scatter(a_here,b_here,8,'filled',MarkerEdgeAlpha=0,MarkerFaceColor=[1,0,0],MarkerFaceAlpha=0.5)
    p = polyfit(a_here, b_here, 1);
    y_fit = polyval(p, a_here);
    plot(a_here, y_fit, 'r','LineWidth',1);
    text(0.8,-0.8,sprintf('corr=%.02f', corr(a_here',b_here','type','Spearman')),Color=[1,0,0],FontWeight="bold")
    xlabel('EBA units response (a.u.)')
    ylabel('FBA units response (a.u.)')
    xlim([-1,2]);ylim([-1,2])
    set_font

    switch interested_area{aa}(1)
        case 'M'
            subplot(2,5,4);
            hold on
            histogram(pos(intersect(eba_neuron,find(ses_id==ses_id(2)))),20,EdgeAlpha=0,Normalization="pdf",FaceColor=cm(3,:));
            histogram(pos(intersect(fba_neuron,find(ses_id==ses_id(2)))),20,EdgeAlpha=0,Normalization="pdf",FaceColor=cm(4,:));
        case 'A'
            subplot(2,5,9);
            hold on
            histogram(pos(intersect(eba_neuron,find(ses_id==56))),20,EdgeAlpha=0,Normalization="pdf",FaceColor=cm(3,:));
            histogram(pos(intersect(fba_neuron,find(ses_id==56))),20,EdgeAlpha=0,Normalization="pdf",FaceColor=cm(4,:));
    end
    xlabel('Depth (um)')
    ylabel('Frequency')
    legend({'EBA unit','FBA unit'},'Box','off',Location='best')
    set_font
end
%
subplot(2,5,[5,10])
a = mean([eba_neuron_rsp{1}]);
[~,a] = sort(a,'descend');
b = mean([fba_neuron_rsp{1}]);
[~,b] = sort(b,'descend');
img=[];
for ee = 1:5
    img = [img; [img_pool{a(ee)},ones([227,100,3])*255,img_pool{b(ee)}]];
end
imshow(img)
title('EBA unit       FBA unit')

figsave(fullfile(root_dir,'Figs/F4/'), 'F4_para')