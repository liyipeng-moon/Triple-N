clear
root_dir = 'C:\Users\moonl\Desktop\NNN';
cd(root_dir)
addpath(genpath(pwd));
[proc_dir,raw_dir] = gen_dirs(root_dir);
reliability_thres = 0.4;
manual_data = readtable("exclude_area.xls");
area_array=[];
r_array = [];
%
resolution_array = [200:200:1600];
rsp_pool={};
for k = 1:length(resolution_array)+1
    rsp_pool{k}=[];
end
    for area_now = 1:102
    if(~strcmp(manual_data.Area{area_now},'IT'))
        continue
    end
    ses_idx = manual_data.SesIdx(area_now);
    proc1_file_name = dir(fullfile(proc_dir,sprintf('Processed_ses%02d*', ses_idx)));
    proc1_file_name = proc1_file_name.name;
    pro1_data = load(fullfile(proc_dir,proc1_file_name));
    x1 = manual_data.y1(area_now);
    x2 = manual_data.y2(area_now);
    pos_now = pro1_data.pos;
    unit_here = find(pos_now>x1 & pos_now<x2 & pro1_data.reliability_best>reliability_thres);
    pos_now = pos_now(unit_here);
    [~, unit_order] = sort(pos_now,'ascend');
    pos_now = pos_now(unit_order);
    rsp_now = pro1_data.response_best(unit_here(unit_order),1:1000);

    for resolution_idx = 1:length(resolution_array)
        re_here = resolution_array(resolution_idx);
        bound_array = floor(min(pos_now):re_here:max(pos_now));
        for bound_idx = 1:length(bound_array)-1
            start_with = bound_array(bound_idx);
            unit_within_bound = find(pos_now>start_with & pos_now<start_with+re_here);
            fprintf('resolution: %d um, %d unit for %d bound \n',re_here, length(unit_within_bound),length(bound_array))
            rsp_pool{resolution_idx} = [rsp_pool{resolution_idx}; sum(rsp_now(unit_within_bound,:),1)];
        end
    end
    sprintf('%d', area_now)
end
%
[alls,nm] = load_embedding;

s{1} = alls{6};
s{2} = alls{14};

%
nm={};
for resolution_idx = 1:length(resolution_array)
    nm{resolution_idx} = sprintf('%d um',resolution_array(resolution_idx));
end
%
close all
figure
set(gcf,'Position',[50 550 800 450])
start_parfor
rng(1009)
for level = 1:8
    rsp_matrix = rsp_pool{level};
    pred_r_array = [];
    pred_r2_array = [];
    parfor unit = 1:size(rsp_matrix,1)
        rsp_now = zscore(double(rsp_matrix(unit,:)));
        [pred_r_array(unit,:),pred_r2_array(unit,:),best_numbers] = compare_encoders(s, rsp_now, 10);
    end
    
    visual_best = pred_r_array(:,2);
    LLM = pred_r_array(:,1);
    good_neuron = find(~isnan(visual_best+LLM));
    visual_best = visual_best(good_neuron);
    LLM = LLM(good_neuron);
    [beta_estimate(level),beta_CI(level,:)] = demingRegression(visual_best, LLM);

    if(level==1 || level==8)
        nexttile;hold on
        bin_scatter(visual_best,LLM,0:0.05:1,0:0.05:1);
        plot([0,1],[0,beta_estimate(level)],'LineWidth',1,Color='k')
        xlabel('Visual performance (r)')
        ylabel('Language performance (r)')
        title(sprintf('Bin size: %s, \\sigma = 1\n LVR=%.02f',nm{level},beta_estimate(level)))
        xtickangle(0)
    end
    drawnow
end

%% 
SE = (beta_CI(:,2)-beta_CI(:,1))/2;
nexttile; hold on
errorbar(resolution_array, beta_estimate, SE','LineWidth',1,Color=[0.5,0.5,0.5]);
ylim([0.5,1]);
scatter(0,0.74,"filled","MarkerEdgeAlpha",0,'MarkerFaceColor',[0,0,0])
xticks(resolution_array)
xlim([-100,resolution_array(end)+100])
xlabel('Spatial bin (mm)')
ylabel('LVR')
legend({'Binned MUA','Unbinned Unit'},Box="off",Location="best")
colorbar
xtickangle(90)
set_font


%%
clear beta_estimate beta_CI
rng(1009)
sigmalevel = 0.1:0.1:0.9;

nm={};

for s_idx = 1:length(sigmalevel)
    nm{s_idx} = sprintf('\\sigma = %.01f',sigmalevel(s_idx));
end
% do compress for 200 um
rsp_matrix = rsp_pool{1};
for level = 1:length(sigmalevel)
    pred_r_array = [];
    pred_r2_array = [];
    parfor unit = 1:size(rsp_matrix,1)
        rsp_now = zscore(double(rsp_matrix(unit,:) .^ sigmalevel(level))) ;
        [pred_r_array(unit,:),pred_r2_array(unit,:),best_numbers] = compare_encoders(s, rsp_now, 10);
    end
    visual_best = pred_r_array(:,2);
    LLM = pred_r_array(:,1);
    good_neuron = find(~isnan(visual_best+LLM));
    visual_best = visual_best(good_neuron);
    LLM = LLM(good_neuron);
    [beta_estimate(level),beta_CI(level,:)] = demingRegression(visual_best, LLM);
    if(level==1 || level==9)
        nexttile;hold on
        bin_scatter(visual_best,LLM,0:0.05:1,0:0.05:1);
        plot([0,1],[0,beta_estimate(level)],'LineWidth',1,Color='k')
        xlabel('Visual performance (r)')
        ylabel('Language performance (r)')
        title(sprintf('Bin size: 200um, %s\n LVR=%.02f',nm{level},beta_estimate(level)))
        xtickangle(0)
    end
end


SE = (beta_CI(:,2)-beta_CI(:,1))/2;
nexttile; hold on
errorbar(sigmalevel, beta_estimate, SE','LineWidth',1,Color=[0.5,0.5,0.5]);
ylim([0.5,1]);
scatter(1,beta_estimate(end),"filled","MarkerEdgeAlpha",0,'MarkerFaceColor',[0,0,0])
xticks(0.1:0.2:0.9)
xlim([0,1])
xlabel('sigma')
ylabel('LVR')
legend({'compressed MUA','MUA'},Box="off",Location="best")
colorbar
set_font

figsave(fullfile(root_dir,"Figs/F5"),'FS5')
