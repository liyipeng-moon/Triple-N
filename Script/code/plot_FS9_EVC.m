clear;clc
root_dir = 'C:\Users\moonl\Desktop\NNN';
cd(root_dir)
addpath(genpath(pwd));
[proc_dir,raw_dir] = gen_dirs(root_dir);
%%
img = load("img_pool.mat").img_pool;
Contrast_Map = zeros([1000,50,50]);
Luminance_Map = zeros([1000,50,50]);
parfor i = 1:1000
    [Contrast_Map(i,:,:), Luminance_Map(i,:,:)] = localContrastMap(img{i});
end
Contrast_Map = Contrast_Map(:, 2:end-1, 2:end-1);
Luminance_Map = Luminance_Map(:,2:end-1,2:end-1);
%%
close all

for interested_ses = [78, 79]
    clear rf_all params fit_result gof
    proc_file = dir(fullfile(proc_dir,sprintf('Processed_ses%02d*',interested_ses)));
    load(fullfile(proc_dir,proc_file.name))
    rf_all = zeros([length(F_SI),48,48]);
    unit_n = length(F_SI);
    tic
    parfor unit = 1:unit_n
        rsp = response_best(unit,1:1000);
        random_data = zeros(48);
        for i = 1:1000
            rf_all(unit,:,:) = rf_all(unit,:,:) + rsp(i) .* Contrast_Map(i,:,:);
        end

        boot_time = 1000;
        for i = 1:1000
            for t = 1:boot_time
                random_data = random_data +  rsp(randi(1000)) .* squeeze(Contrast_Map(i,:,:));
            end
        end
        rf_all(unit,:,:) = squeeze(rf_all(unit,:,:)) - (random_data/boot_time);

        RF_Here = squeeze(rf_all(unit,:,:));
        [params{unit}, fit_result{unit}, gof(unit)] = gaussian2D_fit_symmetric(RF_Here);
    end
    toc
    %
    figure;
    sgtitle(proc_file.name,Interpreter="none")
    set(gcf,'Position',[100 25 660 250])
    [val,b] = sort(gof,'descend');

    subplot(1,2,1); hold on
    unit = b(3);
    
    show_rf(fit_result{unit},squeeze(rf_all(unit,:,:)))
    title('Example Unit')
    set(gca, 'YDir', 'reverse');

    subplot(1,2,2);hold on
    for example = 1:length(F_SI)
        if(gof(example)>0.45)
            plot_circle(fit_result{example},0)
        end
    end
    plot_circle(fit_result{unit},1)
    axis equal tight
    xlim([-5.5,5.5]); ylim([-5.5,5.5]); xline(0); yline(0);
    title('Population')
    set(gca, 'YDir', 'reverse'); % note that in matlab, larger y value corresponding to lower field, and images goes with y axis
    drawnow
    figsave(fullfile(root_dir,"Figs/EVC/"),sprintf('%d',interested_ses))

    
end

%% 
clear;clc
root_dir = 'C:\Users\moonl\Desktop\NNN';
cd(root_dir)
addpath(genpath(pwd));
[proc_dir,raw_dir] = gen_dirs(root_dir);
manual_data = readtable("exclude_area.xls");

[s,s_name] = load_embedding;
alex_array = 9:16;
xlabels={'Conv1','Conv2','Conv3','Conv4','Conv5','FC6','FC7','FC8'};
ylabels={'V1','V2','V4','IT'};
interesed_case = {'V1','V2','V4'};

BigR=[];
for case_idx = 1:length(interesed_case)
    Latency{case_idx} = []; R_array = []; best_layer{case_idx} = [];
    for area_idx = 1:102
        if(strcmp(manual_data.AREALABEL{area_idx}(1:length(interesed_case{case_idx})),interesed_case{case_idx}))
            proc_file = dir(fullfile(proc_dir,sprintf('Processed_ses%02d*',manual_data.SesIdx(area_idx))));
            proc_data = load(fullfile(proc_dir,proc_file.name));
            unit_here = find(proc_data.pos>manual_data.y1(area_idx) & proc_data.pos<manual_data.y2(area_idx) & proc_data.reliability_best>0.4);
            unit_here = find(proc_data.reliability_best>0.4);
            [~, lat] = max(proc_data.mean_psth(unit_here,51:350)');
            Latency{case_idx} = [Latency{case_idx}, lat];
            encoder_data = load(sprintf('s5_encode_ses_%03d.mat', area_idx));
            reliability_here = load(sprintf('s5_encode_ses_%03d.mat', area_idx)).r_here;
            for unit = 1:length(reliability_here)
                r = encoder_data.pred_r_array(unit,:);
                r2 = encoder_data.pred_r2_array(unit,:);
                if(max(r2)>0)
                    tmp = r(alex_array);
                    tmp = tmp-min(tmp);
                    tmp = tmp/max(tmp);
                    R_array = [R_array;tmp];
%                     R_array = [R_array; r(alex_array)./max(r(alex_array))];
                    [~,best_layer{case_idx}(end+1)] = max(r(alex_array));
                end
            end
        end
    end
    BigR = [BigR; mean(R_array)];
end

% Load for IT
case_idx=4;
Latency{case_idx} = []; R_array = []; best_layer{case_idx} = [];
for area_idx = 1:102
    if(strcmp(manual_data.Area{area_idx},'IT'))
        proc_file = dir(fullfile(proc_dir,sprintf('Processed_ses%02d*',manual_data.SesIdx(area_idx))));
        proc_data = load(fullfile(proc_dir,proc_file.name));
        unit_here = find(proc_data.pos>manual_data.y1(area_idx) & proc_data.pos<manual_data.y2(area_idx) & proc_data.reliability_best>0.4);
        [~, lat] = max(proc_data.mean_psth(unit_here,51:350)');
        Latency{case_idx} = [Latency{case_idx}, lat];

        encoder_data = load(sprintf('s5_encode_ses_%03d.mat', area_idx));
        reliability_here = load(sprintf('s5_encode_ses_%03d.mat', area_idx)).r_here;
        for unit = 1:length(reliability_here)
            r = encoder_data.pred_r_array(unit,:);
            r2 = encoder_data.pred_r2_array(unit,:);
            if(max(r2)>0)
                tmp = r(alex_array);
                tmp = tmp-min(tmp);
                tmp = tmp/max(tmp);
                R_array = [R_array;tmp];
%                 R_array = [R_array; r(alex_array)./max(r(alex_array))];
                [~,best_layer{case_idx}(end+1)] = max(r(alex_array));
            end
        end
    end
end

BigR = [BigR; mean(R_array)];

for k = 1:size(BigR,1)
    tmp = BigR(k,:);
    tmp = tmp-min(tmp);
    tmp = tmp/max(tmp);
    BigR(k,:)=tmp;
end
%
cm_here = colormap_matplotlib('Set3',4);
figure; set(gcf,'Position',[200 250 580 400])
subplot(4,3,[1,2,4,5])
hold on
data = []; group = [];
for i = 1:numel(Latency)
    data = [data, Latency{i}];
    length(Latency{i})
    group = [group, i * ones(1, numel(Latency{i}))];
end
boxplot(data, group,'Symbol','','Whisker',0.5,'orientation','horizontal');
pretty_box(cm_here)
title('Latency (ms)')
xlim([25,225]);yticks(1:4);yticklabels(ylabels);ylim([0.5,4.5]);
set_font
%
subplot(4,3,[7,8,10,11])
imagesc(BigR)
yticks(1:4); xticks(1:9); xticklabels(xlabels); yticklabels(ylabels);xtickangle(90);clim([0,1])
title('Normalized Accuracy')
colormap(colormap_matplotlib('plasma'))
axis tight
set_font
%
for a = 1:4
    subplot(4,3,3*a); hold on
    for l = 1:8
        p(l) = sum(best_layer{a}==l)./length(best_layer{a});
    end
    bar(1:8,p,'EdgeAlpha',0,FaceColor=[0.5,0.5,0.5])
    ylabel('Unit %');title(ylabels{a});
    xticks([]);
    set_font
end
xticks([1:8]);
xlabel('Layers')
figsave(fullfile(root_dir,'Figs/EVC/'),'ComparingEVC')