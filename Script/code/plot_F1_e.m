clear
root_dir = 'C:\Users\moonl\Desktop\NNN';
cd(root_dir)
addpath(genpath(pwd));
[proc_dir,raw_dir] = gen_dirs(root_dir);
mkdir Figs\F1
interested_ses = 1;
all_data = dir(fullfile(proc_dir,sprintf("Processed_ses%02d_*", interested_ses)));
load(all_data.name)
all_data = dir(fullfile(raw_dir,sprintf("ses%02d_*info*", interested_ses)));
load(all_data.name)
trial_data = meta_data.trial_valid_idx;
onset_time = meta_data.onset_time_ms;
load img_pool.mat;
% Play with 10 trials, in 
step = 8;
for a_start = 4008:7:5000
    a = a_start;
    b=a+step;
    time_start = onset_time(a);
    time_end = onset_time(b)+300;
    img_idx = trial_data(a:b);

    if(max(img_idx)<1000 && min(img_idx)>0 && max(diff(onset_time(a:b)))<350)
        img = [];onset_bar = [];

        for ii = 1:length(img_idx)
            img = [img, img_pool{img_idx(ii)}];
            img = [img, 128*ones(size(img_pool{1}))];
            onset_bar = [onset_bar, 0.8*ones([50,227,3])];
            onset_bar = [onset_bar, 0.5*ones([50,227,3])];
        end
    else
        continue
    end


    figure;
    set(gcf,'Position',[750 200 900 700])


    subplot(5,1,1)
    illus_dot_size = 6;
    for img_idx = 114:227:size(img,2)
        for xx = -illus_dot_size:illus_dot_size
            for yy = -illus_dot_size:illus_dot_size
                img(114+xx,img_idx+yy,:)=[255,0,0];
            end
        end
    end
    imshow([img])
    subplot(5,1,[2,3])
    hold on
    set(gca,'TickDir','none')
    LOC = 0;
    good_units = find(reliability_best>0.4 & B_SI>0.2);
    cm_here = colormap_matplotlib('plasma');
    cm_here = cm_here(1+floor(200*(1:length(good_units))/length(good_units)),:);
    data_all = zeros([1, time_end-time_start]);
    for units = 1:length(good_units)
        unit_here = good_units(units);
        raster_raw = GoodUnitStrc(unit_here).spiketime_ms;
        raster_time = raster_raw(raster_raw>time_start &raster_raw<time_end);
        scatter(raster_time, repmat(LOC, [1,length(raster_time)]),4,'filled',Marker='square',MarkerFaceColor=cm_here(units,:),MarkerEdgeColor=cm_here(units,:));
        LOC = LOC+1;
        time_here = floor(raster_time-time_start)+1;
        d0 = zeros([1, time_end-time_start]);
        d0(time_here) = d0(time_here)+1;
        data_all = data_all+conv(d0,ones([1,30]),'same');
    end

    xline(onset_time(a+1:b),Color=[0.5,0.5,0.5],LineStyle="--",LineWidth=1)
    % xline((onset_time(a:b)+150),Color=[0.5,0.5,0.5],LineStyle="--")
    ax=gca;

    xlim([time_start,time_end])
    ylim([-2, LOC+2]);
    xticks([])
    ylabel('# Units')
    set_font

    subplot(5,1,5)
    imshow(onset_bar)

    subplot(5,1,4)

    plot(time_start:time_end-1,data_all./length(good_units),'LineWidth',1,Color=[0.2,0.2,0.2])
    xlim([time_start,time_end])
    xticks(time_start:300:time_end);
    xticklabels(0:300:step*300+300)
    box off
    set(gca,'TickDir','none')
    xline(onset_time(a+1:b),Color=[0.5,0.5,0.5],LineStyle="--",LineWidth=1)
    ylabel('% Unit recruited')
    ylim([0, 0.8])
    yticks([0, 0.8])
    xlabel('Time (ms)')

    set_font
    if(max(data_all)/length(good_units)>0.3)
        break
    end
    close all
end

set(gcf,'renderer','painters');
figsave(fullfile(root_dir,"Figs/F1"),sprintf('F1'))
return
% %%
% figure
% best_units = find(reliability_best==max(reliability_best) & UnitType<3);
% [rsp_here,rsp_order] = sort(response_best(best_units,:),'descend');
% i1=2;i2=940;
% trials = [find(trial_data==rsp_order(i1));  find(trial_data==rsp_order(i2))];
% img_here = [add_edge(img_pool{rsp_order(i1)},[155,0,0],10),add_edge(img_pool{rsp_order(i2)},[0,0,155],10)];
% nexttile
% imshow(img_here)
% all_onset_time = onset_time(trials);
% raster_raw = GoodUnitStrc(best_units).spiketime_ms;
% nexttile; hold on
% LOC = 0;
% trials = trials(:);
% trials = trials(1:4:end);
% for tt = 1:4
%     time_start = onset_time(tt);
%     time_end = time_start+300;
%     LOC = 0;
%     for units = 1:length(good_units)
%         unit_here = good_units(units);
%         raster_raw = GoodUnitStrc(unit_here).spiketime_ms;
%         raster_time = raster_raw(raster_raw>time_start & raster_raw<time_end);
%         raster_time = raster_time-time_start;
%         scatter(raster_time+(tt*350), repmat(LOC, [1,length(raster_time)]),12,'filled',Marker='square',MarkerFaceColor=cm_here(units,:),MarkerEdgeAlpha=0);
%         LOC = LOC+1;
% %         time_here = floor(raster_time-time_start)+1;
% %         d0 = zeros([1, time_end-time_start]);
% %         d0(time_here) = d0(time_here)+1;
% %         data_all = data_all+conv(d0,ones([1,20]),'same');
%     drawnow
%     end
% end

%%

step = 50;
for a_start = 2000:7:5000
    a = a_start;
    b=a+step;
    time_start = onset_time(a);
    time_end = onset_time(b)+300;
    img_idx = trial_data(a:b);

    if(max(img_idx)<1000 && min(img_idx)>0 && max(diff(onset_time(a:b)))<350)
        img = [];onset_bar = [];
        for ii = 1:length(img_idx)
            img = [img, img_pool{img_idx(ii)}];
            onset_bar = [onset_bar, 0.8*ones([50,227,3])];
            onset_bar = [onset_bar, 0.5*ones([50,227,3])];
        end
    else
        continue
    end


    figure;
    set(gcf,'Position',[100 100 1700 900])
    subplot(5,1,[1,2])
    imshow(img)
    subplot(5,1,[3,4])
    hold on
    set(gca,'TickDir','none')
    LOC = 0;
    good_units = find(reliability_best>0.4);
    cm_here = colormap_matplotlib('plasma');
    cm_here = cm_here(1+floor(200*(1:length(good_units))/length(good_units)),:);
    data_all = zeros([1, time_end-time_start]);
    for units = 1:length(good_units)
        unit_here = good_units(units);
        raster_raw = GoodUnitStrc(unit_here).spiketime_ms;
        raster_time = raster_raw(raster_raw>time_start &raster_raw<time_end);
        scatter(raster_time, repmat(LOC, [1,length(raster_time)]),8,'filled',Marker='square',MarkerFaceColor=cm_here(units,:),MarkerEdgeAlpha=0);
        LOC = LOC+1;
        time_here = floor(raster_time-time_start)+1;
        d0 = zeros([1, time_end-time_start]);
        d0(time_here) = d0(time_here)+1;
        data_all = data_all+conv(d0,ones([1,30]),'same');
    end
    data_all = data_all./length(good_units);
    xline(onset_time(a+1:b),Color=[0.2,0.2,0.2],LineStyle="--",LineWidth=1)
    ylim([-2, LOC+2]);
    xticks([])
    ylabel('# Units')
    set_font

    subplot(5,1,5)
    plot(1:(time_end-time_start),data_all,'LineWidth',2,Color=[0.2,0.2,0.2])
    box off
    set(gca,'TickDir','none')
    xline(onset_time(a+1:b)-time_start,Color=[0.5,0.5,0.5],LineStyle="--",LineWidth=1)
    xticks(0:300:(time_end-time_start))
    xticklabels(num2cell(0:300:(time_end-time_start)))
    ylabel('% Unit recruited')
    sprintf('%d %d', a_start, max(data_all))
    xlabel('Time (ms)')
    ylim([0,max(data_all)])
    set_font
    break
end
%
drawnow
time_step_ms_here = 5;
time_one_window = 2400;
time_winsow_all = time_start:time_step_ms_here:time_end-time_one_window;
for tt = 1:length(time_winsow_all)
    set(gcf,'color','w')
    hText = findall(gcf, 'Type', 'text');
    for i = 1:length(hText)
        set(hText(i), 'FontName', 'Arial');
        set(hText(i), 'FontSize', 12);
    end
    ax = gca;
    ax.FontName = 'Arial';
    ax.FontSize = 12;
    tt_here = time_winsow_all(tt);
    tt_multi = 1+(tt-1)*time_step_ms_here;
    subplot(5,1,[1,2])
    img_t0 = floor(227*(tt_here-time_start)/300);
    img_t1 = img_t0+227*8;
    xlim([img_t0,img_t1])
    subplot(5,1,[3,4])
    xlim([tt_here,tt_here+time_one_window])
    subplot(5,1,5)
    xlim([tt_here,tt_here+time_one_window]-time_start)
    output = getframe(gcf);
    output   = frame2im(output);
    [imind, cm] = rgb2ind(output, 256);
    if(tt==1)
        imwrite(imind, cm, fullfile("Figs/F1/",'Illustration.gif'), 'gif', 'Loopcount', inf, 'DelayTime', 0.0167);
    else
        imwrite(imind, cm, fullfile("Figs/F1/",'Illustration.gif'), 'gif', 'WriteMode', 'append', 'DelayTime', 0.0167);
    end
end