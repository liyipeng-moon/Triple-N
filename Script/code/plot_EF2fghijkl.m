clear
clc
close all
root_dir = 'C:\Users\moonl\Desktop\NNN';
[proc_dir,raw_dir] = gen_dirs(root_dir);
addpath(genpath(root_dir));
mkdir Figs\
mkdir Figs\F2S\
H5Folder = 'C:\Users\moonl\Desktop\NNN\NNN_Data\Raw\H5FILES';
ProcFolder = 'C:\Users\moonl\Desktop\NNN\NNN_Data\Processed';
IT_session = [1:70, 88];
% Spike Quality, for IT data.
reliability_thres = 0.4;
UnitType = [];
fr = [];
isi_vio = [];
amp = [];
snr = [];
rpv = [];
wf = {};
rr_array = [];
wf_max = [];
wf_min = [];
Area_Array = [];
Unthreshold_Area = [];
for ses_idx = 1:90
    if(sum(ses_idx==IT_session))
        A=1;
    else
        A=2;
    end
    proc1_file_name = dir(fullfile(ProcFolder,sprintf('Processed_ses%02d*', ses_idx)));
    proc1_file_name = proc1_file_name.name;
    pro1_data = load(fullfile(ProcFolder,proc1_file_name));
    metaname_here = dir(fullfile(H5Folder,sprintf('ses%02d*mat',ses_idx)));
    metaname_here = metaname_here.name;
    meta_data = load(fullfile(H5Folder,metaname_here));
    rr = pro1_data.reliability_best;
    unit_here = find(rr>reliability_thres);
    UnitType = [UnitType, pro1_data.UnitType(unit_here)];
    for uu = 1:length(unit_here)
        wf{end+1} = meta_data.GoodUnitStrc(unit_here(uu)).waveformBC;
        fr(end+1) = length(meta_data.GoodUnitStrc(unit_here(uu)).spiketime_ms)./(meta_data.meta_data.onset_time_ms(end)/1000);
        isi_vio(end+1) = meta_data.GoodUnitStrc(unit_here(uu)).qm.fractionRPVs_estimatedTauR;
        wf_here = wf{end};
        [a,b] = find(abs(wf_here)==max(abs(wf_here(:))));
        wf_now = wf_here(a,:);
        amp(end+1) = max(wf_now(:))-min(wf_now(:));
        snr(end+1) = meta_data.GoodUnitStrc(unit_here(uu)).qm.signalToNoiseRatio;
        wf_min(end+1) = min(wf_now(:));
        wf_max(end+1) = max(wf_now(:));
        Area_Array(end+1)=A;
    end
    rr_array = [rr_array, rr];
    Unthreshold_Area = [Unthreshold_Area, A*ones(size(rr))];
    ses_idx
end
%% Further curate
UTC = UnitType;
wf_pool = zeros([61, length(wf)])';
for uu = 1:length(wf)
    wf_here = wf{uu};
    [a,b] = find(abs(wf_here)==max(abs(wf_here(:))));
    wf_single = wf_here(a,:);
    if(max(wf_single)>-min(wf_single))
        UTC(uu)=3;
    end
    wf_pool(uu,:) = wf_single./max(abs(wf_single));
end
%%
IT_r = rr_array(Unthreshold_Area==1);
fprintf('For IT: \n%d out of %d visual responsive units with reliability over 0.4\n',sum(IT_r>reliability_thres),length(IT_r))
IT_UTC = UTC(Area_Array==1);
fprintf('Among them %d is SU, %d MU, %d is NonSomSU\n', sum(IT_UTC==1),sum(IT_UTC==2),sum(IT_UTC==3))

EVC_r = rr_array(Unthreshold_Area==2);
fprintf('For EVC: \n%d out of %d visual responsive units with reliability over 0.4\n',sum(EVC_r>reliability_thres),length(EVC_r))
EVC_UTC = UTC(Area_Array==2);
fprintf('Among them %d is SU, %d MU, %d is NonSomSU\n', sum(EVC_UTC==1),sum(EVC_UTC==2),sum(EVC_UTC==3))

UnitType = UTC;
%%
rng(1009)
close all
figure(10)
set(gcf,'Position',[400 200 1200 700])
gain_to_uv = 1000*1000*meta_data.meta_data.IMEC_AP_META.imAiRangeMax/meta_data.meta_data.IMEC_AP_META.imMaxInt/500;
vals = {fr,isi_vio,amp,snr};
val_name = {'Firing rate (Hz)','Refractory period violations (%)','Amplitude (uv)','Signal to noise ratio'};
edge_val = {-2:1:25, -0.1:0.05:2,-2:8:250,-5:8:300};
for uu = 1:3
    for pp = 1:length(vals)
        data = vals{pp}(UnitType==uu);
        subplot(5,2+length(val_name),(uu-1)*(2+length(val_name))+2+pp); hold on
        histogram(data,edge_val{pp},Normalization="count",FaceColor=[0.1,0.1,0.1],EdgeAlpha=0)
        xlim(edge_val{pp}([1,end]))
        if(uu==3)
            xlabel(val_name{pp})
        end
        if(pp==1)
            ylabel('# Units')
        end
        set(gca,'TickDir','out')
        box on
        set_font
    end

    order = [];
    thisUnit = find(UnitType==uu);
    if(uu<3)
        [a,b]=sort(wf_min(thisUnit),'ascend');
        order = randsample(find(abs(a)<350,1)+[1:10],3);
    else
        [a,b]=sort(wf_max(thisUnit),'descend');
        order = randsample(find(abs(a)<400,1)+[1:10],3);
    end

    sub_val = (uu-1)*(2+length(val_name))+1;
    subplot(5,2+length(val_name),sub_val:sub_val+1);

    for x = 1:length(order)
        hold on
        wf_here = wf{thisUnit(b(order(x)))};
        plot_waveform(wf_here, x)
    end

    switch uu
        case 1
            ylabel(sprintf('SU\n n = %d\n Voltage(uv)',length(thisUnit)))
            set(gca,'XColor',[1,1,1])
            title(' Example Waveform')
        case 2
            ylabel(sprintf('MUA\n n = %d\n Voltage(uv)',length(thisUnit)))
            set(gca,'XColor',[1,1,1])
        case 3
            ylabel(sprintf('NonSomatic\n n = %d\n Voltage(uv)',length(thisUnit)))
    end
    xlim([2,2.1*(1+length(order))])
    set_font
end

% If you are interested, there are also some spike-cluster, such as Fast spiking and regular spiing neuron
% This is just a simple illustration
% please refer to Wavemap paper (Kenji Lee et al, 2021, elife) for better method
% wf_all = [];
% for wf_array = 1:length(UTC)
%     if(UTC(wf_array)<2)
%         wf_now = wf{wf_array};
%         [a,b] = find(wf_now==max(abs(wf_now(:))) | wf_now==-max(abs(wf_now(:))));
%         wf_all = [wf_all; wf_now(a,:)];
%     end
% end
% wf_all = wf_all';
% wf_all = - wf_all ./ min(wf_all);
% cc = kmeans(wf_all',4); 
% xt = (1:61)/30000;
% figure; 
% for c = 1:4
%     nexttile; hold on
%     data_now = wf_all(:,cc==c);
%     plot(xt,data_now,Color=[0.5,0.5,0.5,0.03])
%     plot(xt,median(data_now'),Color='k',LineWidth=1)
%     ylim([-1,1]);
% end
%
root_dir = 'C:\Users\moonl\Desktop\NNN';
cd(root_dir)
addpath(genpath(pwd));
[proc_dir,raw_dir] = gen_dirs(root_dir);
stats_dir = fullfile(root_dir,"Figs/stats/");
manual_data = readtable("exclude_area.xls");
coor_data = readtable('AreaXYZ.xlsx');
all_proc_data = dir(fullfile(proc_dir,'Proces*'));
AP_range = max(coor_data.A)-min(coor_data.A);
all_cm = colormap_matplotlib('plasma');
all_cm = all_cm(1:200,:);
subplot(5,2+length(val_name),[3*(2+length(val_name))+[4,5],4*(2+length(val_name))+[4,5]]);
hold on
xlim([35,150]);
ylim([150,360]);
rr_pool = [];
r_basic_pool = [];
r_find_pool = [];
for ses_idx = 1:90
    proc1_file_name = dir(fullfile(proc_dir,sprintf('Processed_ses%02d*', ses_idx)));
    proc1_file_name = proc1_file_name.name;
    pro1_data = load(fullfile(proc_dir,proc1_file_name));
    unit_here = find(pro1_data.reliability_best>0.4);
    d1 = pro1_data.best_r_time1(unit_here);
    d2 = pro1_data.best_r_time2(unit_here);
    rr_pool = [rr_pool;pro1_data.r_search_pool(unit_here,:)];
    r_basic_pool = [r_basic_pool, pro1_data.reliability_basic(unit_here)];
    r_find_pool = [r_find_pool, pro1_data.reliability_best(unit_here)];
    ses_row = find(manual_data.SesIdx==ses_idx);
    ap_here = coor_data.A(manual_data.RoiIndex(ses_row(1)));
    ap_percent_here = (ap_here-min(coor_data.A))/AP_range;
    CC = floor(ap_percent_here*199)+1;
    cm_here = all_cm(CC,:);
    errorbar(mean(d1),mean(d2), std(d2)./sqrt(length(d1)),Color=cm_here,LineWidth=2, CapSize=0)
    errorbar(mean(d1),mean(d2), std(d1)./sqrt(length(d1)),'horizontal',Color=cm_here,LineWidth=2, CapSize=0)
    scatter(mean(d1),mean(d2),12,'filled',Marker='square',MarkerFaceColor=cm_here,MarkerEdgeAlpha=1)
end

xlabel('Time Start (ms)')
ylabel('Time End (ms)')
cc=colorbar;
cc.Label.String='Posterior                  Anterior';
cc.TickLabels=[];
colormap(gca,all_cm)
set_font

all_cm = colormap_matplotlib('plasma');
all_cm = all_cm(1:200,:);
rr_pool = [];
r_basic_pool = [];
r_find_pool = [];
for ses_idx = 1:length(all_proc_data)
    proc1_file_name = dir(fullfile(proc_dir,sprintf('Processed_ses%02d*', ses_idx)));
    proc1_file_name = proc1_file_name.name;
    pro1_data = load(fullfile(proc_dir,proc1_file_name));
    d1 = pro1_data.best_r_time1;
    d2 = pro1_data.best_r_time2;
    rr_pool = [rr_pool;pro1_data.r_search_pool];
    r_basic_pool = [r_basic_pool, pro1_data.reliability_basic];
    r_find_pool = [r_find_pool, pro1_data.reliability_best];
end
bin_1 = 20:10:200;
bin_2 = 90:10:390;
window_summary = [];
for b1 = 1:length(bin_1)
    for b2 = 1:length(bin_2)
        if(bin_1(b1)<bin_2(b2))
            window_summary(1,end+1)=bin_1(b1);
            window_summary(2,end)=bin_2(b2);
        else
            window_summary(1,end+1)=bin_2(b2);
            window_summary(2,end)=bin_1(b1);
        end
    end
end

plot_unit = find(~isnan(r_basic_pool)&~isnan(r_find_pool)&r_basic_pool>0);
subplot(5,2+length(val_name),[3*(2+length(val_name))+[1:2],4*(2+length(val_name))+[1:2]]);
hold on
bin_scatter(r_basic_pool(plot_unit),r_find_pool(plot_unit),0:0.01:1,0:0.01:1)
[h,p,ci,stats] = ttest(r_basic_pool,r_find_pool);
fid=fopen(fullfile(root_dir,"Figs/stats",'F2S.txt'),'w');
fprintf(fid,'reliability difference, t(%d)=%.02f,log(p)=%d',stats.df,stats.tstat,log10(p));
xlabel('Reliability (70-220ms)')
ylabel('Reliability (best)')
set_font

figsave(fullfile(root_dir,'Figs/FS1/'),sprintf('FS1_2'))