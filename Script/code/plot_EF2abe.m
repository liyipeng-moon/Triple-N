close all
clear
clear;clc
root_dir = 'C:\Users\moonl\Desktop\NNN';
cd(root_dir)
addpath(genpath(pwd));
[proc_dir,raw_dir] = gen_dirs(root_dir);
all_info = dir(fullfile(raw_dir,'\*info*'));
showrange = -5:0.2:5;
cm_here = colormap_matplotlib('set1',9);
MK_array = 1:5;
%%
type_pool ={};
for ii = 1:length(all_info)
    info_name = all_info(ii).name;
    MK_IDX(ii) = str2num(info_name(15));
    % Load about unit type
    GoodUnit = load(fullfile(raw_dir,info_name)).GoodUnitStrc;
    type_array = [];
    for unit = 1:length(GoodUnit)
        wf_here = GoodUnit(unit).waveformBC;
        abs_val = abs(wf_here);
        [a,b] = find(abs_val==max(abs_val(:)));
        waveform = wf_here(a,:);
        if(max(waveform)>-min(waveform))
            type_array(unit)=3;
        else
            type_array(unit)=GoodUnit(unit).unittype;
        end
    end
    type_array(type_array==4)=3;
    type_pool{ii} = type_array;
    % Load about eye
    meta_data = load(fullfile(raw_dir,info_name)).meta_data;
    trial_success = find(meta_data.trial_valid_idx);
    trial_fail = find(~meta_data.trial_valid_idx);
    eye_data = meta_data.eye_matrix;
    eye_data = squeeze(mean(eye_data,3));
    [CIx(ii,:)] = std(eye_data(1,trial_success));
    [CIy(ii,:)] = std(eye_data(2,trial_success));
    eye_data_pool{ii} = eye_data(:,trial_success);
    fprintf('Loading quality for ses %02d out of %d \n',ii,length(all_info))
end
%% Plot result
close all
figure
set(gcf,'Position',[1 30 1800 405])
subplot(2,6,7);hold on;
bin_scatter_eye(eye_data_pool{5}(1,:),eye_data_pool{5}(2,:),showrange,showrange);
xlabel('Horizontal position (deg)')
ylabel('Vertical position (deg)')
title('Example session 5')

subplot(2,6,8);hold on;
bin_scatter_eye(eye_data_pool{40}(1,:),eye_data_pool{40}(2,:),showrange,showrange);
xlabel('Horizontal position (deg)')
ylabel('Vertical position (deg)')
title('Example session 40')

subplot(2,6,9);hold on;
for mk = 1:5
    ss = find(MK_IDX==mk);
    DataX = CIx(ss);
    DataY = CIy(ss);
    scatter(DataX,DataY,12,'filled',MarkerFaceColor=cm_here(mk,:),MarkerFaceAlpha=0.4,DisplayName=['M' num2str(mk)])
end

s = 40;
scatter(CIx(s),CIy(s),24,Marker='*',MarkerEdgeColor=[0,0,0],MarkerFaceColor=[0,0,0],HandleVisibility='off')
s = 5;
scatter(CIx(s),CIy(s),24,Marker='*',MarkerEdgeColor=[0,0,0],MarkerFaceColor=[0,0,0],HandleVisibility='off')
legend(Box="off",Location="best")
xlabel('Horizontal SD (deg)')
ylabel('Vertical SD (deg)')
xlim([0,2])
ylim([0,2])
colorbar

subplot(2,6,3)
histogram(CIx,0:0.1:2,EdgeAlpha=1,FaceColor=[0.5,0.5,0.5])
l = ylim;ylim([l(1),l(2)*4]);xlim([0,2])
colorbar; axis off

subplot(2,6,10)
histogram(CIx,0:0.1:2,EdgeAlpha=1,FaceColor=[0.5,0.5,0.5],Orientation="horizontal")
ylim([0,2]);l = xlim;xlim([l(1),l(2)*4])
colorbar; axis off


subplot(2,6,11);hold on
for ii = 1:length(all_info)
    array = type_pool{ii};
    SO_ratio(ii) = sum(array<3)./length(array);
    SU_ratio(ii) = sum(array==1)./length(array);
end

for mk = 1:5
    ss = find(MK_IDX==mk);
    scatter(SO_ratio(ss),SU_ratio(ss),12,'filled',MarkerFaceColor=cm_here(mk,:),MarkerFaceAlpha=0.4,DisplayName=['M' num2str(mk)])
end
ldg = legend(Box="off",Location="best",FontSize=8);
xlabel('Percentage of somatic unit')
ylabel('Percentage of single unit')
xlim([0.4,0.8])
ylim([0,0.4])
colorbar

subplot(2,6,5)
histogram(SO_ratio,0.4:0.02:0.8,EdgeAlpha=1,FaceColor=[0.5,0.5,0.5])
l = ylim;ylim([l(1),l(2)*4]);xlim([0.4,0.8])
colorbar; axis off

subplot(2,6,12)
histogram(SU_ratio,0:0.02:0.4,EdgeAlpha=1,FaceColor=[0.5,0.5,0.5],Orientation="horizontal")
l = xlim;xlim([l(1),l(2)*4]);ylim([0,0.4])
colorbar; axis off


% Monkey Wise
bar_val = [];
for mk = 1:5
    ss_this_monkey = find(MK_IDX==mk);
    array_this_monkey = [];
    for s = ss_this_monkey
        array_this_monkey = [array_this_monkey, type_pool{s}];
    end
    for t = 1:3
        bar_val(mk,t) = mean(array_this_monkey==t);
    end
end
bar_val = bar_val';
subplot(2,6,6); hold on
h = bar(bar_val,'FaceAlpha',0.4,EdgeAlpha=0);
for mk = 1:5
    h(mk).FaceColor = cm_here(mk,:);
end
xticks([1,2,3]);
xticklabels({'SU','MUA','NonSomatic'});
ylabel('Percentage')
colorbar

subplot(2,6,4); hold on
h = bar(bar_val','stacked','FaceAlpha',0.4,EdgeAlpha=0);
xticks([1:5]);
xlim([0,8]);
yticks([0:0.5:1])
legend({'SU','MUA','NonSomatic'},'box','off');
ylabel('Proportion')
xlabel('Subject')
    
MTX = [SU_ratio; SO_ratio-SU_ratio; 1-(SO_ratio)];
[a,b] = sort(SU_ratio);
MTX = MTX';
subplot(2,6,[1,2])
bar(MTX,'stacked','FaceAlpha',0.4,EdgeAlpha=0);
xticks([]);
yticks([0:0.5:1])
ylabel('Proportion')
xlabel('Session')
figsave(fullfile(root_dir,"Figs/","FS1"),'FS_Fixation_UnitType')