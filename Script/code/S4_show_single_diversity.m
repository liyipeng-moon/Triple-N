clear;clc
root_dir = 'C:\Users\moonl\Desktop\NNN';
cd(root_dir)
addpath(genpath(pwd));
[proc_dir,raw_dir] = gen_dirs(root_dir);
reliability_thres = 0.4;
manual_data = readtable("exclude_area.xls");
stats_dir = fullfile(root_dir,"Figs/stats/");
mkdir Figs\
mkdir Figs\F4\
load img_pool.mat

order = randperm(12);
face_idx = 1000+[1:5,9:15,6,7,8,16:24];
face_idx = face_idx([order,order+12]);
body_idx = 1000+[26:31,43:48,50:61];
body_idx = body_idx([order,order+12]);

obj_idx = [1032:1042,1049,1025,1062:1072];
obj_idx = obj_idx([order,order+12]);
loc_idx = [face_idx, body_idx, obj_idx];

for interested_area =  {'MB3'}
    avg_rsp = [];
    B_si = [];
    ses_id = [];
    fob_rsp = [];
    reliab = [];
    wf = [];
    for aa = 1:length(interested_area)
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
                x1=manual_data.y1(search_area); x2=manual_data.y2(search_area);
                SI = pro1_data.B_SI;
                good_neuron_idx = find(pro1_data.pos>x1 & pro1_data.pos<x2 & pro1_data.reliability_best>0.4 & SI>0.5);
                avg_rsp = [avg_rsp; pro1_data.response_best(good_neuron_idx,1:1000)];
                
                B_si = [B_si, pro1_data.B_SI(good_neuron_idx)];
                ses_id = [ses_id, search_area*ones([1,length(good_neuron_idx)])];
                
                reliab = [reliab, pro1_data.reliability_best(good_neuron_idx)];
                wf = [wf, pro1_data.mean_psth(good_neuron_idx,:)];

                fob_rsp = [fob_rsp;pro1_data.response_best(good_neuron_idx,loc_idx)];
                return
            end
        end
    end
end
%%
fMRI_dataset_path = 'C:\Users\moonl\Desktop\NNN\NNN_Data\FMRI';
interested_subject=5;
hemi_here = 'lh';
fMRI_data = load(fullfile(fMRI_dataset_path,sprintf('S%d_%s_Rsp.mat',interested_subject,hemi_here)));
fMRI_data = fMRI_data.mean_brain_data;
fMRI_data = double(fMRI_data)./300;
lh_data = fMRI_data;
hemi_here = 'rh';
fMRI_data = load(fullfile(fMRI_dataset_path,sprintf('S%d_%s_Rsp.mat',interested_subject,hemi_here)));
fMRI_data = fMRI_data.mean_brain_data;
fMRI_data = double(fMRI_data)./300;
rh_data = fMRI_data;

rdm_lh = zeros([length(B_si), length(lh_data)]);
rdm_rh = zeros([length(B_si), length(rh_data)]);
pval_lh = zeros([length(B_si), length(lh_data)]);
pval_rh = zeros([length(B_si), length(rh_data)]);
start_parfor
tic
parfor unit_idx = 1:length(B_si)
    mean_rsp_unit = avg_rsp(unit_idx,:);
    [rdm_lh(unit_idx,:),pval_lh(unit_idx,:)] = corr(lh_data', mean_rsp_unit',type='Pearson');
    [rdm_rh(unit_idx,:),pval_rh(unit_idx,:)] = corr(rh_data', mean_rsp_unit',type='Pearson');
    rdm_lh(unit_idx,:) = rdm_lh(unit_idx,:)./reliab(unit_idx);
    rdm_rh(unit_idx,:) = rdm_rh(unit_idx,:)./reliab(unit_idx);
end
toc
clear rh_data lh_data fMRI_data
save(fullfile(proc_dir,sprintf('RDM_MSB_ZZ.mat')))
%% look at each unit  and select some unit with different corelation profile
close all
img_example = [];
for cc = 1:6
    for ins = 1
        idx = loc_idx(12*(cc-1)+[ins]);

            img_example = [img_example, img_pool{idx}];

    end
end
figure
set(gcf,'Position',[750 200 1150 410])
example_unit = [101,12,148,210];

for example_now = 1:length(example_unit)
    rsp_here = fob_rsp(example_unit(example_now),:);
    rsp_here = zscore(rsp_here);
    subplot(2, length(example_unit), example_now); hold on
    for cc = 1:6
        idx = 12*(cc-1)+[1:12];
        mm(cc) = mean(rsp_here(idx));
        ee(cc) = std(rsp_here(idx))./sqrt(12);
    end
    bar(1:6, mm,'FaceAlpha',0.1,'FaceColor',[1,0,0])
    errorbar(1:6, mm,ee,'LineStyle','none','Color','k','LineWidth',1,'CapSize',4)
    xticks([1:6])
    xticklabels({'Monkey Face','Human Face','Monkey Body','Animal','Natural Object','Manmade Object'})
    title(sprintf('Example Unit %d\n Body SI = %.02f',example_unit(example_now), B_si(example_unit(example_now))))
    ylabel('Firing rate (Hz)')
    set_font
end

subplot(2,4,[5:8])
imshow(img_example)
figsave(fullfile(root_dir,'Figs/F4/'),'exampleBody')