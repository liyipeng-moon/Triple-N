close all
clear;clc
root_dir = 'C:\Users\moonl\Desktop\NNN';
cd(root_dir)
addpath(genpath(pwd));
[proc_dir,raw_dir] = gen_dirs(root_dir);
reliability_thres = 0.4;


addpath(genpath(pwd))
manual_data = readtable("exclude_area.xls");
interested_time_point = 1:400;
meta_example = load("ses10_240726_M5_3_info.mat"');
load img_pool.mat
load Clus_savee.mat

interested_area={'MO','AO','MF','AF','MB','AB','LP','PI','CL','AM'};
areaname={'M-Object','A-Object','M-Face','A-Face','M-Body','A-Body','Scene1','Scene2','M-Color','A-Color'};
cm_here = colormap_matplotlib('Set1',9);
figure
set(gcf,'Position',[500 200 1200 650])

for aa = 1:length(interested_area)
    for interested_unit = 1:3
        PSTH_CombinedData = [];
        ses_number = [];
        avg_rsp = [];
        for search_area = 1:height(manual_data)
            if(strcmp(manual_data.AREALABEL{search_area}(1:2),interested_area{aa}))
                ThisSes_idx = manual_data.SesIdx(search_area);

                proc1_file_name = dir(fullfile(proc_dir,sprintf('Processed_ses%02d*', ThisSes_idx)));
                proc1_file_name = proc1_file_name.name;
                pro1_data = load(fullfile(proc_dir,proc1_file_name));

                filename_here = dir(fullfile(raw_dir,sprintf('ses%02d*h5',ThisSes_idx)));
                filename_here = filename_here.name;
                metaname_here = dir(fullfile(raw_dir, sprintf('ses%02d*mat',ThisSes_idx)));
                metaname_here = metaname_here.name;
                meta_data = load(metaname_here);
                PSTHData = h5read(filename_here, '/response_matrix_img');


                x1=manual_data.y1(search_area);x2=manual_data.y2(search_area);
                good_neuron_idx = find(pro1_data.pos>x1 & pro1_data.pos<x2 & pro1_data.reliability_best>0.4);
                idx = clus_save{search_area}(find(clus_save{search_area}));

                good_neuron_idx = good_neuron_idx(idx==interested_unit);
                data_here = PSTHData(good_neuron_idx,:, meta_example.global_params.pre_onset+interested_time_point);
                PSTH_CombinedData = [PSTH_CombinedData; data_here];
                ses_number = [ses_number, ones([1, length(good_neuron_idx)])*search_area];
                switch interested_area{aa}(2)
                    case 'F'
                        SI = pro1_data.F_SI;
                    case 'O'
                        SI = pro1_data.O_SI;
                    case 'B'
                        SI = pro1_data.B_SI;
                    otherwise
                        SI = pro1_data.reliability_best;
                end
                avg_rsp = [avg_rsp; pro1_data.response_best(pro1_data.reliability_best>0.4 & SI>0.2,1:1000)];
            end
        end

        new_psth_combined = zeros(size(PSTH_CombinedData));



        for uu = 1:size(PSTH_CombinedData)
            data = squeeze(mean(PSTH_CombinedData(uu,:,1:300),2));
            new_psth_combined(uu,:,1:300) = (PSTH_CombinedData(uu,:,1:300)-mean(data))./std(data);
        end
        avg_rsp = zscore(avg_rsp,0,2);
        mean_avg = mean(avg_rsp);
        [a,b] = sort(mean_avg);

        nexttile
        hold on
        for i = [1,20]
            data_here = mean(new_psth_combined(:,b(50*(i-1)+(1:50)),1:300),1);
            m = mean(data_here);
            e = std(data_here)/sqrt(100);
            shadedErrorBar(1:300,m,e,'LineProps',{'Color',cm_here(interested_unit,:)},'patchSaturation',0.5);
        end
        data_here = mean(new_psth_combined(:,:,1:300),1);
        m = mean(data_here);
        e = std(data_here)/sqrt(100);
        shadedErrorBar(1:300,m,e,'LineProps',{'Color',[0.5,0.5,0.5]},'patchSaturation',0.5);
        xlim([0,300])
        if(interested_unit==1)
            ylabel(sprintf('%s \n Norm. firing rate', areaname{aa}))
        end
        if(aa<3)
            title(sprintf('Cluster %d',interested_unit))
        end
        if(aa>8)
            xlabel('Time (ms)')
        end
        set_font
        drawnow

        title(sprintf('n=%d',size(new_psth_combined,1)))
    end
    if(mod(aa,2)==1)
        nexttile
        axis off
    end
end

figsave(fullfile(root_dir,'Figs/F3'),sprintf('F3SS'))