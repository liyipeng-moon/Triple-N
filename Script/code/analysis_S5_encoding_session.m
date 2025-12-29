function S5_Encoding_Session(area_now)
root_dir = 'C:\Users\moonl\Desktop\NNN';
cd(root_dir)
addpath(genpath(pwd));
[proc_dir,raw_dir] = gen_dirs(root_dir);
reliability_thres = 0.4;
manual_data = readtable("exclude_area.xls");
interested_time_bin = 1:350;
ses_idx = manual_data.SesIdx(area_now);
proc1_file_name = dir(fullfile(proc_dir,sprintf('Processed_ses%02d*', ses_idx)));
proc1_file_name = proc1_file_name.name;
pro1_data = load(fullfile(proc_dir,proc1_file_name));
x1 = manual_data.y1(area_now); x2 = manual_data.y2(area_now);
unit_here = find(pro1_data.pos>x1 & pro1_data.pos<x2 & pro1_data.reliability_best > reliability_thres);
r_here = pro1_data.reliability_best(unit_here);
all_neuron_rsp = pro1_data.response_best(unit_here,1:1000);
metaname_here = dir(fullfile(raw_dir,sprintf('ses%02d*mat',ses_idx)));
metaname_here = metaname_here.name;
meta_data = load(fullfile(raw_dir,metaname_here));
filename_here = dir(fullfile(raw_dir,sprintf('ses%02d*h5',ses_idx)));
filename_here = filename_here.name;
PSTHData = h5read(fullfile(raw_dir,filename_here), '/response_matrix_img');
PSTHData = PSTHData(unit_here,1:1000,meta_data.global_params.pre_onset+interested_time_bin);
[s,model_name] = load_embedding;
%%
n_space = length(s);
pred_r2_array_t = zeros([length(r_here),n_space, 350]);
pred_r_array = zeros([length(r_here),n_space]);
pred_r2_array = zeros([length(r_here),n_space]);
start_parfor
tic
parfor neuron_idx = 1:length(r_here)
    warning('off', 'all');
    single_neuron_data = squeeze(PSTHData(neuron_idx,:,:));
    rsp_now = zscore(all_neuron_rsp(neuron_idx,:));
    [pred_r_array(neuron_idx,:),pred_r2_array(neuron_idx,:),best_numbers] = compare_encoders(s, rsp_now, 10);
    for l = 1:n_space
        pred_r2_array_t(neuron_idx,l,:) = lscov_r_multi_y(s{l}, single_neuron_data, best_numbers(l));
    end
end
toc
save(fullfile(proc_dir,sprintf('s5_encode_ses_%03d.mat', area_now)),"pred_r2_array_t","pred_r_array","pred_r2_array","r_here","unit_here")
end