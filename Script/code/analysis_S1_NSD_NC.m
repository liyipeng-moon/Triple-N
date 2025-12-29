clear
root_dir = 'C:\Users\moonl\Desktop\NNN';
H5_dir = 'C:\Users\moonl\Desktop\NNN\NNN_Data\Raw\H5FILES';
prep_dir = 'C:\Users\moonl\Desktop\NNN\NNN_Data\processed';
cd(root_dir)
addpath(genpath(pwd));
%%
% GLOBAL PARA
face_idx = 1001:1024;
body_idx = 1000+[26:31,43:48,50:61];
obj_idx = setdiff(1025:1072, body_idx);
boot_times = 5;
basic_time_bin = 70:220;
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
bin_size = size(window_summary, 2);
random_order = randperm(1000);
find_sample = random_order(1:500);
test_sample = random_order(501:1000);
%%
cd(H5_dir)
parpool(20)
for ses_number_now = 1:90
    tic
    filename_here = dir(sprintf('ses%02d*h5',ses_number_now));
    filename_here = filename_here.name;
    metaname_here = dir(sprintf('ses%02d*mat',ses_number_now));
    metaname_here = metaname_here.name;
    meta_data = load(metaname_here);
    RasterData = h5read(filename_here, '/raster_matrix_img');
    PSTHData = h5read(filename_here, '/response_matrix_img');
    pre_onset = meta_data.global_params.pre_onset;
    psth_range = meta_data.global_params.PsthRange;
    % params
    UnitNum = size(PSTHData,1);
    row_example = zeros([1, UnitNum]);
    UnitType = row_example;
    mean_psth = zeros([UnitNum, length(psth_range)]);
    response_basic = zeros([UnitNum, 1072]);
    response_best = zeros([UnitNum, 1072]);
    reliability_basic = row_example;
    pos = row_example;
    F_SI = row_example;
    B_SI = row_example;
    O_SI = row_example;
    snr = row_example;
    snrmax = row_example;
    r_search_pool = zeros([UnitNum, bin_size]);
    reliability_best = row_example;
    reliability_find_testset = row_example;
    best_r_time1 = row_example;
    best_r_time2 = row_example;
    
    figure;title(ses_number_now)
    imagesc(squeeze(mean(RasterData,1)))
    title(filename_here,Interpreter="none")
    drawnow

    parfor unit = 1:UnitNum
        single_r_pool = zeros([1, bin_size]);
        raster_this_neuron = squeeze(RasterData(unit,:,:));
        pos(unit)=meta_data.GoodUnitStrc(unit).spikepos(2);
        UnitType(unit)=meta_data.GoodUnitStrc(unit).unittype;
        basic_raster = mean(raster_this_neuron(:, pre_onset+basic_time_bin),2);
        reliability_basic(unit) = give_me_nc(basic_raster, meta_data.img_idx, test_sample,boot_times);
        
        response_basic(unit, :) = mean(PSTHData(unit,1:1072, pre_onset+basic_time_bin),3);
        for bin_idx = 1:bin_size
            t1 = window_summary(1,bin_idx);
            t2 = window_summary(2,bin_idx);
            selected_raster = mean(raster_this_neuron(:, pre_onset+(t1:t2)),2);
            single_r_pool(bin_idx)=give_me_nc(selected_raster, meta_data.img_idx, find_sample,boot_times);
        end
        
        r_search_pool(unit,:)=single_r_pool;
        [~,b]=max(single_r_pool);

        best_r_time1(unit)=[window_summary(1, b)];
        best_r_time2(unit)=[window_summary(2, b)];

        selected_raster = mean(raster_this_neuron(:, pre_onset+(window_summary(1, b):window_summary(2, b))),2);
        reliability_find_testset(unit) = give_me_nc(selected_raster, meta_data.img_idx, test_sample,boot_times);
        reliability_best(unit) = give_me_nc(selected_raster, meta_data.img_idx, 1:1000,boot_times);

        best_time_now = (best_r_time1(unit):best_r_time2(unit))+pre_onset;
        response_best(unit,:)=mean(PSTHData(unit, :, best_time_now),3);

        data_face = mean(PSTHData(unit, face_idx, best_time_now),3);
        data_body = mean(PSTHData(unit, body_idx, best_time_now),3);
        data_obj = mean(PSTHData(unit, obj_idx, best_time_now),3);

        F_SI(unit) = CalcSI(data_face, [data_body,data_obj] );
        B_SI(unit) = CalcSI(data_body, [data_face,data_obj]);
        O_SI(unit) = CalcSI(data_obj, [data_body,data_face]);

        psth_this_neuron = squeeze(PSTHData(unit, 1:1000,:));
        mean_psth(unit, :)= squeeze(mean(psth_this_neuron));
        tmp_mean = squeeze(mean(psth_this_neuron));

        [~,time_here] = sort(tmp_mean(pre_onset+(1:300)),'descend');
        peak_rsp = mean(tmp_mean(pre_onset+time_here(1)));
        bsl_rsp = psth_this_neuron(:,pre_onset+[-20:20]);
        bsl_rsp = mean(bsl_rsp,2);
        std_baseline = std(bsl_rsp);
        bsl_rsp = mean(bsl_rsp);
        max_peak = max(psth_this_neuron(:, pre_onset+time_here(1)));
        snr(unit) = (peak_rsp-bsl_rsp)./std_baseline;
        snrmax(unit) = (max_peak-bsl_rsp)./std_baseline;
    end
    clear RasterData PSTHData ppm data_face data_body data_obj best_time_now selected_raster
    clear b t1 t2 bin_idx basic_raster single_r_pool raster_this_neuron
    
    file_name_save = fullfile(prep_dir,sprintf('Processed_%s.mat',filename_here(1:end-3)));

    save(file_name_save, 'UnitNum', ...
        "F_SI", 'B_SI','O_SI', ...
        'response_best', 'response_basic',...
        "UnitType",'mean_psth','pos', ...
        'r_search_pool','best_r_time1','best_r_time2', ...
        'reliability_best','reliability_find_testset',"reliability_basic","snr","snrmax" ...
        );
    fprintf('ses %d used %.02f mins \n',ses_number_now, toc/60)

    close all
end
