clear
clc
load DIRS.mat
GoodUnit_Dir=fullfile(root_dir,'Data','Raw','GoodUnit');
GoodLFP_Dir=fullfile(root_dir,'Data','Raw','GoodLFP');

all_good_units_data = dir(fullfile(GoodUnit_Dir,'GoodUnit_*mat'));
addpath(genpath(root_dir))
%%
for ses_idx = 1:90
    filename_goodunit = all_good_units_data(ses_idx).name;
    [day,subject,gnumber] = parse_name(filename_goodunit);
    filename_with_dir = fullfile(GoodUnit_Dir,filename_goodunit);
    fprintf('Converting %s \n', filename_with_dir)
    data_here = load(filename_with_dir);
    GoodUnitStrc = load(filename_with_dir).GoodUnitStrc;
    size_matrix = size(GoodUnitStrc(1).response_matrix_img);
    response_matrix_img = single(zeros([length(GoodUnitStrc), size_matrix(1), size_matrix(2)]));
    for uu = 1:length(GoodUnitStrc)
        response_matrix_img(uu,:,:) = GoodUnitStrc(uu).response_matrix_img;
    end
    size_matrix = size(GoodUnitStrc(1).Raster);
    raster_matrix_img = uint8(zeros([length(GoodUnitStrc), size_matrix(1), size_matrix(2)]));
    for uu = 1:length(GoodUnitStrc)
        raster_matrix_img(uu,:,:) = GoodUnitStrc(uu).Raster;
    end
    save_name = fullfile(H5_dir, sprintf('ses%02d_%s_%s_%s',ses_idx,day,subject,gnumber));
    filename_h5 = sprintf('%s.h5',save_name);
    h5create(filename_h5, '/response_matrix_img', size(response_matrix_img),'Datatype','single');
    h5write(filename_h5, '/response_matrix_img', response_matrix_img);
    h5create(filename_h5, '/raster_matrix_img', size(raster_matrix_img),'Datatype','uint8');
    h5write(filename_h5, '/raster_matrix_img', raster_matrix_img);
    global_params = data_here.global_params;
    meta_data = data_here.meta_data;
    trial_ML = data_here.trial_ML;
    img_idx = data_here.meta_data.trial_valid_idx(data_here.meta_data.trial_valid_idx~=0);
    GoodUnitStrc = rmfield(GoodUnitStrc,'response_matrix_img');
    GoodUnitStrc = rmfield(GoodUnitStrc,'Raster');

    
    filename_LFP = fullfile(GoodLFP_Dir, sprintf('GoodLFP_%s',filename_goodunit(10:end)));
    LFP_file_data = load(filename_LFP);
    LFP_META = LFP_file_data.LFP_META;
    LFP_Data = single(LFP_file_data.LFP_data_trial_wise(find(data_here.meta_data.trial_valid_idx~=0),:,:));
    h5create(filename_h5, '/LFP_Data', size(LFP_Data));
    h5write(filename_h5, '/LFP_Data', LFP_Data);
    filename_h5 = sprintf('%s_info.mat',save_name);
    save(filename_h5,"img_idx","global_params","trial_ML","GoodUnitStrc","meta_data","LFP_META")
end