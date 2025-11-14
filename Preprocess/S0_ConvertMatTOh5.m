clear
clc
root_dir = 'C:\Users\moonl\Desktop\NNN';
GU_folder = 'C:\Users\moonl\Desktop\NNN\NNN_Data\Raw\GoodStruct'; 
all_good_units_data = dir(fullfile(GU_folder,'GoodUnit*mat'));
addpath(genpath(root_dir))
%%
for ses_idx = 1:66
    file_name_here = all_good_units_data(ses_idx).name;
    [day,subject,gnumber] = parse_name(file_name_here);

    save_name = fullfile(root_dir,"Data","Raw","H5FILES", sprintf('ses%02d_%s_%s_%s',ses_idx,day,subject,gnumber));
    save_name = fullfile('C:\Users\moonl\Desktop\NNN\NNN_Data\H5', sprintf('ses%02d_%s_%s_%s',ses_idx,day,subject,gnumber));

    file_name =fullfile(GU_folder,file_name_here);
    
    fprintf('Converting %s \n', file_name)
    data_here = load(file_name);
    
    GoodUnitStrc = load(file_name).GoodUnitStrc;
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
    filename = sprintf('%s.h5',save_name);
    h5create(filename, '/response_matrix_img', size(response_matrix_img),'Datatype','single');
    h5write(filename, '/response_matrix_img', response_matrix_img);
    h5create(filename, '/raster_matrix_img', size(raster_matrix_img),'Datatype','uint8');
    h5write(filename, '/raster_matrix_img', raster_matrix_img);
    global_params = data_here.global_params;
    meta_data = data_here.meta_data;
    trial_ML = data_here.trial_ML;
    img_idx = data_here.meta_data.trial_valid_idx(data_here.meta_data.trial_valid_idx~=0);
    GoodUnitStrc = rmfield(GoodUnitStrc,'response_matrix_img');
    GoodUnitStrc = rmfield(GoodUnitStrc,'Raster');
    LFP_file_name = fullfile(GU_folder, sprintf('GoodLFP_%s',file_name_here(10:end)));
    LFP_file_data = load(LFP_file_name);
    LFP_META = LFP_file_data.LFP_META;
    LFP_Data = single(LFP_file_data.LFP_data_trial_wise(find(data_here.meta_data.trial_valid_idx~=0),:,:));
    h5create(filename, '/LFP_Data', size(LFP_Data));
    h5write(filename, '/LFP_Data', LFP_Data);
    filename = sprintf('%s_info.mat',save_name);
    save(filename,"img_idx","global_params","trial_ML","GoodUnitStrc","meta_data","LFP_META")
end