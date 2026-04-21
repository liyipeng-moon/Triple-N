clear
load('DIRS.mat')

% Load Brain Data for VTC ROI
load(fullfile(FMRI_DIR,'ROI_data.mat'))
ROI_Here = [1:6];
compare_vec_brain = struct;
roi_data = [];
voxel_info_1 = [];voxel_info_2 = [];voxel_info_3 = [];
vv = [];
for ROI_idx = ROI_Here
    for hh = [1,2]
        for ss = [1,2,5,7]
            brain_data = ROI_data{ss,ROI_idx,hh};
            vertex_num = size(brain_data,1);
            roi_data = [roi_data; brain_data];
            voxel_info_1 = [voxel_info_1, ROI_idx*ones([1, vertex_num])]; % which ROI
            voxel_info_2 = [voxel_info_2, hh*ones([1, vertex_num])]; % which hemisphere
            voxel_info_3 = [voxel_info_3, ss*ones([1, vertex_num])]; % which subject
            vv(end+1) = vertex_num;
        end
    end
end

% Load Brain Data for EVC
load(fullfile(FMRI_DIR,'ROI_info.mat'))
HH = ['l','r'];
for h = [1, 2]
    for s= [1,2,5,7]
        Vertex_now = find(eval(sprintf('S%d_%sh_EVC',s,HH(h))));
        vertex_num = length(Vertex_now);
        vertex_num
        brain_data = load(fullfile(FMRI_DIR,sprintf('S%d_%sh_Rsp.mat',s,HH(h)))).mean_brain_data;
        brain_data = brain_data(Vertex_now,:);
        roi_data = [roi_data; brain_data];
        voxel_info_1 = [voxel_info_1, 0*ones([1, vertex_num])]; % which ROI, zero is EVC
        voxel_info_2 = [voxel_info_2, h*ones([1, vertex_num])]; % which hemisphere
        voxel_info_3 = [voxel_info_3, s*ones([1, vertex_num])]; % which subject
    end
end

roi_data = double(roi_data)./300;
roi_data = zscore(roi_data,0,2);
%%
[s,s_name] = load_embedding;
n_space = length(s);
r = zeros([length(voxel_info_3),n_space]);
r2 = zeros([length(voxel_info_3),n_space]);
%%
start_parfor
tic
parfor bb = 1:length(voxel_info_3)
    [r(bb,:), r2(bb,:)] = compare_encoders(s,roi_data(bb,:), 10);
end
toc
%%
save_dir=fullfile(prep_dir,'fMRI_result.mat');
save(save_dir,'r','r2','voxel_info_1','voxel_info_2','voxel_info_3','s_name')