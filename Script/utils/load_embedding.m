function [s,model_nm] = load_embedding
% Loads and organizes feature embeddings from multiple models (LLMs and vision networks)
% applies PCA to reduce dimensionality

model_nm = {};
s={};
load("DIRS.mat")

Model_dir = fullfile(root_dir,'Data','others','ModelFeature','LLM');
all_LLM = dir(fullfile(Model_dir,'*mat'));
for L_idx = 1:length(all_LLM)
    LLM = load(fullfile(Model_dir,all_LLM(L_idx).name));
    llm_ebd = zeros([1000, size(LLM.embeddings,2)]);
    loc = 0;
    for img = 1:1000
        loc = find(LLM.image_id==LLM.image_id(loc+1));
        llm_ebd(img,:) = mean(LLM.embeddings(loc,:),1);
        loc = loc(end);
    end
    [~,s{L_idx},~,~,~,~] = pca(llm_ebd,NumComponents=100);
    model_nm{L_idx} = LLM.Model_name;
end
load(fullfile(root_dir,'Data','others','ModelFeature',"alexnet_layer_rsp.mat"))

for ll = 1:8
    s{end+1} = score{ll};
    model_nm{end+1} = sprintf('Alex%02d',ll);
end

visual_dir = fullfile(root_dir,'Data','others','ModelFeature');
res_name = {'R50_DINO','R50_IN1k'};
for res_case =  1:length(res_name)
    data = dir(fullfile(visual_dir,res_name{res_case},'*layer*'));
    data(end+1) = dir(fullfile(visual_dir,res_name{res_case},'*avg*'));
    data(end+1) = dir(fullfile(visual_dir,res_name{res_case},'*fc*'));
    for ll = 1:length(data)
        s{end+1} = load(fullfile(visual_dir,res_name{res_case},data(ll).name)).feat_pca;
        model_nm{end+1} = sprintf('%s_%02d',res_name{res_case},ll);
    end
end

md_name = 'ViTB16';
data = dir(fullfile(visual_dir,md_name,'*l*'));
for ll = 1:length(data)
    s{end+1} = load(fullfile(visual_dir,md_name,data(ll).name)).feat_pca;
    model_nm{end+1} = sprintf('%s_%02d',md_name,ll);
end

md_name = 'InceptionV3';
data = dir(fullfile(visual_dir,md_name,'*v3*'));
for ll = 1:length(data)
    s{end+1} = load(fullfile(visual_dir,md_name,data(ll).name)).feat_pca;
    model_nm{end+1} = sprintf('%s_%02d',md_name,ll);
end

end

