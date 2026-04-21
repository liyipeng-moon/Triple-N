function CheckEnv()
% Checks whether required MATLAB toolboxes and external functions are installed
missing_toolboxes = {};
required_toolboxes = {
    'Signal Processing Toolbox'
    'Optimization Toolbox'
    'Image Processing Toolbox'
    'Parallel Computing Toolbox'
    };
v = ver;
installed_names = {v.Name};
for i = 1:length(required_toolboxes)
    tb = required_toolboxes{i};
    if any(strcmpi(tb, installed_names))
        fprintf('%s is installed.\n', tb);
    else
        fprintf('%s is NOT installed.\n', tb);
        missing_toolboxes{end+1} = tb;
    end
end

items = {
    'colormap_matplotlib'
    'shadedErrorBar'
    'readNPY'
    'violinplot'
    };
for i = 1:size(items, 1)
    func_name = items{i, 1};
    if exist(func_name, 'file') == 2
        fprintf('%s is available.\n', func_name);
    else
        fprintf('%s is NOT available.\n', func_name);
        missing_toolboxes{end+1} = func_name;
    end
end

if ~isempty(missing_toolboxes)
    fprintf('\n️ Missing required toolboxes:\n');
    for i = 1:length(missing_toolboxes)
        fprintf(' - %s\n', missing_toolboxes{i});
    end
    error('Please install the missing toolboxes before running the code.');
else
    fprintf('\n All required toolboxes are available!\n');
end

fprintf('Looking for other data...\n')
load DIRS.mat
ROI_data = dir(fullfile(FMRI_DIR,'ROI_data.mat'));
if(isempty(ROI_data))
    warning('Extract /others/fmri.zip and paste its path in main code! \n')
else
    fprintf('FMRI path good! \n')
end

FeatureData = dir(fullfile(Model_Dir,'alexnet_layer_rsp.mat'));
if(isempty(FeatureData))
    warning('Extract /others/ModelFeature.zip and paste its path in main code! \n')
else
    fprintf('Model Feature path good! \n')
end
end