function PostProcess_function_raw(data_path)
   
cd(data_path)

meta_file = dir('processed/META*');
load(fullfile(pwd,'processed',meta_file(1).name));
meta_data = load(fullfile(pwd,'processed',meta_file(1).name));
ML_FILE = dir("processed\ML*");
trial_ML = load(fullfile('processed',ML_FILE(1).name)).trial_ML;

[qMetric, unitType] = run_bc(data_path);
[UnitStrc] = load_KS4_output('./kilosort_def_5block_97',IMEC_AP_META,SyncLine);


for trial_idx = 1:length(trial_ML)
    trial_ML(trial_idx).AnalogData.Mouse=[];
    trial_ML(trial_idx).AnalogData.KeyInput=[];
end

if(strcmp(pwd,'F:\NSD_FULL_RAW\240709'))
    UnitStrc(638)=[];
end
if(strcmp(pwd,'F:\NSD_FULL_RAW\240924'))
    UnitStrc(328)=[];
end

file_name_LOCAL = fullfile('processed',sprintf('GoodUnitRaw_%s_g%s.mat',meta_file(1).name(6:end-4), meta_data.g_number));
save(file_name_LOCAL, "UnitStrc", "trial_ML",'meta_data','qMetric', 'unitType','-v7.3')
end