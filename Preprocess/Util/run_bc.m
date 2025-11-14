function [qMetric, unitType]=run_bc(data_path)
npx_data = dir('NPX_*');
ephysKilosortPath = fullfile(data_path,'kilosort_def_5block_97');
npx_probe_data = dir(fullfile(data_path,npx_data.name,"*imec0"));
ephysRawDir = dir(fullfile(data_path, npx_data.name, npx_probe_data.name, "*ap*.*bin"));
ephysMetaDir = dir(fullfile(data_path, npx_data.name, npx_probe_data.name, "*ap*.*meta"));
savePath = fullfile(data_path, "processed","BC"); 
mkdir(savePath)
decompressDataLocal = '';
gain_to_uV = NaN;
kilosortVersion = 4;

%% load data

[spikeTimes_samples, spikeTemplates, templateWaveforms, templateAmplitudes, pcFeatures, ...
    pcFeatureIdx, channelPositions] = bc.load.loadEphysData(ephysKilosortPath,savePath);
%% detect whether data is compressed, decompress locally if necessary

rawFile = bc.dcomp.manageDataCompression(ephysRawDir, decompressDataLocal);
%% which quality metric parameters to extract and thresholds
param = bc.qm.qualityParamValues(ephysMetaDir, rawFile, ephysKilosortPath, gain_to_uV, kilosortVersion); %for unitmatch, run this:

%% compute quality metrics
[qMetric, unitType] = bc.qm.runAllQualityMetrics(param, spikeTimes_samples, spikeTemplates, ...
    templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, savePath);
end