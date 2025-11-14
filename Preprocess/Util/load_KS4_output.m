function [ks_output] = load_KS4_output(ks_path, IMEC_AP_META,SyncLine)
loading_data = {'spike_times', 'spike_templates','templates','spike_positions','amplitudes'};
for d_idx = 1:length(loading_data)
    data_now = loading_data{d_idx};
    eval(sprintf("%s = readNPY(fullfile(ks_path, '%s.npy'));",data_now,data_now))
end
spike_times = 1000*double(spike_times)/IMEC_AP_META.imSampRate;
sync_spike_times = interp1(SyncLine.imec_time, SyncLine.NI_time, spike_times, 'linear', 'extrap');

spike_templates = spike_templates+1;

filename = fullfile(ks_path,'cluster_KSLabel.tsv');
fileID = fopen(filename, 'r');
headerLine = fgetl(fileID);
headers = strsplit(headerLine, '\t');
data = textscan(fileID, '%d%s%d', 'Delimiter', '\t', 'HeaderLines',0);
KS_LABEL = data{2};

example_unit.waveform=[];
example_unit.spiketime_ms=[];
example_unit.spikepos=[];
example_unit.amplitudes=[];
strc_unit = repmat(example_unit, [1, max(spike_templates)]);


for spike_idx = 1:max(spike_templates)
    example_unit.waveform = squeeze(templates(spike_idx,:,:));
    example_unit.spiketime_ms = sync_spike_times(spike_templates==spike_idx);
    example_unit.spikepos = mean(spike_positions((spike_templates==spike_idx),:));
    example_unit.amplitudes = amplitudes(spike_templates==spike_idx);
    strc_unit(spike_idx)=example_unit;
    fprintf('Organiza KS output for unit %d %d\n',spike_idx, max(spike_templates))
end
ks_output=strc_unit;
end
