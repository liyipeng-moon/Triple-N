function [channels,wf_near_site] = prune_wf(input_wf)
amp = max(input_wf')-min(input_wf');
[~,peak_channel] = max(amp);
steps = -6:2:6;
channels = peak_channel+steps;
channels = intersect(channels, 1:384);
wf_near_site = input_wf(channels, :);
