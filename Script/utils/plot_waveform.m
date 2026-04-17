function plot_waveform(input_here,x)
% Plots smoothed neural waveforms for channels with maximal absolute 
% deflection, aligning them in time (tt) and applying Gaussian smoothing
% for clarity.

tt = [1:61]*(1000/30000)+2.1*(x);
[a,b] = find(input_here==max(abs(input_here(:))) | input_here==-max(abs(input_here(:))));
interested_channel = a;
LOCs = 0;

for channs_idx = 1:length(interested_channel)
    channel_now = interested_channel(channs_idx);
    LOC_now = LOCs(channs_idx);
    wf_here = LOC_now+squeeze((input_here(channel_now,:)));
    wf_here = smoothdata(wf_here,1,'gaussian',5);
    plot(tt, wf_here,'k','LineWidth',1.5);
end


