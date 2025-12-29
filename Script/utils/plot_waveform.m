function plot_waveform(input_here,x)

tt = [1:61]*(1000/30000)+2.1*(x);
[a,b] = find(input_here==max(abs(input_here(:))) | input_here==-max(abs(input_here(:))));

interested_channel = a;
% LOCs = 0.7*max(abs(input_here(:)))*[1:length(interested_channel)];
LOCs = 0;


for channs_idx = 1:length(interested_channel)
    channel_now = interested_channel(channs_idx);
    LOC_now = LOCs(channs_idx);
    wf_here = LOC_now+squeeze((input_here(channel_now,:)));
    xxx = 1;
    wf_here = smoothdata(wf_here,1,'gaussian',5);
    plot(tt, wf_here,'k','LineWidth',1.5);
end



% % axis off
% yl_here = ylim;
% yl_here = yl_here(1);
% plot([tt(5),tt(20)], [yl_here,yl_here],'k','LineWidth',1);
% text(tt(10),yl_here-10,1,'0.5ms')
% plot([tt(5),tt(5)],[yl_here,yl_here+50],'k','LineWidth',1);
% text(tt(2),yl_here,1,'50uv','Rotation',90)
% 
