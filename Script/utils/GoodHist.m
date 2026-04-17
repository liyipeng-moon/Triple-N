function [binc,cc]=GoodHist(data,edge_here,color_here)
% Computes histogram counts of data using specified bin edges (edge_here) 
% and visualizes them as a step plot with the given color. 

[cc,ee] = histcounts(data,edge_here,Normalization="count");
cc = [0,cc,0];
binc = ee(1:end-1);
binc = [binc(1),binc,binc(end)];
stairs(binc,cc,'LineWidth',1.5,Color=color_here)