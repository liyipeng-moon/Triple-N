function MTX = Organize_Trial_Data(Ras,Idx)
% Reorganizes trial-wise neural raster data into a standardized 3D matrix
% (unit × image × repetition) by equalizing trial counts across images
% , and applies z-score normalization per unit. 
% 
% This is used for tuning raster data  format used for GSN analysis.

rep = 999;
for img = 1:1000
    rep = min(rep, sum(Idx==img));
end

MTX = zeros([size(Ras,1), 1000, rep]);
for img = 1:1000
    LOC = find(Idx==img);
    LOC = LOC(1:rep);
    MTX(:,img,:)=Ras(:,LOC);
end

for unit = 1:size(MTX,1)
    d = MTX(unit,:,:);
    m = mean(d(:));
    e = std(d(:));
    MTX(unit,:,:) = (d-m)./e;
end