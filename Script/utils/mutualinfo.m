function MI = mutualinfo(x, y,nBins, ifplot)
% Calculate the mutual information between x and y using the simple histogram method.
% Jeffrey-Perks law is assumed (adding 0.5 to all cells before MI estimation).
% 'nBins' : number of bins for each dimension. 
% mutual informaiton: I(X;Y) = sum(sum(p(x,y)*log(p(x,y)/p(x)p(y)))
% See "Information Theory" (Shannon, 1948)
% 
% This code is inspired by:
%     Iskarous, K., Mooshammer, C., Hoole, P., Recasens, 
%     D., Shadle, C. H., Saltzman, E., and Whalen, D. H. (2013). 
%     "The coarticulation/invariance scale: Mutual information 
%     as a measure of coarticulation resistance, motor synergy,
%     and articulatory invariance," J. Acoust. Soc. Am. 134, 1271-1282.
% 
% Requirement: "hist3" in Statistics and Machine Learning Toolbox  
% Weirong Chen   March-15-2015

if nargin<3 || isempty(nBins), nBins=6;end;
if nargin<4 || isempty(ifplot), ifplot=0;end;
[bins,binCenters] =hist3([x y],[3 nBins]); % Collecting bins by histogram method. 
bins=bins+0.5; % Jeffrey-Perks law (adding 0.5 to all bins);
totalCount=sum(sum(bins)); 
pxy=bins/totalCount; % pxy = p(x,y) =  measured joint probability.
px=sum(bins,2)/totalCount;  % px = p(x) 
py=sum(bins,1)/totalCount;  % py = p(y)
p_x_y=repmat(px,1,length(py)) .* repmat(py,length(px),1); % p_x_y = p(x)*p(y)   = joint probalility assuming p(x) and p(y) are independent.
MI=sum(sum( pxy .* log2  (pxy ./ p_x_y) )); % Mutual information matrix before summation
if ifplot, close all; myplotpdf(bins, pxy, p_x_y,binCenters); end;
end % end of main function

function H=myplotpdf(bins, pxy, p_x_y,binCenters)
xticks=roundn(binCenters{1},-1); yticks=roundn(binCenters{2},-1); 
H=figure; 
subplot(1,3,1);
bar3( bins);xlabel('y'); ylabel('x');zlabel('Count');
set(gca,'xticklabel',yticks, 'yticklabel',xticks);
subplot(1,3,2);
bar3(pxy);xlabel('y'); ylabel('x');zlabel('Joint Prob');
set(gca,'xticklabel',yticks, 'yticklabel',xticks);
subplot(1,3,3);
bar3(p_x_y);xlabel('y'); ylabel('x');zlabel('Ind Joint Prob');
set(gca,'xticklabel',yticks, 'yticklabel',xticks);
end %function myplotpdf