function [contrastMap,luminanceMap] = localContrastMap(img)

gridSize = 50;
outputSize = [400, 400];
padFraction = 0.5;

imgGray = rgb2gray(img);
imgGray = im2double(imgGray);
imgGray = imresize(imgGray, outputSize,'bilinear');

[H, W] = size(imgGray);
gridPixH = H / gridSize;
gridPixW = W / gridSize;
padH = round(gridPixH * padFraction);
padW = round(gridPixW * padFraction);
imgGray = padarray(imgGray, [padH, padW], 0.5, 'both');
imgGray = imgGray.^2;

contrastMap = zeros(gridSize, gridSize);
luminanceMap = zeros(gridSize, gridSize);
for i = 1:gridSize
    for j = 1:gridSize
        rowStart = round((i-1)*gridPixH + 1);
        rowEnd   = round(i*gridPixH + padH*2);
        colStart = round((j-1)*gridPixW + 1);
        colEnd   = round(j*gridPixW + padW*2);
        block = imgGray(rowStart:rowEnd, colStart:colEnd);
        contrastMap(i,j) = std(block(:));
        luminanceMap(i,j) = mean(block(:));
    end
end

end
