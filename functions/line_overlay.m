function outMovie = line_overlay(imPath, startX, endX, startY, endY, col)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Read image
im = imread(imPath);
nFrames = endX-startX;

%% Allocate memory
outMovie = repmat(im, 1, 1, 1,  nFrames);

%% Color
if ~exist('col') || isempty(col)
    col = [255, 0, 0];
end

%% Add line
for iFrame = 1:endX-startX
    for i = 1:3
        outMovie(startY:endY, (startX+iFrame-1), i, iFrame) = col(i);
    end
end


end

