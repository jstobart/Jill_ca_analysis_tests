function distance=minDistance(Mask1, Mask2) %calculates the closest distance between the edges of 2 ROIs

% add together the 2 ROI masks
Mask = Mask1 + Mask2;

% converts to black and white image
Mask=im2bw(Mask);

% Define boundaries of each ROI
boundaries = bwboundaries(Mask);
numberOfBoundaries = size(boundaries, 1);

%Pythaogrean theorem method
% Calculate the distance of each point on ROI 1 bondary to points on ROI 2
% boundary.  Then determine the minimum (i.e. closest) distance.
if numberOfBoundaries==1  % ROIs are touching
    distance = 0;
elseif numberOfBoundaries>1
    boundary1 = boundaries{1};
    boundary2 = boundaries{2};
    boundary1x = boundary1(:, 2);
    boundary1y = boundary1(:, 1);
    minDistance = [length(boundary2) 1];
    for k = 1 : length(boundary2)
        boundary2x = boundary2(k, 2);
        boundary2y = boundary2(k, 1);
        allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2);
        % Find closest point, min distance.
        [minDistance(k), ~] = min(allDistances);
    end
    % Find the overall min distance
    distance = min(minDistance);
end
