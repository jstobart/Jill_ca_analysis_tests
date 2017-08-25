function distance=minDistance(Mask1, Mask2) %calculates the closest distance between the edges of 2 ROIs

Mask = Mask1 + Mask2;
Mask=im2bw(Mask);

%Pythaogrean theorem method

% Define object boundaries
boundaries = bwboundaries(Mask);
numberOfBoundaries = size(boundaries, 1);
if numberOfBoundaries==1
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
