% % First you need to get a binary image with the center filled in. 
% I suggest you run around the border of the image and if there are any white pixels 
% closer than 10 pixels or so, connect them. Then call imfill(binaryImage, 'holes') 
% to get three solid branches.
% 
% Then call bwdist to get the distance transform on the binary image.
% 
% Next call bwmorph(binaryImage, 'skel', inf) to get the skeleton. The skeleton 
% runs all the way out to the border and you don't want that since you'll get distance 
% values of zero so you want to clip the skeleton by a few pixels with the 'spur' option of 
% bwmorph. Think about it and you'll realize why. If not, plot the skeleton over the
% binary image and then you'll see why.
% 
% The distance transform has distances from every point to the border, which is not 
% what we want. We want only those distances along the centerline, the skeleton. So 
% multiply the skeleton by the distance transform to get only the radii. Double it to 
% get the diameters.
% 
% Then extract all non-zero pixels and find the min. Then you're done.


%create a binary image with each pair of ROIs
% for the first ROI
if isdouble(CorrData{1,11})
    Image1=zeros(128,128);
    Image1(CorrData{1,11})=1;
    %Image1=im2bw(Image1);
elseif islogical(CorrData{1,11})
   Image1= double(CorrData{1,11}); 
else
    Image1=[];
end

% for the second ROI
if isdouble(CorrData{1,16})
    Image2=zeros(128,128);
    Image2(CorrData{1,16})=1;
    Image2=im2bw(Image2);
elseif islogical(CorrData{1,16})
   Image2= double(CorrData{1,16}); 
else
    Image2=[];
end

Mask=Image1+Image2;
Mask=im2bw(Mask);

%Pythaogrean theorem method

% Define object boundaries
boundaries = bwboundaries(Mask);
numberOfBoundaries = size(boundaries, 1);
boundary1 = boundaries{1};
boundary2 = boundaries{2};
boundary1x = boundary1(:, 2);
boundary1y = boundary1(:, 1);
for k = 1 : length(boundary2)
    boundary2x = boundary2(k, 2);
    boundary2y = boundary2(k, 1);	
    % For this blob, compute distances from boundaries to edge.
    allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2);
    % Find closest point, min distance.
    [minDistance(k), indexOfMin] = min(allDistances);
end
% Find the overall min distance
CorrData{1,23} = min(minDistance);
