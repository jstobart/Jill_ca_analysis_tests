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