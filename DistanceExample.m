 
% Determine the distance between 2 ROIs


% step 1: extract the ROI mask for each ROI from the CellScan

% For FLIKA ROIs: we need the pixel indices and we will calculate the ROI
% shape 

% For ImageJ ROIs: the Mask is saved as a logical, so we must convert it to
% a double


% EXAMPLE
% ExampleCS is the name of the CellScan

% create 2 masks:

% if ROI 1  is a FLIKA ROI
numROI=1;  % the ROI number
Mask1= zeros(128,128);   % The size of your images!  
Mask1(ExampleCS.calcFindROIs.data.roiIdxs{numROI,1})=1;  

% OR

% if ROI 1 is an ImageJ ROI
Mask1= double(ExampleCS.calcFindROIs.data.roiMask(:,:,numROI));



% if ROI 2  is a FLIKA ROI
numROI=2;  % the ROI number
Mask2= zeros(128,128);   % The size of your images!  
Mask2(ExampleCS.calcFindROIs.data.roiIdxs{numROI,1})=1;  

% OR

% if ROI 2 is an ImageJ ROI
Mask2= double(ExampleCS.calcFindROIs.data.roiMask(:,:,numROI));


% run the distance function
distance= minDistance(Mask1, Mask2);


