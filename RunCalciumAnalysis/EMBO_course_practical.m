
%% make a mixing matrix for spectral unmixing of RCaMP and GCaMP

mixedImg = SCIM_Tif();

[unmixedImg, RCaMP_GCaMP_Matrix] = unmix_chs(mixedImg); %returns the mixing matrix to be used for all imaging

%% calcium analysis

% frame by frame- FLIKA

%load data images
ImgArray =  SCIM_Tif();


%spectral unmix
ImgArray= ImgArray.unmix_chs(false, [], cell2mat(RCaMP_GCaMP_Matrix));
%ImgArray.plot()

% motion correct 
refImg = squeeze(mean(ImgArray(1,1).rawdata(:,:,2, 5:10),4));

%2D convolution
ImgArray_MC1=ImgArray.motion_correct( 'refImg', refImg,'ch', 2,'maxShift', 10,'minCorr', 0.4);
%ImgArray_MC1.plot()


% load configs
BL_frames=50;
configCS=loadConfigCS(BL_frames);



% find ROIs
% FLIKA
FLIKA_2D_cytoA = CellScan([], ImgArray, configCS{1,1}, 1);  % cyto astrocyte 2D FLIKA
FLIKA_2D_cytoA.process();
FLIKA_2D_cytoA.plot();

FLIKA_2D_membA = CellScan([], ImgArray, configCS{1,2}, 1);  % membrane astrocyte 2D FLIKA
FLIKA_2D_membA.process();
FLIKA_2D_membA.plot();

FLIKA_2DN = CellScan([], ImgArray, configCS{1,4}, 2);  % neuronal 2D FLIKA
FLIKA_2DN.process();
FLIKA_2DN.plot();

% 
Hand_A = CellScan([], ImgArray, configCS{1,3}, 2);  % astrocyte hand selected
Hand_A.process();
Hand_A.plot();

Hand_N = CellScan([], ImgArray, configCS{1,5}, 2);  % neuronal hand selected
Hand_N.process();
Hand_N.plot();

% FLIKA3D = CellScan([], ImgArray, configCS{1,6}, 1);
% FLIKA3D.process();
% FLIKA3D.plot();
% FLIKA3D.plot('video');


%% Blood flow- diameter, velocity, flux

Velocity=  LineScanVel();
Velocity.process();
Velocity.plot();

Diameter= LineScanDiam();
Diameter.process();
Diameter.plot();
