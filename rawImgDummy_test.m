%%%
% Specify the full path to the raw image object
fnRID003 = fullfile(utils.CHIPS_rootdir, 'tests', 'res', ...
    'ArtificialData_SNR_0.1_FOV_250.tif');

% Specify the channels relevant for this raw image
channels003 = struct('cellular_signal', 1);

% Specify some data about the image acquisition
nRows003 = 1024;
nCols003 = 1024;
acq003 = struct('isBiDi', false, 'lineTime', 0.5, 'zoom', 0.9, ...
    'nLinesPerFrameOrig', nRows003, 'nPixelsPerLineOrig', nCols003);

% Create the RawImgDummy object without any interaction
rid003 = RawImgDummy(fnRID003, channels003, calibration, acq003);

% View the RawImgDummy object
rid003.plot()
%%%
