imPath = 'J:\_Group\Projects\Astrocyte Calcium\Current Milestones\IP3 knockouts\Imaging\IPRG2\2017_11_16\spot1_linescan_Stim\lowres_spot1_linescan139.tif';
calPath = 'J:\_Group\Software\CHIPS\calibration_20x.mat';

%load E:\matlab\ca-analysis\Jill_ca_analysis_tests\RunCalciumAnalysis\ConfigCellScanLS1D.mat

load E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\RCaMP_mGCaMP_Matrix_4Ch.mat

calibration = CalibrationPixelSize.load(calPath);
channels = struct('Ca_Memb_Astro', 1, 'Ca_Neuron', 2);



rawImg = SCIM_Tif(imPath, channels, calibration);

rawImg= rawImg.unmix_chs(false, [], cell2mat(RCaMP_mGCaMP_Matrix));

rawImg.plot();

% Astrocyte channel
load E:\matlab\ca-analysis\Jill_ca_analysis_tests\RunCalciumAnalysis\ConfigCellScanLS1D_LckGC.mat
testCS_Ch1 = CellScan('', rawImg, confObj, 1);

testCS_Ch1.process

testCS_Ch1.plot

%testCS_Ch1.opt_config()

% neuronal channel
testCS_Ch2 = CellScan('', rawImg, confObj, 2);

testCS_Ch2.process

testCS_Ch2.plot()

%testCS_Ch2.opt_config()
