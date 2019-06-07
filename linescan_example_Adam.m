
%% Velocity Line Scan

% Load the velocity image using the new helper function
rid_LSV = load_prairie_line();

% Create the LineScanVel using this object (which should have the correct
% line time etc)
testLSV = LineScanVel('', rid_LSV);

% Do other things as normal...
testLSV.process();
testLSV.plot();
testLSV.opt_config()

% This is the parameter set I used to produce the figure I sent you
colsToUseLSV = [33, 77];
confVel_Radon = ConfigVelocityRadon('windowTime', 250, 'thresholdSNR', 2);
testLSV_matt = LineScanVel('', rid_LSV, confVel_Radon, [], colsToUseLSV);
testLSV_matt.process();
testLSV_matt.plot();

% The LSPIV approach also seems to work quite well in this case
confVel_LSPIV = ConfigVelocityLSPIV('windowTime', 250, 'shiftAmt', 10);
testLSV_matt2 = LineScanVel('', rid_LSV, confVel_LSPIV, [], colsToUseLSV);
testLSV_matt2.process();
testLSV_matt2.plot();

%% Diameter Line Scan

% Load the diameter image using the new helper function
rid_LSD = load_prairie_line();

% Create the LineScanDiam using this object (which should have the correct
% line time etc)
testLSD = LineScanDiam('', rid_LSD);

% Do other things as normal...
testLSD.process();
testLSD.plot();
testLSD.opt_config()

% This is the parameter set I used to produce the figure I sent you
colsToUseLSD = [9, 59];
confDiam = ConfigDiameterFWHM('maxRate', 5);
testLSD_matt = LineScanDiam('', rid_LSD, confDiam, colsToUseLSD);
testLSD_matt.process();
testLSD_matt.plot();

%% Frame Scan

% Create the frame scan image using the normal BioFormats approach (or also
% as previously using FrameScan directly)
ri_FS = BioFormats();
testFS = FrameScan('', ri_FS);

% Do other things as normal...
testFS.process();
testFS.plot();
testFS.opt_config()

% This is the parameter set I used to produce the figure I sent you
colsToUseFS_Vel = [27, 98];
rowsToUseFS_Vel = [40, 63];
colsToUseFS_Diam = [2, 93];
confFS = ConfigFrameScan(ConfigVelocityRadon());
testFS_matt = FrameScan('', ri_FS, confFS, [], colsToUseFS_Vel, ...
    rowsToUseFS_Vel, colsToUseFS_Diam);
testFS_matt.process();
testFS_matt.plot();
