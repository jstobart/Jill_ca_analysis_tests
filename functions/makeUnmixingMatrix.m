calibration ='E:\matlab\CalibrationFiles\calibration_20x.mat';
CalFile = CalibrationPixelSize.load(calibration);

% load example high res image
    unmixFile = 'E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\WT_LR1\2018_01_25_atropine_Kirk\spot4_Stim\highres_spot4062.tif';
    
    channel = struct('Ca_Memb_Astro',1,'Ca_Neuron',2);
    
    unmixImg = SCIM_Tif(unmixFile, channel, CalFile);
    
    [unmixImg, RCaMP_mGCaMP_Matrix] = unmix_chs(unmixImg); %returns the mixing matrix to be used for all imaging
    
    % save the mixing matrix for loading later
    cd(fullfile(Settings.MainDir, 'Results'));
    % write matrix to created file
    save('RCaMP_mGCaMP_Matrix_2KCh.mat', 'RCaMP_mGCaMP_Matrix');