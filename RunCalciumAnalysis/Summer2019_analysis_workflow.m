%% WORKFLOW

% RCaMP mice

% All images: unmix RCaMP and FITC dextran
ImgArray= ImgArray.unmix_chs(false, [], cell2mat(RCaMP_mGCaMP_Matrix));
ImgArray = ImgArray.motion_correct('refImg',refImg, 'ch',2,'maxShift', 10,'minCorr', 0.4);%,'doPlot',true);


% measure velocity line scan
% measure diameter line scan
% measure calcium from diameter line scan
% measure calcium from frame scan
%  - motion correction
%  - hand selected ROIs
%  - FLIKA ROIs


%  flikaTest.rawImg.plot()
