%% nomr
ki = 0;
deltaX = [];
zooms      = [1 2 3 4 6 8 10 12 16 18 20];
movement   = [60 60 30 60 60 40 40 40 20 20 10];
for i= zooms
    ki = ki + 1;
    fnames = dir(['zoom0' num2str(i) '*']);
    if i > 9, fnames = dir(['zoom' num2str(i) '*']);end
    A = imread(fnames(1).name);
    B = imread(fnames(2).name);
    C = xcorr2(A,B);
    [mx, idx] = max(C(:));
    deltaX(ki) = abs(idx-floor(idx/511)*511-256);
end
%%

umProPixel = movement./deltaX;
figure,
plot(zooms, umProPixel, 'LineWidth',2);
title('256x256, 2ms linetime, 20x, Calibration ScanImage: 12.03.2012')
ylabel('um pro pixel')
xlabel('zomm')
saveas(gcf,'CalibrationScanImage','fig');
saveas(gcf,'CalibrationScanImage','png');
%%
figure,
plot(zooms, umProPixel*256, 'LineWidth',2);
title('256x256, 2ms linetime, 20x, Calibration ScanImage: 12.03.2012')
ylabel('field of view in um')
xlabel('zomm')
saveas(gcf,'FieldOfViewScanImage','fig');
saveas(gcf,'FieldOfViewScanImage','png');
%%
save('DataCalibrationScanImage','umProPixel','zooms','deltaX')

%% 
slider = [0.001 5 10 15 20 25 30 35 40 50 60 70 80 90 100];
HelioPower = [0.35 0.95 2.54 5.1 8.5 12.8 18.0 23.0 30.0 45.0 60 77.0 90 105 120];
ScanImagePower = [0.35 0.77 2.2 4.6 7.8 12.0 17.0 22.5 28.2 41.0 55 70 78 85 88];
figure,
plot(slider,HelioPower,'-r','DisplayName','Helio','LineWidth',2); hold on
plot(slider,ScanImagePower,'-b','DisplayName','ScanImage','LineWidth',2); 
title('Power over Slider position');
xlabel('Slider position in %');
ylabel('Power in a.u. (mW)');
legend
saveas(gcf,'PowerComparison','fig');
saveas(gcf,'PowerComparison','png');
