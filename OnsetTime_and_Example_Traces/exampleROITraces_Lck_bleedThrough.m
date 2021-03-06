%%
% load necessary example CellScans

stimlength = 1;
%stimwindow=round(50*11.84);
stimwindow=round(35*11.84);


% % Lck LongStim Trial 1
% % extract traces,
shortstim_neur1=CSArray_Ch2_Hand(1,3).calcMeasureROIs.data.tracesNorm(:,3);
shortstim_proc2=CSArray_Ch1_Hand(1,3).calcMeasureROIs.data.tracesNorm(:,3);

shortstim_ACtraces=shortstim_proc2;
shortstim_Ntraces=shortstim_neur1;
%% extract ROI masks
% % long stim

maps_neur=CSArray_Ch2_Hand(1,3).calcFindROIs.data.roiMask(:,:,3);


NeuroMap2=sum(double(maps_neur),3);
NeuroMask=im2bw(NeuroMap2);
Neuro_B=bwboundaries(NeuroMask);

figure();
imshow(zeros(128,128)); hold on
for k=1:length(Neuro_B)
    border=Neuro_B{k};
    plot(border(:,2),border(:,1),'m','linewidth',1.5);
end



%%
grey = [0.8,0.8,0.8];

TimeX=(1:592)/11.84;

figure('name','GCaMP and RCaMP examples')
hold on
axis off
for ii=1:size(shortstim_ACtraces,2)
    tempY1=smooth(shortstim_ACtraces(:,ii),3);
    tempY1=tempY1(1:stimwindow);
plot(TimeX(1:stimwindow),tempY1'+(3*(ii-1)),'g')%'LineWidth',1);
   
end
for ii=1:size(shortstim_Ntraces,2)
    tempY2=smooth(shortstim_Ntraces(:,ii),3);
   tempY2=tempY2(1:stimwindow);

plot(TimeX(1:stimwindow),tempY2'+(3*(ii-1)),'r')
rectangle('Position', [5 -1 1 8.5])
%plot([-1 -1],[-1 1], 'k','LineWidth', 1)
%plot([-5 -1],[-1 -1], 'k','LineWidth', 1)
end
%%
% figure('name','specific')
% hold on
% axis off
% plot(TimeX(1:stimwindow),smooth(DeltaF_R(1:stimwindow,6),5),'r')
% plot(TimeX(1:stimwindow),(smooth(DeltaF_R(1:stimwindow,8),5))+1,'r')
% plot(TimeX(1:stimwindow),(smooth(DeltaF_G(1:stimwindow,3),5))+3,'g')
% plot(TimeX(1:stimwindow),(smooth(DeltaF_G(1:stimwindow,1),5))+5,'g')
% plot(TimeX(1:stimwindow),(smooth(DeltaF_G(1:stimwindow,4),5))+7,'g')
% plot([5 13],[-1 -1], 'k','LineWidth', 2)
% plot([-1 -1],[0 1], 'k','LineWidth', 1)
% 
% 
% 
% 
% 
% %% Output ROIs for examples
%  Image1=zeros(128,128);
%  %Image1(data_traces{68,10})=1;
% %Image1(extraEarly{67,10})=1;
% Image1(earlyGC{5,10})=1;
%  Image1=im2bw(Image1);
%  figure;imshow(Image1)
% 
% %%  
% filename='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\ExampleTraces';
% 
% write_tiff_stacks(CSArray_Ch1_FLIKA, filename)
% 
% 
% %% to look for bleed through
% BL=1:59;
% 
% for iG=1:size(GC,2)
%     base=GC(BL,iG);
%     meanbase=mean(base);
%     %deltaF/F
%     DeltaF_G(:,iG)=(GC(:,iG)-meanbase)/meanbase;
% end
% 
% for iR=1:size(RC,2)
%     base=RC(BL,iR);
%     meanbase=mean(base);
%     %deltaF/F
%     DeltaF_R(:,iR)=(RC(:,iR)-meanbase)/meanbase;
% end
% 
% 
% %%
% TimeX=(1:1065)/11.84;
% stimwindow=round(30*11.84);
% figure('name','GCaMP and RCaMP examples')
% hold on
% axis off
% for ii=1:size(GC,2)
%     tempY1=smooth(DeltaF_G(:,ii),5);
%     tempY1=tempY1(1:stimwindow);
% plot(TimeX(1:stimwindow),tempY1'+(0.5*(ii-1)),'g')%'LineWidth',1);
%    
% end
% for ii=16:20 %ii=1:size(RC,2)
%     tempY2=smooth(DeltaF_R(:,ii),5);
%    tempY2=tempY2(1:stimwindow);
% 
% plot(TimeX(1:stimwindow),tempY2'+(0.5*(ii-1)),'r')
% plot([5 5],[-1 5], 'k--','LineWidth', 0.5)
% plot([5 13],[-1 -1], 'k','LineWidth', 2)
% end
% 
% figure('name','specific')
% hold on
% axis off
% plot(TimeX(1:stimwindow),smooth(DeltaF_R(1:stimwindow,6),5),'r')
% plot(TimeX(1:stimwindow),(smooth(DeltaF_R(1:stimwindow,8),5))+1,'r')
% plot(TimeX(1:stimwindow),(smooth(DeltaF_G(1:stimwindow,3),5))+3,'g')
% plot(TimeX(1:stimwindow),(smooth(DeltaF_G(1:stimwindow,1),5))+5,'g')
% plot(TimeX(1:stimwindow),(smooth(DeltaF_G(1:stimwindow,4),5))+7,'g')
% plot([5 13],[-1 -1], 'k','LineWidth', 2)
% plot([-1 -1],[0 1], 'k','LineWidth', 1)
% 
% 
