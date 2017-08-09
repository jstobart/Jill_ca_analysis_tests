%%
% load necessary example CellScans

stimlength = 8;
%stimwindow=round(50*11.84);
stimwindow=round(25*11.84);

% Fast Lck traces
% long stim RG17, 2016_08_03 trial 8
shortstim_proc1=CSArray_Ch1_FLIKA_2D.calcMeasureROIs.data.tracesNorm(:,9);
shortstim_dend1=CSArray_Ch2_FLIKA_2D.calcMeasureROIs.data.tracesNorm(:,5);
shortstim_neur1=CSArray_Ch2_Hand.calcMeasureROIs.data.tracesNorm(:,3);

shortstim_ACtraces=cat(2,shortstim_proc1);
shortstim_Ntraces=cat(2, shortstim_dend1, shortstim_neur1);


%% extract ROI masks

 % fast Lck long stim
 maps_proc=CSArray_Ch1_FLIKA_2D.calcFindROIs.data.puffIdxs(9,1);
% 
 maps_dend=CSArray_Ch2_FLIKA_2D.calcFindROIs.data.puffIdxs(5,1);
 maps_neur=CSArray_Ch2_Hand.calcFindROIs.data.roiMask(:,:,3);

ProcMap1=zeros(127,128);
for imaps=1:length(maps_proc)
    Image1=zeros(127,128);
    Image1(maps_proc{imaps,1})=1;
    ProcMap1=ProcMap1+Image1;
end
ProcMap2=zeros(1,128);
ProcMaps=[ProcMap1;ProcMap2];

ACMask=im2bw(ProcMaps);
AC_B=bwboundaries(ACMask);

figure();
imshow(zeros(128,128)); hold on
for k=1:length(AC_B)
    border=AC_B{k};
    plot(border(:,2),border(:,1),'g','linewidth',1.5);
end

DendMap1=zeros(127,128);
for imaps=1:length(maps_dend)
    Image1=zeros(127,128);
    Image1(maps_dend{imaps,1})=1;
    DendMap1=DendMap1+Image1;
end
DendMap2=zeros(1,128);
DendMaps=[DendMap1;DendMap2];

DendMask=im2bw(DendMaps);
Dend_B=bwboundaries(DendMask);

figure();
imshow(zeros(128,128)); hold on
for k=1:length(Dend_B)
    border=Dend_B{k};
    plot(border(:,2),border(:,1),'m','linewidth',1.5);
end


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
rectangle('Position', [5 -1 8 20])
plot([-1 -1],[-1 1], 'k','LineWidth', 1)
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
