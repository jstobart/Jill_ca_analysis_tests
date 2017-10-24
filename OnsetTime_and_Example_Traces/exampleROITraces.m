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
