%%
% load necessary example CellScans

% % extract traces,
stim_traces=Data.calcMeasureROIs.data.tracesNorm(:,3);

%% extract ROI masks
% % long stim

maps=Data.calcFindROIs.data.roiMask(:,:,3);


NeuroMap2=sum(double(maps),3);
NeuroMask=im2bw(NeuroMap2);
Neuro_B=bwboundaries(NeuroMask);

figure();
imshow(zeros(128,128)); hold on
for k=1:length(Neuro_B)
    border=Neuro_B{k};
    plot(border(:,2),border(:,1),'m','linewidth',1.5);
end



%%

TimeX=(1:592)/11.84;

figure('name','GCaMP examples')
hold on
axis off
for ii=1:size(stim_traces,2)
    tempY1=smooth(stim_traces(:,ii),3);
plot(TimeX,tempY1'+(3*(ii-1)),'g')%'LineWidth',1);
   
end
