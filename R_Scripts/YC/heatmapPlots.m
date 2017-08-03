% import csv file

% find unique ROI names
listROI=unique(beforestimROI(:,3));

colourScale=

figure('name','amplitudes')
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iROI=1:length(listROI)
   ROIIdx=find(~cellfun(@isempty, regexp(beforestimROI(:,3), listROI{iROI})));
   amplitudes=beforestimROI(ROIIdx,9);
   for iTrial=1:length(ROIIdx)
       scatter(iTrial, iROI, 5, 'filled','s')
   end
   
end


for iComp=1:length(AvsN_SpaceOnsetsD_F)
    
    xlim([-15 15])
end
plot([0 0],[0 length(AvsN_SpaceOnsetsD_F)], 'k--','LineWidth', 1)
xlabel('time from neuronal event')
