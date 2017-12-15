% plot all pericyte ROI traces overlaid on the same graph

% load the traces file saved from the CellScan
load('E:\Data\Pericyte_project\Two-photon-data\Calcium\Results\test_traces_18_08_2017.mat')

FrameRate=11.84;
nframes=1065;
Grouped_traces=[];

figure ('name', 'Overlaid traces: All pericyte ROIs')
hold on
%axis off
for xROI= 1:size(All_traces,1)
    trace=All_traces{xROI,10}; % each individual trace
    %nframes=length(trace);
    TimeX(1:nframes) = (1:nframes)/FrameRate;  % time series for plotting X
    
    grey = [0.8,0.8,0.8];
    plot(TimeX,trace(1:nframes),'Color',grey,'LineWidth',0.01);  
    Grouped_traces= horzcat(Grouped_traces, trace);
end
MeanTrace=nanmean(Grouped_traces');
plot(TimeX, MeanTrace(1:nframes), 'Color', 'k','LineWidth',1);
