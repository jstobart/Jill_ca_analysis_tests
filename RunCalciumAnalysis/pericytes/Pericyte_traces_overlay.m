% plot all pericyte ROI traces overlaid on the same graph

% load the traces file saved from the CellScan
AllTraces=load();

FrameRate=11.84;

Grouped_traces=[];

figure ('name', 'Overlaid traces: All pericyte ROIs')
hold on
axis off
for xROI= 1:size(AllTraces,2)
    trace=AllTraces{xROI,10};
    nframes=length(trace);
    TimeX(1:nframes) = (1:nframes)/FrameRate;

    grey = [0.8,0.8,0.8];
    plot(TimeX,trace,'Color',grey,'LineWidth',0.01);
    Grouped_traces= horzcat(Grouped_traces, trace);
    end
MeanTrace=mean(Grouped_traces');
plot(TimeX, Grouped_traces, 'Color', 'k','LineWidth',1);
