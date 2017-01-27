%% Plot all traces for the entire trial

% individual traces in grey, mean in red
% plot neurons and astrocytes separately
% plot each ROItype separately (neurons, neuropil, AC processes, somata,
% endfeet)
% plot mean neuron and astrocyte on the same graph (with error bars filled
% in)

nROIs=length(All_traces);
nframes=length(All_traces{1,8}); 

FrameRate = 11.84;
TimeX(1:nframes) = (1:nframes)/FrameRate;


figure ('name', 'traces for entire trial')
hold on
for x = 1:6
    subplot(2, 3, x)
    hold on
    set(gca,'Xlim',[-20 80],'Ylim',[-0.5 5.5],'XTick', 0:20:200)
    xlabel('Time (s)')
    ylabel('dF/F')
    grey = [0.7,0.7,0.7];
    rectangle('Position',[0,-0.5,5,6],'EdgeColor',grey,'FaceColor',grey)
end


figure ('name', 'All ROIs')
hold on

for xROI: 1:nROIs
        tempY = All_traces{xROI,8};NostimTrace2(:,xROI,itrial);
        grey = [0.7,0.7,0.7];
    plot(TimeX,tempY,'Color',grey,'LineWidth',0.5);
    plot([0 60],[18.5 18.5], 'k','LineWidth', 2)
        plot([0 0],[18.5 19.5], 'k','LineWidth', 2)
end

for x = 1:6
    subplot(2, 3, x)
    hold on
    set(gca,'Xlim',[-20 80],'Ylim',[-0.5 5.5],'XTick', 0:20:200)
    xlabel('Time (s)')
    ylabel('dF/F')
    grey = [0.7,0.7,0.7];
    rectangle('Position',[0,-0.5,5,6],'EdgeColor',grey,'FaceColor',grey)
end


figure('name','all traces')
hold on
axis off
        tempY = All_traces{1,8};
    plot(TimeX,tempY,'Color','b','LineWidth',1);


  


%% Plot only the window around the stimulus

% individual traces in grey, mean in red
% plot neurons and astrocytes separately
% plot each ROItype separately
% plot mean neuron and astrocyte on the same graph (with error bars filled
% in)


%% Plot only the responding neurons and astrocytes from the same field of view

%% Plot only responding neurons and responding astrocytes from the same field of view

