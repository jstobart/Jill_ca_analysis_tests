function tOnset = find_first_onset_time(tt, traces, nSDs, nMA, varargin)

% tt= time series (negative numbers are the baseline, time=0 is the onset
% of stimulus, etc)

% traces= the dF/F trace for each ROI (from the CellScan)

% nSDs= number of standard deviations from the baseline for the threshold
% (I used 2.5 for my astrocyte data)

% nMA= the width of the moving average filter for smoothing the traces 

% if varargin is true a plot will be created for each ROI

% Determine the threshold each trace must exceed
blFrames = find(tt < 0, 1, 'last');
blStd = std(traces(1:blFrames, :), 0, 1);
blThresh = blStd.*nSDs;

% Filter the traces with a centered moving average
tracesFilt = utils.moving_average(traces, nMA, 1);

% Find where the traces exceed the threshold
maskOver = bsxfun(@gt, tracesFilt, blThresh);

% Setup for the loop
nROIs = size(traces, 2);
idxOnset = nan(nROIs, 1);
tOnset = idxOnset;

% Find the first index after the baseline that exceeds the threshold
for iROI = 1:nROIs
    try %#ok<TRYNC>
        idxOnset(iROI) = find(maskOver(blFrames+1:end, iROI), 1, 'first');
    end
end

% Add back the baseline frames to the index
idxOnset = idxOnset+blFrames;

% Work out the actual onset time
maskOK = ~isnan(idxOnset);
tOnset(maskOK) = tt(idxOnset(maskOK));

% Plot a debugging figure, if desired
if nargin > 4
    plotROI = varargin{1};
else
    plotROI = [];
end
for ii = 1:numel(plotROI)
    plot_onset(plotROI(ii))
end

% Create a nested function to produce a debugging plot
function plot_onset(plotROI)
    figure, hold on
    plot(tt, traces(:, plotROI), 'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.5)
    plot(tt, tracesFilt(:, plotROI), 'k', 'LineWidth', 2)
    plot(tt([1, end])', blThresh(plotROI).*ones(1, 2), 'k--')
    if ~isnan(idxOnset(plotROI))
        plot(tt(idxOnset(plotROI)), tracesFilt(idxOnset(plotROI), plotROI), 'r+')
    end
    hold off
    axis tight
end

end
