function tDecay = find_decay_time(tt, traces, nSDs, nMA, varargin)

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
idxLast = nan(nROIs, 1);
tDecay = idxLast;

% Find the first index after the baseline that exceeds the threshold
for iROI = 1:nROIs
    try %#ok<TRYNC>

        % find the crossing point of the peak as it drops
        ii = [0, diff(maskOver(blFrames+1:end, iROI)')==0,0];
        i1 = strfind(ii,[0 1]);  % peaks start (frame number)
        i2 = strfind(ii,[1 0]);  % peaks end (frame number)
        out = [i1(:),i2(:)];
        idxLast(iROI) = out (2,2);
        
        %find the peak max
        FirstpeakTrace= tracesFilt(out(2,1):out(2,2),iROI);
        maxId=find(FirstpeakTrace==max(FirstpeakTrace));
        
        maxIdx(iROI)=maxId+out(2,1);
    end
end

% Add back the baseline frames to the index
idxLast = idxLast+blFrames;
maxIdx = maxIdx+blFrames;

% number of frames during decay
tDecay=idxLast-maxIdx;

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
    if ~isnan(idxLast(plotROI))
        plot(tt(idxLast(plotROI)), tracesFilt(idxLast(plotROI), plotROI), 'r+')
    end
    hold off
    axis tight
end

end
