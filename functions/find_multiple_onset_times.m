function tOnset = find_multiple_onset_times(tt, trace, nSDs, nMA, varargin)

% Determine the threshold each trace must exceed
blFrames = find(tt < 0, 1, 'last');
blStd = std(trace(1:blFrames, :), 0, 1);
blThresh = blStd.*nSDs;

% Filter the traces with a centered moving average
tracesFilt = utils.moving_average(trace, nMA, 1);

% Find where the traces exceed the threshold
maskOver = bsxfun(@gt, tracesFilt, blThresh);

% Setup for the loop
%nROIs = size(traces, 2);
%idxOnset = nan(nROIs, 1);
%tOnset = idxOnset;

% Find the first index after the baseline that exceeds the threshold
%for iROI = 1:nROIs
    try %#ok<TRYNC>
        idxOnset = find(maskOver(blFrames+1:end, 1));
        if isempty(idxOnset)
            idxOnset=nan(1,1);
        end
    end
%end

% Add back the baseline frames to the index
idxOnset = idxOnset+blFrames;

% find first of consecutive numbers 
Idiff = @(idxOnset) [1; find(diff(idxOnset)-1)+1; length(idxOnset)];

q=Idiff(idxOnset); % index of first consecutive number
idxOnset2=idxOnset(q);  % this gives all the pairs of consecutive numbers
idxOnset2=idxOnset2(1:end-1); %ignore the redundant last entry

%find numbers that are close together (oscillating peaks that cross the
%threshold with ~400 ms)
i1 = 1;
C{i1}=idxOnset2(i1);
for j1 = 2:numel(idxOnset2)
    t = idxOnset2(j1)-idxOnset2(j1-1);
    if t <= 5
        C{i1} = idxOnset2(j1-1);
    else
        i1  = i1 + 1;
        C{i1} = idxOnset2(j1);
    end
end

idxOnset3=cell2mat(C);

% Work out the actual onset time
maskOK = ~isnan(idxOnset3);
% p=find(diff(a)==1)
% q=[p;p+1];
% a(q)  % this gives all the pairs of consecutive numbers

tOnset(maskOK) = tt(idxOnset3(maskOK));

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
    plot(tt, trace(:, plotROI), 'Color', [0.7, 0.7, 0.7], 'LineWidth', 0.5)
    plot(tt, tracesFilt(:, plotROI), 'k', 'LineWidth', 2)
    plot(tt([1, end])', blThresh(plotROI).*ones(1, 2), 'k--')
%     if ~isnan(idxOnset(plotROI))
%         for iPeak=1:length(idxOnset)
%         plot(tt(idxOnset(plotROI,iPeak)), tracesFilt(idxOnset(plotROI,iPeak), plotROI), 'r+')
%         end
%     end
    hold off
    axis tight
end

end
