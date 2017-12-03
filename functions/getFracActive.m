function [signalN, varargout] = getFracActive(cs, varargin)

% Define the allowed optional arguments and default values, and create a
% default parameters structure
pnames = {'backgroundLevel', 'baselineFrames', 'nanMask'};
dflts  = {cs.calcMeasureROIs.config.backgroundLevel, ...
    cs.calcFindROIs.config.baselineFrames, 0};
params = cell2struct(dflts, pnames, 2);

% Parse function input arguments
params = utils.parsepropval(params, varargin{:});

% Create temporary mean image
tempImg = mean(cs.rawImg.rawdata(:,:,cs.channelToUse,params.baselineFrames),4);

% Resize mask if necessary
params.nanMask = imresize(params.nanMask, size(tempImg));

% Remove pixels from count
tempImg(params.nanMask) = NaN;

% Threshold image
signalIdx = (tempImg > prctile(tempImg(:),params.backgroundLevel));

% Count pixels above threshold
signalN = nansum(signalIdx(:));

if nargout > 1
    % Count ROI mask pixels
    activeN = sum(cs.calcFindROIs.data.roiMask(:));
    varargout{1} = activeN;
end

if nargout > 2
    % Total pixels
    totalN = numel(cs.calcFindROIs.data.roiMask(:)) - sum(params.nanMask(:));
    varargout{2} = totalN;
end

end


