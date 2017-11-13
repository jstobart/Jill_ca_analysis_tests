function [signalN, varargout] = getFracActive(cs, varargin)
    
    % get Config or user-defined values
    if nargin == 1
        backgroundLevel = cs.calcMeasureROIs.config.backgroundLevel;
        baselineFrames = cs.calcFindROIs.config.baselineFrames;
    elseif nargin == 2
        backgroundLevel = varargin{1};
    elseif nargin == 3
        baselineFrames = varargin{2};
    else 
        error('Too many input arguments! Maximum is 3, given are %i', nargin)
    end
    
    % Create temporary mean image
    tempImg = mean(cs.rawImg.rawdata(:,:,channelToUse,baselineFrames),4);
    
    % Threshold image
    signalIdx = (tempImg > prctile(tempImg(:),backgroundLevel));
    
    % Count pixels above threshold
    signalN = sum(signalIdx(:));
    
    if nargout > 1
        % Count ROI mask pixels
        activeN = sum(cs.calcFindROIs.data.roiMask(:));
        varargout{1} = activeN;
    end
    
    if nargout > 2
        % Total pixels
        totalN = numel(cs.calcFindROIs.data.roiMask(:));
        varargout{2} = totalN;
    end

end


