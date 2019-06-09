function [cfun zoom pixelSize gof] = microscopecalibration(type,varargin)

%MICROSCOPECALIBRATION Generates calibration curve for two-photon scope
%
%   CFUN = MICROSCOPECALIBRATION(TYPE) produces the calibration curve CFUN
%   (a curve fitting toolbox cfit object) for the calibration TYPE, where
%   TYPE is a string defining a particular calibration.
%
%   CFUN = MICROSCOPECALIBRATION(TYPE,DOPLOT) explicitly states whether or
%   not to display a figure of the calibration curve.  DOPLOT must be a 
%   logical scalar.  If unspecified, the default value is true.
%
%   CFUN = MICROSCOPECALIBRATION(TYPE,DOPLOT,FRAMESIZE) specifies which
%   frame size to produce the calibration curve for, where more than one
%   option exists.  FRAMESIZE must be a numeric scalar for which there is a
%   valid calibration curve.  If unspecified, the default value is 256.
%
%   [CFUN ZOOM] = MICROSCOPECALIBRATION(...) and 
%   [CFUN ZOOM PIXELSIZE] = MICROSCOPECALIBRATION(...) ouput the raw values
%   of ZOOM and PIXELSIZE which were used to produce the calibration curve
%
%   [CFUN ZOOM PIXELSIZE GOF] = MICROSCOPECALIBRATION(...) ouputs the 
%   structure GOF, which contains information about the goodness-of-fit of
%   the calibration.  See the curve fitting toolbox for further
%   documentation.

%   Last updated 15 Feb 2012
%
%   Written by Matthew Barrett
%   Auckland Bioengineering Institute
%   University of Auckland
%   New Zealand
%   matthewjpbarrett@gmail.com

% =========================================================================
% %%%%%%%%%%%%%%%%%%%%%% Process/validate arguments %%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================

DEFAULT_DOPLOT = true;

% Assign input arguments
if nargin > 1
    doPlot = varargin{1};
    if nargin>2
        frameSize = varargin{2};
    end
end

% Check type is assigned correctly
validType = exist('type','var') && ~isempty(type) && ischar(type);
if ~validType
    while true
        strInput = sprintf('\tEmpty or invalid calibration.\n\t\tWhich lens was used (16 or 20)? ');
        inputLens = input(strInput);
        if isempty(inputLens), inputLens = 0; end
        switch inputLens
            case 16
                type = 'renaud__16x__28_01_2011';
                fprintf('\t\tAssuming %s\n',type)
                break
            case 20
                type = 'matthias__20x__10_02_2012';
                fprintf('\t\tAssuming %s\n',type)
                break
            otherwise
                fprintf('\t\tInvalid lens.  Enter "16" or "20"\n');
        end
    end
end

% Check doPlot is assigned correctly
validDoPlot = exist('doPlot','var') && ~isempty(doPlot) && ...
                islogical(doPlot) && isscalar(doPlot);
if ~validDoPlot
    doPlot = DEFAULT_DOPLOT;
end

% Create function to check frameSize is assigned correctly.  Do the actual
% checking only where relevant
function checkframesize(VALID_SIZES,DEFAULT_FRAMESIZE)
    
    % Check for validity
    validFrameSize = exist('frameSize','var') && ~isempty(frameSize) && ...
                        isnumeric(frameSize) && isscalar(frameSize) && ...
                        ismember(frameSize,VALID_SIZES);
    if ~validFrameSize
        frameSize = DEFAULT_FRAMESIZE;
    end

end


% =========================================================================
% %%%%%%%%%%%%%%%%%%%%%%% Type-specific parameters %%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================

switch type
    case 'brunoOriginal__16x__unknown'
        
        % Ensure we're calling for a valid frameSize
        validSizes = 256;
        defaultFrameSize = validSizes(1);
        checkframesize(validSizes,defaultFrameSize)
        
        % WARNING: Don't actually know what frameSize this was captured at,
        % but it's a pretty safe assumption that it was 256 because of the
        % way it's used in the pixelSize calculation below
        
        % Data from calibration
        data = [

        6 107240 107003
        8 107213 107030
        10 107196 107049
        12 107184 107059
        14 107176 107059
        16 107169 107076
        18 107162 107081
        20 107159 107086

        ];

        % Format data
        zoom = data(:,1);
        pixelSize = ( data(:,2) - data(:,3) ) / 256;
        
    case 'renaud__16x__28_01_2011'
        
        % WARNING: this calibration is artificially constructed from the
        % information sent below.  The ZOOM, PIXELSIZE, and GOF information 
        % is meaningless, although the calibration curve itself is valid.
        
        % 16x water immersion, checked that the ruler was horizontal and 
        % the objective vertical. Ruler is aligned on the x-axis direction.
        % I imaged pollen grains with 256 x 256 pixels.
        % 
        % I got: 
        % 
        % pixel size [um] = P1/(zoom factor)+P2
        %
        % with  P1 = 5.31979 pm 0.02848, 
        %       P2 = 0.03464 pm 0.00937, and 
        %       R2 = 0.99966.
        
        % Ensure we're calling for a valid frameSize
        validSizes = 256;
        defaultFrameSize = validSizes(1);
        checkframesize(validSizes,defaultFrameSize)
        
        % artificially constructed from what Renaud sent
        zoom = (20:-2:2)';
        
        % params
        P1 = 5.31979;
        P2 = 0.03464;
        
        % setup fit function
        funPixelSize = @(zoomFactor) P1./zoomFactor + P2;
        
        % create artificial data
        pixelSize = funPixelSize(zoom)*16/20/1.1;
        
    case 'matthew__20x__13_05_2011'
        
        % Ensure we're calling for a valid frameSize
        validSizes = [128 256 512];
        defaultFrameSize = validSizes(2);
        checkframesize(validSizes,defaultFrameSize)
        
        switch frameSize
            case 128, row = 1;
            case 256, row = 2;
            case 512, row = 3;
            otherwise, % unknown frameSize
        end   
        
        % The range of zoom values I calibrated for.
        zoom = (20:-2:2)';
        
        % How far I moved the stage at each zoom level (to maximize the
        % number of pixels the pollen had moved).  This was calculated by
        % taking the difference between the stage positions displayed on 
        % the labview control panel for the microscope.
        stageShift = [  27 30 34 36 49 63 88 127 203 432    ]; % microns
        
        % How many pixels the grain of pollen moved (calculated in a
        % seperate matlab function using cross-correlation).
        pixelShiftArray = ...
                     [  60 67 74 78 85 91 100 106 111 117   ;     % 128^2
                        120 134 148 156 171 183 200 212 222 235 ; % 256^2
                        242 268 298 313 344 367 402 426 447 472]; % 512^2       
                    
        % The pixel size for a given zoom.
        pixelSize = (stageShift./pixelShiftArray(row,:))';
        
    case 'matthias__20x__10_02_2012'
        
        % Ensure we're calling for a valid frameSize
        validSizes = 256;
        defaultFrameSize = validSizes(1);
        checkframesize(validSizes,defaultFrameSize)
        
        % The range of zoom values I calibrated for.
        zoom = (2:2:20)';
        
        % How far the stage moved at each zoom level
        stageShift = [ 122 120 113 76 59 44 34 27 22 14 ]; % microns
        
        % How many pixels the grain of pollen moved (calculated in a
        % seperate matlab function using cross-correlation).
        pixelShiftArray = [ 64 127 179 159 153 133 113 105 100 88 ]; % pixels
                    
        % The pixel size for a given zoom.
        pixelSize = (stageShift./pixelShiftArray)';
     case 'Johannes_ScanImage_20x__09_03_2012'
        
        % Ensure we're calling for a valid frameSize
        validSizes = 256;
        defaultFrameSize = validSizes(1);
        checkframesize(validSizes,defaultFrameSize)
        
        % The range of zoom values I calibrated for.
        zoom = [1,2,3,4,5,8,11]';
        
        % How far the stage moved at each zoom level
        stageShift = [60,60,40,40,40,20,10]; % microns

        % How many pixels the grain of pollen moved (calculated in a
        % seperate matlab function using cross-correlation).
        pixelShiftArray = [34,67,66,83,107,90,46]; % pixels
                    
        % The pixel size for a given zoom.
        pixelSize = (stageShift./pixelShiftArray)';    
    case 'Johannes_ScanImage_20x__12_03_2012'
        
        % Ensure we're calling for a valid frameSize
        validSizes = 256;
        defaultFrameSize = validSizes(1);
        checkframesize(validSizes,defaultFrameSize)
        
        % The range of zoom values I calibrated for.
        zoom = [1,2,3,4,6,8,10,12,16,18,20]';
        
        % How far the stage moved at each zoom level
        stageShift = [60    60    30    60    60    40    40    40    20    20    10]; % microns
        
        % How many pixels the grain of pollen moved (calculated in a
        % seperate matlab function using cross-correlation).
        pixelShiftArray = [16    32    24    63    95    83   105   123    82    89    50]; % pixels
        
        % The pixel size for a given zoom.
        pixelSize = (stageShift./pixelShiftArray)';
        
    otherwise
        
        % Unknown calibration type
        fprintf('\tUnknown calibration type "%s". Trying again.\n',type)
        varargout = microscopecalibration([],varargin);
        return
        
end

% =========================================================================
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Peform fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================

% Setup fit
options = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0 0],...
               'Upper',[10, 2*min(pixelSize)],...
               'Startpoint',[1, mean(pixelSize)]);

% Use a hyperbolic shape.  This function is superior to y = a*exp(bx)
ffun = fittype('a/x + b');

% Perform fit
[cfun,gof] = fit(zoom, pixelSize, ffun, options);

% =========================================================================
% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Peform plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =========================================================================

if doPlot
    figure, hold on
    title(['Frame Size ' num2str(frameSize) 'x' num2str(frameSize)])
    plot(cfun,'b-',zoom,pixelSize,'ro','predfunc')
    xlabel('Zoom'), ylabel('Pixelsize (micrometers)');
    hold off
end

end

