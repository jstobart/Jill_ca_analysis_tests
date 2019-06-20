function rid = load_prairie_line(varargin)
%load_prairie_line - Load image data and metadata from Prairie line scan
%
%   OBJ = load_prairie_line() prompts for all required information and
%   loads the image data and metadata from a Prairie Technologies line
%   scan.  OBJ is a RawImgDummy object that can be used like a BioFormats
%   object for all processing in CHIPS.
%
%   OBJ = load_prairie_line(FN_TIF) specifies the filename of the OME-TIF
%   image to load.  The function will attempt to locate the XML file
%   containing the metadata in the same folder.  If it cannot find it (or
%   is confused by multiple XML files, it will prompt to choose the
%   appropriate XML file.
%
%   OBJ = load_prairie_line(FN_TIF, FN_XML) specifies the filename of the
%   XML file containing the relevant metadata.  If the function cannot
%   locate the relevant data in the XML file, it will prompt to specify the
%   required values.
%
%   OBJ = load_prairie_line(FN_TIF, FN_XML, CHS, CAL) specifies the image
%   channels and calibration. CHS and CAL must be scalar and in the format
%   expected by the Metadata constructor (see link below).  If either of
%   these arguments are empty, the constructor will prompt for the required
%   information.
%
%   This function requires the additional function xml2struct, which is
%   available from the <a href="matlab:web('https://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct', '-browser')">MATLAB file exchange(2012)</a>.
%
%   See also BioFormats, RawImgDummy, Metadata, xml2struct

%   Copyright (C) 2017  Matthew J.P. Barrett
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%   
%   You should have received a copy of the GNU General Public License 
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Initial variables + parseing

% Declare the 'guess' filename as a persistent variable so we
% can remember where to start from next time.
persistent fnGuessTif fnGuessXML

% Parse the optional arguments
[fnTif, fnXML, chsIn, calIn] = utils.parse_opt_args(...
    {'', '', [], []}, varargin);

%% Sort out the TIF file

extTif = '.ome.tif';
if isempty(fnTif)
    
    % Choose the tif file, if one is not specified
    filterSpecTif = {['*', extTif], 'Prairie OME-TIF'};
    strTitleTif = 'Select an image file';
    [fTif, pTif] = uigetfile(filterSpecTif, strTitleTif, fnGuessTif);
    fnTif = check_is_cancelled(fTif, pTif, 'TIF');
    fnGuessTif = fnTif;
else
    
    % Otherwise check the supplied image
    utils.checks.equal(fnTif(end-7:end), extTif, ...
        'TIF file extension', extTif)
    utils.checks.file_exists(fnTif);
    
end

%% Sort out the XML file

extXML = '.xml';
if isempty(fnXML)
    
    
    % First, try to find if there's an xml file in the image directory,
    % excluding any XML files that contain 'VoltageOut'
    [pTif, fTif, eTif] = fileparts(fnTif);
    strDir = dir(pTif);
    maskDir = [strDir(:).isdir];
    fnList = {strDir(:).name};
    fnList = fnList(~maskDir);
    maskAllXML = ~cellfun(@isempty, strfind(fnList, extXML));
    maskVoltage = ~cellfun(@isempty, strfind(fnList, 'VoltageOut'));
    maskXML = maskAllXML & ~maskVoltage;
    
    % If that works, good, otherwise prompt the user to choose an XML file
    hasOnlyOne = sum(maskXML) == 1;
    if hasOnlyOne
        fnXML = fullfile(pTif, fnList{maskXML});
    else
        filterSpecXML = {['*', extXML], 'Prairie XML'};
        strTitleXML = sprintf('Select the XML file for %s.%s', fTif, eTif);
        [fXML, pXML] = uigetfile(filterSpecXML, strTitleXML, fnGuessXML);
        fnXML = check_is_cancelled(fXML, pXML, 'XML');
        fnGuessXML = fnXML;
    end
    
else
    
    % Otherwise check the supplied XML file
    utils.checks.equal(fnXML(end-4:end), extXML, ...
        'XML file extension', extXML)
    utils.checks.file_exists(fnXML);
    
end

%% Read in the data from the XML File

% Load the XML file
xmlStruct = xml2struct(fnXML);

% Check the file looks like we expect
try
    strKey = xmlStruct.PVScan.Sequence.Frame{1}.PVStateShard.PVStateValue{...
        4}.Attributes.key;
catch
    strKey = 'ERROR!';
end

% Try to read in the line time from the XML file, otherwise prompt the user
% for the value
isCorrectKey = strcmpi(strKey, 'scanLinePeriod');
if isCorrectKey
    try
        lineTime = str2double(xmlStruct.PVScan.Sequence.Frame{...
            1}.PVStateShard.PVStateValue{4}.Attributes.value)*1000; % ms
    catch
        lineTime = prompt_for_lineTime();
    end
else
    lineTime = prompt_for_lineTime();
end

%% Read in the rest of the image data and create the RawImgDummy object

% Read in the original image using BioFormats
riBF = BioFormats(fnTif, chsIn, calIn);

% Adjust the metadata to the correct lineTime
acq = riBF.metadata.get_acq();
acq.lineTime = lineTime;
acq.pixelTime = lineTime/size(riBF.rawdata, 2);

% Create the new RawImgDummy object!
<<<<<<< HEAD
rid = RawImgDummy(riBF.name, riBF.rawdata,...
=======
rid = RawImgDummy(riBF.name, riBF.rawdata,... 
>>>>>>> d65e00cda651017ee9b33ab750acfee87c330068
    riBF.metadata.channels, riBF.metadata.calibration, acq);

end

% ----------------------------------------------------------------------- %

function fnFull = check_is_cancelled(fn, pn, ftype)

% Throw an error if the user cancelled, otherwise return the filename
hasCancelled = ~(ischar(fn) || iscell(fn)) && (fn == 0) && (pn == 0);
if ~hasCancelled
    fnFull = fullfile(pn, fn);
else
    error('LoadPrairieLine:DidNotChooseFile', ['You must select ' ...
        'one or more %s files to load.'], ftype)
end
            
end

% ----------------------------------------------------------------------- %

function lineTime = prompt_for_lineTime()

% Give a warning
warning('LoadPrairieLine:NoLineTime', ['The line time could not be ' ...
        'automatically extracted.'])

% Prompt the user for the lineTime
while true
    lineTime = input('Please enter the line time [in ms]: ', 's');
    ME = utils.checks.prfs(str2double(lineTime), 'lineTime');
    isOk = isempty(ME);
    if isOk
        break
    else
        warning('LoadPrairieLine:BadLineTime', ['The lineTime ' ...
            'must be a positive real finite scalar number.'])
    end
end
    
end
