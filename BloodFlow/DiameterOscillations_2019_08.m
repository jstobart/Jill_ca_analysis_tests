Alldata=[]; Alldata2=[]; DataTable=[];
Settings.Drug ='9months';
Settings.Animal ='MJ';
Settings.SaveFile = {'E:\Jill\Data\Winnipeg\Pericytes\In vivo 2P\RCaMP Mice\Results\Bloodflow Results\RCaMP 59721 MJ\Diameter_Oscillations\DiameterOscillationsAnalysis_9months.csv'};

Settings.FileNames = {  % folder names where images and roiSet.zip is found
    'E:\Jill\Data\Winnipeg\Pericytes\In vivo 2P\RCaMP Mice\Results\Bloodflow Results\RCaMP 59721 MJ\Diameter\2019_11_07',...
    %'J:\PericyteData\In vivo 2P\GCaMP Mice\Results\Bloodflow Results\GCaMP 79335\Diameter\2019_12_06',...
    %'J:\PericyteData\In vivo 2P\GCaMP Mice\Results\Bloodflow Results\GCaMP 79335\Diameter\2019_12_06',...
    };


for iDrug= 1:length(Settings.FileNames)
    
    % loop through each Drug
    SpotRoot= Settings.FileNames{iDrug};
     
    % Get a list of all csv in this folder.
    % different trials
    FileNames = dir(SpotRoot);  % look for folders
    FileNames(ismember( {FileNames.name}, {'.', '..'})) = [];
    
    for iFile= 1:length(FileNames)
        % Get the path from the folder
        names   = {FileNames.name};
        file=names{iFile};
        SpotId=file(12:19);  %vessel name

        testRoot =fullfile(SpotRoot, names{iFile});
        
        data= xlsread(testRoot);

        Trace=data(1:150,2);
        time= data(1:150,1);
        %%%
        %figure();
        %plot(time, smooth(BL(1:150,2),15))
        
        [peaks,locs,w,p]=findpeaks(smooth(Trace,15),time,'MinPeakProminence',0.05, 'WidthReference','halfheight');
        
        % mean time between peaks
        meanCycle = mean(diff(locs)); 
        
         % full width at half max
        meanWidth = mean(w);
       
        % mean prominence
        meanProm = mean(p); 
        
        % number of peaks per minute
        freq = (length(peaks)/time(end))*60;  
        
        
        % output the data
        output=[];
        output{1,1}=Settings.Animal;
        output{1,2} =Settings.Drug;
        output{1,3}=SpotId;
        output{1,4}= freq;
        output{1,5} =meanProm;
        output{1,6}=meanWidth;
        output{1,7}=meanCycle;
        
        %All data
         Alldata= vertcat(Alldata,output);
        
         clearvars output
        
    end
Alldata2= vertcat(Alldata2,Alldata);
clearvars Alldata
end

% %% Save all data for R analysis
      % append peak data
        dataNames={'Animal','Drug','Vessel','peaksPerMin','peakProminence','peakWidth','CycleTime'};
        DataTable= vertcat(DataTable,Alldata2);
        DataTable=vertcat(dataNames,DataTable);

cell2csv(Settings.SaveFile, DataTable);