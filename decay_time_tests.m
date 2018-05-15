nROI=1:length(All_traces); 

for iROI=1:nROI   
    trace=All_traces{iROI,8};
    FrameRate=All_traces{iROI,11};
    nframes=length(trace);
    TimeX(1:nframes) = (1:nframes)/FrameRate;
    BL_time=All_traces{iROI, 7}/All_traces{iROI,11};  % baseline time
    baselineCorrectedTime=TimeX-BL_time;  % time vector with negative numbers for baseline
    
    % decay time:  peak time - lowest time of the trace (above the
    % threshold)
    DecayFrames=find_decay_time(baselineCorrectedTime(10:end), trace(10:end),2.5,2,1);
    
    if isempty(DecayFrames)
        Decay=nan(1,1);
    else
        Decay=DecayFrames/FrameRate;   % time in s
    end
    Shortstim{iROI, 17}= Decay;   % add to table
 end
    
