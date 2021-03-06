% all ROIs  correlation
clearvars

doplot=0;

%FileSave{1,1}='E:\Data\Two_Photon_Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\LongStim_Correlations.mat';
%FileSave{1,2}='E:\Data\Two_Photon_Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\LongStim_Correlations.csv';

FileSave{1,1}='D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LongStim_Correlations.mat';
FileSave{1,2}='D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LongStim_Correlations.csv';


% load data traces
load('D:\Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\S&LStim_cGC&RC_traces_17_02_2017.mat');

Long_Short=All_traces;

load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\LStim_cGC&RC_traces_17_02_2017.mat');

Long=All_traces;

All_traces=vertcat(Long_Short,Long);

FrameRate=11.84;

searchString='Stim';

for iROI=1:length(All_traces)
    Stim_str{iROI} = strfind(All_traces{iROI,6},searchString);
    str_idx(iROI)=~isempty(Stim_str{iROI});
    
    %find ROITypes
    N_str= strfind(All_traces{iROI, 1},'N');
    S_str= strfind(All_traces{iROI, 1},'S');
    EF_str= strfind(All_traces{iROI, 1},'E');
    P_str= strfind(All_traces{iROI, 1},'r');
    NP_str= strfind(All_traces{iROI, 1},'np');
    
    if ~isempty(N_str)
        All_traces{iROI,13}='Neuron';
    elseif ~isempty(S_str)
        All_traces{iROI,13}='Soma';
    elseif ~isempty(EF_str)
        All_traces{iROI,13}='Endfeet';
    elseif ~isempty(P_str)
        All_traces{iROI,13}='Process';
    elseif ~isempty(NP_str)
        All_traces{iROI,13}='Neuropil';
    end
    
    % make new names
    All_traces{iROI,14}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,1});
    All_traces{iROI,15}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,1});
    All_traces{iROI,16}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,1},'_',All_traces{iROI,6});
    All_traces{iROI,17}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,6});
    
end

data_traces = All_traces(str_idx',:);

for iROI=1:length(data_traces)
    %if a process ROI is overlapping with a soma or process, exclude it
    %from data
    %All_traces{iROI,11}=0;
    nonOverlapIdx(iROI)=~ischar(data_traces{iROI,12});
    
end

% remove overlapping processes
data_traces = data_traces(nonOverlapIdx',:);

%% Correlation of each ROI with all ROIs in the same trial
% CONSIDER WHOLE TRIAL and 25 sec stim window

% only compare ROIs from the same trial together
trial_condition=unique(data_traces(:,17));
All_Corrs=[];
for itrial= 1:length(trial_condition)
    CurrentTrial=trial_condition{itrial};
    
    for iTrace= 1:size(data_traces,1)
        TrialIdx(iTrace)=strcmp(data_traces{iTrace,17}, CurrentTrial);
    end
    
    CurrentData=data_traces(TrialIdx,:);
    
    for iROI=1:size(CurrentData,1)
        CorrData(1,1:6)=CurrentData(iROI,2:7);
        CorrData(1,7)=CurrentData(iROI,17); 
        
        CorrData(1,8)= CurrentData(iROI,1);
        CorrData(1,9)= CurrentData(iROI,13);
        CorrData(1,10)= CurrentData(iROI,9);
        CorrData(1,11)= CurrentData(iROI,10);
        
        %whole trace
        TraceX= CurrentData{1,8};
        % stim window
        TraceX_short=TraceX(1:(25*FrameRate),:);
        
        for kROI = (iROI+1):size(CurrentData,1)
            CorrData(1,12)= CurrentData(kROI,3);
            CorrData(1,13)= CurrentData(kROI,1);
            CorrData(1,14)= CurrentData(kROI,13);
            CorrData(1,15)= CurrentData(kROI,9);
            CorrData(1,16)= CurrentData(kROI,10);
            
            %whole trace
            TraceY= CurrentData{kROI,8};
            % stim window
            TraceY_short=TraceY(1:(25*FrameRate),:);
            
            % whole trial correlation
            [R,Ps] = corrcoef(TraceX, TraceY);
            CorrData{1,17}= R(1,2); % correlation coeff
            CorrData{1,18}= Ps(1,2); % p value
            
            % stim window correlation
            [R2,Ps2]= corrcoef(TraceX_short, TraceY_short);
            CorrData{1,19}= R2(1,2); % correlation coeff
            CorrData{1,20}= Ps2(1,2); % p value
            
            % stim window cross correlation
            [Corr_Sequence,lag]  = xcorr(TraceX_short, TraceY_short,'coeff');
            [~,I] = max(abs(Corr_Sequence));
            CorrData{1,21}=  max(abs(Corr_Sequence)); % max of correlation sequence
            CorrData{1,22}= (lag(I))/FrameRate; % lag time
            
            % correlation plots
            %             subplot(3,1,1); plot(TraceX_short); title('s1');
            %             subplot(3,1,2); plot(TraceY_short); title('s2');
            %             subplot(3,1,3); plot(lag,Corr_Sequence);
            %             title('Cross-correlation between s1 and s2')
            
            %% Distance Calculations
            % find the minimium distance between the edges of the two ROIs
            
            %create a binary image with each pair of ROIs
            % for the first ROI
            if isnumeric(CorrData{1,11})
                Image1=zeros(128,128);
                Image1(CorrData{1,11})=1;
                %Image1=im2bw(Image1);
            elseif islogical(CorrData{1,11})
                Image1= double(CorrData{1,11});
            else
                Image1=[];
            end
            
            % for the second ROI
            if isnumeric(CorrData{1,16})
                Image2=zeros(128,128);
                Image2(CorrData{1,16})=1;
                Image2=im2bw(Image2);
            elseif islogical(CorrData{1,16})
                Image2= double(CorrData{1,16});
            else
                Image2=[];
            end
            
            Mask=Image1+Image2;
            Mask=im2bw(Mask);
            
            %Pythaogrean theorem method
            
            % Define object boundaries
            boundaries = bwboundaries(Mask);
            numberOfBoundaries = size(boundaries, 1);
            if numberOfBoundaries==1
                CorrData{1,23} = 0;
            elseif numberOfBoundaries>1
                boundary1 = boundaries{1};
                boundary2 = boundaries{2};
                boundary1x = boundary1(:, 2);
                boundary1y = boundary1(:, 1);
                for k = 1 : length(boundary2)
                    boundary2x = boundary2(k, 2);
                    boundary2y = boundary2(k, 1);
                    % For this blob, compute distances from boundaries to edge.
                    allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2);
                    % Find closest point, min distance.
                    [minDistance(k), indexOfMin] = min(allDistances);
                end
                % Find the overall min distance
                CorrData{1,23} = (min(minDistance)*cell2mat(CurrentData(kROI,11)));
            else
                CorrData{1,23} = [];
            end
            CorrData{1,24} = numberOfBoundaries; % number of ROIs in the mask
            
            All_Corrs=vertcat(All_Corrs,CorrData); % concatenate all data into a big matrix
            
            clear boundaries boundary1 boundary2 boundary1x boundary1y boundary2x boundary2y minDistance indexofMin
        end
        
    end
    
end




%% save as a CSV
names = {'Trial','ChannelX','Spot','Animal','Condition','Depth','TrialName','ROI_X',...
    'ROI_X_type','ChannelY','ROI_Y','ROI_Y_type','Long_Corr','LongPvalue',...
    'Short_Corr','ShortPvalue','xCorr','Lag','MinDistance','numROIs'};

CorrelationData=All_Corrs;
CorrelationData(:,10) = [];
CorrelationData(:,10) = [];
CorrelationData(:,13) = [];
CorrelationData(:,13) = [];

AllData=[names;CorrelationData];

% write date to created file
cell2csv(FileSave{1,2}, AllData);
save(FileSave{1,1}, 'All_Corrs','-v7.3');

%%

if doplot
    figure()
    title('Whole trial', 'FontSize', 14); % set title
    hold on; axis off
    %for iCond = 1:mCond
    ActTrials= [];
    for kTrial = 1:mtrials
        % mean correlation plot across trials- Stim
        % find "active" trials
        %                     ActiveTrials = find(AllTrials(iCond,:,kTrial)>0);
        %if ~isempty(ActiveTrials)  % trials with at least 1 peak in first 30 sec
        temp = cell2mat(Long_Rcor(:,:,kTrial,iCond));
        subplot((mtrials+1),mCond,(mCond*(kTrial-1))+iCond)
        imagesc(squeeze(temp));
        set(gca, 'XTick', 1:numROIs); % center x-axis ticks on bins
        set(gca, 'YTick', 1:numROIs); % center y-axis ticks on bins
        set(gca, 'XTickLabel', ROI_stuff(:,1)); % set x-axis labels
        set(gca, 'YTickLabel', ROI_stuff(:,1)); % set y-axis labels
        colormap('jet'); % set the colorscheme
        colorbar; % enable colorbar
        %                         ActTrials = cat(3,ActTrials,temp);
        %                     %end
    end
    % mean correlation plot across trials- NO Stim
    %tempMean= mean(ActTrials,3);
    hold on
    subplot((mtrials+1),mCond,(mCond*(mtrials))+iCond)
    imagesc(squeeze(Mean_Rcor(:,:,iCond)));
    set(gca, 'XTick', 1:numROIs); % center x-axis ticks on bins
    set(gca, 'YTick', 1:numROIs); % center y-axis ticks on bins
    set(gca, 'XTickLabel', ROI_stuff(:,1)); % set x-axis labels
    set(gca, 'YTickLabel', ROI_stuff(:,1)); % set y-axis labels
    title(strcat('Mean Correlation Matrix',Spots{iSpot}), 'FontSize', 14); % set title
    colormap('jet'); % set the colorscheme
    colorbar; % enable colorbar
    % end
end

%% DSP4

% all ROIs  correlation
clearvars

doplot=0;

FileSave{1,1}='E:\Data\Two_Photon_Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\DSP4_LongStim_Correlations.mat';
FileSave{1,2}='E:\Data\Two_Photon_Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\DSP4_LongStim_Correlations.csv';

% load data traces
load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\DSP4_cGC&RC_traces_17_02_2017.mat');

FrameRate=11.84;

searchString='Stim';

for iROI=1:length(All_traces)
    Stim_str{iROI} = strfind(All_traces{iROI,6},searchString);
    str_idx(iROI)=~isempty(Stim_str{iROI});
    
    %find ROITypes
    N_str= strfind(All_traces{iROI, 1},'N');
    S_str= strfind(All_traces{iROI, 1},'S');
    EF_str= strfind(All_traces{iROI, 1},'E');
    P_str= strfind(All_traces{iROI, 1},'r');
    NP_str= strfind(All_traces{iROI, 1},'np');
    
    if ~isempty(N_str)
        All_traces{iROI,13}='Neuron';
    elseif ~isempty(S_str)
        All_traces{iROI,13}='Soma';
    elseif ~isempty(EF_str)
        All_traces{iROI,13}='Endfeet';
    elseif ~isempty(P_str)
        All_traces{iROI,13}='Process';
    elseif ~isempty(NP_str)
        All_traces{iROI,13}='Neuropil';
    end
    
    % make new names
    All_traces{iROI,14}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,1});
    All_traces{iROI,15}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,1});
    All_traces{iROI,16}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,1},'_',All_traces{iROI,6});
    All_traces{iROI,17}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,6});
    
end

data_traces = All_traces(str_idx',:);

for iROI=1:length(data_traces)
    %if a process ROI is overlapping with a soma or process, exclude it
    %from data
    %All_traces{iROI,11}=0;
    nonOverlapIdx(iROI)=~ischar(data_traces{iROI,12});
    
end

% remove overlapping processes
data_traces = data_traces(nonOverlapIdx',:);

%% Correlation of each ROI with all ROIs in the same trial
% CONSIDER WHOLE TRIAL and 25 sec stim window

% only compare ROIs from the same trial together
trial_condition=unique(data_traces(:,17));
All_Corrs=[];
for itrial= 1:length(trial_condition)
    CurrentTrial=trial_condition{itrial};
    
    for iTrace= 1:size(data_traces,1)
        TrialIdx(iTrace)=strcmp(data_traces{iTrace,17}, CurrentTrial);
    end
    
    CurrentData=data_traces(TrialIdx,:);
    
    for iROI=1:size(CurrentData,1)
        CorrData(1,1:6)=CurrentData(iROI,2:7);
        CorrData(1,7)=CurrentData(iROI,17);
        
        CorrData(1,8)= CurrentData(iROI,1);
        CorrData(1,9)= CurrentData(iROI,13);
        CorrData(1,10)= CurrentData(iROI,9);
        CorrData(1,11)= CurrentData(iROI,10);
        
        %whole trace
        TraceX= CurrentData{1,8};
        % stim window
        TraceX_short=TraceX(1:(25*FrameRate),:);
        
        for kROI = (iROI+1):size(CurrentData,1)
            CorrData(1,12)= CurrentData(kROI,3);
            CorrData(1,13)= CurrentData(kROI,1);
            CorrData(1,14)= CurrentData(kROI,13);
            CorrData(1,15)= CurrentData(kROI,9);
            CorrData(1,16)= CurrentData(kROI,10);
            
            %whole trace
            TraceY= CurrentData{kROI,8};
            % stim window
            TraceY_short=TraceY(1:(25*FrameRate),:);
            
            % whole trial correlation
            [R,Ps] = corrcoef(TraceX, TraceY);
            CorrData{1,17}= R(1,2); % correlation coeff
            CorrData{1,18}= Ps(1,2); % p value
            
            % stim window correlation
            [R2,Ps2]= corrcoef(TraceX_short, TraceY_short);
            CorrData{1,19}= R2(1,2); % correlation coeff
            CorrData{1,20}= Ps2(1,2); % p value
            
            % stim window cross correlation
            [Corr_Sequence,lag]  = xcorr(TraceX_short, TraceY_short,'coeff');
            [~,I] = max(abs(Corr_Sequence));
            CorrData{1,21}=  max(abs(Corr_Sequence)); % max of correlation sequence
            CorrData{1,22}= (lag(I))/FrameRate; % lag time
            
            % correlation plots
            %             subplot(3,1,1); plot(TraceX_short); title('s1');
            %             subplot(3,1,2); plot(TraceY_short); title('s2');
            %             subplot(3,1,3); plot(lag,Corr_Sequence);
            %             title('Cross-correlation between s1 and s2')
            
            %% Distance Calculations
            % find the minimium distance between the edges of the two ROIs
            
            %create a binary image with each pair of ROIs
            % for the first ROI
            if isnumeric(CorrData{1,11})
                Image1=zeros(128,128);
                Image1(CorrData{1,11})=1;
                %Image1=im2bw(Image1);
            elseif islogical(CorrData{1,11})
                Image1= double(CorrData{1,11});
            else
                Image1=[];
            end
            
            % for the second ROI
            if isnumeric(CorrData{1,16})
                Image2=zeros(128,128);
                Image2(CorrData{1,16})=1;
                Image2=im2bw(Image2);
            elseif islogical(CorrData{1,16})
                Image2= double(CorrData{1,16});
            else
                Image2=[];
            end
            
            Mask=Image1+Image2;
            Mask=im2bw(Mask);
            
            %Pythaogrean theorem method
            
            % Define object boundaries
            boundaries = bwboundaries(Mask);
            numberOfBoundaries = size(boundaries, 1);
            if numberOfBoundaries==1
                CorrData{1,23} = 0;
            elseif numberOfBoundaries>1
                boundary1 = boundaries{1};
                boundary2 = boundaries{2};
                boundary1x = boundary1(:, 2);
                boundary1y = boundary1(:, 1);
                for k = 1 : length(boundary2)
                    boundary2x = boundary2(k, 2);
                    boundary2y = boundary2(k, 1);
                    % For this blob, compute distances from boundaries to edge.
                    allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2);
                    % Find closest point, min distance.
                    [minDistance(k), indexOfMin] = min(allDistances);
                end
                % Find the overall min distance
                CorrData{1,23} = (min(minDistance)*cell2mat(CurrentData(kROI,11)));
            else
                CorrData{1,23} = [];
            end
            CorrData{1,24} = numberOfBoundaries; % number of ROIs in the mask
            
            All_Corrs=vertcat(All_Corrs,CorrData); % concatenate all data into a big matrix
            
            clear boundaries boundary1 boundary2 boundary1x boundary1y boundary2x boundary2y minDistance indexofMin
        end
        
    end
    
end


%% save as a CSV
names = {'Trial','ChannelX','Spot','Animal','Condition','Depth','TrialName','ROI_X',...
    'ROI_X_type','ChannelY','ROI_Y','ROI_Y_type','Long_Corr','LongPvalue',...
    'Short_Corr','ShortPvalue','xCorr','Lag','MinDistance','numROIs'};

CorrelationData=All_Corrs;
CorrelationData(:,10) = [];
CorrelationData(:,10) = [];
CorrelationData(:,13) = [];
CorrelationData(:,13) = [];

AllData=[names;CorrelationData];

% write date to created file
cell2csv(FileSave{1,2}, AllData);
save(FileSave{1,1}, 'All_Corrs','-v7.3');

%% Lck

% all ROIs  correlation
clearvars

doplot=0;

FileSave{1,1}='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LongStim_Correlations.mat';
FileSave{1,2}='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LongStim_Correlations.csv';

% load data traces
load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\S&LStim_LckGC&RC_traces_17_02_2017.mat');

Long_Short=All_traces;

load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LStim_LckGC&RC_traces_17_02_2017.mat');

Long=All_traces;

All_traces=vertcat(Long_Short,Long);

FrameRate=11.84;

searchString='Stim';

for iROI=1:length(All_traces)
    Stim_str{iROI} = strfind(All_traces{iROI,6},searchString);
    str_idx(iROI)=~isempty(Stim_str{iROI});
    
    %find ROITypes
    N_str= strfind(All_traces{iROI, 1},'N');
    S_str= strfind(All_traces{iROI, 1},'S');
    EF_str= strfind(All_traces{iROI, 1},'E');
    P_str= strfind(All_traces{iROI, 1},'r');
    NP_str= strfind(All_traces{iROI, 1},'np');
    
    if ~isempty(N_str)
        All_traces{iROI,13}='Neuron';
    elseif ~isempty(S_str)
        All_traces{iROI,13}='Soma';
    elseif ~isempty(EF_str)
        All_traces{iROI,13}='Endfeet';
    elseif ~isempty(P_str)
        All_traces{iROI,13}='Process';
    elseif ~isempty(NP_str)
        All_traces{iROI,13}='Neuropil';
    end
    
    % make new names
    All_traces{iROI,14}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,1});
    All_traces{iROI,15}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,1});
    All_traces{iROI,16}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,1},'_',All_traces{iROI,6});
    All_traces{iROI,17}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,6});
    
end

data_traces = All_traces(str_idx',:);

for iROI=1:length(data_traces)
    %if a process ROI is overlapping with a soma or process, exclude it
    %from data
    %All_traces{iROI,11}=0;
    nonOverlapIdx(iROI)=~ischar(data_traces{iROI,12});
    
end

% remove overlapping processes
data_traces = data_traces(nonOverlapIdx',:);

%% Correlation of each ROI with all ROIs in the same trial
% CONSIDER WHOLE TRIAL and 25 sec stim window

% only compare ROIs from the same trial together
trial_condition=unique(data_traces(:,17));
All_Corrs=[];
for itrial= 1:length(trial_condition)
    CurrentTrial=trial_condition{itrial};
    
    for iTrace= 1:size(data_traces,1)
        TrialIdx(iTrace)=strcmp(data_traces{iTrace,17}, CurrentTrial);
    end
    
    CurrentData=data_traces(TrialIdx,:);
    
    for iROI=1:size(CurrentData,1)
        CorrData(1,1:6)=CurrentData(iROI,2:7);
        CorrData(1,7)=CurrentData(iROI,17);
        
        CorrData(1,8)= CurrentData(iROI,1);
        CorrData(1,9)= CurrentData(iROI,13);
        CorrData(1,10)= CurrentData(iROI,9);
        CorrData(1,11)= CurrentData(iROI,10);
        
        %whole trace
        TraceX= CurrentData{1,8};
        % stim window
        TraceX_short=TraceX(1:(25*FrameRate),:);
        
        for kROI = (iROI+1):size(CurrentData,1)
            CorrData(1,12)= CurrentData(kROI,3);
            CorrData(1,13)= CurrentData(kROI,1);
            CorrData(1,14)= CurrentData(kROI,13);
            CorrData(1,15)= CurrentData(kROI,9);
            CorrData(1,16)= CurrentData(kROI,10);
            
            %whole trace
            TraceY= CurrentData{kROI,8};
            % stim window
            TraceY_short=TraceY(1:(25*FrameRate),:);
            
            % whole trial correlation
            [R,Ps] = corrcoef(TraceX, TraceY);
            CorrData{1,17}= R(1,2); % correlation coeff
            CorrData{1,18}= Ps(1,2); % p value
            
            % stim window correlation
            [R2,Ps2]= corrcoef(TraceX_short, TraceY_short);
            CorrData{1,19}= R2(1,2); % correlation coeff
            CorrData{1,20}= Ps2(1,2); % p value
            
            % stim window cross correlation
            [Corr_Sequence,lag]  = xcorr(TraceX_short, TraceY_short,'coeff');
            [~,I] = max(abs(Corr_Sequence));
            CorrData{1,21}=  max(abs(Corr_Sequence)); % max of correlation sequence
            CorrData{1,22}= (lag(I))/FrameRate; % lag time
            
            % correlation plots
            %             subplot(3,1,1); plot(TraceX_short); title('s1');
            %             subplot(3,1,2); plot(TraceY_short); title('s2');
            %             subplot(3,1,3); plot(lag,Corr_Sequence);
            %             title('Cross-correlation between s1 and s2')
            
            %% Distance Calculations
            % find the minimium distance between the edges of the two ROIs
            
            %create a binary image with each pair of ROIs
            % for the first ROI
            if isnumeric(CorrData{1,11})
                Image1=zeros(128,128);
                Image1(CorrData{1,11})=1;
                %Image1=im2bw(Image1);
            elseif islogical(CorrData{1,11})
                Image1= double(CorrData{1,11});
            else
                Image1=[];
            end
            
            % for the second ROI
            if isnumeric(CorrData{1,16})
                Image2=zeros(128,128);
                Image2(CorrData{1,16})=1;
                Image2=im2bw(Image2);
            elseif islogical(CorrData{1,16})
                Image2= double(CorrData{1,16});
            else
                Image2=[];
            end
            
            Mask=Image1+Image2;
            Mask=im2bw(Mask);
            
            %Pythaogrean theorem method
            
            % Define object boundaries
            boundaries = bwboundaries(Mask);
            numberOfBoundaries = size(boundaries, 1);
            if numberOfBoundaries==1
                CorrData{1,23} = 0;
            elseif numberOfBoundaries>1
                boundary1 = boundaries{1};
                boundary2 = boundaries{2};
                boundary1x = boundary1(:, 2);
                boundary1y = boundary1(:, 1);
                for k = 1 : length(boundary2)
                    boundary2x = boundary2(k, 2);
                    boundary2y = boundary2(k, 1);
                    % For this blob, compute distances from boundaries to edge.
                    allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2);
                    % Find closest point, min distance.
                    [minDistance(k), indexOfMin] = min(allDistances);
                end
                % Find the overall min distance
                CorrData{1,23} = (min(minDistance)*cell2mat(CurrentData(kROI,11)));
            else
                CorrData{1,23} = [];
            end
            CorrData{1,24} = numberOfBoundaries; % number of ROIs in the mask
            
            All_Corrs=vertcat(All_Corrs,CorrData); % concatenate all data into a big matrix
            
            clear boundaries boundary1 boundary2 boundary1x boundary1y boundary2x boundary2y minDistance indexofMin
        end
        
    end
    
end




%% save as a CSV
names = {'Trial','ChannelX','Spot','Animal','Condition','Depth','TrialName','ROI_X',...
    'ROI_X_type','ChannelY','ROI_Y','ROI_Y_type','Long_Corr','LongPvalue',...
    'Short_Corr','ShortPvalue','xCorr','Lag','MinDistance','numROIs'};

CorrelationData=All_Corrs;
CorrelationData(:,10) = [];
CorrelationData(:,10) = [];
CorrelationData(:,13) = [];
CorrelationData(:,13) = [];

AllData=[names;CorrelationData];

% write date to created file
cell2csv(FileSave{1,2}, AllData);
save(FileSave{1,1}, 'All_Corrs','-v7.3');


%% Lck- short stim

% all ROIs  correlation
clearvars

doplot=0;

FileSave{1,1}='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\ShortStim_Correlations.mat';
FileSave{1,2}='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\ShortStim_Correlations.csv';

% load data traces
load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\S&LStim_LckGC&RC_traces_17_02_2017.mat');

Long_Short=All_traces;

load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LStim_LckGC&RC_traces_17_02_2017.mat');

Long=All_traces;

All_traces=vertcat(Long_Short,Long);

FrameRate=11.84;

searchString='shortstim';

for iROI=1:length(All_traces)
    Stim_str{iROI} = strfind(All_traces{iROI,6},searchString);
    str_idx(iROI)=~isempty(Stim_str{iROI});
    
    %find ROITypes
    N_str= strfind(All_traces{iROI, 1},'N');
    S_str= strfind(All_traces{iROI, 1},'S');
    EF_str= strfind(All_traces{iROI, 1},'E');
    P_str= strfind(All_traces{iROI, 1},'r');
    NP_str= strfind(All_traces{iROI, 1},'np');
    
    if ~isempty(N_str)
        All_traces{iROI,13}='Neuron';
    elseif ~isempty(S_str)
        All_traces{iROI,13}='Soma';
    elseif ~isempty(EF_str)
        All_traces{iROI,13}='Endfeet';
    elseif ~isempty(P_str)
        All_traces{iROI,13}='Process';
    elseif ~isempty(NP_str)
        All_traces{iROI,13}='Neuropil';
    end
    
    % make new names
    All_traces{iROI,14}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,1});
    All_traces{iROI,15}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,1});
    All_traces{iROI,16}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,1},'_',All_traces{iROI,6});
    All_traces{iROI,17}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,6});
    
end

data_traces = All_traces(str_idx',:);

for iROI=1:length(data_traces)
    %if a process ROI is overlapping with a soma or process, exclude it
    %from data
    %All_traces{iROI,11}=0;
    nonOverlapIdx(iROI)=~ischar(data_traces{iROI,12});
    
end

% remove overlapping processes
data_traces = data_traces(nonOverlapIdx',:);

%% Correlation of each ROI with all ROIs in the same trial
% CONSIDER WHOLE TRIAL and 25 sec stim window

% only compare ROIs from the same trial together
trial_condition=unique(data_traces(:,17));
All_Corrs=[];
for itrial= 1:length(trial_condition)
    CurrentTrial=trial_condition{itrial};
    
    for iTrace= 1:size(data_traces,1)
        TrialIdx(iTrace)=strcmp(data_traces{iTrace,17}, CurrentTrial);
    end
    
    CurrentData=data_traces(TrialIdx,:);
    if length(CurrentData{1,8})>200
        
        for iROI=1:size(CurrentData,1)
            CorrData(1,1:6)=CurrentData(iROI,2:7);
            CorrData(1,7)=CurrentData(iROI,17);
            
            CorrData(1,8)= CurrentData(iROI,1);
            CorrData(1,9)= CurrentData(iROI,13);
            CorrData(1,10)= CurrentData(iROI,9);
            CorrData(1,11)= CurrentData(iROI,10);
            
            %whole trace
            TraceX= CurrentData{1,8};
            % stim window
            TraceX_short=TraceX(1:(25*FrameRate),:);
            
            for kROI = (iROI+1):size(CurrentData,1)
                CorrData(1,12)= CurrentData(kROI,3);
                CorrData(1,13)= CurrentData(kROI,1);
                CorrData(1,14)= CurrentData(kROI,13);
                CorrData(1,15)= CurrentData(kROI,9);
                CorrData(1,16)= CurrentData(kROI,10);
                
                %whole trace
                TraceY= CurrentData{kROI,8};
                % stim window
                TraceY_short=TraceY(1:(25*FrameRate),:);
                
                % whole trial correlation
                [R,Ps] = corrcoef(TraceX, TraceY);
                CorrData{1,17}= R(1,2); % correlation coeff
                CorrData{1,18}= Ps(1,2); % p value
                
                % stim window correlation
                [R2,Ps2]= corrcoef(TraceX_short, TraceY_short);
                CorrData{1,19}= R2(1,2); % correlation coeff
                CorrData{1,20}= Ps2(1,2); % p value
                
                % stim window cross correlation
                [Corr_Sequence,lag]  = xcorr(TraceX_short, TraceY_short,'coeff');
                [~,I] = max(abs(Corr_Sequence));
                CorrData{1,21}=  max(abs(Corr_Sequence)); % max of correlation sequence
                CorrData{1,22}= (lag(I))/FrameRate; % lag time
                
                % correlation plots
                %             subplot(3,1,1); plot(TraceX_short); title('s1');
                %             subplot(3,1,2); plot(TraceY_short); title('s2');
                %             subplot(3,1,3); plot(lag,Corr_Sequence);
                %             title('Cross-correlation between s1 and s2')
                
                %% Distance Calculations
                % find the minimium distance between the edges of the two ROIs
                
                %create a binary image with each pair of ROIs
                % for the first ROI
                if isnumeric(CorrData{1,11})
                    Image1=zeros(128,128);
                    Image1(CorrData{1,11})=1;
                    %Image1=im2bw(Image1);
                elseif islogical(CorrData{1,11})
                    Image1= double(CorrData{1,11});
                else
                    Image1=[];
                end
                
                % for the second ROI
                if isnumeric(CorrData{1,16})
                    Image2=zeros(128,128);
                    Image2(CorrData{1,16})=1;
                    Image2=im2bw(Image2);
                elseif islogical(CorrData{1,16})
                    Image2= double(CorrData{1,16});
                else
                    Image2=[];
                end
                
                Mask=Image1+Image2;
                Mask=im2bw(Mask);
                
                %Pythaogrean theorem method
                
                % Define object boundaries
                boundaries = bwboundaries(Mask);
                numberOfBoundaries = size(boundaries, 1);
                if numberOfBoundaries==1
                    CorrData{1,23} = 0;
                elseif numberOfBoundaries>1
                    boundary1 = boundaries{1};
                    boundary2 = boundaries{2};
                    boundary1x = boundary1(:, 2);
                    boundary1y = boundary1(:, 1);
                    for k = 1 : length(boundary2)
                        boundary2x = boundary2(k, 2);
                        boundary2y = boundary2(k, 1);
                        % For this blob, compute distances from boundaries to edge.
                        allDistances = sqrt((boundary1x - boundary2x).^2 + (boundary1y - boundary2y).^2);
                        % Find closest point, min distance.
                        [minDistance(k), indexOfMin] = min(allDistances);
                    end
                    % Find the overall min distance
                    CorrData{1,23} = (min(minDistance)*cell2mat(CurrentData(kROI,11)));
                else
                    CorrData{1,23} = [];
                end
                CorrData{1,24} = numberOfBoundaries; % number of ROIs in the mask
                
                All_Corrs=vertcat(All_Corrs,CorrData); % concatenate all data into a big matrix
                
                clear boundaries boundary1 boundary2 boundary1x boundary1y boundary2x boundary2y minDistance indexofMin
            end
            
        end
    end
    
end




%% save as a CSV
names = {'Trial','ChannelX','Spot','Animal','Condition','Depth','TrialName','ROI_X',...
    'ROI_X_type','ChannelY','ROI_Y','ROI_Y_type','Long_Corr','LongPvalue',...
    'Short_Corr','ShortPvalue','xCorr','Lag','MinDistance','numROIs'};

CorrelationData=All_Corrs;
CorrelationData(:,10) = [];
CorrelationData(:,10) = [];
CorrelationData(:,13) = [];
CorrelationData(:,13) = [];

AllData=[names;CorrelationData];

% write date to created file
cell2csv(FileSave{1,2}, AllData);
save(FileSave{1,1}, 'All_Corrs','-v7.3');