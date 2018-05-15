%
% clearvars
% close all

%% Load data

% save files names
saveFiles1='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\time_comparisons.mat';
saveFiles2='E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\time_comparisons.csv';

%peak data
% load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\Lck_GCaMP6f\Results\LckGC&RC_2D_longstim_28_04_2017.csv');
% 


%Peaks=AllData2(2:end,:);
Peaks=datarawshortnoGeph;
for iROI=1:length(Peaks)
    % make new unique trial names
    Peaks{iROI,19}=strcat(Peaks{iROI,18},'_',Peaks{iROI,1});
    % make new Condition/treatment names
    Peaks{iROI,20}=strcat(Peaks{iROI,13},'_',num2str(Peaks{iROI,14}));
end

%% compare peak times for each ROI

if ~exist(saveFiles1, 'file')
    tic
    TimeComparisons=[];
    Trials= unique(Peaks(:,19));
    CondTreat= unique(Peaks(:,20));
    
    % peak comparisons for each field of view and each trial
    for iTreat=1:length(CondTreat)
        CurrentTreat=CondTreat(iTreat);
        
        % Find the idx of matching trials
        matchingTreatIdx = find(~cellfun(@isempty, strfind(Peaks(:,20), CurrentTreat)));
        TreatData = Peaks(matchingTreatIdx,:);
        
        % peak comparisons for each field of view and each trial
        for itrial=1:length(Trials)
            CurrentTrial=Trials(itrial);
            
            % Find the idx of matching trials
            matchingTrialIdx = find(~cellfun(@isempty, str(TreatData(:,19), CurrentTrial)));
            TrialData = TreatData(matchingTrialIdx,:);
            
            numROIs = size(TrialData,1);
            
            for kROI = 1:numROIs
                for iROI = (kROI+1):numROIs
                    %build giant data table
                    temp1{1,1:5} = TrialData(kROI,[13,14,16,17,18]); %condition,treatment,ROIname1, etc.
                    temp1{1,6} = TrialData(kROI,15); % peak time for ROI1
                    temp1{1,7} = TrialData(iROI,16); % ROIname2
                    temp1{1,8} = TrialData(iROI,15); % peak time for ROI2
                    temp1{1,9} = (TrialData(kROI,15)-TrialData(iROI,15)); % peak time difference
                end
            end
            
            TimeComparisons=vertcat(TimeComparisons,temp1); % concatenate all data into a big matrix
            
            clear TrialData
        end
    end
    
    save(saveFiles1, 'TimeComparisons','-v7.3');
    
    names={'condition','treatment','ROI1','areaU','sessionU','peakTime1','ROI2','peakTime2','TimeDiff'};
    TimeComp2=vertcat(names,TimeComparisons);
    cell2csv(saveFiles2, TimeComp2);
    
    toc
else
    load(saveFiles1);
end

%% histograms
figure('name', 'histogram: peakTime time differences')
histogram(cell2mat(TimeComparisons(:,9)))


%% Raster plot

figure('name', 'Raster plot time differences')
hold on
set(gca,'ytick',[])
set(gca,'YColor',get(gcf,'Color'))
for iComp=1:length(TimeComparisons)
    
    scatter(cell2mat(TimeComparisons(iComp,9)), iComp, 5, 'filled','k')
end
plot([0 0],[0 length(TimeComparisons)], 'r--','LineWidth', 1)


