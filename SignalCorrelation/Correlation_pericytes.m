% all ROIs  correlation
clearvars

FileSave{1,1}='E:\Data\Pericyte_project\Two-photon-data\Calcium\Results\Pericyte_Correlations.csv';

% load data traces
load('E:\Data\Pericyte_project\Two-photon-data\Calcium\Results\test_traces_21_03_2017.mat');

FrameRate=11.84;


for iROI=1:length(All_traces)   
    %find ROITypes
    S_str= strfind(All_traces{iROI, 1},'S');
    P_str= strfind(All_traces{iROI, 1},'P');
    
    if ~isempty(S_str)
        All_traces{iROI,10}='Soma';
    elseif ~isempty(P_str)
        All_traces{iROI,10}='Process';
    end
    
    % make new ROI names
    All_traces{iROI,11}=strcat(All_traces{iROI,1},'_',All_traces{iROI,3},'_',All_traces{iROI,4}); % unique ROI name
   % make new ROI names including trial name
    All_traces{iROI,12}=strcat(All_traces{iROI,1},'_',All_traces{iROI,2},'_',All_traces{iROI,3},'_',All_traces{iROI,4}); % unique ROI name
      % make new trial names
    All_traces{iROI,13}=strcat(All_traces{iROI,2},'_',All_traces{iROI,3},'_',All_traces{iROI,4}); % unique ROI name

end


%% Correlation of each ROI with all ROIs in the same trial
% CONSIDER WHOLE TRIAL 

% only compare ROIs from the same trial together
trials=unique(All_traces(:,13));
All_Corrs=[];

for itrial= 1:length(trials)
    CurrentTrial=trials{itrial};
    
    for iTrace= 1:size(All_traces,1)
        TrialIdx(iTrace)=strcmp(All_traces{iTrace,13}, CurrentTrial);
    end
    
    CurrentData=All_traces(TrialIdx,:);
    
    for iROI=1:size(CurrentData,1)
        CorrData(1,1:5)=CurrentData(iROI,2:6); %info about mouse, trial, etc.
        CorrData(1,6)=CurrentData(iROI,13); % trial name
        
        CorrData(1,7)= CurrentData(iROI,1); % first ROI
        CorrData(1,8)= CurrentData(iROI,10); % ROI type
        CorrData(1,9)= CurrentData(iROI,8); % ROI mask
        
        %whole trace
        TraceX= CurrentData{1,7};
        
        for kROI = (iROI+1):size(CurrentData,1)
            CorrData(1,10)= CurrentData(kROI,1); % second ROI
            CorrData(1,11)= CurrentData(kROI,10); % ROI type 2
            CorrData(1,12)= CurrentData(kROI,8); % ROI mask 2
            
            %whole trace
            TraceY= CurrentData{kROI,7};
            
            % whole trial linear correlation with pearsons correlation
            % coefficient
            [R,Ps] = corrcoef(TraceX, TraceY);
            CorrData{1,13}= R(1,2); % correlation coeff
            CorrData{1,14}= Ps(1,2); % p value
            
            % cross correlation
            [Corr_Sequence,lag]  = xcorr(TraceX, TraceY,'coeff');
            [~,I] = max(abs(Corr_Sequence));
            CorrData{1,15}=  max(abs(Corr_Sequence)); % max of correlation sequence
            CorrData{1,16}= (lag(I))/FrameRate; % lag time
            
            % correlation plots
            %             subplot(3,1,1); plot(TraceX); title('s1');
            %             subplot(3,1,2); plot(TraceY); title('s2');
            %             subplot(3,1,3); plot(lag,Corr_Sequence);
            %             title('Cross-correlation between s1 and s2')
            
            %% Distance Calculations
            % find the minimium distance between the edges of the two ROIs
            
            %create a binary image with each pair of ROIs
            % for the first ROI
            if isnumeric(CorrData{1,9})
                Image1=zeros(128,128);
                Image1(CorrData{1,9})=1;
                %Image1=im2bw(Image1);
            elseif islogical(CorrData{1,9})
                Image1= double(CorrData{1,9});
            else
                Image1=[];
            end
            
            % for the second ROI
            if isnumeric(CorrData{1,12})
                Image2=zeros(128,128);
                Image2(CorrData{1,12})=1;
                Image2=im2bw(Image2);
            elseif islogical(CorrData{1,12})
                Image2= double(CorrData{1,12});
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
                CorrData{1,17} = 0;
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
                % Find the overall min distance in um using the pixel size
                CorrData{1,17} = (min(minDistance)*cell2mat(CurrentData(kROI,9)));
            else
                CorrData{1,17} = [];
            end
            CorrData{1,18} = numberOfBoundaries; % number of ROIs in the mask
            
            All_Corrs=vertcat(All_Corrs,CorrData); % concatenate all data into a big matrix
            
            clear boundaries boundary1 boundary2 boundary1x boundary1y boundary2x boundary2y minDistance indexofMin
        end
        
    end
    
end




%% save as a CSV
names = {'Trial','Spot','Animal','CellType','Depth','TrialName','ROI_X',...
    'ROI_X_type','ROI_Y','ROI_Y_type','LinCorr','LinPvalue',...
    'xCorr','Lag','MinDistance','numROIs'};

CorrelationData=All_Corrs;
CorrelationData(:,9) = [];
CorrelationData(:,11) = [];


AllData=[names;CorrelationData];

% write date to created file
cell2csv(FileSave{1,1}, AllData);
%save(FileSave{1,1}, 'All_Corrs','-v7.3');

