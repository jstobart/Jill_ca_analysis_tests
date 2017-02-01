% all ROIs  correlation
clearvars


% load data traces
load('E:\Data\Two_Photon_Data\GCaMP_RCaMP\cyto_GCaMP6s\Results\S&LStim_cGC&RC_traces_01_30_2017.mat');

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
        All_traces{iROI,12}='Neuron';
    elseif ~isempty(S_str)
        All_traces{iROI,12}='Soma';
    elseif ~isempty(EF_str)
        All_traces{iROI,12}='Endfeet';
    elseif ~isempty(P_str)
        All_traces{iROI,12}='Process';
    elseif ~isempty(NP_str)
        All_traces{iROI,12}='Neuropil';
    end
    
    % make new names
    All_traces{iROI,13}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,1});
    All_traces{iROI,14}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,1});
    All_traces{iROI,15}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,1},'_',All_traces{iROI,6});
    All_traces{iROI,16}=strcat(All_traces{iROI,5},'_',All_traces{iROI,4},'_',All_traces{iROI,2}, '_',All_traces{iROI,6});
    
    %if a process ROI is overlapping with a soma or process, exclude it
    %from data
    %All_traces{iROI,11}=0;
    nonOverlapIdx(iROI)=~ischar(All_traces{iROI,11});
    
end


data_traces = All_traces(str_idx',:);

% remove overlapping processes
All_traces = All_traces(nonOverlapIdx',:);

%% Correlation of each ROI with all ROIs in the same trial
% CONSIDER WHOLE TRIAL and 25 sec stim window

% only compare ROIs from the same trial together
trial_condition=unique(All_traces(:,16));

for itrial= 1:length(trial_condition)
    CurrentTrial=trial_condition{itrial};
    
    for iTrace= 1:size(All_traces,1)
        TrialIdx(iTrace)=strcmp(All_traces{iTrace,16}, CurrentTrial);
    end
    
    CurrentData=All_traces(TrialIdx,:);
    
    for iROI=1:size(CurrentData,1)
        CorrData{iROI,1:6}=CurrentData{iROI,2:7};
        CorrData{iROI,7}=CurrentData{iROI,16};
        
        CorrData{iROI,8}= CurrentData{iROI,1};
        CorrData{iROI,9}= CurrentData{iROI,12};
        
        %whole trace
        TraceX= CurrentData{iROI,8};
        % stim window
        TraceX_short=TraceX(1:(25*FrameRate),:);
        
        for kROI = (iROI+1):size(CurrentData,1)
            CorrData{iROI,10}= CurrentData{kROI,1};
            CorrData{iROI,11}= CurrentData{kROI,12};
            
            %whole trace
            TraceY= CurrentData{kROI,8};
            % stim window
            TraceY_short=TraceY(1:(25*FrameRate),:);
            
            % whole trial correlation
            [R,Ps] = corrcoef(TraceX, TraceY);
            CorrData{iROI,12}= R(1,2); % correlation coeff
            CorrData{iROI,13}= Ps(1,2); % p value
            
            % stim window correlation
            [R2,Ps2]= corrcoef(TraceX_short, TraceY_short);
            CorrData{iROI,14}= R2(1,2); % correlation coeff
            CorrData{iROI,15}= Ps2(1,2); % p value
            
            % stim window cross correlation
            [Corr_Sequence,lag]  = xcorr(TraceX_short, TraceY_short,'coeff');
            [~,I] = max(abs(Corr_Sequence));
            timeDiff=(lag(I))*FrameRate;
            CorrData{iROI,16}=  max(abs(Corr_Sequence)); % max of correlation sequence
            CorrData{iROI,17}= timeDiff; % lag time
        end
        
    end
    
    % CORRELATE SIGNALS
    
    % Preallocation
    Long_Rcor = cell(size(ROI_stuff,1),size(ROI_stuff,1),mtrials,mCond);
    Long_P = cell(size(ROI_stuff,1),size(ROI_stuff,1),mtrials,mCond);
    ROILabel = cell(mtrials,mCond);
    dist= cell(mtrials,size(ROI_stuff,1),size(ROI_stuff,1),mCond);
    Mean_Rcor = [];
    
    % Build a ROI matrix
    if ~exist('tempNames2','var')
        tempNames2=tempNames;
    end
    ROI_stuff(:,1) = tempNames2;
    for iROI = 1:length(tempNames2)
        ROI_indx{iROI} = ~cellfun(@isempty,regexp(ROI_info,tempNames2{iROI}));
        %x_info
    end
    for iROI = 1:length(tempNames2)
        xpos{iROI} = x_info(ROI_indx{iROI});
        ypos{iROI} = y_info(ROI_indx{iROI});
        ROI_stuff(iROI,2) = num2cell(xpos{iROI}(1,1));
        ROI_stuff(iROI,3) = num2cell(ypos{iROI}(1,1));
    end
    
    clear matches mat CompROIs matchROIs tempNames2 ROI_info
    
    
    % correlate each ROI with all other ROIs
    for iCond = 1:mCond
        for kTrial = 1:mtrials
            if ~exist('ROInum','var')
                for kROI = 1:numROIs
                    %iROI =
                    tempX = tempAll(:,kROI,kTrial,iCond); % load X trace data
                    % normalize by calculating change in fluorescence
                    BLX = mean(tempX(1:10));
                    dffX = (tempX-BLX)/BLX;
                    for iROI = 1:numROIs % load other somal data
                        tempY = tempAll(:,iROI,kTrial,iCond); % load Y trace data
                        % normalize by calculating change in fluorescence
                        BLY = mean(tempY(1:10));
                        dffY = (tempY-BLY)/BLY;
                        % correlation of the whole trial
                        [Long_R,Long_Ps] = corrcoef(dffX, dffY);
                        Long_Rcor{kROI,iROI,kTrial,iCond}= Long_R(1,2); % correlation coeff
                        Long_P{kROI,iROI,kTrial,iCond} = Long_Ps(1,2); % p value
                    end
                end
            elseif exist('ROInum','var') && length(ROInum)==1;
                for xx = 1:size(ROI_stuff,1)
                    iROI = [1:(ROInum(1)-1),(ROInum(1)+1):numROIs];
                    kROI = iROI(xx);
                    tempX = tempAll(:,kROI,kTrial,iCond); % load X trace data
                    % normalize by calculating change in fluorescence
                    BLX = mean(tempX(1:10));
                    dffX = (tempX-BLX)/BLX;
                    for jj = 1:size(ROI_stuff,1) % load other somal data
                        iROI = iROI(jj);
                        tempY = tempAll(:,iROI,kTrial,iCond); % load Y trace data
                        % normalize by calculating change in fluorescence
                        BLY = mean(tempY(1:10));
                        dffY = (tempY-BLY)/BLY;
                        % correlation of the whole trial
                        [Long_R,Long_Ps] = corrcoef(dffX, dffY);
                        Long_Rcor{xx,jj,kTrial,iCond}= Long_R(1,2); % correlation coeff
                        Long_P{xx,jj,kTrial,iCond} = Long_Ps(1,2); % p value
                    end
                end
            elseif exist('ROInum','var') && length(ROInum)==2;
                for xx = 1:size(ROI_stuff,1)
                    if ROInum(2)-ROInum(1)==0
                        iROI =  [1:(ROInum(1)-1),(ROInum(1)+2):numROIs];
                    else
                        iROI =  [1:(ROInum(1)-1),(ROInum(1)+1):(ROInum(2)-1),(ROInum(2)+1):numROIs];
                    end
                    kROI = iROI(xx);
                    tempX = tempAll(:,kROI,kTrial,iCond); % load X trace data
                    % normalize by calculating change in fluorescence
                    BLX = mean(tempX(1:10));
                    dffX = (tempX-BLX)/BLX;
                    for jj = 1:size(ROI_stuff,1) % load other somal data
                        iROI = iROI(jj);
                        tempY = tempAll(:,iROI,kTrial,iCond); % load Y trace data
                        % normalize by calculating change in fluorescence
                        BLY = mean(tempY(1:10));
                        dffY = (tempY-BLY)/BLY;
                        % correlation of the whole trial
                        [Long_R,Long_Ps] = corrcoef(dffX, dffY);
                        Long_Rcor{xx,jj,kTrial,iCond}= Long_R(1,2); % correlation coeff
                        Long_P{xx,jj,kTrial,iCond} = Long_Ps(1,2); % p value
                    end
                end
            elseif exist('ROInum','var') && length(ROInum)==3;
                for xx = 1:size(ROI_stuff,1)
                    if ROInum(2)-ROInum(1)==0
                        iROI =  [1:(ROInum(1)-1),(ROInum(1)+2):(ROInum(3)-1),(ROInum(3)+1):numROIs];
                    elseif ROInum(3)-ROInum(2)==0
                        iROI =  [1:(ROInum(1)-1),(ROInum(1)+1):(ROInum(2)-1),(ROInum(2)+2):numROIs];
                    else
                        iROI =  [1:(ROInum(1)-1),(ROInum(1)+1):(ROInum(2)-1),(ROInum(2)+1):(ROInum(3)-1),(ROInum(3)+1):numROIs];
                    end
                    kROI = iROI(xx);
                    tempX = tempAll(:,kROI,kTrial,iCond); % load X trace data
                    % normalize by calculating change in fluorescence
                    BLX = mean(tempX(1:10));
                    dffX = (tempX-BLX)/BLX;
                    for jj = 1:size(ROI_stuff,1) % load other somal data
                        iROI = iROI(jj);
                        tempY = tempAll(:,iROI,kTrial,iCond); % load Y trace data
                        % normalize by calculating change in fluorescence
                        BLY = mean(tempY(1:10));
                        dffY = (tempY-BLY)/BLY;
                        % correlation of the whole trial
                        [Long_R,Long_Ps] = corrcoef(dffX, dffY);
                        Long_Rcor{xx,jj,kTrial,iCond}= Long_R(1,2); % correlation coeff
                        Long_P{xx,jj,kTrial,iCond} = Long_Ps(1,2); % p value
                    end
                end
            else
                Error('Check ROI comparisons- 3 or more matching ROIs')
            end
        end
        Mean_Rcor(:,:,iCond) = mean(cell2mat(Long_Rcor(:,:,:,iCond)),3);
    end
    
    % distance calculations
    for iCond = 1:mCond
        for kTrial = 1:mtrials
            for kROI = 1:size(ROI_stuff,1)
                Loc1 = [ROI_stuff{kROI,2},ROI_stuff{kROI,3}];  % x and y ROI position
                for iROI = 1:size(ROI_stuff,1) % load other somal data
                    Loc2 = [ROI_stuff{iROI,2},ROI_stuff{iROI,3}];
                    
                    % distance calculations
                    Loc_df = Loc2-Loc1;
                    % convert from pixels to um, x = zoom
                    %(formula: y = 3.7347/x+0.0121)
                    if strcmp(CurrentAnimal,'GC62')
                        x = 8;
                    elseif strcmp(Spots{iSpot},'2013_5_15_spot1A')||strcmp(Spots{iSpot},'2013_6_11_spot1')||strcmp(Spots{iSpot},'2013_6_11_spot2')
                        x = 3.5;
                    elseif strcmp(Spots{iSpot},'2013_6_10_spot4')
                        x = 5;
                    elseif strcmp(Spots{iSpot},'2013_5_15_spot1B')
                        x = 6.5;
                    else
                        x = 4;
                    end
                    y = 3.7347/x +0.0121;
                    dist{kTrial,kROI,iROI,iCond} = norm(Loc_df*y);
                end
            end
        end
    end
    
    numROIs2 = size(Long_P,1);
    for  iCond =1:mCond;
        if iCond ==1;
            for kTrial = 1:mtrials
                for kROI = 1:numROIs2
                    for iROI = (kROI+1):numROIs2
                        %build giant data table
                        temp1{1,1} = CurrentAnimal;
                        temp1{1,2} = Spots{iSpot};
                        temp1{1,3} = Layer{1,1};
                        temp1{1,4} = 'Nostim';
                        temp1{1,5} = strcat('Trial',num2str(kTrial,'%02d'));
                        temp1{1,6} = cell2mat(Long_Rcor(kROI,iROI,kTrial,iCond));
                        temp1{1,7} = cell2mat(dist(kTrial,kROI,iROI,iCond));
                        temp1{1,8} = cell2mat(ROI_stuff(kROI,1));
                        temp1{1,9} = cell2mat(ROI_stuff(iROI,1));
                        temp1{1,10} = cell2mat(Long_P(kROI,iROI,kTrial,iCond));
                        %                                 if (~isempty(strfind(temp1{1,8},'E')) && ~isempty(strfind(temp1{1,9},'RO')))
                        %                                     EfvsP_NS = vertcat(EfvsP_NS,temp1(1,:));
                        %                                 elseif (~isempty(strfind(temp1{1,8},'S')) && ~isempty(strfind(temp1{1,9},'RO')));
                        %                                     SvsP_NS = vertcat(SvsP_NS,temp1(1,:));
                        %                                 elseif (~isempty(strfind(temp1{1,8},'E')) && ~isempty(strfind(temp1{1,9},'S')));
                        %                                     EfvsS_NS = vertcat(EfvsS_NS,temp1(1,:));
                        %                                 elseif (~isempty(strfind(temp1{1,8},'S')) && ~isempty(strfind(temp1{1,9},'S')));
                        %                                     SvsS_NS = vertcat(SvsS_NS,temp1(1,:));
                        %                                 elseif (~isempty(strfind(temp1{1,8},'RO')) && ~isempty(strfind(temp1{1,9},'RO')));
                        %                                     PvsP_NS = vertcat(PvsP_NS,temp1(1,:));
                        %                                 elseif(~isempty(strfind(temp1{1,8},'E')) && ~isempty(strfind(temp1{1,9},'E')));
                        %                                     EfvsEf_NS = vertcat(EfvsEf_NS,temp1(1,:));
                        if ~isempty(strfind(temp1{1,8},'E'))
                            Endfeet_NS = vertcat(Endfeet_NS,temp1(1,:));
                        elseif ~isempty(strfind(temp1{1,8},'S'))
                            Somas_NS = vertcat(Somas_NS,temp1(1,:));
                        else %~isempty(strfind(temp1{1,8},'RO'))
                            Processes_NS = vertcat(Processes_NS,temp1(1,:));
                        end
                    end
                end
            end
        else
            for kTrial = 1:mtrials
                for kROI = 1:numROIs2
                    for iROI = (kROI+1):numROIs2
                        %build giant data table
                        temp1{1,1} = CurrentAnimal;
                        temp1{1,2} = Spots{iSpot};
                        temp1{1,3} = Layer{1,1};
                        temp1{1,4} = 'Stim';
                        temp1{1,5} = strcat('Trial',num2str(kTrial,'%02d'));
                        temp1{1,6} = cell2mat(Long_Rcor(kROI,iROI,kTrial,iCond));
                        temp1{1,7} = cell2mat(dist(kTrial,kROI,iROI,iCond));
                        temp1{1,8} = cell2mat(ROI_stuff(kROI,1));
                        temp1{1,9} = cell2mat(ROI_stuff(iROI,1));
                        temp1{1,10} = cell2mat(Long_P(kROI,iROI,kTrial,iCond));
                        %                                 if (~isempty(strfind(temp1{1,8},'E')) && ~isempty(strfind(temp1{1,9},'RO')))
                        %                                     EfvsP_S = vertcat(EfvsP_S,temp1(1,:));
                        %                                 elseif (~isempty(strfind(temp1{1,8},'S')) && ~isempty(strfind(temp1{1,9},'RO')));
                        %                                     SvsP_S = vertcat(SvsP_S,temp1(1,:));
                        %                                 elseif (~isempty(strfind(temp1{1,8},'E')) && ~isempty(strfind(temp1{1,9},'S')));
                        %                                     EfvsS_S = vertcat(EfvsS_S,temp1(1,:));
                        %                                 elseif (~isempty(strfind(temp1{1,8},'S')) && ~isempty(strfind(temp1{1,9},'S')));
                        %                                     SvsS_S = vertcat(SvsS_S,temp1(1,:));
                        %                                 elseif (~isempty(strfind(temp1{1,8},'RO')) && ~isempty(strfind(temp1{1,9},'RO')));
                        %                                     PvsP_S = vertcat(PvsP_S,temp1(1,:));
                        %                                 elseif(~isempty(strfind(temp1{1,8},'E')) && ~isempty(strfind(temp1{1,9},'E')));
                        %                                     EfvsEf_S = vertcat(EfvsEf_S,temp1(1,:));
                        %else
                        % if none of the conditions is true '
                        %   fprintf('None of the values are matching\n');
                        if ~isempty(strfind(temp1{1,8},'E'))
                            Endfeet_S = vertcat(Endfeet_S,temp1(1,:));
                        elseif ~isempty(strfind(temp1{1,8},'S'))
                            Somas_S = vertcat(Somas_S,temp1(1,:));
                        else %~isempty(strfind(temp1{1,8},'RO'))
                            Processes_S = vertcat(Processes_S,temp1(1,:));
                        end
                    end
                end
            end
        end
    end
    
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
    
    clear tempSig1 tempSig2 tempAll Labelnames Spot_info ROI_info x_info y_info xpos ypos ROI_indx
    clear Mean_Rcor temp1 ROI_stuff ROInum SpotROIs numROIs numROIs2
    
end
end
clear Spots

end

%% save as a CSV
names = {'Animal','Spot', 'Layer','Condition','Trial','Corr','Distance','ROI_X','ROI_Y','Pvalue'};

Endfeet_S = vertcat(names,Endfeet_S);
Endfeet_NS = vertcat(names,Endfeet_NS);
Somas_S = vertcat(names,Somas_S);
Somas_NS = vertcat(names,Somas_NS);
Processes_S = vertcat(names,Processes_S);
Processes_NS = vertcat(names,Processes_NS);


cd(fullfile(MainDir, 'Results'));
% write date to created file
cell2csv('Nostim_Endfeet_Corr.csv', Endfeet_NS);
cell2csv('Nostim_Somas_Corr.csv', Somas_NS);
cell2csv('Nostim_Processes_Corr.csv', Processes_NS);
cell2csv('Stim_Endfeet_Corr.csv', Endfeet_S);
cell2csv('Stim_Somas_Corr.csv', Somas_S);
cell2csv('Stim_Processes_Corr.csv', Processes_S);
%%  Plots (scatter and histograms)
%         figure
%         subplot(1,2,1)
%         hold on
%         title('Endfeet vs Processes')
%         scatter(cell2mat(EfvsP_NS(:,8)),cell2mat(EfvsP_NS(:,7)),'b')
%         scatter(cell2mat(EfvsP(:,8)),cell2mat(EfvsP(:,7)),'r')
%         hold off
%         subplot(1,2,2)
%         hold on
%         histogram(cell2mat(EfvsP_NS(:,5)),'BinWidth',0.025,'FaceColor','b')
%         histogram(cell2mat(EfvsP(:,5)),'BinWidth',0.025,'FaceColor','r')
%         hold off
%
%         figure
%         subplot(1,2,1)
%         hold on
%         title('Somas vs Processes')
%         scatter(cell2mat(SvsP_NS(:,8)),cell2mat(SvsP_NS(:,7)))
%         scatter(cell2mat(SvsP(:,8)),cell2mat(SvsP(:,7)))
%         hold off
%         subplot(1,2,2)
%         hold on
%         title('Somas vs Processes')
%         histogram(cell2mat(SvsP_NS(:,5)),'BinWidth',0.025,'FaceColor','b')
%         histogram(cell2mat(SvsP(:,5)),'BinWidth',0.025,'FaceColor','r')
%         hold off
%
%         figure
%         subplot(1,2,1)
%         hold on
%         title('Endfeet vs Somas')
%         scatter(cell2mat(EfvsS_NS(:,8)),cell2mat(EfvsS_NS(:,7)))
%         scatter(cell2mat(EfvsS(:,8)),cell2mat(EfvsS(:,7)))
%         hold off
%         subplot(1,2,2)
%         hold on
%         title('Endfeet vs Somas')
%         histogram(cell2mat(EfvsS_NS(:,5)),'BinWidth',0.025,'FaceColor','b')
%         histogram(cell2mat(EfvsS(:,5)),'BinWidth',0.025,'FaceColor','r')
%         hold off
%
% %%
% % save('crossCor_11_6_15.mat','Rxcor','-v7.3')
% % save('crossCor_lag_11_6_15.mat','Rxcor_lags','-v7.3')
% % save('RCor_11_6_15.mat','Rcor','-v7.3')
% % save('RCor_mean_11_6_15.mat','Rcor_mean','-v7.3')
% % save('RCor_pvalues_11_6_15.mat','Pvalue','-v7.3')
