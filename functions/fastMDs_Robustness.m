%% Robustness of fast MDs

load('D:\Data\GCaMP_RCaMP\Revision\Lck_GCaMP\FilesforMatlab\Control_untreated\Traces_1stCohort_Lck_nostim_vs_longstim_12_2017.mat')
All_traces(:,14)=strcat(All_traces(:,5),{'_'},All_traces(:,4),{'_'},All_traces(:,2),...
    {'_'},All_traces(:,1),{'_'},All_traces(:,6));

Cohort1=All_traces;
Cohort1(:,7)= [];
Cohort1(:,7)= [];

load('D:\Data\GCaMP_RCaMP\Revision\Lck_GCaMP\FilesforMatlab\Control_untreated\Traces_2ndCohort_Lck_nostim_vs_longstim_01_2018.mat')
All_traces(:,18)=strcat(All_traces(:,5),{'_'},All_traces(:,4),{'_'},All_traces(:,2),...
    {'_'},All_traces(:,1),{'_'},All_traces(:,6));

All_traces(:,7)= [];
All_traces(:,7)= [];
All_traces(:,7)= [];
All_traces(:,12)= [];
All_traces(:,12)= [];
All_traces(:,12)= [];

AllData=vertcat(Cohort1, All_traces);



respondingData=load('D:\Data\GCaMP_RCaMP\Revision\Lck_GCaMP\FilesforMatlab\Control_untreated\Astrocyte_respondingROIs.mat');
respondingROIs=respondingData.x.ROIs_Cond;
respondingROIs(:,2)=respondingData.x.Group;

DelayedMatch=ismember(respondingROIs(:,2), {'delayed'});
DelayedROIs=respondingROIs(DelayedMatch,1);
FastROIs=respondingROIs(~DelayedMatch,1);


Match1=ismember(AllData(:,12), DelayedROIs);
DelayedTraces=AllData(Match1,:);

Match2=ismember(AllData(:,12), FastROIs);
FastTraces=AllData(Match2,:);


clearvars AllData All_traces Cohort1

%% FAST ROIs
bigtable = cell2table(FastTraces);
bigtable = bigtable(:, [1:6, 8,12]);

% 127x128 px images
nRows = 127;
nCols = 128;

% Only FLIKA ROIs
hcIdx = regexp(bigtable{:,1}, regexptranslate('wildcard', 'roi*'));
hcIdx = cellfun(@isempty, hcIdx);
bigtable = bigtable(~hcIdx,:);

% Response probability in n trials, divided by number of active px
stimFracAc = [];
nostimFracAc = [];

% Total number of individual FOV images
masks = unique([bigtable(:,2), bigtable(:,4), bigtable(:,5), bigtable(:,6)]);
masks = [masks, table(cell(size(masks,1),1), 'VariableName', {'trialmask'}), table(cell(size(masks,1),1), ...
    'VariableName', {'score'}), table(cell(size(masks,1),1), 'VariableName', {'FOVmask'})];

% Create ROI image for single trial
animals = unique(bigtable{:,5});
for iAnimal = 1:numel(animals)
    animalIdx = strcmp(bigtable{:,5}, animals{iAnimal});
    animalTable = bigtable(animalIdx, :);
    
    conditions = unique(animalTable{:,6});
    for iCond = 1:numel(conditions)
        condIdx = strcmp(animalTable{:,6}, conditions{iCond});
        condTable = animalTable(condIdx, :);
        
        spots = unique(condTable{:,4});
        for iSpot = 1: numel(spots)
            spotIdx = strcmp(condTable{:,4}, spots{iSpot});
            spotTable = condTable(spotIdx, :);
            
            trials = unique(spotTable{:,2});
            trialmasks = [];
            for iTrial = 1:numel(trials)
                trialIdx = strcmp(spotTable{:,2}, trials{iTrial});
                trialTable = spotTable(trialIdx, :);
                
                rois = unique(trialTable{:,1});
                mask = zeros(nRows, nCols);
                
                for iROI = 1:numel(rois)
                    roiIdx = strcmp(trialTable{:,1}, rois{iROI});
                    roiPx = trialTable{roiIdx, 7}{1};
                    
                    mask(roiPx) = true;
                end
                
                % Store mask
                for iLen = 1:size(masks, 1)
                    if all(strcmp([masks{iLen,1:4}], [trials(iTrial), spots(iSpot), animals(iAnimal), conditions(iCond)]))
                        masks(iLen, 5) = cell2table({mask});
                        break
                    end
                end
                
                trialmasks(:,:,iTrial) = mask;
            end
            
            
            fracImg = sum(trialmasks,3)./numel(trials);
            
            for iLen = 1:size(masks, 1)
                if ~isempty(masks{iLen, 5}) && all(strcmp([masks{iLen,2:4}], [spots(iSpot), animals(iAnimal), conditions(iCond)]))
                    masks{iLen, 7} = {fracImg};
                    break
                end
            end
            
            %% method with all active pixels 
%             activePx = sum(sum(sumImg>0));
%             score = sum(sum(fracImg)) / activePx;
            
            %% alternative method with threshold
            thresh = 1/numel(trials); % or 1/numel(trials)
            activePxIdx = fracImg > 0;
            activePx = sum(activePxIdx(:));
            fracImg(fracImg <= thresh) = NaN;
            score = nansum(nansum(fracImg)) / activePx;
            
            if strcmp(conditions{iCond}, 'Stim')
                stimFracAc = [stimFracAc, activePx/(127*128)];
            elseif strcmp(conditions{iCond}, 'Nostim')
                nostimFracAc = [nostimFracAc, activePx/(127*128)];
            else
                warning('Uncaught condition');
            end
            
%             % Store score
            for iLen = 1:size(masks, 1)
                if ~isempty(masks{iLen, 5}) && all(strcmp([masks{iLen,2:4}], [spots(iSpot), animals(iAnimal), conditions(iCond)]))
                    masks{iLen, 6} = {score};
                    break
                end
            end
            
        end        
    end
end

% Remove duplicates
badIdx = cellfun(@isempty, masks{:,5});
masks = masks(~badIdx, :);

%% Create main table "scores"
scores = masks(:,2:end);
scores = scores(~cellfun(@isempty, scores{:,end}), :);

% Extract values for some simple tests
allscores = scores{:,5};
allscores=cell2mat(allscores);
%nostim = allscores(1:2:size(allscores,1));
%stim = allscores(2:2:size(allscores,1));

% Simple ttest on stim/nostim scores
%[h,p]=ttest(nostim, stim, 'Alpha', 0.001);

%% Plot random sample

% Show figure
figure;
imagesc(scores{36,6}{1});
axis square
%title(sprintf('fast; score = %f', scores{3,1}{1}))
caxis([0 1])
colormap('parula')
colorbar



