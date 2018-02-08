%% Robustness of FOV

load('P:\Kim\Traces_2ndCohort_Lck_nostim_vs_longstim_01_2018.mat')
bigtable = cell2table(All_traces);
bigtable = bigtable(:, [1:6, 11]);

% 127x128 px images
nRows = 127;
nCols = 128;

% Only GCaMP ROIs
gIdx = strcmp(bigtable{:,3}, 'GCaMP');
bigtable = bigtable(gIdx, :);

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
                        masks(iLen, 5) = {mask};
                        break
                    end
                end
                
                trialmasks(:,:,iTrial) = mask;
            end
            
            
            fracImg = sum(trialmasks,3)./numel(trials);
            
            for iLen = 1:size(masks, 1)
                if ~isempty(masks{iLen, 5}) && all(strcmp([masks{iLen,1:4}], [{'trial01'}, spots(iSpot), animals(iAnimal), conditions(iCond)]))
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
            
            % Store score
            for iLen = 1:size(masks, 1)
                if ~isempty(masks{iLen, 5}) && all(strcmp([masks{iLen,1:4}], [{'trial01'}, spots(iSpot), animals(iAnimal), conditions(iCond)]))
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
nostim = allscores(1:2:size(allscores,1));
stim = allscores(2:2:size(allscores,1));

% Simple ttest on stim/nostim scores
[h,p]=ttest(nostim, stim, 'Alpha', 0.001)

%% Plot random sample

% remove entries that don't have a matching stim or no stim
% cytosolic data
%scores([5,72,75],:) = [];

% Lck data  1st cohort
%scores([1,32,33,44,53,78],:) = [];


% sampIdx = randi(size(scores,1));
% if mod(sampIdx,2)
%     sampIdx = [sampIdx, sampIdx+1];
% else
%     sampIdx = [sampIdx-1, sampIdx];
% end
samp = scores(sampIdx, :);

% Show figure
figure;
subplot(1,2,1)
imagesc(samp{1,5}{1});
axis square
title(sprintf('Nostim; score = %f', samp{1,4}{1}))
caxis([0 1])
colormap('parula')
colorbar

subplot(1,2,2)
imagesc(samp{2,5}{1});
axis square
title(sprintf('Stim; score = %f', samp{2,4}{1}))
caxis([0 1])
colormap('jet')
colorbar

% cyto example:  RG14-18_02_2016.... not great example
% lck example: RG14 08_03_2016.... ok

%RG14- 16_03_10_spot1
% RG17 -{'16_02_24_spot1'}
% RG14 -{'16_02_24_spot1'}


% alternative Lck examples
% RG14 {'16_02_26_spot2'}
% IPRG2 {'17_11_22_spot5'}
% RG14 {'16_03_04_spot1'}
% RG16 {'16_03_11_spot1'}


%WT_LR4.... %spot3

%% pharmacology traces
load('D:\Data\GCaMP_RCaMP\Revision\Lck_GCaMP\FilesforMatlab\Robustness\Results_pharmacology_Lck_nostim_vs_longstim_12_2017.mat');
pharmacology=scores;
DMSOIdx = regexp(pharmacology{:,4}, regexptranslate('wildcard', 'DMSO*'));
DMSOIdx = cellfun(@isempty, DMSOIdx);
pharmacology = pharmacology(DMSOIdx, :);

gIdx2 = strcmp(pharmacology{:,2}, 'WT_LR1');
pharma_WT_LR1 = pharmacology(gIdx2, :);

%load('D:\Data\GCaMP_RCaMP\Revision\Lck_GCaMP\FilesforMatlab\Robustness\Results_1stCohort_Lck_nostim_vs_longstim_12_2017.mat');
load('D:\Data\GCaMP_RCaMP\Revision\Lck_GCaMP\FilesforMatlab\Robustness\Results_2ndCohort_Lck_nostim_vs_longstim_01_2018.mat');
control=scores;
gIdx = strcmp(control{:,2}, 'WT_LR1');
control_WT_LR1 = control(gIdx, :);


%%
% spot 1
% Show figure
figure;
subplot(2,3,1)
imagesc(control_WT_LR1{10,5}{1});
axis square
title(sprintf('control; score = %f', control_WT_LR1{2,4}{1}))
caxis([0 1])
colormap('parula')
colorbar

subplot(2,3,2)
imagesc(pharma_WT_LR1{29,7}{1});
axis square
title(sprintf('atropine; score = %f', pharma_WT_LR1{13,6}{1}))
caxis([0 1])
colormap('parula')
colorbar

subplot(2,3,3)
imagesc(pharma_WT_LR1{30,7}{1});
axis square
title(sprintf('metergoline; score = %f', pharma_WT_LR1{6,6}{1}))
caxis([0 1])
colormap('jet')
colorbar

subplot(2,3,4)
imagesc(pharma_WT_LR1{31,7}{1});
axis square
title(sprintf('prazosin; score = %f', pharma_WT_LR1{7,6}{1}))
caxis([0 1])
colormap('jet')
colorbar

subplot(2,3,5)
imagesc(pharma_WT_LR1{32,7}{1});
axis square
title(sprintf('trazodone; score = %f', pharma_WT_LR1{8,6}{1}))
caxis([0 1])
colormap('jet')
colorbar