ExperimentFolder = 'E:\Data\Two_Photon_Data\P2Y1_Mice\';

%ZipFolder = 'P:\_Group\Projects\Astrocyte Calcium\Current Milestones\Astrocyte Calcium Imaging During Plasticity\Imaging Data\ImagesforClicking\Controls';

AnimalNames = {
    'PY05',...    %control
    };

ScoreSheetNames = {
    'PY05_Scoresheet_Imaging_WT.xls',...
    };

spot1Path = {
    };

spot2Path = {
    'E:\Data\Two_Photon_Data\P2Y1_Mice\PY05\2015_09_28\Lck_spot2_WT_nostim\highres_spot2_stim020.tif',...
    };


ScoreSheetFolder = ...
    'E:\Data\Two_Photon_Data\P2Y1_Mice\Scoresheets';

ScoreSheetFolders = fullfile(ScoreSheetFolder, ScoreSheetNames);

useParallel = false;
useHandROIs = false;

for iAnimal = 1:numel(AnimalNames)
    
    savepath = fullfile(ExperimentFolder, AnimalNames{iAnimal}, ...
        '151006_ResultsTable_PY05.mat');
    
    if exist(savepath, 'file')
        continue
    end
    
    [~,spot1] = scim.scim_openTif([ExperimentFolder, ...
        spot1Path{iAnimal}]);
    [~,spot2] = scim.scim_openTif([ExperimentFolder, ...
        spot2Path{iAnimal}]);

    spot1Ref = mean(squeeze(spot1(:,:,1,:)),3);
    spot2Ref = mean(squeeze(spot2(:,:,1,:)),3);
    
    refImg = cat(3, spot1Ref, spot2Ref);
    
    AnimalObj = Animal(ExperimentFolder, ...
                'WT', ...
                AnimalNames(iAnimal), ...
                ScoreSheetFolders(iAnimal));
    
    
    currROIFolder = fullfile(ZipFolder, AnimalNames{iAnimal});
    
    if ~strcmp(AnimalObj.state, 'preprocessed')
        AnimalObj.preprocess(refImg, useHandROIs, currROIFolder);
    end
    
    nBaseSessions = 1;
    
    AnimalObj.process(nBaseSessions, useParallel, useHandROIs);
    
    % Save data, because of RAM problem
    dataTableTemp = AnimalObj.output_data();
    save(savepath, 'dataTableTemp', '-v7.3');
    
    clearvars -except AnimalNames ExperimentFolder ScoreSheetNames spot1Path ...
        spot2Path ScoreSheetFolders useParallel useHandROIs ZipFolder
end

newTable = table();
dataTable = table();
for iAnimal = 1:numel(AnimalNames)
    
    savepath = fullfile(ExperimentFolder, AnimalNames{iAnimal}, ...
        '150902_ResultsTable_handclicked.mat');
    load(savepath)
    
    Animal = repmat(AnimalNames(iAnimal), size(dataTableTemp, 1), 1);
    animalCol = table(Animal);
    dataTable = [animalCol, dataTableTemp];
    newTable = [newTable; dataTable];
    
end

writetable(newTable, fullfile(ExperimentFolder, 'Results', ...
    'Table_PY05.csv'), 'Delimiter', ',');