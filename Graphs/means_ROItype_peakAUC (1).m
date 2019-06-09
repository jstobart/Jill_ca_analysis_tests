%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds the largest peak in every trial and calculates the average area %
% under curve over all trials in one ROI. Then calculates average area  %
% under the curve for peaks in each ROI type.                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1st Version 29.08.13 KF


clear all
%open data files
DirName = '/Users/kimferrari1/Documents/Uni/FS13/Forschungspraktikum/Plasticity/Analysis_Results';
cd(DirName)
flist=dir(fullfile(DirName, '*.mat')); % will give you list of specifed files
for ii=1:length(flist)
     load(flist(ii).name) % load files in workspace
end

%% day0-25 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for x = [0,4,7,11,18,25]
    t = eval(['Spared_Day' num2str(x)]);
    name = ['Spared_Lay1_Day' num2str(x)];
    
%% Layer 1 Trimmed

RowInfo = ['AP_stim  ';'AP_nostim';'AC_stim  ';'AC_nostim';'BR_stim  ';'BR_nostim';'EF_stim  ';'EF_nostim'];
MEAN.(name) = cellstr(RowInfo);

%% Average area under curve for all AP ROIs in stimulation trials

y = [];
%y = cells(10,length(t.AP_Lay1_8secstim)) %preallocation for y
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.AP_Lay1_8secstim)  % No. of ROIs
    for itrial = 1:length(t.AP_Lay1_8secstim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.AP_Lay1_8secstim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){1,2} = mean(nanmean(cellfun(@mean, m)));                 %mean over all trials

%MEAN.(name){1,2} = mean(nanmean(cellfun(@mean, m)));;
   


%% Average area under curve for all AP ROIs in no-stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.AP_Lay1_Nostim)  % No. of ROIs
    for itrial = 1:length(t.AP_Lay1_Nostim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.AP_Lay1_Nostim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){2,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all AC ROIs in stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.AC_Lay1_8secstim)  % No. of ROIs
    for itrial = 1:length(t.AC_Lay1_8secstim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.AC_Lay1_8secstim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){3,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all AC ROIs in no-stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.AC_Lay1_Nostim)  % No. of ROIs
    for itrial = 1:length(t.AC_Lay1_Nostim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.AC_Lay1_Nostim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){4,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all BR ROIs in stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.BR_Lay1_8secstim)  % No. of ROIs
    for itrial = 1:length(t.BR_Lay1_8secstim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.BR_Lay1_8secstim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){5,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all BR ROIs in no-stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.BR_Lay1_Nostim)  % No. of ROIs
    for itrial = 1:length(t.BR_Lay1_Nostim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.BR_Lay1_Nostim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){6,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all EF ROIs in stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.EF_Lay1_8secstim)  % No. of ROIs
    for itrial = 1:length(t.EF_Lay1_8secstim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.EF_Lay1_8secstim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){7,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all EF ROIs in no-stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.EF_Lay1_Nostim)  % No. of ROIs
    for itrial = 1:length(t.EF_Lay1_Nostim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.EF_Lay1_Nostim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){8,2} = mean(nanmean(cellfun(@mean, m)));
end
for x = [0,4,7,11,18,25]
    t = eval(['Trimmed_Day' num2str(x)]);
    name = ['Trimmed_Lay1_Day' num2str(x)];
    
%% Layer 1 Spared

RowInfo = ['AP_stim  ';'AP_nostim';'AC_stim  ';'AC_nostim';'BR_stim  ';'BR_nostim';'EF_stim  ';'EF_nostim'];
MEAN.(name) = cellstr(RowInfo);

%% Average area under curve for all AP ROIs in stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.AP_Lay1_8secstim)  % No. of ROIs
    for itrial = 1:length(t.AP_Lay1_8secstim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.AP_Lay1_8secstim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){1,2} = mean(nanmean(cellfun(@mean, m)));
   


%% Average area under curve for all AP ROIs in no-stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.AP_Lay1_Nostim)  % No. of ROIs
    for itrial = 1:length(t.AP_Lay1_Nostim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.AP_Lay1_Nostim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){2,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all AC ROIs in stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.AC_Lay1_8secstim)  % No. of ROIs
    for itrial = 1:length(t.AC_Lay1_8secstim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.AC_Lay1_8secstim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){3,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all AC ROIs in no-stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.AC_Lay1_Nostim)  % No. of ROIs
    for itrial = 1:length(t.AC_Lay1_Nostim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.AC_Lay1_Nostim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){4,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all BR ROIs in stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.BR_Lay1_8secstim)  % No. of ROIs
    for itrial = 1:length(t.BR_Lay1_8secstim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.BR_Lay1_8secstim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){5,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all BR ROIs in no-stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.BR_Lay1_Nostim)  % No. of ROIs
    for itrial = 1:length(t.BR_Lay1_Nostim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.BR_Lay1_Nostim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){6,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all EF ROIs in stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.EF_Lay1_8secstim)  % No. of ROIs
    for itrial = 1:length(t.EF_Lay1_8secstim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.EF_Lay1_8secstim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){7,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all EF ROIs in no-stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.EF_Lay1_Nostim)  % No. of ROIs
    for itrial = 1:length(t.EF_Lay1_Nostim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.EF_Lay1_Nostim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){8,2} = mean(nanmean(cellfun(@mean, m)));
end
for x = [0,4,7,11,18,25]
    t = eval(['Spared_Day' num2str(x)]);
    name = ['Spared_Lay2_Day' num2str(x)];
    
%% Layer 1 Trimmed

RowInfo = ['AP_stim  ';'AP_nostim';'AC_stim  ';'AC_nostim';'BR_stim  ';'BR_nostim';'EF_stim  ';'EF_nostim'];
MEAN.(name) = cellstr(RowInfo);

%% Average area under curve for all AP ROIs in stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.AP_Lay2_8secstim)  % No. of ROIs
    for itrial = 1:length(t.AP_Lay2_8secstim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.AP_Lay2_8secstim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){1,2} = mean(nanmean(cellfun(@mean, m)));
   


%% Average area under curve for all AP ROIs in no-stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.AP_Lay2_Nostim)  % No. of ROIs
    for itrial = 1:length(t.AP_Lay2_Nostim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.AP_Lay2_Nostim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){2,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all AC ROIs in stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.AC_Lay2_8secstim)  % No. of ROIs
    for itrial = 1:length(t.AC_Lay2_8secstim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.AC_Lay2_8secstim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){3,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all AC ROIs in no-stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.AC_Lay2_Nostim)  % No. of ROIs
    for itrial = 1:length(t.AC_Lay2_Nostim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.AC_Lay2_Nostim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){4,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all BR ROIs in stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.BR_Lay2_8secstim)  % No. of ROIs
    for itrial = 1:length(t.BR_Lay2_8secstim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.BR_Lay2_8secstim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){5,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all BR ROIs in no-stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.BR_Lay2_Nostim)  % No. of ROIs
    for itrial = 1:length(t.BR_Lay2_Nostim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.BR_Lay2_Nostim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){6,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all EF ROIs in stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.EF_Lay2_8secstim)  % No. of ROIs
    for itrial = 1:length(t.EF_Lay2_8secstim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.EF_Lay2_8secstim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){7,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all EF ROIs in no-stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.EF_Lay2_Nostim)  % No. of ROIs
    for itrial = 1:length(t.EF_Lay2_Nostim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.EF_Lay2_Nostim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){8,2} = mean(nanmean(cellfun(@mean, m)));
end
for x = [0,4,7,11,18,25]
    t = eval(['Trimmed_Day' num2str(x)]);
    name = ['Trimmed_Lay2_Day' num2str(x)];
    
%% Layer 1 Spared

RowInfo = ['AP_stim  ';'AP_nostim';'AC_stim  ';'AC_nostim';'BR_stim  ';'BR_nostim';'EF_stim  ';'EF_nostim'];
MEAN.(name) = cellstr(RowInfo);

%% Average area under curve for all AP ROIs in stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.AP_Lay2_8secstim)  % No. of ROIs
    for itrial = 1:length(t.AP_Lay2_8secstim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.AP_Lay2_8secstim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){1,2} = mean(nanmean(cellfun(@mean, m)));
   


%% Average area under curve for all AP ROIs in no-stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.AP_Lay2_Nostim)  % No. of ROIs
    for itrial = 1:length(t.AP_Lay2_Nostim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.AP_Lay2_Nostim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){2,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all AC ROIs in stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.AC_Lay2_8secstim)  % No. of ROIs
    for itrial = 1:length(t.AC_Lay2_8secstim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.AC_Lay2_8secstim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){3,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all AC ROIs in no-stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.AC_Lay2_Nostim)  % No. of ROIs
    for itrial = 1:length(t.AC_Lay2_Nostim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.AC_Lay2_Nostim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){4,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all BR ROIs in stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.BR_Lay2_8secstim)  % No. of ROIs
    for itrial = 1:length(t.BR_Lay2_8secstim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.BR_Lay2_8secstim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){5,2} = mean(nanmean(cellfun(@mean, m)));

%% Average area under curve for all BR ROIs in no-stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.BR_Lay2_Nostim)  % No. of ROIs
    for itrial = 1:length(t.BR_Lay2_Nostim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.BR_Lay2_Nostim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){6,2} = mean(nanmean(cellfun(@mean, m)));;

%% Average area under curve for all EF ROIs in stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.EF_Lay2_8secstim)  % No. of ROIs
    for itrial = 1:length(t.EF_Lay2_8secstim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.EF_Lay2_8secstim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){7,2} = mean(nanmean(cellfun(@mean, m)));;

%% Average area under curve for all EF ROIs in no-stimulation trials

y = [];
m = [];

% find largest peaks per trial 
for iROI = 1:length(t.EF_Lay2_Nostim)  % No. of ROIs
    for itrial = 1:length(t.EF_Lay2_Nostim{1,iROI}{14,2})         % No. of trials per ROI
            y{itrial,iROI} = t.EF_Lay2_Nostim{1,iROI}{14,2}{1,itrial};
            m{itrial,iROI}= max(cell2mat(y{itrial,iROI}));
    end
end

MEAN.(name){8,2} = mean(nanmean(cellfun(@mean, m)));;
end

%% Create Output folder
file = '/Users/kimferrari1/Documents/Uni/FS13/Forschungspraktikum/Plasticity/Analysis_Results';
outputFolder = fullfile(file, 'Means');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

cd(outputFolder);

save('Means_peakAUC.mat', 'MEAN');
