library("lme4")
library("lmerTest")
#library("lattice")
library("plyr")
library("dplyr")
library("ggplot2")
#library("ggpubr")
library("lsmeans")
library("Rmisc")
#library("MASS")
library("multcomp")
library("reshape2")
library("tidyr")
#library("data.table")
library("Hmisc")
library("stringr")

########################

# theme for plots
max.theme <- theme_classic() + 
  theme(
    text=element_text(size=12),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14, face="bold"),
    axis.title.y=element_text(vjust=1),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14, face="bold"),
    axis.line.x = element_line(colour = 'black', size=1, linetype='solid'),
    axis.line.y = element_line(colour = 'black', size=1, linetype='solid')) 

###########
# NOTES
# COLOUR BLIND FRIENDLY PALETTE FOR PLOTS
# The palette with black:
cbbPalette <- c("#000000","#D55E00","#009E73","#E69F00","#56B4E9","#CC79A7","#F0E442")

########################
# load data

lck.peaks1 <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/Peaks_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
lck.OT1<-read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/OnsetTimes_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
lck.peaks2 <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/Peaks_2ndCohort_Lck_nostim_vs_longstim_01_2018.csv", header=TRUE, sep = ",")
lck.OT2<-read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/OnsetTimes_2ndCohort_Lck_nostim_vs_longstim_01_2018.csv", header=TRUE, sep = ",")


all.cyto.peaks <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/FilesforR/Peaks_allMice_cyto_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
all.cyto.OT<-read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/FilesforR/OnsetTimes_allMice_cyto_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")

# 2s before and after stim
longstim.lck.OT.2s <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/First_Neuron_Submission/Lck_longstim_onset_2sbeforeStim.csv", header=TRUE, sep = ",")
nostim.lck.OT.2s <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/First_Neuron_Submission/Lck_nostim_onset_2sbeforeStim.csv", header=TRUE, sep = ",")
longstim.cyto.OT.2s <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/First_Neuron_Submission/cyto_longstim_onset_2sbeforeStim.csv", header=TRUE, sep = ",")
nostim.cyto.OT.2s <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/First_Neuron_Submission/cyto_nostim_onset_2sbeforeStim.csv", header=TRUE, sep = ",")


######
#OLD DATA
#long.lck <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_longstim_28_04_2017.csv", header=TRUE, sep = ",")
#short.lck <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_shortstim_28_04_2017.csv", header=TRUE, sep = ",")
#nostim.lck <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_nostim_28_04_2017.csv", header=TRUE, sep = ",")

long.cyto<- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cytoGC&RC_2D_longstim_28_04_2017.csv", header=TRUE, sep = ",")
short.cyto <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cytoGC&RC_2D_shortstim_28_04_2017.csv", header=TRUE, sep = ",")
nostim.cyto <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cytoGC&RC_2D_nostim_28_04_2017.csv", header=TRUE, sep = ",")


# onset time comparisons for nostim data
#longstim.OT.comp <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/longstim_firstonset_comparisons.csv", header=TRUE, sep = ",")
#shortstim.OT.comp <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/shortstim_firstonset_comparisons.csv", header=TRUE, sep = ",")
#nostim.OT.comp <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/nostim_firstonset_comparisons.csv", header=TRUE, sep = ",")

#longstim.lck.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/Lck_longstim_onset&AUC.csv", header=TRUE, sep = ",")
#shortstim.lck.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/shortstim_onset&AUC.csv", header=TRUE, sep = ",")
#nostim.lck.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/Lck_nostim_onset&AUC.csv", header=TRUE, sep = ",")

longstim.lck.OT.2s <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/First_Neuron_Submission/Lck_longstim_onset_2sbeforeStim.csv", header=TRUE, sep = ",")
nostim.lck.OT.2s <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/First_Neuron_Submission/Lck_nostim_onset_2sbeforeStim.csv", header=TRUE, sep = ",")
longstim.cyto.OT.2s <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_longstim_onset_2sbeforeStim.csv", header=TRUE, sep = ",")
nostim.cyto.OT.2s <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_nostim_onset_2sbeforeStim.csv", header=TRUE, sep = ",")

longstim.cyto.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_longstim_onset&AUC.csv", header=TRUE, sep = ",")
shortstim.cyto.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_shortstim_onset&AUC.csv", header=TRUE, sep = ",")
nostim.cyto.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_nostim_firstonset&AUC.csv", header=TRUE, sep = ",")

##### 
#home files
#long.lck <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_longstim_28_04_2017.csv", header=TRUE, sep = ",")
#short.lck <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_shortstim_28_04_2017.csv", header=TRUE, sep = ",")
#nostim.lck <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_nostim_28_04_2017.csv", header=TRUE, sep = ",")

lck.peaks1 <- read.table("D:/Data/GCaMP_RCaMP/Revision/Lck_GCaMP/FilesforR/Control_untreated/Peaks_1stCohort_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
lck.OT1<-read.table("D:/Data/GCaMP_RCaMP/Revision/Lck_GCaMP/FilesforR/Control_untreated/OnsetTimes_1stCohort_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
lck.peaks2 <- read.table("D:/Data/GCaMP_RCaMP/Revision/Lck_GCaMP/FilesforR/Control_untreated/Peaks_2ndCohort_Lck_nostim_vs_longstim_01_2018.csv", header=TRUE, sep = ",")
lck.OT2<-read.table("D:/Data/GCaMP_RCaMP/Revision/Lck_GCaMP/FilesforR/Control_untreated/OnsetTimes_2ndCohort_Lck_nostim_vs_longstim_01_2018.csv", header=TRUE, sep = ",")


all.cyto.peaks <- read.table("D:/Data/GCaMP_RCaMP/Revision/cytoGCaMP/FilesforR/Peaks_allMice_cyto_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
all.cyto.OT<-read.table("D:/Data/GCaMP_RCaMP/Revision/cytoGCaMP/FilesforR/OnsetTimes_allMice_cyto_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")


##########
# reshape data

lsm.options(pbkrtest.limit = 100000)


# remove genotype information 
lck.peaks2$Genotype<-NULL
lck.OT2$Genotype<-NULL
lck.OT1$Spot_trial<-NULL
lck.OT1$ROIs_trial<-NULL
lck.OT1$ROIType<-NULL
names(lck.OT1)[names(lck.OT1)=="Overlap"] <- "overlap"

# join data sets
all.lck.peaks<-rbind(lck.peaks1,lck.peaks2)
all.lck.OT<-rbind(lck.OT1, lck.OT2)
all.lck.OT.2s<-rbind(nostim.lck.OT.2s,longstim.lck.OT.2s)

# remove the data frames that are combined
rm(lck.peaks1,lck.peaks2,lck.OT1, lck.OT2, nostim.lck.OT.2s, longstim.lck.OT.2s)

# only consider wildtype mice, remove IP3R2KO
all.lck.peaks<-all.lck.peaks[!(grepl("IP",all.lck.peaks$Animal)),]
all.lck.OT<-all.lck.OT[!(grepl("IP",all.lck.OT$Animal)),]

#############
# ONSET Times
# unique ROI names
all.lck.OT$ROIs_trial<-paste(all.lck.OT$Animal, all.lck.OT$Spot, all.lck.OT$Trial, all.lck.OT$ROI, sep= "_")
all.lck.OT$Spot_trial<-paste(all.lck.OT$Animal, all.lck.OT$Spot, all.lck.OT$Trial, sep= "_")
all.lck.OT$Spot_trial_Cond<-paste(all.lck.OT$Spot_trial, all.lck.OT$Condition, sep="_")
all.lck.OT$ROIs_Cond<-paste(all.lck.OT$ROIs_trial, all.lck.OT$Condition, sep="_")
all.lck.OT$Animal_Spot<-paste(all.lck.OT$Animal, all.lck.OT$Spot, sep="_")

all.lck.OT.2s$Spot_trial_Cond<-paste(all.lck.OT.2s$Spot_trial, all.lck.OT.2s$Condition, sep="_")
all.lck.OT.2s$ROIs_Cond<-paste(all.lck.OT.2s$ROIs_trial, all.lck.OT.2s$Condition, sep="_")

# ROI Types
all.lck.OT$ROIType= "none"
all.lck.OTA<- subset(all.lck.OT, Channel=="GCaMP")
all.lck.OTB<- subset(all.lck.OT, Channel=="RCaMP")

# ROITypes
all.lck.OTA$ROIType[grepl("r",all.lck.OTA$ROI)]="Process"
all.lck.OTA$ROIType[grepl("E",all.lck.OTA$ROI)]="Endfoot"
all.lck.OTB$ROIType[grepl("r",all.lck.OTB$ROI)]="Dendrite"
all.lck.OTB$ROIType[grepl("D",all.lck.OTB$ROI)]="Dendrite"
all.lck.OTB$ROIType[grepl("N",all.lck.OTB$ROI)]="Neuron"
#all.lck.OTB$ROIType[grepl("S",all.lck.OTB$ROI)]="SmoothMuscle"

all.lck.OT<-rbind(all.lck.OTA, all.lck.OTB)
all.lck.OT$ROIType<- as.factor(all.lck.OT$ROIType)

# exclude the neuropil ROIs, because they were hand selected and not necessary
all.lck.OT<-all.lck.OT[!(all.lck.OT$ROIType=="none"),]

# REMOVE duplicate entries from onset time and peak time data frames 
# only the first entry will be used
all.lck.OT2=all.lck.OT[order(all.lck.OT$OnsetTime),] # sort by ascending onset time
all.lck.OT2.2s=all.lck.OT.2s[order(all.lck.OT.2s$OnsetTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.lck.OT<-distinct(all.lck.OT2, ROIs_Cond, .keep_all = TRUE)
all.lck.OT.2s<-distinct(all.lck.OT2.2s, ROIs_Cond, .keep_all = TRUE)

#adjust onset time for the data with 2 s before stimulation included:

all.lck.OT.2s$OnsetTimeAdjust<-all.lck.OT.2s$OnsetTime-2

rm(all.lck.OT2,all.lck.OT2.2s, all.lck.OTA, all.lck.OTB)


##############
# PEAK DATA
all.lck.peaks$ROIType= "none"
all.lck.peaksA<- subset(all.lck.peaks, Channel=="GCaMP")
all.lck.peaksB<- subset(all.lck.peaks, Channel=="RCaMP")

# ROITypes
all.lck.peaksA$ROIType[grepl("r",all.lck.peaksA$roiName)]="Process"
all.lck.peaksA$ROIType[grepl("E",all.lck.peaksA$roiName)]="Endfoot"
all.lck.peaksB$ROIType[grepl("r",all.lck.peaksB$roiName)]="Dendrite"
all.lck.peaksB$ROIType[grepl("D",all.lck.peaksB$roiName)]="Dendrite"
all.lck.peaksB$ROIType[grepl("N",all.lck.peaksB$roiName)]="Neuron"
#all.lck.peaksB$ROIType[grepl("S",all.lck.peaksB$roiName)]="SmoothMuscle"

all.lck.peaks<-rbind(all.lck.peaksA, all.lck.peaksB)
all.lck.peaks$ROIType<- as.factor(all.lck.peaks$ROIType)

#unique ROI names
all.lck.peaks$ROIs_trial<-paste(all.lck.peaks$Animal, all.lck.peaks$Spot, all.lck.peaks$Trial,all.lck.peaks$roiName, sep= "_")
all.lck.peaks$trials<-paste(all.lck.peaks$Animal, all.lck.peaks$Spot, all.lck.peaks$Trial, sep= "_")
all.lck.peaks$trials_Cond<-paste(all.lck.peaks$trials, all.lck.peaks$Condition, sep= "_")

all.lck.peaks$ROIs_Cond<-paste(all.lck.peaks$ROIs_trial, all.lck.peaks$Condition, sep= "_")
all.lck.peaks$Animal_Spot<-paste(all.lck.peaks$Animal, all.lck.peaks$Spot, sep= "_")


# count number of trials per spot
Spot.lck.ntrials<-ddply(all.lck.peaks, c("Animal","Spot","Condition"), summarise, nTrials=length(unique(Trial)))
Spot.lck.ntrials$Ani_Spot_Cond<-paste(Spot.lck.ntrials$Animal, Spot.lck.ntrials$Spot, Spot.lck.ntrials$Condition, sep="_")

Spot.TotalROIs<- ddply(all.lck.peaks[all.lck.peaks$Condition=="Stim",], c("trials","trials_Cond", "Channel"), summarise, nROIs=length(unique(ROIs_Cond)))

# if there were no ROIs with a signal during the trial, make it zero
noPeaks.stim<-subset(all.lck.peaks, peakType=="NoPeak" & Condition=="Stim")
Spot.noPeaks<- ddply(noPeaks.stim, c("trials_Cond", "Channel"), summarise, nROIs=length(unique(ROIs_Cond)))


# remove ROIs with no peaks
all.lck.peaks<-all.lck.peaks[!(all.lck.peaks$peakType=="NoPeak"),]

# exclude the neuropil ROIs, because they were hand selected and not necessary
all.lck.peaks<-all.lck.peaks[!(all.lck.peaks$ROIType=="none"),]

# remove matching astrocyte process and soma ROIs
Overlap= all.lck.peaks$overlap!=0
all.lck.peaks<-all.lck.peaks[!Overlap,]
#OverlapROIs<-unique(nostim$ROIs_trial[Overlap])


#add baseline time to peaks table
all.lck.OT$BL_time<-all.lck.OT$baseline/all.lck.OT$FrameRate
all.lck.peaks<-merge(all.lck.peaks, all.lck.OT[, c("ROIs_Cond", "BL_time")], by="ROIs_Cond", all.x=TRUE)

# adjust peak time and duration
all.lck.peaks$peakTime<- all.lck.peaks$peakTime-all.lck.peaks$BL_time
all.lck.peaks$peakStart<- all.lck.peaks$peakStart-all.lck.peaks$BL_time
all.lck.peaks$peakStartHalf<- all.lck.peaks$peakStartHalf-all.lck.peaks$BL_time
all.lck.peaks$Duration<- all.lck.peaks$halfWidth*2



# drop peaks that occur before the start of stimulation
all.lck.peaks2<-subset(all.lck.peaks,peakTime>0)

# only the first entry will be used
all.lck.peaks3<-all.lck.peaks2[order(all.lck.peaks2$peakTime),] # sort by ascending peak time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.lck.peaks<-distinct(all.lck.peaks3, ROIs_Cond,.keep_all = TRUE)


rm(all.lck.peaks2,all.lck.peaks3, all.lck.peaksA, all.lck.peaksB)







################### 
# cytoGCaMP6s data
all.cyto.peaks$Genotype<-NULL
all.cyto.OT$Genotype<- NULL
all.cyto.OT.2s<-rbind(nostim.cyto.OT.2s,longstim.cyto.OT.2s)

# remove the data frames that are combined
rm(nostim.cyto.OT.2s, longstim.cyto.OT.2s)

## ONSET TIMES

# unique ROI names
all.cyto.OT$ROIs_trial<-paste(all.cyto.OT$Animal, all.cyto.OT$Spot, all.cyto.OT$Trial, all.cyto.OT$ROI, sep= "_")
all.cyto.OT$Spot_trial<-paste(all.cyto.OT$Animal, all.cyto.OT$Spot, all.cyto.OT$Trial, sep= "_")
all.cyto.OT$Spot_trial_Cond<-paste(all.cyto.OT$Spot_trial, all.cyto.OT$Condition, sep="_")
all.cyto.OT$ROIs_Cond<-paste(all.cyto.OT$ROIs_trial, all.cyto.OT$Condition, sep="_")
all.cyto.OT$Animal_Spot<-paste(all.cyto.OT$Animal, all.cyto.OT$Spot, sep="_")

all.cyto.OT.2s$Spot_trial_Cond<-paste(all.cyto.OT.2s$Spot_trial, all.cyto.OT.2s$Condition, sep="_")
all.cyto.OT.2s$ROIs_Cond<-paste(all.cyto.OT.2s$ROIs_trial, all.cyto.OT.2s$Condition, sep="_")

# ROI Types
all.cyto.OT$ROIType= "none"
all.cyto.OTA<- subset(all.cyto.OT, Channel=="GCaMP")
all.cyto.OTB<- subset(all.cyto.OT, Channel=="RCaMP")

# ROITypes
all.cyto.OTA$ROIType[grepl("r",all.cyto.OTA$ROI)]="Process"
all.cyto.OTA$ROIType[grepl("E",all.cyto.OTA$ROI)]="Endfoot"
all.cyto.OTB$ROIType[grepl("r",all.cyto.OTB$ROI)]="Dendrite"
all.cyto.OTB$ROIType[grepl("D",all.cyto.OTB$ROI)]="Dendrite"
all.cyto.OTB$ROIType[grepl("N",all.cyto.OTB$ROI)]="Neuron"
#all.cyto.OTB$ROIType[grepl("S",all.cyto.OTB$ROI)]="SmoothMuscle"

all.cyto.OT<-rbind(all.cyto.OTA, all.cyto.OTB)
all.cyto.OT$ROIType<- as.factor(all.cyto.OT$ROIType)

# REMOVE duplicate entries from onset time and peak time data frames 
# only the first entry will be used
all.cyto.OT2=all.cyto.OT[order(all.cyto.OT$OnsetTime),] # sort by ascending onset time
all.cyto.OT2.2s=all.cyto.OT.2s[order(all.cyto.OT.2s$OnsetTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.cyto.OT<-distinct(all.cyto.OT2, ROIs_Cond, .keep_all = TRUE)
all.cyto.OT.2s<-distinct(all.cyto.OT2.2s, ROIs_Cond, .keep_all = TRUE)

#adjust onset time for the data with 2 s before stimulation included:

all.cyto.OT.2s$OnsetTimeAdjust<-all.cyto.OT.2s$OnsetTime-2

rm(all.cyto.OT2,all.cyto.OT2.2s, all.cyto.OTA, all.cyto.OTB)




## PEAK DATA

# exclude the neuropil ROIs, because they were hand selected and not necessary
all.cyto.peaks<-all.cyto.peaks[!(all.cyto.peaks$roiName=="np"),]
all.cyto.peaks<-all.cyto.peaks[!(all.cyto.peaks$peakType=="NoPeak"),]

# no stim peak data
all.cyto.peaks$ROIType= 0
all.cyto.peaksA<- subset(all.cyto.peaks, Channel=="GCaMP")
all.cyto.peaksB<- subset(all.cyto.peaks, Channel=="RCaMP")

# ROITypes
all.cyto.peaksA$ROIType[grepl("r",all.cyto.peaksA$roiName)]="Process"
all.cyto.peaksA$ROIType[grepl("E",all.cyto.peaksA$roiName)]="Endfoot"
all.cyto.peaksA$ROIType[grepl("S",all.cyto.peaksA$roiName)]="Soma"
all.cyto.peaksB$ROIType[grepl("r",all.cyto.peaksB$roiName)]="Dendrite"
all.cyto.peaksB$ROIType[grepl("D",all.cyto.peaksB$roiName)]="Dendrite"
all.cyto.peaksB$ROIType[grepl("N",all.cyto.peaksB$roiName)]="Neuron"

all.cyto.peaks<-rbind(all.cyto.peaksA, all.cyto.peaksB)
all.cyto.peaks$ROIType<- as.factor(all.cyto.peaks$ROIType)

#unique ROI names
all.cyto.peaks$ROIs_trial<-paste(all.cyto.peaks$Animal, all.cyto.peaks$Spot, all.cyto.peaks$Trial,all.cyto.peaks$roiName, sep= "_")
all.cyto.peaks$trials<-paste(all.cyto.peaks$Animal, all.cyto.peaks$Spot, all.cyto.peaks$Trial, sep= "_")
all.cyto.peaks$trials_Cond<-paste(all.cyto.peaks$trials, all.cyto.peaks$Condition, sep= "_")
all.cyto.peaks$ROIs_Cond<-paste(all.cyto.peaks$ROIs_trial, all.cyto.peaks$Condition, sep= "_")


# count number of trials per spot
Spot.cyto.ntrials<-ddply(all.cyto.peaks, c("Animal","Spot","Condition"), summarise, nTrials=length(unique(Trial)))
Spot.cyto.ntrials$Ani_Spot_Cond<-paste(Spot.cyto.ntrials$Animal, Spot.cyto.ntrials$Spot, Spot.cyto.ntrials$Condition, sep="_")


# remove matching astrocyte process and soma ROIs
Overlap= all.cyto.peaks$overlap!=0
all.cyto.peaks<-all.cyto.peaks[!Overlap,]
#OverlapROIs<-unique(nostim$ROIs_trial[Overlap])



# adjust peak time and duration
all.cyto.OT$BL_time<-all.cyto.OT$baseline/all.cyto.OT$FrameRate
all.cyto.peaks<-merge(all.cyto.peaks, all.cyto.OT[, c("ROIs_Cond", "BL_time")], by="ROIs_Cond", all.x=TRUE)

all.cyto.peaks$peakTime<- all.cyto.peaks$peakTime-all.cyto.peaks$BL_time
all.cyto.peaks$peakStart<- all.cyto.peaks$peakStart-all.cyto.peaks$BL_time
all.cyto.peaks$peakStartHalf<- all.cyto.peaks$peakStartHalf-all.cyto.peaks$BL_time
all.cyto.peaks$Duration<- all.cyto.peaks$halfWidth*2


# drop peaks that occur before the start of stimulation
all.cyto.peaks2<-subset(all.cyto.peaks,peakTime>0)

# only the first entry will be used
all.cyto.peaks3<-all.cyto.peaks2[order(all.cyto.peaks2$peakTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.cyto.peaks<-distinct(all.cyto.peaks3, ROIs_Cond,.keep_all = TRUE)

rm(all.cyto.peaks2,all.cyto.peaks3, all.cyto.peaksA, all.cyto.peaksB)



################
# NEURONAL PERCENTILE OF ONSET

# combine neuronal responses from cyto-GCaMP and Lck-GCaMP trials
# neuronal responses to stimulation
NeuronalStim1<-subset(all.lck.OT, Channel=="RCaMP" & Condition=="Stim")
NeuronalStim2<-subset(all.cyto.OT, Channel=="RCaMP" & Condition=="Stim")
all.NeuronalStim<-rbind(NeuronalStim1, NeuronalStim2)

ggplot(all.NeuronalStim,aes(x=OnsetTime)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("all RCaMP onset times from stim trials")+
  max.theme

Neuron95Onset<-quantile(all.NeuronalStim$OnsetTime[all.NeuronalStim$OnsetTime<8], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
print(Neuron95Onset)

NeuronOT25<-Neuron95Onset[[6]]
NeuronOT50<-Neuron95Onset[[11]]
NeuronOT75<-Neuron95Onset[[16]]

# should have an onset time in 8 s stimulus

ggplot(all.NeuronalStim[all.NeuronalStim$OnsetTime<8,],aes(x=OnsetTime)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP onset times between 0 and 8 s from stim trials")+
  geom_vline(xintercept=NeuronOT25) +
  geom_vline(xintercept=NeuronOT50) + 
  geom_vline(xintercept=NeuronOT75) + 
  max.theme


rm(all.NeuronalStim, NeuronalStim1, NeuronalStim2)

###
# NEURONAL PERCENTILE OF Peak TIMES

# combine neuronal responses from cyto-GCaMP and Lck-GCaMP trials
# neuronal responses to stimulation
NeuronalStim1<-subset(all.lck.peaks, Channel=="RCaMP" & Condition=="Stim")
NeuronalStim2<-subset(all.cyto.peaks, Channel=="RCaMP" & Condition=="Stim")
all.NeuronalStim.peaks<-rbind(NeuronalStim1, NeuronalStim2)

ggplot(all.NeuronalStim.peaks[all.NeuronalStim.peaks$peakTime<40,],aes(x=peakTime)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("all RCaMP peak times from stim trials")+
  max.theme

Neuron95PeakTime<-quantile(all.NeuronalStim.peaks$peakTime[all.NeuronalStim.peaks$peakTime<40], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
print(Neuron95PeakTime)

NeuronPT25<-Neuron95Onset[[6]]
NeuronPT50<-Neuron95Onset[[11]]
NeuronPT75<-Neuron95Onset[[16]]

# should have a peak time in 8 s stimulus

ggplot(all.NeuronalStim.peaks[all.NeuronalStim.peaks$peakTime<8,],aes(x=peakTime)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP onset times between 0 and 8 s from stim trials")+
  max.theme

rm(all.NeuronalStim.peaks)

#######
# ASTROCYTE PERCENTILE ONSET
# Lck Astrocytes

# subset data to include only astrocytes, stim trials, with an onset greater than 0 and less than 40
AstroStim.Lck<-subset(all.lck.OT, Channel=="GCaMP" & Condition=="Stim" & OnsetTime<40 & OnsetTime>0)

ggplot(AstroStim.Lck,aes(x=OnsetTime, group=Condition, fill=Condition)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("0-15 s Lck.GCaMP onset times from stim trials")+
  max.theme

# data with onset time= 0s has been excluded because it screws up the log transformation
AstroStim.Lck$Log_OT<-log(AstroStim.Lck$OnsetTime)

ggplot(AstroStim.Lck,aes(x=Log_OT, group=Condition, fill=Condition)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("0-15 s Lck.GCaMP onset times from stim trials- log scale")+
  max.theme

stD.logOT<-sd(AstroStim.Lck$Log_OT,na.rm = TRUE)
mean.logOT<-mean(AstroStim.Lck$Log_OT,na.rm = TRUE)

threshold.Lck.OT<-exp(mean.logOT)+(2.33*exp(stD.logOT))



# ASTROCYTE PERCENTILE PEAK TIME
# Lck Astrocytes
AstroStim.Lck.PT<-subset(all.lck.peaks, Channel=="GCaMP" & Condition=="Stim" & peakTime<40 & peakTime>0)

ggplot(AstroStim.Lck.PT,aes(x=peakTime)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("0-20 s Lck.GCaMP peak times from stim trials")+
  max.theme

AstroStim.Lck.PT$Log_PT<-log(AstroStim.Lck.PT$peakTime)

ggplot(AstroStim.Lck.PT,aes(x=Log_PT, group=Condition, fill=Condition)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("0-15 s Lck.GCaMP peak times from stim trials- log scale")+
  max.theme

stD.logPT<-sd(AstroStim.Lck.PT$Log_PT,na.rm = TRUE)
mean.logPT<-mean(AstroStim.Lck.PT$Log_PT,na.rm = TRUE)

threshold.Lck.PT<-exp(mean.logPT)+(2.33*exp(stD.logPT))



# combine Onset times and peak times to allow plotting on the same histogram
AstroStim.Lck$TimeType="Onset"
AstroStim.Lck.PT$TimeType="PeakMax"

Lck.Onsettime<-AstroStim.Lck[,c(12,22)]
Lck.Peaktime<-AstroStim.Lck.PT[,c(6,31)]

names(Lck.Onsettime)[names(Lck.Onsettime)=="OnsetTime"] <- "time"
names(Lck.Peaktime)[names(Lck.Peaktime)=="peakTime"] <- "time"

Lck.combinedTimes<-rbind(Lck.Onsettime,Lck.Peaktime)

ggplot(Lck.combinedTimes[Lck.combinedTimes$time<40,],aes(x=time, fill=TimeType)) +
  geom_histogram(binwidth=(0.084*5), position="dodge") +
  ggtitle("Lck Onset and peak time distributions from stim trials, 0.01 log transformed Threshold")+
  geom_vline(xintercept=threshold.Lck.OT) +
  geom_vline(xintercept=threshold.Lck.PT) + 
  max.theme

rm(AstroStim.Lck.PT, AstroStim.Lck, Lck.Onsettime,Lck.Peaktime,Lck.combinedTimes)


##################
# cyto Astrocytes
# subset data to include only astrocytes, stim trials, with an onset greater than 0 and less than 40
AstroStim.cyto<-subset(all.cyto.OT, Channel=="GCaMP" & Condition=="Stim" & OnsetTime<40 & OnsetTime>0)

ggplot(AstroStim.cyto,aes(x=OnsetTime, group=Condition, fill=Condition)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("0-15 s cyto.GCaMP onset times from stim trials")+
  max.theme

# data with onset time= 0s has been excluded because it screws up the log transformation
AstroStim.cyto$Log_OT<-log(AstroStim.cyto$OnsetTime)

ggplot(AstroStim.cyto,aes(x=Log_OT, group=Condition, fill=Condition)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("0-15 s cyto.GCaMP onset times from stim trials- log scale")+
  max.theme

stD.logOT<-sd(AstroStim.cyto$Log_OT,na.rm = TRUE)
mean.logOT<-mean(AstroStim.cyto$Log_OT,na.rm = TRUE)

threshold.cyto.OT<-exp(mean.logOT)+(2.33*exp(stD.logOT))


# ASTROCYTE PERCENTILE PEAK TIME
# cyto Astrocytes
AstroStim.cyto.PT<-subset(all.cyto.peaks, Channel=="GCaMP" & Condition=="Stim" & peakTime<40 & peakTime>0)

ggplot(AstroStim.cyto.PT,aes(x=peakTime)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("0-20 s cyto.GCaMP peak times from stim trials")+
  max.theme

AstroStim.cyto.PT$Log_PT<-log(AstroStim.cyto.PT$peakTime)

ggplot(AstroStim.cyto.PT,aes(x=Log_PT, group=Condition, fill=Condition)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("0-15 s cyto.GCaMP peak times from stim trials- log scale")+
  max.theme

stD.logPT<-sd(AstroStim.cyto.PT$Log_PT,na.rm = TRUE)
mean.logPT<-mean(AstroStim.cyto.PT$Log_PT,na.rm = TRUE)

threshold.cyto.PT<-exp(mean.logPT)+(2.33*exp(stD.logPT))


# combine Onset times and peak times to allow plotting on the same histogram
AstroStim.cyto$TimeType="Onset"
AstroStim.cyto.PT$TimeType="PeakMax"

cyto.Onsettime<-AstroStim.cyto[,c(12,22)]
cyto.Peaktime<-AstroStim.cyto.PT[,c(6,31)]

names(cyto.Onsettime)[names(cyto.Onsettime)=="OnsetTime"] <- "time"
names(cyto.Peaktime)[names(cyto.Peaktime)=="peakTime"] <- "time"

cyto.combinedTimes<-rbind(cyto.Onsettime,cyto.Peaktime)

ggplot(cyto.combinedTimes[cyto.combinedTimes$time<40,],aes(x=time, fill=TimeType)) +
  geom_histogram(binwidth=(0.084*5), position="dodge") +
  ggtitle("cyto Onset and peak time distributions from stim trials, 80th per")+
  geom_vline(xintercept=threshold.cyto.OT) +
  geom_vline(xintercept=threshold.cyto.PT) + 
  max.theme

rm(AstroStim.cyto.PT, AstroStim.cyto, cyto.Onsettime, cyto.Peaktime, cyto.combinedTimes)
######

## NOTE: define a stimulation window
# for distributions: histograms

stimwindow=15    # chosen because most of the astrocyte responses were done by this point, only for distributions


stim.lck.OT.dist<-subset(all.lck.OT,Condition!="shortstim" & OnsetTime<stimwindow)
stim.lck.OT.dist.2s<-subset(all.lck.OT.2s,Condition!="shortstim" & OnsetTime<stimwindow)

stim.cyto.OT.dist<-subset(all.cyto.OT,Condition!="shortstim" & OnsetTime<stimwindow)
stim.cyto.OT.dist.2s<-subset(all.cyto.OT.2s,Condition!="shortstim" & OnsetTime<stimwindow)


# LCK DATA
# Onset time histograms- normalized to the number of trials

# long stim vs no stim
ntrials.lck.OT.long.R.dis<- ddply(stim.lck.OT.dist[stim.lck.OT.dist$Channel=="RCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))
ntrials.lck.OT.long.G.dis<- ddply(stim.lck.OT.dist[stim.lck.OT.dist$Channel=="GCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))


#histogram bins
histseq= seq(0,15,0.5)
Nostim.N=0
Stim.N=0
Nostim.A=0
Stim.A=0
zeroRow<-data.frame(cbind(Nostim.N, Stim.N,Nostim.A, Stim.A))

# neuronal lck onset histogram
# counts for each condition in the histogram
Nostim.N=hist(stim.lck.OT.dist$OnsetTime[(stim.lck.OT.dist$Channel=="RCaMP" & stim.lck.OT.dist$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim.N=hist(stim.lck.OT.dist$OnsetTime[(stim.lck.OT.dist$Channel=="RCaMP" & stim.lck.OT.dist$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts

Nostim.A=hist(stim.lck.OT.dist$OnsetTime[(stim.lck.OT.dist$Channel=="GCaMP" & stim.lck.OT.dist$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim.A=hist(stim.lck.OT.dist$OnsetTime[(stim.lck.OT.dist$Channel=="GCaMP" & stim.lck.OT.dist$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim.N=Nostim.N/ntrials.lck.OT.long.R.dis$ntrials[(ntrials.lck.OT.long.R.dis$Condition=="Nostim")]
Stim.N=Stim.N/ntrials.lck.OT.long.R.dis$ntrials[(ntrials.lck.OT.long.R.dis$Condition=="Stim")]

Nostim.A=Nostim.A/ntrials.lck.OT.long.G.dis$ntrials[(ntrials.lck.OT.long.G.dis$Condition=="Nostim")]
Stim.A=Stim.A/ntrials.lck.OT.long.G.dis$ntrials[(ntrials.lck.OT.long.G.dis$Condition=="Stim")]

#Shortstim=hist(all.lck.OT$OnsetTime[(all.lck.OT$Channel=="RCaMP" & all.lck.OT$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

#Shortstim=Shortstim/ntrials.lck.OT$ntrials[(ntrials.lck.OT$Condition=="shortstim")]

#make a data frame for plotting
lck.long.histo <- data.frame(cbind(Nostim.N, Stim.N,Nostim.A, Stim.A))
lck.long.histo2<-rbind(zeroRow,lck.long.histo)
lck.long.histo2$time<-histseq

ggplot(NULL, aes(x=time))+
  geom_line(data=lck.long.histo2, aes(y=Nostim.N, color="Nostim.N")) +
  geom_line(data=lck.long.histo2, aes(y=Stim.N, color="Stim.N")) +
  geom_line(data=lck.long.histo2, aes(y=Nostim.A*7, color="Nostim.A")) +
  geom_line(data=lck.long.histo2, aes(y=Stim.A*7, color="Stim.A")) +
  scale_y_continuous(sec.axis = sec_axis(~./7, name = "Astrocyte peaks/trial")) + # secondary axis
  ggtitle("stim- neurons vs astrocytes-lck data") + 
  xlab("Onset Time (s)") + 
  ylab("Neuron peaks/trial") + 
  xlim(-0.5,15) +
  max.theme



#compare distributions (what is plotted in the figure)

# ks test- astrocytes
OT.lck.GC.NSvsS.kstest<- ks.test(Nostim.A,Stim.A)
print(OT.lck.GC.NSvsS.kstest)

# neurons
OT.lck.RC.NSvsS.kstest<- ks.test(Nostim.N,Stim.N)
print(OT.lck.RC.NSvsS.kstest)



# Anderson Darling
library("kSamples")
#library("SuppDist")

OT.lck.GC.NSvsS.adtest<- ad.test(Nostim.A,Stim.A)
print(OT.lck.GC.NSvsS.adtest)

# neurons
OT.lck.RC.NSvsS.adtest<- ad.test(Nostim.N,Stim.N)
print(OT.lck.RC.NSvsS.adtest)

rm(lck.long.histo, lck.long.histo2)

######
# Lck peak time distributions

stimwindow=15
histseq= seq(-5,15,0.5)

stim.lck.PT.dist<-subset(all.lck.peaks,Condition!="shortstim" & peakTime<stimwindow)
short.lck.PT.dist<-subset(all.lck.peaks,Condition!="Stim" & peakTime<stimwindow)

# long stim vs no stim
ntrials.lck.PT.long.R.dis<- ddply(stim.lck.PT.dist[stim.lck.PT.dist$Channel=="RCaMP",], c("Condition"), summarise, ntrials=length(unique(trials)))
ntrials.lck.PT.long.G.dis<- ddply(stim.lck.PT.dist[stim.lck.PT.dist$Channel=="GCaMP",], c("Condition"), summarise, ntrials=length(unique(trials)))

# neuronal lck onset histogram
# counts for each condition in the histogram
Nostim.N=hist(stim.lck.PT.dist$peakTime[(stim.lck.PT.dist$Channel=="RCaMP" & stim.lck.PT.dist$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim.N=hist(stim.lck.PT.dist$peakTime[(stim.lck.PT.dist$Channel=="RCaMP" & stim.lck.PT.dist$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts

Nostim.A=hist(stim.lck.PT.dist$peakTime[(stim.lck.PT.dist$Channel=="GCaMP" & stim.lck.PT.dist$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim.A=hist(stim.lck.PT.dist$peakTime[(stim.lck.PT.dist$Channel=="GCaMP" & stim.lck.PT.dist$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim.N=Nostim.N/ntrials.lck.PT.long.R.dis$ntrials[(ntrials.lck.PT.long.R.dis$Condition=="Nostim")]
Stim.N=Stim.N/ntrials.lck.PT.long.R.dis$ntrials[(ntrials.lck.PT.long.R.dis$Condition=="Stim")]

Nostim.A=Nostim.A/ntrials.lck.PT.long.G.dis$ntrials[(ntrials.lck.PT.long.G.dis$Condition=="Nostim")]
Stim.A=Stim.A/ntrials.lck.PT.long.G.dis$ntrials[(ntrials.lck.PT.long.G.dis$Condition=="Stim")]

#Shortstim=hist(all.lck.OT$OnsetTime[(all.lck.OT$Channel=="RCaMP" & all.lck.OT$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

#Shortstim=Shortstim/ntrials.lck.OT$ntrials[(ntrials.lck.OT$Condition=="shortstim")]

#make a data frame for plotting
lck.long.PT.histo <- data.frame(cbind(Nostim.N, Stim.N,Nostim.A, Stim.A))
lck.long.PT.histo2<-rbind(zeroRow,lck.long.PT.histo)
lck.long.PT.histo2$time<-histseq

ggplot(NULL, aes(x=time))+
  geom_line(data=lck.long.PT.histo2, aes(y=Nostim.N, color="Nostim.N")) +
  geom_line(data=lck.long.PT.histo2, aes(y=Stim.N, color="Stim.N")) +
  geom_line(data=lck.long.PT.histo2, aes(y=Nostim.A*2.5, color="Nostim.A")) +
  geom_line(data=lck.long.PT.histo2, aes(y=Stim.A*2.5, color="Stim.A")) +
  scale_y_continuous(sec.axis = sec_axis(~./2.5, name = "Astrocyte peaks/trial")) + # secondary axis
  ggtitle("stim- neurons vs astrocytes-lck data") + 
  xlab("Peak Max Time (s)") + 
  ylab("Neuron peaks/trial") + 
  xlim(-0.5,15) +
  max.theme


#compare distributions (what is plotted in the figure)

# ks test- astrocytes
pT.lck.GC.NSvsS.kstest<- ks.test(Nostim.A,Stim.A)
print(pT.lck.GC.NSvsS.kstest)

# neurons
pT.lck.RC.NSvsS.kstest<- ks.test(Nostim.N,Stim.N)
print(pT.lck.RC.NSvsS.kstest)



# Anderson Darling
library("kSamples")
#library("SuppDist")

pT.lck.GC.NSvsS.adtest<- ad.test(Nostim.A,Stim.A)
print(pT.lck.GC.NSvsS.adtest)

# neurons
pT.lck.RC.NSvsS.adtest<- ad.test(Nostim.N,Stim.N)
print(pT.lck.RC.NSvsS.adtest)

rm(lck.long.PT.histo, lck.long.PT.histo2)
##################### 
histseq= seq(0,15,0.5)
#cyto data  ONSET time distributions
ntrials.cyto.OT.long.R.dis<- ddply(stim.cyto.OT.dist[stim.cyto.OT.dist$Channel=="RCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))
ntrials.cyto.OT.long.G.dis<- ddply(stim.cyto.OT.dist[stim.cyto.OT.dist$Channel=="GCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))

# neuronal cyto onset histogram
# counts for each condition in the histogram
Nostim.N=hist(stim.cyto.OT.dist$OnsetTime[(stim.cyto.OT.dist$Channel=="RCaMP" & stim.cyto.OT.dist$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim.N=hist(stim.cyto.OT.dist$OnsetTime[(stim.cyto.OT.dist$Channel=="RCaMP" & stim.cyto.OT.dist$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts

Nostim.A=hist(stim.cyto.OT.dist$OnsetTime[(stim.cyto.OT.dist$Channel=="GCaMP" & stim.cyto.OT.dist$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim.A=hist(stim.cyto.OT.dist$OnsetTime[(stim.cyto.OT.dist$Channel=="GCaMP" & stim.cyto.OT.dist$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim.N=Nostim.N/ntrials.cyto.OT.long.R.dis$ntrials[(ntrials.cyto.OT.long.R.dis$Condition=="Nostim")]
Stim.N=Stim.N/ntrials.cyto.OT.long.R.dis$ntrials[(ntrials.cyto.OT.long.R.dis$Condition=="Stim")]

Nostim.A=Nostim.A/ntrials.cyto.OT.long.G.dis$ntrials[(ntrials.cyto.OT.long.G.dis$Condition=="Nostim")]
Stim.A=Stim.A/ntrials.cyto.OT.long.G.dis$ntrials[(ntrials.cyto.OT.long.G.dis$Condition=="Stim")]

#Shortstim=hist(all.lck.OT$OnsetTime[(all.lck.OT$Channel=="RCaMP" & all.lck.OT$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts
#Shortstim=Shortstim/ntrials.lck.OT$ntrials[(ntrials.lck.OT$Condition=="shortstim")]

#make a data frame for plotting
cyto.long.histo <- data.frame(cbind(Nostim.N, Stim.N,Nostim.A, Stim.A))
cyto.long.histo2<-rbind(zeroRow,cyto.long.histo)
cyto.long.histo2$time<-histseq

ggplot(NULL, aes(x=time))+
  geom_line(data=cyto.long.histo2, aes(y=Nostim.N, color="Nostim.N")) +
  geom_line(data=cyto.long.histo2, aes(y=Stim.N, color="Stim.N")) +
  geom_line(data=cyto.long.histo2, aes(y=Nostim.A*1.5, color="Nostim.A")) +
  geom_line(data=cyto.long.histo2, aes(y=Stim.A*1.5, color="Stim.A")) +
  scale_y_continuous(sec.axis = sec_axis(~./1.5, name = "Astrocyte peaks/trial")) + # secondary axis
  ggtitle("stim- neurons vs astrocytes-cyto data") + 
  xlab("Onset Time (s)") + 
  ylab("Neuron peaks/trial") + 
  xlim(-0.5,15) +
  max.theme



#compare distributions (what is plotted in the figure)

# ks test- astrocytes
OT.cyto.GC.NSvsS.kstest<- ks.test(Nostim.A,Stim.A)
print(OT.cyto.GC.NSvsS.kstest)

# neurons
OT.cyto.RC.NSvsS.kstest<- ks.test(Nostim.N,Stim.N)
print(OT.cyto.RC.NSvsS.kstest)



# Anderson Darling
library("kSamples")
#library("SuppDist")

OT.cyto.GC.NSvsS.adtest<- ad.test(Nostim.A,Stim.A)
print(OT.cyto.GC.NSvsS.adtest)

# neurons
OT.cyto.RC.NSvsS.adtest<- ad.test(Nostim.N,Stim.N)
print(OT.cyto.RC.NSvsS.adtest)

rm(cyto.long.histo, cyto.long.histo2)
#####

# cyto peak time distributions

stimwindow=15
histseq= seq(-5,15,0.5)

stim.cyto.PT.dist<-subset(all.cyto.peaks,Condition!="shortstim" & peakTime<stimwindow)
short.cyto.PT.dist<-subset(all.cyto.peaks,Condition!="Stim" & peakTime<stimwindow)

# long stim vs no stim
ntrials.cyto.PT.long.R.dis<- ddply(stim.cyto.PT.dist[stim.cyto.PT.dist$Channel=="RCaMP",], c("Condition"), summarise, ntrials=length(unique(trials)))
ntrials.cyto.PT.long.G.dis<- ddply(stim.cyto.PT.dist[stim.cyto.PT.dist$Channel=="GCaMP",], c("Condition"), summarise, ntrials=length(unique(trials)))

# neuronal cyto onset histogram
# counts for each condition in the histogram
Nostim.N=hist(stim.cyto.PT.dist$peakTime[(stim.cyto.PT.dist$Channel=="RCaMP" & stim.cyto.PT.dist$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim.N=hist(stim.cyto.PT.dist$peakTime[(stim.cyto.PT.dist$Channel=="RCaMP" & stim.cyto.PT.dist$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts

Nostim.A=hist(stim.cyto.PT.dist$peakTime[(stim.cyto.PT.dist$Channel=="GCaMP" & stim.cyto.PT.dist$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim.A=hist(stim.cyto.PT.dist$peakTime[(stim.cyto.PT.dist$Channel=="GCaMP" & stim.cyto.PT.dist$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim.N=Nostim.N/ntrials.cyto.PT.long.R.dis$ntrials[(ntrials.cyto.PT.long.R.dis$Condition=="Nostim")]
Stim.N=Stim.N/ntrials.cyto.PT.long.R.dis$ntrials[(ntrials.cyto.PT.long.R.dis$Condition=="Stim")]

Nostim.A=Nostim.A/ntrials.cyto.PT.long.G.dis$ntrials[(ntrials.cyto.PT.long.G.dis$Condition=="Nostim")]
Stim.A=Stim.A/ntrials.cyto.PT.long.G.dis$ntrials[(ntrials.cyto.PT.long.G.dis$Condition=="Stim")]

#Shortstim=hist(all.cyto.OT$OnsetTime[(all.cyto.OT$Channel=="RCaMP" & all.cyto.OT$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

#Shortstim=Shortstim/ntrials.cyto.OT$ntrials[(ntrials.cyto.OT$Condition=="shortstim")]

#make a data frame for plotting
cyto.long.PT.histo <- data.frame(cbind(Nostim.N, Stim.N,Nostim.A, Stim.A))
cyto.long.PT.histo2<-rbind(zeroRow,cyto.long.PT.histo)
cyto.long.PT.histo2$time<-histseq

ggplot(NULL, aes(x=time))+
  geom_line(data=cyto.long.PT.histo2, aes(y=Nostim.N, color="Nostim.N")) +
  geom_line(data=cyto.long.PT.histo2, aes(y=Stim.N, color="Stim.N")) +
  geom_line(data=cyto.long.PT.histo2, aes(y=Nostim.A*2, color="Nostim.A")) +
  geom_line(data=cyto.long.PT.histo2, aes(y=Stim.A*2, color="Stim.A")) +
  scale_y_continuous(sec.axis = sec_axis(~./2, name = "Astrocyte peaks/trial")) + # secondary axis
  ggtitle("stim- neurons vs astrocytes-cyto data") + 
  xlab("Peak Max Time (s)") + 
  ylab("Neuron peaks/trial") + 
  xlim(-0.5,15) +
  max.theme


#compare distributions (what is plotted in the figure)

# ks test- astrocytes
pT.cyto.GC.NSvsS.kstest<- ks.test(Nostim.A,Stim.A)
print(pT.cyto.GC.NSvsS.kstest)

# neurons
pT.cyto.RC.NSvsS.kstest<- ks.test(Nostim.N,Stim.N)
print(pT.cyto.RC.NSvsS.kstest)



# Anderson Darling
library("kSamples")
#library("SuppDist")

pT.cyto.GC.NSvsS.adtest<- ad.test(Nostim.A,Stim.A)
print(pT.cyto.GC.NSvsS.adtest)

# neurons
pT.cyto.RC.NSvsS.adtest<- ad.test(Nostim.N,Stim.N)
print(pT.cyto.RC.NSvsS.adtest)

rm(cyto.long.PT.histo, cyto.long.PT.histo2, ntrials.cyto.OT.long.G.dis,ntrials.cyto.OT.long.R.dis,ntrials.cyto.PT.long.G.dis,
   ntrials.cyto.PT.long.R.dis, ntrials.lck.OT.long.G.dis, ntrials.lck.OT.long.R.dis, ntrials.lck.PT.long.G.dis, ntrials.lck.PT.long.R.dis)

##################
######
# number of ROIs in each trial for each field of view (across the whole trial)
lck.OT.8strial<-all.lck.OT[all.lck.OT$OnsetTime<8,]

lck.OT.8strial$Channel <- factor(lck.OT.8strial$Channel, levels = c("RCaMP","GCaMP"))

ROInum.lck.8strial<-ddply(lck.OT.8strial, c("Animal","Spot","Condition","Channel"), summarise, nROIs=length(OnsetTime))

ROInum.lck.8strial<-ROInum.lck.8strial[complete.cases(ROInum.lck.8strial),]
ROInum.lck.8strial$Ani_Spot<- paste(ROInum.lck.8strial$Animal, ROInum.lck.8strial$Spot, sep="_")

# add in number of trials
ROInum.lck.8strial$Ani_Spot_Cond<-paste(ROInum.lck.8strial$Animal, ROInum.lck.8strial$Spot, ROInum.lck.8strial$Condition, sep="_")
ROInum.lck.8strial<-merge(ROInum.lck.8strial, Spot.lck.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)
ROInum.lck.8strial$ROIsPerTrial<-ROInum.lck.8strial$nROIs/ROInum.lck.8strial$nTrials

# mean
df.lck.ROInum.8strial<-summarySE(ROInum.lck.8strial, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition"))

# plots
#boxplot
ggplot(ROInum.lck.8strial, aes(x=Channel,y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("Lck- ROIs per FOV") +
  max.theme

# paired line plots
df.lck.ROInum.8strial$ROIsPerTrialMean<-df.lck.ROInum.8strial$ROIsPerTrial
df.lck.ROInum.8strial$Chan_Cond<-interaction(df.lck.ROInum.8strial$Channel, df.lck.ROInum.8strial$Condition)
ROInum.lck.8strial$Chan_Cond<-interaction(ROInum.lck.8strial$Channel, ROInum.lck.8strial$Condition)

ROInum.lck.8strial<-merge(ROInum.lck.8strial, df.lck.ROInum.8strial[, c("Chan_Cond", "ROIsPerTrialMean","se")], by="Chan_Cond", all.x=TRUE)

ggplot(ROInum.lck.8strial[ROInum.lck.8strial$Channel=="RCaMP",], aes(x=Condition, y = ROIsPerTrial)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(group=Ani_Spot), colour="#b5b5b5")+
  ylim(0, 50)+
  geom_point(aes(x=Condition, y=ROIsPerTrialMean), size = 5, colour="#7b3294")+
  geom_line(aes(x=Condition, y=ROIsPerTrialMean,group=Ani_Spot), size=1.5, colour="#7b3294")+
  geom_errorbar(aes(x=Condition,ymin=ROIsPerTrialMean-se, ymax=ROIsPerTrialMean+se), colour="#7b3294", width=0.2,  size=1.5,position=position_dodge(.9)) +
    max.theme

ggplot(ROInum.lck.8strial[ROInum.lck.8strial$Channel=="GCaMP",], aes(x=Condition, y = ROIsPerTrial)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(group=Ani_Spot), colour="#b5b5b5")+
  ylim(0, 50)+
  geom_point(aes(x=Condition, y=ROIsPerTrialMean), size = 5, colour="#008837")+
  geom_line(aes(x=Condition, y=ROIsPerTrialMean,group=Ani_Spot), size=1.5, colour="#008837")+
  geom_errorbar(aes(x=Condition,ymin=ROIsPerTrialMean-se, ymax=ROIsPerTrialMean+se), colour="#008837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  max.theme

ggplot(df.lck.ROInum.8strial, aes(x=Channel,y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  max.theme

Condition_Channel2= interaction(ROInum.lck.8strial$Condition,ROInum.lck.8strial$Channel)
nROI.lck.stim.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.lck.8strial,REML=FALSE)
nROI.lck.stim.model1 = lmer(ROIsPerTrial~ Channel + (1|Animal), ROInum.lck.8strial,REML=FALSE)
nROI.lck.stim.model2 = lmer(ROIsPerTrial ~ Condition + (1|Animal), ROInum.lck.8strial,REML=FALSE)
nROI.lck.stim.model3 = lmer(ROIsPerTrial ~ Condition_Channel2 + (1|Animal), ROInum.lck.8strial,REML=FALSE)
nROI.lck.stim.anova <- anova(nROI.lck.stim.null, nROI.lck.stim.model1,nROI.lck.stim.model2,nROI.lck.stim.model3)
print(nROI.lck.stim.anova)

nROI.lck.stim.Cond_Channel<- glht(nROI.lck.stim.model3, mcp(Condition_Channel2= "Tukey"))
summary(nROI.lck.stim.Cond_Channel)

################
# number of ROIs in each trial for each field of view (across the whole trial)
cyto.OT.8strial<-all.cyto.OT[all.cyto.OT$OnsetTime<8,]

cyto.OT.8strial$Channel <- factor(cyto.OT.8strial$Channel, levels = c("RCaMP","GCaMP"))

ROInum.cyto.8strial<-ddply(cyto.OT.8strial, c("Animal","Spot","Condition","Channel"), summarise, nROIs=length(OnsetTime))

ROInum.cyto.8strial<-ROInum.cyto.8strial[complete.cases(ROInum.cyto.8strial),]
ROInum.cyto.8strial$Ani_Spot<- paste(ROInum.cyto.8strial$Animal, ROInum.cyto.8strial$Spot, sep="_")

# add in number of trials
ROInum.cyto.8strial$Ani_Spot_Cond<-paste(ROInum.cyto.8strial$Animal, ROInum.cyto.8strial$Spot, ROInum.cyto.8strial$Condition, sep="_")
ROInum.cyto.8strial<-merge(ROInum.cyto.8strial, Spot.cyto.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)
ROInum.cyto.8strial$ROIsPerTrial<-ROInum.cyto.8strial$nROIs/ROInum.cyto.8strial$nTrials

# mean
df.cyto.ROInum.8strial<-summarySE(ROInum.cyto.8strial, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition"))

# plots
#boxplot
ggplot(ROInum.cyto.8strial, aes(x=Channel,y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("cyto- ROIs per FOV") +
  max.theme

# paired line plots
df.cyto.ROInum.8strial$ROIsPerTrialMean<-df.cyto.ROInum.8strial$ROIsPerTrial
df.cyto.ROInum.8strial$Chan_Cond<-interaction(df.cyto.ROInum.8strial$Channel, df.cyto.ROInum.8strial$Condition)
ROInum.cyto.8strial$Chan_Cond<-interaction(ROInum.cyto.8strial$Channel, ROInum.cyto.8strial$Condition)

ROInum.cyto.8strial<-merge(ROInum.cyto.8strial, df.cyto.ROInum.8strial[, c("Chan_Cond", "ROIsPerTrialMean","se")], by="Chan_Cond", all.x=TRUE)

ggplot(ROInum.cyto.8strial[ROInum.cyto.8strial$Channel=="RCaMP",], aes(x=Condition, y = ROIsPerTrial)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(group=Ani_Spot), colour="#b5b5b5")+
  ylim(0, 25)+
  geom_point(aes(x=Condition, y=ROIsPerTrialMean), size = 5, colour="#7b3294")+
  geom_line(aes(x=Condition, y=ROIsPerTrialMean,group=Ani_Spot), size=1.5, colour="#7b3294")+
  geom_errorbar(aes(x=Condition,ymin=ROIsPerTrialMean-se, ymax=ROIsPerTrialMean+se), colour="#7b3294", width=0.2,  size=1.5,position=position_dodge(.9)) +
  max.theme

ggplot(ROInum.cyto.8strial[ROInum.cyto.8strial$Channel=="GCaMP",], aes(x=Condition, y = ROIsPerTrial)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(group=Ani_Spot), colour="#b5b5b5")+
  ylim(0, 25)+
  geom_point(aes(x=Condition, y=ROIsPerTrialMean), size = 5, colour="#008837")+
  geom_line(aes(x=Condition, y=ROIsPerTrialMean,group=Ani_Spot), size=1.5, colour="#008837")+
  geom_errorbar(aes(x=Condition,ymin=ROIsPerTrialMean-se, ymax=ROIsPerTrialMean+se), colour="#008837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  max.theme

ggplot(df.cyto.ROInum.8strial, aes(x=Channel,y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("cyto num ROIs/trial per field of view during 8s stim") +
  max.theme

Condition_Channel2= interaction(ROInum.cyto.8strial$Condition,ROInum.cyto.8strial$Channel)
nROI.cyto.stim.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.cyto.8strial,REML=FALSE)
nROI.cyto.stim.model1 = lmer(ROIsPerTrial~ Channel + (1|Animal), ROInum.cyto.8strial,REML=FALSE)
nROI.cyto.stim.model2 = lmer(ROIsPerTrial ~ Condition + (1|Animal), ROInum.cyto.8strial,REML=FALSE)
nROI.cyto.stim.model3 = lmer(ROIsPerTrial ~ Condition_Channel2 + (1|Animal), ROInum.cyto.8strial,REML=FALSE)
nROI.cyto.stim.anova <- anova(nROI.cyto.stim.null, nROI.cyto.stim.model1,nROI.cyto.stim.model2,nROI.cyto.stim.model3)
print(nROI.cyto.stim.anova)

nROI.cyto.stim.Cond_Channel<- glht(nROI.cyto.stim.model3, mcp(Condition_Channel2= "Tukey"))
summary(nROI.cyto.stim.Cond_Channel)


######
 # Lck
# considering 1 s before and 1 s after stimulation 

lck.eitherside.stim2<-subset(stim.lck.OT.dist.2s, OnsetTimeAdjust<1 & OnsetTimeAdjust>-1)
lck.eitherside.stim2$either1s="other"
lck.eitherside.stim2$either1s[lck.eitherside.stim2$OnsetTimeAdjust<0]="before"
lck.eitherside.stim2$either1s[lck.eitherside.stim2$OnsetTimeAdjust>0]="after"

###
# ROI number before and after

# greater number ROIs with fast onset than Nostim?
ROInum.eitherside.1s<-ddply(lck.eitherside.stim2, c("Animal","Spot","Condition","Channel","either1s"), summarise, nROIs=length(OnsetTime))

# add in number of trials
ROInum.eitherside.1s$Ani_Spot_Cond<-paste(ROInum.eitherside.1s$Animal, ROInum.eitherside.1s$Spot, ROInum.eitherside.1s$Condition, sep="_")

ROInum.eitherside.1s<-merge(ROInum.eitherside.1s, Spot.lck.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)

#ROInum.eitherside.1s$nTrials=Spot.lck.ntrials$nTrials[match(ROInum.eitherside.1s$Ani_Spot_Cond,Spot.lck.ntrials$Ani_Spot_Cond)]

ROInum.eitherside.1s$ROIsPerTrial<-ROInum.eitherside.1s$nROIs/ROInum.eitherside.1s$nTrials

ROInum.eitherside.1s$either1s<-factor(ROInum.eitherside.1s$either1s, levels=c("before","after"))

# only consider GCaMP
ROInum.eitherside.1s.GC<-subset(ROInum.eitherside.1s, Channel=="GCaMP")

df.ROInum.either1s.GC<-summarySE(ROInum.eitherside.1s.GC, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition","either1s"))
df.ROInum.either1s.GC2<-summarySE(ROInum.eitherside.1s.GC[ROInum.eitherside.1s.GC$Condition=="Stim",], measurevar = "ROIsPerTrial", groupvars = c("either1s"))

ggplot(df.ROInum.either1s.GC, aes(x=Condition,y=ROIsPerTrial, fill= either1s)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ggtitle("num ROIs per field of view- 1 s before vs 1 s after") +
  max.theme

ggplot(df.ROInum.either1s.GC2, aes(x=either1s,y=ROIsPerTrial, fill= either1s)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ggtitle("num ROIs per field of view- 1 s before vs 1 s after") +
  max.theme

ROInum.eitherside.1s.GC$either1s<-as.factor(ROInum.eitherside.1s.GC$either1s)
Condition_either1s= interaction(ROInum.eitherside.1s.GC$Condition,ROInum.eitherside.1s.GC$either1s)

nROI.either1s.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.eitherside.1s.GC,REML=FALSE)
nROI.either1s.model1 = lmer(ROIsPerTrial ~ either1s + (1|Animal), ROInum.eitherside.1s.GC,REML=FALSE)
nROI.either1s.model2 = lmer(ROIsPerTrial ~ Condition + (1|Animal), ROInum.eitherside.1s.GC,REML=FALSE)
nROI.either1s.model3 = lmer(ROIsPerTrial ~ Condition_either1s + (1|Animal), ROInum.eitherside.1s.GC,REML=FALSE)
nROI.either1s.anova <- anova(nROI.either1s.null,nROI.either1s.model1,nROI.either1s.model2,nROI.either1s.model3)
print(nROI.either1s.anova)

nROI.either1s.Cond_Gr<- glht(nROI.either1s.model3, mcp(Condition_either1s= "Tukey"))
summary(nROI.either1s.Cond_Gr)

nROI.either1s.pv<- glht(nROI.either1s.model1, mcp(either1s= "Tukey"))
summary(nROI.either1s.pv)

# only consider RCaMP
ROInum.eitherside.1s.RC<-subset(ROInum.eitherside.1s, Channel=="RCaMP")

df.ROInum.either1s.RC<-summarySE(ROInum.eitherside.1s.RC, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition","either1s"))

ggplot(df.ROInum.either1s.RC, aes(x=Condition,y=ROIsPerTrial, fill= either1s)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ggtitle("num ROIs per field of view- 1 s before vs 1 s after") +
  max.theme

ROInum.eitherside.1s.RC$either1s<-as.factor(ROInum.eitherside.1s.RC$either1s)
Condition_either1s= interaction(ROInum.eitherside.1s.RC$Condition,ROInum.eitherside.1s.RC$either1s)

nROI.either1s.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.eitherside.1s.RC,REML=FALSE)
nROI.either1s.model1 = lmer(ROIsPerTrial ~ either1s + (1|Animal), ROInum.eitherside.1s.RC,REML=FALSE)
nROI.either1s.model2 = lmer(ROIsPerTrial ~ Condition + (1|Animal), ROInum.eitherside.1s.RC,REML=FALSE)
nROI.either1s.model3 = lmer(ROIsPerTrial ~ Condition_either1s + (1|Animal), ROInum.eitherside.1s.RC,REML=FALSE)
nROI.either1s.anova <- anova(nROI.either1s.null,nROI.either1s.model1,nROI.either1s.model2,nROI.either1s.model3)
print(nROI.either1s.anova)

nROI.either1s.Cond_Gr<- glht(nROI.either1s.model3, mcp(Condition_either1s= "Tukey"))
summary(nROI.either1s.Cond_Gr)






 ##############
# cyto
# considering 1 s before and 1 s after stimulation 

cyto.eitherside.stim2<-subset(stim.cyto.OT.dist.2s, OnsetTimeAdjust<1 & OnsetTimeAdjust>-1)
cyto.eitherside.stim2$either1s="other"
cyto.eitherside.stim2$either1s[cyto.eitherside.stim2$OnsetTimeAdjust<0]="before"
cyto.eitherside.stim2$either1s[cyto.eitherside.stim2$OnsetTimeAdjust>0]="after"

####
# ROI number before and after

# greater number ROIs with fast onset than Nostim?
cyto.ROInum.eitherside.1s<-ddply(cyto.eitherside.stim2, c("Animal","Spot","Condition","Channel","either1s"), summarise, nROIs=length(OnsetTime))

# add in number of trials
cyto.ROInum.eitherside.1s$Ani_Spot_Cond<-paste(cyto.ROInum.eitherside.1s$Animal, cyto.ROInum.eitherside.1s$Spot, cyto.ROInum.eitherside.1s$Condition, sep="_")

cyto.ROInum.eitherside.1s<-merge(cyto.ROInum.eitherside.1s, Spot.cyto.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)

#ROInum.eitherside.1s$nTrials=Spot.cyto.ntrials$nTrials[match(ROInum.eitherside.1s$Ani_Spot_Cond,Spot.cyto.ntrials$Ani_Spot_Cond)]

cyto.ROInum.eitherside.1s$ROIsPerTrial<-cyto.ROInum.eitherside.1s$nROIs/cyto.ROInum.eitherside.1s$nTrials

cyto.ROInum.eitherside.1s$either1s<-factor(cyto.ROInum.eitherside.1s$either1s, levels=c("before","after"))

# only consider GCaMP
cyto.ROInum.eitherside.1s.GC<-subset(cyto.ROInum.eitherside.1s, Channel=="GCaMP")

df.cyto.ROInum.either1s.GC<-summarySE(cyto.ROInum.eitherside.1s.GC, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition","either1s"))
df.cyto.ROInum.either1s.GC2<-summarySE(cyto.ROInum.eitherside.1s.GC[cyto.ROInum.eitherside.1s.GC$Condition=="Stim",], measurevar = "ROIsPerTrial", groupvars = c("either1s"))

ggplot(df.cyto.ROInum.either1s.GC, aes(x=Condition,y=ROIsPerTrial, fill= either1s)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ggtitle("num ROIs per field of view- 1 s before vs 1 s after") +
  max.theme

ggplot(df.cyto.ROInum.either1s.GC2, aes(x=either1s,y=ROIsPerTrial, fill= either1s)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ggtitle("num ROIs per field of view- 1 s before vs 1 s after") +
  max.theme

cyto.ROInum.eitherside.1s.GC$either1s<-as.factor(cyto.ROInum.eitherside.1s.GC$either1s)
Condition_either1s= interaction(cyto.ROInum.eitherside.1s.GC$Condition,cyto.ROInum.eitherside.1s.GC$either1s)

cyto.nROI.either1s.null = lmer(ROIsPerTrial ~ (1|Animal), cyto.ROInum.eitherside.1s.GC,REML=FALSE)
cyto.nROI.either1s.model1 = lmer(ROIsPerTrial ~ either1s + (1|Animal), cyto.ROInum.eitherside.1s.GC,REML=FALSE)
cyto.nROI.either1s.model2 = lmer(ROIsPerTrial ~ Condition + (1|Animal), cyto.ROInum.eitherside.1s.GC,REML=FALSE)
cyto.nROI.either1s.model3 = lmer(ROIsPerTrial ~ Condition_either1s + (1|Animal), cyto.ROInum.eitherside.1s.GC,REML=FALSE)
cyto.nROI.either1s.anova <- anova(cyto.nROI.either1s.null,cyto.nROI.either1s.model1,cyto.nROI.either1s.model2,cyto.nROI.either1s.model3)
print(cyto.nROI.either1s.anova)

cyto.nROI.either1s.Cond_Gr<- glht(cyto.nROI.either1s.model3, mcp(Condition_either1s= "Tukey"))
summary(cyto.nROI.either1s.Cond_Gr)

# only consider RCaMP
cyto.ROInum.eitherside.1s.RC<-subset(cyto.ROInum.eitherside.1s, Channel=="RCaMP")

df.cyto.ROInum.either1s.RC<-summarySE(cyto.ROInum.eitherside.1s.RC, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition","either1s"))

ggplot(df.cyto.ROInum.either1s.RC, aes(x=Condition,y=ROIsPerTrial, fill= either1s)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ggtitle("num ROIs per field of view- 1 s before vs 1 s after") +
  max.theme

cyto.ROInum.eitherside.1s.RC$either1s<-as.factor(cyto.ROInum.eitherside.1s.RC$either1s)
Condition_either1s= interaction(cyto.ROInum.eitherside.1s.RC$Condition,cyto.ROInum.eitherside.1s.RC$either1s)

cyto.nROI.either1s.null = lmer(ROIsPerTrial ~ (1|Animal), cyto.ROInum.eitherside.1s.RC,REML=FALSE)
cyto.nROI.either1s.model1 = lmer(ROIsPerTrial ~ either1s + (1|Animal), cyto.ROInum.eitherside.1s.RC,REML=FALSE)
cyto.nROI.either1s.model2 = lmer(ROIsPerTrial ~ Condition + (1|Animal), cyto.ROInum.eitherside.1s.RC,REML=FALSE)
cyto.nROI.either1s.model3 = lmer(ROIsPerTrial ~ Condition_either1s + (1|Animal), cyto.ROInum.eitherside.1s.RC,REML=FALSE)
cyto.nROI.either1s.anova <- anova(cyto.nROI.either1s.null,cyto.nROI.either1s.model1,cyto.nROI.either1s.model2,cyto.nROI.either1s.model3)
print(cyto.nROI.either1s.anova)

cyto.nROI.either1s.Cond_Gr<- glht(cyto.nROI.either1s.model3, mcp(Condition_either1s= "Tukey"))
summary(cyto.nROI.either1s.Cond_Gr)


###########################################

#################

# for median and mean calculations

# use thresholds calculated from the 2.33 SD of the log transformed astrocyte data


LongN_PTwind=8  # 8 seconds
LongN_OTwind=8  # 8 seconds

# test for normal distribution of astrocyte onset times
astro.LCk.stim<-subset(all.lck.OT, Channel=="GCaMP" & Condition=="Stim")
shapiro.test(astro.LCk.stim$OnsetTime)

# remove data that is outside the above windows
# lck
stim.lck.OT.window<-subset(all.lck.OT, OnsetTime<=threshold.Lck.OT)

stim.lck.OT.R<-subset(all.lck.OT, Channel=="RCaMP" & OnsetTime<=LongN_OTwind)
stim.lck.OT.G<-subset(all.lck.OT, Channel=="GCaMP" & OnsetTime<=threshold.Lck.OT)

# neurons have a different window that astrocytes
stim.lck.OT.window2<-rbind(stim.lck.OT.R, stim.lck.OT.G)

# test for normal distribution of astrocyte onset times
shapiro.test(stim.lck.OT.G$OnsetTime[stim.lck.OT.G$Condition=="Stim"])


#cyto
stim.cyto.OT.window<-subset(all.cyto.OT,OnsetTime<=threshold.cyto.OT)

stim.cyto.OT.R<-subset(all.cyto.OT, Channel=="RCaMP" & OnsetTime<=LongN_OTwind)
stim.cyto.OT.G<-subset(all.cyto.OT, Channel=="GCaMP" & OnsetTime<=threshold.cyto.OT)

stim.cyto.OT.window2<-rbind(stim.cyto.OT.R, stim.cyto.OT.G)



# peak times
stim.lck.PT.window<-subset(all.lck.peaks, peakTime<=threshold.Lck.PT)

stim.lck.PT.R<-subset(all.lck.peaks, Channel=="RCaMP" & peakTime<=LongN_PTwind)
stim.lck.PT.G<-subset(all.lck.peaks, Channel=="GCaMP" & peakTime<=threshold.Lck.PT)

# neurons have a different window that astrocytes
stim.lck.PT.window2<-rbind(stim.lck.PT.R, stim.lck.PT.G)


#cyto
stim.cyto.PT.window<-subset(all.cyto.peaks,peakTime<=threshold.cyto.PT)

stim.cyto.PT.R<-subset(all.cyto.peaks, Channel=="RCaMP" & peakTime<=LongN_PTwind)
stim.cyto.PT.G<-subset(all.cyto.peaks, Channel=="GCaMP" & peakTime<=threshold.cyto.PT)

stim.cyto.PT.window2<-rbind(stim.cyto.PT.R, stim.cyto.PT.G)

rm(stim.cyto.OT.G, stim.cyto.OT.R, stim.lck.OT.G, stim.lck.OT.R, stim.cyto.PT.G, stim.cyto.PT.R, stim.lck.PT.G, stim.lck.PT.R)


# no need to combine lck and cyto data tables, as we never make direct comparisons between the two sensors

# onset times- based on neurons and astrocytes from the same astrocyte time window (see threshold above)

stim.lck.OT.window$Channel=factor(stim.lck.OT.window$Channel, levels = c("RCaMP","GCaMP"))
stim.cyto.OT.window$Channel=factor(stim.cyto.OT.window$Channel, levels = c("RCaMP","GCaMP"))


df.lck.OT.mean<-summarySE(stim.lck.OT.window, measurevar = "OnsetTime", groupvars = c("Channel", "Condition"))
df.lck.OT.mean.stim<-summarySE(stim.lck.OT.window[stim.lck.OT.window$Condition=="Stim",], measurevar = "OnsetTime", groupvars = c("Channel", "Condition"))
df.lck.OT.mean.stim2<-summarySE(stim.lck.alldata[stim.lck.alldata$Condition=="Stim",], measurevar = "OnsetTime", groupvars = c("Channel", "Condition"))

df.cyto.OT.mean<-summarySE(stim.cyto.OT.window, measurevar = "OnsetTime", groupvars = c("Channel", "Condition"))
df.cyto.OT.mean.stim<-summarySE(stim.cyto.OT.window[stim.cyto.OT.window$Condition=="Stim",], measurevar = "OnsetTime", groupvars = c("Channel", "Condition"))


ggplot(df.lck.OT.mean.stim, aes(x=Channel,y=OnsetTime, fill= Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("lck Mean Onset Time (s)") +
  max.theme

ggplot(df.cyto.OT.mean.stim, aes(x=Channel,y=OnsetTime, fill= Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("cyto Mean Onset Time (s)") +
  max.theme


# means stats Lck
#only stim
OT.lck.onlystim.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial_Cond), stim.lck.OT.window[stim.lck.OT.window$Condition=="Stim",],REML=FALSE)
OT.lck.onlystim.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Spot_trial_Cond), stim.lck.OT.window[stim.lck.OT.window$Condition=="Stim",],REML=FALSE)
OT.lck.onlystim.anova <- anova(OT.lck.onlystim.null, OT.lck.onlystim.model1)
print(OT.lck.onlystim.anova)

OT.lck.onlystim.Channel<- glht(OT.lck.onlystim.model1, mcp(Channel= "Tukey"))
summary(OT.lck.onlystim.Channel)

# means stats cyto
#only stim
OT.cyto.onlystim.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial), stim.cyto.OT.window[stim.cyto.OT.window$Condition=="Stim",],REML=FALSE)
OT.cyto.onlystim.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.cyto.OT.window[stim.cyto.OT.window$Condition=="Stim",],REML=FALSE)
OT.cyto.onlystim.anova <- anova(OT.cyto.onlystim.null, OT.cyto.onlystim.model1)
print(OT.cyto.onlystim.anova)

OT.cyto.onlystim.Channel<- glht(OT.cyto.onlystim.model1, mcp(Channel= "Tukey"))
summary(OT.cyto.onlystim.Channel)


####
# peak times

stim.lck.PT.window$Channel=factor(stim.lck.PT.window$Channel, levels = c("RCaMP","GCaMP"))
stim.cyto.PT.window$Channel=factor(stim.cyto.PT.window$Channel, levels = c("RCaMP","GCaMP"))


df.lck.PT.mean<-summarySE(stim.lck.PT.window, measurevar = "peakTime", groupvars = c("Channel", "Condition"))
df.lck.PT.mean.stim<-summarySE(stim.lck.PT.window[stim.lck.PT.window$Condition=="Stim",], measurevar = "peakTime", groupvars = c("Channel", "Condition"))

df.cyto.PT.mean<-summarySE(stim.cyto.PT.window, measurevar = "peakTime", groupvars = c("Channel", "Condition"))
df.cyto.PT.mean.stim<-summarySE(stim.cyto.PT.window[stim.cyto.PT.window$Condition=="Stim",], measurevar = "peakTime", groupvars = c("Channel", "Condition"))


ggplot(df.lck.PT.mean.stim, aes(x=Channel,y=peakTime, fill= Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("lck Mean Peak Time (s)") +
  max.theme

ggplot(df.cyto.PT.mean.stim, aes(x=Channel,y=peakTime, fill= Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("cyto Mean Peak Time (s)") +
  ylim(0,10)+
  max.theme


# means stats Lck
#only stim
PT.lck.onlystim.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.PT.window[stim.lck.PT.window$Condition=="Stim",],REML=FALSE)
PT.lck.onlystim.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.PT.window[stim.lck.PT.window$Condition=="Stim",],REML=FALSE)
PT.lck.onlystim.anova <- anova(PT.lck.onlystim.null, PT.lck.onlystim.model1)
print(PT.lck.onlystim.anova)

PT.lck.onlystim.Channel<- glht(PT.lck.onlystim.model1, mcp(Channel= "Tukey"))
summary(PT.lck.onlystim.Channel)

# means stats cyto
#only stim
PT.cyto.onlystim.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.cyto.PT.window[stim.cyto.PT.window$Condition=="Stim",],REML=FALSE)
PT.cyto.onlystim.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.cyto.PT.window[stim.cyto.PT.window$Condition=="Stim",],REML=FALSE)
PT.cyto.onlystim.anova <- anova(PT.cyto.onlystim.null, PT.cyto.onlystim.model1)
print(PT.cyto.onlystim.anova)

PT.cyto.onlystim.Channel<- glht(PT.cyto.onlystim.model1, mcp(Channel= "Tukey"))
summary(PT.cyto.onlystim.Channel)


######
# combine peak times and onset time tables so we only consider peaks with an onset time

stim.lck.alldata<-merge(stim.lck.PT.window, stim.lck.OT.window[, c("ROIs_Cond", "OnsetTime","TraceAUC1","TraceAUC10")], by="ROIs_Cond")
stim.cyto.alldata<-merge(stim.cyto.PT.window, stim.cyto.OT.window[, c("ROIs_Cond", "OnsetTime","TraceAUC1","TraceAUC10")], by="ROIs_Cond")


# identify "FAST" astrocytes
stim.lck.alldata$Group<-0
stim.lck.alldata$Group[stim.lck.alldata$OnsetTime<=NeuronOT50]<-"fast"
stim.lck.alldata$Group[stim.lck.alldata$OnsetTime>NeuronOT50]<-"delayed"


#####
# number of fast ROIs in each trial for each field of view 
# based on neuronal percentile onset time

stim.lck.alldata$Channel <- factor(stim.lck.alldata$Channel, levels = c("RCaMP","GCaMP"))
lck.stim<-subset(stim.lck.alldata, Channel=="GCaMP" & Condition=="Stim")
lck.stim$Ani_Spot_Cond<-paste(lck.stim$Animal, lck.stim$Spot, lck.stim$Condition, sep="_")

ROInum.lck.group<-ddply(lck.stim, c("Animal","Spot","Condition","Group","Channel","Animal_Spot", "Ani_Spot_Cond"), summarise, nROIs=length(OnsetTime))

# add in number of trials
ROInum.lck.group<-merge(ROInum.lck.group, Spot.lck.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)

ROInum.lck.group$ROIsPerTrial<-ROInum.lck.group$nROIs/ROInum.lck.group$nTrials


# mean
df.lck.ROInum.group.mean<-summarySE(ROInum.lck.group, measurevar = "ROIsPerTrial", groupvars = c("Group"))


ggplot(df.lck.ROInum.group.mean, aes(x=Group,y=ROIsPerTrial, fill= Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view") +
  max.theme

######
# fast vs delayed

stim.lck.alldata$Channel_Group<-interaction(stim.lck.alldata$Channel, stim.lck.alldata$Group)
stim.lck.alldata$Channel_Group<-as.factor(stim.lck.alldata$Channel_Group)
stim.lck.alldata$Channel_Group<-factor(stim.lck.alldata$Channel_Group, levels = c("RCaMP.fast","GCaMP.fast","GCaMP.delayed"))

stim.lck.compdata<-stim.lck.alldata[!(stim.lck.alldata$Channel=="RCaMP"& stim.lck.alldata$Group=="delayed"),]

# take out the effect of Condition
# we are only interested in stim case
stim.lck.compdata.STIM<-subset(stim.lck.compdata, Condition=="Stim")

# mean onset times
df.OT2<- summarySE(stim.lck.alldata, measurevar = "OnsetTime", groupvars = c("Channel_Group","Condition"))
df.OT3<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim",], measurevar = "OnsetTime", groupvars = c("Channel_Group"))
df.OT4<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim"& stim.lck.compdata$Channel=="GCaMP",], measurevar = "OnsetTime", groupvars = c("ROIType","Channel_Group"))
df.OT5<- summarySE(stim.lck.compdata.STIM, measurevar = "OnsetTime", groupvars = c("ROIType","Channel_Group"))

median.Lck<-group_by(stim.lck.compdata, Channel_Group, Condition)
summarise(median.Lck, Median=median(OnsetTime))

ggplot(df.OT2, aes(x=Channel_Group,y=OnsetTime, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.OT3, aes(x=Channel_Group,y=OnsetTime, fill=Channel_Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.OT4, aes(x=Channel_Group,y=OnsetTime, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme


ggplot(stim.lck.alldata, aes(x=Channel_Group,y=OnsetTime, fill=Channel_Group)) +
  geom_boxplot(notch=TRUE)+
  ylab("Onset Time (s)") +
  ggtitle("lck data- fast vs delayed")+
  max.theme

ggplot(stim.lck.compdata.STIM, aes(x=Channel_Group,y=OnsetTime, fill=Channel_Group)) +
  geom_boxplot()+
  ylab("Onset Time (s)") +
  ggtitle("lck data- fast vs delayed")+
  max.theme

# stats
Group_Channel_Cond_Type=interaction(stim.lck.alldata$Group,stim.lck.alldata$Channel,stim.lck.alldata$Condition,stim.lck.alldata$ROIType)
Group_Channel_Cond=interaction(stim.lck.alldata$Channel_Group,stim.lck.alldata$Condition)

# stats for onset times- neurons vs astrocytes
OT.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model2 = lmer(OnsetTime ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.alldata,REML=FALSE)
OT.model3 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model4 = lmer(OnsetTime ~ Channel_Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model5 = lmer(OnsetTime ~ Group_Channel_Cond + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model6 = lmer(OnsetTime ~ Group_Channel_Cond_Type + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.anova <- anova(OT.null, OT.model1,OT.model2,OT.model3,OT.model4,OT.model5,OT.model6)
print(OT.anova)

OT.Group_channel<- glht(OT.model5, mcp(Group_Channel_Cond= "Tukey"))
summary(OT.Group_channel)

#OT.Group_channel_ty<- glht(OT.model6, mcp(Group_Channel_Cond_Type= "Tukey"))
#summary(OT.Group_channel_ty)


Group_Channel_Type=interaction(stim.lck.compdata.STIM$Group,stim.lck.compdata.STIM$Channel,stim.lck.compdata.STIM$ROIType)
Group_Channel=interaction(stim.lck.compdata.STIM$Channel, stim.lck.compdata.STIM$Group)

# stats for onset times- neurons vs astrocytes
OT.stim.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.model3 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.model4 = lmer(OnsetTime ~ Group_Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.model6 = lmer(OnsetTime ~ Group_Channel_Type + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.anova <- anova(OT.stim.null, OT.stim.model1,OT.stim.model3,OT.stim.model4,OT.stim.model6)
print(OT.stim.anova)

OT.stim.Group_channel<- glht(OT.stim.model4, mcp(Group_Channel= "Tukey"))
summary(OT.stim.Group_channel)

OT.stim.Group_channel_type<- glht(OT.stim.model6, mcp(Group_Channel_Type= "Tukey"))
summary(OT.stim.Group_channel_type)

summary(OT.stim.model4)

# check residuals for linearity
plot(fitted(OT.stim.model4), residuals(OT.stim.model4),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(OT.stim.model4), residuals(OT.stim.model4)), col=46, lwd=2.5)


#median test 

median.test <- function(x, y){
  z <- c(x, y)
  g <- rep(1:2, c(length(x), length(y)))
  m <- median(z)
  fisher.test(z < m, g)$p.value
}

median.test(stim.lck.compdata.STIM$OnsetTime[stim.lck.compdata.STIM$Channel_Group=="RCaMP.fast"], 
            stim.lck.compdata.STIM$OnsetTime[stim.lck.compdata.STIM$Channel_Group=="GCaMP.fast"])

wilcox.test(stim.lck.compdata.STIM$OnsetTime[stim.lck.compdata.STIM$Channel_Group=="RCaMP.fast"], 
            stim.lck.compdata.STIM$OnsetTime[stim.lck.compdata.STIM$Channel_Group=="GCaMP.fast"])



######
# peak time
df.pT1<-summarySE(stim.cyto.alldata, measurevar = "peakTime", groupvars = c("Channel", "Condition"))
df.pT2<-summarySE(stim.both.alldata, measurevar = "peakTime", groupvars = c("Channel", "Group","Condition"))
df.pT3<- summarySE(stim.lck.alldata, measurevar = "peakTime", groupvars = c("Channel_Group","Condition"))
df.pT4<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim",], measurevar = "peakTime", groupvars = c("Channel_Group"))
df.pT5<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim"& stim.lck.compdata$Channel=="GCaMP",], measurevar = "peakTime", groupvars = c("ROIType","Channel_Group"))
df.pT6<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim"& stim.lck.compdata$Channel=="GCaMP",], measurevar = "peakTime", groupvars = c("Channel_Group"))

ggplot(df.pT1, aes(x=Channel,y=peakTime, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakTime") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.pT3, aes(x=Channel_Group,y=peakTime, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakTime") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.pT4, aes(x=Channel_Group,y=peakTime, fill= Channel_Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakTime") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.pT5, aes(x=Channel_Group,y=peakTime, fill= ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakTime") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.pT6, aes(x=Channel_Group,y=peakTime, fill= Channel_Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakTime") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(stim.lck.compdata, aes(x=Channel_Group,y=peakTime, fill=Channel_Group)) +
  geom_boxplot()+
  ylab("peakTime") +
  ggtitle("lck data- fast vs delayed")+
  max.theme

ggplot(stim.lck.compdata, aes(x=peakTime, y=..density..,fill=Channel_Group)) +
  geom_histogram(binwidth=0.5, position="dodge")+
  ylab("density") +
  xlim(-2,15)+
  ggtitle("lck data- fast vs delayed")+
  max.theme

ggplot(stim.lck.alldata[stim.lck.alldata$Channel_Group=="lck_GCaMP.fast",], aes(x=peakTime, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.5, position="dodge")+
  ylab("density") +
  xlim(-2,15)+
  ggtitle("lck data- fast no stim vs stim")+
  max.theme

# stats for duration- neurons vs astrocytes
pT.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
pT.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
pT.model2 = lmer(peakTime ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.alldata,REML=FALSE)
pT.model3 = lmer(peakTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
pT.model4 = lmer(peakTime ~ Channel_Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
pT.model5 = lmer(peakTime ~ Group_Channel_Cond + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
pT.model6 = lmer(peakTime ~ Group_Channel_Cond_Type + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
pT.anova <- anova(pT.null, pT.model1,pT.model2,pT.model3,pT.model4,pT.model5,pT.model6)
print(pT.anova)

pT.Group_channel<- glht(pT.model5, mcp(Group_Channel_Cond= "Tukey"))
summary(pT.Group_channel)

pT.Group_channel_ty<- glht(pT.model6, mcp(Group_Channel_Cond_Type= "Tukey"))
summary(pT.Group_channel_ty)



## only for STIM case

# stats for onset times- neurons vs astrocytes
pT.stim.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
pT.stim.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
pT.stim.model3 = lmer(peakTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
pT.stim.model4 = lmer(peakTime ~ Group_Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
pT.stim.model6 = lmer(peakTime ~ Group_Channel_Type + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
pT.stim.anova <- anova(pT.stim.null, pT.stim.model1,pT.stim.model3,pT.stim.model4,pT.stim.model6)
print(pT.stim.anova)

pT.stim.Group_channel<- glht(pT.stim.model4, mcp(Group_Channel= "Tukey"))
summary(pT.stim.Group_channel)

pT.stim.Group_channel_type<- glht(pT.stim.model6, mcp(Group_Channel_Type= "Tukey"))
summary(pT.stim.Group_channel_type)

########
#amplitude

#aggregate Lck & cyto data, so we have the mean per FOV

lck.GC.spot<-ddply(stim.lck.alldata[stim.lck.alldata$Channel=="GCaMP",], c("Animal","Spot","Animal_Spot","Condition"), summarise, 
                meanAmp=mean(amplitude), meanDur=mean(Duration))
lck.GC.spot$Ani_Cond<-paste(lck.GC.spot$Animal_Spot, lck.GC.spot$Condition, sep="_")
lck.GC.spot.proc<-ddply(stim.lck.alldata[stim.lck.alldata$Channel=="GCaMP" & stim.lck.alldata$ROIType=="Process",], c("Animal","Spot","Animal_Spot","Condition"), summarise, 
                   meanProcArea=mean(area))
lck.GC.spot.proc$Ani_Cond<-paste(lck.GC.spot.proc$Animal_Spot, lck.GC.spot.proc$Condition, sep="_")

lck.GC.spot<-merge(lck.GC.spot, lck.GC.spot.proc[, c("Ani_Cond", "meanProcArea")], by="Ani_Cond", all.x=TRUE)

#cyto
stim.cyto.alldata$Animal_Spot<- paste(stim.cyto.alldata$Animal, stim.cyto.alldata$Spot, sep="_")

cyto.GC.spot<-ddply(stim.cyto.alldata[stim.cyto.alldata$Channel=="GCaMP",], c("Animal","Spot","Animal_Spot","Condition"), summarise, 
                   meanAmp=mean(amplitude), meanDur=mean(Duration))
cyto.GC.spot$Ani_Cond<-paste(cyto.GC.spot$Animal_Spot, cyto.GC.spot$Condition, sep="_")
cyto.GC.spot.proc<-ddply(stim.cyto.alldata[stim.cyto.alldata$Channel=="GCaMP" & stim.cyto.alldata$ROIType=="Process",], c("Animal","Spot","Animal_Spot","Condition"), summarise, 
                        meanProcArea=mean(area))
cyto.GC.spot.proc$Ani_Cond<-paste(cyto.GC.spot.proc$Animal_Spot, cyto.GC.spot.proc$Condition, sep="_")

cyto.GC.spot<-merge(cyto.GC.spot, cyto.GC.spot.proc[, c("Ani_Cond", "meanProcArea")], by="Ani_Cond", all.x=TRUE)


# means for plots
df.amp1A1<-summarySE(stim.cyto.alldata[stim.cyto.alldata$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("Channel", "Condition"))
df.amp1A2<-summarySE(cyto.GC.spot, measurevar = "meanAmp", groupvars = c("Condition"))
# paired line plots
df.amp1A2$Mean_meanAmp<-df.amp1A2$meanAmp
cyto.GC.spot<-merge(cyto.GC.spot, df.amp1A2[, c("Condition", "Mean_meanAmp","se")], by="Condition", all.x=TRUE)

df.amp1B1<-summarySE(stim.lck.alldata[stim.lck.alldata$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("Channel", "Condition"))
df.amp1B2<-summarySE(lck.GC.spot, measurevar = "meanAmp", groupvars = c("Condition"))
# paired line plots
df.amp1B2$Mean_meanAmp<-df.amp1B2$meanAmp
lck.GC.spot<-merge(lck.GC.spot, df.amp1B2[, c("Condition", "Mean_meanAmp","se")], by="Condition", all.x=TRUE)

df.amp3<- summarySE(stim.lck.compdata, measurevar = "amplitude", groupvars = c("Channel_Group","Condition"))
df.amp4<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim",], measurevar = "amplitude", groupvars = c("Channel_Group"))
df.amp5<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim"& stim.lck.compdata$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("ROIType","Channel_Group"))
df.amp6<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="RCaMP",], measurevar = "amplitude", groupvars = c("Channel_Group"))


# plots
ggplot(df.amp1A1, aes(x=Condition,y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("cyto amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp1A2, aes(x=Condition,y=meanAmp, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanAmp-se, ymax=meanAmp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("cyto amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

# paired plot
ggplot(cyto.GC.spot, aes(x=Condition, y = meanAmp)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(group=Animal_Spot), colour="#b5b5b5")+
  geom_point(aes(x=Condition, y=Mean_meanAmp), size = 5, colour="#1b7837")+
  geom_line(aes(x=Condition, y=Mean_meanAmp,group=Animal_Spot), size=1.5, colour="#1b7837")+
  geom_errorbar(aes(x=Condition,ymin=Mean_meanAmp-se, ymax=Mean_meanAmp+se), colour="#1b7837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  ggtitle("cyto amp per FOV")+
  max.theme



ggplot(df.amp1B1, aes(x=Condition,y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("lck amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp1B2, aes(x=Condition,y=meanAmp, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanAmp-se, ymax=meanAmp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("lck amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(lck.GC.spot, aes(x=Condition, y = meanAmp)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(group=Animal_Spot), colour="#b5b5b5")+
  geom_point(aes(x=Condition, y=Mean_meanAmp), size = 5, colour="#1b7837")+
  geom_line(aes(x=Condition, y=Mean_meanAmp,group=Animal_Spot), size=1.5, colour="#1b7837")+
  geom_errorbar(aes(x=Condition,ymin=Mean_meanAmp-se, ymax=Mean_meanAmp+se), colour="#1b7837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  ggtitle("lck amp per FOV")+
max.theme


ggplot(df.amp3, aes(x=Channel_Group,y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp4, aes(x=Channel_Group,y=amplitude, fill= Channel_Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp5, aes(x=Channel_Group,y=amplitude, fill= ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp6, aes(x=Channel_Group,y=amplitude, fill= Channel_Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(stim.lck.compdata, aes(x=Channel_Group,y=amplitude, fill=Channel_Group)) +
  geom_boxplot(notch=TRUE)+
  ylab("amplitude") +
  ggtitle("lck data- fast vs delayed")+
  max.theme

ggplot(stim.lck.compdata, aes(x=amplitude, y=..density..,fill=Channel_Group)) +
  geom_histogram(binwidth=0.5, position="dodge")+
  ylab("density") +
  xlim(-2,15)+
  ggtitle("lck data- fast vs delayed")+
  max.theme

ggplot(stim.lck.alldata[stim.lck.alldata$Channel_Group=="GCaMP.fast",], aes(x=amplitude, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.5, position="dodge")+
  ylab("density") +
  xlim(-2,15)+
  ggtitle("lck data- fast no stim vs stim")+
  max.theme


# supplementary figures

#cyto GC
Cond_Channel=interaction(stim.cyto.alldata$Condition,stim.cyto.alldata$Channel)
amp.cyto.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), stim.cyto.alldata,REML=FALSE)
amp.cyto.model1 = lmer(amplitude ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.cyto.alldata,REML=FALSE)
amp.cyto.model2 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.cyto.alldata,REML=FALSE)
amp.cyto.model3 = lmer(amplitude ~ Cond_Channel + (1|Animal) + (1|Spot) + (1|trials), stim.cyto.alldata,REML=FALSE)

amp.cyto.anova <- anova(amp.cyto.null, amp.cyto.model1,amp.cyto.model2,amp.cyto.model3)
print(amp.cyto.anova)

amp.cyto.Cond_channel<- glht(amp.cyto.model3, mcp(Cond_Channel= "Tukey"))
summary(amp.cyto.Cond_channel)

#cyto GC spots
amp.cyto.spot.null = lmer(meanAmp ~ (1|Animal), cyto.GC.spot,REML=FALSE)
amp.cyto.spot.model1 = lmer(meanAmp ~ Condition + (1|Animal), cyto.GC.spot,REML=FALSE)
amp.cyto.spot.anova <- anova(amp.cyto.spot.null, amp.cyto.spot.model1)
print(amp.cyto.spot.anova)

amp.cyto.spot.pv<- glht(amp.cyto.spot.model1, mcp(Condition= "Tukey"))
summary(amp.cyto.spot.pv)

#lck GC
Cond_Channel=interaction(stim.lck.alldata$Condition,stim.lck.alldata$Channel)
amp.lck.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
amp.lck.model1 = lmer(amplitude ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
amp.lck.model2 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.alldata,REML=FALSE)
amp.lck.model3 = lmer(amplitude ~ Cond_Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)

amp.lck.anova <- anova(amp.lck.null, amp.lck.model1,amp.lck.model2,amp.lck.model3)
print(amp.lck.anova)

amp.lck.Cond_channel<- glht(amp.lck.model3, mcp(Cond_Channel= "Tukey"))
summary(amp.lck.Cond_channel)

#lck GC spots
amp.lck.spot.null = lmer(meanAmp ~ (1|Animal), lck.GC.spot,REML=FALSE)
amp.lck.spot.model1 = lmer(meanAmp ~ Condition + (1|Animal), lck.GC.spot,REML=FALSE)
amp.lck.spot.anova <- anova(amp.lck.spot.null, amp.lck.spot.model1)
print(amp.lck.spot.anova)

amp.lck.spot.pv<- glht(amp.lck.spot.model1, mcp(Condition= "Tukey"))
summary(amp.lck.spot.pv)

## lck data and group (fast, delayed)
comp_Channel_Group=interaction(stim.lck.compdata$Group,stim.lck.compdata$Channel)
comp_Channel_Group_Cond=interaction(stim.lck.compdata$Group,stim.lck.compdata$Channel,stim.lck.compdata$Condition)
comp_Channel_Group_Cond_Type=interaction(stim.lck.compdata$Group,stim.lck.compdata$Channel,
                                         stim.lck.compdata$Condition,stim.lck.compdata$ROIType)


# stats for onset times- neurons vs astrocytes
amp.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata,REML=FALSE)
amp.model1 = lmer(amplitude ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata,REML=FALSE)
amp.model2 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata,REML=FALSE)
amp.model3 = lmer(amplitude ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata,REML=FALSE)
amp.model4 = lmer(amplitude ~ comp_Channel_Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata,REML=FALSE)
amp.model5 = lmer(amplitude ~ comp_Channel_Group_Cond + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata,REML=FALSE)
amp.model6 = lmer(amplitude ~ comp_Channel_Group_Cond_Type + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata,REML=FALSE)
amp.anova <- anova(amp.null, amp.model1,amp.model2,amp.model3,amp.model4,amp.model5,amp.model6)
print(amp.anova)

amp.Group_channel<- glht(amp.model5, mcp(comp_Channel_Group_Cond= "Tukey"))
summary(amp.Group_channel)

amp.Group_channel_ty<- glht(amp.model6, mcp(comp_Channel_Group_Cond_Type= "Tukey"))
summary(amp.Group_channel_ty)

# because we used different sensors, separate RCaMP and GCaMP

#lck-GCaMP
Group_Cond_GC=interaction(stim.lck.compdata$Group[stim.lck.compdata$Channel=="GCaMP"],stim.lck.compdata$Condition[stim.lck.compdata$Channel=="GCaMP"])
Group_Cond_Type_GC=interaction(stim.lck.compdata$Group[stim.lck.compdata$Channel=="GCaMP"],
                               stim.lck.compdata$Condition[stim.lck.compdata$Channel=="GCaMP"],
                               stim.lck.compdata$ROIType[stim.lck.compdata$Channel=="GCaMP"])

amp.null.GC = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="GCaMP",],REML=FALSE)
amp.model2.GC  = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata[stim.lck.compdata$Channel=="GCaMP",],REML=FALSE)
amp.model3.GC  = lmer(amplitude ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="GCaMP",],REML=FALSE)
amp.model5.GC  = lmer(amplitude ~ Group_Cond_GC + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="GCaMP",],REML=FALSE)
amp.model6.GC  = lmer(amplitude ~ Group_Cond_Type_GC + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="GCaMP",],REML=FALSE)
amp.anova.GC  <- anova(amp.null.GC, amp.model2.GC,amp.model3.GC,amp.model5.GC,amp.model6.GC)
print(amp.anova.GC)

amp.Group_Cond.GC<- glht(amp.model5.GC, mcp(Group_Cond_GC= "Tukey"))
summary(amp.Group_Cond.GC)

amp.Group_channel_ty.GC<- glht(amp.model6.GC, mcp(Group_Cond_Type_GC= "Tukey"))
summary(amp.Group_channel_ty.GC)

#RCaMP
amp.null.RC = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="lck_RCaMP",],REML=FALSE)
amp.model2.RC  = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata[stim.lck.compdata$Channel=="lck_RCaMP",],REML=FALSE)
amp.anova.RC  <- anova(amp.null.RC, amp.model2.RC)
print(amp.anova.RC)

amp.Cond.RC<- glht(amp.model2.RC, mcp(Condition= "Tukey"))
summary(amp.Cond.RC)


#only consider STIM case
amp.null.stim = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="RCaMP",],REML=FALSE)
amp.model1.stim  = lmer(amplitude ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="RCaMP",],REML=FALSE)

amp.anova.stim  <- anova(amp.null.stim, amp.model1.stim)
print(amp.anova.stim)

amp.stim.group<- glht(amp.model1.stim, mcp(Group= "Tukey"))
summary(amp.stim.group)

########
#duration
# means for plots
df.Dur1A1<-summarySE(stim.cyto.alldata[stim.cyto.alldata$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Channel", "Condition"))
df.Dur1A2<-summarySE(cyto.GC.spot, measurevar = "meanDur", groupvars = c("Condition"))
# paired line plots
df.Dur1A2$Mean_meanDur<-df.Dur1A2$meanDur
df.Dur1A2$Dur_se<-df.Dur1A2$se
cyto.GC.spot<-merge(cyto.GC.spot, df.Dur1A2[, c("Condition", "Mean_meanDur","Dur_se")], by="Condition", all.x=TRUE)

df.Dur1B1<-summarySE(stim.lck.alldata[stim.lck.alldata$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Channel", "Condition"))
df.Dur1B2<-summarySE(lck.GC.spot, measurevar = "meanDur", groupvars = c("Condition"))
# paired line plots
df.Dur1B2$Mean_meanDur<-df.Dur1B2$meanDur
df.Dur1B2$Dur_se<-df.Dur1B2$se
lck.GC.spot<-merge(lck.GC.spot, df.Dur1B2[, c("Condition", "Mean_meanDur","Dur_se")], by="Condition", all.x=TRUE)

df.Dur3<- summarySE(stim.lck.compdata, measurevar = "Duration", groupvars = c("Channel_Group","Condition"))
df.Dur4<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim",], measurevar = "Duration", groupvars = c("Channel_Group"))
df.Dur5<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim"& stim.lck.compdata$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("ROIType","Channel_Group"))
df.Dur6<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="RCaMP",], measurevar = "Duration", groupvars = c("Channel_Group"))



# plots
ggplot(df.Dur1A1, aes(x=Condition,y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("cyto Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.Dur1A2, aes(x=Condition,y=meanDur, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanDur-se, ymax=meanDur+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("cyto Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

# paired plot
ggplot(cyto.GC.spot, aes(x=Condition, y = meanDur)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(group=Animal_Spot), colour="#b5b5b5")+
  geom_point(aes(x=Condition, y=Mean_meanDur), size = 5, colour="#1b7837")+
  geom_line(aes(x=Condition, y=Mean_meanDur,group=Animal_Spot), size=1.5, colour="#1b7837")+
  geom_errorbar(aes(x=Condition,ymin=Mean_meanDur-se, ymax=Mean_meanDur+se), colour="#1b7837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  ggtitle("cyto Dur per FOV")+
  max.theme



ggplot(df.Dur1B1, aes(x=Condition,y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("lck Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.Dur1B2, aes(x=Condition,y=meanDur, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanDur-se, ymax=meanDur+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("lck Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(lck.GC.spot, aes(x=Condition, y = meanDur)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(group=Animal_Spot), colour="#b5b5b5")+
  geom_point(aes(x=Condition, y=Mean_meanDur), size = 5, colour="#1b7837")+
  geom_line(aes(x=Condition, y=Mean_meanDur,group=Animal_Spot), size=1.5, colour="#1b7837")+
  geom_errorbar(aes(x=Condition,ymin=Mean_meanDur-se, ymax=Mean_meanDur+se), colour="#1b7837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  ggtitle("lck Dur per FOV")+
  max.theme

# NOT SIG DIFFERENT!!



ggplot(df.Dur3, aes(x=Channel_Group,y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.Dur4, aes(x=Channel_Group,y=Duration, fill= Channel_Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.Dur5, aes(x=Channel_Group,y=Duration, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.Dur6, aes(x=Channel_Group,y=Duration, fill=Channel_Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(stim.lck.compdata, aes(x=Channel_Group,y=Duration, fill=Channel_Group)) +
  geom_boxplot(notch=TRUE)+
  ylab("Duration") +
  ggtitle("lck data- fast vs delayed")+
  max.theme

ggplot(stim.lck.compdata, aes(x=Duration, y=..density..,fill=Channel_Group)) +
  geom_histogram(binwidth=0.5, position="dodge")+
  ylab("density") +
  xlim(-2,15)+
  ggtitle("lck data- fast vs delayed")+
  max.theme

ggplot(stim.lck.alldata[stim.lck.alldata$Channel_Group=="GCaMP.fast",], aes(x=Duration, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.5, position="dodge")+
  ylab("density") +
  xlim(-2,15)+
  ggtitle("lck data- fast no stim vs stim")+
  max.theme

# for supplementary figure
ggplot(stim.both.alldata[stim.both.alldata$Channel=="cyto_GCaMP",], aes(x=Duration, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=1, position="dodge")+
  ylab("density") +
  xlim(0,50)+
  ggtitle("cyto data- no stim vs stim")+
  max.theme

ggplot(stim.both.alldata[stim.both.alldata$Channel=="lck_GCaMP",], aes(x=Duration, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=1, position="dodge")+
  ylab("density") +
  xlim(0,50)+
  ggtitle("lck data- no stim vs stim")+
  max.theme


# supplementary figures

#cyto GC
Cond_Channel=interaction(stim.cyto.alldata$Condition,stim.cyto.alldata$Channel)
dur.cyto.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.cyto.alldata,REML=FALSE)
dur.cyto.model1 = lmer(Duration ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.cyto.alldata,REML=FALSE)
dur.cyto.model2 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.cyto.alldata,REML=FALSE)
dur.cyto.model3 = lmer(Duration ~ Cond_Channel + (1|Animal) + (1|Spot) + (1|trials), stim.cyto.alldata,REML=FALSE)

dur.cyto.anova <- anova(dur.cyto.null, dur.cyto.model1,dur.cyto.model2,dur.cyto.model3)
print(dur.cyto.anova)

dur.cyto.Cond_channel<- glht(dur.cyto.model3, mcp(Cond_Channel= "Tukey"))
summary(dur.cyto.Cond_channel)

#cyto GC spots
dur.cyto.spot.null = lmer(meanDur ~ (1|Animal), cyto.GC.spot,REML=FALSE)
dur.cyto.spot.model1 = lmer(meanDur ~ Condition + (1|Animal), cyto.GC.spot,REML=FALSE)
dur.cyto.spot.anova <- anova(dur.cyto.spot.null, dur.cyto.spot.model1)
print(dur.cyto.spot.anova)

dur.cyto.spot.pv<- glht(dur.cyto.spot.model1, mcp(Condition= "Tukey"))
summary(dur.cyto.spot.pv)

#lck GC
Cond_Channel=interaction(stim.lck.alldata$Condition,stim.lck.alldata$Channel)
dur.lck.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
dur.lck.model1 = lmer(Duration ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
dur.lck.model2 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.alldata,REML=FALSE)
dur.lck.model3 = lmer(Duration ~ Cond_Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)

dur.lck.anova <- anova(dur.lck.null, dur.lck.model1,dur.lck.model2,dur.lck.model3)
print(dur.lck.anova)

dur.lck.Cond_channel<- glht(dur.lck.model3, mcp(Cond_Channel= "Tukey"))
summary(dur.lck.Cond_channel)

#lck GC spots
dur.lck.spot.null = lmer(meanDur ~ (1|Animal), lck.GC.spot,REML=FALSE)
dur.lck.spot.model1 = lmer(meanDur ~ Condition + (1|Animal), lck.GC.spot,REML=FALSE)
dur.lck.spot.anova <- anova(dur.lck.spot.null, dur.lck.spot.model1)
print(dur.lck.spot.anova)

dur.lck.spot.pv<- glht(dur.lck.spot.model1, mcp(Condition= "Tukey"))
summary(dur.lck.spot.pv)


# stats for duration- neurons vs astrocytes
dur.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
dur.model1 = lmer(Duration ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
dur.model2 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.alldata,REML=FALSE)
dur.model3 = lmer(Duration ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
dur.model4 = lmer(Duration ~ Channel_Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
dur.model5 = lmer(Duration ~ Group_Channel_Cond + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
dur.model6 = lmer(Duration ~ Group_Channel_Cond_Type + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
dur.anova <- anova(dur.null, dur.model1,dur.model2,dur.model3,dur.model4,dur.model5,dur.model6)
print(dur.anova)

dur.Group_channel<- glht(dur.model5, mcp(Group_Channel_Cond= "Tukey"))
summary(dur.Group_channel)

dur.Group_channel_ty<- glht(dur.model6, mcp(Group_Channel_Cond_Type= "Tukey"))
summary(dur.Group_channel_ty)



# because we used different sensors, separate RCaMP and GCaMP

#lck-GCaMP
dur.null.GC = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="GCaMP",],REML=FALSE)
dur.model2.GC  = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata[stim.lck.compdata$Channel=="GCaMP",],REML=FALSE)
dur.model3.GC  = lmer(Duration ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="GCaMP",],REML=FALSE)
dur.model5.GC  = lmer(Duration ~ Group_Cond_GC + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="GCaMP",],REML=FALSE)
dur.model6.GC  = lmer(Duration ~ Group_Cond_Type_GC + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="GCaMP",],REML=FALSE)
dur.anova.GC  <- anova(dur.null.GC, dur.model2.GC,dur.model3.GC,dur.model5.GC,dur.model6.GC)
print(dur.anova.GC)

dur.Group_Cond.GC<- glht(dur.model5.GC, mcp(Group_Cond_GC= "Tukey"))
summary(dur.Group_Cond.GC)

dur.Group_channel_ty.GC<- glht(dur.model6.GC, mcp(Group_Cond_Type_GC= "Tukey"))
summary(dur.Group_channel_ty.GC)

#RCaMP
dur.null.RC = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="RCaMP",],REML=FALSE)
dur.model2.RC  = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata[stim.lck.compdata$Channel=="RCaMP",],REML=FALSE)
dur.anova.RC  <- anova(dur.null.RC, dur.model2.RC)
print(dur.anova.RC)

dur.Cond.RC<- glht(dur.model2.RC, mcp(Condition= "Tukey"))
summary(dur.Cond.RC)


#only consider STIM case
dur.null.stim = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="RCaMP",],REML=FALSE)
dur.model1.stim  = lmer(Duration ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="RCaMP",],REML=FALSE)

dur.anova.stim  <- anova(dur.null.stim, dur.model1.stim)
print(dur.anova.stim)

dur.stim.group<- glht(dur.model1.stim, mcp(Group= "Tukey"))
summary(dur.stim.group)



#######

# Process ROI area
df.Rarea1<-summarySE(stim.cyto.alldata[stim.cyto.alldata$ROIType=="Process",], measurevar = "area", groupvars = c("Channel", "Condition"))
df.Rarea2<-summarySE(stim.lck.alldata[stim.lck.alldata$ROIType=="Process",], measurevar = "area", groupvars = c("Channel", "Condition"))
df.Rarea3<- summarySE(stim.lck.compdata[stim.lck.compdata$ROIType=="Process",], measurevar = "area", groupvars = c("Channel_Group","Condition"))
df.Rarea4<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim"& stim.lck.compdata$ROIType=="Process",], measurevar = "area", groupvars = c("Channel_Group"))

ggplot(df.Rarea1, aes(x=Channel,y=area, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab(" cyto area") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.Rarea2, aes(x=Channel,y=area, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("lck area") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.Rarea3, aes(x=Channel_Group,y=area, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("area") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.Rarea4, aes(x=Channel_Group,y=area, fill= Channel_Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("area") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(stim.lck.compdata, aes(x=Channel_Group,y=area, fill=Channel_Group)) +
  geom_boxplot(notch=TRUE)+
  ylab("area") +
  ggtitle("lck data- fast vs delayed")+
  max.theme

ggplot(stim.lck.compdata, aes(x=area, y=..density..,fill=Channel_Group)) +
  geom_histogram(binwidth=10, position="dodge")+
  ylab("density") +
  xlim(-2,200)+
  ggtitle("lck data- fast vs delayed")+
  max.theme


#only consider STIM case
area.null.stim = lmer(area ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="RCaMP" & stim.lck.compdata.STIM$ROIType=="Process",],REML=FALSE)
area.model1.stim  = lmer(area ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="RCaMP" & stim.lck.compdata.STIM$ROIType=="Process",],REML=FALSE)

area.anova.stim  <- anova(area.null.stim, area.model1.stim)
print(area.anova.stim)

area.stim.group<- glht(area.model1.stim, mcp(Group= "Tukey"))
summary(area.stim.group)


######
# fraction of active pixels from all the astrocyte pixels

# only consider active based ROIs (i.e. processes) with peaks that have a peak time and onset time within the time windows

Lck.processes<- subset(stim.lck.alldata, ROIType=="Process" & Channel=="GCaMP")

# only ROIS with a signal near the stimulus
StimROI.lck.proc.Pixels<-ddply(Lck.processes, c("ROIs_Cond","trials","Animal","Spot","Animal_Spot","Condition","pixelsize", "nFluoPix"), summarise, 
                              nActivePix=sum(nActivePix))

# separate also the fast or delayed groups
groups.lck.proc.Pixels<-ddply(Lck.processes, c("ROIs_Cond","trials","Animal","Spot","Animal_Spot","Condition","Group","pixelsize", "nFluoPix"), summarise, 
                               nActivePix=sum(nActivePix))

# summarize per trial and spot
Trial.lck.proc.StimPixels<-ddply(StimROI.lck.proc.Pixels, c("trials","Animal","Spot","Animal_Spot","Condition","pixelsize", "nFluoPix"), summarise, nROIs=length(unique(ROIs_Cond)), 
                                nActivePix=sum(nActivePix))
Trial.lck.proc.StimPixels$FracActive=Trial.lck.proc.StimPixels$nActivePix/Trial.lck.proc.StimPixels$nFluoPix
Trial.lck.proc.StimPixels$FracActivePerROI=Trial.lck.proc.StimPixels$FracActive/Trial.lck.proc.StimPixels$nROIs


Trial.lck.group.StimPixels<-ddply(groups.lck.proc.Pixels, c("trials","Animal","Spot","Animal_Spot","Condition","Group","pixelsize", "nFluoPix"), summarise, nROIs=length(unique(ROIs_Cond)), 
                                nActivePix=sum(nActivePix))
Trial.lck.group.StimPixels$FracActive=Trial.lck.group.StimPixels$nActivePix/Trial.lck.group.StimPixels$nFluoPix
Trial.lck.group.StimPixels$FracActivePerROI=Trial.lck.group.StimPixels$FracActive/Trial.lck.group.StimPixels$nROIs


Spot.lck.StimPixels<-ddply(Trial.lck.proc.StimPixels, c("Animal","Spot","Animal_Spot","Condition","pixelsize"), summarise, nTrials=length(unique(trials)),
                           meanTotalPix=mean(nFluoPix), MeanFracActive=mean(FracActive), MeanFracActPerROI=mean(FracActivePerROI))

Spot.lck.group.StimPixels<-ddply(Trial.lck.group.StimPixels, c("Animal","Spot","Animal_Spot","Condition","Group","pixelsize"), summarise, nTrials=length(unique(trials)),
                                 meanTotalPix=mean(nFluoPix), MeanFracActive=mean(FracActive), MeanFracActPerROI=mean(FracActivePerROI))

#fracActive per trial
df.FracActive.lck<- summarySE(Trial.lck.proc.StimPixels, measurevar = "FracActive", groupvars = c("Condition"))
df.FracActive.lck.group<- summarySE(Trial.lck.group.StimPixels, measurevar = "FracActive", groupvars = c("Condition","Group"))

#fracActive per ROI per trial
df.FracActivePerROI.lck<- summarySE(Trial.lck.proc.StimPixels, measurevar = "FracActivePerROI", groupvars = c("Condition"))
df.FracActivePerROI.lck.group<- summarySE(Trial.lck.group.StimPixels, measurevar = "FracActivePerROI", groupvars = c("Condition","Group"))

#fracActive per Spot
df.FracActive.Spot<- summarySE(Spot.lck.StimPixels, measurevar = "MeanFracActive", groupvars = c("Condition"))
df.FracActivePerROI.Spot<- summarySE(Spot.lck.StimPixels, measurevar = "MeanFracActPerROI", groupvars = c("Condition"))
df.FracActive.Spot2<- summarySE(Spot.lck.group.StimPixels, measurevar = "MeanFracActive", groupvars = c("Condition","Group"))

#fraction active  (includes ROI size and number)
ggplot(data=df.FracActive.lck, aes(x=Condition, y= FracActive, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=FracActive-se, ymax=FracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Fraction of Active Astrocyte Process Area") +
  ggtitle("frac act per trial") +
  scale_fill_manual(values=c("#b4e2a8", "#1b7837")) + 
  max.theme

ggplot(data=df.FracActive.Spot, aes(x=Condition, y= MeanFracActive, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=MeanFracActive-se, ymax=MeanFracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Fraction of Active Astrocyte Process Area") +
  ggtitle("MeanFracActive per spot") +
  scale_fill_manual(values=c("#b4e2a8", "#1b7837")) + 
  max.theme

#fraction active per ROI  (only a factor of ROI size)
ggplot(data=df.FracActivePerROI.lck, aes(x=Condition, y= FracActivePerROI, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=FracActivePerROI-se, ymax=FracActivePerROI+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActivePerROI") +
  ggtitle("FracActivePerROI per trial") +
  scale_fill_manual(values=c("#b4e2a8", "#1b7837")) + 
  max.theme

ggplot(data=df.FracActivePerROI.Spot, aes(x=Condition, y= MeanFracActPerROI, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=MeanFracActPerROI-se, ymax=MeanFracActPerROI+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Fraction of Active Astrocyte Process Area") +
  ggtitle("MeanFracActPerROI per spot") +
  scale_fill_manual(values=c("#b4e2a8", "#1b7837")) + 
  max.theme

# fraction active by group
ggplot(data=df.FracActive.lck.group, aes(x=Group, y= FracActive, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=FracActive-se, ymax=FracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive") +
  ggtitle("FracActive per trial") +
  scale_fill_manual(values=c("#b4e2a8", "#1b7837")) + 
  max.theme

ggplot(data=df.FracActive.Spot2, aes(x=Group, y= MeanFracActive, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=MeanFracActive-se, ymax=MeanFracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive") +
  ggtitle("FracActive per trial") +
  scale_fill_manual(values=c("#b4e2a8", "#1b7837")) + 
  max.theme



# paired line plots
df.FracActive.Spot$MeanFracActiveTotal<-df.FracActive.Spot$MeanFracActive
Spot.lck.StimPixels<-merge(Spot.lck.StimPixels, df.FracActive.Spot[, c("Condition", "MeanFracActiveTotal","se")], by="Condition", all.x=TRUE)


# paired plot
ggplot(Spot.lck.StimPixels, aes(x=Condition, y = MeanFracActive)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(group=Animal_Spot), colour="#b5b5b5")+
  geom_point(aes(x=Condition, y=MeanFracActiveTotal), size = 5, colour="#1b7837")+
  geom_line(aes(x=Condition, y=MeanFracActiveTotal,group=Animal_Spot), size=1.5, colour="#1b7837")+
  geom_errorbar(aes(x=Condition,ymin=MeanFracActiveTotal-se, ymax=MeanFracActiveTotal+se), colour="#1b7837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  max.theme


## Stats

#lck

fracActive.lck.stim.null = lmer(MeanFracActive ~ (1|Animal) , Spot.lck.StimPixels,REML=FALSE)
fracActive.lck.stim.model1 = lmer(MeanFracActive ~ Condition + (1|Animal) , Spot.lck.StimPixels,REML=FALSE)
fracActive.lck.stim.anova  <- anova(fracActive.lck.stim.null, fracActive.lck.stim.model1)
print(fracActive.lck.stim.anova)

fracActive.lck.stim.pvalue<- glht(fracActive.lck.stim.model1, mcp(Condition= "Tukey"))
summary(fracActive.lck.stim.pvalue)


# groups
group_Cond= interaction(Spot.lck.group.StimPixels$Group, Spot.lck.group.StimPixels$Condition)
fracActive.lck.group.null = lmer(MeanFracActive ~ (1|Animal) , Spot.lck.group.StimPixels,REML=FALSE)
fracActive.lck.group.model1 = lmer(MeanFracActive ~ Condition + (1|Animal) , Spot.lck.group.StimPixels,REML=FALSE)
fracActive.lck.group.model2 = lmer(MeanFracActive ~ Group + (1|Animal) , Spot.lck.group.StimPixels,REML=FALSE)
fracActive.lck.group.model3 = lmer(MeanFracActive ~ group_Cond + (1|Animal) , Spot.lck.group.StimPixels,REML=FALSE)
fracActive.lck.group.anova  <- anova(fracActive.lck.group.null, fracActive.lck.group.model1,
                                     fracActive.lck.group.model2, fracActive.lck.group.model3)
print(fracActive.lck.group.anova)

fracActive.lck.group.pvalue<- glht(fracActive.lck.group.model3, mcp(group_Cond= "Tukey"))
summary(fracActive.lck.group.pvalue)



# cyto GCaMP
# only consider active based ROIs (i.e. processes) with peaks that have a peak time and onset time within the time windows

cyto.processes<- subset(stim.cyto.alldata, ROIType=="Process" & Channel=="GCaMP")

# only ROIS with a signal near the stimulus
StimROI.cyto.proc.Pixels<-ddply(cyto.processes, c("ROIs_Cond","trials","Animal","Spot","Animal_Spot","Condition","pixelsize", "nFluoPix"), summarise, 
                               nActivePix=sum(nActivePix))

# summarize per trial and spot
Trial.cyto.proc.StimPixels<-ddply(StimROI.cyto.proc.Pixels, c("trials","Animal","Spot","Animal_Spot","Condition","pixelsize", "nFluoPix"), summarise, nROIs=length(unique(ROIs_Cond)), 
                                 nActivePix=sum(nActivePix))
Trial.cyto.proc.StimPixels$FracActive=Trial.cyto.proc.StimPixels$nActivePix/Trial.cyto.proc.StimPixels$nFluoPix
Trial.cyto.proc.StimPixels$FracActivePerROI=Trial.cyto.proc.StimPixels$FracActive/Trial.cyto.proc.StimPixels$nROIs


Spot.cyto.StimPixels<-ddply(Trial.cyto.proc.StimPixels, c("Animal","Spot","Animal_Spot","Condition","pixelsize"), summarise, nTrials=length(unique(trials)),
                            meanTotalPix=mean(nFluoPix),MeanFracActive=mean(FracActive), MeanFracActPerROI=mean(FracActivePerROI))


#fracActive per trial
df.FracActive.cyto<- summarySE(Trial.cyto.proc.StimPixels, measurevar = "FracActive", groupvars = c("Condition"))


#fracActive per ROI per trial
df.FracActivePerROI.cyto<- summarySE(Trial.cyto.proc.StimPixels, measurevar = "FracActivePerROI", groupvars = c("Condition"))

#fracActive per Spot
df.FracActive.cyto.Spot<- summarySE(Spot.cyto.StimPixels, measurevar = "MeanFracActive", groupvars = c("Condition"))
df.FracActivePerROI.cyto.Spot<- summarySE(Spot.cyto.StimPixels, measurevar = "MeanFracActPerROI", groupvars = c("Condition"))


#fraction active  (includes ROI size and number)
ggplot(data=df.FracActive.cyto, aes(x=Condition, y= FracActive, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=FracActive-se, ymax=FracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Fraction of Active Astrocyte Process Area") +
  ggtitle("frac act per trial") +
  scale_fill_manual(values=c("#b4e2a8", "#1b7837")) + 
  max.theme

ggplot(data=df.FracActive.cyto.Spot, aes(x=Condition, y= MeanFracActive, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=MeanFracActive-se, ymax=MeanFracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Fraction of Active Astrocyte Process Area") +
  ggtitle("MeanFracActive per spot") +
  scale_fill_manual(values=c("#b4e2a8", "#1b7837")) + 
  max.theme

#fraction active per ROI  (only a factor of ROI size)
ggplot(data=df.FracActivePerROI.cyto, aes(x=Condition, y= FracActivePerROI, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=FracActivePerROI-se, ymax=FracActivePerROI+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActivePerROI") +
  ggtitle("FracActivePerROI per trial") +
  scale_fill_manual(values=c("#b4e2a8", "#1b7837")) + 
  max.theme

ggplot(data=df.FracActivePerROI.cyto.Spot, aes(x=Condition, y= MeanFracActPerROI, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=MeanFracActPerROI-se, ymax=MeanFracActPerROI+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Fraction of Active Astrocyte Process Area") +
  ggtitle("MeanFracActPerROI per spot") +
  scale_fill_manual(values=c("#b4e2a8", "#1b7837")) + 
  max.theme


# paired line plots
df.FracActive.cyto.Spot$MeanFracActiveTotal<-df.FracActive.cyto.Spot$MeanFracActive
Spot.cyto.StimPixels<-merge(Spot.cyto.StimPixels, df.FracActive.cyto.Spot[, c("Condition", "MeanFracActiveTotal","se")], by="Condition", all.x=TRUE)


# paired plot
ggplot(Spot.cyto.StimPixels, aes(x=Condition, y = MeanFracActive)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(group=Animal_Spot), colour="#b5b5b5")+
  geom_point(aes(x=Condition, y=MeanFracActiveTotal), size = 5, colour="#1b7837")+
  geom_line(aes(x=Condition, y=MeanFracActiveTotal,group=Animal_Spot), size=1.5, colour="#1b7837")+
  geom_errorbar(aes(x=Condition,ymin=MeanFracActiveTotal-se, ymax=MeanFracActiveTotal+se), colour="#1b7837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  max.theme

## Stats

#cyto

fracActive.cyto.stim.null = lmer(MeanFracActive ~ (1|Animal) , Spot.cyto.StimPixels,REML=FALSE)
fracActive.cyto.stim.model1 = lmer(MeanFracActive ~ Condition + (1|Animal) , Spot.cyto.StimPixels,REML=FALSE)
fracActive.cyto.stim.anova  <- anova(fracActive.cyto.stim.null, fracActive.cyto.stim.model1)
print(fracActive.cyto.stim.anova)

fracActive.cyto.stim.pvalue<- glht(fracActive.cyto.stim.model1, mcp(Condition= "Tukey"))
summary(fracActive.cyto.stim.pvalue)



################
#Robustness scores

lck.robustness<-read.table("D:/Data/GCaMP_RCaMP/Revision/Lck_GCaMP/FilesforR/Lck_robustnessScores.csv", header=TRUE, sep = ",")
lck.robustness$Animal_Spot<- paste(lck.robustness$All_traces5, lck.robustness$All_traces4, sep="_")
lck.robustness$Condition<- lck.robustness$All_traces6
lck.robustness$Ani_Cond<-paste(lck.robustness$Animal_Spot, lck.robustness$Condition, sep="_")

#remove IP3KOs and IP3WTs
KOs<-c("IPRG1","IPRG4","IPRG5","IPRG7")
WTs<-c("IPRG2","IPRG3","IPRG6")
lck.robustness<-subset(lck.robustness, !(All_traces5 %in% KOs))
lck.robustness<-subset(lck.robustness, !(All_traces5 %in% WTs))

# incorporate total number of astrocyte pixels for each FOV
Spot.lck.StimPixels$Ani_Cond<-paste(Spot.lck.StimPixels$Animal_Spot, Spot.lck.StimPixels$Condition, sep="_")

lck.robustness<-merge(lck.robustness, Spot.lck.StimPixels[, c("Ani_Cond", "meanTotalPix")], by="Ani_Cond", all.x=TRUE)
lck.robustness$FracPixAboveThres<-lck.robustness$nPxAboveThresh/lck.robustness$meanTotalPix

df.lck.robust1<-summarySE(lck.robustness, measurevar = "score", groupvars = c("Condition"))
df.lck.robust2<-summarySE(lck.robustness, measurevar = "FracPixAboveThres", groupvars = c("Condition"))


# paired line plots
df.lck.robust1$MeanFracofROIarea<-df.lck.robust1$score
df.lck.robust1$MeanFracofROIarea_se<-df.lck.robust1$se
lck.robustness<-merge(lck.robustness, df.lck.robust1[, c("Condition", "MeanFracofROIarea","MeanFracofROIarea_se")], by="Condition", all.x=TRUE)
df.lck.robust2$MeanFracofACarea<-df.lck.robust2$FracPixAboveThres
df.lck.robust2$MeanFracofACarea_se<-df.lck.robust2$se
lck.robustness<-merge(lck.robustness, df.lck.robust2[, c("Condition", "MeanFracofACarea","MeanFracofACarea_se")], by="Condition", all.x=TRUE)


# paired plot
ggplot(lck.robustness, aes(x=Condition, y = score)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(group=Animal_Spot), colour="#b5b5b5")+
  geom_point(aes(x=Condition, y=MeanFracofROIarea), size = 5, colour="#1b7837")+
  geom_line(aes(x=Condition, y=MeanFracofROIarea,group=Animal_Spot), size=1.5, colour="#1b7837")+
  geom_errorbar(aes(x=Condition,ymin=MeanFracofROIarea-MeanFracofROIarea_se, ymax=MeanFracofROIarea+MeanFracofROIarea_se), colour="#1b7837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  ggtitle("Fraction of Repeated Pixels per Total ROI area")+
  max.theme

ggplot(lck.robustness, aes(x=Condition, y = FracPixAboveThres)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(group=Animal_Spot), colour="#b5b5b5")+
  geom_point(aes(x=Condition, y=MeanFracofACarea), size = 5, colour="#1b7837")+
  geom_line(aes(x=Condition, y=MeanFracofACarea,group=Animal_Spot), size=1.5, colour="#1b7837")+
  geom_errorbar(aes(x=Condition,ymin=MeanFracofACarea-MeanFracofACarea_se, ymax=MeanFracofACarea+MeanFracofACarea_se), colour="#1b7837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  ggtitle("Fraction of Repeated Pixels per Total lck astrocyte area")+
  max.theme



# stats
# fraction of pixels vs total ROI area
lck.robust.null = lmer(score ~ (1|All_traces5) , lck.robustness, REML=FALSE)
lck.robust.model1 = lmer(score ~ Condition + (1|All_traces5) , lck.robustness,REML=FALSE)
lck.robust.anova  <- anova(lck.robust.null, lck.robust.model1)
print(lck.robust.anova)

lck.robust.pvalue<- glht(lck.robust.model1, mcp(Condition= "Tukey"))
summary(lck.robust.pvalue)

# fraction of pixels vs total astrocyte area
lck.robust.totalAC.null = lmer(FracPixAboveThres ~ (1|All_traces5) , lck.robustness, REML=FALSE)
lck.robust.totalAC.model1 = lmer(FracPixAboveThres ~ Condition + (1|All_traces5) , lck.robustness,REML=FALSE)
lck.robust.totalAC.anova  <- anova(lck.robust.totalAC.null, lck.robust.totalAC.model1)
print(lck.robust.totalAC.anova)

lck.robust.totalAC.pvalue<- glht(lck.robust.totalAC.model1, mcp(Condition= "Tukey"))
summary(lck.robust.totalAC.pvalue)




# cyto
cyto.robustness<-read.table("D:/Data/GCaMP_RCaMP/Revision/cytoGCaMP/FilesforR/cyto_robustnessScores.csv", header=TRUE, sep = ",")

cyto.robustness$Animal_Spot<- paste(cyto.robustness$All_traces5, cyto.robustness$All_traces4, sep="_")
cyto.robustness$Condition<- cyto.robustness$All_traces6

# incorporate total number of astrocyte pixels for each FOV
Spot.cyto.StimPixels$Ani_Cond<-paste(Spot.cyto.StimPixels$Animal_Spot, Spot.cyto.StimPixels$Condition, sep="_")
cyto.robustness$Ani_Cond<-paste(cyto.robustness$Animal_Spot, cyto.robustness$Condition, sep="_")
cyto.robustness<-merge(cyto.robustness, Spot.cyto.StimPixels[, c("Ani_Cond", "meanTotalPix")], by="Ani_Cond", all.x=TRUE)

cyto.robustness$FracPixAboveThres<-cyto.robustness$nPxAboveThresh/cyto.robustness$meanTotalPix

df.cyto.robust1<-summarySE(cyto.robustness, measurevar = "score", groupvars = c("Condition"))
df.cyto.robust2<-summarySE(cyto.robustness, measurevar = "FracPixAboveThres", groupvars = c("Condition"))

# paired line plots
df.cyto.robust1$MeanFracofROIarea<-df.cyto.robust1$score
df.cyto.robust1$MeanFracofROIarea_se<-df.cyto.robust1$se
cyto.robustness<-merge(cyto.robustness, df.cyto.robust1[, c("Condition", "MeanFracofROIarea","MeanFracofROIarea_se")], by="Condition", all.x=TRUE)
df.cyto.robust2$MeanFracofACarea<-df.cyto.robust2$FracPixAboveThres
df.cyto.robust2$MeanFracofACarea_se<-df.cyto.robust2$se
cyto.robustness<-merge(cyto.robustness, df.cyto.robust2[, c("Condition", "MeanFracofACarea","MeanFracofACarea_se")], by="Condition", all.x=TRUE)


# paired plot
ggplot(cyto.robustness, aes(x=Condition, y = score)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(group=Animal_Spot), colour="#b5b5b5")+
  geom_point(aes(x=Condition, y=MeanFracofROIarea), size = 5, colour="#1b7837")+
  geom_line(aes(x=Condition, y=MeanFracofROIarea,group=Animal_Spot), size=1.5, colour="#1b7837")+
  geom_errorbar(aes(x=Condition,ymin=MeanFracofROIarea-MeanFracofROIarea_se, ymax=MeanFracofROIarea+MeanFracofROIarea_se), colour="#1b7837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  ggtitle("Fraction of Repeated Pixels per Total ROI area")+
  max.theme

ggplot(cyto.robustness, aes(x=Condition, y = FracPixAboveThres)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(group=Animal_Spot), colour="#b5b5b5")+
  geom_point(aes(x=Condition, y=MeanFracofACarea), size = 5, colour="#1b7837")+
  geom_line(aes(x=Condition, y=MeanFracofACarea,group=Animal_Spot), size=1.5, colour="#1b7837")+
  geom_errorbar(aes(x=Condition,ymin=MeanFracofACarea-MeanFracofACarea_se, ymax=MeanFracofACarea+MeanFracofACarea_se), colour="#1b7837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  ggtitle("Fraction of Repeated Pixels per Total cyto astrocyte area")+
  max.theme

# stats
# fraction of pixels above threshold for total ROI area
cyto.robust.null = lmer(score ~ (1|All_traces5) , cyto.robustness, REML=FALSE)
cyto.robust.model1 = lmer(score ~ Condition + (1|All_traces5) , cyto.robustness,REML=FALSE)
cyto.robust.anova  <- anova(cyto.robust.null, cyto.robust.model1)
print(cyto.robust.anova)

cyto.robust.pvalue<- glht(cyto.robust.model1, mcp(Condition= "Tukey"))
summary(cyto.robust.pvalue)


# fraction of pixels vs total astrocyte area
cyto.robust.totalAC.null = lmer(FracPixAboveThres ~ (1|All_traces5) , cyto.robustness, REML=FALSE)
cyto.robust.totalAC.model1 = lmer(FracPixAboveThres ~ Condition + (1|All_traces5) , cyto.robustness,REML=FALSE)
cyto.robust.totalAC.anova  <- anova(cyto.robust.totalAC.null, cyto.robust.totalAC.model1)
print(cyto.robust.totalAC.anova)

cyto.robust.totalAC.pvalue<- glht(cyto.robust.totalAC.model1, mcp(Condition= "Tukey"))
summary(cyto.robust.totalAC.pvalue)

##### 
# correlation data 

CorrData <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LongStim_Correlations.csv", header=TRUE, sep = ",")
CorrData <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LongStim_Correlations_fixedDis.csv", header=TRUE, sep = ",")


CorrData$ROI_Y_type<-as.character(CorrData$ROI_Y_type)
CorrData$ROI_Y_type[grepl("D",CorrData$ROI_Y)]="Dendrite"
CorrData$ROI_Y_type<-as.factor(CorrData$ROI_Y_type)
CorrData$CompChannel<-paste(CorrData$ChannelX, CorrData$ChannelY, sep= "_")
CorrData$CompType<-paste(CorrData$ROI_X_type, CorrData$ROI_Y_type, sep= "_")

# GCaMP and RCaMP ROIs
GCaMP_RCaMP<-subset(CorrData, CompChannel=="GCaMP_RCaMP")
GCaMP_RCaMP$ROIs_trial<-paste(GCaMP_RCaMP$Animal, GCaMP_RCaMP$Spot, GCaMP_RCaMP$Trial,GCaMP_RCaMP$ROI_X, sep= "_")
GCaMP_RCaMP$RCaMP_ROIs<-paste(GCaMP_RCaMP$Animal, GCaMP_RCaMP$Spot, GCaMP_RCaMP$Trial,GCaMP_RCaMP$ROI_Y, sep= "_")

#remove random spots with gcamp neurons
ntypes=c("Neuron_Neuron","Neuron_Neuropil","Neuropil_Neuron","Neuropil_Neuropil")
GCaMP_RCaMP<-subset(GCaMP_RCaMP, !(CompType %in% ntypes))



# consider only comparisons between astrocyte and neuronal ROIs that respond to stimulation 

respGC<-subset(stim.both.alldata, Channel=="lck_GCaMP" & Condition=="Stim") 
respRC<-subset(stim.both.alldata, Channel=="lck_RCaMP" & Condition=="Stim") 

# list of responding astrocytes and their corresponding group
respGC2<-unique(respGC[c("ROIs_trial","Group")])
respRC2<-unique(respRC[c("ROIs_trial","Group")])

# only correlations from responding astrocytes and neurons
GCaMP_RCaMP<-subset(GCaMP_RCaMP, ROIs_trial %in% respGC2$ROIs_trial)
GCaMP_RCaMP<-subset(GCaMP_RCaMP, RCaMP_ROIs %in% unique(respRC$ROIs_trial))

# put group information into correlation table ("fast or delayed")
GCaMP_RCaMP.group<-merge(GCaMP_RCaMP, respGC2[, c("ROIs_trial", "Group")], by="ROIs_trial")


ggplot(GCaMP_RCaMP.group, aes(x=Short_Corr, fill=Group)) + geom_density(alpha=0.3)+
  ggtitle("density-") + max.theme

ggplot(GCaMP_RCaMP.group, aes(x=MinDistance, y=Short_Corr, colour=Group))+
  geom_point(alpha=0.3) +
  max.theme

ggplot(GCaMP_RCaMP.group, aes(x=MinDistance, y=Short_Corr, color=Group)) + 
  geom_point(alpha=0.3) + 
  geom_density2d()+
  max.theme

ggplot(GCaMP_RCaMP.group[GCaMP_RCaMP.group$Group=="fast",], aes(x=MinDistance, y=Short_Corr)) + 
  geom_point() + 
  geom_smooth(method='lm')+
  max.theme


#  means considering all correlations
GCaMP_RCaMP.group$Group <- factor(GCaMP_RCaMP.group$Group, levels = c("fast","delayed"))

#########
# significant correlation (and respond to stimulation)

# comparisons with significant P values
GR_Corr_Sig<-subset(GCaMP_RCaMP.group, ShortPvalue<0.05)

ggplot(GR_Corr_Sig, aes(x=Short_Corr, fill=Group)) + geom_density(alpha=0.3)+
  ggtitle("density-") + max.theme

ggplot(GR_Corr_Sig, aes(x=MinDistance, y=Short_Corr, colour=Group))+
  geom_point(alpha=0.3) +
  max.theme


df5A1<- summarySE(GCaMP_RCaMP.group, measurevar="Short_Corr", groupvars=c("Group"))
df5A2<- summarySE(GCaMP_RCaMP.group, measurevar="Short_Corr", groupvars=c("CompType"))
df5A3<- summarySE(GCaMP_RCaMP.group, measurevar="Short_Corr", groupvars=c("Group","CompType"))

df5B1<- summarySE(GR_Corr_Sig, measurevar="Short_Corr", groupvars=c("Group"))
df5B2<- summarySE(GR_Corr_Sig, measurevar="Short_Corr", groupvars=c("CompType"))
df5B3<- summarySE(GR_Corr_Sig, measurevar="Short_Corr", groupvars=c("Group","CompType"))

ggplot(data=df5A1, aes(x=Group, y=Short_Corr, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Short_Corr for all responding GCaMP") +
  max.theme

ggplot(data=df5A2, aes(x=CompType, y=Short_Corr, fill=CompType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Short_Corr for all responding GCaMP") +
  max.theme

ggplot(data=df5A3, aes(x=CompType, y=Short_Corr, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Short_Corr for all responding GCaMP") +
  max.theme

ggplot(data=df5B1, aes(x=Group, y=Short_Corr, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Short_Corr for sig responding GCaMP") +
  max.theme

ggplot(data=df5B2, aes(x=CompType, y=Short_Corr, fill=CompType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Short_Corr for sig responding GCaMP") +
  max.theme

ggplot(data=df5B3, aes(x=CompType, y=Short_Corr, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Short_Corr for sig responding GCaMP") +
  max.theme


df6A1<- summarySE(GCaMP_RCaMP.group, measurevar="xCorr", groupvars=c("Group"))
df6A2<- summarySE(GCaMP_RCaMP.group, measurevar="xCorr", groupvars=c("CompType"))
df6A3<- summarySE(GCaMP_RCaMP.group, measurevar="xCorr", groupvars=c("Group","CompType"))

df6B1<- summarySE(GR_Corr_Sig, measurevar="xCorr", groupvars=c("Group"))
df6B2<- summarySE(GR_Corr_Sig, measurevar="xCorr", groupvars=c("CompType"))
df6B3<- summarySE(GR_Corr_Sig, measurevar="xCorr", groupvars=c("Group","CompType"))

ggplot(data=df6A1, aes(x=Group, y=xCorr, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("xCorr") +
  ggtitle("xCorr for all responding GCaMP") +
  max.theme

ggplot(data=df6A2, aes(x=CompType, y=xCorr, fill=CompType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("xCorr") +
  ggtitle("xCorr for all responding GCaMP") +
  max.theme

ggplot(data=df6A3, aes(x=CompType, y=xCorr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("xCorr") +
  ggtitle("xCorr for all responding GCaMP") +
  max.theme

ggplot(data=df6B1, aes(x=Group, y=xCorr, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("xCorr") +
  ggtitle("xCorr for sig responding GCaMP") +
  max.theme

ggplot(data=df6B2, aes(x=CompType, y=xCorr, fill=CompType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("xCorr") +
  ggtitle("xCorr for sig responding GCaMP") +
  max.theme

ggplot(data=df6B3, aes(x=CompType, y=xCorr, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("xCorr") +
  ggtitle("xCorr for sig responding GCaMP") +
  max.theme


df7A1<- summarySE(GCaMP_RCaMP.group, measurevar="Lag", groupvars=c("Group"))
df7A2<- summarySE(GCaMP_RCaMP.group, measurevar="Lag", groupvars=c("CompType"))
df7A3<- summarySE(GCaMP_RCaMP.group, measurevar="Lag", groupvars=c("Group","CompType"))

df7B1<- summarySE(GR_Corr_Sig, measurevar="Lag", groupvars=c("Group"))
df7B2<- summarySE(GR_Corr_Sig, measurevar="Lag", groupvars=c("CompType"))
df7B3<- summarySE(GR_Corr_Sig, measurevar="Lag", groupvars=c("Group","CompType"))

ggplot(data=df7A1, aes(x=Group, y=Lag, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Lag") +
  ggtitle("Lag for all responding GCaMP") +
  max.theme

ggplot(data=df7A2, aes(x=CompType, y=Lag, fill=CompType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Lag") +
  ggtitle("Lag for all responding GCaMP") +
  max.theme

ggplot(data=df7A3, aes(x=CompType, y=Lag, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Lag") +
  ggtitle("Lag for all responding GCaMP") +
  max.theme

ggplot(data=df7B1, aes(x=GroupX, y=Lag, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Lag") +
  ggtitle("Lag for sig responding GCaMP") +
  max.theme

ggplot(data=df7B2, aes(x=CompType, y=Lag, fill=CompType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Lag") +
  ggtitle("Lag for sig responding GCaMP") +
  max.theme

ggplot(data=df7B3, aes(x=CompType, y=Lag, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Lag") +
  ggtitle("Lag for sig responding GCaMP") +
  max.theme

#######
# mean distance
df8A1<- summarySE(GCaMP_RCaMP.group, measurevar="MinDistance", groupvars=c("Group"))
df8A2<- summarySE(GCaMP_RCaMP.group, measurevar="MinDistance", groupvars=c("CompType"))
df8A3<- summarySE(GCaMP_RCaMP.group, measurevar="MinDistance", groupvars=c("Group","CompType"))

df8B1<- summarySE(GR_Corr_Sig, measurevar="MinDistance", groupvars=c("Group"))
df8B2<- summarySE(GR_Corr_Sig, measurevar="MinDistance", groupvars=c("CompType"))
df8B3<- summarySE(GR_Corr_Sig, measurevar="MinDistance", groupvars=c("Group","CompType"))

ggplot(data=df8A1, aes(x=Group, y=MinDistance, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=MinDistance-se, ymax=MinDistance+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("MinDistance") +
  ggtitle("MinDistance for all responding GCaMP") +
  max.theme

ggplot(data=df8B1, aes(x=Group, y=MinDistance, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=MinDistance-se, ymax=MinDistance+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("MinDistance") +
  ggtitle("MinDistance for all sig GCaMP") +
  max.theme


# stats 
minDis.lck.null = lmer(MinDistance ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP_RCaMP.group,REML=FALSE)
minDis.lck.model1 = lmer(MinDistance ~ Group + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP_RCaMP.group,REML=FALSE)
minDis.lck.anova <- anova(minDis.lck.null, minDis.lck.model1)
print(minDis.lck.anova)

#sig ROIs
minDis2.lck.null = lmer(MinDistance ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), GR_Corr_Sig,REML=FALSE)
minDis2.lck.model1 = lmer(MinDistance ~ Group + (1|Animal) + (1|Spot) + (1|ROIs_trial), GR_Corr_Sig,REML=FALSE)
minDis2.lck.anova <- anova(minDis2.lck.null, minDis2.lck.model1)
print(minDis2.lck.anova)

minDis.lck.group<- glht(minDis.lck.model1, mcp(Group= "Tukey"))
summary(minDis.lck.group)


######
# distance to the nearest active neuron
# NOTE: no effect (see below)- distance is similar between groups and ROI types

respGC$minDistance=NA
for (ii in 1:nrow(respGC))
{
  ROIx=respGC$ROIs_trial[ii]
  subset1=subset(GCaMP_RCaMP.groups, ROIs_trial==ROIx)
  if (nrow(subset1)>0)
  {
    respGC$minDistance[ii]=min(subset1$MinDistance)
  }
}

dfDis<- summarySE(respGC, measurevar="minDistance", groupvars=c("Group"),na.rm = TRUE)
dfDis2<- summarySE(respGC, measurevar="minDistance", groupvars=c("ROIType"),na.rm = TRUE)
dfDis3<- summarySE(respGC, measurevar="minDistance", groupvars=c("ROIType","Group"),na.rm = TRUE)


ggplot(stim.lck.alldata[stim.lck.alldata$Channel_Group=="lck_GCaMP.fast",], aes(x=area, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=5, position="dodge")+
  ylab("density") +
  xlim(-2,200)+
  ggtitle("lck data- fast no stim vs stim")+
  max.theme

Group_ROIType=interaction(respGC$Group, respGC$ROIType)
# stats for onset times- neurons vs astrocytes
minDis.lck.type.null = lmer(minDistance ~ (1|Animal) + (1|Spot) + (1|trials), respGC,REML=FALSE)
minDis.lck.type.model1 = lmer(minDistance ~ Group + (1|Animal) + (1|Spot) + (1|trials), respGC,REML=FALSE)
minDis.lck.type.model2 = lmer(minDistance ~ ROIType + (1|Animal) + (1|Spot) + (1|trials), respGC,REML=FALSE)
minDis.lck.type.model3 = lmer(minDistance ~ Group_ROIType + (1|Animal) + (1|Spot) + (1|trials),respGC,REML=FALSE)
minDis.lck.type.anova <- anova(minDis.lck.type.null, minDis.lck.type.model1, 
                               minDis.lck.type.model2,minDis.lck.type.model3)
print(minDis.lck.type.anova)



########
#mean peak amplitude for each ROI with peaks from the stim window
# use this to calculate high, mid and low responding neuronal somas

ROIwise<- ddply(all.lck.peaks, c("Animal", "Spot","Trial","Channel","roiName","Condition","trials", "ROIs_trial" ,"ROIType"), 
                summarise, meanAmp=mean(amplitude), meanDur= mean(Duration), meanProm=mean(prominence),
                meanPAUC=mean(peakAUC), nPeaks= length(amplitude))

ROIwise$roiNameUnique<-paste(ROIwise$Animal,ROIwise$Spot,ROIwise$roiName, sep="_")
ROIwise.Neurons<-subset(ROIwise, Condition!="Nostim" & ROIType=="Neuron")
ROIwise.Dendrites<-subset(ROIwise, Condition!="Nostim" & ROIType=="Dendrite")

# MEAN info for each ROI for trials of short OR long stim
NeuronSomas<-ddply(ROIwise.Neurons, c("Animal", "Spot", "roiNameUnique"),
                   summarise, meanAmp2=mean(meanAmp), meanProm2=mean(meanProm),
                   meanDur2=mean(meanDur), meanPAUC2=mean(meanPAUC))
NeuronDendrites<-ROIwise.Dendrites

# histograms for neurons
ggplot(NeuronSomas, aes(x=meanAmp2)) + geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("neuronal soma mean amplitude- both stim types averaged") +
  max.theme

ggplot(NeuronDendrites, aes(x=meanAmp)) + geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("neuronal dendrites mean amplitude- both stim types averaged") +
  max.theme


# calculate responders
neuronS_percentiles<-quantile(NeuronSomas$meanAmp2, prob = seq(0, 1, length = 21), type = 5)
neuronD_percentiles<-quantile(NeuronDendrites$meanAmp, prob = seq(0, 1, length = 21), type = 5)

NeuronSomas$responders=0
NeuronSomas$responders[NeuronSomas$meanAmp2>neuronS_percentiles[20]]="high"
NeuronSomas$responders[NeuronSomas$meanAmp2<neuronS_percentiles[11]]="low"
NeuronSomas$responders[NeuronSomas$meanAmp2<=neuronS_percentiles[20] & NeuronSomas$meanAmp2>=neuronS_percentiles[11]]="mid"

NeuronDendrites$responders=0
NeuronDendrites$responders[NeuronDendrites$meanAmp>neuronD_percentiles[20]]="high"
NeuronDendrites$responders[NeuronDendrites$meanAmp<neuronD_percentiles[11]]="low"
NeuronDendrites$responders[NeuronDendrites$meanAmp<=neuronD_percentiles[20] & NeuronDendrites$meanAmp>=neuronD_percentiles[11]]="mid"

#dendrites vs somas

df.amp.somas<- summarySE(data=NeuronSomas, measurevar = "meanAmp2", groupvars = c("responders"))
df.amp.dend<- summarySE(data=NeuronDendrites, measurevar = "meanAmp", groupvars = c("responders"))


#Add the responder type info to the peak table

# enter responding ROI information- high, mid or low responding neurons
all.lck.peaks.group$Nresponders=0
all.lck.peaks.Nresp=data.frame()
for (ii in 1:nrow(NeuronDendrites))
{
  Dend=NeuronDendrites$ROIs_trial[ii]
  RespType=NeuronDendrites$responders[ii]
  subset1=subset(all.lck.peaks.group, ROIs_trial==Dend)
  if (nrow(subset1)>0)
  {
    subset1$Nresponders=RespType
    all.lck.peaks.Nresp<-rbind(all.lck.peaks.Nresp, subset1)
  }
}

all.lck.peaks.group$roiNameUnique<-paste(all.lck.peaks.group$Animal,all.lck.peaks.group$Spot,all.lck.peaks.group$roiName, sep="_")
#GCaMP_RCaMP.groups$Nresponders=0
for (ii in 1:nrow(NeuronSomas))
{
  Soma=NeuronSomas$roiNameUnique[ii]
  RespType=NeuronSomas$responders[ii]
  subset1=subset(all.lck.peaks.group, roiNameUnique==Soma)
  if (nrow(subset1)>0)
  {
    subset1$Nresponders=RespType
    subset1$roiNameUnique=NULL
    all.lck.peaks.Nresp<-rbind(all.lck.peaks.Nresp, subset1)
  }
}

all.lck.peaks.group$roiNameUnique=NULL
all.lck.peaks.Nresp<-rbind(all.lck.peaks.Nresp, all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",])

######

# subsets to look for example ROIs for movies

lck.fast<-subset(stim.lck.compdata.STIM, Channel_Group=="lck_GCaMP.fast")
lck.fast.corr<-subset(GCaMP_RCaMP.groups, GroupX=="fast A")
lck.fast.corr.close<-subset(lck.fast.corr, MinDistance<15)


lck.fast.corr.close<-merge(lck.fast.corr.close, lck.fast[, c("ROIs_trial", "area")], by="ROIs_trial", all.x=TRUE)


#sort by correlation
lck.fast.corr.close<-lck.fast.corr.close[order(lck.fast.corr.close$Short_Corr),] 





########

# histograms for onset time comparison RASTERs
AvsN.Space.Nostim <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/OnsetTimeComps/Nostim/Nostim_AvsN_SpaceOnsets.csv", header=TRUE, sep = ",")
AvsN.Time.Nostim <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/OnsetTimeComps/Nostim/Nostim_AvsN_TimeOnsets.csv", header=TRUE, sep = ",")
NvsA.Space.Nostim <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/OnsetTimeComps/Nostim/Nostim_NvsA_SpaceOnsets.csv", header=TRUE, sep = ",")
NvsA.Time.Nostim <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/OnsetTimeComps/Nostim/Nostim_NvsA_TimeOnsets.csv", header=TRUE, sep = ",")

AvsN.Space.Nostim <- read.table("D:/Data/GCaMP_RCaMP/Manuscript/Figures/OnsetTime_histograms/OnsetTimeComps/Nostim/Nostim_AvsN_SpaceOnsets.csv", header=TRUE, sep = ",")
AvsN.Time.Nostim <- read.table("D:/Data/GCaMP_RCaMP/Manuscript/Figures/OnsetTime_histograms/OnsetTimeComps/Nostim/Nostim_AvsN_TimeOnsets.csv", header=TRUE, sep = ",")
NvsA.Space.Nostim <- read.table("D:/Data/GCaMP_RCaMP/Manuscript/Figures/OnsetTime_histograms/OnsetTimeComps/Nostim/Nostim_NvsA_SpaceOnsets.csv", header=TRUE, sep = ",")
NvsA.Time.Nostim <- read.table("D:/Data/GCaMP_RCaMP/Manuscript/Figures/OnsetTime_histograms/OnsetTimeComps/Nostim/Nostim_NvsA_TimeOnsets.csv", header=TRUE, sep = ",")

AvsN.Space.Nostim<-AvsN.Space.Nostim[-1,]

AvsN.Space.Stim <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/OnsetTimeComps/Stim/Stim_AvsN_SpaceOnsets.csv", header=TRUE, sep = ",")
AvsN.Time.Stim <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/OnsetTimeComps/Stim/Stim_AvsN_TimeOnsets.csv", header=TRUE, sep = ",")
NvsA.Space.Stim <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/OnsetTimeComps/Stim/Stim_NvsA_SpaceOnsets.csv", header=TRUE, sep = ",")
NvsA.Time.Stim <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/OnsetTimeComps/Stim/Stim_NvsA_TimeOnsets.csv", header=TRUE, sep = ",")

AvsN.Space.Stim <- read.table("D:/Data/GCaMP_RCaMP/Manuscript/Figures/OnsetTime_histograms/OnsetTimeComps/Stim/Stim_AvsN_SpaceOnsets.csv", header=TRUE, sep = ",")
AvsN.Time.Stim <- read.table("D:/Data/GCaMP_RCaMP/Manuscript/Figures/OnsetTime_histograms/OnsetTimeComps/Stim/Stim_AvsN_TimeOnsets.csv", header=TRUE, sep = ",")
NvsA.Space.Stim <- read.table("D:/Data/GCaMP_RCaMP/Manuscript/Figures/OnsetTime_histograms/OnsetTimeComps/Stim/Stim_NvsA_SpaceOnsets.csv", header=TRUE, sep = ",")
NvsA.Time.Stim <- read.table("D:/Data/GCaMP_RCaMP/Manuscript/Figures/OnsetTime_histograms/OnsetTimeComps/Stim/Stim_NvsA_TimeOnsets.csv", header=TRUE, sep = ",")

AvsN.Space.Nostim$Group<-"other"
AvsN.Space.Nostim$AOnset<-as.numeric(as.character(AvsN.Space.Nostim$AOnset))
AvsN.Space.Nostim$TimeDiff<-as.numeric(as.character(AvsN.Space.Nostim$TimeDiff))
AvsN.Space.Nostim$Group[AvsN.Space.Nostim$AOnset<1]<-"fast"
AvsN.Space.Nostim$Group[(AvsN.Space.Nostim$AOnset>=1 & AvsN.Space.Nostim$AOnset<12)]<-"delayed"

AvsN.Time.Nostim$Group<-"other"
AvsN.Time.Nostim$AOnset<-as.numeric(as.character(AvsN.Time.Nostim$AOnset))
AvsN.Time.Nostim$TimeDiff<-as.numeric(as.character(AvsN.Time.Nostim$TimeDiff))
AvsN.Time.Nostim$Group[AvsN.Time.Nostim$AOnset<1]<-"fast"
AvsN.Time.Nostim$Group[(AvsN.Time.Nostim$AOnset>=1 & AvsN.Time.Nostim$AOnset<12)]<-"delayed"


NvsA.Space.Nostim$Group<-"other"
NvsA.Space.Nostim$AOnset<-as.numeric(as.character(NvsA.Space.Nostim$AOnset))
NvsA.Space.Nostim$TimeDiff<-as.numeric(as.character(NvsA.Space.Nostim$TimeDiff))
NvsA.Space.Nostim$Group[NvsA.Space.Nostim$AOnset<1]<-"fast"
NvsA.Space.Nostim$Group[(NvsA.Space.Nostim$AOnset>=1 & NvsA.Space.Nostim$AOnset<12)]<-"delayed"


NvsA.Time.Nostim$Group<-"other"
NvsA.Time.Nostim$AOnset<-as.numeric(as.character(NvsA.Time.Nostim$AOnset))
NvsA.Time.Nostim$TimeDiff<-as.numeric(as.character(NvsA.Time.Nostim$TimeDiff))
NvsA.Time.Nostim$Group[NvsA.Time.Nostim$AOnset<1]<-"fast"
NvsA.Time.Nostim$Group[(NvsA.Time.Nostim$AOnset>=1 & NvsA.Time.Nostim$AOnset<12)]<-"delayed"


AvsN.Space.Stim$Group<-"other"
AvsN.Space.Stim$AOnset<-as.numeric(as.character(AvsN.Space.Stim$AOnset))
AvsN.Space.Stim$TimeDiff<-as.numeric(as.character(AvsN.Space.Stim$TimeDiff))
AvsN.Space.Stim$Group[AvsN.Space.Stim$AOnset<1]<-"fast"
AvsN.Space.Stim$Group[(AvsN.Space.Stim$AOnset>=1 & AvsN.Space.Stim$AOnset<12)]<-"delayed"

AvsN.Time.Stim$Group<-"other"
AvsN.Time.Stim$AOnset<-as.numeric(as.character(AvsN.Time.Stim$AOnset))
AvsN.Time.Stim$TimeDiff<-as.numeric(as.character(AvsN.Time.Stim$TimeDiff))
AvsN.Time.Stim$Group[AvsN.Time.Stim$AOnset<1]<-"fast"
AvsN.Time.Stim$Group[(AvsN.Time.Stim$AOnset>=1 & AvsN.Time.Stim$AOnset<12)]<-"delayed"


NvsA.Space.Stim$Group<-"other"
NvsA.Space.Stim$AOnset<-as.numeric(as.character(NvsA.Space.Stim$AOnset))
NvsA.Space.Stim$TimeDiff<-as.numeric(as.character(NvsA.Space.Stim$TimeDiff))
NvsA.Space.Stim$Group[NvsA.Space.Stim$AOnset<1]<-"fast"
NvsA.Space.Stim$Group[(NvsA.Space.Stim$AOnset>=1 & NvsA.Space.Stim$AOnset<12)]<-"delayed"


NvsA.Time.Stim$Group<-"other"
NvsA.Time.Stim$AOnset<-as.numeric(as.character(NvsA.Time.Stim$AOnset))
NvsA.Time.Stim$TimeDiff<-as.numeric(as.character(NvsA.Time.Stim$TimeDiff))
NvsA.Time.Stim$Group[NvsA.Time.Stim$AOnset<1]<-"fast"
NvsA.Time.Stim$Group[(NvsA.Time.Stim$AOnset>=1 & NvsA.Time.Stim$AOnset<12)]<-"delayed"


# histograms
ggplot(AvsN.Space.Nostim[AvsN.Space.Nostim$Group!="other",], aes(x=TimeDiff, fill=Group)) + #y=..density..,
  geom_histogram(aes(y = ..density..),binwidth=1, position="dodge")+
  #geom_density(aes(y = ..density..*(164*0.1)))+
  ggtitle("AvsN Space Nostim")+
  ylim(0,0.4)+
  xlim(-15,15)+
  max.theme

ggplot(AvsN.Space.Stim[AvsN.Space.Stim$Group!="other",], aes(x=TimeDiff, fill=Group)) + #y=..density..,
  geom_histogram(aes(y = ..density..),binwidth=1, position="dodge")+
  ggtitle("AvsN Space Stim")+
  ylim(0,0.4)+
  xlim(-15,15)+
  max.theme

# no colour
ggplot(AvsN.Space.Nostim[AvsN.Space.Nostim$Group!="other",], aes(x=TimeDiff, fill="blue")) + #y=..density..,
  geom_histogram(binwidth=1, position="dodge")+
  #geom_density(aes(y = ..density..*(164*0.1)))+
  ggtitle("AvsN Space Nostim")+
  ylim(0,350)+
  xlim(-12,12)+
  max.theme

ggplot(AvsN.Space.Stim[AvsN.Space.Stim$Group!="other",], aes(x=TimeDiff, fill="blue")) + #y=..density..,
  geom_histogram(binwidth=1, position="dodge")+
  ggtitle("AvsN Space Stim")+
  ylim(0,350)+
  xlim(-12,12)+
  max.theme


# compare distributions
AvsN.Space.kstest<- ks.test(AvsN.Space.Stim$TimeDiff[AvsN.Space.Stim$Group!="other"],
                            AvsN.Space.Nostim$TimeDiff[AvsN.Space.Nostim$Group!="other"])
print(AvsN.Space.kstest)



ggplot(AvsN.Time.Nostim[AvsN.Time.Nostim$Group!="other",], aes(x=TimeDiff, fill=Group)) + #y=..density..,
  geom_histogram(aes(y = ..density..),binwidth=1, position="dodge")+
  ggtitle("AvsN Time Nostim")+
  ylim(0,0.8)+
  xlim(-15,15)+
  max.theme

ggplot(AvsN.Time.Stim[AvsN.Time.Stim$Group!="other",], aes(x=TimeDiff, fill=Group)) + #y=..density..,
  geom_histogram(aes(y = ..density..),binwidth=1, position="dodge")+
  ggtitle("AvsN Time Stim")+
  ylim(0,0.8)+
  xlim(-15,15)+
  max.theme

ggplot(NvsA.Time.Nostim[NvsA.Time.Nostim$Group!="other",], aes(x=TimeDiff, fill=Group)) + #y=..density..,
  geom_histogram(aes(y = ..density..),binwidth=1, position="dodge")+
  ggtitle("NvsA Time Nostim")+
  ylim(0,0.7)+
  xlim(-15,15)+
  max.theme

ggplot(NvsA.Time.Stim[NvsA.Time.Stim$Group!="other",], aes(x=TimeDiff, fill=Group)) + #y=..density..,
  geom_histogram(aes(y = ..density..),binwidth=1, position="dodge")+
  ggtitle("NvsA Time Stim")+
  ylim(0,0.7)+
  xlim(-15,15)+
  max.theme

ggplot(NvsA.Space.Nostim[NvsA.Space.Nostim$Group!="other",], aes(x=TimeDiff, fill=Group)) + #y=..density..,
  geom_histogram(aes(y = ..density..),binwidth=1, position="dodge")+
  ggtitle("NvsA Space Nostim")+
  ylim(0,0.4)+
  xlim(-15,15)+
  max.theme

ggplot(NvsA.Space.Stim[NvsA.Space.Stim$Group!="other",], aes(x=TimeDiff, fill=Group)) + #y=..density..,
  geom_histogram(aes(y = ..density..),binwidth=1, position="dodge")+
  ggtitle("NvsA Space Stim")+
  ylim(0,0.4)+
  xlim(-15,15)+
  max.theme






###### 
# proportion of ROIs that are fast or delayed

#what about proportion of each group based on ROI type??

GCaMP<-subset(all.lck.peaks.Nresp,Channel=="GCaMP")
RCaMP<-subset(all.lck.peaks.Nresp,Channel=="RCaMP")

GCaMP$Group<-factor(GCaMP$Group,levels= c("fast","delayed"))

fastROIs<-unique(GCaMP$ROIs_trial[GCaMP$Group=="fast"])

allROIs.all<-length(unique(GCaMP$ROIs_trial))
fastROIs.all<-length(unique(GCaMP$ROIs_trial[GCaMP$Group=="fast"]))
delayedROIs.all<-length(unique(GCaMP$ROIs_trial[GCaMP$Group=="delayed"]))

longstim.GC<-subset(GCaMP, Condition=="Stim")
allROIs.L<-length(unique(longstim.GC$ROIs_trial))
fastROIs.L<-length(unique(longstim.GC$ROIs_trial[longstim.GC$Group=="fast"]))
delayedROIs.L<-length(unique(longstim.GC$ROIs_trial[longstim.GC$Group=="delayed"]))

shortstim.GC<-subset(GCaMP, Condition=="shortstim")
allROIs.S<-length(unique(shortstim.GC$ROIs_trial))
fastROIs.S<-length(unique(shortstim.GC$ROIs_trial[shortstim.GC$Group=="fast"]))
delayedROIs.S<-length(unique(shortstim.GC$ROIs_trial[shortstim.GC$Group=="delayed"]))

#proportion
propfast.all=fastROIs.all/allROIs.all
propdelayed.all=delayedROIs.all/allROIs.all

propfast.L=fastROIs.L/allROIs.L
propdelayed.L=delayedROIs.L/allROIs.L

propfast.S=fastROIs.S/allROIs.S
propdelayed.S=delayedROIs.S/allROIs.S


#ef vs proc
allROIs.ef<-length(unique(GCaMP$ROIs_trial[GCaMP$ROIType=="Endfoot"]))
fastROIs.ef<-length(unique(GCaMP$ROIs_trial[GCaMP$Group=="fast"&GCaMP$ROIType=="Endfoot"]))
delayedROIs.ef<-length(unique(GCaMP$ROIs_trial[GCaMP$Group=="delayed"&GCaMP$ROIType=="Endfoot"]))

allROIs.proc<-length(unique(GCaMP$ROIs_trial[GCaMP$ROIType=="Process"]))
fastROIs.proc<-length(unique(GCaMP$ROIs_trial[GCaMP$Group=="fast"&GCaMP$ROIType=="Process"]))
delayedROIs.proc<-length(unique(GCaMP$ROIs_trial[GCaMP$Group=="delayed"&GCaMP$ROIType=="Process"]))


propfast.ef=fastROIs.ef/allROIs.ef
propdelayed.ef=delayedROIs.ef/allROIs.ef

propfast.proc=fastROIs.proc/allROIs.proc
propdelayed.proc=delayedROIs.proc/allROIs.proc




#################
# conditional probabilty of astrocyte response given a nearby neuronal response


# find total number of astrocyte ROIs or neuronal ROIs in a field of view
# find number of responding astrocytes or responding neuronal per trial
# find the number of fast astrocytes per field of view


# number of ROIs in each trial for each field of view (across the time window (2 s for neurons, 15 s for AC))

stim.lck.OT.window2$Channel <- factor(stim.lck.OT.window2$Channel, levels = c("RCaMP","GCaMP"))

ROInum.lck.stim<-ddply(stim.lck.OT.window2, c("Animal","Spot","Condition","Channel"), summarise, nROIs=length(OnsetTime))


# number of ROIs

Spot.TotalROIs<-merge(Spot.TotalROIs2, Spot.noPeaks[, c("trials_Cond", "nROIs")], by="trials_Cond", all.x=TRUE)



#remove peaks with zero peakAUC

count.peaks <- ddply (post30.sig.peaks, c("Animal", "Spot", "Layer", "Condition","ROI","ROIType","peakType"), summarise,
                      Peaks = length(numPeaks))

#divide by the number of trials for each ROI
count.peaks$peakPertrial <-0

#Nostim
#NS.peaks <- subset(count.peaks,Condition=="Nostim")
#NS.ROI.p30 <- subset(post30.ROI,Condition=="Nostim")
sigroiNamesNS <-as.character(unique(post30.ROI$ROI))
count.peaks2 <-data.frame()
for (ii in 1:length(sigroiNamesNS))
{
  name =sigroiNamesNS[ii]
  subset1 = subset(post30.ROI, ROI == name)
  subset2 = subset(count.peaks, ROI == name)
  for (iii in 1:length(subset2$ROI))
  {
    #subset2$peakPertrial<- subset2$Peaks[iii]/subset1$N    
    subset2$Trial<- subset1$N[1]   
  }
  count.peaks2<- rbind(count.peaks2,subset2)
}


# if no peak of a particular type exists for individual ROI, it becomes zero
ROI = as.character(unique(count.peaks$ROI))
count.peaks3 <- data.frame()

for (ii in 1:length(ROI))
{
  name = ROI[ii]
  ROIubset.NS = subset(count.peaks2 , ROI == name & Condition == "Nostim")
  ROIubset.S = subset(count.peaks2 , ROI == name & Condition == "Stim")
  # use data from stim if there are no peaks during no stim
  if ((nrow(ROIubset.NS) ==0)==TRUE)
  {ROIubset.NS<- head(ROIubset.S,1)
  ROIubset.NS$Condition ="Nostim"
  ROIubset.NS$peakType ="NaN"
  }
  # fill in zero ROI for each type
  if ((nrow(ROIubset.NS) ==3)==TRUE)
  {count.peaks3<- rbind(count.peaks3,ROIubset.NS)
  } else {
    # find peakTypes with matches
    single.NS = grepl("Singlepeak",ROIubset.NS$peakType)
    multi.NS = grepl("Multipeak",ROIubset.NS$peakType)
    plateau.NS = grepl("Plateau",ROIubset.NS$peakType)
    
    if (nrow(ROIubset.NS[single.NS,])>0)
    { count.peaks3<- rbind(count.peaks3,ROIubset.NS[single.NS,])
    } else {
      zeros <- head(ROIubset.NS,1)
      zeros$peakType = "Singlepeak"
      zeros$Peaks = 0
      zeros$peakPertrial = 0
      count.peaks3<- rbind(count.peaks3,zeros)  
    }
    if (nrow(ROIubset.NS[multi.NS,])>0)
    { count.peaks3<- rbind(count.peaks3,ROIubset.NS[multi.NS,])
    } else {
      zeros <- head(ROIubset.NS,1)
      zeros$peakType = "Multipeak"
      zeros$Peaks = 0
      zeros$peakPertrial = 0
      count.peaks3<- rbind(count.peaks3,zeros)  
    }
    if (nrow(ROIubset.NS[plateau.NS,])>0)
    { count.peaks3<- rbind(count.peaks3,ROIubset.NS[plateau.NS,])
    } else {
      zeros <- head(ROIubset.NS,1)
      zeros$peakType = "Plateau"
      zeros$Peaks = 0
      zeros$peakPertrial = 0
      count.peaks3<- rbind(count.peaks3,zeros)  
    }  
  }
  
  if ((nrow(ROIubset.S) ==3)==TRUE)
  {count.peaks3<- rbind(count.peaks3,ROIubset.S)
  } else {
    # find peakTypes with matches
    single.S = grepl("Singlepeak",ROIubset.S$peakType)
    multi.S = grepl("Multipeak",ROIubset.S$peakType)
    plateau.S = grepl("Plateau",ROIubset.S$peakType)
    
    if (nrow(ROIubset.S[single.S,])>0)
    { count.peaks3<- rbind(count.peaks3,ROIubset.S[single.S,])
    } else {
      zeros <- head(ROIubset.S,1)
      zeros$peakType = "Singlepeak"
      zeros$Peaks = 0
      zeros$peakPertrial = 0
      count.peaks3<- rbind(count.peaks3,zeros)  
    }
    if (nrow(ROIubset.S[multi.S,])>0)
    { count.peaks3<- rbind(count.peaks3,ROIubset.S[multi.S,])
    } else {
      zeros <- head(ROIubset.S,1)
      zeros$peakType = "Multipeak"
      zeros$Peaks = 0
      zeros$peakPertrial = 0
      count.peaks3<- rbind(count.peaks3,zeros)  
    }
    if (nrow(ROIubset.S[plateau.S,])>0)
    { count.peaks3<- rbind(count.peaks3,ROIubset.S[plateau.S,])
    } else {
      zeros <- head(ROIubset.S,1)
      zeros$peakType = "Plateau"
      zeros$Peaks = 0
      zeros$peakPertrial = 0
      count.peaks3<- rbind(count.peaks3,zeros)  
    }  
  }
}


count.peaks3$peakPertrial<- count.peaks3$Peaks/count.peaks3$Trial 
