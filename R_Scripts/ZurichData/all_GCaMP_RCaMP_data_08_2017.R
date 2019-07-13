library("lme4")
library("lmerTest")
library("plyr")
library("dplyr")
library("ggplot2")
library("lsmeans")
library("Rmisc")
library("multcomp")
library("reshape2")
library("tidyr")
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

long.lck <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_longstim_28_04_2017.csv", header=TRUE, sep = ",")
short.lck <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_shortstim_28_04_2017.csv", header=TRUE, sep = ",")
nostim.lck <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_nostim_28_04_2017.csv", header=TRUE, sep = ",")

long.cyto<- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cytoGC&RC_2D_longstim_28_04_2017.csv", header=TRUE, sep = ",")
short.cyto <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cytoGC&RC_2D_shortstim_28_04_2017.csv", header=TRUE, sep = ",")
nostim.cyto <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cytoGC&RC_2D_nostim_28_04_2017.csv", header=TRUE, sep = ",")

long.cyto.DSP4 <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cytoGC&RC_2D_longstim_28_04_2017.csv", header=TRUE, sep = ",")
short.cyto.DSP4<- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cytoGC&RC_2D_shortstim_28_04_2017.csv", header=TRUE, sep = ",")
nostim.cyto.DSP4<- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cytoGC&RC_2D_nostim_28_04_2017.csv", header=TRUE, sep = ",")

# onset time comparisons for nostim data
#longstim.OT.comp <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/longstim_firstonset_comparisons.csv", header=TRUE, sep = ",")
#shortstim.OT.comp <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/shortstim_firstonset_comparisons.csv", header=TRUE, sep = ",")
#nostim.OT.comp <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/nostim_firstonset_comparisons.csv", header=TRUE, sep = ",")

longstim.lck.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/Lck_longstim_onset&AUC.csv", header=TRUE, sep = ",")
shortstim.lck.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/shortstim_onset&AUC.csv", header=TRUE, sep = ",")
nostim.lck.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/Lck_nostim_onset&AUC.csv", header=TRUE, sep = ",")

longstim.lck.OT.2s <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/Lck_longstim_onset_2sbeforeStim.csv", header=TRUE, sep = ",")
nostim.lck.OT.2s <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/Lck_nostim_onset_2sbeforeStim.csv", header=TRUE, sep = ",")
longstim.cyto.OT.2s <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_longstim_onset_2sbeforeStim.csv", header=TRUE, sep = ",")
nostim.cyto.OT.2s <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_nostim_onset_2sbeforeStim.csv", header=TRUE, sep = ",")

longstim.cyto.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_longstim_onset&AUC.csv", header=TRUE, sep = ",")
shortstim.cyto.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_shortstim_onset&AUC.csv", header=TRUE, sep = ",")
nostim.cyto.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_nostim_firstonset&AUC.csv", header=TRUE, sep = ",")

longstim.cyto.DSP4.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cyto_longstim_firstonset&AUC.csv", header=TRUE, sep = ",")
shortstim.cyto.DSP4.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cyto_shortstim_firstonset&AUC.csv", header=TRUE, sep = ",")
nostim.cyto.DSP4.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cyto_nostim_firstonset&AUC.csv", header=TRUE, sep = ",")

##### 
#home files
long.lck <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_longstim_28_04_2017.csv", header=TRUE, sep = ",")
short.lck <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_shortstim_28_04_2017.csv", header=TRUE, sep = ",")
nostim.lck <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_nostim_28_04_2017.csv", header=TRUE, sep = ",")

long.cyto<- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cytoGC&RC_2D_longstim_28_04_2017.csv", header=TRUE, sep = ",")
short.cyto <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cytoGC&RC_2D_shortstim_28_04_2017.csv", header=TRUE, sep = ",")
nostim.cyto <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cytoGC&RC_2D_nostim_28_04_2017.csv", header=TRUE, sep = ",")

long.cyto.DSP4 <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cytoGC&RC_2D_longstim_28_04_2017.csv", header=TRUE, sep = ",")
short.cyto.DSP4<- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cytoGC&RC_2D_shortstim_28_04_2017.csv", header=TRUE, sep = ",")
nostim.cyto.DSP4<- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cytoGC&RC_2D_nostim_28_04_2017.csv", header=TRUE, sep = ",")

longstim.lck.OT <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/longstim_onset&AUC.csv", header=TRUE, sep = ",")
shortstim.lck.OT <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/shortstim_onset&AUC.csv", header=TRUE, sep = ",")
nostim.lck.OT <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/nostim_onset&AUC.csv", header=TRUE, sep = ",")

longstim.cyto.OT <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_longstim_onset&AUC.csv", header=TRUE, sep = ",")
shortstim.cyto.OT <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_shortstim_onset&AUC.csv", header=TRUE, sep = ",")
nostim.cyto.OT <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_nostim_firstonset&AUC.csv", header=TRUE, sep = ",")

longstim.cyto.DSP4.OT <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cyto_longstim_firstonset&AUC.csv", header=TRUE, sep = ",")
shortstim.cyto.DSP4.OT <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cyto_shortstim_firstonset&AUC.csv", header=TRUE, sep = ",")
nostim.cyto.DSP4.OT <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cyto_nostim_firstonset&AUC.csv", header=TRUE, sep = ",")


##########

lsm.options(pbkrtest.limit = 100000)

all.lck.peaks<-rbind(nostim.lck,long.lck,short.lck)
all.lck.OT<-rbind(nostim.lck.OT,longstim.lck.OT,shortstim.lck.OT)
all.lck.OT.2s<-rbind(nostim.lck.OT.2s,longstim.lck.OT.2s)
all.lck.OT$Spot_trial_Cond<-paste(all.lck.OT$Spot_trial, all.lck.OT$Condition, sep="_")
all.lck.OT$ROIs_Cond<-paste(all.lck.OT$ROIs_trial, all.lck.OT$Condition, sep="_")
all.lck.OT.2s$Spot_trial_Cond<-paste(all.lck.OT.2s$Spot_trial, all.lck.OT.2s$Condition, sep="_")
all.lck.OT.2s$ROIs_Cond<-paste(all.lck.OT.2s$ROIs_trial, all.lck.OT.2s$Condition, sep="_")



# remove trials with no ROIs
all.lck.peaks<-all.lck.peaks[!(all.lck.peaks$ROIname=="none"),]


all.lck.peaks$ROIType= 0
all.lck.peaksA<- subset(all.lck.peaks, Channel=="GCaMP")
all.lck.peaksB<- subset(all.lck.peaks, Channel=="RCaMP")

# ROITypes
all.lck.peaksA$ROIType[grepl("r",all.lck.peaksA$ROIname)]="Process"
all.lck.peaksA$ROIType[grepl("E",all.lck.peaksA$ROIname)]="Endfoot"
all.lck.peaksB$ROIType[grepl("r",all.lck.peaksB$ROIname)]="Dendrite"
all.lck.peaksB$ROIType[grepl("D",all.lck.peaksB$ROIname)]="Dendrite"
all.lck.peaksB$ROIType[grepl("N",all.lck.peaksB$ROIname)]="Neuron"

all.lck.peaks<-rbind(all.lck.peaksA, all.lck.peaksB)
all.lck.peaks$ROIType<- as.factor(all.lck.peaks$ROIType)

#unique ROI names
all.lck.peaks$ROIs_trial<-paste(all.lck.peaks$Animal, all.lck.peaks$Spot, all.lck.peaks$Trial,all.lck.peaks$ROIname, sep= "_")
all.lck.peaks$trials<-paste(all.lck.peaks$Animal, all.lck.peaks$Spot, all.lck.peaks$Trial, sep= "_")
all.lck.peaks$trials_Cond<-paste(all.lck.peaks$trials, all.lck.peaks$Condition, sep= "_")

all.lck.peaks$ROIs_Cond<-paste(all.lck.peaks$ROIs_trial, all.lck.peaks$Condition, sep= "_")


# count number of trials per spot
Spot.lck.ntrials<-ddply(all.lck.peaks, c("Animal","Spot","Condition"), summarise, nTrials=length(unique(Trial)))
Spot.lck.ntrials$Ani_Spot_Cond<-paste(Spot.lck.ntrials$Animal, Spot.lck.ntrials$Spot, Spot.lck.ntrials$Condition, sep="_")
                     
# exclude the neuropil ROIs, because they were hand selected and not necessary
all.lck.peaks<-all.lck.peaks[!(all.lck.peaks$ROIname=="np"),]

# remove matching astrocyte process and soma ROIs
Overlap= all.lck.peaks$overlap!=0
all.lck.peaks<-all.lck.peaks[!Overlap,]
#OverlapROIs<-unique(nostim$ROIs_trial[Overlap])

# adjust peak time and duration
all.lck.peaks$peakTime<- all.lck.peaks$peakTime-5
all.lck.peaks$peakStart<- all.lck.peaks$peakStart-5
all.lck.peaks$peakStartHalf<- all.lck.peaks$peakStartHalf-5
all.lck.peaks$Duration<- all.lck.peaks$halfWidth*2

# drop peaks that occur before the start of stimulation
all.lck.peaks2<-subset(all.lck.peaks,peakTime>0)

# REMOVE duplicate entries from onset time and peak time data frames 
# only the first entry will be used
all.lck.OT2=all.lck.OT[order(all.lck.OT$OnsetTime),] # sort by ascending onset time
all.lck.OT2.2s=all.lck.OT.2s[order(all.lck.OT.2s$OnsetTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.lck.OT<-distinct(all.lck.OT2, ROIs_Cond, .keep_all = TRUE)
all.lck.OT.2s<-distinct(all.lck.OT2.2s, ROIs_Cond, .keep_all = TRUE)

# only the first entry will be used
all.lck.peaks3<-all.lck.peaks2[order(all.lck.peaks2$peakTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.lck.peaks<-distinct(all.lck.peaks3, ROIs_Cond,.keep_all = TRUE)


#adjust onset time for the data with 2 s before stimulation included:

all.lck.OT.2s$OnsetTimeAdjust<-all.lck.OT.2s$OnsetTime-2


##### 
# cytoGCaMP6s data

all.cyto.peaks<-rbind(nostim.cyto,long.cyto,short.cyto)
all.cyto.OT<-rbind(nostim.cyto.OT,longstim.cyto.OT,shortstim.cyto.OT)
all.cyto.OT$Spot_trial_Cond<-paste(all.cyto.OT$Spot_trial, all.cyto.OT$Condition, sep="_")
all.cyto.OT$ROIs_Cond<-paste(all.cyto.OT$ROIs_trial, all.cyto.OT$Condition, sep="_")

all.cyto.OT.2s<-rbind(nostim.cyto.OT.2s,longstim.cyto.OT.2s)
all.cyto.OT.2s$Spot_trial_Cond<-paste(all.cyto.OT.2s$Spot_trial, all.cyto.OT.2s$Condition, sep="_")
all.cyto.OT.2s$ROIs_Cond<-paste(all.cyto.OT.2s$ROIs_trial, all.cyto.OT.2s$Condition, sep="_")


# exclude the neuropil ROIs, because they were hand selected and not necessary
all.cyto.peaks<-all.cyto.peaks[!(all.cyto.peaks$ROIname=="np"),]
all.cyto.peaks<-all.cyto.peaks[!(all.cyto.peaks$ROIname=="none"),]

# no stim peak data
all.cyto.peaks$ROIType= 0
all.cyto.peaksA<- subset(all.cyto.peaks, Channel=="GCaMP")
all.cyto.peaksB<- subset(all.cyto.peaks, Channel=="RCaMP")

# ROITypes
all.cyto.peaksA$ROIType[grepl("r",all.cyto.peaksA$ROIname)]="Process"
all.cyto.peaksA$ROIType[grepl("E",all.cyto.peaksA$ROIname)]="Endfoot"
all.cyto.peaksA$ROIType[grepl("S",all.cyto.peaksA$ROIname)]="Soma"
all.cyto.peaksB$ROIType[grepl("r",all.cyto.peaksB$ROIname)]="Dendrite"
all.cyto.peaksB$ROIType[grepl("D",all.cyto.peaksB$ROIname)]="Dendrite"
all.cyto.peaksB$ROIType[grepl("N",all.cyto.peaksB$ROIname)]="Neuron"

all.cyto.peaks<-rbind(all.cyto.peaksA, all.cyto.peaksB)
all.cyto.peaks$ROIType<- as.factor(all.cyto.peaks$ROIType)

#unique ROI names
all.cyto.peaks$ROIs_trial<-paste(all.cyto.peaks$Animal, all.cyto.peaks$Spot, all.cyto.peaks$Trial,all.cyto.peaks$ROIname, sep= "_")
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
all.cyto.peaks$peakTime<- all.cyto.peaks$peakTime-5
all.cyto.peaks$peakStart<- all.cyto.peaks$peakStart-5
all.cyto.peaks$peakStartHalf<- all.cyto.peaks$peakStartHalf-5
all.cyto.peaks$Duration<- all.cyto.peaks$halfWidth*2


# drop peaks that occur before the start of stimulation
all.cyto.peaks2<-subset(all.cyto.peaks,peakTime>0)

# REMOVE duplicate entries from onset time and peak time data frames 
# only the first entry will be used
all.cyto.OT2=all.cyto.OT[order(all.cyto.OT$OnsetTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.cyto.OT<-distinct(all.cyto.OT2, ROIs_Cond,.keep_all = TRUE)

# only the first entry will be used
all.cyto.peaks3<-all.cyto.peaks2[order(all.cyto.peaks2$peakTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.cyto.peaks<-distinct(all.cyto.peaks3, ROIs_Cond,.keep_all = TRUE)


# REMOVE duplicate entries from onset time and peak time data frames 
# only the first entry will be used
all.cyto.OT2=all.cyto.OT[order(all.cyto.OT$OnsetTime),] # sort by ascending onset time
all.cyto.OT2.2s=all.cyto.OT.2s[order(all.cyto.OT.2s$OnsetTime),] # sort by ascending onset time

#adjust onset time for the data with 2 s before stimulation included:

all.cyto.OT.2s$OnsetTimeAdjust<-all.cyto.OT.2s$OnsetTime-2

########
# cytoGCaMP6s DSP4 data

all.cyto.DSP4.peaks<-rbind(nostim.cyto.DSP4,long.cyto.DSP4,short.cyto.DSP4)
all.cyto.DSP4.OT<-rbind(nostim.cyto.DSP4.OT,longstim.cyto.DSP4.OT,shortstim.cyto.DSP4.OT)
all.cyto.DSP4.OT$Spot_trial_Cond<-paste(all.cyto.DSP4.OT$Spot_trial, all.cyto.DSP4.OT$Condition, sep="_")
all.cyto.DSP4.OT$ROIs_Cond<-paste(all.cyto.DSP4.OT$ROIs_trial, all.cyto.DSP4.OT$Condition, sep="_")

# exclude the neuropil ROIs, because they were hand selected and not necessary
all.cyto.DSP4.peaks<-all.cyto.DSP4.peaks[!(all.cyto.DSP4.peaks$ROIname=="np"),]
all.cyto.DSP4.peaks<-all.cyto.DSP4.peaks[!(all.cyto.DSP4.peaks$ROIname=="none"),]

# no stim peak data
all.cyto.DSP4.peaks$ROIType= 0
all.cyto.DSP4.peaksA<- subset(all.cyto.DSP4.peaks, Channel=="GCaMP")
all.cyto.DSP4.peaksB<- subset(all.cyto.DSP4.peaks, Channel=="RCaMP")

# ROITypes
all.cyto.DSP4.peaksA$ROIType[grepl("r",all.cyto.DSP4.peaksA$ROIname)]="Process"
all.cyto.DSP4.peaksA$ROIType[grepl("E",all.cyto.DSP4.peaksA$ROIname)]="Endfoot"
all.cyto.DSP4.peaksA$ROIType[grepl("S",all.cyto.DSP4.peaksA$ROIname)]="Soma"
all.cyto.DSP4.peaksB$ROIType[grepl("r",all.cyto.DSP4.peaksB$ROIname)]="Dendrite"
all.cyto.DSP4.peaksB$ROIType[grepl("D",all.cyto.DSP4.peaksB$ROIname)]="Dendrite"
all.cyto.DSP4.peaksB$ROIType[grepl("N",all.cyto.DSP4.peaksB$ROIname)]="Neuron"

all.cyto.DSP4.peaks<-rbind(all.cyto.DSP4.peaksA, all.cyto.DSP4.peaksB)
all.cyto.DSP4.peaks$ROIType<- as.factor(all.cyto.DSP4.peaks$ROIType)

#unique ROI names
all.cyto.DSP4.peaks$ROIs_trial<-paste(all.cyto.DSP4.peaks$Animal, all.cyto.DSP4.peaks$Spot, all.cyto.DSP4.peaks$Trial,all.cyto.DSP4.peaks$ROIname, sep= "_")
all.cyto.DSP4.peaks$trials<-paste(all.cyto.DSP4.peaks$Animal, all.cyto.DSP4.peaks$Spot, all.cyto.DSP4.peaks$Trial, sep= "_")
all.cyto.DSP4.peaks$trials_Cond<-paste(all.cyto.DSP4.peaks$trials, all.cyto.DSP4.peaks$Condition, sep= "_")
all.cyto.DSP4.peaks$ROIs_Cond<-paste(all.cyto.DSP4.peaks$ROIs_trial, all.cyto.DSP4.peaks$Condition, sep= "_")

# remove matching astrocyte process and soma ROIs
Overlap= all.cyto.DSP4.peaks$overlap!=0
all.cyto.DSP4.peaks<-all.cyto.DSP4.peaks[!Overlap,]
#OverlapROIs<-unique(nostim$ROIs_trial[Overlap])

# adjust peak time and duration
all.cyto.DSP4.peaks$peakTime<- all.cyto.DSP4.peaks$peakTime-5
all.cyto.DSP4.peaks$peakStart<- all.cyto.DSP4.peaks$peakStart-5
all.cyto.DSP4.peaks$peakStartHalf<- all.cyto.DSP4.peaks$peakStartHalf-5
all.cyto.DSP4.peaks$Duration<- all.cyto.DSP4.peaks$halfWidth*2


# drop peaks that occur before the start of stimulation
all.cyto.DSP4.peaks2<-subset(all.cyto.DSP4.peaks,peakTime>0)

# REMOVE duplicate entries from onset time and peak time data frames 
# only the first entry will be used
all.cyto.DSP4.OT2=all.cytoDSP4..OT[order(all.cyto.DSP4.OT$OnsetTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.cyto.DSP4.OT<-distinct(all.cyto.DSP4.OT2, ROIs_Cond,.keep_all = TRUE)

# only the first entry will be used
all.cyto.DSP4.peaks3<-all.cyto.DSP4.peaks2[order(all.cyto.DSP4.peaks2$peakTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.cyto.DSP4.peaks<-distinct(all.cyto.DSP4.peaks3, ROIs_Cond,.keep_all = TRUE)


######

# distributions of onset times and peak times

# both Lck and cyto GCaMP

# all signals within 15 s of stimulus start

stimwindow=15

stim.lck.OT.dist<-subset(all.lck.OT,Condition!="shortstim" & OnsetTime<stimwindow)
short.lck.OT.dist<-subset(all.lck.OT,Condition!="Stim" & OnsetTime<stimwindow)
stim.lck.OT.dist.2s<-subset(all.lck.OT.2s,Condition!="shortstim" & OnsetTime<stimwindow)

stim.cyto.OT.dist<-subset(all.cyto.OT,Condition!="shortstim" & OnsetTime<stimwindow)
short.cyto.OT.dist<-subset(all.cyto.OT,Condition!="Stim" & OnsetTime<stimwindow)
stim.cyto.OT.dist.2s<-subset(all.cyto.OT.2s,Condition!="shortstim" & OnsetTime<stimwindow)

stim.cyto.DSP4.OT.dist<-subset(all.cyto.DSP4.OT,Condition!="shortstim" & OnsetTime<stimwindow)
short.cyto.DSP4.OT.dist<-subset(all.cyto.DSP4.OT,Condition!="Stim" & OnsetTime<stimwindow)

#########
# LCK DATA Onset time histograms- normalized to the number of trials

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


##

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

#####################
#CYTO data onset times
histseq= seq(0,15,0.5)

ntrials.cyto.OT.long.R.dis<- ddply(stim.cyto.OT.dist[stim.cyto.OT.dist$Channel=="RCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))
ntrials.cyto.OT.long.G.dis<- ddply(stim.cyto.OT.dist[stim.cyto.OT.dist$Channel=="GCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))

ntrials.cyto.DSP4.OT.long.R<- ddply(stim.cyto.DSP4.OT.dist[stim.cyto.DSP4.OT.dist$Channel=="RCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))
ntrials.cyto.DSP4.OT.long.G<- ddply(stim.cyto.DSP4.OT.dist[stim.cyto.DSP4.OT.dist$Channel=="GCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))

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

##
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



#######
#LCK data peak time distributions

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


##

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

##########
#CYTO data peak times

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


##

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




#######
# number of ROIs per field of view (in 15s time window after stim start- for both neurons and astrocytes)

# use data from onset time distributions


# LCK data
stim.lck.OT.dist$Channel <- factor(stim.lck.OT.dist$Channel, levels = c("RCaMP","GCaMP"))

ROInum.lck.stim<-ddply(stim.lck.OT.dist, c("Animal","Spot","Condition","Channel"), summarise, nROIs=length(OnsetTime))

# add in number of trials
ROInum.lck.stim$Ani_Spot_Cond<-paste(ROInum.lck.stim$Animal, ROInum.lck.stim$Spot, ROInum.lck.stim$Condition, sep="_")

ROInum.lck.stim<-merge(ROInum.lck.stim, Spot.lck.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)

ROInum.lck.stim$ROIsPerTrial<-ROInum.lck.stim$nROIs/ROInum.lck.stim$nTrials


# mean
df.lck.ROInum.mean<-summarySE(ROInum.lck.stim, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition"))


ggplot(df.lck.ROInum.mean, aes(x=Channel,y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view") +
  max.theme

Condition_Channel2= interaction(ROInum.lck.stim$Condition,ROInum.lck.stim$Channel)
nROI.lck.stim.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.lck.stim,REML=FALSE)
nROI.lck.stim.model1 = lmer(ROIsPerTrial~ Channel + (1|Animal), ROInum.lck.stim,REML=FALSE)
nROI.lck.stim.model2 = lmer(ROIsPerTrial ~ Condition + (1|Animal), ROInum.lck.stim,REML=FALSE)
nROI.lck.stim.model3 = lmer(ROIsPerTrial ~ Condition_Channel2 + (1|Animal), ROInum.lck.stim,REML=FALSE)
nROI.lck.stim.anova <- anova(nROI.lck.stim.null, nROI.lck.stim.model1,nROI.lck.stim.model2,nROI.lck.stim.model3)
print(nROI.lck.stim.anova)

nROI.lck.stim.Cond_Channel<- glht(nROI.lck.stim.model3, mcp(Condition_Channel2= "Tukey"))
summary(nROI.lck.stim.Cond_Channel)


##
# CYTO data

stim.cyto.OT.dist$Channel <- factor(stim.cyto.OT.dist$Channel, levels = c("RCaMP","GCaMP"))

ROInum.cyto.stim<-ddply(stim.cyto.OT.dist, c("Animal","Spot","Condition","Channel"), summarise, nROIs=length(OnsetTime))

# add in number of trials
ROInum.cyto.stim$Ani_Spot_Cond<-paste(ROInum.cyto.stim$Animal, ROInum.cyto.stim$Spot, ROInum.cyto.stim$Condition, sep="_")

ROInum.cyto.stim<-merge(ROInum.cyto.stim, Spot.cyto.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)

ROInum.cyto.stim$ROIsPerTrial<-ROInum.cyto.stim$nROIs/ROInum.cyto.stim$nTrials


# mean
df.cyto.ROInum.mean<-summarySE(ROInum.cyto.stim, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition"))


ggplot(df.cyto.ROInum.mean, aes(x=Channel,y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view") +
  max.theme

Condition_Channel2= interaction(ROInum.cyto.stim$Condition,ROInum.cyto.stim$Channel)
nROI.cyto.stim.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.cyto.stim,REML=FALSE)
nROI.cyto.stim.model1 = lmer(ROIsPerTrial~ Channel + (1|Animal), ROInum.cyto.stim,REML=FALSE)
nROI.cyto.stim.model2 = lmer(ROIsPerTrial ~ Condition + (1|Animal), ROInum.cyto.stim,REML=FALSE)
nROI.cyto.stim.model3 = lmer(ROIsPerTrial ~ Condition_Channel2 + (1|Animal), ROInum.cyto.stim,REML=FALSE)
nROI.cyto.stim.anova <- anova(nROI.cyto.stim.null, nROI.cyto.stim.model1,nROI.cyto.stim.model2,nROI.cyto.stim.model3)
print(nROI.cyto.stim.anova)

nROI.cyto.stim.Cond_Channel<- glht(nROI.cyto.stim.model3, mcp(Condition_Channel2= "Tukey"))
summary(nROI.cyto.stim.Cond_Channel)


#####
# mean onset times for all sensors

# first 10s after stimulation
Long_OTwind=10

# lck
stim.lck.OT<-subset(all.lck.OT,Condition!="shortstim")
stim.lck.OT.window<-subset(stim.lck.OT, OnsetTime<=Long_OTwind)

#cyto
stim.cyto.OT<-subset(all.cyto.OT,Condition!="shortstim")
stim.cyto.OT.window<-subset(stim.cyto.OT, OnsetTime<=Long_OTwind)


# combine cyto and cyto data sets for onset times
stim.lck.OT.window$sensor<- "lck"
stim.cyto.OT.window$sensor<-"cyto"
stim.lck.OT.window$Channel<-paste(stim.lck.OT.window$sensor, stim.lck.OT.window$Channel, sep="_")
stim.cyto.OT.window$Channel<-paste(stim.cyto.OT.window$sensor, stim.cyto.OT.window$Channel, sep="_")

stim.both.OT.window<-rbind(stim.lck.OT.window,stim.cyto.OT.window)

stim.both.OT.window$Channel <- factor(stim.both.OT.window$Channel, levels = c("cyto_RCaMP","cyto_GCaMP","lck_RCaMP","lck_GCaMP"))


#means
df.both.OT.mean<-summarySE(stim.both.OT.window, measurevar = "OnsetTime", groupvars = c("Channel", "Condition"))
df.both.OT.mean.stim<-summarySE(stim.both.OT.window[stim.both.OT.window$Condition=="Stim",], measurevar = "OnsetTime", groupvars = c("Channel", "Condition"))

ggplot(df.both.OT.mean, aes(x=Channel,y=OnsetTime, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme

ggplot(df.both.OT.mean.stim, aes(x=Channel,y=OnsetTime, fill= Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme


ggplot(stim.both.OT.window, aes(x=Channel,y=OnsetTime, fill= Condition)) +
  geom_boxplot(notch=TRUE)+
  ylab("Onset Time (s)") +
  ggtitle("notched")+
  max.theme

# means stats
Condition_Channel=interaction(stim.both.OT.window$Condition,stim.both.OT.window$Channel)

OT.both.stim.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial_Cond), stim.both.OT.window,REML=FALSE)
OT.both.stim.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Spot_trial_Cond), stim.both.OT.window,REML=FALSE)
OT.both.stim.model2 = lmer(OnsetTime ~ Condition + (1|Animal) + (1|Spot) + (1|Spot_trial_Cond) , stim.both.OT.window,REML=FALSE)
OT.both.stim.model3 = lmer(OnsetTime ~ Condition_Channel + (1|Animal) + (1|Spot) + (1|Spot_trial_Cond), stim.both.OT.window,REML=FALSE)
OT.both.stim.anova <- anova(OT.both.stim.null, OT.both.stim.model1,OT.both.stim.model2,OT.both.stim.model3)
print(OT.both.stim.anova)

OT.both.stim.Cond_Channel<- glht(OT.both.stim.model3, mcp(Condition_Channel= "Tukey"))
summary(OT.both.stim.Cond_Channel)


#only stim
OT.both.onlystim.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial_Cond), stim.both.OT.window[stim.both.OT.window$Condition=="Stim",],REML=FALSE)
OT.both.onlystim.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Spot_trial_Cond), stim.both.OT.window[stim.both.OT.window$Condition=="Stim",],REML=FALSE)
OT.both.onlystim.anova <- anova(OT.both.onlystim.null, OT.both.onlystim.model1)
print(OT.both.onlystim.anova)

OT.both.onlystim.Channel<- glht(OT.both.onlystim.model1, mcp(Channel= "Tukey"))
summary(OT.both.onlystim.Channel)


# compare onset time boxplots and their variances

# ks test
OT.lck.NSvsS.kstest<- ks.test(stim.both.OT.window$OnsetTime[stim.both.OT.window$Channel=="lck_GCaMP" & stim.both.OT.window$Condition=="Nostim"],
                              stim.both.OT.window$OnsetTime[stim.both.OT.window$Channel=="lck_GCaMP" & stim.both.OT.window$Condition=="Stim"])
print(OT.lck.NSvsS.kstest)

OT.cyto.NSvsS.kstest<- ks.test(stim.both.OT.window$OnsetTime[stim.both.OT.window$Channel=="cyto_GCaMP" & stim.both.OT.window$Condition=="Nostim"],
                               stim.both.OT.window$OnsetTime[stim.both.OT.window$Channel=="cyto_GCaMP" & stim.both.OT.window$Condition=="Stim"])
print(OT.cyto.NSvsS.kstest)


# lck stim vs cyto stim
OT.both.SvsS.kstest<- ks.test(stim.both.OT.window$OnsetTime[stim.both.OT.window$Channel=="lck_GCaMP" & stim.both.OT.window$Condition=="Stim"],
                              stim.both.OT.window$OnsetTime[stim.both.OT.window$Channel=="cyto_GCaMP" & stim.both.OT.window$Condition=="Stim"])
print(OT.both.SvsS.kstest)


# Anderson Darling
library("kSamples")
#library("SuppDist")

#lck no stim vs stim
OT.lck.NSvsS.adtest<- ad.test(stim.both.OT.window$OnsetTime[stim.both.OT.window$Channel=="lck_GCaMP" & stim.both.OT.window$Condition=="Nostim"],
                              stim.both.OT.window$OnsetTime[stim.both.OT.window$Channel=="lck_GCaMP" & stim.both.OT.window$Condition=="Stim"])
print(OT.lck.NSvsS.adtest)


#cyto no stim vs stim
OT.cyto.NSvsS.adtest<- ad.test(stim.both.OT.window$OnsetTime[stim.both.OT.window$Channel=="cyto_GCaMP" & stim.both.OT.window$Condition=="Nostim"],
                               stim.both.OT.window$OnsetTime[stim.both.OT.window$Channel=="cyto_GCaMP" & stim.both.OT.window$Condition=="Stim"])
print(OT.cyto.NSvsS.adtest)

# lck stim vs cyto stim
OT.both.SvsS.adtest<- ad.test(stim.both.OT.window$OnsetTime[stim.both.OT.window$Channel=="lck_GCaMP" & stim.both.OT.window$Condition=="Stim"],
                              stim.both.OT.window$OnsetTime[stim.both.OT.window$Channel=="cyto_GCaMP" & stim.both.OT.window$Condition=="Stim"])
print(OT.both.SvsS.adtest)

# lck no stim vs cyto no stim 
OT.both.NSvsNS.adtest<- ad.test(stim.both.OT.window$OnsetTime[stim.both.OT.window$Channel=="lck_GCaMP" & stim.both.OT.window$Condition=="Nostim"],
                                stim.both.OT.window$OnsetTime[stim.both.OT.window$Channel=="cyto_GCaMP" & stim.both.OT.window$Condition=="Nostim"])
print(OT.both.NSvsNS.adtest)


########
# mean peak time calculations

# in first 12 s after stimulation

Long_PTwind=12


# lck
stim.lck.PT<-subset(all.lck.peaks,Condition!="shortstim")
stim.lck.peaks.window<-subset(stim.lck.PT, peakTime<=Long_PTwind & peakTime>=0 & Duration<45)


stim.cyto.PT<-subset(all.cyto.peaks,Condition!="shortstim")
stim.cyto.peaks.window<-subset(stim.cyto.PT, peakTime<=Long_PTwind & peakTime>=0 & Duration<55)


# combine lck and cyto data sets for peak times

stim.lck.peaks.window$sensor<- "lck"
stim.cyto.peaks.window$sensor<-"cyto"
stim.lck.peaks.window$Channel<-paste(stim.lck.peaks.window$sensor, stim.lck.peaks.window$Channel, sep="_")
stim.cyto.peaks.window$Channel<-paste(stim.cyto.peaks.window$sensor, stim.cyto.peaks.window$Channel, sep="_")

stim.both.peaks.window<-rbind(stim.lck.peaks.window,stim.cyto.peaks.window)

stim.both.peaks.window$Channel <- factor(stim.both.peaks.window$Channel, levels = c("cyto_RCaMP","cyto_GCaMP","lck_RCaMP","lck_GCaMP"))


# means

df.both.PT.mean<-summarySE(stim.both.peaks.window, measurevar = "peakTime", groupvars = c("Channel", "Condition"))
df.both.PT.mean.stim<-summarySE(stim.both.peaks.window[stim.both.peaks.window$Condition=="Stim",], measurevar = "peakTime", groupvars = c("Channel", "Condition"))


ggplot(df.both.PT.mean, aes(x=Channel,y=peakTime, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean peak Time (s)") +
  max.theme

ggplot(df.both.PT.mean.stim, aes(x=Channel,y=peakTime, fill= Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Peak Time (s)") +
  max.theme

ggplot(stim.both.peaks.window, aes(x=Channel,y=peakTime, fill= Condition)) +
  geom_boxplot(notch=TRUE)+
  ylab("peak Time (s)") +
  ggtitle("notched")+
  max.theme


# means stats
Condition_Channel2=interaction(stim.both.peaks.window$Condition,stim.both.peaks.window$Channel)

PT.both.stim.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials_Cond), stim.both.peaks.window,REML=FALSE)
PT.both.stim.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials_Cond), stim.both.peaks.window,REML=FALSE)
PT.both.stim.model2 = lmer(peakTime ~ Condition + (1|Animal) + (1|Spot) + (1|trials_Cond) , stim.both.peaks.window,REML=FALSE)
PT.both.stim.model3 = lmer(peakTime ~ Condition_Channel2 + (1|Animal) + (1|Spot) + (1|trials_Cond), stim.both.peaks.window,REML=FALSE)
PT.both.stim.anova <- anova(PT.both.stim.null, PT.both.stim.model1,PT.both.stim.model2,PT.both.stim.model3)
print(PT.both.stim.anova)

PT.both.stim.Cond_Channel<- glht(PT.both.stim.model3, mcp(Condition_Channel2= "Tukey"))
summary(PT.both.stim.Cond_Channel)


#only stim
PT.both.onlystim.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials_Cond), stim.both.peaks.window[stim.both.peaks.window$Condition=="Stim",],REML=FALSE)
PT.both.onlystim.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials_Cond), stim.both.peaks.window[stim.both.peaks.window$Condition=="Stim",],REML=FALSE)
PT.both.onlystim.anova <- anova(PT.both.onlystim.null, PT.both.onlystim.model1)
print(PT.both.onlystim.anova)

PT.both.onlystim.Channel<- glht(PT.both.onlystim.model1, mcp(Channel= "Tukey"))
summary(PT.both.onlystim.Channel)

# compare onset time distributions and their variances

# Anderson Darling
#library("kSamples")
#library("SuppDist")

#lck no stim vs stim
PT.lck.NSvsS.adtest<- ad.test(stim.both.peaks.window$peakTime[stim.both.peaks.window$Channel=="lck_GCaMP" & stim.both.peaks.window$Condition=="Nostim"],
                              stim.both.peaks.window$peakTime[stim.both.peaks.window$Channel=="lck_GCaMP" & stim.both.peaks.window$Condition=="Stim"])
print(PT.lck.NSvsS.adtest)

#cyto no stim vs stim
PT.cyto.NSvsS.adtest<- ad.test(stim.both.peaks.window$peakTime[stim.both.peaks.window$Channel=="cyto_GCaMP" & stim.both.peaks.window$Condition=="Nostim"],
                               stim.both.peaks.window$peakTime[stim.both.peaks.window$Channel=="cyto_GCaMP" & stim.both.peaks.window$Condition=="Stim"])
print(PT.cyto.NSvsS.adtest)

# lck stim vs cyto stim
PT.both.SvsS.adtest<- ad.test(stim.both.peaks.window$peakTime[stim.both.peaks.window$Channel=="lck_GCaMP" & stim.both.peaks.window$Condition=="Stim"],
                              stim.both.peaks.window$peakTime[stim.both.peaks.window$Channel=="cyto_GCaMP" & stim.both.peaks.window$Condition=="Stim"])
print(PT.both.SvsS.adtest)

# lck no stim vs cyto no stim 
PT.both.NSvsNS.adtest<- ad.test(stim.both.peaks.window$peakTime[stim.both.peaks.window$Channel=="lck_GCaMP" & stim.both.peaks.window$Condition=="Nostim"],
                                stim.both.peaks.window$peakTime[stim.both.peaks.window$Channel=="cyto_GCaMP" & stim.both.peaks.window$Condition=="Nostim"])
print(PT.both.NSvsNS.adtest)


########
# ROI number 1s before and 1s after stimulation


# Lck data 

lck.eitherside.stim2<-subset(stim.lck.OT.dist.2s, OnsetTimeAdjust<1 & OnsetTimeAdjust>-1)
lck.eitherside.stim2$either1s="other"
lck.eitherside.stim2$either1s[lck.eitherside.stim2$OnsetTimeAdjust<0]="before"
lck.eitherside.stim2$either1s[lck.eitherside.stim2$OnsetTimeAdjust>0]="after"

# sum number of ROIs per field of view
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

#####
#cyto data 

# considering 1 s before and 1 s after stimulation 

cyto.eitherside.stim2<-subset(stim.cyto.OT.dist.2s, OnsetTimeAdjust<1 & OnsetTimeAdjust>-1)
cyto.eitherside.stim2$either1s="other"
cyto.eitherside.stim2$either1s[cyto.eitherside.stim2$OnsetTimeAdjust<0]="before"
cyto.eitherside.stim2$either1s[cyto.eitherside.stim2$OnsetTimeAdjust>0]="after"


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



#########
# FAST vs DELAYED groups

# ONLY Lck DATA!

LongN_PTwind=8  # neurons with a peak within 8s (fast neurons)
LongAC_PTwind=15 # astrocytes with a peak within 15s (likely delayed AC)

LongN_OTwind=1  # neurons with an onset within 1s (fast neurons)
LongAC_OTwind=12  # astrocytes with an onset within 12s (likely delayed AC)


# onset
stim.lck.OT.R<-subset(stim.lck.OT, Channel=="RCaMP" & OnsetTime<=LongN_OTwind)
stim.lck.OT.G<-subset(stim.lck.OT, Channel=="GCaMP" & OnsetTime<=LongAC_OTwind)
stim.lck.OT.F_D<-rbind(stim.lck.OT.R, stim.lck.OT.G)

stim.cyto.OT.R<-subset(stim.cyto.OT, Channel=="RCaMP" & OnsetTime<=LongN_OTwind)
stim.cyto.OT.G<-subset(stim.cyto.OT, Channel=="GCaMP" & OnsetTime<=LongAC_OTwind)
stim.cyto.OT.F_D<-rbind(stim.cyto.OT.R, stim.cyto.OT.G)

# peak times
stim.lck.PT.R<-subset(stim.lck.PT, Channel=="RCaMP" & peakTime<=LongN_PTwind & peakTime>=0 & Duration<45)
stim.lck.PT.G<-subset(stim.lck.PT, Channel=="GCaMP" & peakTime<=LongAC_PTwind & peakTime>=0 & Duration<45)
stim.lck.peaks.F_D<-rbind(stim.lck.PT.R, stim.lck.PT.G)

stim.cyto.PT.R<-subset(stim.cyto.PT, Channel=="RCaMP" & peakTime<=LongN_PTwind & peakTime>=0 & Duration<55)
stim.cyto.PT.G<-subset(stim.cyto.PT, Channel=="GCaMP" & peakTime<=LongAC_PTwind & peakTime>=0 & Duration<55)
stim.cyto.peaks.F_D<-rbind(stim.cyto.PT.R, stim.cyto.PT.G)

# peak times only for ROIs with onset times
stim.lck.F_D<-merge(stim.lck.peaks.F_D, stim.lck.OT.F_D[, c("ROIs_Cond", "OnsetTime","TraceAUC1","TraceAUC10")], by="ROIs_Cond")
stim.cyto.F_D<-merge(stim.cyto.peaks.F_D, stim.cyto.OT.F_D[, c("ROIs_Cond", "OnsetTime","TraceAUC1","TraceAUC10")], by="ROIs_Cond")


# identify "FAST" astrocytes- only in Lck data
stim.lck.F_D$Group<-0
stim.lck.F_D$Group[stim.lck.F_D$OnsetTime<1]<-"fast"
stim.lck.F_D$Group[stim.lck.F_D$OnsetTime>=1]<-"delayed"

stim.lck.F_D$Group <- factor(stim.lck.F_D$Group, levels = c("fast","delayed"))
stim.lck.F_D$Channel <- factor(stim.lck.F_D$Channel, levels = c("RCaMP","GCaMP"))


stim.lck.F_D$Channel_Group<-interaction(stim.lck.F_D$Channel, stim.lck.F_D$Group)
stim.lck.F_D$Channel_Group<-as.factor(stim.lck.F_D$Channel_Group)

# take out the effect of Condition
# we are only interested in stim case
stim.lck.F_D.STIM<-subset(stim.lck.F_D, Condition=="Stim")

#######
# mean onset times
df.OT1<- summarySE(stim.lck.F_D, measurevar = "OnsetTime", groupvars = c("Channel_Group","Condition"))
df.OT2<- summarySE(stim.lck.F_D.STIM, measurevar = "OnsetTime", groupvars = c("Channel_Group"))
df.OT3<- summarySE(stim.lck.F_D.STIM[stim.lck.F_D.STIM$Channel=="GCaMP",], measurevar = "OnsetTime", groupvars = c("ROIType","Channel_Group"))


ggplot(df.OT1, aes(x=Channel_Group,y=OnsetTime, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.OT2, aes(x=Channel_Group,y=OnsetTime, fill=Channel_Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.OT3, aes(x=Channel_Group,y=OnsetTime, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

# stats
Group_Channel_Cond_Type=interaction(stim.lck.F_D$Group,stim.lck.F_D$Channel,stim.lck.F_D$Condition,stim.lck.F_D$ROIType)
Group_Channel_Cond=interaction(stim.lck.F_D$Channel_Group,stim.lck.F_D$Condition)

# stats for onset times- nostim and stim
OT.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
OT.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
OT.model2 = lmer(OnsetTime ~ Condition + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
OT.model3 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
OT.model4 = lmer(OnsetTime ~ Channel_Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
OT.model5 = lmer(OnsetTime ~ Group_Channel_Cond + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
OT.model6 = lmer(OnsetTime ~ Group_Channel_Cond_Type + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
OT.anova <- anova(OT.null, OT.model1,OT.model2,OT.model3,OT.model4,OT.model5,OT.model6)
print(OT.anova)

OT.Group_channel<- glht(OT.model5, mcp(Group_Channel_Cond= "Tukey"))
summary(OT.Group_channel)


# stim only
Group_Channel_Type=interaction(stim.lck.F_D.STIM$Group,stim.lck.F_D.STIM$Channel,stim.lck.F_D.STIM$ROIType)
Group_Channel=interaction(stim.lck.F_D.STIM$Channel, stim.lck.F_D.STIM$Group)

# stats for onset times- neurons vs astrocytes
OT.stim.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM,REML=FALSE)
OT.stim.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM,REML=FALSE)
OT.stim.model3 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM,REML=FALSE)
OT.stim.model4 = lmer(OnsetTime ~ Group_Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM,REML=FALSE)
OT.stim.model6 = lmer(OnsetTime ~ Group_Channel_Type + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM,REML=FALSE)
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

######
# mean peak times
df.PT1<- summarySE(stim.lck.F_D, measurevar = "peakTime", groupvars = c("Channel_Group","Condition"))
df.PT2<- summarySE(stim.lck.F_D.STIM, measurevar = "peakTime", groupvars = c("Channel_Group"))
df.PT3<- summarySE(stim.lck.F_D.STIM[stim.lck.F_D.STIM$Channel=="GCaMP",], measurevar = "peakTime", groupvars = c("ROIType","Channel_Group"))


ggplot(df.PT1, aes(x=Channel_Group,y=peakTime, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.PT2, aes(x=Channel_Group,y=peakTime, fill=Channel_Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.PT3, aes(x=Channel_Group,y=peakTime, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

# stats

# stats for onset times- nostim and stim
PT.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
PT.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
PT.model2 = lmer(peakTime ~ Condition + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
PT.model3 = lmer(peakTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
PT.model4 = lmer(peakTime ~ Channel_Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
PT.model5 = lmer(peakTime ~ Group_Channel_Cond + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
PT.model6 = lmer(peakTime ~ Group_Channel_Cond_Type + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
PT.anova <- anova(PT.null, PT.model1,PT.model2,PT.model3,PT.model4,PT.model5,PT.model6)
print(PT.anova)

PT.Group_channel<- glht(PT.model5, mcp(Group_Channel_Cond= "Tukey"))
summary(PT.Group_channel)


# stim only

# stats for onset times- neurons vs astrocytes
PT.stim.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM,REML=FALSE)
PT.stim.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM,REML=FALSE)
PT.stim.model3 = lmer(peakTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM,REML=FALSE)
PT.stim.model4 = lmer(peakTime ~ Group_Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM,REML=FALSE)
PT.stim.model6 = lmer(peakTime ~ Group_Channel_Type + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM,REML=FALSE)
PT.stim.anova <- anova(PT.stim.null, PT.stim.model1,PT.stim.model3,PT.stim.model4,PT.stim.model6)
print(PT.stim.anova)

PT.stim.Group_channel<- glht(PT.stim.model4, mcp(Group_Channel= "Tukey"))
summary(PT.stim.Group_channel)

PT.stim.Group_channel_type<- glht(PT.stim.model6, mcp(Group_Channel_Type= "Tukey"))
summary(PT.stim.Group_channel_type)

summary(PT.stim.model4)

# check residuals for linearity
plot(fitted(PT.stim.model4), residuals(PT.stim.model4),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(PT.stim.model4), residuals(PT.stim.model4)), col=46, lwd=2.5)


########
#amplitude

# comparing stim and no stim
df.amp1.lck<- summarySE(stim.lck.F_D[stim.lck.F_D$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("Condition"))
df.amp1.cyto<- summarySE(stim.cyto.F_D[stim.cyto.F_D$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("Condition"))

ggplot(df.amp1.lck, aes(x=Condition,y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp1.cyto, aes(x=Condition,y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

# for supplementary figure histogram
ggplot(stim.lck.F_D[stim.lck.F_D$Channel=="GCaMP",], aes(x=amplitude, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.5, position="dodge")+
  ylab("density") +
  xlim(0,12)+
  ggtitle("lck data- no stim vs stim")+
  max.theme

ggplot(stim.cyto.F_D[stim.cyto.F_D$Channel=="GCaMP",], aes(x=amplitude, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.5, position="dodge")+
  ylab("density") +
  xlim(0,12)+
  ggtitle("lck data- no stim vs stim")+
  max.theme


# lck no stim vs stim
amp.lck.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D[stim.lck.F_D$Channel=="GCaMP",],REML=FALSE)
amp.lck.model1 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.F_D[stim.lck.F_D$Channel=="GCaMP",],REML=FALSE)

amp.lck.anova <- anova(amp.lck.null, amp.lck.model1)
print(amp.lck.anova)

amp.lck.cond<- glht(amp.lck.model1, mcp(Condition= "Tukey"))
summary(amp.lck.cond)

# cyto no stim vs stim
amp.cyto.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), stim.cyto.F_D[stim.cyto.F_D$Channel=="GCaMP",],REML=FALSE)
amp.cyto.model1 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.cyto.F_D[stim.cyto.F_D$Channel=="GCaMP",],REML=FALSE)

amp.cyto.anova <- anova(amp.cyto.null, amp.cyto.model1)
print(amp.cyto.anova)

amp.cyto.cond<- glht(amp.cyto.model1, mcp(Condition= "Tukey"))
summary(amp.cyto.cond)


##
# comparing astrocyte groups
df.amp2<- summarySE(stim.lck.F_D, measurevar = "amplitude", groupvars = c("Channel_Group","Condition"))
df.amp3<- summarySE(stim.lck.F_D.STIM[stim.lck.F_D.STIM$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("Channel_Group"))
df.amp4<- summarySE(stim.lck.F_D.STIM[stim.lck.F_D.STIM$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("ROIType","Channel_Group"))

ggplot(df.amp2, aes(x=Channel_Group,y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp3, aes(x=Channel_Group,y=amplitude, fill= Channel_Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp4, aes(x=Channel_Group,y=amplitude, fill= ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

## lck data and group (fast, delayed)
comp_Channel_Group=interaction(stim.lck.F_D$Group,stim.lck.F_D$Channel)
comp_Channel_Group_Cond=interaction(stim.lck.F_D$Group,stim.lck.F_D$Channel,stim.lck.F_D$Condition)
comp_Channel_Group_Cond_Type=interaction(stim.lck.F_D$Group,stim.lck.F_D$Channel,
                                         stim.lck.F_D$Condition,stim.lck.F_D$ROIType)


# stats for onset times- neurons vs astrocytes
amp.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
amp.model1 = lmer(amplitude ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
amp.model2 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.F_D,REML=FALSE)
amp.model3 = lmer(amplitude ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
amp.model4 = lmer(amplitude ~ comp_Channel_Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
amp.model5 = lmer(amplitude ~ comp_Channel_Group_Cond + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
amp.model6 = lmer(amplitude ~ comp_Channel_Group_Cond_Type + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
amp.anova <- anova(amp.null, amp.model1,amp.model2,amp.model3,amp.model4,amp.model5,amp.model6)
print(amp.anova)

amp.Group_channel<- glht(amp.model5, mcp(comp_Channel_Group_Cond= "Tukey"))
summary(amp.Group_channel)

amp.Group_channel_ty<- glht(amp.model6, mcp(comp_Channel_Group_Cond_Type= "Tukey"))
summary(amp.Group_channel_ty)



# consider only the Lck data and STIM case

#lck-GCaMP
Group_Type_GC=interaction(stim.lck.F_D.STIM$Group[stim.lck.F_D.STIM$Channel=="GCaMP"],stim.lck.F_D.STIM$ROIType[stim.lck.F_D.STIM$Channel=="GCaMP"])

amp.null.GC = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM[stim.lck.F_D.STIM$Channel=="GCaMP",],REML=FALSE)
amp.model1.GC  = lmer(amplitude ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM[stim.lck.F_D.STIM$Channel=="GCaMP",],REML=FALSE)
amp.model2.GC  = lmer(amplitude ~ Group_Type_GC + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM[stim.lck.F_D.STIM$Channel=="GCaMP",],REML=FALSE)
amp.anova.GC  <- anova(amp.null.GC, amp.model1.GC,amp.model2.GC)
print(amp.anova.GC)

amp.Group.GC<- glht(amp.model1.GC, mcp(Group= "Tukey"))
summary(amp.Group.GC)

amp.Group.type.GC<- glht(amp.model2.GC, mcp(Group_Type_GC= "Tukey"))
summary(amp.Group.type.GC)


plot(fitted(amp.model2.GC), residuals(amp.model2.GC),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(amp.model2.GC), residuals(amp.model2.GC)), col=46, lwd=2.5)

########
#duration
# comparing stim and no stim
df.dur1.lck<- summarySE(stim.lck.F_D[stim.lck.F_D$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Condition"))
df.dur1.cyto<- summarySE(stim.cyto.F_D[stim.cyto.F_D$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Condition"))

ggplot(df.dur1.lck, aes(x=Condition,y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.dur1.cyto, aes(x=Condition,y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

# for supplementary figure histogram
ggplot(stim.lck.F_D[stim.lck.F_D$Channel=="GCaMP",], aes(x=Duration, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.5, position="dodge")+
  ylab("density") +
  ggtitle("lck data- no stim vs stim")+
  max.theme

ggplot(stim.cyto.F_D[stim.cyto.F_D$Channel=="GCaMP",], aes(x=Duration, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.5, position="dodge")+
  ylab("density") +
  ggtitle("cyto data- no stim vs stim")+
  max.theme


# lck no stim vs stim
dur.lck.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D[stim.lck.F_D$Channel=="GCaMP",],REML=FALSE)
dur.lck.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.F_D[stim.lck.F_D$Channel=="GCaMP",],REML=FALSE)

dur.lck.anova <- anova(dur.lck.null, dur.lck.model1)
print(dur.lck.anova)

dur.lck.cond<- glht(dur.lck.model1, mcp(Condition= "Tukey"))
summary(dur.lck.cond)

# cyto no stim vs stim
dur.cyto.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.cyto.F_D[stim.cyto.F_D$Channel=="GCaMP",],REML=FALSE)
dur.cyto.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.cyto.F_D[stim.cyto.F_D$Channel=="GCaMP",],REML=FALSE)

dur.cyto.anova <- anova(dur.cyto.null, dur.cyto.model1)
print(dur.cyto.anova)

dur.cyto.cond<- glht(dur.cyto.model1, mcp(Condition= "Tukey"))
summary(dur.cyto.cond)


##
# comparing astrocyte groups
df.dur2<- summarySE(stim.lck.F_D, measurevar = "Duration", groupvars = c("Channel_Group","Condition"))
df.dur3<- summarySE(stim.lck.F_D.STIM[stim.lck.F_D.STIM$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Channel_Group"))
df.dur4<- summarySE(stim.lck.F_D.STIM[stim.lck.F_D.STIM$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("ROIType","Channel_Group"))

ggplot(df.dur2, aes(x=Channel_Group,y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.dur3, aes(x=Channel_Group,y=Duration, fill= Channel_Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.dur4, aes(x=Channel_Group,y=Duration, fill= ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

## lck data and group (fast, delayed)

# stats for onset times- neurons vs astrocytes
dur.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
dur.model1 = lmer(Duration ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
dur.model2 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.F_D,REML=FALSE)
dur.model3 = lmer(Duration ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
dur.model4 = lmer(Duration ~ comp_Channel_Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
dur.model5 = lmer(Duration ~ comp_Channel_Group_Cond + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
dur.model6 = lmer(Duration ~ comp_Channel_Group_Cond_Type + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D,REML=FALSE)
dur.anova <- anova(dur.null, dur.model1,dur.model2,dur.model3,dur.model4,dur.model5,dur.model6)
print(dur.anova)

dur.Group_channel<- glht(dur.model5, mcp(comp_Channel_Group_Cond= "Tukey"))
summary(dur.Group_channel)

dur.Group_channel_ty<- glht(dur.model6, mcp(comp_Channel_Group_Cond_Type= "Tukey"))
summary(dur.Group_channel_ty)



# consider only the Lck data and STIM case

#lck-GCaMP
dur.null.GC = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM[stim.lck.F_D.STIM$Channel=="GCaMP",],REML=FALSE)
dur.model1.GC  = lmer(Duration ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM[stim.lck.F_D.STIM$Channel=="GCaMP",],REML=FALSE)
dur.model2.GC  = lmer(Duration ~ Group_Type_GC + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM[stim.lck.F_D.STIM$Channel=="GCaMP",],REML=FALSE)
dur.anova.GC  <- anova(dur.null.GC, dur.model1.GC,dur.model2.GC)
print(dur.anova.GC)

dur.Group.GC<- glht(dur.model1.GC, mcp(Group= "Tukey"))
summary(dur.Group.GC)

dur.Group.type.GC<- glht(dur.model2.GC, mcp(Group_Type_GC= "Tukey"))
summary(dur.Group.type.GC)


######

# Process ROI area
df.Rarea1<- summarySE(stim.lck.F_D[stim.lck.F_D$ROIType=="Process",], measurevar = "area", groupvars = c("Channel_Group","Condition"))
df.Rarea2<- summarySE(stim.lck.F_D.STIM[stim.lck.F_D.STIM$Channel=="GCaMP"& stim.lck.F_D.STIM$ROIType=="Process",], measurevar = "area", groupvars = c("Channel_Group"))


ggplot(df.Rarea1, aes(x=Channel_Group,y=area, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("area") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.Rarea2, aes(x=Channel_Group,y=area, fill= Channel_Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("area") +
  scale_fill_manual(values=cbbPalette)+
  max.theme


#only consider STIM case
area.null.stim = lmer(area ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM[stim.lck.F_D.STIM$Channel!="lck_RCaMP" & stim.lck.F_D.STIM$ROIType=="Process",],REML=FALSE)
area.model1.stim  = lmer(area ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.F_D.STIM[stim.lck.F_D.STIM$Channel!="lck_RCaMP" & stim.lck.F_D.STIM$ROIType=="Process",],REML=FALSE)

area.anova.stim  <- anova(area.null.stim, area.model1.stim)
print(area.anova.stim)

area.stim.group<- glht(area.model1.stim, mcp(Group= "Tukey"))
summary(area.stim.group)



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

respGC<-subset(stim.lck.F_D.STIM, Channel=="GCaMP") 
respRC<-subset(stim.lck.F_D.STIM, Channel=="RCaMP") 

# list of responding astrocytes and their corresponding group
respGC2<-unique(respGC[c("ROIs_trial","Group")])
respRC2<-unique(respRC[c("ROIs_trial","Group")])

# only correlations from responding astrocytes and neurons
GCaMP_RCaMP<-subset(GCaMP_RCaMP, ROIs_trial %in% respGC2$ROIs_trial)
GCaMP_RCaMP<-subset(GCaMP_RCaMP, RCaMP_ROIs %in% unique(respRC$ROIs_trial))

# put group information into correlation table ("fast or delayed")
GCaMP_RCaMP.group<-merge(GCaMP_RCaMP, stim.lck.F_D.STIM[, c("ROIs_trial", "Group")], by="ROIs_trial")


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

# linear correlation 
shortcorr.lck.null = lmer(Short_Corr ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP_RCaMP.group,REML=FALSE)
shortcorr.lck.model1 = lmer(Short_Corr ~ Group + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP_RCaMP.group,REML=FALSE)
shortcorr.lck.anova <- anova(shortcorr.lck.null, shortcorr.lck.model1)
print(shortcorr.lck.anova)

shortcorr.lck.group<- glht(shortcorr.lck.model1, mcp(Group= "Tukey"))
summary(shortcorr.lck.group)

#####
# cross correlation
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

ggplot(data=df6A3, aes(x=CompType, y=xCorr, fill=Group)) +
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

#cross correlation 
xcorr.lck.null = lmer(xCorr ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP_RCaMP.group,REML=FALSE)
xcorr.lck.model1 = lmer(xCorr ~ Group + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP_RCaMP.group,REML=FALSE)
xcorr.lck.anova <- anova(xcorr.lck.null, xcorr.lck.model1)
print(xcorr.lck.anova)

xcorr.lck.group<- glht(xcorr.lck.model1, mcp(Group= "Tukey"))
summary(xcorr.lck.group)

#####
# lag

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

ggplot(data=df7B1, aes(x=Group, y=Lag, fill=Group)) +
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

#cross correlation 
lag.lck.null = lmer(Lag ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP_RCaMP.group,REML=FALSE)
lag.lck.model1 = lmer(Lag ~ Group + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP_RCaMP.group,REML=FALSE)
lag.lck.anova <- anova(lag.lck.null, lag.lck.model1)
print(lag.lck.anova)

lag.lck.group<- glht(lag.lck.model1, mcp(Group= "Tukey"))
summary(lag.lck.group)


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
AvsN.Space.Nostim<-AvsN.Space.Nostim[AvsN.Space.Nostim$Group!="other",]
AvsN.Space.Nostim<-subset(AvsN.Space.Nostim, TimeDiff>-12 & TimeDiff<12)

AvsN.Space.Stim<-AvsN.Space.Stim[AvsN.Space.Stim$Group!="other",]
AvsN.Space.Stim<-subset(AvsN.Space.Stim, TimeDiff>-12 & TimeDiff<12)

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


NS.AvsN_space.median<-quantile(AvsN.Space.Nostim$TimeDiff,
                           prob = seq(0, 1, length = 5),type=5, na.rm=TRUE)
print(NS.AvsN_space.median)

S.AvsN_space.median<-quantile(AvsN.Space.Stim$TimeDiff,
                               prob = seq(0, 1, length = 5),type=5, na.rm=TRUE)
print(S.AvsN_space.median)

# compare distributions
AvsN.Space.kstest<- ks.test(AvsN.Space.Stim$TimeDiff,
                            AvsN.Space.Nostim$TimeDiff)
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


############



#####
# cyto DSP4 peak time data

# find spot 1 and spot 2
all.cyto.DSP4.peaks$SpotTime<-0
all.cyto.DSP4.peaks$SpotTime[grepl("spot1",all.cyto.DSP4.peaks$Spot)]="Spot1"
all.cyto.DSP4.peaks$SpotTime[grepl("spot2",all.cyto.DSP4.peaks$Spot)]="Spot2"
all.cyto.DSP4.peaks$SpotTime<-as.factor(all.cyto.DSP4.peaks$SpotTime)

# control in the first 2 days
all.cyto.DSP4.peaks$Drug="DSP4"
all.cyto.DSP4.peaks$Drug[grepl("04_06",all.cyto.DSP4.peaks$Spot)]="control"
all.cyto.DSP4.peaks$Drug[grepl("04_08",all.cyto.DSP4.peaks$Spot)]="control"


######
# astrocyte cyto DSP4 

# COMPARE CONTORL AND TREATMENT
# mean duration times
all.cyto.DSP4.peaks$Condition<- factor(all.cyto.DSP4.peaks$Condition, levels = c("Nostim","shortstim","Stim"))

df.dur.cyto4<- summarySE(all.cyto.DSP4.peaks, measurevar = "Duration", groupvars = c("Channel", "Condition","Drug"))



ggplot(df.dur.cyto4[df.dur.cyto4$Channel!="RCaMP",], aes(x=Condition,y=Duration, fill= Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration (s)") +
  ggtitle("astrocytes") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.dur.cyto4[df.dur.cyto4$Channel=="RCaMP",], aes(x=Condition,y=Duration, fill= Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration (s)") +
  ggtitle("neurons") +
  scale_fill_manual(values=cbbPalette)+
  max.theme


#astrocytes with before or after for the same spot
dur.cytovsDSP2.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
dur.cytovsDSP2.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
dur.cytovsDSP2.model2 = lmer(Duration ~ Drug+ (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
dur.cytovsDSP2.model3 = lmer(Duration ~ Condition+Drug + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
dur.cytovsDSP2.model4 = lmer(Duration ~ Drug_Condition + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
dur.cytovsDSP2.anova <- anova(dur.cytovsDSP2.null, dur.cytovsDSP2.model1, dur.cytovsDSP2.model2, 
                              dur.cytovsDSP2.model3, dur.cytovsDSP2.model4)
print(dur.cytovsDSP2.anova)

dur.cytovsDSP.treatment_Condition2<- glht(dur.cytovsDSP2.model4, mcp(Drug_Condition= "Tukey"))
summary(dur.cytovsDSP.treatment_Condition2)


dur.cytovsDSP2.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="RCaMP",],REML=FALSE)
dur.cytovsDSP2.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="RCaMP",],REML=FALSE)
dur.cytovsDSP2.model2 = lmer(Duration ~ Drug+ (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="RCaMP",],REML=FALSE)
dur.cytovsDSP2.model3 = lmer(Duration ~ Condition+Drug + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="RCaMP",],REML=FALSE)
dur.cytovsDSP2.model4 = lmer(Duration ~ Drug_Condition + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="RCaMP",],REML=FALSE)
dur.cytovsDSP2.anova <- anova(dur.cytovsDSP2.null, dur.cytovsDSP2.model1, dur.cytovsDSP2.model2, 
                              dur.cytovsDSP2.model3, dur.cytovsDSP2.model4)
print(dur.cytovsDSP2.anova)

dur.cytovsDSP.treatment_Condition3<- glht(dur.cytovsDSP2.model4, mcp(Drug_Condition= "Tukey"))
summary(dur.cytovsDSP.treatment_Condition3)




######
# astrocyte cyto DSP4 

# COMPARE CONTORL AND TREATMENT

# comparing to early time point controls
df.amp.cyto5<- summarySE(all.cyto.DSP4.peaks, measurevar = "amplitude", groupvars = c("Channel", "Condition","Drug"))
df.amp.cyto6<- summarySE(all.cyto.DSP4.peaks, measurevar = "amplitude", groupvars = c("Channel", "Condition","Drug","ROIType"))
df.amp.cyto7<- summarySE(all.cyto.DSP4.peaks, measurevar = "amplitude", groupvars = c("Channel", "Condition","Drug","SpotTime"))


ggplot(df.amp.cyto5[df.amp.cyto5$Channel!="RCaMP",], aes(x=Condition,y=amplitude, fill= Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("astrocytes") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp.cyto5[df.amp.cyto5$Channel=="RCaMP",], aes(x=Condition,y=amplitude, fill= Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("neurons") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp.cyto6[df.amp.cyto6$Channel!="RCaMP",], aes(x=interaction(Condition,ROIType),y=amplitude, fill= Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("astrocytes") +
  scale_fill_manual(values=cbbPalette)+
  max.theme



#astrocytes with before or after
amp.cytovsDSP2.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
amp.cytovsDSP2.model1 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
amp.cytovsDSP2.model2 = lmer(amplitude ~ Drug+ (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
amp.cytovsDSP2.model3 = lmer(amplitude ~ Condition+Drug + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
amp.cytovsDSP2.model4 = lmer(amplitude ~ Drug_Condition + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
amp.cytovsDSP2.anova <- anova(amp.cytovsDSP2.null, amp.cytovsDSP2.model1, amp.cytovsDSP2.model2, 
                              amp.cytovsDSP2.model3, amp.cytovsDSP2.model4)
print(amp.cytovsDSP2.anova)

amp.cytovsDSP.treatment_Condition2<- glht(amp.cytovsDSP2.model4, mcp(Drug_Condition= "Tukey"))
summary(amp.cytovsDSP.treatment_Condition2)


amp.cytovsDSP2.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="RCaMP",],REML=FALSE)
amp.cytovsDSP2.model1 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="RCaMP",],REML=FALSE)
amp.cytovsDSP2.model2 = lmer(amplitude ~ Drug+ (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="RCaMP",],REML=FALSE)
amp.cytovsDSP2.model3 = lmer(amplitude ~ Condition+Drug + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="RCaMP",],REML=FALSE)
amp.cytovsDSP2.model4 = lmer(amplitude ~ Drug_Condition + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="RCaMP",],REML=FALSE)
amp.cytovsDSP2.anova <- anova(amp.cytovsDSP2.null, amp.cytovsDSP2.model1, amp.cytovsDSP2.model2, 
                              amp.cytovsDSP2.model3, amp.cytovsDSP2.model4)
print(amp.cytovsDSP2.anova)

amp.cytovsDSP.treatment_Condition3<- glht(amp.cytovsDSP2.model4, mcp(Drug_Condition= "Tukey"))
summary(amp.cytovsDSP.treatment_Condition3)


######

# cyto and DSP4

# peaks per trial
all.cyto.DSP4.peaks<-subset(all.cyto.DSP4.peaks, peakTime<15)

all.cyto.DSP4.trials<-ddply(all.cyto.DSP4.peaks, c("Animal","Spot","trials","SpotTime","Condition","Channel","Drug"), summarise, nPeaks=length(amplitude))
all.cyto.DSP4.trials.type<-ddply(all.cyto.DSP4.peaks, c("Animal","Spot","trials","SpotTime","Condition","Channel","ROIType","Drug"), summarise, nPeaks=length(amplitude))


df.cyto.numPeaks4<-summarySE(all.cyto.DSP4.trials, measurevar="nPeaks", groupvars=c("Channel","Condition","Drug"))
df.cyto.numPeaks5<-summarySE(all.cyto.DSP4.trials.type, measurevar="nPeaks", groupvars=c("ROIType","Condition","Channel","Drug"))



ggplot(data=df.cyto.numPeaks4[df.cyto.numPeaks4$Channel=="GCaMP",], aes(x=Condition, y= nPeaks, fill=Drug)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=nPeaks-se, ymax=nPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("nPeaks/trial") +
  ggtitle("astrocytes")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme
  
  ggplot(data=df.cyto.numPeaks4[df.cyto.numPeaks4$Channel=="RCaMP",], aes(x=Condition, y= nPeaks, fill=Drug)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=nPeaks-se, ymax=nPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
    ylab("nPeaks/trial") +
    ggtitle("neurons")+
  scale_fill_manual(values=cbbPalette) + 
    max.theme
  
  
  ggplot(data=df.cyto.numPeaks5[df.cyto.numPeaks5$Channel=="GCaMP",], aes(x=interaction(Condition,ROIType), y= nPeaks, fill=Drug)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=nPeaks-se, ymax=nPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
    ylab("nPeaks/trial") +
    ggtitle("astrocytes")+
  scale_fill_manual(values=cbbPalette) + 
    max.theme
  
  ggplot(data=df.cyto.numPeaks5[df.cyto.numPeaks5$Channel=="RCaMP",], aes(x=interaction(Condition,ROIType), y= nPeaks, fill=Drug)) + 
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    geom_errorbar(aes(ymin=nPeaks-se, ymax=nPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
    ylab("nPeaks/trial") +
    ggtitle("neurons")+
  scale_fill_manual(values=cbbPalette) + 
    max.theme

#########
all.cyto.DSP4.trials$Drug_Condition=interaction(all.cyto.DSP4.trials$Drug, all.cyto.DSP4.trials$Condition)
all.cyto.DSP4.trials.type$Drug_Condition= interaction(all.cyto.DSP4.trials.type$Drug,
                                                      all.cyto.DSP4.trials.type$Condition)
all.cyto.DSP4.trials.type$Drug_Cond_type= interaction(all.cyto.DSP4.trials.type$Drug,
                                                      all.cyto.DSP4.trials.type$Condition,
                                                      all.cyto.DSP4.trials.type$ROIType)
#before and after dsp4
nP.cytovsDSP2.null = lmer(nPeaks ~ (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials[all.cyto.DSP4.trials$Channel!="RCaMP",],REML=FALSE)
nP.cytovsDSP2.model1 = lmer(nPeaks~ Condition + (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials[all.cyto.DSP4.trials$Channel!="RCaMP",],REML=FALSE)
nP.cytovsDSP2.model2 = lmer(nPeaks ~ Drug+ (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials[all.cyto.DSP4.trials$Channel!="RCaMP",],REML=FALSE)
nP.cytovsDSP2.model3 = lmer(nPeaks ~ Condition+Drug + (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials[all.cyto.DSP4.trials$Channel!="RCaMP",],REML=FALSE)
nP.cytovsDSP2.model4 = lmer(nPeaks ~ Drug_Condition + (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials[all.cyto.DSP4.trials$Channel!="RCaMP",],REML=FALSE)
nP.cytovsDSP2.anova <- anova(nP.cytovsDSP2.null, nP.cytovsDSP2.model1, nP.cytovsDSP2.model2, 
                             nP.cytovsDSP2.model3, nP.cytovsDSP2.model4)
print(nP.cytovsDSP2.anova)

nP.cytovsDSP.treatment_Condition2<- glht(nP.cytovsDSP2.model4, mcp(Drug_Condition= "Tukey"))
summary(nP.cytovsDSP.treatment_Condition2)


nP.cytovsDSP2.null = lmer(nPeaks ~ (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials[all.cyto.DSP4.trials$Channel=="RCaMP",],REML=FALSE)
nP.cytovsDSP2.model1 = lmer(nPeaks~ Condition + (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials[all.cyto.DSP4.trials$Channel=="RCaMP",],REML=FALSE)
nP.cytovsDSP2.model2 = lmer(nPeaks ~ Drug+ (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials[all.cyto.DSP4.trials$Channel=="RCaMP",],REML=FALSE)
nP.cytovsDSP2.model3 = lmer(nPeaks ~ Condition+Drug + (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials[all.cyto.DSP4.trials$Channel=="RCaMP",],REML=FALSE)
nP.cytovsDSP2.model4 = lmer(nPeaks ~ Drug_Condition + (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials[all.cyto.DSP4.trials$Channel=="RCaMP",],REML=FALSE)
nP.cytovsDSP2.anova <- anova(nP.cytovsDSP2.null, nP.cytovsDSP2.model1, nP.cytovsDSP2.model2, 
                             nP.cytovsDSP2.model3, nP.cytovsDSP2.model4)
print(nP.cytovsDSP2.anova)

nP.cytovsDSP.treatment_Condition3<- glht(nP.cytovsDSP2.model4, mcp(Drug_Condition= "Tukey"))
summary(nP.cytovsDSP.treatment_Condition3)


###
# DSP4 nPeaks and ROItype

nP.type.cytovsDSP2.null = lmer(nPeaks ~ (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials.type[all.cyto.DSP4.trials.type$Channel!="RCaMP",],REML=FALSE)
nP.type.cytovsDSP2.model1 = lmer(nPeaks~ Condition + (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials.type[all.cyto.DSP4.trials.type$Channel!="RCaMP",],REML=FALSE)
nP.type.cytovsDSP2.model2 = lmer(nPeaks ~ Drug+ (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials.type[all.cyto.DSP4.trials.type$Channel!="RCaMP",],REML=FALSE)
nP.type.cytovsDSP2.model3 = lmer(nPeaks ~ Condition+Drug + (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials.type[all.cyto.DSP4.trials.type$Channel!="RCaMP",],REML=FALSE)
nP.type.cytovsDSP2.model4 = lmer(nPeaks ~ Drug_Condition + (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials.type[all.cyto.DSP4.trials.type$Channel!="RCaMP",],REML=FALSE)
nP.type.cytovsDSP2.model5 = lmer(nPeaks ~ Condition+Drug+ROIType + (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials.type[all.cyto.DSP4.trials.type$Channel!="RCaMP",],REML=FALSE)
nP.type.cytovsDSP2.model6 = lmer(nPeaks ~ Drug_Cond_type + (1|Animal) + (1|Spot) + (1|SpotTime), all.cyto.DSP4.trials.type[all.cyto.DSP4.trials.type$Channel!="RCaMP",],REML=FALSE)

nP.type.cytovsDSP2.anova <- anova(nP.type.cytovsDSP2.null, nP.type.cytovsDSP2.model1, nP.type.cytovsDSP2.model2, 
                             nP.type.cytovsDSP2.model3, nP.type.cytovsDSP2.model4,
                             nP.type.cytovsDSP2.model5, nP.type.cytovsDSP2.model6)
print(nP.type.cytovsDSP2.anova)

nP.type.cytovsDSP.treatment_Condition2<- glht(nP.type.cytovsDSP2.model4, mcp(Drug_Condition= "Tukey"))
summary(nP.type.cytovsDSP.treatment_Condition2)

nP.type.cytovsDSP.treatment_Cond_type2<- glht(nP.type.cytovsDSP2.model6, mcp(Drug_Cond_type= "Tukey"))
summary(nP.type.cytovsDSP.treatment_Cond_type2)

########
#mean peak amplitude for each ROI with peaks from the stim window
# use this to calculate high, mid and low responding neuronal somas

ROIwise<- ddply(all.lck.peaks, c("Animal", "Spot","Trial","Channel","ROIname","Condition","trials", "ROIs_trial" ,"ROIType"), 
                summarise, meanAmp=mean(amplitude), meanDur= mean(Duration), meanProm=mean(prominence),
                meanPAUC=mean(peakAUC), nPeaks= length(amplitude))

ROIwise$ROInameUnique<-paste(ROIwise$Animal,ROIwise$Spot,ROIwise$ROIname, sep="_")
ROIwise.Neurons<-subset(ROIwise, Condition!="Nostim" & ROIType=="Neuron")
ROIwise.Dendrites<-subset(ROIwise, Condition!="Nostim" & ROIType=="Dendrite")

# MEAN info for each ROI for trials of short OR long stim
NeuronSomas<-ddply(ROIwise.Neurons, c("Animal", "Spot", "ROInameUnique"),
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

all.lck.peaks.group$ROInameUnique<-paste(all.lck.peaks.group$Animal,all.lck.peaks.group$Spot,all.lck.peaks.group$ROIname, sep="_")
#GCaMP_RCaMP.groups$Nresponders=0
for (ii in 1:nrow(NeuronSomas))
{
  Soma=NeuronSomas$ROInameUnique[ii]
  RespType=NeuronSomas$responders[ii]
  subset1=subset(all.lck.peaks.group, ROInameUnique==Soma)
  if (nrow(subset1)>0)
  {
    subset1$Nresponders=RespType
    subset1$ROInameUnique=NULL
    all.lck.peaks.Nresp<-rbind(all.lck.peaks.Nresp, subset1)
  }
}

all.lck.peaks.group$ROInameUnique=NULL
all.lck.peaks.Nresp<-rbind(all.lck.peaks.Nresp, all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",])

