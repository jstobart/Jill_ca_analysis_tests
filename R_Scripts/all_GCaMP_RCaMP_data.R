
library("lme4")
library("lmerTest")
#library("lattice")
library("plyr")
library("ggplot2")
#library("gplots")
library("lsmeans")
library("Rmisc")
#library("MASS")
library("multcomp")
library("reshape2")
library("tidyr")
library("data.table")
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

longstim.lck.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/longstim_onset&AUC.csv", header=TRUE, sep = ",")
shortstim.lck.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/shortstim_onset&AUC.csv", header=TRUE, sep = ",")
nostim.lck.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/nostim_onset&AUC.csv", header=TRUE, sep = ",")

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

# exclude the neuropil ROIs, because they were hand selected and not necessary
all.lck.peaks<-all.lck.peaks[!(all.lck.peaks$ROIname=="np"),]
all.lck.peaks<-all.lck.peaks[!(all.lck.peaks$ROIname=="none"),] # remove trials with NO Peaks

# no stim peak data
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


# remove matching astrocyte process and soma ROIs
Overlap= all.lck.peaks$overlap!=0
all.lck.peaks<-all.lck.peaks[!Overlap,]
#OverlapROIs<-unique(nostim$ROIs_trial[Overlap])

# adjust peak time and duration
all.lck.peaks$peakTime<- all.lck.peaks$peakTime-5
all.lck.peaks$Duration<- all.lck.peaks$halfWidth*2

##### 
# cytoGCaMP6s data

all.cyto.peaks<-rbind(nostim.cyto,long.cyto,short.cyto)
all.cyto.OT<-rbind(nostim.cyto.OT,longstim.cyto.OT,shortstim.cyto.OT)

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


# remove matching astrocyte process and soma ROIs
Overlap= all.cyto.peaks$overlap!=0
all.cyto.peaks<-all.cyto.peaks[!Overlap,]
#OverlapROIs<-unique(nostim$ROIs_trial[Overlap])

# adjust peak time and duration
all.cyto.peaks$peakTime<- all.cyto.peaks$peakTime-5
all.cyto.peaks$Duration<- all.cyto.peaks$halfWidth*2
########
# cytoGCaMP6s DSP4 data

all.cyto.DSP4.peaks<-rbind(nostim.cyto.DSP4,long.cyto.DSP4,short.cyto.DSP4)
all.cyto.DSP4.OT<-rbind(nostim.cyto.DSP4.OT,longstim.cyto.DSP4.OT,shortstim.cyto.DSP4.OT)

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


# remove matching astrocyte process and soma ROIs
Overlap= all.cyto.DSP4.peaks$overlap!=0
all.cyto.DSP4.peaks<-all.cyto.DSP4.peaks[!Overlap,]
#OverlapROIs<-unique(nostim$ROIs_trial[Overlap])

# adjust peak time and duration
all.cyto.DSP4.peaks$peakTime<- all.cyto.DSP4.peaks$peakTime-5
all.cyto.DSP4.peaks$Duration<- all.cyto.DSP4.peaks$halfWidth*2

######

## NOTE: define a stimulation window
# for distributions: histograms

stimwindow=15

stim.lck.OT.dist<-subset(all.lck.OT,Condition!="shortstim" & OnsetTime<stimwindow)
short.lck.OT.dist<-subset(all.lck.OT,Condition!="Stim" & OnsetTime<stimwindow)
stim.cyto.OT.dist<-subset(all.cyto.OT,Condition!="shortstim" & OnsetTime<stimwindow)
short.cyto.OT.dist<-subset(all.cyto.OT,Condition!="Stim" & OnsetTime<stimwindow)
stim.cyto.DSP4.OT.dist<-subset(all.cyto.DSP4.OT,Condition!="shortstim" & OnsetTime<stimwindow)
short.cyto.DSP4.OT.dist<-subset(all.cyto.DSP4.OT,Condition!="Stim" & OnsetTime<stimwindow)

#####
# for median and mean calculations

# no stim vs 8 s stim- neuronal window=9s, AC window= 15 s for peak time, neuronal window=2s, AC window=12 s for onset
# no stim vs 1 s stim- neuronal window=2s, AC window= 10 s for peak time, neuronal window=2s, AC window=8 s for onset

LongN_PTwind=9
LongAC_PTwind=15

LongN_OTwind=2
LongAC_OTwind=12

ShortN_PTwind=2
ShortAC_PTwind=10

ShortN_OTwind=2
ShortAC_OTwind=8

# remove data that is outside the above windows
# lck
stim.lck.OT<-subset(all.lck.OT,Condition!="shortstim")
stim.lck.OT.R<-subset(stim.lck.OT, Channel=="RCaMP" & OnsetTime<=LongN_OTwind)
stim.lck.OT.G<-subset(stim.lck.OT, Channel=="GCaMP" & OnsetTime<=LongAC_OTwind)

stim.lck.OT.window<-rbind(stim.lck.OT.R, stim.lck.OT.G)

short.lck.OT<-subset(all.lck.OT,Condition!="Stim")
short.lck.OT.R<-subset(short.lck.OT, Channel=="RCaMP" & OnsetTime<=ShortN_OTwind)
short.lck.OT.G<-subset(short.lck.OT, Channel=="GCaMP" & OnsetTime<=ShortAC_OTwind)

short.lck.OT.window<-rbind(short.lck.OT.R, short.lck.OT.G)


#cyto
stim.cyto.OT<-subset(all.cyto.OT,Condition!="shortstim")
stim.cyto.OT.R<-subset(stim.cyto.OT, Channel=="RCaMP" & OnsetTime<=LongN_OTwind)
stim.cyto.OT.G<-subset(stim.cyto.OT, Channel=="GCaMP" & OnsetTime<=LongAC_OTwind)

stim.cyto.OT.window<-rbind(stim.cyto.OT.R, stim.cyto.OT.G)

short.cyto.OT<-subset(all.cyto.OT,Condition!="Stim")
short.cyto.OT.R<-subset(short.cyto.OT, Channel=="RCaMP" & OnsetTime<=ShortN_OTwind)
short.cyto.OT.G<-subset(short.cyto.OT, Channel=="GCaMP" & OnsetTime<=ShortAC_OTwind)

short.cyto.OT.window<-rbind(short.cyto.OT.R, short.cyto.OT.G)


#cyto DSP4
stim.cyto.DSP4.OT<-subset(all.cyto.DSP4.OT,Condition!="shortstim")
stim.cyto.DSP4.OT.R<-subset(stim.cyto.DSP4.OT, Channel=="RCaMP" & OnsetTime<=LongN_OTwind)
stim.cyto.DSP4.OT.G<-subset(stim.cyto.DSP4.OT, Channel=="GCaMP" & OnsetTime<=LongAC_OTwind)

stim.cyto.DSP4.OT.window<-rbind(stim.cyto.DSP4.OT.R, stim.cyto.DSP4.OT.G)

short.cyto.DSP4.OT<-subset(all.cyto.DSP4.OT,Condition!="Stim")
short.cyto.DSP4.OT.R<-subset(short.cyto.DSP4.OT, Channel=="RCaMP" & OnsetTime<=ShortN_OTwind)
short.cyto.DSP4.OT.G<-subset(short.cyto.DSP4.OT, Channel=="GCaMP" & OnsetTime<=ShortAC_OTwind)

short.cyto.DSP4.OT.window<-rbind(short.cyto.DSP4.OT.R, short.cyto.DSP4.OT.G)


# peak times
# lck
stim.lck.PT<-subset(all.lck.peaks,Condition!="shortstim")
stim.lck.PT.R<-subset(stim.lck.PT, Channel=="RCaMP" & peakTime<=LongN_PTwind)
stim.lck.PT.G<-subset(stim.lck.PT, Channel=="GCaMP" & peakTime<=LongAC_PTwind)

stim.lck.peaks.window<-rbind(stim.lck.PT.R, stim.lck.PT.G)

short.lck.PT<-subset(all.lck.peaks,Condition!="Stim")
short.lck.PT.R<-subset(short.lck.PT, Channel=="RCaMP" & peakTime<=ShortN_PTwind)
short.lck.PT.G<-subset(short.lck.PT, Channel=="GCaMP" & peakTime<=ShortAC_PTwind)

short.lck.peaks.window<-rbind(short.lck.PT.R, short.lck.PT.G)


#cyto
stim.cyto.PT<-subset(all.cyto.peaks,Condition!="shortstim")
stim.cyto.PT.R<-subset(stim.cyto.PT, Channel=="RCaMP" & peakTime<=LongN_PTwind)
stim.cyto.PT.G<-subset(stim.cyto.PT, Channel=="GCaMP" & peakTime<=LongAC_PTwind)

stim.cyto.peaks.window<-rbind(stim.cyto.PT.R, stim.cyto.PT.G)

short.cyto.PT<-subset(all.cyto.peaks,Condition!="Stim")
short.cyto.PT.R<-subset(short.cyto.PT, Channel=="RCaMP" & peakTime<=ShortN_PTwind)
short.cyto.PT.G<-subset(short.cyto.PT, Channel=="GCaMP" & peakTime<=ShortAC_PTwind)

short.cyto.peaks.window<-rbind(short.cyto.PT.R, short.cyto.PT.G)


#cyto DSP4
stim.cyto.DSP4.PT<-subset(all.cyto.DSP4.peaks,Condition!="shortstim")
stim.cyto.DSP4.PT.R<-subset(stim.cyto.DSP4.PT, Channel=="RCaMP" & peakTime<=LongN_PTwind)
stim.cyto.DSP4.PT.G<-subset(stim.cyto.DSP4.PT, Channel=="GCaMP" & peakTime<=LongAC_PTwind)

stim.cyto.DSP4.peaks.window<-rbind(stim.cyto.DSP4.PT.R, stim.cyto.DSP4.PT.G)

short.cyto.DSP4.PT<-subset(all.cyto.DSP4.peaks,Condition!="Stim")
short.cyto.DSP4.PT.R<-subset(short.cyto.DSP4.PT, Channel=="RCaMP" & peakTime<=ShortN_PTwind)
short.cyto.DSP4.PT.G<-subset(short.cyto.DSP4.PT, Channel=="GCaMP" & peakTime<=ShortAC_PTwind)

short.cyto.DSP4.peaks.window<-rbind(short.cyto.DSP4.PT.R, short.cyto.DSP4.PT.G)


######
# LCK DATA
# Onset time histograms- normalized to the number of trials

# long stim vs no stim
ntrials.lck.OT.long.R.dis<- ddply(stim.lck.OT.dist[stim.lck.OT.dist$Channel=="RCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))
ntrials.lck.OT.long.G.dis<- ddply(stim.lck.OT.dist[stim.lck.OT.dist$Channel=="GCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))


#histogram bins
histseq= seq(0,15,1)
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



# means, medians, modes
# novel function for mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

stim.lck.OT.mean1<-summarySE(stim.lck.OT.dist, measurevar = "OnsetTime", groupvars = c("Channel", "Condition"))
stim.lck.OT.mean2<-summarySE(stim.lck.OT.window, measurevar = "OnsetTime", groupvars = c("Channel", "Condition"))


stim.lck.OT.stats1<-ddply(stim.lck.OT.dist, c("Condition","Channel"), summarise, mean=mean(OnsetTime),
                       median=median(OnsetTime), mode=getmode(OnsetTime))

stim.lck.OT.stats2<-ddply(stim.lck.OT.window, c("Condition","Channel"), summarise, mean=mean(OnsetTime),
                       median=median(OnsetTime), mode=getmode(OnsetTime))

# stats
Condition_Channel=interaction(stim.lck.OT.window$Condition,stim.lck.OT.window$Channel)

# stats for onset times- neurons vs astrocytes
OT.lck.stim.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial), stim.lck.OT.window,REML=FALSE)
OT.lck.stim.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.lck.OT.window,REML=FALSE)
OT.lck.stim.model2 = lmer(OnsetTime ~ Condition + (1|Animal) + (1|Spot) + (1|Spot_trial) , stim.lck.OT.window,REML=FALSE)
OT.lck.stim.model3 = lmer(OnsetTime ~ Condition_Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.lck.OT.window,REML=FALSE)
OT.lck.stim.anova <- anova(OT.lck.stim.null, OT.lck.stim.model1,OT.lck.stim.model2,OT.lck.stim.model3)
print(OT.lck.stim.anova)

OT.lck.stim.Cond_Channel<- glht(OT.lck.stim.model3, mcp(Condition_Channel= "Tukey"))
summary(OT.lck.stim.Cond_Channel)


ggplot(stim.lck.OT.mean2, aes(x=Channel,y=OnsetTime, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme

ggplot(stim.lck.OT.stats2, aes(x=Channel,y=median, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  ylab("Median Onset Time (s)") +
  max.theme
 
######
# number of ROIs per trial per field of view?
ROInum.lck.stim<-ddply(stim.lck.OT.window, c("Animal","Spot","Spot_trial","Condition","Channel"), summarise, nROIs=length(OnsetTime))

stim.lck.ROInum.mean<-summarySE(ROInum.lck.stim, measurevar = "nROIs", groupvars = c("Channel", "Condition"))

ggplot(stim.lck.ROInum.mean, aes(x=Channel,y=nROIs, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=nROIs-se, ymax=nROIs+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view") +
  max.theme

Condition_Channel2= interaction(ROInum.lck.stim$Condition,ROInum.lck.stim$Channel)
nROI.lck.stim.null = lmer(nROIs ~ (1|Animal) + (1|Spot), ROInum.lck.stim,REML=FALSE)
nROI.lck.stim.model1 = lmer(nROIs~ Channel + (1|Animal) + (1|Spot), ROInum.lck.stim,REML=FALSE)
nROI.lck.stim.model2 = lmer(nROIs ~ Condition + (1|Animal) + (1|Spot), ROInum.lck.stim,REML=FALSE)
nROI.lck.stim.model3 = lmer(nROIs ~ Condition_Channel2 + (1|Animal) + (1|Spot), ROInum.lck.stim,REML=FALSE)
nROI.lck.stim.anova <- anova(nROI.lck.stim.null, nROI.lck.stim.model1,nROI.lck.stim.model2,nROI.lck.stim.model3)
print(nROI.lck.stim.anova)

nROI.lck.stim.Cond_Channel<- glht(nROI.lck.stim.model3, mcp(Condition_Channel2= "Tukey"))
summary(nROI.lck.stim.Cond_Channel)

######
#peak times

stimwindow=15
histseq= seq(-5,15,1)

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
  geom_line(data=lck.long.PT.histo2, aes(y=Nostim.A*2, color="Nostim.A")) +
  geom_line(data=lck.long.PT.histo2, aes(y=Stim.A*2, color="Stim.A")) +
  scale_y_continuous(sec.axis = sec_axis(~./2, name = "Astrocyte peaks/trial")) + # secondary axis
  ggtitle("stim- neurons vs astrocytes-lck data") + 
  xlab("Peak Max Time (s)") + 
  ylab("Neuron peaks/trial") + 
  xlim(-0.5,15) +
  max.theme



##################### 
#cyto data
ntrials.cyto.OT.long.R.dis<- ddply(stim.cyto.OT.window[stim.cyto.OT.window$Channel=="RCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))
ntrials.cyto.OT.long.G.dis<- ddply(stim.cyto.OT.window[stim.cyto.OT.window$Channel=="GCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))

ntrials.cyto.DSP4.OT.long.R<- ddply(stim.cyto.DSP4.OT.window[stim.cyto.DSP4.OT.window$Channel=="RCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))
ntrials.cyto.DSP4.OT.long.G<- ddply(stim.cyto.DSP4.OT.window[stim.cyto.DSP4.OT.window$Channel=="GCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))

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



# means, medians, modes

stim.cyto.OT.mean1<-summarySE(stim.cyto.OT.dist, measurevar = "OnsetTime", groupvars = c("Channel", "Condition"))
stim.cyto.OT.mean2<-summarySE(stim.cyto.OT.window, measurevar = "OnsetTime", groupvars = c("Channel", "Condition"))


stim.cyto.OT.stats1<-ddply(stim.cyto.OT.dist, c("Condition","Channel"), summarise, mean=mean(OnsetTime),
                          median=median(OnsetTime), mode=getmode(OnsetTime))

stim.cyto.OT.stats2<-ddply(stim.cyto.OT.window, c("Condition","Channel"), summarise, mean=mean(OnsetTime),
                          median=median(OnsetTime), mode=getmode(OnsetTime))

# stats
Condition_Channel=interaction(stim.cyto.OT.window$Condition,stim.cyto.OT.window$Channel)

# stats for onset times- neurons vs astrocytes
OT.cyto.stim.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial), stim.cyto.OT.window,REML=FALSE)
OT.cyto.stim.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.cyto.OT.window,REML=FALSE)
OT.cyto.stim.model2 = lmer(OnsetTime ~ Condition + (1|Animal) + (1|Spot) + (1|Spot_trial) , stim.cyto.OT.window,REML=FALSE)
OT.cyto.stim.model3 = lmer(OnsetTime ~ Condition_Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.cyto.OT.window,REML=FALSE)
OT.cyto.stim.anova <- anova(OT.cyto.stim.null, OT.cyto.stim.model1,OT.cyto.stim.model2,OT.cyto.stim.model3)
print(OT.cyto.stim.anova)

OT.cyto.stim.Cond_Channel<- glht(OT.cyto.stim.model3, mcp(Condition_Channel= "Tukey"))
summary(OT.cyto.stim.Cond_Channel)


ggplot(stim.cyto.OT.mean2, aes(x=Channel,y=OnsetTime, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme

ggplot(stim.cyto.OT.stats2, aes(x=Channel,y=median, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  ylab("Median Onset Time (s)") +
  max.theme

 ##############
# number of ROIs per trial per field of view?
ROInum.cyto.stim<-ddply(stim.cyto.OT.window, c("Animal","Spot","Spot_trial","Condition","Channel"), summarise, nROIs=length(OnsetTime))

stim.cyto.ROInum.mean<-summarySE(ROInum.cyto.stim, measurevar = "nROIs", groupvars = c("Channel", "Condition"))

ggplot(stim.cyto.ROInum.mean, aes(x=Channel,y=nROIs, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=nROIs-se, ymax=nROIs+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view") +
  max.theme

Condition_Channel2= interaction(ROInum.cyto.stim$Condition,ROInum.cyto.stim$Channel)
nROI.cyto.stim.null = lmer(nROIs ~ (1|Animal) + (1|Spot), ROInum.cyto.stim,REML=FALSE)
nROI.cyto.stim.model1 = lmer(nROIs~ Channel + (1|Animal) + (1|Spot), ROInum.cyto.stim,REML=FALSE)
nROI.cyto.stim.model2 = lmer(nROIs ~ Condition + (1|Animal) + (1|Spot), ROInum.cyto.stim,REML=FALSE)
nROI.cyto.stim.model3 = lmer(nROIs ~ Condition_Channel2 + (1|Animal) + (1|Spot), ROInum.cyto.stim,REML=FALSE)
nROI.cyto.stim.anova <- anova(nROI.cyto.stim.null, nROI.cyto.stim.model1,nROI.cyto.stim.model2,nROI.cyto.stim.model3)
print(nROI.cyto.stim.anova)

nROI.cyto.stim.Cond_Channel<- glht(nROI.cyto.stim.model3, mcp(Condition_Channel2= "Tukey"))
summary(nROI.cyto.stim.Cond_Channel)





#####
# neurons and astrocytes together by density to account for different number of ROIs

ggplot(NULL, aes(x=OnsetTime))+
  geom_freqpoly(data= stim.lck.OT.window[stim.lck.OT.window$Channel=="RCaMP",], binwidth = 0.5, lwd=1,
                                         aes(y=..density.., colour = Channel, linetype=Condition)) +
  geom_freqpoly(data= stim.lck.OT.window[stim.lck.OT.window$Channel=="GCaMP",], binwidth = 0.5, lwd=1,
                                         aes(y=(..density..)*5, colour = Channel, linetype=Condition)) +
  scale_y_continuous(sec.axis = sec_axis(~./5, name = "Astrocyte density")) + # secondary axis
  ggtitle("stim- neurons vs astrocytes-lck data") + 
  xlab("Onset Time (s)") + 
  ylab("Neuron density") + 
  xlim(-0.5,LongAC_OTwind) +
  max.theme  
  

#####

# zoomed in plot of onset times- identify "FAST" astrocytes
stim.lck.OT.dist<-subset(stim.lck.OT.dist,Condition!="Nostim")
stim.lck.OT.dist$Group<-0
stim.lck.OT.dist$Group[stim.lck.OT.dist$OnsetTime<1]<-"fast"
stim.lck.OT.dist$Group[stim.lck.OT.dist$OnsetTime>=1]<-"delayed"

ggplot(stim.lck.OT.dist[stim.lck.OT.dist$Group=="fast"& stim.lck.OT.dist$Channel=="GCaMP",], aes(x=OnsetTime, y=..density.., colour = interaction(ROIType,Group), linetype=Condition)) +
  geom_freqpoly(binwidth = (0.0845*3), lwd=1)+
  ggtitle("EF vs P-lck data") + 
  xlab("Onset Time (s)") + 
  #scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

ggplot(stim.lck.OT[stim.lck.OT$Group=="fast",], aes(x=OnsetTime, y=..density.., colour = interaction(Channel,Group), linetype=Condition)) +
  geom_freqpoly(binwidth = 0.0845, lwd=1)+
  ggtitle("stimwindow- neurons vs astrocytes-lck data") + 
  xlab("Onset Time (s)") + 
  xlim(0,1)+
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

ggplot(stim.lck.OT[stim.lck.OT$Group=="fast",], aes(x=OnsetTime, y=..density.., colour = interaction(Channel,Group))) +
  geom_freqpoly(binwidth = 0.0845, lwd=1)+
  ggtitle("stimwindow- neurons vs astrocytes-lck data") + 
  xlab("Onset Time (s)") + 
  xlim(0,1)+
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

ggplot(stim.lck.OT[stim.lck.OT$Group=="fast",], aes(x=OnsetTime, colour = interaction(Channel,Group), linetype=Condition)) + 
  stat_ecdf(lwd=1) +
  ggtitle("cdf- neurons vs astrocytes-fast ROIs") + 
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

# mean onset times
df.OT1<- summarySE(stim.lck.OT, measurevar = "OnsetTime", groupvars = c("Channel", "Group"))
df.OT1B<- summarySE(stim.lck.OT, measurevar = "OnsetTime", groupvars = c("ROIType", "Group"))
df.OT1C<- summarySE(all.lck.OT, measurevar = "OnsetTime", groupvars = c("ROIType", "Condition"))

#df.OT1$OnsetTime<-df.OT1$OnsetTime*1000  # convert to ms
#df.OT1$se<-df.OT1$se*1000 # convert to ms

df.OT1$Channel <- factor(df.OT1$Channel, levels = c("RCaMP","GCaMP"))
#df.OT1$Condition <- factor(df.OT1$Condition, levels = c("shortstim","Stim"))
df.OT1$Group <- factor(df.OT1$Group, levels = c("fast","delayed"))

df.OT1 = df.OT1[!(df.OT1$Channel=="RCaMP"&df.OT1$Group=="delayed"),]

ggplot(df.OT1, aes(x=interaction(Channel,Group),y=OnsetTime, fill= interaction(Channel,Group))) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

stim.lck.OT$Group<-as.factor(stim.lck.OT$Group)
Condition_Channel=interaction(stim.lck.OT$Condition,stim.lck.OT$Channel)
Group_Channel=interaction(stim.lck.OT$Group,stim.lck.OT$Channel)
Group_Channel_Type=interaction(stim.lck.OT$Group,stim.lck.OT$Channel,stim.lck.OT$ROIType)
# stats for onset times- neurons vs astrocytes
OT.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial), stim.lck.OT,REML=FALSE)
OT.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.lck.OT,REML=FALSE)
OT.model2 = lmer(OnsetTime ~ Condition + (1|Animal) + (1|Spot) + (1|Spot_trial) , stim.lck.OT,REML=FALSE)
OT.model3 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.lck.OT,REML=FALSE)
OT.model4 = lmer(OnsetTime ~ Group_Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.lck.OT,REML=FALSE)
OT.model5 = lmer(OnsetTime ~ Condition_Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.lck.OT,REML=FALSE)
OT.model6 = lmer(OnsetTime ~ Group_Channel_Type + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.lck.OT,REML=FALSE)
OT.anova <- anova(OT.null, OT.model1,OT.model2,OT.model3,OT.model4,OT.model5,OT.model6)
print(OT.anova)

OT.Group_channel<- glht(OT.model4, mcp(Group_Channel= "Tukey"))
summary(OT.Group_channel)

OT.Group_channel_ty<- glht(OT.model6, mcp(Group_Channel_Type= "Tukey"))
summary(OT.Group_channel_ty)


Condition_ROIType=interaction(all.lck.OT$Condition,all.lck.OT$ROIType)
# stats for onset times- neurons vs astrocytes
OT.lck.type.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial), all.lck.OT,REML=FALSE)
OT.lck.type.model1 = lmer(OnsetTime ~ ROIType + (1|Animal) + (1|Spot) + (1|Spot_trial), all.lck.OT,REML=FALSE)
OT.lck.type.model2 = lmer(OnsetTime ~ Condition + (1|Animal) + (1|Spot) + (1|Spot_trial), all.lck.OT,REML=FALSE)
OT.lck.type.model4 = lmer(OnsetTime ~ ROIType+Condition + (1|Animal) + (1|Spot) + (1|Spot_trial), all.lck.OT,REML=FALSE)
OT.lck.type.model5 = lmer(OnsetTime ~ Condition_ROIType + (1|Animal) + (1|Spot) + (1|Spot_trial), all.lck.OT,REML=FALSE)
OT.lck.type.anova <- anova(OT.lck.type.null, OT.lck.type.model1, OT.lck.type.model2, 
                           OT.lck.type.model4, OT.lck.type.model5)
print(OT.lck.type.anova)

OT.lck.Cond_ROIType<- glht(OT.lck.type.model5, mcp(Condition_ROIType= "Tukey"))
summary(OT.lck.Cond_ROIType)


#AUC vs onset time scatter plot
ggplot(stim.lck.OT, aes(x=OnsetTime,y=TraceAUC10, colour=Group)) +
  geom_point() +
  scale_colour_manual(values=cbbPalette)+
  max.theme



###########################

# mean onset times
df.OT2<- summarySE(stim.cyto.OT, measurevar = "OnsetTime", groupvars = c("Channel", "Group"))
df.OT2C<- summarySE(all.cyto.OT, measurevar = "OnsetTime", groupvars = c("ROIType", "Condition"))

#df.OT1$OnsetTime<-df.OT1$OnsetTime*1000  # convert to ms
#df.OT1$se<-df.OT1$se*1000 # convert to ms

df.OT2$Channel <- factor(df.OT2$Channel, levels = c("RCaMP","GCaMP"))
#df.OT2$Condition <- factor(df.OT2$Condition, levels = c("shortstim","Stim"))
df.OT2$Group <- factor(df.OT2$Group, levels = c("fast","delayed"))

df.OT2 = df.OT2[!(df.OT2$Channel=="RCaMP"&df.OT2$Group=="delayed"),]

ggplot(df.OT2, aes(x=interaction(Channel,Group),y=OnsetTime, fill= interaction(Channel,Group))) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

stim.cyto.OT$Group<-as.factor(stim.cyto.OT$Group)
Condition_Channel=interaction(stim.cyto.OT$Condition,stim.cyto.OT$Channel)
Group_Channel=interaction(stim.cyto.OT$Group,stim.cyto.OT$Channel)
# stats for onset times- neurons vs astrocytes
OT.cyto.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial), stim.cyto.OT,REML=FALSE)
OT.cyto.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.cyto.OT,REML=FALSE)
OT.cyto.model2 = lmer(OnsetTime ~ Condition + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.cyto.OT,REML=FALSE)
OT.cyto.model3 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.cyto.OT,REML=FALSE)
OT.cyto.model4 = lmer(OnsetTime ~ Group_Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.cyto.OT,REML=FALSE)
OT.cyto.model5 = lmer(OnsetTime ~ Condition_Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.cyto.OT,REML=FALSE)
OT.cyto.anova <- anova(OT.cyto.null, OT.cyto.model1, OT.cyto.model2, OT.cyto.model3, OT.cyto.model4, OT.cyto.model5)
print(OT.cyto.anova)

OT.cyto.Group_channel<- glht(OT.cyto.model4, mcp(Group_Channel= "Tukey"))
summary(OT.cyto.Group_channel)


Condition_ROIType=interaction(all.cyto.OT$Condition,all.cyto.OT$ROIType)
# stats for onset times- neurons vs astrocytes
OT.cyto.type.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial), all.cyto.OT,REML=FALSE)
OT.cyto.type.model1 = lmer(OnsetTime ~ ROIType + (1|Animal) + (1|Spot) + (1|Spot_trial), all.cyto.OT,REML=FALSE)
OT.cyto.type.model2 = lmer(OnsetTime ~ Condition + (1|Animal) + (1|Spot) + (1|Spot_trial), all.cyto.OT,REML=FALSE)
OT.cyto.type.model4 = lmer(OnsetTime ~ ROIType+Condition + (1|Animal) + (1|Spot) + (1|Spot_trial), all.cyto.OT,REML=FALSE)
OT.cyto.type.model5 = lmer(OnsetTime ~ Condition_ROIType + (1|Animal) + (1|Spot) + (1|Spot_trial), all.cyto.OT,REML=FALSE)
OT.cyto.type.anova <- anova(OT.cyto.type.null, OT.cyto.type.model1, OT.cyto.type.model2, 
                            OT.cyto.type.model4, OT.cyto.type.model5)
print(OT.cyto.type.anova)

OT.cyto.Cond_ROIType<- glht(OT.cyto.type.model5, mcp(Condition_ROIType= "Tukey"))
summary(OT.cyto.Cond_ROIType)


######
# cyto DSP4 data

#cyto data
# neuronal onset histogram
#histogram bins
histseq= seq(0,stimwindow,1)
# counts for each condition in the histogram
Nostim=hist(all.cyto.DSP4.OT$OnsetTime[(all.cyto.DSP4.OT$Channel=="RCaMP" & all.cyto.DSP4.OT$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim=hist(all.cyto.DSP4.OT$OnsetTime[(all.cyto.DSP4.OT$Channel=="RCaMP" & all.cyto.DSP4.OT$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts
Shortstim=hist(all.cyto.DSP4.OT$OnsetTime[(all.cyto.DSP4.OT$Channel=="RCaMP" & all.cyto.DSP4.OT$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim=Nostim/ntrials.cyto.DSP4.OT$ntrials[(ntrials.cyto.DSP4.OT$Condition=="Nostim")]
Stim=Stim/ntrials.cyto.DSP4.OT$ntrials[(ntrials.cyto.DSP4.OT$Condition=="Stim")]
Shortstim=Shortstim/ntrials.cyto.DSP4.OT$ntrials[(ntrials.cyto.DSP4.OT$Condition=="shortstim")]

#make a data frame for plotting
Neuron.cyto.DSP4.histo <- data.frame(cbind(Nostim, Stim, Shortstim))
Neuron.cyto.DSP4.histo$time<-histseq[2:length(histseq)]
Neuron.cyto.DSP4.histo<-melt(Neuron.cyto.DSP4.histo,id="time")

ggplot(data=Neuron.cyto.DSP4.histo, aes(x=time, y= value, colour=variable)) + 
  geom_line(size=1)+
  xlab("Onset time (s)") +
  ylab("number of ROIs per trial")+
  ggtitle("Neurons- onset time- ROIs per trial- cyto.DSP4 data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

# astrocyte lck onset histogram
#histogram bins
histseq= seq(0,stimwindow,1)

# counts for each condition in the histogram
Nostim=hist(all.cyto.DSP4.OT$OnsetTime[(all.cyto.DSP4.OT$Channel=="GCaMP" & all.cyto.DSP4.OT$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim=hist(all.cyto.DSP4.OT$OnsetTime[(all.cyto.DSP4.OT$Channel=="GCaMP" & all.cyto.DSP4.OT$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts
Shortstim=hist(all.cyto.DSP4.OT$OnsetTime[(all.cyto.DSP4.OT$Channel=="GCaMP" & all.cyto.DSP4.OT$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim=Nostim/ntrials.cyto.DSP4.OT$ntrials[(ntrials.cyto.DSP4.OT$Condition=="Nostim")]
Stim=Stim/ntrials.cyto.DSP4.OT$ntrials[(ntrials.cyto.DSP4.OT$Condition=="Stim")]
Shortstim=Shortstim/ntrials.cyto.DSP4.OT$ntrials[(ntrials.cyto.DSP4.OT$Condition=="shortstim")]

#make a data frame for plotting
AC.cyto.DSP4.histo<- data.frame(cbind(Nostim, Stim, Shortstim))
AC.cyto.DSP4.histo$time<-histseq[2:length(histseq)]
AC.cyto.DSP4.histo<-melt(AC.cyto.DSP4.histo,id="time")

ggplot(data=AC.cyto.DSP4.histo, aes(x=time, y= value, colour=variable)) + 
  geom_line(size=1)+
  xlab("Onset time (s)") +
  ylab("number of ROIs per trial")+
  ggtitle("Astrocytes- onset time- ROIs per trial- cyto.DSP4 data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

###
# neurons and astrocytes together by density to account for different number of ROIs
ggplot(all.cyto.DSP4.OT[(all.cyto.DSP4.OT$Condition!="Nostim"),], aes(x=OnsetTime, y=..density.., colour = Channel, linetype=Condition)) +
  geom_freqpoly(binwidth = 0.5, lwd=1)+
  ggtitle("stim- neurons vs astrocytes-cyto.DSP4 data") + 
  xlab("Onset Time (s)") + 
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

ggplot(all.cyto.DSP4.OT[(all.cyto.DSP4.OT$Condition!="Nostim"),], aes(x=OnsetTime, colour = Channel, linetype=Condition)) + 
  stat_ecdf() +
  ggtitle("cdf- neurons vs astrocytes-cyto.DSP4 data- stim types") + 
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

ggplot(all.cyto.DSP4.OT[(all.cyto.DSP4.OT$Condition!="Nostim"),], aes(x=OnsetTime, colour = Channel)) + 
  stat_ecdf(lwd=1) +
  ggtitle("cdf- neurons vs astrocytes-cyto.DSP4 data") + 
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme



########
# COMPARE CONTROL VS DSP4 Onset times
all.cyto.DSP4.OT$treatment="DSP4"
all.cyto.OT$treatment="control"

control.vs.DSP4.OT<-rbind(all.cyto.DSP4.OT,all.cyto.OT)

ggplot(control.vs.DSP4.OT[(control.vs.DSP4.OT$Channel!="RCaMP"),], aes(x=OnsetTime, colour = Condition, linetype=treatment)) + 
  stat_ecdf(lwd=1) +
  ggtitle("cdf- astrocytes-control vs DSP4- stim types") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

ggplot(control.vs.DSP4.OT[(control.vs.DSP4.OT$Channel!="GCaMP"),], aes(x=OnsetTime, colour = Condition, linetype=treatment)) + 
  stat_ecdf(lwd=1) +
  ggtitle("cdf- neurons-control vs DSP4- stim types") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

# mean onset times
df.OT3<- summarySE(control.vs.DSP4.OT, measurevar = "OnsetTime", groupvars = c("Channel", "treatment"))

df.OT3$Channel <- factor(df.OT3$Channel, levels = c("RCaMP","GCaMP"))
df.OT3$treatment <- factor(df.OT3$treatment, levels = c("fast","delayed"))

#df.OT1 = df.OT1[!(df.OT1$Channel=="RCaMP"&df.OT1$Group=="delayed"),]

ggplot(df.OT3, aes(x=Channel,y=OnsetTime, fill= interaction(Channel,treatment))) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme



#######
# peak time histograms

# NOTE: sometimes a peak is detected for the whole trial so exclude peaks with really long durations!

all.lck.peaks<-subset(all.lck.peaks,peakTime<stimwindow & peakTime>0 & Duration<45)
all.cyto.peaks<-subset(all.cyto.peaks,peakTime<stimwindow & peakTime>0 & Duration<80)
all.cyto.DSP4.peaks<-subset(all.cyto.DSP4.peaks,peakTime<stimwindow & peakTime>0 & Duration<80)

ntrials.lck.peaks<- ddply(all.lck.peaks, c("Condition"), summarise, ntrials=length(unique(trials)))
ntrials.cyto.peaks<- ddply(all.cyto.peaks, c("Condition"), summarise, ntrials=length(unique(trials)))
ntrials.cyto.DSP4.peaks<- ddply(all.cyto.DSP4.peaks, c("Condition"), summarise, ntrials=length(unique(trials)))


#####
# lck peak time data

#histogram bins
histseq= seq(0,stimwindow, 1)

# neuronal peak time histogram
# counts for each condition in the histogram
Nostim=hist(all.lck.peaks$peakTime[(all.lck.peaks$Channel=="RCaMP" & all.lck.peaks$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim=hist(all.lck.peaks$peakTime[(all.lck.peaks$Channel=="RCaMP" & all.lck.peaks$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts
Shortstim=hist(all.lck.peaks$peakTime[(all.lck.peaks$Channel=="RCaMP" & all.lck.peaks$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim=Nostim/ntrials.lck.peaks$ntrials[(ntrials.lck.peaks$Condition=="Nostim")]
Stim=Stim/ntrials.lck.peaks$ntrials[(ntrials.lck.peaks$Condition=="Stim")]
Shortstim=Shortstim/ntrials.lck.peaks$ntrials[(ntrials.lck.peaks$Condition=="shortstim")]

#make a data frame for plotting
Neuron.lck.histo <- data.frame(cbind(Nostim, Stim, Shortstim))
Neuron.lck.histo$time<-histseq[2:length(histseq)]
Neuron.lck.histo<-melt(Neuron.lck.histo,id="time")

ggplot(data=Neuron.lck.histo, aes(x=time, y= value, colour=variable)) + 
  geom_line(size=1)+
  xlab("Time to Peak Max (s)") +
  ylab("number of peaks per trial")+
  ggtitle("Neurons- peak time- ROIs per trial- lck data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

# astrocyte lck onset histogram
#histogram bins
histseq= seq(0,stimwindow, 1)

# counts for each condition in the histogram
Nostim=hist(all.lck.peaks$peakTime[(all.lck.peaks$Channel=="GCaMP" & all.lck.peaks$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim=hist(all.lck.peaks$peakTime[(all.lck.peaks$Channel=="GCaMP" & all.lck.peaks$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts
Shortstim=hist(all.lck.peaks$peakTime[(all.lck.peaks$Channel=="GCaMP" & all.lck.peaks$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim=Nostim/ntrials.lck.peaks$ntrials[(ntrials.lck.peaks$Condition=="Nostim")]
Stim=Stim/ntrials.lck.peaks$ntrials[(ntrials.lck.peaks$Condition=="Stim")]
Shortstim=Shortstim/ntrials.lck.peaks$ntrials[(ntrials.lck.peaks$Condition=="shortstim")]

#make a data frame for plotting
AC.lck.histo<- data.frame(cbind(Nostim, Stim, Shortstim))
AC.lck.histo$time<-histseq[2:length(histseq)]
AC.lck.histo<-melt(AC.lck.histo,id="time")

ggplot(data=AC.lck.histo, aes(x=time, y= value, colour=variable)) + 
  geom_line(size=1)+
  xlab("Time to Peak Max (s)") +
  xlim(0,stimwindow) +
  ylab("number of peaks per trial")+
  ggtitle("Astrocytes- peak time- ROIs per trial- lck data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme


#DEFINE FAST and DELAYED
all.lck.peaks$Group<-"NaN"
all.lck.peaksA=subset(all.lck.peaks, Condition!="Nostim" & ROIs_trial %in% stim.lck.OT$ROIs_trial[stim.lck.OT$OnsetTime<1])
all.lck.peaksA$Group="fast"
all.lck.peaksB=subset(all.lck.peaks, Condition!="Nostim" & ROIs_trial %in% stim.lck.OT$ROIs_trial[stim.lck.OT$OnsetTime>1])
all.lck.peaksB$Group="delayed"
all.lck.peaks.group<-rbind(all.lck.peaksA, all.lck.peaksB)


ggplot(all.lck.peaks.group, aes(x=peakTime, colour = interaction(Channel,Group)))+ #, linetype=Condition)) +
  stat_ecdf(lwd=1) +
  ggtitle("cdf- lck- neurons vs astrocytes-fast ROIs") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

# mean peak times
df.pT1A<- summarySE(all.lck.peaks.group, measurevar = "peakTime", groupvars = c("Channel", "Group","Condition"))
df.pT1B<- summarySE(all.lck.peaks, measurevar = "peakTime", groupvars = c("ROIType", "Condition"))

ggplot(all.lck.peaks[(all.lck.peaks$Condition!="shortstim"),], aes(x=peakTime, y=..density.., fill = Channel)) +
  geom_histogram(binwidth = 0.5, position=position_dodge())+
  ggtitle("long stim- neurons vs astrocytes-lck data") + 
  xlab("Peak Time (s)") + 
  xlim(-0.5,20)+
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme



#print median values
all.lck.peaks.short<-subset(all.lck.peaks, peakTime<15)

lck.pt.median<-quantile(all.lck.peaks.short$peakTime[(all.lck.peaks.short$Condition!="shortstim"&all.lck.peaks.short$Channel=="GCaMP")])
lck.pt.median=lck.pt.median[3]

rcamp1.pt.median<-quantile(all.lck.peaks.short$peakTime[(all.lck.peaks.short$Condition!="shortstim"&all.lck.peaks.short$Channel=="RCaMP")],
                           prob = seq(0, 1, length = 5), type = 5)
rcamp1.pt.median<-rcamp1.pt.median[3]
#summary(all.lck.peaks.short$peakTime[(all.lck.peaks.short$Condition!="shortstim"&all.lck.peaks.short$Channel=="GCaMP")])
#summary(all.lck.peaks.short$peakTime[(all.lck.peaks.short$Condition!="shortstim"&all.lck.peaks.short$Channel=="RCaMP")])

ggplot(all.lck.peaks.short[(all.lck.peaks.short$Condition!="shortstim"),], aes(x=peakTime, y=..density.., fill = Channel)) +
  geom_histogram(binwidth = 0.5, position=position_dodge())+
  ggtitle("long stim- neurons vs astrocytes-lck data") + 
  geom_vline(xintercept=lck.pt.median, linetype="dashed", size=1) +
  geom_vline(xintercept=rcamp1.pt.median, linetype="dashed", size=1) +
  xlab("Peak Time (s)") + 
  xlim(-0.5,15)+
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme



Condition_ROIType=interaction(all.lck.peaks$Condition,all.lck.peaks$ROIType)
# stats for onset times- neurons vs astrocytes
pT.lck.type.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), all.lck.peaks,REML=FALSE)
pT.lck.type.model1 = lmer(peakTime ~ ROIType + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.lck.peaks,REML=FALSE)
pT.lck.type.model2 = lmer(peakTime ~ Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.lck.peaks,REML=FALSE)
pT.lck.type.model4 = lmer(peakTime ~ ROIType+Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.lck.peaks,REML=FALSE)
pT.lck.type.model5 = lmer(peakTime ~ Condition_ROIType + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.lck.peaks,REML=FALSE)
pT.lck.type.anova <- anova(pT.lck.type.null, pT.lck.type.model1, pT.lck.type.model2, 
                            pT.lck.type.model4, pT.lck.type.model5)
print(pT.lck.type.anova)

pT.lck.Cond_ROIType<- glht(pT.lck.type.model5, mcp(Condition_ROIType= "Tukey"))
summary(pT.lck.Cond_ROIType)



######
# cyto peak time data

#histogram bins
histseq= seq(0,stimwindow, 1)

# neuronal peak time histogram
# counts for each condition in the histogram
Nostim=hist(all.cyto.peaks$peakTime[(all.cyto.peaks$Channel=="RCaMP" & all.cyto.peaks$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim=hist(all.cyto.peaks$peakTime[(all.cyto.peaks$Channel=="RCaMP" & all.cyto.peaks$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts
Shortstim=hist(all.cyto.peaks$peakTime[(all.cyto.peaks$Channel=="RCaMP" & all.cyto.peaks$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim=Nostim/ntrials.cyto.peaks$ntrials[(ntrials.cyto.peaks$Condition=="Nostim")]
Stim=Stim/ntrials.cyto.peaks$ntrials[(ntrials.cyto.peaks$Condition=="Stim")]
Shortstim=Shortstim/ntrials.cyto.peaks$ntrials[(ntrials.cyto.peaks$Condition=="shortstim")]

#make a data frame for plotting
Neuron.cyto.histo <- data.frame(cbind(Nostim, Stim, Shortstim))
Neuron.cyto.histo$time<-histseq[2:length(histseq)]
Neuron.cyto.histo<-melt(Neuron.cyto.histo,id="time")

ggplot(data=Neuron.cyto.histo, aes(x=time, y= value, colour=variable)) + 
  geom_line(size=1)+
  xlab("Time to Peak Max (s)") +
  ylab("number of peaks per trial")+
  ggtitle("Neurons- peak time- ROIs per trial- cyto data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

# astrocyte lck onset histogram
#histogram bins
histseq= seq(0,stimwindow, 1)

# counts for each condition in the histogram
Nostim=hist(all.cyto.peaks$peakTime[(all.cyto.peaks$Channel=="GCaMP" & all.cyto.peaks$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim=hist(all.cyto.peaks$peakTime[(all.cyto.peaks$Channel=="GCaMP" & all.cyto.peaks$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts
Shortstim=hist(all.cyto.peaks$peakTime[(all.cyto.peaks$Channel=="GCaMP" & all.cyto.peaks$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim=Nostim/ntrials.cyto.peaks$ntrials[(ntrials.cyto.peaks$Condition=="Nostim")]
Stim=Stim/ntrials.cyto.peaks$ntrials[(ntrials.cyto.peaks$Condition=="Stim")]
Shortstim=Shortstim/ntrials.cyto.peaks$ntrials[(ntrials.cyto.peaks$Condition=="shortstim")]

#make a data frame for plotting
AC.cyto.histo<- data.frame(cbind(Nostim, Stim, Shortstim))
AC.cyto.histo$time<-histseq[2:length(histseq)]
AC.cyto.histo<-melt(AC.cyto.histo,id="time")

ggplot(data=AC.cyto.histo, aes(x=time, y= value, colour=variable)) + 
  geom_line(size=1)+
  xlab("Peak time (s)") +
  ylab("number of peaks per trial")+
  ggtitle("Astrocytes- peak time- ROIs per trial- cyto data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

ggplot(all.cyto.peaks[(all.cyto.peaks$Condition=="Stim"),], aes(x=peakTime, y=..density.., fill = Channel)) +
  geom_histogram(binwidth = 0.5, position=position_dodge())+
  ggtitle("long stim- neurons vs astrocytes-cyto data") + 
  xlab("Peak Time (s)") + 
  xlim(-0.5,20)+
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme


df.pT2A<- summarySE(all.cyto.peaks, measurevar = "peakTime", groupvars = c("Condition"))
df.pT2B<- summarySE(all.cyto.peaks, measurevar = "peakTime", groupvars = c("ROIType", "Condition"))

Condition_ROIType=interaction(all.cyto.peaks$Condition,all.cyto.peaks$ROIType)
# stats for onset times- neurons vs astrocytes
pT.cyto.type.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), all.cyto.peaks,REML=FALSE)
pT.cyto.type.model1 = lmer(peakTime ~ ROIType + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.cyto.peaks,REML=FALSE)
pT.cyto.type.model2 = lmer(peakTime ~ Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.cyto.peaks,REML=FALSE)
pT.cyto.type.model4 = lmer(peakTime ~ ROIType+Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.cyto.peaks,REML=FALSE)
pT.cyto.type.model5 = lmer(peakTime ~ Condition_ROIType + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.cyto.peaks,REML=FALSE)
pT.cyto.type.anova <- anova(pT.cyto.type.null, pT.cyto.type.model1, pT.cyto.type.model2, 
                            pT.cyto.type.model4, pT.cyto.type.model5)
print(pT.cyto.type.anova)

pT.cyto.Cond_ROIType<- glht(pT.cyto.type.model5, mcp(Condition_ROIType= "Tukey"))
summary(pT.cyto.Cond_ROIType)


ggplot(all.cyto.peaks[(all.cyto.peaks$Condition!="shortstim"),], aes(x=peakTime, y=..density.., fill = Channel)) +
  geom_histogram(binwidth = 0.5, position=position_dodge())+
  ggtitle("long stim- neurons vs astrocytes-cyto data") + 
  xlab("Peak Time (s)") + 
  xlim(-0.5,15)+
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

ggplot(all.cyto.peaks[(all.cyto.peaks$Condition!="shortstim"),], aes(x=peakTime, y=..density.., colour = ROIType)) +
  geom_freqpoly(binwidth = 0.5, lwd=1)+
  ggtitle("long stim- neurons vs astrocytes-cyto data") + 
  xlab("Peak Time (s)") + 
  xlim(-0.5,15)+
  max.theme

#print median values
all.cyto.peaks.short<-subset(all.cyto.peaks, peakTime<15)

cyto.pt.median<-quantile(all.cyto.peaks.short$peakTime[(all.cyto.peaks.short$Condition!="shortstim"&all.cyto.peaks.short$Channel=="GCaMP")])
cyto.pt.median=cyto.pt.median[3]

rcamp2.pt.median<-quantile(all.cyto.peaks.short$peakTime[(all.cyto.peaks.short$Condition!="shortstim"&all.cyto.peaks.short$Channel=="RCaMP")],
                           prob = seq(0, 1, length = 5), type = 5)
rcamp2.pt.median<-rcamp2.pt.median[3]
#summary(all.lck.peaks.short$peakTime[(all.lck.peaks.short$Condition!="shortstim"&all.lck.peaks.short$Channel=="GCaMP")])
#summary(all.lck.peaks.short$peakTime[(all.lck.peaks.short$Condition!="shortstim"&all.lck.peaks.short$Channel=="RCaMP")])

ggplot(all.cyto.peaks.short[(all.cyto.peaks.short$Condition!="shortstim"),], aes(x=peakTime, y=..density.., fill = Channel)) +
  geom_histogram(binwidth = 0.5, position=position_dodge())+
  ggtitle("long stim- neurons vs astrocytes-cyto data") + 
  geom_vline(xintercept=cyto.pt.median, linetype="dashed", size=1) +
  geom_vline(xintercept=rcamp1.pt.median, linetype="dashed", size=1) +
  xlab("Peak Time (s)") + 
  xlim(-0.5,15)+
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

####
#compare cyto to lck peak times

cyto<-subset(all.cyto.peaks.short, Channel=="GCaMP" & Condition!="Nostim")
lck<-subset(all.lck.peaks.short, Channel=="GCaMP"& Condition!="Nostim")

cyto$Channel="cyto"
lck$Channel="lck"
lck$Group=NULL
GCaMP.pT<-rbind(cyto,lck)

df.pT2C1<- summarySE(GCaMP.pT, measurevar = "peakTime", groupvars = c("Channel"))
df.pT2C2<- summarySE(GCaMP.pT, measurevar = "peakTime", groupvars = c("Channel", "Condition"))
df.pT2D<- summarySE(GCaMP.pT, measurevar = "peakTime", groupvars = c("Channel", "Condition","ROIType"))

Condition_Channel=interaction(GCaMP.pT$Condition,GCaMP.pT$Channel)
# stats for onset times- neurons vs astrocytes
pT.cytovslck.type.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP.pT,REML=FALSE)
pT.cytovslck.type.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP.pT,REML=FALSE)
pT.cytovslck.type.model2 = lmer(peakTime ~ Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP.pT,REML=FALSE)
pT.cytovslck.type.model4 = lmer(peakTime ~ Channel+Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP.pT,REML=FALSE)
pT.cytovslck.type.model5 = lmer(peakTime ~ Condition_Channel + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP.pT,REML=FALSE)
pT.cytovslck.type.anova <- anova(pT.cytovslck.type.null, pT.cytovslck.type.model1, pT.cytovslck.type.model2, 
                            pT.cytovslck.type.model4, pT.cytovslck.type.model5)
print(pT.cytovslck.type.anova)

pT.cytovslck.channel<- glht(pT.cytovslck.type.model1, mcp(Channel= "Tukey"))
summary(pT.cytovslck.channel)

pT.cytovslck.Cond_Ch<- glht(pT.cytovslck.type.model5, mcp(Condition_Channel= "Tukey"))
summary(pT.cytovslck.Cond_Ch)

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

#histogram bins
histseq= seq(0,stimwindow, 1)

# neuronal peak time histogram
# counts for each condition in the histogram
Nostim=hist(all.cyto.DSP4.peaks$peakTime[(all.cyto.DSP4.peaks$Channel=="RCaMP" & all.cyto.DSP4.peaks$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim=hist(all.cyto.DSP4.peaks$peakTime[(all.cyto.DSP4.peaks$Channel=="RCaMP" & all.cyto.DSP4.peaks$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts
Shortstim=hist(all.cyto.DSP4.peaks$peakTime[(all.cyto.DSP4.peaks$Channel=="RCaMP" & all.cyto.DSP4.peaks$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim=Nostim/ntrials.cyto.DSP4.peaks$ntrials[(ntrials.cyto.DSP4.peaks$Condition=="Nostim")]
Stim=Stim/ntrials.cyto.DSP4.peaks$ntrials[(ntrials.cyto.DSP4.peaks$Condition=="Stim")]
Shortstim=Shortstim/ntrials.cyto.DSP4.peaks$ntrials[(ntrials.cyto.DSP4.peaks$Condition=="shortstim")]

#make a data frame for plotting
Neuron.cyto.DSP4.histo <- data.frame(cbind(Nostim, Stim, Shortstim))
Neuron.cyto.DSP4.histo$time<-histseq[2:length(histseq)]
Neuron.cyto.DSP4.histo<-melt(Neuron.cyto.DSP4.histo,id="time")

ggplot(data=Neuron.cyto.DSP4.histo, aes(x=time, y= value, colour=variable)) + 
  geom_line(size=1)+
  xlab("Time to peak max (s)") +
  ylab("number of peaks per trial")+
  ggtitle("Neurons- duration- ROIs per trial- cyto.DSP4 data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

# astrocyte cyto dsp4 onset histogram
#histogram bins
histseq= seq(0,stimwindow, 1)

# counts for each condition in the histogram
Nostim=hist(all.cyto.DSP4.peaks$peakTime[(all.cyto.DSP4.peaks$Channel=="GCaMP" & all.cyto.DSP4.peaks$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim=hist(all.cyto.DSP4.peaks$peakTime[(all.cyto.DSP4.peaks$Channel=="GCaMP" & all.cyto.DSP4.peaks$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts
Shortstim=hist(all.cyto.DSP4.peaks$peakTime[(all.cyto.DSP4.peaks$Channel=="GCaMP" & all.cyto.DSP4.peaks$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim=Nostim/ntrials.cyto.DSP4.peaks$ntrials[(ntrials.cyto.DSP4.peaks$Condition=="Nostim")]
Stim=Stim/ntrials.cyto.DSP4.peaks$ntrials[(ntrials.cyto.DSP4.peaks$Condition=="Stim")]
Shortstim=Shortstim/ntrials.cyto.DSP4.peaks$ntrials[(ntrials.cyto.DSP4.peaks$Condition=="shortstim")]

#make a data frame for plotting
AC.cyto.DSP4.histo<- data.frame(cbind(Nostim, Stim, Shortstim))
AC.cyto.DSP4.histo$time<-histseq[2:length(histseq)]
AC.cyto.DSP4.histo<-melt(AC.cyto.DSP4.histo,id="time")

ggplot(data=AC.cyto.DSP4.histo, aes(x=time, y= value, colour=variable)) + 
  geom_line(size=1)+
  xlab("Time to peak max (s)") +
  ylab("number of peaks per trial")+
  ggtitle("Astrocytes- duration- ROIs per trial- cyto.DSP4 data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme


##########
# compare to other control trials
# COMPARE CONTROL VS DSP4 peak times
all.cyto.DSP4.peaks$treatment="DSP4"
all.cyto.peaks$treatment="control"

control.vs.DSP4.peaks<-rbind(all.cyto.DSP4.peaks,all.cyto.peaks)

ggplot(control.vs.DSP4.peaks[(control.vs.DSP4.peaks$Channel!="RCaMP" & control.vs.DSP4.peaks$Condition!="Nostim"),], aes(x=peakTime, colour = Condition, linetype=treatment)) + 
  stat_ecdf(lwd=1) +
  ggtitle("cdf- astrocytes-control vs DSP4- peak time") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

ggplot(control.vs.DSP4.peaks[(control.vs.DSP4.peaks$Channel!="GCaMP"),], aes(x=peakTime, colour = Condition, linetype=treatment)) + 
  stat_ecdf(lwd=1) +
  ggtitle("cdf- neurons-control vs DSP4- stim types") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

# mean onset times
df.peaks3<- summarySE(control.vs.DSP4.peaks, measurevar = "peakTime", groupvars = c("Channel", "treatment"))

df.peaks3$Channel <- factor(df.peaks3$Channel, levels = c("RCaMP","GCaMP"))
df.peaks.dsp4<- summarySE(all.cyto.DSP4.peaks, measurevar = "peakTime", groupvars = c("Channel", "Condition","Drug"))
df.peaks.dsp5<- summarySE(all.cyto.DSP4.peaks, measurevar = "peakTime", groupvars = c("Channel", "Condition","Drug","ROIType"))

df.peaks.dsp4$Condition <- factor(df.peaks.dsp4$Condition, levels = c("Nostim","shortstim","Stim"))
#df.OT3$treatment <- factor(df.OT3$treatment, levels = c("fast","delayed"))

#df.OT1 = df.OT1[!(df.OT1$Channel=="RCaMP"&df.OT1$Group=="delayed"),]

ggplot(df.peaks3, aes(x=Channel,y=peakTime, fill= interaction(Channel,treatment))) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean peak Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.peaks.dsp4[df.peaks.dsp4$Channel=="GCaMP",], aes(x=interaction(Channel,Condition),y=peakTime, fill= Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean peak Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.peaks.dsp5[df.peaks.dsp5$Channel=="GCaMP",], aes(x=interaction(ROIType,Condition),y=peakTime, fill= Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean peak Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.peaks.dsp4[df.peaks.dsp4$Channel=="RCaMP",], aes(x=interaction(Channel,Condition),y=peakTime, fill= Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean peak Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

# delayed peak time with dsp4 treatment?

control.vs.DSP4.peaks$treatment<-as.factor(control.vs.DSP4.peaks$treatment)
control.vs.DSP4.peaks$treatment_Condition=interaction(control.vs.DSP4.peaks$treatment,control.vs.DSP4.peaks$Condition)
# stats for peak times- control vs DSP4
pT.cytovsDSP.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials), control.vs.DSP4.peaks[control.vs.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
pT.cytovsDSP.model1 = lmer(peakTime ~ Condition + (1|Animal) + (1|Spot) + (1|trials), control.vs.DSP4.peaks[control.vs.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
pT.cytovsDSP.model2 = lmer(peakTime ~ treatment + (1|Animal) + (1|Spot) + (1|trials), control.vs.DSP4.peaks[control.vs.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
pT.cytovsDSP.model3 = lmer(peakTime ~ Condition+treatment + (1|Animal) + (1|Spot) + (1|trials), control.vs.DSP4.peaks[control.vs.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
pT.cytovsDSP.model4 = lmer(peakTime ~ treatment_Condition + (1|Animal) + (1|Spot) + (1|trials), control.vs.DSP4.peaks[control.vs.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
pT.cytovsDSP.anova <- anova(pT.cytovsDSP.null, pT.cytovsDSP.model1, pT.cytovsDSP.model2, pT.cytovsDSP.model3, pT.cytovsDSP.model4)
print(pT.cytovsDSP.anova)

pT.cytovsDSP.Condition<- glht(pT.cytovsDSP.model1, mcp(Condition= "Tukey"))
summary(pT.cytovsDSP.Condition)

pT.cytovsDSP.treatment_Condition<- glht(pT.cytovsDSP.model4, mcp(treatment_Condition= "Tukey"))
summary(pT.cytovsDSP.treatment_Condition)


all.cyto.DSP4.peaks$Drug<-as.factor(all.cyto.DSP4.peaks$Drug)
all.cyto.DSP4.peaks$SpotTime<-as.factor(all.cyto.DSP4.peaks$SpotTime)
all.cyto.DSP4.peaks$Drug_Condition=interaction(all.cyto.DSP4.peaks$Drug,all.cyto.DSP4.peaks$Condition)

# stats for peak times- control vs DSP4
pT.cytovsDSP2.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
pT.cytovsDSP2.model1 = lmer(peakTime ~ Condition + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
pT.cytovsDSP2.model2 = lmer(peakTime ~ Drug+ (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
pT.cytovsDSP2.model3 = lmer(peakTime ~ Condition+Drug + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
pT.cytovsDSP2.model4 = lmer(peakTime ~ Drug_Condition + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
pT.cytovsDSP2.anova <- anova(pT.cytovsDSP2.null, pT.cytovsDSP2.model1, pT.cytovsDSP2.model2, pT.cytovsDSP2.model3, 
                             pT.cytovsDSP2.model4)
print(pT.cytovsDSP2.anova)

pT.cytovsDSP.treatment_Condition2<- glht(pT.cytovsDSP2.model4, mcp(Drug_Condition= "Tukey"))
summary(pT.cytovsDSP.treatment_Condition2)


pT.cytovsDSP2.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="RCaMP",],REML=FALSE)
pT.cytovsDSP2.model1 = lmer(peakTime ~ Condition + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="RCaMP",],REML=FALSE)
pT.cytovsDSP2.model2 = lmer(peakTime ~ Drug+ (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="RCaMP",],REML=FALSE)
pT.cytovsDSP2.model3 = lmer(peakTime ~ Condition+Drug + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="RCaMP",],REML=FALSE)
pT.cytovsDSP2.model4 = lmer(peakTime ~ Drug_Condition + (1|Animal) + (1|Spot) + (1|SpotTime) + (1|trials), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="RCaMP",],REML=FALSE)
pT.cytovsDSP2.anova <- anova(pT.cytovsDSP2.null, pT.cytovsDSP2.model1, pT.cytovsDSP2.model2, pT.cytovsDSP2.model3, 
                             pT.cytovsDSP2.model4)
print(pT.cytovsDSP2.anova)

pT.cytovsDSP.treatment_Condition2<- glht(pT.cytovsDSP2.model4, mcp(Drug_Condition= "Tukey"))
summary(pT.cytovsDSP.treatment_Condition2)

#Results: not significantly different between control and DSP4 for astrocyte peak times!




######
# Astrocyte Duration histograms

# lck data

# astrocyte lck duration histogram
#histogram bins
histseq= seq(0,45, 1)

# counts for each condition in the histogram
Nostim=hist(all.lck.peaks$Duration[(all.lck.peaks$Channel=="GCaMP" & all.lck.peaks$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim=hist(all.lck.peaks$Duration[(all.lck.peaks$Channel=="GCaMP" & all.lck.peaks$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts
Shortstim=hist(all.lck.peaks$Duration[(all.lck.peaks$Channel=="GCaMP" & all.lck.peaks$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim=Nostim/ntrials.lck.peaks$ntrials[(ntrials.lck.peaks$Condition=="Nostim")]
Stim=Stim/ntrials.lck.peaks$ntrials[(ntrials.lck.peaks$Condition=="Stim")]
Shortstim=Shortstim/ntrials.lck.peaks$ntrials[(ntrials.lck.peaks$Condition=="shortstim")]

#make a data frame for plotting
AC.lck.histo<- data.frame(cbind(Nostim, Stim, Shortstim))
AC.lck.histo$time<-histseq[2:length(histseq)]
AC.lck.histo<-melt(AC.lck.histo,id="time")

ggplot(data=AC.lck.histo, aes(x=time, y= value, colour=variable)) + 
  geom_line(size=1)+
  xlab("Duration (s)") +
  ylab("number of peaks per trial")+
  ggtitle("Astrocytes- duration- ROIs per trial- lck data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

# density histogram
ggplot(data=all.lck.peaks[(all.lck.peaks$Channel=="GCaMP" & all.lck.peaks$Condition!="Nostim"),], aes(x=Duration, y=..density..))+ 
  geom_histogram(size=1, binwidth = 1)+
  xlab("Duration (s)") +
  ylab("density")+
  ggtitle("Astrocytes- duration after stim- density- lck data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

summary(all.lck.peaks$Duration[all.lck.peaks$Channel=="GCaMP" & all.lck.peaks$Condition!="Nostim"])

all.lck.peaks.group$Group<-as.factor(all.lck.peaks.group$Group)
all.lck.peaks.group$Group<-factor(all.lck.peaks.group$Group, levels=c("fast","delayed"))

# bar graph
df.dur.lck1<-summarySE(data=all.lck.peaks[all.lck.peaks$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Condition"))

df.dur.lck2<-summarySE(data=all.lck.peaks[all.lck.peaks$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Condition","ROIType"))

df.dur.lck3<-summarySE(data=all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Group"))
df.dur.lck4<-summarySE(data=all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Group","ROIType"))

ggplot(data=df.dur.lck1, aes(x=Condition, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
    max.theme

ggplot(data=df.dur.lck2, aes(x=ROIType, y=Duration, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("ROIType") +
  ylab("Duration (s)") +
  ggtitle("astrocyte lck duration- ROIType") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

ggplot(data=df.dur.lck3, aes(x=Group, y=Duration, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Group") +
  ylab("Duration (s)") +
  ggtitle("astrocyte lck duration- ROIType") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

all.lck.peaks$Condition_type<-interaction(all.lck.peaks$Condition,all.lck.peaks$ROIType)
# stats for astrocyte durations
dur.lck.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],REML=FALSE)
dur.lck.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],REML=FALSE)
dur.lck.model2= lmer(Duration ~ Condition_type + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],REML=FALSE)
dur.lck.model3= lmer(Duration ~ Group + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],REML=FALSE)
dur.lck.anova <- anova(dur.lck.null, dur.lck.model1,dur.lck.model2,dur.lck.model3)
print(dur.lck.anova)

dur.lck.Condition<- glht(dur.lck.model1, mcp(Condition= "Tukey"))
summary(dur.lck.Condition)

dur.lck.Condition_type<- glht(dur.lck.model2, mcp(Condition_type= "Tukey"))
summary(dur.lck.Condition_type)


#compare fast and delayed astrocytes
ggplot(all.lck.peaks.group[(all.lck.peaks.group$Channel!="GCaMP"),], aes(x=Duration, colour = Condition, linetype=Group)) + 
  stat_ecdf(lwd=1) +
  ggtitle("cdf- astrocytes-fast vs delyed-lck") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

ggplot(data=all.lck.peaks.group[(all.lck.peaks.group$Channel!="Nostim" & all.lck.peaks.group$Condition=="Stim"),], aes(x=Duration, y=..density.., fill=Group))+ 
  geom_histogram(size=1, binwidth = 1)+
  xlab("Duration (s)") +
  ylab("density")+
  ggtitle("Astrocytes- duration after stim- density- lck data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

df.dur.lck3<-summarySE(data=all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Condition","Group"))

ggplot(data=df.dur.lck3, aes(x=Condition, y=Duration, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("ROIType") +
  ylab("Duration (s)") +
  ggtitle("astrocyte lck duration- fast vs delayed") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

all.lck.peaks.group$Condition_group<-interaction(all.lck.peaks.group$Condition,all.lck.peaks.group$Group)
# stats for astrocyte durations
dur.lck.group.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",],REML=FALSE)
dur.lck.group.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",],REML=FALSE)
dur.lck.group.model2= lmer(Duration ~ Group + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",],REML=FALSE)
dur.lck.group.model3= lmer(Duration ~ Condition_group + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",],REML=FALSE)
dur.lck.group.anova <- anova(dur.lck.group.null, dur.lck.group.model1,dur.lck.group.model2,dur.lck.group.model3)
print(dur.lck.group.anova)

dur.lck.group.Condition<- glht(dur.lck.group.model1, mcp(Condition= "Tukey"))
summary(dur.lck.group.Condition)

dur.lck.group.Condition_group<- glht(dur.lck.group.model3, mcp(Condition_group= "Tukey"))
summary(dur.lck.group.Condition_group)

dur.lck.group.pv<- glht(dur.lck.group.model2, mcp(Group= "Tukey"))
summary(dur.lck.group.pv)

######
# astrocyte cyto duration

#histogram bins
histseq= seq(0,80, 1)

# counts for each condition in the histogram
Nostim=hist(all.cyto.peaks$Duration[(all.cyto.peaks$Channel=="GCaMP" & all.cyto.peaks$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim=hist(all.cyto.peaks$Duration[(all.cyto.peaks$Channel=="GCaMP" & all.cyto.peaks$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts
Shortstim=hist(all.cyto.peaks$Duration[(all.cyto.peaks$Channel=="GCaMP" & all.cyto.peaks$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim=Nostim/ntrials.cyto.peaks$ntrials[(ntrials.cyto.peaks$Condition=="Nostim")]
Stim=Stim/ntrials.cyto.peaks$ntrials[(ntrials.cyto.peaks$Condition=="Stim")]
Shortstim=Shortstim/ntrials.cyto.peaks$ntrials[(ntrials.cyto.peaks$Condition=="shortstim")]

#make a data frame for plotting
AC.cyto.histo<- data.frame(cbind(Nostim, Stim, Shortstim))
AC.cyto.histo$time<-histseq[2:length(histseq)]
AC.cyto.histo<-melt(AC.cyto.histo,id="time")

ggplot(data=AC.cyto.histo, aes(x=time, y= value, colour=variable)) + 
  geom_line(size=1)+
  xlab("Duration (s)") +
  ylab("number of peaks per trial")+
  ggtitle("Astrocytes- duration- ROIs per trial- cyto data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

# denstiy histogram
ggplot(data=all.cyto.peaks[(all.cyto.peaks$Channel=="GCaMP" & all.cyto.peaks$Condition!="Nostim"),], aes(x=Duration, y=..density..))+ 
  geom_histogram(size=1, binwidth = 2)+
  xlab("Duration (s)") +
  ylab("density")+
  ggtitle("Astrocytes- duration after stim- density- cyto data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

summary(all.cyto.peaks$Duration[all.cyto.peaks$Channel=="GCaMP" & all.cyto.peaks$Condition!="Nostim"])

# boxplot
all.cyto.peaks$Condition <- factor(all.cyto.peaks$Condition, levels = c("Nostim","shortstim","Stim"))

ggplot(data=all.cyto.peaks[all.cyto.peaks$Channel=="GCaMP",], aes(x=Condition, y= Duration, colour=Condition)) + 
  geom_boxplot(size=1)+
  ylab("Duration (s)")+
  ggtitle("Astrocytes- duration - cyto data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme


df.dur.cyto1<-summarySE(data=all.cyto.peaks[all.cyto.peaks$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Condition"))
df.dur.cyto2<-summarySE(data=all.cyto.peaks[all.cyto.peaks$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Condition","ROIType"))

ggplot(data=df.dur.cyto1, aes(x=Condition, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.dur.cyto2, aes(x=ROIType, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

all.cyto.peaks$Condition_type<-interaction(all.cyto.peaks$Condition,all.cyto.peaks$ROIType)
# stats for astrocyte durations
dur.cyto.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.cyto.peaks[all.cyto.peaks$Channel=="GCaMP",],REML=FALSE)
dur.cyto.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.cyto.peaks[all.cyto.peaks$Channel=="GCaMP",],REML=FALSE)
dur.cyto.model2= lmer(Duration ~ Condition_type + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.cyto.peaks[all.cyto.peaks$Channel=="GCaMP",],REML=FALSE)
dur.cyto.anova <- anova(dur.cyto.null, dur.cyto.model1,dur.cyto.model2)
print(dur.cyto.anova)

dur.cyto.Condition<- glht(dur.cyto.model1, mcp(Condition= "Tukey"))
summary(dur.cyto.Condition)

dur.cyto.Condition_type<- glht(dur.cyto.model2, mcp(Condition_type= "Tukey"))
summary(dur.cyto.Condition_type)

######
# astrocyte cyto DSP4 

# COMPARE CONTORL AND TREATMENT
# mean duration times
control.vs.DSP4.peaks$Condition<- factor(control.vs.DSP4.peaks$Condition, levels = c("Nostim","shortstim","Stim"))
all.cyto.DSP4.peaks$Condition<- factor(all.cyto.DSP4.peaks$Condition, levels = c("Nostim","shortstim","Stim"))


df.dur.cyto3<- summarySE(control.vs.DSP4.peaks, measurevar = "Duration", groupvars = c("Channel", "Condition","treatment"))
df.dur.cyto4<- summarySE(all.cyto.DSP4.peaks, measurevar = "Duration", groupvars = c("Channel", "Condition","Drug"))

ggplot(df.dur.cyto3[df.dur.cyto3$Channel!="RCaMP",], aes(x=Condition,y=Duration, fill= treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration (s)") +
  ggtitle("astrocytes") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.dur.cyto3[df.dur.cyto3$Channel!="GCaMP",], aes(x=Condition,y=Duration, fill= treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration (s)") +
  ggtitle("neurons") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

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

# stats for duration- control vs DSP4
dur.cytovsDSP.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), control.vs.DSP4.peaks[control.vs.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
dur.cytovsDSP.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials), control.vs.DSP4.peaks[control.vs.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
dur.cytovsDSP.model2 = lmer(Duration ~ treatment + (1|Animal) + (1|Spot) + (1|trials), control.vs.DSP4.peaks[control.vs.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
dur.cytovsDSP.model3 = lmer(Duration ~ Condition+treatment + (1|Animal) + (1|Spot) + (1|trials), control.vs.DSP4.peaks[control.vs.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
dur.cytovsDSP.model4 = lmer(Duration ~ treatment_Condition + (1|Animal) + (1|Spot) + (1|trials), control.vs.DSP4.peaks[control.vs.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
dur.cytovsDSP.anova <- anova(dur.cytovsDSP.null, dur.cytovsDSP.model1, dur.cytovsDSP.model2, dur.cytovsDSP.model3, dur.cytovsDSP.model4)
print(dur.cytovsDSP.anova)

dur.cytovsDSP.Condition<- glht(dur.cytovsDSP.model1, mcp(Condition= "Tukey"))
summary(dur.cytovsDSP.Condition)

dur.cytovsDSP.treatment_Condition<- glht(dur.cytovsDSP.model4, mcp(treatment_Condition= "Tukey"))
summary(dur.cytovsDSP.treatment_Condition)

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

# amplitude histograms

# astrocyte lck amplitude histogram
#histogram bins
histseq= seq(-1,15, 0.1)

# counts for each condition in the histogram
Nostim=hist(all.lck.peaks$amplitude[(all.lck.peaks$Channel=="GCaMP" & all.lck.peaks$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim=hist(all.lck.peaks$amplitude[(all.lck.peaks$Channel=="GCaMP" & all.lck.peaks$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts
Shortstim=hist(all.lck.peaks$amplitude[(all.lck.peaks$Channel=="GCaMP" & all.lck.peaks$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim=Nostim/ntrials.lck.peaks$ntrials[(ntrials.lck.peaks$Condition=="Nostim")]
Stim=Stim/ntrials.lck.peaks$ntrials[(ntrials.lck.peaks$Condition=="Stim")]
Shortstim=Shortstim/ntrials.lck.peaks$ntrials[(ntrials.lck.peaks$Condition=="shortstim")]

#make a data frame for plotting
AC.lck.histo<- data.frame(cbind(Nostim, Stim, Shortstim))
AC.lck.histo$time<-histseq[2:length(histseq)]
AC.lck.histo<-melt(AC.lck.histo,id="time")

ggplot(data=AC.lck.histo, aes(x=time, y= value, colour=variable)) + 
  geom_line(size=1)+
  xlab("amplitude dF/F") +
  ylab("number of peaks per trial")+
  ggtitle("Astrocytes- peak time- ROIs per trial- lck data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme


ggplot(all.lck.peaks[(all.lck.peaks$Channel=="GCaMP"&all.lck.peaks$Condition!="Nostim"),], aes(x=amplitude, y=..density.., colour = ROIType)) +
  geom_freqpoly(binwidth = 0.1, lwd=1)+
  ggtitle("stim- astrocytes-lck") + 
  xlab("amplitude") + 
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

ggplot(all.lck.peaks[(all.lck.peaks$Channel=="GCaMP"&all.lck.peaks$Condition!="Nostim"),], aes(x=amplitude, y=..density..)) +
  geom_histogram(binwidth = 0.3, lwd=1)+
  ggtitle("stim- astrocytes-lck") + 
  xlab("amplitude") + 
  scale_fill_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

summary(all.lck.peaks$amplitude[all.lck.peaks$Channel=="GCaMP" & all.lck.peaks$Condition!="Nostim"])

# boxplot
levels(all.lck.peaks$Condition) <- c("Nostim","shortstim","Stim")
ggplot(data=all.lck.peaks[all.lck.peaks$Channel=="GCaMP",], aes(x=Condition, y= amplitude, colour=Condition)) + 
  geom_boxplot(size=1)+
  ylab("amplitude dF/F")+
  ggtitle("Astrocytes- amplitude- - lck data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme


df.amp.lck1<-summarySE(data=all.lck.peaks[all.lck.peaks$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("Condition"))

df.amp.lck2<-summarySE(data=all.lck.peaks[all.lck.peaks$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("Condition","ROIType"))


ggplot(data=df.amp.lck1, aes(x=Condition, y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amp") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp.lck2, aes(x=ROIType, y=amplitude, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("ROIType") +
  ylab("amplitude") +
  ggtitle("astrocyte lck amplitude- ROIType") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

# stats for astrocyte durations
amp.lck.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],REML=FALSE)
amp.lck.model1 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],REML=FALSE)
amp.lck.model2= lmer(amplitude ~ Condition_type + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],REML=FALSE)
amp.lck.anova <- anova(amp.lck.null, amp.lck.model1,amp.lck.model2)
print(amp.lck.anova)

amp.lck.Condition<- glht(amp.lck.model1, mcp(Condition= "Tukey"))
summary(amp.lck.Condition)

amp.lck.Condition_type<- glht(amp.lck.model2, mcp(Condition_type= "Tukey"))
summary(amp.lck.Condition_type)


# amplitude of fast vs delayed

df.amp.lck3<-summarySE(data=all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("Condition","Group"))
df.amp.lck4<-summarySE(data=all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("Group"))

ggplot(data=df.amp.lck3, aes(x=Condition, y=amplitude, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("ROIType") +
  ylab("amplitude") +
  ggtitle("astrocyte lck duration- fast vs delayed") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

ggplot(data=df.amp.lck4, aes(x=Group, y=amplitude, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("ROIType") +
  ylab("amplitude") +
  ggtitle("astrocyte lck duration- fast vs delayed") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

all.lck.peaks.group$Condition_group<-interaction(all.lck.peaks.group$Condition,all.lck.peaks.group$Group)
# stats for astrocyte durations
amp.lck.group.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",],REML=FALSE)
amp.lck.group.model1 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",],REML=FALSE)
amp.lck.group.model2= lmer(amplitude ~ Group + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",],REML=FALSE)
amp.lck.group.model3= lmer(amplitude ~ Condition_group + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",],REML=FALSE)
amp.lck.group.anova <- anova(amp.lck.group.null, amp.lck.group.model1,amp.lck.group.model2,amp.lck.group.model3)
print(amp.lck.group.anova)

amp.lck.group.Condition<- glht(amp.lck.group.model1, mcp(Condition= "Tukey"))
summary(amp.lck.group.Condition)

amp.lck.group.Condition_group<- glht(amp.lck.group.model3, mcp(Condition_group= "Tukey"))
summary(amp.lck.group.Condition_group)

amp.lck.group.pv<- glht(amp.lck.group.model2, mcp(Group= "Tukey"))
summary(amp.lck.group.pv)



######
# cyto amplitude

# astrocyte cyto duration

all.cyto.peaks2<-subset(all.cyto.peaks, amplitude<30)
#histogram bins
histseq= seq(-5,31, 1)

# counts for each condition in the histogram
Nostim=hist(all.cyto.peaks2$amplitude[(all.cyto.peaks2$Channel=="GCaMP" & all.cyto.peaks2$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim=hist(all.cyto.peaks2$amplitude[(all.cyto.peaks2$Channel=="GCaMP" & all.cyto.peaks2$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts
Shortstim=hist(all.cyto.peaks2$amplitude[(all.cyto.peaks2$Channel=="GCaMP" & all.cyto.peaks2$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim=Nostim/ntrials.cyto.peaks$ntrials[(ntrials.cyto.peaks$Condition=="Nostim")]
Stim=Stim/ntrials.cyto.peaks$ntrials[(ntrials.cyto.peaks$Condition=="Stim")]
Shortstim=Shortstim/ntrials.cyto.peaks$ntrials[(ntrials.cyto.peaks$Condition=="shortstim")]

#make a data frame for plotting
AC.cyto.histo<- data.frame(cbind(Nostim, Stim, Shortstim))
AC.cyto.histo$time<-histseq[2:length(histseq)]
AC.cyto.histo<-melt(AC.cyto.histo,id="time")

ggplot(data=AC.cyto.histo, aes(x=time, y= value, colour=variable)) + 
  geom_line(size=1)+
  xlab("amplitude (s)") +
  ylab("number of peaks per trial")+
  ggtitle("Astrocytes- amplitude- ROIs per trial- cyto data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

# denstiy histogram
ggplot(data=all.cyto.peaks2[(all.cyto.peaks2$Channel=="GCaMP" & all.cyto.peaks2$Condition!="Nostim"),], aes(x=amplitude, y=..density..))+ 
  geom_histogram(size=1, binwidth = 0.5)+
  xlab("amplitude (s)") +
  ylab("density")+
  ggtitle("Astrocytes- amplitude after stim- density- cyto data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

summary(all.cyto.peaks2$amplitude[all.cyto.peaks2$Channel=="GCaMP" & all.cyto.peaks2$Condition!="Nostim"])

# boxplot
all.cyto.peaks2$Condition <- factor(all.cyto.peaks2$Condition, levels = c("Nostim","shortstim","Stim"))

ggplot(data=all.cyto.peaks2[all.cyto.peaks2$Channel=="GCaMP",], aes(x=Condition, y= amplitude, colour=Condition)) + 
  geom_boxplot(size=1)+
  ylab("amplitude (s)")+
  ggtitle("Astrocytes- amplitude - cyto data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme


df.amp2<-summarySE(data=all.cyto.peaks[all.cyto.peaks$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("Condition"))

ggplot(data=df.amp2, aes(x=Condition, y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

all.cyto.peaks$Condition_type<-interaction(all.cyto.peaks$Condition,all.cyto.peaks$ROIType)
# stats for astrocyte amplitudes
amp.cyto.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.cyto.peaks[all.cyto.peaks$Channel=="GCaMP",],REML=FALSE)
amp.cyto.model1 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.cyto.peaks[all.cyto.peaks$Channel=="GCaMP",],REML=FALSE)
amp.cyto.model2= lmer(amplitude ~ Condition_type + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.cyto.peaks[all.cyto.peaks$Channel=="GCaMP",],REML=FALSE)
amp.cyto.anova <- anova(amp.cyto.null, amp.cyto.model1,amp.cyto.model2)
print(amp.cyto.anova)

amp.cyto.Condition<- glht(amp.cyto.model1, mcp(Condition= "Tukey"))
summary(amp.cyto.Condition)

amp.cyto.Condition_type<- glht(amp.cyto.model2, mcp(Condition_type= "Tukey"))
summary(amp.cyto.Condition_type)


######
# astrocyte cyto DSP4 

# COMPARE CONTORL AND TREATMENT

# mean amp
control.vs.DSP4.peaks$Condition<- factor(control.vs.DSP4.peaks$Condition, levels = c("Nostim","shortstim","Stim"))
df.amp.cyto3<- summarySE(control.vs.DSP4.peaks, measurevar = "amplitude", groupvars = c("Channel", "Condition","treatment"))
df.amp.cyto4<- summarySE(control.vs.DSP4.peaks, measurevar = "amplitude", groupvars = c("Channel", "Condition","treatment","ROIType"))

# comparing to early time point controls
df.amp.cyto5<- summarySE(all.cyto.DSP4.peaks, measurevar = "amplitude", groupvars = c("Channel", "Condition","Drug"))
df.amp.cyto6<- summarySE(all.cyto.DSP4.peaks, measurevar = "amplitude", groupvars = c("Channel", "Condition","Drug","ROIType"))
df.amp.cyto7<- summarySE(all.cyto.DSP4.peaks, measurevar = "amplitude", groupvars = c("Channel", "Condition","Drug","SpotTime"))

ggplot(df.amp.cyto3[df.amp.cyto3$Channel!="GCaMP",], aes(x=Condition,y=amplitude, fill= treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("neurons") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp.cyto3[df.amp.cyto3$Channel!="RCaMP",], aes(x=Condition,y=amplitude, fill= treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("astrocytes") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

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



# stats for peak times- control vs DSP4
amp.cytovsDSP.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), control.vs.DSP4.peaks[control.vs.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
amp.cytovsDSP.model1 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials), control.vs.DSP4.peaks[control.vs.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
amp.cytovsDSP.model2 = lmer(amplitude ~ treatment + (1|Animal) + (1|Spot) + (1|trials), control.vs.DSP4.peaks[control.vs.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
amp.cytovsDSP.model3 = lmer(amplitude ~ Condition+treatment + (1|Animal) + (1|Spot) + (1|trials), control.vs.DSP4.peaks[control.vs.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
amp.cytovsDSP.model4 = lmer(amplitude ~ treatment_Condition + (1|Animal) + (1|Spot) + (1|trials), control.vs.DSP4.peaks[control.vs.DSP4.peaks$Channel!="RCaMP",],REML=FALSE)
amp.cytovsDSP.anova <- anova(amp.cytovsDSP.null, amp.cytovsDSP.model1, amp.cytovsDSP.model2, amp.cytovsDSP.model3, amp.cytovsDSP.model4)
print(amp.cytovsDSP.anova)


amp.cytovsDSP.treatment_Condition<- glht(amp.cytovsDSP.model4, mcp(treatment_Condition= "Tukey"))
summary(amp.cytovsDSP.treatment_Condition)


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
#  Lck ROI area

lck.ROIarea<-subset(all.lck.peaks, (ROIType=="Dendrite" |ROIType=="Process"))
lck.ROIarea.group<-subset(all.lck.peaks.group, ROIType=="Dendrite" |ROIType=="Process")

lck.ROIarea$Condition<- factor(lck.ROIarea$Condition, levels=c("Nostim","shortstim","Stim"))
lck.ROIarea$Group<- factor(lck.ROIarea$Condition, levels=c("fast","delayed"))

#boxplot
ggplot(data=lck.ROIarea, aes(x=Condition, y= area, colour=ROIType)) + 
  geom_boxplot(size=1)+
  ylab("area (sq m)")+
  ggtitle("Astrocytes- roi area- - lck data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

df.area1<- summarySE(lck.ROIarea, measurevar="area", groupvars=c("Condition","ROIType"))
df.area2<- summarySE(lck.ROIarea.group, measurevar="area", groupvars=c("Condition","ROIType", "Group"))
df.area3<- summarySE(lck.ROIarea.group[lck.ROIarea.group$Condition!="Nostim",], measurevar="area", groupvars=c("ROIType", "Group"))

ggplot(data=df.area1, aes(x=ROIType, y= area, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIarea") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.area2, aes(x=interaction(ROIType,Group), y= area, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIarea") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.area3, aes(x=ROIType, y= area, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIarea") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.area3[df.area3$ROIType!="Dendrite",], aes(x=Group, y= area, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIarea") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

lck.ROIarea$Condition_type<-interaction(lck.ROIarea$Condition,lck.ROIarea$ROIType)
# stats for astrocyte amplitudes
lck.area.null = lmer(area ~ (1|Animal) + (1|Spot) + (1|trials), lck.ROIarea,REML=FALSE)
lck.area.model1 = lmer(area ~ Condition + (1|Animal) + (1|Spot) + (1|trials),lck.ROIarea ,REML=FALSE)
lck.area.model2= lmer(area~ Condition_type + (1|Animal) + (1|Spot) + (1|trials), lck.ROIarea, REML=FALSE)
lck.area.anova <- anova(lck.area.null, lck.area.model1, lck.area.model2)
print(lck.area.anova)

lck.area.Condition_type<- glht(lck.area.model2, mcp(Condition_type= "Tukey"))
summary(lck.area.Condition_type)

lck.ROIarea.group$Condition_type<-interaction(lck.ROIarea.group$Condition,lck.ROIarea.group$ROIType)
lck.ROIarea.group$Condition_type_group<-interaction(lck.ROIarea.group$Condition,lck.ROIarea.group$ROIType, lck.ROIarea.group$Group)
lck.ROIarea.group$type_group<-interaction(lck.ROIarea.group$ROIType, lck.ROIarea.group$Group)

lck.area.group.null = lmer(area ~ (1|Animal) + (1|Spot) + (1|trials), lck.ROIarea.group,REML=FALSE)
lck.area.group.model1 = lmer(area ~ Condition + (1|Animal) + (1|Spot) + (1|trials),lck.ROIarea.group ,REML=FALSE)
lck.area.group.model2= lmer(area~ Condition_type + (1|Animal) + (1|Spot) + (1|trials), lck.ROIarea.group, REML=FALSE)
lck.area.group.model3= lmer(area~ type_group + (1|Animal) + (1|Spot) + (1|trials), lck.ROIarea.group, REML=FALSE)
lck.area.group.model4= lmer(area~ Condition_type_group + (1|Animal) + (1|Spot) + (1|trials), lck.ROIarea.group, REML=FALSE)
lck.area.group.anova <- anova(lck.area.group.null, lck.area.group.model1, 
                              lck.area.group.model2,lck.area.group.model3,
                              lck.area.group.model4)
print(lck.area.group.anova)

lck.area.group.type<- glht(lck.area.group.model3, mcp(type_group= "Tukey"))
summary(lck.area.group.type)

lck.area.group.Condition_type<- glht(lck.area.group.model4, mcp(Condition_type_group= "Tukey"))
summary(lck.area.group.Condition_type)


######
# cyto and DSP4 ROI area

cyto.vs.dsp4.ROIarea<-subset(control.vs.DSP4.peaks, (ROIType=="Dendrite" |ROIType=="Process"))
cyto.vs.dsp4.ROIarea$Condition<- factor(cyto.vs.dsp4.ROIarea$Condition, levels=c("Nostim","shortstim","Stim"))

#boxplot
ggplot(data=cyto.vs.dsp4.ROIarea, aes(x=Condition, y= area, colour=interaction(ROIType,treatment))) + 
  geom_boxplot(size=1)+
  ylab("area (sq m)")+
  ggtitle("Astrocytes- roi area- - cyto data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

df.area1<- summarySE(cyto.vs.dsp4.ROIarea, measurevar="area", groupvars=c("Condition","ROIType","treatment"))
df.area2<- summarySE(cyto.vs.dsp4.ROIarea, measurevar="area", groupvars=c("treatment","ROIType"))
df.area3<- summarySE(cyto.vs.dsp4.ROIarea[cyto.vs.dsp4.ROIarea$treatment=="control",], measurevar="area", groupvars=c("Condition","ROIType"))

ggplot(data=df.area1[df.area1$treatment=="control",], aes(x=ROIType, y= area, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIarea") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.area1, aes(x=interaction(ROIType, Condition), y= area, fill=treatment)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIarea") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.area2, aes(x=ROIType, y= area, fill=treatment)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIarea") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.area3, aes(x=ROIType, y= area, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIarea") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

cyto.vs.dsp4.ROIarea$Condition_type<-interaction(cyto.vs.dsp4.ROIarea$Condition,cyto.vs.dsp4.ROIarea$ROIType)
cyto.vs.dsp4.ROIarea$treatment_type<-interaction(cyto.vs.dsp4.ROIarea$treatment,cyto.vs.dsp4.ROIarea$ROIType)
cyto.vs.dsp4.ROIarea$Condition_type_treatment<-interaction(cyto.vs.dsp4.ROIarea$Condition,cyto.vs.dsp4.ROIarea$ROIType,cyto.vs.dsp4.ROIarea$treatment)

# stats for astrocyte amplitudes
cyto.vs.dsp4.ROIarea.null = lmer(area ~ (1|Animal) + (1|Spot) + (1|trials), cyto.vs.dsp4.ROIarea,REML=FALSE)
cyto.vs.dsp4.ROIarea.model1 = lmer(area ~ Condition + (1|Animal) + (1|Spot) + (1|trials), cyto.vs.dsp4.ROIarea,REML=FALSE)
cyto.vs.dsp4.ROIarea.model2= lmer(area~ Condition_type + (1|Animal) + (1|Spot) + (1|trials), cyto.vs.dsp4.ROIarea, REML=FALSE)
cyto.vs.dsp4.ROIarea.model3= lmer(area~ treatment_type + (1|Animal) + (1|Spot) + (1|trials), cyto.vs.dsp4.ROIarea, REML=FALSE)
cyto.vs.dsp4.ROIarea.model4= lmer(area~ Condition_type_treatment + (1|Animal) + (1|Spot) + (1|trials), cyto.vs.dsp4.ROIarea, REML=FALSE)
cyto.vs.dsp4.ROIarea.anova <- anova(cyto.vs.dsp4.ROIarea.null, cyto.vs.dsp4.ROIarea.model1, cyto.vs.dsp4.ROIarea.model2,
                                    cyto.vs.dsp4.ROIarea.model3,cyto.vs.dsp4.ROIarea.model4)
print(cyto.vs.dsp4.ROIarea.anova)

cyto.vs.dsp4.ROIarea.treatment_type<- glht(cyto.vs.dsp4.ROIarea.model3, mcp(treatment_type= "Tukey"))
summary(cyto.vs.dsp4.ROIarea.treatment_type)

cyto.vs.dsp4.ROIarea.Condition_type_treatment <- glht(cyto.vs.dsp4.ROIarea.model4, mcp(Condition_type_treatment = "Tukey"))
summary(cyto.vs.dsp4.ROIarea.Condition_type_treatment )

######

# peaks per trial
lck.trials<-ddply(all.lck.peaks, c("Animal","Spot","trials","Condition","Channel"), summarise, nPeaks=length(amplitude))
lck.trials.type<-ddply(all.lck.peaks, c("Animal","Spot","trials","Condition","Channel","ROIType"), summarise, nPeaks=length(amplitude))

df.lck.numPeaks1<-summarySE(lck.trials, measurevar="nPeaks", groupvars=c("Condition","Channel"))

ggplot(data=df.lck.numPeaks1, aes(x=Channel, y= nPeaks, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=nPeaks-se, ymax=nPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("nPeaks") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

lck.trials$Condition_channel<-interaction(lck.trials$Condition,lck.trials$Channel)
# stats for astrocyte amplitudes
lck.nPeak.null = lmer(nPeaks ~ (1|Animal) + (1|Spot), lck.trials,REML=FALSE)
lck.nPeak.model1 = lmer(nPeaks ~ Condition + (1|Animal) + (1|Spot),lck.trials,REML=FALSE)
lck.nPeak.model2= lmer(nPeaks~ Condition_channel + (1|Animal) + (1|Spot),lck.trials, REML=FALSE)
lck.nPeak.anova <- anova(lck.nPeak.null, lck.nPeak.model1, lck.nPeak.model2)
print(lck.nPeak.anova)

lck.nPeak.Condition_channel<- glht(lck.nPeak.model2, mcp(Condition_channel= "Tukey"))
summary(lck.nPeak.Condition_channel)

#ROITypes
df.lck.numPeaks2<-summarySE(lck.trials.type, measurevar="nPeaks", groupvars=c("Condition","Channel","ROIType"))

ggplot(data=df.lck.numPeaks2, aes(x=interaction(Channel,ROIType), y= nPeaks, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=nPeaks-se, ymax=nPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("nPeaks") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


######

# cyto and DSP4

# peaks per trial
control.vs.DSP4.trials<-ddply(control.vs.DSP4.peaks, c("Animal","Spot","trials","Condition","Channel","treatment"), summarise, nPeaks=length(amplitude))
control.vs.DSP4.trials.type<-ddply(control.vs.DSP4.peaks, c("Animal","Spot","trials","Condition","Channel","ROIType","treatment"), summarise, nPeaks=length(amplitude))

all.cyto.DSP4.peaks<-subset(all.cyto.DSP4.peaks, peakTime<15)

all.cyto.DSP4.trials<-ddply(all.cyto.DSP4.peaks, c("Animal","Spot","trials","SpotTime","Condition","Channel","Drug"), summarise, nPeaks=length(amplitude))
all.cyto.DSP4.trials.type<-ddply(all.cyto.DSP4.peaks, c("Animal","Spot","trials","SpotTime","Condition","Channel","ROIType","Drug"), summarise, nPeaks=length(amplitude))


df.cyto.numPeaks1<-summarySE(control.vs.DSP4.trials, measurevar="nPeaks", groupvars=c("Channel","treatment"))
df.cyto.numPeaks2<-summarySE(control.vs.DSP4.trials, measurevar="nPeaks", groupvars=c("Condition","Channel","treatment"))
df.cyto.numPeaks3<-summarySE(control.vs.DSP4.trials.type, measurevar="nPeaks", groupvars=c("ROIType","Condition","Channel","treatment"))

df.cyto.numPeaks4<-summarySE(all.cyto.DSP4.trials, measurevar="nPeaks", groupvars=c("Channel","Condition","Drug"))
df.cyto.numPeaks5<-summarySE(all.cyto.DSP4.trials.type, measurevar="nPeaks", groupvars=c("ROIType","Condition","Channel","Drug"))


ggplot(data=df.cyto.numPeaks1, aes(x=Channel, y= nPeaks, fill=treatment)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=nPeaks-se, ymax=nPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("nPeaks/trial") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.cyto.numPeaks2, aes(x=interaction(Channel,Condition), y= nPeaks, fill=treatment)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=nPeaks-se, ymax=nPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("nPeaks/trial") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.cyto.numPeaks3, aes(x=interaction(Channel,interaction(ROIType,Condition)), y= nPeaks, fill=treatment)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=nPeaks-se, ymax=nPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("nPeaks/trial") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

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

#######
control.vs.DSP4.trials.type$Condition_channel<-interaction(control.vs.DSP4.trials.type$Condition,control.vs.DSP4.trials.type$Channel)
control.vs.DSP4.trials.type$treatment_channel<-interaction(control.vs.DSP4.trials.type$treatment,control.vs.DSP4.trials.type$Channel)
control.vs.DSP4.trials.type$Condition_channel_treatment<-interaction(control.vs.DSP4.trials.type$Condition,control.vs.DSP4.trials.type$Channel,
                                                                     control.vs.DSP4.trials.type$treatment)
control.vs.DSP4.trials.type$Condition_channel_treatment_type<-interaction(control.vs.DSP4.trials.type$Condition,control.vs.DSP4.trials.type$Channel,
                                                                          control.vs.DSP4.trials.type$treatment,control.vs.DSP4.trials.type$ROIType)
# stats for astrocyte amplitudes
cyto.nPeak.null = lmer(nPeaks ~ (1|Animal) + (1|Spot), control.vs.DSP4.trials.type,REML=FALSE)
cyto.nPeak.model1 = lmer(nPeaks ~ Condition + (1|Animal) + (1|Spot), control.vs.DSP4.trials.type,REML=FALSE)
cyto.nPeak.model2A= lmer(nPeaks~ Condition_channel + (1|Animal) + (1|Spot), control.vs.DSP4.trials.type, REML=FALSE)
cyto.nPeak.model2B= lmer(nPeaks~ treatment_channel + (1|Animal) + (1|Spot), control.vs.DSP4.trials.type, REML=FALSE)
cyto.nPeak.model3= lmer(nPeaks~ Condition_channel_treatment + (1|Animal) + (1|Spot), control.vs.DSP4.trials.type, REML=FALSE)
cyto.nPeak.model4= lmer(nPeaks~ Condition_channel_treatment_type + (1|Animal) + (1|Spot), control.vs.DSP4.trials.type, REML=FALSE)
cyto.nPeak.anova <- anova(cyto.nPeak.null, cyto.nPeak.model1, cyto.nPeak.model2A,cyto.nPeak.model2B,
                          cyto.nPeak.model3,cyto.nPeak.model4)
print(cyto.nPeak.anova)

cyto.nPeak.Condition_channel<- glht(cyto.nPeak.model2A, mcp(Condition_channel= "Tukey"))
summary(cyto.nPeak.Condition_channel)

cyto.nPeak.treatment_channel<- glht(cyto.nPeak.model2B, mcp(treatment_channel= "Tukey"))
summary(cyto.nPeak.treatment_channel)

cyto.nPeak.Condition_channel_t<- glht(cyto.nPeak.model3, mcp(Condition_channel_treatment= "Tukey"))
summary(cyto.nPeak.Condition_channel_t)

cyto.nPeak.Condition_channel_t_t<- glht(cyto.nPeak.model4, mcp(Condition_channel_treatment_type= "Tukey"))
summary(cyto.nPeak.Condition_channel_t_t)

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


#####

GCaMP$Condition <- factor(GCaMP$Condition, levels = c("shortstim","Stim"))

df.fast.amp1<- summarySE(GCaMP, measurevar = "amplitude", groupvars = c("Group"))
df.fast.amp2<- summarySE(GCaMP, measurevar = "amplitude", groupvars = c("Group","Condition"))
df.fast.amp3<- summarySE(GCaMP, measurevar = "amplitude", groupvars = c("Group","ROIType"))
df.fast.amp4<- summarySE(GCaMP, measurevar = "amplitude", groupvars = c("Group","Condition","ROIType"))


ggplot(data=df.fast.amp1, aes(x=Group, y= amplitude, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.fast.amp2, aes(x=Condition, y= amplitude, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.fast.amp4, aes(x=interaction(ROIType,Condition), y= amplitude, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


ggplot(data=df.fast.amp3, aes(x=ROIType, y= amplitude, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


GCaMP$Group_type<-interaction(GCaMP$Group, GCaMP$ROIType)
GCaMP$Group_type_Cond<-interaction(GCaMP$Group, GCaMP$ROIType, GCaMP$Condition)

group.amp.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP,REML=FALSE)
group.amp.model1 = lmer(amplitude  ~ Group + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP,REML=FALSE)
group.amp.model2= lmer(amplitude ~ Condition_group + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP, REML=FALSE)
group.amp.model3= lmer(amplitude ~ Group_type + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP, REML=FALSE)
group.amp.model4= lmer(amplitude ~ Group_type_Cond + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP, REML=FALSE)
group.amp.anova <- anova(group.amp.null, group.amp.model1, 
                         group.amp.model2,group.amp.model3,
                         group.amp.model4)
print(group.amp.anova)

group.amp.pv<- glht(group.amp.model1, mcp(Group= "Tukey"))
summary(group.amp.pv)

group.amp.cond.pv<- glht(group.amp.model2, mcp(Condition_group= "Tukey"))
summary(group.amp.cond.pv)

group.amp.Group_type.pv<- glht(group.amp.model3, mcp(Group_type= "Tukey"))
summary(group.amp.Group_type.pv)



#####
#duration

df.fast.dur1<- summarySE(GCaMP, measurevar = "Duration", groupvars = c("Group"))
df.fast.dur2<- summarySE(GCaMP, measurevar = "Duration", groupvars = c("Group","Condition"))
df.fast.dur3<- summarySE(GCaMP, measurevar = "Duration", groupvars = c("Group","ROIType"))
df.fast.dur4<- summarySE(GCaMP, measurevar = "Duration", groupvars = c("Group","Condition","ROIType"))


ggplot(data=df.fast.dur1, aes(x=Group, y= Duration, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.fast.dur2, aes(x=Condition, y= Duration, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.fast.dur4, aes(x=interaction(ROIType,Condition), y= Duration, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


ggplot(data=df.fast.dur3, aes(x=ROIType, y= Duration, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


group.dur.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP,REML=FALSE)
group.dur.model1 = lmer(Duration  ~ Group + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP,REML=FALSE)
group.dur.model2= lmer(Duration ~ Condition_group + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP, REML=FALSE)
group.dur.model3= lmer(Duration ~ Group_type + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP, REML=FALSE)
group.dur.model4= lmer(Duration ~ Group_type_Cond + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP, REML=FALSE)
group.dur.anova <- anova(group.dur.null, group.dur.model1, 
                         group.dur.model2,group.dur.model3,
                         group.dur.model4)
print(group.dur.anova)

group.dur.pv<- glht(group.dur.model1, mcp(Group= "Tukey"))
summary(group.dur.pv)

group.dur.cond.pv<- glht(group.dur.model2, mcp(Condition_group= "Tukey"))
summary(group.dur.cond.pv)

group.dur.Group_type.pv<- glht(group.dur.model3, mcp(Group_type= "Tukey"))
summary(group.dur.Group_type.pv)
