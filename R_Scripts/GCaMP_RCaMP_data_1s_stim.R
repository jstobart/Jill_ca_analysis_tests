
library("lme4")
library("lmerTest")
#library("lattice")
library("plyr")
library("dplyr")
library("ggplot2")
#library("gplots")
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

long.lck <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_longstim_28_04_2017.csv", header=TRUE, sep = ",")
short.lck <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_shortstim_28_04_2017.csv", header=TRUE, sep = ",")
nostim.lck <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_nostim_28_04_2017.csv", header=TRUE, sep = ",")

long.cyto<- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cytoGC&RC_2D_longstim_28_04_2017.csv", header=TRUE, sep = ",")
short.cyto <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cytoGC&RC_2D_shortstim_28_04_2017.csv", header=TRUE, sep = ",")
nostim.cyto <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cytoGC&RC_2D_nostim_28_04_2017.csv", header=TRUE, sep = ",")

longstim.lck.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/longstim_onset&AUC.csv", header=TRUE, sep = ",")
shortstim.lck.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/shortstim_onset&AUC.csv", header=TRUE, sep = ",")
nostim.lck.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/nostim_onset&AUC.csv", header=TRUE, sep = ",")

longstim.cyto.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_longstim_onset&AUC.csv", header=TRUE, sep = ",")
shortstim.cyto.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_shortstim_onset&AUC.csv", header=TRUE, sep = ",")
nostim.cyto.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_nostim_firstonset&AUC.csv", header=TRUE, sep = ",")

##### 
#home files
long.lck <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_longstim_28_04_2017.csv", header=TRUE, sep = ",")
short.lck <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_shortstim_28_04_2017.csv", header=TRUE, sep = ",")
nostim.lck <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_nostim_28_04_2017.csv", header=TRUE, sep = ",")

long.cyto<- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cytoGC&RC_2D_longstim_28_04_2017.csv", header=TRUE, sep = ",")
short.cyto <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cytoGC&RC_2D_shortstim_28_04_2017.csv", header=TRUE, sep = ",")
nostim.cyto <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cytoGC&RC_2D_nostim_28_04_2017.csv", header=TRUE, sep = ",")

longstim.lck.OT <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/longstim_onset&AUC.csv", header=TRUE, sep = ",")
shortstim.lck.OT <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/shortstim_onset&AUC.csv", header=TRUE, sep = ",")
nostim.lck.OT <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/nostim_onset&AUC.csv", header=TRUE, sep = ",")

longstim.cyto.OT <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_longstim_onset&AUC.csv", header=TRUE, sep = ",")
shortstim.cyto.OT <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_shortstim_onset&AUC.csv", header=TRUE, sep = ",")
nostim.cyto.OT <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/cyto_nostim_firstonset&AUC.csv", header=TRUE, sep = ",")


##########

lsm.options(pbkrtest.limit = 100000)

all.lck.peaks<-rbind(nostim.lck,long.lck,short.lck)
all.lck.OT<-rbind(nostim.lck.OT,longstim.lck.OT,shortstim.lck.OT)
all.lck.OT$Spot_trial_Cond<-paste(all.lck.OT$Spot_trial, all.lck.OT$Condition, sep="_")
all.lck.OT$ROIs_Cond<-paste(all.lck.OT$ROIs_trial, all.lck.OT$Condition, sep="_")

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
all.lck.peaks$trials_Cond<-paste(all.lck.peaks$trials, all.lck.peaks$Condition, sep= "_")

all.lck.peaks$ROIs_Cond<-paste(all.lck.peaks$ROIs_trial, all.lck.peaks$Condition, sep= "_")

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

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.lck.OT<-distinct(all.lck.OT2, ROIs_Cond, .keep_all = TRUE)

# only the first entry will be used
all.lck.peaks3<-all.lck.peaks2[order(all.lck.peaks2$peakTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.lck.peaks<-distinct(all.lck.peaks3, ROIs_Cond,.keep_all = TRUE)

##### 
# cytoGCaMP6s data

all.cyto.peaks<-rbind(nostim.cyto,long.cyto,short.cyto)
all.cyto.OT<-rbind(nostim.cyto.OT,longstim.cyto.OT,shortstim.cyto.OT)
all.cyto.OT$Spot_trial_Cond<-paste(all.cyto.OT$Spot_trial, all.cyto.OT$Condition, sep="_")
all.cyto.OT$ROIs_Cond<-paste(all.cyto.OT$ROIs_trial, all.cyto.OT$Condition, sep="_")

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


########

## NOTE: define a stimulation window
# for distributions: histograms

stimwindow=15

short.lck.OT.dist<-subset(all.lck.OT,Condition!="Stim" & OnsetTime<stimwindow)
short.cyto.OT.dist<-subset(all.cyto.OT,Condition!="Stim" & OnsetTime<stimwindow)


#####
# for median and mean calculations

# no stim vs 8 s stim- neuronal window=9s, AC window= 15 s for peak time, neuronal window=2s, AC window=12 s for onset
# no stim vs 1 s stim- neuronal window=2s, AC window= 10 s for peak time, neuronal window=2s, AC window=8 s for onset


ShortN_PTwind=2
ShortAC_PTwind=10

ShortN_OTwind=2
ShortAC_OTwind=8

# remove data that is outside the above windows
# lck

short.lck.OT<-subset(all.lck.OT,Condition!="Stim")
short.lck.OT.R<-subset(short.lck.OT, Channel=="RCaMP" & OnsetTime<=ShortN_OTwind)
short.lck.OT.G<-subset(short.lck.OT, Channel=="GCaMP" & OnsetTime<=ShortAC_OTwind)

short.lck.OT.window<-rbind(short.lck.OT.R, short.lck.OT.G)


#cyto

short.cyto.OT<-subset(all.cyto.OT,Condition!="Stim")
short.cyto.OT.R<-subset(short.cyto.OT, Channel=="RCaMP" & OnsetTime<=ShortN_OTwind)
short.cyto.OT.G<-subset(short.cyto.OT, Channel=="GCaMP" & OnsetTime<=ShortAC_OTwind)

short.cyto.OT.window<-rbind(short.cyto.OT.R, short.cyto.OT.G)



# peak times
# lck

short.lck.PT<-subset(all.lck.peaks,Condition!="Stim")
short.lck.PT.R<-subset(short.lck.PT, Channel=="RCaMP" & peakTime<=ShortN_PTwind & peakTime>=0)
short.lck.PT.G<-subset(short.lck.PT, Channel=="GCaMP" & peakTime<=ShortAC_PTwind & peakTime>=0)

short.lck.peaks.window<-rbind(short.lck.PT.R, short.lck.PT.G)


#cyto
short.cyto.PT<-subset(all.cyto.peaks,Condition!="Stim")
short.cyto.PT.R<-subset(short.cyto.PT, Channel=="RCaMP" & peakTime<=ShortN_PTwind & peakTime>=0)
short.cyto.PT.G<-subset(short.cyto.PT, Channel=="GCaMP" & peakTime<=ShortAC_PTwind & peakTime>=0)

short.cyto.peaks.window<-rbind(short.cyto.PT.R, short.cyto.PT.G)



######
# LCK DATA
# Onset time histograms- normalized to the number of trials

# long stim vs no stim
ntrials.lck.OT.long.R.dis<- ddply(short.lck.OT.dist[short.lck.OT.dist$Channel=="RCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))
ntrials.lck.OT.long.G.dis<- ddply(short.lck.OT.dist[short.lck.OT.dist$Channel=="GCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))


#histogram bins
histseq= seq(0,15,1)
Nostim.N=0
Stim.N=0
Nostim.A=0
Stim.A=0
zeroRow<-data.frame(cbind(Nostim.N, Stim.N,Nostim.A, Stim.A))

# neuronal lck onset histogram
# counts for each condition in the histogram
Nostim.N=hist(short.lck.OT.dist$OnsetTime[(short.lck.OT.dist$Channel=="RCaMP" & short.lck.OT.dist$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim.N=hist(short.lck.OT.dist$OnsetTime[(short.lck.OT.dist$Channel=="RCaMP" & short.lck.OT.dist$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

Nostim.A=hist(short.lck.OT.dist$OnsetTime[(short.lck.OT.dist$Channel=="GCaMP" & short.lck.OT.dist$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim.A=hist(short.lck.OT.dist$OnsetTime[(short.lck.OT.dist$Channel=="GCaMP" & short.lck.OT.dist$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim.N=Nostim.N/ntrials.lck.OT.long.R.dis$ntrials[(ntrials.lck.OT.long.R.dis$Condition=="Nostim")]
Stim.N=Stim.N/ntrials.lck.OT.long.R.dis$ntrials[(ntrials.lck.OT.long.R.dis$Condition=="shortstim")]

Nostim.A=Nostim.A/ntrials.lck.OT.long.G.dis$ntrials[(ntrials.lck.OT.long.G.dis$Condition=="Nostim")]
Stim.A=Stim.A/ntrials.lck.OT.long.G.dis$ntrials[(ntrials.lck.OT.long.G.dis$Condition=="shortstim")]

#make a data frame for plotting
lck.long.histo <- data.frame(cbind(Nostim.N, Stim.N,Nostim.A, Stim.A))
lck.long.histo2<-rbind(zeroRow,lck.long.histo)
lck.long.histo2$time<-histseq

ggplot(NULL, aes(x=time))+
  geom_line(data=lck.long.histo2, aes(y=Nostim.N, color="Nostim.N")) +
  geom_line(data=lck.long.histo2, aes(y=Stim.N, color="Stim.N")) +
  geom_line(data=lck.long.histo2, aes(y=Nostim.A*10, color="Nostim.A")) +
  geom_line(data=lck.long.histo2, aes(y=Stim.A*10, color="Stim.A")) +
  scale_y_continuous(sec.axis = sec_axis(~./10, name = "Astrocyte peaks/trial")) + # secondary axis
  ggtitle("shortstim- neurons vs astrocytes-lck data") + 
  xlab("Onset Time (s)") + 
  ylab("Neuron peaks/trial") + 
  xlim(-0.5,15) +
  max.theme



######
# number of ROIs per trial per field of view?

short.lck.OT.window$Channel <- factor(short.lck.OT.window$Channel, levels = c("RCaMP","GCaMP"))

ROInum.lck.short<-ddply(short.lck.OT.window, c("Animal","Spot","Spot_trial","Condition","Channel"), summarise, nROIs=length(OnsetTime))

short.lck.ROInum.mean<-summarySE(ROInum.lck.short, measurevar = "nROIs", groupvars = c("Channel", "Condition"))

ggplot(short.lck.ROInum.mean, aes(x=Channel,y=nROIs, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=nROIs-se, ymax=nROIs+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view") +
  max.theme

Condition_Channel2= interaction(ROInum.lck.short$Condition,ROInum.lck.short$Channel)
nROI.lck.short.null = lmer(nROIs ~ (1|Animal) + (1|Spot), ROInum.lck.short,REML=FALSE)
nROI.lck.short.model1 = lmer(nROIs~ Channel + (1|Animal) + (1|Spot), ROInum.lck.short,REML=FALSE)
nROI.lck.short.model2 = lmer(nROIs ~ Condition + (1|Animal) + (1|Spot), ROInum.lck.short,REML=FALSE)
nROI.lck.short.model3 = lmer(nROIs ~ Condition_Channel2 + (1|Animal) + (1|Spot), ROInum.lck.short,REML=FALSE)
nROI.lck.short.anova <- anova(nROI.lck.short.null, nROI.lck.short.model1,nROI.lck.short.model2,nROI.lck.short.model3)
print(nROI.lck.short.anova)

nROI.lck.short.Cond_Channel<- glht(nROI.lck.short.model3, mcp(Condition_Channel2= "Tukey"))
summary(nROI.lck.short.Cond_Channel)

######
#peak times

stimwindow=15
histseq= seq(-5,15,1)

short.lck.PT.dist<-subset(all.lck.peaks,Condition!="Stim" & peakTime<stimwindow)

# long stim vs no stim
ntrials.lck.PT.long.R.dis<- ddply(short.lck.PT.dist[short.lck.PT.dist$Channel=="RCaMP",], c("Condition"), summarise, ntrials=length(unique(trials)))
ntrials.lck.PT.long.G.dis<- ddply(short.lck.PT.dist[short.lck.PT.dist$Channel=="GCaMP",], c("Condition"), summarise, ntrials=length(unique(trials)))

# neuronal lck onset histogram
# counts for each condition in the histogram
Nostim.N=hist(short.lck.PT.dist$peakTime[(short.lck.PT.dist$Channel=="RCaMP" & short.lck.PT.dist$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim.N=hist(short.lck.PT.dist$peakTime[(short.lck.PT.dist$Channel=="RCaMP" & short.lck.PT.dist$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

Nostim.A=hist(short.lck.PT.dist$peakTime[(short.lck.PT.dist$Channel=="GCaMP" & short.lck.PT.dist$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim.A=hist(short.lck.PT.dist$peakTime[(short.lck.PT.dist$Channel=="GCaMP" & short.lck.PT.dist$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim.N=Nostim.N/ntrials.lck.PT.long.R.dis$ntrials[(ntrials.lck.PT.long.R.dis$Condition=="Nostim")]
Stim.N=Stim.N/ntrials.lck.PT.long.R.dis$ntrials[(ntrials.lck.PT.long.R.dis$Condition=="shortstim")]

Nostim.A=Nostim.A/ntrials.lck.PT.long.G.dis$ntrials[(ntrials.lck.PT.long.G.dis$Condition=="Nostim")]
Stim.A=Stim.A/ntrials.lck.PT.long.G.dis$ntrials[(ntrials.lck.PT.long.G.dis$Condition=="shortstim")]


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
  ggtitle("shortstim- neurons vs astrocytes-lck data") + 
  xlab("Peak Max Time (s)") + 
  ylab("Neuron peaks/trial") + 
  xlim(-0.5,15) +
  max.theme


##################### 
histseq= seq(0,15,1)
#cyto data
ntrials.cyto.OT.long.R.dis<- ddply(short.cyto.OT.window[short.cyto.OT.window$Channel=="RCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))
ntrials.cyto.OT.long.G.dis<- ddply(short.cyto.OT.window[short.cyto.OT.window$Channel=="GCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))

# neuronal cyto onset histogram
# counts for each condition in the histogram
Nostim.N=hist(short.cyto.OT.dist$OnsetTime[(short.cyto.OT.dist$Channel=="RCaMP" & short.cyto.OT.dist$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim.N=hist(short.cyto.OT.dist$OnsetTime[(short.cyto.OT.dist$Channel=="RCaMP" & short.cyto.OT.dist$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

Nostim.A=hist(short.cyto.OT.dist$OnsetTime[(short.cyto.OT.dist$Channel=="GCaMP" & short.cyto.OT.dist$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim.A=hist(short.cyto.OT.dist$OnsetTime[(short.cyto.OT.dist$Channel=="GCaMP" & short.cyto.OT.dist$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim.N=Nostim.N/ntrials.cyto.OT.long.R.dis$ntrials[(ntrials.cyto.OT.long.R.dis$Condition=="Nostim")]
Stim.N=Stim.N/ntrials.cyto.OT.long.R.dis$ntrials[(ntrials.cyto.OT.long.R.dis$Condition=="shortstim")]

Nostim.A=Nostim.A/ntrials.cyto.OT.long.G.dis$ntrials[(ntrials.cyto.OT.long.G.dis$Condition=="Nostim")]
Stim.A=Stim.A/ntrials.cyto.OT.long.G.dis$ntrials[(ntrials.cyto.OT.long.G.dis$Condition=="shortstim")]

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
  ggtitle("shortstim- neurons vs astrocytes-cyto data") + 
  xlab("Onset Time (s)") + 
  ylab("Neuron peaks/trial") + 
  xlim(-0.5,15) +
  max.theme

 ##############
# number of ROIs per trial per field of view?

short.cyto.OT.window$Channel <- factor(short.cyto.OT.window$Channel, levels = c("RCaMP","GCaMP"))

ROInum.cyto.short<-ddply(short.cyto.OT.window, c("Animal","Spot","Spot_trial","Condition","Channel"), summarise, nROIs=length(OnsetTime))

short.cyto.ROInum.mean<-summarySE(ROInum.cyto.short, measurevar = "nROIs", groupvars = c("Channel", "Condition"))

ggplot(short.cyto.ROInum.mean, aes(x=Channel,y=nROIs, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=nROIs-se, ymax=nROIs+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view") +
  max.theme

Condition_Channel2= interaction(ROInum.cyto.short$Condition,ROInum.cyto.short$Channel)
nROI.cyto.short.null = lmer(nROIs ~ (1|Animal) + (1|Spot), ROInum.cyto.short,REML=FALSE)
nROI.cyto.short.model1 = lmer(nROIs~ Channel + (1|Animal) + (1|Spot), ROInum.cyto.short,REML=FALSE)
nROI.cyto.short.model2 = lmer(nROIs ~ Condition + (1|Animal) + (1|Spot), ROInum.cyto.short,REML=FALSE)
nROI.cyto.short.model3 = lmer(nROIs ~ Condition_Channel2 + (1|Animal) + (1|Spot), ROInum.cyto.short,REML=FALSE)
nROI.cyto.short.anova <- anova(nROI.cyto.short.null, nROI.cyto.short.model1,nROI.cyto.short.model2,nROI.cyto.short.model3)
print(nROI.cyto.short.anova)

nROI.cyto.short.Cond_Channel<- glht(nROI.cyto.short.model3, mcp(Condition_Channel2= "Tukey"))
summary(nROI.cyto.short.Cond_Channel)

#####

#peak times

stimwindow=15
histseq= seq(-5,15,1)

short.cyto.PT.dist<-subset(all.cyto.peaks,Condition!="Stim" & peakTime<stimwindow)

# long stim vs no stim
ntrials.cyto.PT.long.R.dis<- ddply(short.cyto.PT.dist[short.cyto.PT.dist$Channel=="RCaMP",], c("Condition"), summarise, ntrials=length(unique(trials)))
ntrials.cyto.PT.long.G.dis<- ddply(short.cyto.PT.dist[short.cyto.PT.dist$Channel=="GCaMP",], c("Condition"), summarise, ntrials=length(unique(trials)))

# neuronal cyto onset histogram
# counts for each condition in the histogram
Nostim.N=hist(short.cyto.PT.dist$peakTime[(short.cyto.PT.dist$Channel=="RCaMP" & short.cyto.PT.dist$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim.N=hist(short.cyto.PT.dist$peakTime[(short.cyto.PT.dist$Channel=="RCaMP" & short.cyto.PT.dist$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

Nostim.A=hist(short.cyto.PT.dist$peakTime[(short.cyto.PT.dist$Channel=="GCaMP" & short.cyto.PT.dist$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim.A=hist(short.cyto.PT.dist$peakTime[(short.cyto.PT.dist$Channel=="GCaMP" & short.cyto.PT.dist$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim.N=Nostim.N/ntrials.cyto.PT.long.R.dis$ntrials[(ntrials.cyto.PT.long.R.dis$Condition=="Nostim")]
Stim.N=Stim.N/ntrials.cyto.PT.long.R.dis$ntrials[(ntrials.cyto.PT.long.R.dis$Condition=="shortstim")]

Nostim.A=Nostim.A/ntrials.cyto.PT.long.G.dis$ntrials[(ntrials.cyto.PT.long.G.dis$Condition=="Nostim")]
Stim.A=Stim.A/ntrials.cyto.PT.long.G.dis$ntrials[(ntrials.cyto.PT.long.G.dis$Condition=="shortstim")]

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
  ggtitle("short stim- neurons vs astrocytes-cyto data") + 
  xlab("Peak Max Time (s)") + 
  ylab("Neuron peaks/trial") + 
  xlim(-0.5,15) +
  max.theme


#########
# combine lck and cyto data sets for onset times
short.lck.OT.window$sensor<- "lck"
short.cyto.OT.window$sensor<-"cyto"
short.lck.OT.window$Channel<-paste(short.lck.OT.window$sensor, short.lck.OT.window$Channel, sep="_")
short.cyto.OT.window$Channel<-paste(short.cyto.OT.window$sensor, short.cyto.OT.window$Channel, sep="_")

short.both.OT.window<-rbind(short.lck.OT.window,short.cyto.OT.window)


# combine lck and cyto data sets for peak times

short.lck.peaks.window$sensor<- "lck"
short.cyto.peaks.window$sensor<-"cyto"
short.lck.peaks.window$Channel<-paste(short.lck.peaks.window$sensor, short.lck.peaks.window$Channel, sep="_")
short.cyto.peaks.window$Channel<-paste(short.cyto.peaks.window$sensor, short.cyto.peaks.window$Channel, sep="_")

short.both.peaks.window<-rbind(short.lck.peaks.window,short.cyto.peaks.window)


short.both.OT.window$Channel <- factor(short.both.OT.window$Channel, levels = c("cyto_RCaMP","cyto_GCaMP","lck_RCaMP","lck_GCaMP"))
short.both.peaks.window$Channel <- factor(short.both.peaks.window$Channel, levels = c("cyto_RCaMP","cyto_GCaMP","lck_RCaMP","lck_GCaMP"))



# means, medians, modes
# novel function for mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


# onset times

short.both.OT.mean<-summarySE(short.both.OT.window, measurevar = "OnsetTime", groupvars = c("Channel", "Condition"))

short.both.OT.stats<-ddply(short.both.OT.window, c("Condition","Channel"), summarise, mean=mean(OnsetTime),
                          median=median(OnsetTime), mode=getmode(OnsetTime))

ggplot(short.both.OT.mean, aes(x=Channel,y=OnsetTime, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme



ggplot(short.both.OT.stats, aes(x=Channel,y=median, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  ylab("Median Onset Time (s)") +
  max.theme

ggplot(short.both.OT.stats, aes(x=Channel,y=mode, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  ylab("Mode Onset Time (s)") +
  max.theme


ggplot(short.both.OT.window, aes(x=Channel,y=OnsetTime, fill= Condition)) +
  geom_boxplot()+
  ylab("Onset Time (s)") +
  max.theme

ggplot(short.both.OT.window, aes(x=Channel,y=OnsetTime, fill= Condition)) +
  geom_boxplot(notch=TRUE)+
  ylab("Onset Time (s)") +
  ggtitle("notched")+
  max.theme

ggplot(short.both.OT.window[short.both.OT.window$Condition=="shortstim",], aes(x=Channel,y=OnsetTime, fill= Condition)) +
  geom_boxplot(notch=TRUE)+
  ylab("Onset Time (s)") +
  ggtitle("notched")+
  max.theme

# means stats
Condition_Channel=interaction(short.both.OT.window$Condition,short.both.OT.window$Channel)

OT.both.short.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial_Cond), short.both.OT.window,REML=FALSE)
OT.both.short.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Spot_trial_Cond), short.both.OT.window,REML=FALSE)
OT.both.short.model2 = lmer(OnsetTime ~ Condition + (1|Animal) + (1|Spot) + (1|Spot_trial_Cond) , short.both.OT.window,REML=FALSE)
OT.both.short.model3 = lmer(OnsetTime ~ Condition_Channel + (1|Animal) + (1|Spot) + (1|Spot_trial_Cond), short.both.OT.window,REML=FALSE)
OT.both.short.anova <- anova(OT.both.short.null, OT.both.short.model1,OT.both.short.model2,OT.both.short.model3)
print(OT.both.short.anova)

OT.both.short.Cond_Channel<- glht(OT.both.short.model3, mcp(Condition_Channel= "Tukey"))
summary(OT.both.short.Cond_Channel)


# compare onset time distributions and their variances

# Anderson Darling
library("kSamples")
#library("SuppDist")

#lck no stim vs stim
OT.lck.NSvsS.adtest<- ad.test(short.both.OT.window$OnsetTime[short.both.OT.window$Channel=="lck_GCaMP" & short.both.OT.window$Condition=="Nostim"],
                              short.both.OT.window$OnsetTime[short.both.OT.window$Channel=="lck_GCaMP" & short.both.OT.window$Condition=="shortstim"])
print(OT.lck.NSvsS.adtest)

#cyto no stim vs stim
OT.cyto.NSvsS.adtest<- ad.test(short.both.OT.window$OnsetTime[short.both.OT.window$Channel=="cyto_GCaMP" & short.both.OT.window$Condition=="Nostim"],
                               short.both.OT.window$OnsetTime[short.both.OT.window$Channel=="cyto_GCaMP" & short.both.OT.window$Condition=="shortstim"])
print(OT.cyto.NSvsS.adtest)

# lck stim vs cyto stim
OT.both.SvsS.adtest<- ad.test(short.both.OT.window$OnsetTime[short.both.OT.window$Channel=="lck_GCaMP" & short.both.OT.window$Condition=="shortstim"],
                              short.both.OT.window$OnsetTime[short.both.OT.window$Channel=="cyto_GCaMP" & short.both.OT.window$Condition=="shortstim"])
print(OT.both.SvsS.adtest)

# lck no stim vs cyto no stim 
OT.both.NSvsNS.adtest<- ad.test(short.both.OT.window$OnsetTime[short.both.OT.window$Channel=="lck_GCaMP" & short.both.OT.window$Condition=="Nostim"],
                                short.both.OT.window$OnsetTime[short.both.OT.window$Channel=="cyto_GCaMP" & short.both.OT.window$Condition=="Nostim"])
print(OT.both.NSvsNS.adtest)



# peak times

short.both.PT.mean<-summarySE(short.both.peaks.window, measurevar = "peakTime", groupvars = c("Channel", "Condition"))

short.both.PT.stats<-ddply(short.both.peaks.window, c("Condition","Channel"), summarise, mean=mean(peakTime),
                          median=median(peakTime), mode=getmode(peakTime))

ggplot(short.both.PT.mean, aes(x=Channel,y=peakTime, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean peak Time (s)") +
  max.theme



ggplot(short.both.PT.stats, aes(x=Channel,y=median, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  ylab("Median peak Time (s)") +
  max.theme

ggplot(short.both.PT.stats, aes(x=Channel,y=mode, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  ylab("Mode peak Time (s)") +
  max.theme


ggplot(short.both.peaks.window, aes(x=Channel,y=peakTime, fill= Condition)) +
  geom_boxplot()+
  ylab("peak Time (s)") +
  max.theme

ggplot(short.both.peaks.window, aes(x=Channel,y=peakTime, fill= Condition)) +
  geom_boxplot(notch=TRUE)+
  ylab("peak Time (s)") +
  ggtitle("notched")+
  max.theme

ggplot(short.both.peaks.window[short.both.peaks.window$Condition=="shortstim",], aes(x=Channel,y=peakTime, fill= Condition)) +
  geom_boxplot(notch=TRUE)+
  ylab("peak Time (s)") +
  ggtitle("notched")+
  max.theme

# means stats
Condition_Channel2=interaction(short.both.peaks.window$Condition,short.both.peaks.window$Channel)

PT.both.short.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials_Cond), short.both.peaks.window,REML=FALSE)
PT.both.short.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials_Cond), short.both.peaks.window,REML=FALSE)
PT.both.short.model2 = lmer(peakTime ~ Condition + (1|Animal) + (1|Spot) + (1|trials_Cond) , short.both.peaks.window,REML=FALSE)
PT.both.short.model3 = lmer(peakTime ~ Condition_Channel2 + (1|Animal) + (1|Spot) + (1|trials_Cond), short.both.peaks.window,REML=FALSE)
PT.both.short.anova <- anova(PT.both.short.null, PT.both.short.model1,PT.both.short.model2,PT.both.short.model3)
print(PT.both.short.anova)

PT.both.short.Cond_Channel<- glht(PT.both.short.model3, mcp(Condition_Channel2= "Tukey"))
summary(PT.both.short.Cond_Channel)



# compare onset time distributions and their variances

# Anderson Darling
#library("kSamples")
#library("SuppDist")

#lck no stim vs stim
PT.lck.NSvsS.adtest<- ad.test(short.both.peaks.window$peakTime[short.both.peaks.window$Channel=="lck_GCaMP" & short.both.peaks.window$Condition=="Nostim"],
                              short.both.peaks.window$peakTime[short.both.peaks.window$Channel=="lck_GCaMP" & short.both.peaks.window$Condition=="shortstim"])
print(PT.lck.NSvsS.adtest)

#cyto no stim vs stim
PT.cyto.NSvsS.adtest<- ad.test(short.both.peaks.window$peakTime[short.both.peaks.window$Channel=="cyto_GCaMP" & short.both.peaks.window$Condition=="Nostim"],
                               short.both.peaks.window$peakTime[short.both.peaks.window$Channel=="cyto_GCaMP" & short.both.peaks.window$Condition=="shortstim"])
print(PT.cyto.NSvsS.adtest)

# lck stim vs cyto stim
PT.both.SvsS.adtest<- ad.test(short.both.peaks.window$peakTime[short.both.peaks.window$Channel=="lck_GCaMP" & short.both.peaks.window$Condition=="shortstim"],
                              short.both.peaks.window$peakTime[short.both.peaks.window$Channel=="cyto_GCaMP" & short.both.peaks.window$Condition=="shortstim"])
print(PT.both.SvsS.adtest)

# lck no stim vs cyto no stim 
PT.both.NSvsNS.adtest<- ad.test(short.both.peaks.window$peakTime[short.both.peaks.window$Channel=="lck_GCaMP" & short.both.peaks.window$Condition=="Nostim"],
                                short.both.peaks.window$peakTime[short.both.peaks.window$Channel=="cyto_GCaMP" & short.both.peaks.window$Condition=="Nostim"])
print(PT.both.NSvsNS.adtest)





