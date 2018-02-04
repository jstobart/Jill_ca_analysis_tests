library("lme4")
library("lmerTest")
library("lattice")
library("plyr")
library("dplyr")
library("ggplot2")
library("lsmeans")
library("Rmisc")
#library("MASS")
library("multcomp")
library("reshape2")
library("tidyr")
#library("data.table")
library("Hmisc")
library("stringr")

#################
#IP3R2KO

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


##### 
#home files

lck.peaks1 <- read.table("D:/Data/GCaMP_RCaMP/Revision/Lck_GCaMP/FilesforR/Control_untreated/Peaks_1stCohort_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
lck.OT1<-read.table("D:/Data/GCaMP_RCaMP/Revision/Lck_GCaMP/FilesforR/Control_untreated/OnsetTimes_1stCohort_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
lck.peaks2 <- read.table("D:/Data/GCaMP_RCaMP/Revision/Lck_GCaMP/FilesforR/Control_untreated/Peaks_2ndCohort_Lck_nostim_vs_longstim_01_2018.csv", header=TRUE, sep = ",")
lck.OT2<-read.table("D:/Data/GCaMP_RCaMP/Revision/Lck_GCaMP/FilesforR/Control_untreated/OnsetTimes_2ndCohort_Lck_nostim_vs_longstim_01_2018.csv", header=TRUE, sep = ",")


##########

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

# remove the data frames that are combined
rm(lck.peaks1,lck.peaks2,lck.OT1, lck.OT2)

# only consider IP3R2KO
all.lck.peaks<-all.lck.peaks[grepl("IP",all.lck.peaks$Animal),]
all.lck.OT<-all.lck.OT[grepl("IP",all.lck.OT$Animal),]

#############
# ONSET Times
# unique ROI names
all.lck.OT$ROIs_trial<-paste(all.lck.OT$Animal, all.lck.OT$Spot, all.lck.OT$Trial, all.lck.OT$ROI, sep= "_")
all.lck.OT$Spot_trial<-paste(all.lck.OT$Animal, all.lck.OT$Spot, all.lck.OT$Trial, sep= "_")
all.lck.OT$Spot_trial_Cond<-paste(all.lck.OT$Spot_trial, all.lck.OT$Condition, sep="_")
all.lck.OT$ROIs_Cond<-paste(all.lck.OT$ROIs_trial, all.lck.OT$Condition, sep="_")

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


# REMOVE duplicate entries from onset time and peak time data frames 
# only the first entry will be used
all.lck.OT2=all.lck.OT[order(all.lck.OT$OnsetTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.lck.OT<-distinct(all.lck.OT2, ROIs_Cond, .keep_all = TRUE)

#adjust onset time for the data with 2 s before stimulation included:

rm(all.lck.OT2,all.lck.OTA, all.lck.OTB)

# Genotype
all.lck.peaks$Genotype="IP3R2_WT"
all.lck.peaks$Genotype[grepl("IPRG1",all.lck.peaks$Animal)]="IP3R2_KO"
all.lck.peaks$Genotype[grepl("IPRG4",all.lck.peaks$Animal)]="IP3R2_KO"
all.lck.peaks$Genotype[grepl("IPRG5",all.lck.peaks$Animal)]="IP3R2_KO"
all.lck.peaks$Genotype[grepl("IPRG7",all.lck.peaks$Animal)]="IP3R2_KO"

all.lck.OT$Genotype="IP3R2_WT"
all.lck.OT$Genotype[grepl("IPRG1",all.lck.OT$Animal)]="IP3R2_KO"
all.lck.OT$Genotype[grepl("IPRG4",all.lck.OT$Animal)]="IP3R2_KO"
all.lck.OT$Genotype[grepl("IPRG5",all.lck.OT$Animal)]="IP3R2_KO"
all.lck.OT$Genotype[grepl("IPRG7",all.lck.OT$Animal)]="IP3R2_KO"

all.lck.peaks$Genotype<-as.factor(all.lck.peaks$Genotype)
all.lck.OT$Genotype<-as.factor(all.lck.OT$Genotype)
all.lck.OT$Genotype<-factor(all.lck.OT$Genotype,levels=c("IP3R2_WT","IP3R2_KO"))
all.lck.peaks$Genotype<-factor(all.lck.peaks$Genotype,levels=c("IP3R2_WT","IP3R2_KO"))

all.lck.peaks$ROIType= "none"
all.lck.peaksA<- subset(all.lck.peaks, Channel=="GCaMP")
all.lck.peaksB<- subset(all.lck.peaks, Channel=="RCaMP")

# ROITypes
all.lck.peaksA$ROIType[grepl("r",all.lck.peaksA$roiName)]="Process"
all.lck.peaksA$ROIType[grepl("E",all.lck.peaksA$roiName)]="Endfoot"
all.lck.peaksB$ROIType[grepl("r",all.lck.peaksB$roiName)]="Dendrite"
all.lck.peaksB$ROIType[grepl("D",all.lck.peaksB$roiName)]="Dendrite"
all.lck.peaksB$ROIType[grepl("N",all.lck.peaksB$roiName)]="Neuron"

all.lck.peaks<-rbind(all.lck.peaksA, all.lck.peaksB)
rm(all.lck.peaksA, all.lck.peaksB)
all.lck.peaks$ROIType<- as.factor(all.lck.peaks$ROIType)

#unique ROI names
all.lck.peaks$ROIs_trial<-paste(all.lck.peaks$Animal, all.lck.peaks$Spot, all.lck.peaks$Trial,all.lck.peaks$roiName, sep= "_")
all.lck.peaks$trials<-paste(all.lck.peaks$Animal, all.lck.peaks$Spot, all.lck.peaks$Trial, sep= "_")
all.lck.peaks$trials_Cond<-paste(all.lck.peaks$trials, all.lck.peaks$Condition, sep= "_")

all.lck.peaks$ROIs_Cond<-paste(all.lck.peaks$ROIs_trial, all.lck.peaks$Condition, sep= "_")


# count number of trials per spot
Spot.lck.ntrials<-ddply(all.lck.peaks, c("Animal","Genotype","Spot","Condition"), summarise, nTrials=length(unique(Trial)))
Spot.lck.ntrials$Ani_Spot_Cond<-paste(Spot.lck.ntrials$Animal, Spot.lck.ntrials$Spot, Spot.lck.ntrials$Condition, sep="_")


# remove ROIs with no peaks
all.lck.peaks<-all.lck.peaks[!(all.lck.peaks$peakType=="NoPeak"),]


# remove matching astrocyte process and soma ROIs
Overlap= all.lck.peaks$overlap!=0
all.lck.peaks<-all.lck.peaks[!Overlap,]
#OverlapROIs<-unique(nostim$ROIs_trial[Overlap])


#add baseline time to peaks table
all.lck.OT$BL_time<-all.lck.OT$baseline/all.lck.OT$FrameRate
all.lck.peaks<-merge(all.lck.peaks, all.lck.OT[, c("ROIs_Cond", "BL_time","OnsetTime", "TraceAUC1", "TraceAUC10")], by="ROIs_Cond", all.x=TRUE)

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

all.lck.peaks$Genotype<-factor(all.lck.peaks$Genotype, levels=c("IP3R2_WT", "IP3R2_KO"))
rm(all.lck.peaks2,all.lck.peaks3, all.lck.OT)


#################
# neuronal responses to stimulation
NeuronalStim<-subset(all.lck.OT, Channel=="RCaMP" & Condition=="Stim" & OnsetTime<8)

# should have an onset time in 8 s stimulus
ggplot(NeuronalStim[NeuronalStim$OnsetTime<10,],aes(x=OnsetTime,y=..density..,fill=Genotype)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP onset times between 0 and 10 s from stim trials")+
  max.theme


Neuron95Onset.KO<-quantile(NeuronalStim$OnsetTime[(NeuronalStim$Genotype=="IP3R2_KO")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
#print(Neuron95Onset.C)
NeuronPT50.KO<-Neuron95Onset.KO[[11]]

Neuron95Onset.WT<-quantile(NeuronalStim$OnsetTime[(NeuronalStim$Genotype=="IP3R2_WT")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
#print(Neuron95Onset.C)
NeuronPT50.WT<-Neuron95Onset.WT[[11]]

print(NeuronPT50.KO)
print(NeuronPT50.WT)

##
# Astrocyte onset times
AstroStim<-subset(all.lck.OT, Channel=="GCaMP" & Condition=="Stim")

# should have an onset time in 8 s stimulus
ggplot(AstroStim[AstroStim$OnsetTime<40,],aes(x=OnsetTime,fill=Genotype)) +
  geom_histogram(binwidth=(0.084), position="dodge") +
  ggtitle("Lck-GCaMP onset times between 0 and 15 s from stim trials")+
  max.theme

Astro95Onset.WT<-quantile(AstroStim$OnsetTime[(AstroStim$OnsetTime<40 & AstroStim$Genotype=="IP3R2_WT")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
print(Astro95Onset.WT)

Astro95Onset.KO<-quantile(AstroStim$OnsetTime[(AstroStim$OnsetTime<40 & AstroStim$Genotype=="IP3R2_KO")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
print(Astro95Onset.KO)

WT.Astro.OT80<-Astro95Onset.WT[[17]]


# Astrocyte peak times
AstroPeaks<-subset(all.lck.peaks, Channel=="GCaMP" & Condition=="Stim")

# should have an onset time in 8 s stimulus
ggplot(AstroPeaks[AstroPeaks$peakTime<40,],aes(x=peakTime,fill=Genotype)) +
  geom_histogram(binwidth=(0.084), position="dodge") +
  ggtitle("Lck-GCaMP peakTime times between 0 and 40 s from stim trials")+
  max.theme

Astro95Peaks.WT<-quantile(AstroPeaks$peakTime[(AstroPeaks$peakTime<40 & AstroPeaks$Genotype=="IP3R2_WT")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
print(Astro95Peaks.WT)

Astro95Peaks.KO<-quantile(AstroPeaks$peakTime[(AstroPeaks$peakTime<40 & AstroPeaks$Genotype=="IP3R2_KO")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
print(Astro95Peaks.KO)

WT.Astro.PT80<-Astro95Peaks.WT[[17]]


############# 
# Onset time histograms- normalized to the number of trials

stimwindow=15

GC.lck.OT.dist<-subset(all.lck.OT,Condition=="Stim" & OnsetTime<stimwindow & Channel=="GCaMP")
GC.lck.OT.NS.dist<-subset(all.lck.OT,Condition=="Nostim" & OnsetTime<stimwindow & Channel=="GCaMP")

# KO vs WT stim
ntrials.GC.KOvsWT.stim<- ddply(GC.lck.OT.dist, c("Genotype"), summarise, ntrials=length(unique(Spot_trial)))

ntrials.GC.KOvsWT.nostim<- ddply(GC.lck.OT.NS.dist, c("Genotype"), summarise, ntrials=length(unique(Spot_trial)))

#histogram bins
histseq= seq(0,15,0.5)
KO.A.S=0
WT.A.S=0
KO.A.NS=0
WT.A.NS=0
zeroRow<-data.frame(cbind(KO.A.S, WT.A.S,KO.A.NS, WT.A.NS))

# neuronal lck onset histogram
# counts for each condition in the histogram
KO.A.S=hist(GC.lck.OT.dist$OnsetTime[GC.lck.OT.dist$Genotype=="IP3R2_KO"], breaks=histseq, plot=FALSE)$counts
WT.A.S=hist(GC.lck.OT.dist$OnsetTime[GC.lck.OT.dist$Genotype=="IP3R2_WT"], breaks=histseq, plot=FALSE)$counts
KO.A.NS=hist(GC.lck.OT.NS.dist$OnsetTime[GC.lck.OT.NS.dist$Genotype=="IP3R2_KO"], breaks=histseq, plot=FALSE)$counts
WT.A.NS=hist(GC.lck.OT.NS.dist$OnsetTime[GC.lck.OT.NS.dist$Genotype=="IP3R2_WT"], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
KO.A.S=KO.A.S/ntrials.GC.KOvsWT.stim$ntrials[(ntrials.GC.KOvsWT.stim$Genotype=="IP3R2_KO")]
WT.A.S=WT.A.S/ntrials.GC.KOvsWT.stim$ntrials[(ntrials.GC.KOvsWT.stim$Genotype=="IP3R2_WT")]
KO.A.NS=KO.A.NS/ntrials.GC.KOvsWT.stim$ntrials[(ntrials.GC.KOvsWT.stim$Genotype=="IP3R2_KO")]
WT.A.NS=WT.A.NS/ntrials.GC.KOvsWT.stim$ntrials[(ntrials.GC.KOvsWT.stim$Genotype=="IP3R2_WT")]

#make a data frame for plotting
lck.histo <- data.frame(cbind(KO.A.S, WT.A.S,KO.A.NS, WT.A.NS))
lck.histo2<-rbind(zeroRow,lck.histo)
lck.histo2$time<-histseq

ggplot(NULL, aes(x=time))+
  geom_line(data=lck.histo2, aes(y=KO.A.S, color="KO.A.S")) +
  geom_line(data=lck.histo2, aes(y=WT.A.S, color="WT.A.S")) +
  geom_line(data=lck.histo2, aes(y=KO.A.NS, color="KO.A.NS")) +
  geom_line(data=lck.histo2, aes(y=WT.A.NS, color="WT.A.NS")) +
  ggtitle("IP3R2 KO vs WT astrocytes-lck data") + 
  xlab("Onset Time (s)") + 
  ylab("Astrocytes peaks/trial") + 
  xlim(-0.5,15) +
  max.theme


ggplot(NULL, aes(x=time))+
  geom_line(data=lck.histo2, aes(y=KO.A.S, color="KO.A.S")) +
  geom_line(data=lck.histo2, aes(y=WT.A.S, color="WT.A.S")) +
  ggtitle("IP3R2 KO vs WT astrocytes-lck data") + 
  xlab("Onset Time (s)") + 
  ylab("Astrocytes peaks/trial") + 
  xlim(-0.5,15) +
  max.theme


#compare distributions (what is plotted in the figure)

# ks test- astrocytes
OT.lck.GC.KOvsWT.kstest<- ks.test(KO.A.S,WT.A.S)
print(OT.lck.GC.KOvsWT.kstest)

############
# Peak time histograms- normalized to the number of trials

stimwindow=25

GC.lck.PT.dist<-subset(all.lck.peaks,Condition=="Stim" & peakTime<stimwindow & Channel=="GCaMP")
GC.lck.PT.NS.dist<-subset(all.lck.peaks,Condition=="Nostim" & peakTime<stimwindow & Channel=="GCaMP")

# KO vs WT stim
ntrials.GC.KOvsWT.stim<- ddply(GC.lck.PT.dist, c("Genotype"), summarise, ntrials=length(unique(trials)))

ntrials.GC.KOvsWT.nostim<- ddply(GC.lck.PT.NS.dist, c("Genotype"), summarise, ntrials=length(unique(trials)))

#histogram bins
histseq= seq(0,25,0.5)
KO.A.S=0
WT.A.S=0
KO.A.NS=0
WT.A.NS=0
zeroRow<-data.frame(cbind(KO.A.S, WT.A.S,KO.A.NS, WT.A.NS))

# neuronal lck onset histogram
# counts for each condition in the histogram
KO.A.S=hist(GC.lck.PT.dist$peakTime[GC.lck.PT.dist$Genotype=="IP3R2_KO"], breaks=histseq, plot=FALSE)$counts
WT.A.S=hist(GC.lck.PT.dist$peakTime[GC.lck.PT.dist$Genotype=="IP3R2_WT"], breaks=histseq, plot=FALSE)$counts
KO.A.NS=hist(GC.lck.PT.NS.dist$peakTime[GC.lck.PT.NS.dist$Genotype=="IP3R2_KO"], breaks=histseq, plot=FALSE)$counts
WT.A.NS=hist(GC.lck.PT.NS.dist$peakTime[GC.lck.PT.NS.dist$Genotype=="IP3R2_WT"], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
KO.A.S=KO.A.S/ntrials.GC.KOvsWT.stim$ntrials[(ntrials.GC.KOvsWT.stim$Genotype=="IP3R2_KO")]
WT.A.S=WT.A.S/ntrials.GC.KOvsWT.stim$ntrials[(ntrials.GC.KOvsWT.stim$Genotype=="IP3R2_WT")]
KO.A.NS=KO.A.NS/ntrials.GC.KOvsWT.stim$ntrials[(ntrials.GC.KOvsWT.stim$Genotype=="IP3R2_KO")]
WT.A.NS=WT.A.NS/ntrials.GC.KOvsWT.stim$ntrials[(ntrials.GC.KOvsWT.stim$Genotype=="IP3R2_WT")]

#make a data frame for plotting
lck.histo <- data.frame(cbind(KO.A.S, WT.A.S,KO.A.NS, WT.A.NS))
lck.histo2<-rbind(zeroRow,lck.histo)
lck.histo2$time<-histseq

ggplot(NULL, aes(x=time))+
  geom_line(data=lck.histo2, aes(y=KO.A.S, color="KO.A.S")) +
  geom_line(data=lck.histo2, aes(y=WT.A.S, color="WT.A.S")) +
  geom_line(data=lck.histo2, aes(y=KO.A.NS, color="KO.A.NS")) +
  geom_line(data=lck.histo2, aes(y=WT.A.NS, color="WT.A.NS")) +
  ggtitle("IP3R2 KO vs WT astrocytes-lck data") + 
  xlab("Peak Time (s)") + 
  ylab("Astrocytes peaks/trial") + 
  xlim(-0.5,25) +
  max.theme


ggplot(NULL, aes(x=time))+
  geom_line(data=lck.histo2, aes(y=KO.A.S, color="KO.A.S")) +
  geom_line(data=lck.histo2, aes(y=WT.A.S, color="WT.A.S")) +
  ggtitle("IP3R2 KO vs WT astrocytes-lck data") + 
  xlab("Peak Time (s)") + 
  ylab("Astrocytes peaks/trial") + 
  xlim(-0.5,25) +
  max.theme


#compare distributions (what is plotted in the figure)

# ks test- astrocytes
OT.lck.GC.KOvsWT.kstest<- ks.test(KO.A.S,WT.A.S)
print(OT.lck.GC.KOvsWT.kstest)



######

LongN_OTwind=8
LongAC_OTwind=12.76
fastTh<-1.09

# remove data that is outside the above windows
# lck
stim.lck.OT.R<-subset(all.lck.peaks, Channel=="RCaMP" & OnsetTime<=LongN_OTwind & peakTime>=0 & Duration<45)
stim.lck.OT.G<-subset(all.lck.peaks, Channel=="GCaMP" & OnsetTime<=LongAC_OTwind & peakTime>=0 & Duration<45)

lck.peaks.window<-rbind(stim.lck.OT.R, stim.lck.OT.G)

rm(stim.lck.OT.G,stim.lck.OT.R)


# onset times
ggplot(lck.peaks.window[(lck.peaks.window$Channel=="RCaMP" & lck.peaks.window$Condition=="Stim"),], aes(x=interaction(Channel, Genotype),y=OnsetTime, fill= Condition)) +
  geom_boxplot()+
  ylab("Onset Time (s)") +
  ggtitle("time window 0-12.76s, rcamp, stim")+
  max.theme

ggplot(lck.peaks.window[(lck.peaks.window$Channel=="GCaMP" & lck.peaks.window$Condition=="Stim"),], aes(x=interaction(Channel, Genotype),y=OnsetTime, fill= Condition)) +
  geom_boxplot()+
  ylab("Onset Time (s)") +
  ggtitle("time window 0-12.76s, gcamp, stim")+
  max.theme


#peak times
ggplot(lck.peaks.window[(lck.peaks.window$Channel=="RCaMP" & lck.peaks.window$Condition=="Stim"),], aes(x=interaction(Channel, Genotype),y=peakTime, fill= Condition)) +
  geom_boxplot()+
  ylab("peak Time (s)") +
  ggtitle("time window 0-12.76s, rcamp, stim")+
  max.theme

ggplot(lck.peaks.window[(lck.peaks.window$Channel=="RCaMP" & lck.peaks.window$Condition=="Stim"),], aes(x=interaction(Channel, Genotype),y=peakTime, fill= Condition)) +
  geom_boxplot()+
  ylab("peak Time (s)") +
  ggtitle("time window 0-12.76s, rcamp, stim")+
  max.theme

ggplot(lck.peaks.window[(lck.peaks.window$Channel=="GCaMP" & lck.peaks.window$Condition=="Stim"),], aes(x=interaction(Channel, Genotype),y=peakTime, fill= Condition)) +
  geom_boxplot()+
  ylab("peak Time (s)") +
  ggtitle("time window 0-12.76s, gcamp, stim")+
  max.theme



######

# identify "FAST" astrocytes
lck.peaks.window$Group<-0
lck.peaks.window$Group[lck.peaks.window$OnsetTime<fastTh]<-"fast_MDs"
lck.peaks.window$Group[lck.peaks.window$OnsetTime>=fastTh]<-"delayed_MDs"

lck.peaks.window$Group <- factor(lck.peaks.window$Group, levels = c("fast_MDs","delayed_MDs"))
lck.peaks.window$Channel <- factor(lck.peaks.window$Channel, levels = c("RCaMP","GCaMP"))


# fast vs delayed
lck.peaks.window$Channel_Group<-interaction(lck.peaks.window$Channel, lck.peaks.window$Group)
lck.peaks.window$Channel_Group<-as.factor(lck.peaks.window$Channel_Group)

stim.lck.compdata<-lck.peaks.window[!(lck.peaks.window$Channel=="RCaMP"& lck.peaks.window$Group=="delayed_MDs"),]

# take out the effect of Condition
# we are only interested in stim case
stim.lck.compdata.STIM<-subset(stim.lck.compdata, Condition=="Stim")


######

# fraction of active pixels from all the astrocyte pixels
# with peaks near stimulus

# only consider active based ROIs (i.e. processes)

stim.lck.processes<-subset(lck.peaks.window, Channel=="GCaMP" & ROIType=="Process")

# any ROIs with a signal in the trial
AllROI.lck.proc.Pixels<-ddply(stim.lck.processes, c("Genotype","ROIs_Cond","trials","Animal","Spot","Condition","Channel","pixelsize", "nFluoPix"), summarise, 
                              nActivePix=sum(nActivePix))

# only ROIS with a signal near the stimulus
#StimROI.lck.proc.Pixels<-ddply(Lck.processes[Lck.processes$peakTime<8,], c("ROIs_Cond","trials","Animal","Spot","Condition","Channel","pixelsize", "nFluoPix"), summarise, 
 #                              nActivePix=sum(nActivePix))

# summarize per trial and spot
Spot.lck.proc.AllPixels<-ddply(AllROI.lck.proc.Pixels, c("Genotype","trials","Animal","Spot","Condition","Channel","pixelsize", "nFluoPix"), summarise, nROIs=length(unique(ROIs_Cond)), 
                               nActivePix=sum(nActivePix))
Spot.lck.proc.AllPixels$FracActive=Spot.lck.proc.AllPixels$nActivePix/Spot.lck.proc.AllPixels$nFluoPix
Spot.lck.proc.AllPixels$FracActivePerROI=Spot.lck.proc.AllPixels$FracActive/Spot.lck.proc.AllPixels$nROIs


df.FracActive.lck.all<- summarySE(Spot.lck.proc.AllPixels, measurevar = "FracActive", groupvars = c("Condition","Genotype"))
df.FracActive.lck.perROI<- summarySE(Spot.lck.proc.AllPixels, measurevar = "FracActivePerROI", groupvars = c("Condition","Genotype"))

ggplot(data=df.FracActive.lck.all, aes(x=Genotype, y= FracActive, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=FracActive-se, ymax=FracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.FracActive.lck.perROI, aes(x=Genotype, y= FracActivePerROI, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=FracActivePerROI-se, ymax=FracActivePerROI+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive Per ROI") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

#only astrocytes
Cond_Genotype= interaction(Spot.lck.proc.AllPixels$Condition,Spot.lck.proc.AllPixels$Genotype)
FracActive.lck.stim.null = lmer(FracActive ~ (1|Animal), Spot.lck.proc.AllPixels,REML=FALSE)
FracActive.lck.stim.model1 = lmer(FracActive~ Condition+ (1|Animal), Spot.lck.proc.AllPixels,REML=FALSE)
FracActive.lck.stim.model2 = lmer(FracActive ~ Genotype + (1|Animal), Spot.lck.proc.AllPixels,REML=FALSE)
FracActive.lck.stim.model3 = lmer(FracActive ~ Condition + Genotype + (1|Animal), Spot.lck.proc.AllPixels,REML=FALSE)
FracActive.lck.stim.model4 = lmer(FracActive ~ Cond_Genotype + (1|Animal), Spot.lck.proc.AllPixels,REML=FALSE)
FracActive.lck.stim.anova <- anova(FracActive.lck.stim.null, FracActive.lck.stim.model1,FracActive.lck.stim.model2,
                                 FracActive.lck.stim.model3,FracActive.lck.stim.model4)
print(FracActive.lck.stim.anova)

FracActive.lck.stim.Cond_Genotype<- glht(FracActive.lck.stim.model4, mcp(Cond_Genotype= "Tukey"))
summary(FracActive.lck.stim.Cond_Genotype)

# not different between stim and no stim for IP3KOs!

#####################
lck.peaks.window$Animal_Spot<- paste(lck.peaks.window$Animal, lck.peaks.window$Spot, sep="_")

# number of ROIs in each trial for each field of view (across the whole trial)
lck.8strial<-lck.peaks.window[lck.peaks.window$OnsetTime<8,]

lck.8strial$Channel <- factor(lck.8strial$Channel, levels = c("RCaMP","GCaMP"))

ROInum.8strial<-ddply(lck.8strial, c("Animal","Spot","Genotype","Condition","Channel","Animal_Spot"), summarise, nROIs=length(OnsetTime))
ROInum.8strial.group<-ddply(lck.8strial, c("Animal","Spot","Genotype","Condition","Channel","Animal_Spot","Group"), summarise, nROIs=length(OnsetTime))

# add in number of trials
ROInum.8strial$Ani_Spot_Cond_Genotype<-paste(ROInum.8strial$Animal_Spot, ROInum.8strial$Condition,ROInum.8strial$Genotype, sep="_")
ROInum.8strial.group$Ani_Spot_Cond_Genotype<-paste(ROInum.8strial.group$Animal_Spot, ROInum.8strial.group$Condition,ROInum.8strial.group$Genotype, sep="_")

Spot.lck.ntrials$Ani_Spot_Cond_Genotype<-paste(Spot.lck.ntrials$Animal, Spot.lck.ntrials$Spot,Spot.lck.ntrials$Condition,Spot.lck.ntrials$Genotype, sep="_")
ROInum.8strial<-merge(ROInum.8strial, Spot.lck.ntrials[, c("Ani_Spot_Cond_Genotype", "nTrials")], by="Ani_Spot_Cond_Genotype", all.x=TRUE)
ROInum.8strial$ROIsPerTrial<-ROInum.8strial$nROIs/ROInum.8strial$nTrials
ROInum.8strial.group<-merge(ROInum.8strial.group, Spot.lck.ntrials[, c("Ani_Spot_Cond_Genotype", "nTrials")], by="Ani_Spot_Cond_Genotype", all.x=TRUE)
ROInum.8strial.group$ROIsPerTrial<-ROInum.8strial.group$nROIs/ROInum.8strial.group$nTrials

# remove outliers
#outlier <- boxplot.stats(ROInum.8strial$ROIsPerTrial[ROInum.8strial$Channel=="GCaMP"])$out
#ROInum.8strial<-subset(ROInum.8strial, !(ROIsPerTrial %in% outlier))

# mean
df.ROInum.8strial<-summarySE(ROInum.8strial, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition","Genotype"))
df.ROInum.8strial.group<-summarySE(ROInum.8strial.group, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition","Genotype","Group"))


# paired line plots
df.ROInum.8strial$ROIsPerTrialMean<-df.ROInum.8strial$ROIsPerTrial
df.ROInum.8strial$Chan_Cond_Genotype<-interaction(df.ROInum.8strial$Channel, df.ROInum.8strial$Condition,df.ROInum.8strial$Genotype)
ROInum.8strial$Chan_Cond_Genotype<-interaction(ROInum.8strial$Channel, ROInum.8strial$Condition,ROInum.8strial$Genotype)

ROInum.8strial<-merge(ROInum.8strial, df.ROInum.8strial[, c("Chan_Cond_Genotype", "ROIsPerTrialMean","se")], by="Chan_Cond_Genotype", all.x=TRUE)

df.ROInum.8strial.group$ROIsPerTrialMean<-df.ROInum.8strial.group$ROIsPerTrial
df.ROInum.8strial.group$Chan_Cond_Genotype_Group<-interaction(df.ROInum.8strial.group$Channel, 
                                                          df.ROInum.8strial.group$Condition,
                                                          df.ROInum.8strial.group$Genotype,
                                                          df.ROInum.8strial.group$Group)
ROInum.8strial.group$Chan_Cond_Genotype_Group<-interaction(ROInum.8strial.group$Channel, 
                                                       ROInum.8strial.group$Condition,
                                                       ROInum.8strial.group$Genotype,
                                                       ROInum.8strial.group$Group)

ROInum.8strial.group<-merge(ROInum.8strial.group, df.ROInum.8strial.group[, c("Chan_Cond_Genotype_Group", "ROIsPerTrialMean","se")], by="Chan_Cond_Genotype_Group", all.x=TRUE)


# plots
# RCaMP
ggplot(df.ROInum.8strial[df.ROInum.8strial$Channel=="RCaMP",], aes(x=Genotype,y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  max.theme

ggplot(df.ROInum.8strial[(df.ROInum.8strial$Channel=="RCaMP" & df.ROInum.8strial$Condition=="Stim"),], aes(x=Genotype,y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  max.theme

ggplot(df.ROInum.8strial[df.ROInum.8strial$Channel=="RCaMP",],aes(x=Genotype,y=ROIsPerTrial, group=Condition)) +
  geom_point(aes(x=Genotype,y=ROIsPerTrial),stat="identity", size=3)+
  geom_line(aes(x=Genotype, y=ROIsPerTrial, colour=Condition))+
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("RCaMP ROIs")+
  max.theme

ggplot(ROInum.8strial[ROInum.8strial$Channel=="RCaMP",], aes(x=interaction(Genotype, Channel),y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("Lck- ROIs per FOV") +
  max.theme

# paired line plot
ggplot(ROInum.8strial[(ROInum.8strial$Channel=="RCaMP" & ROInum.8strial$Condition=="Stim"),], aes(x=Genotype, y = ROIsPerTrial)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(x=Genotype, y=ROIsPerTrial,group=Animal_Spot), colour="#b5b5b5")+
  #ylim(0, 50)+
  geom_point(aes(x=Genotype, y=ROIsPerTrialMean), size = 5, colour="#7b3294")+
  geom_line(aes(x=Genotype, y=ROIsPerTrialMean,group=Animal_Spot), size=1.5, colour="#7b3294")+
  geom_errorbar(aes(x=Genotype,ymin=ROIsPerTrialMean-se, ymax=ROIsPerTrialMean+se), colour="#7b3294", width=0.2,  size=1.5,position=position_dodge(.9)) +
  max.theme


# gcamp
ggplot(df.ROInum.8strial[df.ROInum.8strial$Channel=="GCaMP",], aes(x=Genotype,y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  max.theme

ggplot(df.ROInum.8strial[(df.ROInum.8strial$Channel=="GCaMP" & df.ROInum.8strial$Condition=="Stim"),], aes(x=Genotype,y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  max.theme

ggplot(df.ROInum.8strial[df.ROInum.8strial$Channel=="GCaMP",],aes(x=Genotype,y=ROIsPerTrial, group=Condition)) +
  geom_point(aes(x=Genotype,y=ROIsPerTrial),stat="identity", size=3)+
  geom_line(aes(x=Genotype, y=ROIsPerTrial, colour=Condition))+
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("GCaMP ROIs") +
  max.theme

ggplot(ROInum.8strial[ROInum.8strial$Channel=="GCaMP",], aes(x=interaction(Genotype, Channel),y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("Lck- ROIs per FOV") +
  max.theme

# paired line plot
ggplot(ROInum.8strial[(ROInum.8strial$Channel=="GCaMP" & ROInum.8strial$Condition=="Stim"),], aes(x=Genotype, y = ROIsPerTrial)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(x=Genotype, y=ROIsPerTrial,group=Animal_Spot), colour="#b5b5b5")+
  #ylim(0, 50)+
  geom_point(aes(x=Genotype, y=ROIsPerTrialMean), size = 5, colour="#008837")+
  geom_line(aes(x=Genotype, y=ROIsPerTrialMean,group=Animal_Spot), size=1.5, colour="#008837")+
  geom_errorbar(aes(x=Genotype,ymin=ROIsPerTrialMean-se, ymax=ROIsPerTrialMean+se), colour="#008837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  max.theme

# only fast GcaMP

fastROInum<- subset(ROInum.8strial.group, Channel=="GCaMP" & Condition=="Stim" & Group=="fast_MDs")

ggplot(fastROInum, aes(x=interaction(Genotype, Channel),y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("Lck- fast ROIs per FOV") +
  max.theme

ggplot(fastROInum,aes(x=Genotype, y = ROIsPerTrial)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(x=Genotype, y=ROIsPerTrial,group=Animal_Spot), colour="#b5b5b5")+
  #ylim(0, 50)+
  geom_point(aes(x=Genotype, y=ROIsPerTrialMean), size = 5, colour="#008837")+
  geom_line(aes(x=Genotype, y=ROIsPerTrialMean, group=Animal_Spot), size=1.5, colour="#008837")+
  geom_errorbar(aes(x=Genotype,ymin=ROIsPerTrialMean-se, ymax=ROIsPerTrialMean+se), colour="#008837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  max.theme


ggplot(df.ROInum.8strial.group[df.ROInum.8strial.group$Channel=="GCaMP" & 
                                 df.ROInum.8strial.group$Group=="fast_MDs" &
                                 df.ROInum.8strial.group$Condition=="Stim",],
       aes(x=Genotype,y=ROIsPerTrial, colour=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("fast GCaMP ROIs") +
  max.theme



#only delayed GCaMP

delayedROInum<- subset(ROInum.8strial.group, Channel=="GCaMP" & Condition=="Stim" & Group=="delayed_MDs")

ggplot(delayedROInum, aes(x=interaction(Genotype, Channel),y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("Lck- delayed ROIs per FOV") +
  max.theme

ggplot(delayedROInum,aes(x=Genotype, y = ROIsPerTrial)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(x=Genotype, y=ROIsPerTrial,group=Animal_Spot), colour="#b5b5b5")+
  #ylim(0, 50)+
  geom_point(aes(x=Genotype, y=ROIsPerTrialMean), size = 5, colour="#008837")+
  geom_line(aes(x=Genotype, y=ROIsPerTrialMean, group=Animal_Spot), size=1.5, colour="#008837")+
  geom_errorbar(aes(x=Genotype,ymin=ROIsPerTrialMean-se, ymax=ROIsPerTrialMean+se), colour="#008837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  max.theme

ggplot(df.ROInum.8strial.group[df.ROInum.8strial.group$Channel=="GCaMP" & 
                                 df.ROInum.8strial.group$Group=="delayed_MDs" &
                                 df.ROInum.8strial.group$Condition=="Stim",],
       aes(x=Genotype,y=ROIsPerTrial, colour=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("delayed GCaMP ROIs") +
  max.theme



#stats
# RCaMP
Condition_Genotype_RC= interaction(ROInum.8strial$Condition[ROInum.8strial$Channel=="RCaMP"],ROInum.8strial$Genotype[ROInum.8strial$Channel=="RCaMP"])
nROI.RC.stim.null = lmer(ROIsPerTrial ~ (1|Animal) + (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="RCaMP",],REML=FALSE)
nROI.RC.stim.model1 = lmer(ROIsPerTrial ~ Condition + (1|Animal)+ (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="RCaMP",],REML=FALSE)
nROI.RC.stim.model2 = lmer(ROIsPerTrial ~ Condition_Genotype_RC + (1|Animal)+ (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="RCaMP",],REML=FALSE)
nROI.RC.stim.anova <- anova(nROI.RC.stim.null, nROI.RC.stim.model1,nROI.RC.stim.model2)
print(nROI.RC.stim.anova)

nROI.RC.stim.Cond_Genotype<- glht(nROI.RC.stim.model2, mcp(Condition_Genotype_RC= "Tukey"))
summary(nROI.RC.stim.Cond_Genotype)

# ALL GenotypeS ARE P<0.01 less than control for stim case
# atropine and metergoline do not have sig difference between no stim and stim

# GCaMP
Condition_Genotype_GC= interaction(ROInum.8strial$Condition[ROInum.8strial$Channel=="GCaMP"],ROInum.8strial$Genotype[ROInum.8strial$Channel=="GCaMP"])
Condition_Genotype_Group_GC= interaction(ROInum.8strial$Condition[ROInum.8strial$Channel=="GCaMP"],
                                         ROInum.8strial$Genotype[ROInum.8strial$Channel=="GCaMP"],
                                         ROInum.8strial$Group[ROInum.8strial$Channel=="GCaMP"])

nROI.GC.stim.null = lmer(ROIsPerTrial ~ (1|Animal) + (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="GCaMP",],REML=FALSE)
nROI.GC.stim.model1 = lmer(ROIsPerTrial ~ Condition + (1|Animal)+ (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="GCaMP",],REML=FALSE)
nROI.GC.stim.model2 = lmer(ROIsPerTrial ~ Condition_Genotype_GC + (1|Animal)+ (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="GCaMP",],REML=FALSE)
nROI.GC.stim.model3 = lmer(ROIsPerTrial ~ Condition_Genotype_Group_GC + (1|Animal)+ (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="GCaMP",],REML=FALSE)
nROI.GC.stim.anova <- anova(nROI.GC.stim.null, nROI.GC.stim.model1,nROI.GC.stim.model2,nROI.GC.stim.model3)
print(nROI.GC.stim.anova)

nROI.GC.stim.Cond_Genotype<- glht(nROI.GC.stim.model2, mcp(Condition_Genotype_GC= "Tukey"))
summary(nROI.GC.stim.Cond_Genotype)

# metergoline, trazodone, prazosin- no sig difference between stim and control

# trazodone, and prazosin- sig fewer ROIs compared to control


# fast GC
nROI.GC.fast.stim.null = lmer(ROIsPerTrial ~ (1|Animal) + (1|Spot), fastROInum, REML=FALSE)
nROI.GC.fast.stim.model2 = lmer(ROIsPerTrial ~ Genotype + (1|Animal) + (1|Spot),fastROInum, REML=FALSE)
nROI.GC.fast.stim.anova <- anova(nROI.GC.fast.stim.null, nROI.GC.fast.stim.model2)
print(nROI.GC.fast.stim.anova)

nROI.GC.fast.stim.Cond_Genotype<- glht(nROI.GC.fast.stim.model2, mcp(Genotype= "Tukey"))
summary(nROI.GC.fast.stim.Cond_Genotype)

# no difference in ROI number across fast AC types


# delayed GC
nROI.GC.delayed.stim.null = lmer(ROIsPerTrial ~ (1|Animal) + (1|Spot), delayedROInum, REML=FALSE)
nROI.GC.delayed.stim.model2 = lmer(ROIsPerTrial ~ Genotype + (1|Animal) + (1|Spot),delayedROInum, REML=FALSE)
nROI.GC.delayed.stim.anova <- anova(nROI.GC.delayed.stim.null, nROI.GC.delayed.stim.model2)
print(nROI.GC.delayed.stim.anova)

nROI.GC.delayed.stim.Cond_Genotype<- glht(nROI.GC.delayed.stim.model2, mcp(Genotype= "Tukey"))
summary(nROI.GC.delayed.stim.Cond_Genotype)




#########
# mean onset times
df.OT1<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim",], measurevar = "OnsetTime", groupvars = c("Channel_Group", "Genotype"))
df.OT2<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",], measurevar = "OnsetTime", groupvars = c("Group", "Genotype"))
df.OT3<- summarySE(stim.lck.compdata.STIM, measurevar = "OnsetTime", groupvars = c("ROIType","Channel_Group", "Genotype"))

ggplot(df.OT1, aes(x=interaction(Genotype,Channel_Group),y=OnsetTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme

ggplot(df.OT2, aes(x=Group,y=OnsetTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme

ggplot(df.OT3, aes(x=interaction(Channel_Group,ROIType),y=OnsetTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme



# stats
Group_Channel_Type_Cond_Gen=interaction(lck.peaks.window$Group,lck.peaks.window$Channel,lck.peaks.window$ROIType, 
                                        lck.peaks.window$Genotype,lck.peaks.window$Condition)
Group_Channel_Cond_Gen=interaction(lck.peaks.window$Group,lck.peaks.window$Channel,lck.peaks.window$Genotype,
                                   lck.peaks.window$Condition)

# stats for onset times- neurons vs astrocytes
OT.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|trials), lck.peaks.window,REML=FALSE)
OT.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), lck.peaks.window,REML=FALSE)
OT.model2 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), lck.peaks.window,REML=FALSE)
OT.model3 = lmer(OnsetTime ~ Group_Channel_Cond_Gen + (1|Animal) + (1|Spot) + (1|trials), lck.peaks.window,REML=FALSE)
OT.model4 = lmer(OnsetTime ~ Group_Channel_Type_Cond_Gen + (1|Animal) + (1|Spot) + (1|trials), lck.peaks.window,REML=FALSE)
OT.anova <- anova(OT.null, OT.model1,OT.model2,OT.model3,OT.model4)
print(OT.anova)

OT.Group_channel_Gen<- glht(OT.model3, mcp(Group_Channel_Cond_Gen= "Tukey"))
summary(OT.Group_channel_Gen)

#OT.Group_Channel_Type_Gen<- glht(OT.model4, mcp(Group_Channel_Type_Cond_Gen= "Tukey"))
#summary(OT.Group_Channel_Type_Gen)


astro.compdata.stim<-subset(stim.lck.compdata.STIM, Channel=="GCaMP")

Group_Type_Gen=interaction(astro.compdata.stim$Group,astro.compdata.stim$ROIType, astro.compdata.stim$Genotype)
Group_Gen=interaction(astro.compdata.stim$Group, astro.compdata.stim$Genotype)

# stats for onset times- astrocytes only during stim
OT.stim.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
OT.stim.model3 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
OT.stim.model4 = lmer(OnsetTime ~ Group_Gen + (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
OT.stim.model6 = lmer(OnsetTime ~ Group_Type_Gen + (1|Animal) + (1|Spot) + (1|trials),astro.compdata.stim,REML=FALSE)
OT.stim.anova <- anova(OT.stim.null, OT.stim.model1,OT.stim.model3,OT.stim.model4,OT.stim.model6)
print(OT.stim.anova)

OT.stim.Group_Gen<- glht(OT.stim.model4, mcp(Group_Gen= "Tukey"))
summary(OT.stim.Group_Gen)

#OT.stim.Group_type_Gen<- glht(OT.stim.model6, mcp(Group_Type_Gen= "Tukey"))
#summary(OT.stim.Group_type_Gen)

summary(OT.stim.model4)

# check residuals for linearity
plot(fitted(OT.stim.model4), residuals(OT.stim.model4),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(OT.stim.model4), residuals(OT.stim.model4)), col=46, lwd=2.5)





#############
# mean peak time
# mean onset times
df.PT1<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim",], measurevar = "peakTime", groupvars = c("Channel_Group", "Genotype"))
df.PT2<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",], measurevar = "peakTime", groupvars = c("Group", "Genotype"))
df.PT3<- summarySE(stim.lck.compdata.STIM, measurevar = "peakTime", groupvars = c("ROIType","Channel_Group", "Genotype"))

ggplot(df.PT1, aes(x=interaction(Genotype,Channel_Group),y=peakTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Peak Max Time (s)") +
  max.theme

ggplot(df.PT2, aes(x=Group,y=peakTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Peak Max Time (s)") +
  max.theme

ggplot(df.PT3, aes(x=interaction(Channel_Group,ROIType),y=peakTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Peak Max Time (s)") +
  max.theme



# stats
# stats for onset times- neurons vs astrocytes
PT.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials), lck.peaks.window,REML=FALSE)
PT.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), lck.peaks.window,REML=FALSE)
PT.model2 = lmer(peakTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), lck.peaks.window,REML=FALSE)
PT.model3 = lmer(peakTime ~ Group_Channel_Cond_Gen + (1|Animal) + (1|Spot) + (1|trials), lck.peaks.window,REML=FALSE)
PT.model4 = lmer(peakTime ~ Group_Channel_Type_Cond_Gen + (1|Animal) + (1|Spot) + (1|trials), lck.peaks.window,REML=FALSE)
PT.anova <- anova(PT.null, PT.model1,PT.model2,PT.model3,PT.model4)
print(PT.anova)

PT.Group_channel_Gen<- glht(PT.model3, mcp(Group_Channel_Cond_Gen= "Tukey"))
summary(PT.Group_channel_Gen)

#PT.Group_Channel_Type_Gen<- glht(PT.model4, mcp(Group_Channel_Type_Cond_Gen= "Tukey"))
#summary(PT.Group_Channel_Type_Gen)


# stats for onset times- only ASTROCYTES
PT.stim.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
PT.stim.model1 = lmer(peakTime ~ Genotype + (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
PT.stim.model3 = lmer(peakTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
PT.stim.model4 = lmer(peakTime ~ Group_Gen + (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
PT.stim.model6 = lmer(peakTime ~ Group_Type_Gen + (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
PT.stim.anova <- anova(PT.stim.null, PT.stim.model1,PT.stim.model3,PT.stim.model4,PT.stim.model6)
print(PT.stim.anova)

PT.stim.Group_Gen<- glht(PT.stim.model4, mcp(Group_Gen= "Tukey"))
summary(PT.stim.Group_Gen)

#PT.stim.Group_type_Gen<- glht(PT.stim.model6, mcp(Group_Type_Gen= "Tukey"))
#summary(PT.stim.Group_type_Gen)


########
#amplitude
df.amp1<-summarySE(lck.peaks.window, measurevar = "amplitude", groupvars = c("Channel","Condition", "Genotype"))
df.amp2<-summarySE(lck.peaks.window, measurevar = "amplitude", groupvars = c("Channel", "ROIType","Condition", "Genotype"))

df.amp3A<- summarySE(stim.lck.compdata, measurevar = "amplitude", groupvars = c("Channel_Group","Genotype","Condition"))
df.amp3B<- summarySE(stim.lck.compdata.STIM, measurevar = "amplitude", groupvars = c("Channel_Group","Genotype"))

df.amp4<- summarySE(astro.compdata.stim, measurevar = "amplitude", groupvars = c("Group","Genotype"))
df.amp5<- summarySE(astro.compdata.stim, measurevar = "amplitude", groupvars = c("Group","ROIType","Genotype"))


ggplot(df.amp1, aes(x=interaction(Channel,Genotype),y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  max.theme

ggplot(df.amp2, aes(x=interaction(Channel,interaction(Genotype, ROIType)),y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  max.theme

ggplot(df.amp3A, aes(x=interaction(Channel_Group,Genotype),y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  max.theme

ggplot(df.amp3B, aes(x=Channel_Group,y=amplitude, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude during stim trials") +
  max.theme

ggplot(df.amp4, aes(x=Group,y=amplitude, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude during stim trials") +
  max.theme


ggplot(df.amp5, aes(x=interaction(Group, ROIType),y=amplitude, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude during stim trials") +
  max.theme





#lck-GCaMP ONLY
amp.null.GC = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
amp.model2.GC  = lmer(amplitude ~ Genotype + (1|Animal) + (1|Spot) + (1|trials) , astro.compdata.stim,REML=FALSE)
amp.model3.GC  = lmer(amplitude ~ Group + (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
amp.model5.GC  = lmer(amplitude ~ Group_Gen + (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
amp.model6.GC  = lmer(amplitude ~ Group_Type_Gen + (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
amp.anova.GC  <- anova(amp.null.GC, amp.model2.GC,amp.model3.GC,amp.model5.GC,amp.model6.GC)
print(amp.anova.GC)

amp.Group_Gen.GC<- glht(amp.model5.GC, mcp(Group_Gen= "Tukey"))
summary(amp.Group_Gen.GC)

amp.Group_gen_ty.GC<- glht(amp.model6.GC, mcp(Group_Type_Gen= "Tukey"))
summary(amp.Group_gen_ty.GC)


# KNOCKOUTS have significantly lower amplitude that WT!
# still a significant difference between fast and delayed WTs



########
#duration
df.Dur1<-summarySE(lck.peaks.window, measurevar = "Duration", groupvars = c("Channel","Condition", "Genotype"))
df.Dur2<-summarySE(lck.peaks.window, measurevar = "Duration", groupvars = c("Channel", "ROIType","Condition", "Genotype"))

df.Dur3A<- summarySE(stim.lck.compdata, measurevar = "Duration", groupvars = c("Channel_Group","Genotype","Condition"))
df.Dur3B<- summarySE(stim.lck.compdata.STIM, measurevar = "Duration", groupvars = c("Channel_Group","Genotype"))

df.Dur4<- summarySE(astro.compdata.stim, measurevar = "Duration", groupvars = c("Group","Genotype"))
df.Dur5<- summarySE(astro.compdata.stim, measurevar = "Duration", groupvars = c("Group","ROIType","Genotype"))


ggplot(df.Dur1, aes(x=interaction(Channel,Genotype),y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  max.theme

ggplot(df.Dur2, aes(x=interaction(Channel,interaction(Genotype, ROIType)),y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  max.theme

ggplot(df.Dur3A, aes(x=interaction(Channel_Group,Genotype),y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  max.theme

ggplot(df.Dur3B, aes(x=Channel_Group,y=Duration, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration during stim trials") +
  max.theme

ggplot(df.Dur4, aes(x=Group,y=Duration, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration during stim trials") +
  max.theme

ggplot(df.Dur5, aes(x=interaction(Group, ROIType),y=Duration, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration during stim trials") +
  max.theme





#lck-GCaMP ONLY
Group_Genotype=interaction(stim.lck.compdata.STIM$Group[stim.lck.compdata.STIM$Channel=="GCaMP"],stim.lck.compdata.STIM$Genotype[stim.lck.compdata.STIM$Channel=="GCaMP"])
Group_Type_Genotype=interaction(stim.lck.compdata.STIM$Group[stim.lck.compdata.STIM$Channel=="GCaMP"],
                                stim.lck.compdata.STIM$Genotype[stim.lck.compdata.STIM$Channel=="GCaMP"],
                                stim.lck.compdata.STIM$ROIType[stim.lck.compdata.STIM$Channel=="GCaMP"])

Dur.null.GC = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
Dur.model2.GC  = lmer(Duration ~ Genotype + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
Dur.model3.GC  = lmer(Duration ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
Dur.model5.GC  = lmer(Duration ~ Group_Genotype + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
Dur.model6.GC  = lmer(Duration ~ Group_Type_Genotype + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
Dur.anova.GC  <- anova(Dur.null.GC, Dur.model2.GC,Dur.model3.GC,Dur.model5.GC,Dur.model6.GC)
print(Dur.anova.GC)

Dur.Group_Gen.GC<- glht(Dur.model5.GC, mcp(Group_Genotype= "Tukey"))
summary(Dur.Group_Gen.GC)

Dur.Group_gen_ty.GC<- glht(Dur.model6.GC, mcp(Group_Type_Genotype= "Tukey"))
summary(Dur.Group_gen_ty.GC)

######
#trace AUc 10s
df.tAUC101<-summarySE(lck.peaks.window, measurevar = "TraceAUC10", groupvars = c("Channel","Condition", "Genotype"), na.rm=TRUE)
df.tAUC102<-summarySE(lck.peaks.window, measurevar = "TraceAUC10", groupvars = c("Channel", "ROIType","Condition", "Genotype"),na.rm=TRUE)

df.tAUC103A<- summarySE(stim.lck.compdata, measurevar = "TraceAUC10", groupvars = c("Channel_Group","Genotype","Condition"),na.rm=TRUE)
df.tAUC103B<- summarySE(stim.lck.compdata.STIM, measurevar = "TraceAUC10", groupvars = c("Channel_Group","Genotype"),na.rm=TRUE)

df.tAUC104<- summarySE(astro.compdata.stim, measurevar = "TraceAUC10", groupvars = c("Group","Genotype"),na.rm=TRUE)
df.tAUC105<- summarySE(astro.compdata.stim, measurevar = "TraceAUC10", groupvars = c("Group","ROIType","Genotype"),na.rm=TRUE)


ggplot(df.tAUC101, aes(x=interaction(Channel,Genotype),y=TraceAUC10, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=TraceAUC10-se, ymax=TraceAUC10+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("TraceAUC10") +
  max.theme

ggplot(df.tAUC102, aes(x=interaction(Channel,interaction(Genotype, ROIType)),y=TraceAUC10, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=TraceAUC10-se, ymax=TraceAUC10+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("TraceAUC10") +
  max.theme

ggplot(df.tAUC103A, aes(x=interaction(Channel_Group,Genotype),y=TraceAUC10, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=TraceAUC10-se, ymax=TraceAUC10+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("TraceAUC10") +
  max.theme

ggplot(df.tAUC103B, aes(x=Channel_Group,y=TraceAUC10, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=TraceAUC10-se, ymax=TraceAUC10+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("TraceAUC10 tAUC10ing stim trials") +
  max.theme

ggplot(df.tAUC104, aes(x=Group,y=TraceAUC10, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=TraceAUC10-se, ymax=TraceAUC10+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("TraceAUC10 during stim trials") +
  max.theme

ggplot(df.tAUC105, aes(x=interaction(Group, ROIType),y=TraceAUC10, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=TraceAUC10-se, ymax=TraceAUC10+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("TraceAUC10 during stim trials") +
  max.theme





#lck-GCaMP ONLY
Group_Genotype=interaction(stim.lck.compdata.STIM$Group[stim.lck.compdata.STIM$Channel=="GCaMP"],stim.lck.compdata.STIM$Genotype[stim.lck.compdata.STIM$Channel=="GCaMP"])
Group_Type_Genotype=interaction(stim.lck.compdata.STIM$Group[stim.lck.compdata.STIM$Channel=="GCaMP"],
                                stim.lck.compdata.STIM$Genotype[stim.lck.compdata.STIM$Channel=="GCaMP"],
                                stim.lck.compdata.STIM$ROIType[stim.lck.compdata.STIM$Channel=="GCaMP"])

tAUC10.null.GC = lmer(TraceAUC10 ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
tAUC10.model2.GC  = lmer(TraceAUC10 ~ Genotype + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
tAUC10.model3.GC  = lmer(TraceAUC10 ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
tAUC10.model5.GC  = lmer(TraceAUC10 ~ Group_Genotype + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
tAUC10.model6.GC  = lmer(TraceAUC10 ~ Group_Type_Genotype + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
tAUC10.anova.GC  <- anova(tAUC10.null.GC, tAUC10.model2.GC,tAUC10.model3.GC,tAUC10.model5.GC,tAUC10.model6.GC)
print(tAUC10.anova.GC)

tAUC10.Group_Gen.GC<- glht(tAUC10.model5.GC, mcp(Group_Genotype= "Tukey"))
summary(tAUC10.Group_Gen.GC)

tAUC10.Group_gen_ty.GC<- glht(tAUC10.model6.GC, mcp(Group_Type_Genotype= "Tukey"))
summary(tAUC10.Group_gen_ty.GC)


######
#trace AUc 1s
df.tAUC11<-summarySE(lck.peaks.window, measurevar = "TraceAUC1", groupvars = c("Channel","Condition", "Genotype"), na.rm=TRUE)
df.tAUC12<-summarySE(lck.peaks.window, measurevar = "TraceAUC1", groupvars = c("Channel", "ROIType","Condition", "Genotype"),na.rm=TRUE)

df.tAUC13A<- summarySE(stim.lck.compdata, measurevar = "TraceAUC1", groupvars = c("Channel_Group","Genotype","Condition"),na.rm=TRUE)
df.tAUC13B<- summarySE(stim.lck.compdata.STIM, measurevar = "TraceAUC1", groupvars = c("Channel_Group","Genotype"),na.rm=TRUE)

df.tAUC14<- summarySE(astro.compdata.stim, measurevar = "TraceAUC1", groupvars = c("Group","Genotype"),na.rm=TRUE)
df.tAUC15<- summarySE(astro.compdata.stim, measurevar = "TraceAUC1", groupvars = c("Group","ROIType","Genotype"),na.rm=TRUE)


ggplot(df.tAUC11, aes(x=interaction(Channel,Genotype),y=TraceAUC1, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=TraceAUC1-se, ymax=TraceAUC1+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("TraceAUC1") +
  max.theme

ggplot(df.tAUC12, aes(x=interaction(Channel,interaction(Genotype, ROIType)),y=TraceAUC1, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=TraceAUC1-se, ymax=TraceAUC1+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("TraceAUC1") +
  max.theme

ggplot(df.tAUC13A, aes(x=interaction(Channel_Group,Genotype),y=TraceAUC1, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=TraceAUC1-se, ymax=TraceAUC1+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("TraceAUC1") +
  max.theme

ggplot(df.tAUC13B, aes(x=Channel_Group,y=TraceAUC1, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=TraceAUC1-se, ymax=TraceAUC1+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("TraceAUC1 tAUC1ing stim trials") +
  max.theme

ggplot(df.tAUC14, aes(x=Group,y=TraceAUC1, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=TraceAUC1-se, ymax=TraceAUC1+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("TraceAUC1 during stim trials") +
  max.theme

ggplot(df.tAUC15, aes(x=interaction(Group, ROIType),y=TraceAUC1, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=TraceAUC1-se, ymax=TraceAUC1+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("TraceAUC1 during stim trials") +
  max.theme





#lck-GCaMP ONLY
Group_Genotype=interaction(stim.lck.compdata.STIM$Group[stim.lck.compdata.STIM$Channel=="GCaMP"],stim.lck.compdata.STIM$Genotype[stim.lck.compdata.STIM$Channel=="GCaMP"])
Group_Type_Genotype=interaction(stim.lck.compdata.STIM$Group[stim.lck.compdata.STIM$Channel=="GCaMP"],
                                stim.lck.compdata.STIM$Genotype[stim.lck.compdata.STIM$Channel=="GCaMP"],
                                stim.lck.compdata.STIM$ROIType[stim.lck.compdata.STIM$Channel=="GCaMP"])

tAUC1.null.GC = lmer(TraceAUC1 ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
tAUC1.model2.GC  = lmer(TraceAUC1 ~ Genotype + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
tAUC1.model3.GC  = lmer(TraceAUC1 ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
tAUC1.model5.GC  = lmer(TraceAUC1 ~ Group_Genotype + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
tAUC1.model6.GC  = lmer(TraceAUC1 ~ Group_Type_Genotype + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
tAUC1.anova.GC  <- anova(tAUC1.null.GC, tAUC1.model2.GC,tAUC1.model3.GC,tAUC1.model5.GC,tAUC1.model6.GC)
print(tAUC1.anova.GC)

tAUC1.Group_Gen.GC<- glht(tAUC1.model5.GC, mcp(Group_Genotype= "Tukey"))
summary(tAUC1.Group_Gen.GC)

tAUC1.Group_gen_ty.GC<- glht(tAUC1.model6.GC, mcp(Group_Type_Genotype= "Tukey"))
summary(tAUC1.Group_gen_ty.GC)


#######
#peakAUC
df.peakAUC1<-summarySE(lck.peaks.window, measurevar = "peakAUC", groupvars = c("Channel","Condition", "Genotype"), na.rm=TRUE)
df.peakAUC2<-summarySE(lck.peaks.window, measurevar = "peakAUC", groupvars = c("Channel", "ROIType","Condition", "Genotype"),na.rm=TRUE)

df.peakAUC3A<- summarySE(stim.lck.compdata, measurevar = "peakAUC", groupvars = c("Channel_Group","Genotype","Condition"),na.rm=TRUE)
df.peakAUC3B<- summarySE(stim.lck.compdata.STIM, measurevar = "peakAUC", groupvars = c("Channel_Group","Genotype"),na.rm=TRUE)

df.peakAUC4<- summarySE(astro.compdata.stim, measurevar = "peakAUC", groupvars = c("Group","Genotype"),na.rm=TRUE)
df.peakAUC5<- summarySE(astro.compdata.stim, measurevar = "peakAUC", groupvars = c("Group","ROIType","Genotype"),na.rm=TRUE)


ggplot(df.peakAUC1, aes(x=interaction(Channel,Genotype),y=peakAUC, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakAUC-se, ymax=peakAUC+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakAUC") +
  max.theme

ggplot(df.peakAUC2, aes(x=interaction(Channel,interaction(Genotype, ROIType)),y=peakAUC, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakAUC-se, ymax=peakAUC+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakAUC") +
  max.theme

ggplot(df.peakAUC3A, aes(x=interaction(Channel_Group,Genotype),y=peakAUC, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakAUC-se, ymax=peakAUC+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakAUC") +
  max.theme

ggplot(df.peakAUC3B, aes(x=Channel_Group,y=peakAUC, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakAUC-se, ymax=peakAUC+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakAUC peakAUCing stim trials") +
  max.theme

ggplot(df.peakAUC4, aes(x=Group,y=peakAUC, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakAUC-se, ymax=peakAUC+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakAUC during stim trials") +
  max.theme

ggplot(df.peakAUC5, aes(x=interaction(Group, ROIType),y=peakAUC, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakAUC-se, ymax=peakAUC+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakAUC during stim trials") +
  max.theme





#lck-GCaMP ONLY
Group_Genotype=interaction(stim.lck.compdata.STIM$Group[stim.lck.compdata.STIM$Channel=="GCaMP"],stim.lck.compdata.STIM$Genotype[stim.lck.compdata.STIM$Channel=="GCaMP"])
Group_Type_Genotype=interaction(stim.lck.compdata.STIM$Group[stim.lck.compdata.STIM$Channel=="GCaMP"],
                                stim.lck.compdata.STIM$Genotype[stim.lck.compdata.STIM$Channel=="GCaMP"],
                                stim.lck.compdata.STIM$ROIType[stim.lck.compdata.STIM$Channel=="GCaMP"])

peakAUC.null.GC = lmer(peakAUC ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
peakAUC.model2.GC  = lmer(peakAUC ~ Genotype + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
peakAUC.model3.GC  = lmer(peakAUC ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
peakAUC.model5.GC  = lmer(peakAUC ~ Group_Genotype + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
peakAUC.model6.GC  = lmer(peakAUC ~ Group_Type_Genotype + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
peakAUC.anova.GC  <- anova(peakAUC.null.GC, peakAUC.model2.GC,peakAUC.model3.GC,peakAUC.model5.GC,peakAUC.model6.GC)
print(peakAUC.anova.GC)

peakAUC.Group_Gen.GC<- glht(peakAUC.model5.GC, mcp(Group_Genotype= "Tukey"))
summary(peakAUC.Group_Gen.GC)

peakAUC.Group_gen_ty.GC<- glht(peakAUC.model6.GC, mcp(Group_Type_Genotype= "Tukey"))
summary(peakAUC.Group_gen_ty.GC)


###### 
# proportion of ROIs that are fast or delayed

#what about proportion of each group based on ROI type??

GCaMP.KO<-subset(lck.peaks.window,Channel=="GCaMP" & Genotype=="IP3R2_KO" & Condition=="Stim")
GCaMP.WT<-subset(lck.peaks.window,Channel=="GCaMP" & Genotype=="IP3R2_WT" & Condition=="Stim")
RCaMP.KO<-subset(lck.peaks.window,Channel=="RCaMP" & Genotype=="IP3R2_KO" & Condition=="Stim")
RCaMP.WT<-subset(lck.peaks.window,Channel=="RCaMP" & Genotype=="IP3R2_WT" & Condition=="Stim")

allROIs.KO<-length(unique(GCaMP.KO$ROIs_trial))
fastROIs.KO<-length(unique(GCaMP.KO$ROIs_trial[GCaMP.KO$Group=="fast_MDs"]))
delayedROIs.KO<-length(unique(GCaMP.KO$ROIs_trial[GCaMP.KO$Group=="delayed_MDs"]))

allROIs.WT<-length(unique(GCaMP.WT$ROIs_trial))
fastROIs.WT<-length(unique(GCaMP.WT$ROIs_trial[GCaMP.WT$Group=="fast_MDs"]))
delayedROIs.WT<-length(unique(GCaMP.WT$ROIs_trial[GCaMP.WT$Group=="delayed_MDs"]))

#proportion
propfast.KO=fastROIs.KO/allROIs.KO
propdelayed.KO=delayedROIs.KO/allROIs.KO

propfast.WT=fastROIs.WT/allROIs.WT
propdelayed.WT=delayedROIs.WT/allROIs.WT


#ef vs proc  wildtypes
allROIs.ef.WT<-length(unique(GCaMP.WT$ROIs_trial[GCaMP.WT$ROIType=="Endfoot"]))
fastROIs.ef.WT<-length(unique(GCaMP.WT$ROIs_trial[GCaMP.WT$Group=="fast_MDs"&GCaMP.WT$ROIType=="Endfoot"]))
delayedROIs.ef.WT<-length(unique(GCaMP.WT$ROIs_trial[GCaMP.WT$Group=="delayed_MDs"&GCaMP.WT$ROIType=="Endfoot"]))

allROIs.proc.WT<-length(unique(GCaMP.WT$ROIs_trial[GCaMP.WT$ROIType=="Process"]))
fastROIs.proc.WT<-length(unique(GCaMP.WT$ROIs_trial[GCaMP.WT$Group=="fast_MDs"&GCaMP.WT$ROIType=="Process"]))
delayedROIs.proc.WT<-length(unique(GCaMP.WT$ROIs_trial[GCaMP.WT$Group=="delayed_MDs"&GCaMP.WT$ROIType=="Process"]))


propfast.ef.WT=fastROIs.ef.WT/allROIs.ef.WT
propdelayed.ef.WT=delayedROIs.ef.WT/allROIs.ef.WT

propfast.proc.WT=fastROIs.proc.WT/allROIs.proc.WT
propdelayed.proc.WT=delayedROIs.proc.WT/allROIs.proc.WT

#ef vs proc knockouts
allROIs.ef.KO<-length(unique(GCaMP.KO$ROIs_trial[GCaMP.KO$ROIType=="Endfoot"]))
fastROIs.ef.KO<-length(unique(GCaMP.KO$ROIs_trial[GCaMP.KO$Group=="fast_MDs"&GCaMP.KO$ROIType=="Endfoot"]))
delayedROIs.ef.KO<-length(unique(GCaMP.KO$ROIs_trial[GCaMP.KO$Group=="delayed_MDs"&GCaMP.KO$ROIType=="Endfoot"]))

allROIs.proc.KO<-length(unique(GCaMP.KO$ROIs_trial[GCaMP.KO$ROIType=="Process"]))
fastROIs.proc.KO<-length(unique(GCaMP.KO$ROIs_trial[GCaMP.KO$Group=="fast_MDs"&GCaMP.KO$ROIType=="Process"]))
delayedROIs.proc.KO<-length(unique(GCaMP.KO$ROIs_trial[GCaMP.KO$Group=="delayed_MDs"&GCaMP.KO$ROIType=="Process"]))


propfast.ef.KO=fastROIs.ef.KO/allROIs.ef.KO
propdelayed.ef.KO=delayedROIs.ef.KO/allROIs.ef.KO

propfast.proc.KO=fastROIs.proc.KO/allROIs.proc.KO
propdelayed.proc.KO=delayedROIs.proc.KO/allROIs.proc.KO



