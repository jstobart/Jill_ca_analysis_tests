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

all.lck.peaks <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/Peaks_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
all.lck.OT<-read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/OnsetTimes_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")


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




#################
# neuronal responses to stimulation
NeuronalStim<-subset(all.lck.OT, Channel=="RCaMP" & Condition=="Stim")

# should have an onset time in 8 s stimulus
ggplot(NeuronalStim[NeuronalStim$OnsetTime<10,],aes(x=OnsetTime,y=..density..,fill=Genotype)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP onset times between 0 and 10 s from stim trials")+
  max.theme


Neuron95Onset<-quantile(NeuronalStim$OnsetTime[NeuronalStim$OnsetTime<8], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
print(Neuron95Onset)


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
# for median and mean calculations

# no stim vs 8 s stim- 
#neuronal window=9s, AC window= 15 s for peak time, 
#neuronal window=2s, AC window=12 s for onset


LongN_PTwind=8
LongAC_PTwind=WT.Astro.PT80

LongN_OTwind=8
LongAC_OTwind=WT.Astro.OT80

# remove data that is outside the above windows
# lck
stim.lck.OT.R<-subset(all.lck.OT, Channel=="RCaMP" & OnsetTime<=LongN_OTwind)
stim.lck.OT.G<-subset(all.lck.OT, Channel=="GCaMP" & OnsetTime<=LongAC_OTwind)

stim.lck.OT.window<-rbind(stim.lck.OT.R, stim.lck.OT.G)


# peak times
# lck
stim.lck.PT.R<-subset(all.lck.peaks, Channel=="RCaMP" & peakTime<=LongN_PTwind & peakTime>=0 & Duration<45)
stim.lck.PT.G<-subset(all.lck.peaks, Channel=="GCaMP" & peakTime<=LongAC_PTwind & peakTime>=0 & Duration<45)

stim.lck.peaks.window<-rbind(stim.lck.PT.R, stim.lck.PT.G)

rm(stim.lck.OT.G,stim.lck.OT.R,stim.lck.PT.G, stim.lck.PT.R)
#######
# peak times only for ROIs with onset times?

stim.lck.alldata<-merge(stim.lck.peaks.window, stim.lck.OT.window[, c("ROIs_Cond", "OnsetTime","TraceAUC1","TraceAUC10")], by="ROIs_Cond")


# onset times
ggplot(stim.lck.alldata, aes(x=interaction(Channel, Genotype),y=OnsetTime, fill= Condition)) +
  geom_boxplot(notch=TRUE)+
  ylab("OnsetTimes Time (s)") +
  ggtitle("matched data= peak times and onset times")+
  max.theme

#peak times
ggplot(stim.lck.alldata, aes(x=interaction(Channel, Genotype),y=peakTime, fill= Condition)) +
  geom_boxplot(notch=TRUE)+
  ylab("peak Time (s)") +
  ggtitle("matched data= peak times and onset times")+
  max.theme

######

# identify "FAST" astrocytes
stim.lck.alldata$Group<-0
stim.lck.alldata$Group[stim.lck.alldata$OnsetTime<1]<-"fast_MDs"
stim.lck.alldata$Group[stim.lck.alldata$OnsetTime>=1]<-"delayed_MDs"

stim.lck.alldata$Group <- factor(stim.lck.alldata$Group, levels = c("fast_MDs","delayed_MDs"))
stim.lck.alldata$Channel <- factor(stim.lck.alldata$Channel, levels = c("RCaMP","GCaMP"))


# fast vs delayed
stim.lck.alldata$Channel_Group<-interaction(stim.lck.alldata$Channel, stim.lck.alldata$Group)
stim.lck.alldata$Channel_Group<-as.factor(stim.lck.alldata$Channel_Group)

stim.lck.compdata<-stim.lck.alldata[!(stim.lck.alldata$Channel=="RCaMP"& stim.lck.alldata$Group=="delayed_MDs"),]

# take out the effect of Condition
# we are only interested in stim case
stim.lck.compdata.STIM<-subset(stim.lck.compdata, Condition=="Stim")

###### 
# proportion of ROIs that are fast or delayed

#what about proportion of each group based on ROI type??

GCaMP.KO<-subset(stim.lck.alldata,Channel=="GCaMP" & Genotype=="IP3R2_KO" & Condition=="Stim")
GCaMP.WT<-subset(stim.lck.alldata,Channel=="GCaMP" & Genotype=="IP3R2_WT" & Condition=="Stim")
RCaMP.KO<-subset(stim.lck.alldata,Channel=="RCaMP" & Genotype=="IP3R2_KO" & Condition=="Stim")
RCaMP.WT<-subset(stim.lck.alldata,Channel=="RCaMP" & Genotype=="IP3R2_WT" & Condition=="Stim")

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

######

# fraction of active pixels from all the astrocyte pixels
# with peaks near stimulus

# only consider active based ROIs (i.e. processes)

stim.lck.processes<-subset(stim.lck.alldata, Channel=="GCaMP" & ROIType=="Process")

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
# number of ROIs in each trial for each field of view (during the 8s stimulus)

ROInum.lck<-ddply(stim.lck.alldata, c("Animal","Spot","Condition","Channel","Genotype","trials"), summarise, nROIs=length(OnsetTime))

# add in number of trials
ROInum.lck$Ani_Spot_Cond<-paste(ROInum.lck$Animal, ROInum.lck$Spot, ROInum.lck$Condition, sep="_")
ROInum.lck<-merge(ROInum.lck, Spot.lck.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)
ROInum.lck$ROIsPerTrial<-ROInum.lck$nROIs/ROInum.lck$nTrials

# mean
df.lck.ROInum<-summarySE(ROInum.lck, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition", "Genotype"))


ggplot(df.lck.ROInum, aes(x=interaction(Channel,Genotype),y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view") +
  max.theme

ggplot(df.lck.ROInum[df.lck.ROInum$Channel=="GCaMP",], aes(x=Genotype,y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view") +
  max.theme

Condition_Channel2= interaction(ROInum.lck$Condition,ROInum.lck$Channel)
Condition_Channel_Geno= interaction(ROInum.lck$Condition,ROInum.lck$Channel,ROInum.lck$Genotype)
nROI.lck.stim.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.lck,REML=FALSE)
nROI.lck.stim.model1 = lmer(ROIsPerTrial~ Channel + (1|Animal), ROInum.lck,REML=FALSE)
nROI.lck.stim.model2A = lmer(ROIsPerTrial ~ Condition + (1|Animal), ROInum.lck,REML=FALSE)
nROI.lck.stim.model2B = lmer(ROIsPerTrial ~ Genotype + (1|Animal), ROInum.lck,REML=FALSE)
nROI.lck.stim.model3 = lmer(ROIsPerTrial ~ Condition_Channel2+ Genotype + (1|Animal), ROInum.lck,REML=FALSE)
nROI.lck.stim.model4 = lmer(ROIsPerTrial ~ Condition_Channel_Geno + (1|Animal), ROInum.lck,REML=FALSE)
nROI.lck.stim.anova <- anova(nROI.lck.stim.null, nROI.lck.stim.model1,nROI.lck.stim.model2A,
                             nROI.lck.stim.model2B,nROI.lck.stim.model3,nROI.lck.stim.model4)
print(nROI.lck.stim.anova)

nROI.lck.stim.Cond_Channel_Geno<- glht(nROI.lck.stim.model4, mcp(Condition_Channel_Geno= "Tukey"))
summary(nROI.lck.stim.Cond_Channel_Geno)


# only GCaMP
Condition_Geno= interaction(ROInum.lck$Condition[ROInum.lck$Channel=="GCaMP"],ROInum.lck$Genotype[ROInum.lck$Channel=="GCaMP"])
nROI.lck.stim.null.gc = lmer(ROIsPerTrial ~ (1|Animal), ROInum.lck[ROInum.lck$Channel=="GCaMP",],REML=FALSE)
nROI.lck.stim.model1.gc = lmer(ROIsPerTrial ~ Condition + (1|Animal), ROInum.lck[ROInum.lck$Channel=="GCaMP",],REML=FALSE)
nROI.lck.stim.model2.gc = lmer(ROIsPerTrial ~ Genotype + (1|Animal), ROInum.lck[ROInum.lck$Channel=="GCaMP",],REML=FALSE)
nROI.lck.stim.model3.gc = lmer(ROIsPerTrial ~ Condition_Geno + (1|Animal), ROInum.lck[ROInum.lck$Channel=="GCaMP",],REML=FALSE)
nROI.lck.stim.anova.gc <- anova(nROI.lck.stim.null.gc, nROI.lck.stim.model1.gc,nROI.lck.stim.model2.gc,nROI.lck.stim.model3.gc)
print(nROI.lck.stim.anova.gc)

nROI.lck.stim.Cond_Geno.gc<- glht(nROI.lck.stim.model3.gc, mcp(Condition_Geno= "Tukey"))
summary(nROI.lck.stim.Cond_Geno.gc)

#################
# number of fastROIs per trial
FastROInum.lck.stim<-ddply(stim.lck.alldata[stim.lck.alldata$Condition=="Stim",], c("Animal","Spot","Genotype","Condition","Channel","Group","Channel_Group"), summarise, nROIs=length(OnsetTime))

# add in number of trials
FastROInum.lck.stim$Ani_Spot_Cond<-paste(FastROInum.lck.stim$Animal, FastROInum.lck.stim$Spot, FastROInum.lck.stim$Condition, sep="_")

FastROInum.lck.stim<-merge(FastROInum.lck.stim, Spot.lck.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)

FastROInum.lck.stim$ROIsPerTrial<-FastROInum.lck.stim$nROIs/FastROInum.lck.stim$nTrials

# mean number of total ROIs per trial
df.lck.FastROInum.mean<-summarySE(FastROInum.lck.stim, measurevar = "ROIsPerTrial", groupvars = c("Channel_Group", "Genotype"))
df.lck.FastROInum.GC<-summarySE(FastROInum.lck.stim[FastROInum.lck.stim$Channel=="GCaMP",], measurevar = "ROIsPerTrial", groupvars = c("Group", "Genotype"))


ggplot(df.lck.FastROInum.mean, aes(x=Channel_Group,y=ROIsPerTrial, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during stimulation") +
  max.theme

ggplot(df.lck.FastROInum.GC, aes(x=Group,y=ROIsPerTrial, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during stimulation") +
  max.theme

#only astrocytes
Group_Genotype= interaction(FastROInum.lck.stim$Group[FastROInum.lck.stim$Channel=="GCaMP"],FastROInum.lck.stim$Genotype[FastROInum.lck.stim$Channel=="GCaMP"])
FastnROI.lck.stim.null = lmer(ROIsPerTrial ~ (1|Animal), FastROInum.lck.stim[FastROInum.lck.stim$Channel=="GCaMP",],REML=FALSE)
FastnROI.lck.stim.model1 = lmer(ROIsPerTrial~ Group + (1|Animal), FastROInum.lck.stim[FastROInum.lck.stim$Channel=="GCaMP",],REML=FALSE)
FastnROI.lck.stim.model2 = lmer(ROIsPerTrial ~ Genotype + (1|Animal), FastROInum.lck.stim[FastROInum.lck.stim$Channel=="GCaMP",],REML=FALSE)
FastnROI.lck.stim.model3 = lmer(ROIsPerTrial ~ Group + Genotype + (1|Animal), FastROInum.lck.stim[FastROInum.lck.stim$Channel=="GCaMP",],REML=FALSE)
FastnROI.lck.stim.model4 = lmer(ROIsPerTrial ~ Group_Genotype + (1|Animal), FastROInum.lck.stim[FastROInum.lck.stim$Channel=="GCaMP",],REML=FALSE)
FastnROI.lck.stim.anova <- anova(FastnROI.lck.stim.null, FastnROI.lck.stim.model1,FastnROI.lck.stim.model2,
                                 FastnROI.lck.stim.model3,FastnROI.lck.stim.model4)
print(FastnROI.lck.stim.anova)

FastnROI.lck.stim.Group_Genotype<- glht(FastnROI.lck.stim.model4, mcp(Group_Genotype= "Tukey"))
summary(FastnROI.lck.stim.Group_Genotype)


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
Group_Channel_Type_Cond_Gen=interaction(stim.lck.alldata$Group,stim.lck.alldata$Channel,stim.lck.alldata$ROIType, 
                                        stim.lck.alldata$Genotype,stim.lck.alldata$Condition)
Group_Channel_Cond_Gen=interaction(stim.lck.alldata$Group,stim.lck.alldata$Channel,stim.lck.alldata$Genotype,
                                   stim.lck.alldata$Condition)

# stats for onset times- neurons vs astrocytes
OT.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model2 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model3 = lmer(OnsetTime ~ Group_Channel_Cond_Gen + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model4 = lmer(OnsetTime ~ Group_Channel_Type_Cond_Gen + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.anova <- anova(OT.null, OT.model1,OT.model2,OT.model3,OT.model4)
print(OT.anova)

OT.Group_channel_Gen<- glht(OT.model3, mcp(Group_Channel_Cond_Gen= "Tukey"))
summary(OT.Group_channel_Gen)

OT.Group_Channel_Type_Gen<- glht(OT.model4, mcp(Group_Channel_Type_Cond_Gen= "Tukey"))
summary(OT.Group_Channel_Type_Gen)


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

OT.stim.Group_type_Gen<- glht(OT.stim.model6, mcp(Group_Type_Gen= "Tukey"))
summary(OT.stim.Group_type_Gen)

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
PT.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
PT.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
PT.model2 = lmer(peakTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
PT.model3 = lmer(peakTime ~ Group_Channel_Cond_Gen + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
PT.model4 = lmer(peakTime ~ Group_Channel_Type_Cond_Gen + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
PT.anova <- anova(PT.null, PT.model1,PT.model2,PT.model3,PT.model4)
print(PT.anova)

PT.Group_channel_Gen<- glht(PT.model3, mcp(Group_Channel_Cond_Gen= "Tukey"))
summary(PT.Group_channel_Gen)

PT.Group_Channel_Type_Gen<- glht(PT.model4, mcp(Group_Channel_Type_Cond_Gen= "Tukey"))
summary(PT.Group_Channel_Type_Gen)


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

PT.stim.Group_ype_Gen<- glht(PT.stim.model6, mcp(Group_Type_Gen= "Tukey"))
summary(PT.stim.Group_ype_Gen)


########
#amplitude
df.amp1<-summarySE(stim.lck.alldata, measurevar = "amplitude", groupvars = c("Channel","Condition", "Genotype"))
df.amp2<-summarySE(stim.lck.alldata, measurevar = "amplitude", groupvars = c("Channel", "ROIType","Condition", "Genotype"))

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
df.Dur1<-summarySE(stim.lck.alldata, measurevar = "Duration", groupvars = c("Channel","Condition", "Genotype"))
df.Dur2<-summarySE(stim.lck.alldata, measurevar = "Duration", groupvars = c("Channel", "ROIType","Condition", "Genotype"))

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
# peak time
df.pT1<-summarySE(stim.both.alldata, measurevar = "peakTime", groupvars = c("Channel", "Condition"))
df.pT2<-summarySE(stim.both.alldata, measurevar = "peakTime", groupvars = c("Channel", "Group","Condition"))
df.pT3<- summarySE(stim.lck.alldata, measurevar = "peakTime", groupvars = c("Channel_Group","Condition"))
df.pT4<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim",], measurevar = "peakTime", groupvars = c("Channel_Group"))
df.pT5<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim"& stim.lck.compdata$Channel=="lck_GCaMP",], measurevar = "peakTime", groupvars = c("ROIType","Channel_Group"))

ggplot(df.pT1, aes(x=Channel,y=peakTime, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakTime") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.pT2, aes(x=interaction(Channel,Group),y=peakTime, fill= Condition)) +
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

ggplot(stim.lck.compdata, aes(x=Channel_Group,y=peakTime, fill=Channel_Group)) +
  geom_boxplot(notch=TRUE)+
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


#######

# Process ROI area
df.Rarea1<-summarySE(stim.both.alldata[stim.both.alldata$ROIType=="Process",], measurevar = "area", groupvars = c("Channel", "Condition"))
df.Rarea2<-summarySE(stim.both.alldata[stim.both.alldata$ROIType=="Process",], measurevar = "area", groupvars = c("Channel", "Group","Condition"))
df.Rarea3<- summarySE(stim.lck.compdata[stim.lck.compdata$ROIType=="Process",], measurevar = "area", groupvars = c("Channel_Group","Condition"))
df.Rarea4<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim"& stim.lck.compdata$ROIType=="Process",], measurevar = "area", groupvars = c("Channel_Group"))
df.Rarea6<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="lck_RCaMP"& stim.lck.compdata.STIM$ROIType=="Process",], measurevar = "area", groupvars = c("Group"))

ggplot(df.Rarea1, aes(x=Channel,y=area, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("area") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.Rarea2, aes(x=interaction(Channel,Group),y=area, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("area") +
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

ggplot(df.Rarea6, aes(x=Group,y=area, fill= Group)) +
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

ggplot(stim.lck.alldata[stim.lck.alldata$Channel_Group=="lck_GCaMP.fast",], aes(x=area, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=5, position="dodge")+
  ylab("density") +
  xlim(-2,200)+
  ggtitle("lck data- fast no stim vs stim")+
  max.theme

#only consider STIM case
area.null.stim = lmer(area ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="lck_RCaMP" & stim.lck.compdata.STIM$ROIType=="Process",],REML=FALSE)
area.model1.stim  = lmer(area ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="lck_RCaMP" & stim.lck.compdata.STIM$ROIType=="Process",],REML=FALSE)

area.anova.stim  <- anova(area.null.stim, area.model1.stim)
print(area.anova.stim)

area.stim.group<- glht(area.model1.stim, mcp(Group= "Tukey"))
summary(area.stim.group)



