
library("lme4")
library("lmerTest")
library("lattice")
library("plyr")
library("ggplot2")
library("gplots")
library("lsmeans")
library("Rmisc")
library("MASS")
library("multcomp")
library("reshape2")
library("data.table")
library("Hmisc")
library("stringr")
library("spatstat")
library('hexbin')

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

# COLOUR BLIND FRIENDLY PALETTE FOR PLOTS
# The palette with black:
cbbPalette <- c("#000000","#D55E00","#009E73","#E69F00","#56B4E9","#CC79A7","#F0E442")


###########
# NOTES
# cyto data- onset time histogram, peak time histogram, duration, amplitude of peaks
# Fast peaks?
# rise time? from onset time to peak time
# ROI area of dendrites and processes?
# distance calculations (could they be close together?- consider the onset time comparisons)

# DSP4- doesn't change the response to stimulation

#Lck- fa

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
longstim.OT.comp <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/longstim_firstonset_comparisons.csv", header=TRUE, sep = ",")
shortstim.OT.comp <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/shortstim_firstonset_comparisons.csv", header=TRUE, sep = ",")
nostim.OT.comp <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/nostim_firstonset_comparisons.csv", header=TRUE, sep = ",")

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

## NOTE: define a stimulation window: 30 s after stim for both sensors

stimwindow=30


# onset time distributions
all.lck.OT<-subset(all.lck.OT,OnsetTime<stimwindow)
all.cyto.OT<-subset(all.cyto.OT,OnsetTime<stimwindow)
all.cyto.DSP4.OT<-subset(all.cyto.DSP4.OT,OnsetTime<stimwindow)

ntrials.lck.OT<- ddply(all.lck.OT, c("Condition"), summarise, ntrials=length(unique(Spot_trial)))
ntrials.cyto.OT<- ddply(all.cyto.OT, c("Condition"), summarise, ntrials=length(unique(Spot_trial)))
ntrials.cyto.DSP4.OT<- ddply(all.cyto.DSP4.OT, c("Condition"), summarise, ntrials=length(unique(Spot_trial)))

#histogram bins
histseq= seq(0,stimwindow,1)

# neuronal lck onset histogram
# counts for each condition in the histogram
Nostim=hist(all.lck.OT$OnsetTime[(all.lck.OT$Channel=="RCaMP" & all.lck.OT$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim=hist(all.lck.OT$OnsetTime[(all.lck.OT$Channel=="RCaMP" & all.lck.OT$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts
Shortstim=hist(all.lck.OT$OnsetTime[(all.lck.OT$Channel=="RCaMP" & all.lck.OT$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim=Nostim/ntrials.lck.OT$ntrials[(ntrials.lck.OT$Condition=="Nostim")]
Stim=Stim/ntrials.lck.OT$ntrials[(ntrials.lck.OT$Condition=="Stim")]
Shortstim=Shortstim/ntrials.lck.OT$ntrials[(ntrials.lck.OT$Condition=="shortstim")]
                     
#make a data frame for plotting
Neuron.lck.histo <- data.frame(cbind(Nostim, Stim, Shortstim))
Neuron.lck.histo$time<-histseq[2:length(histseq)]
Neuron.lck.histo2<-melt(Neuron.lck.histo,id="time")

ggplot(data=Neuron.lck.histo2, aes(x=time, y= value, colour=variable)) + 
  geom_line(size=1)+
  xlab("Onset time (s)") +
  ylab("number of ROIs per trial")+
  ggtitle("Neurons- onset time- ROIs per trial- lck data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

# astrocyte lck onset histogram
#histogram bins
histseq2= seq(0,stimwindow, 1)

# counts for each condition in the histogram
Nostim=hist(all.lck.OT$OnsetTime[(all.lck.OT$Channel=="GCaMP" & all.lck.OT$Condition=="Nostim")], breaks=histseq2, plot=FALSE)$counts
Stim=hist(all.lck.OT$OnsetTime[(all.lck.OT$Channel=="GCaMP" & all.lck.OT$Condition=="Stim")], breaks=histseq2, plot=FALSE)$counts
Shortstim=hist(all.lck.OT$OnsetTime[(all.lck.OT$Channel=="GCaMP" & all.lck.OT$Condition=="shortstim")], breaks=histseq2, plot=FALSE)$counts
  
# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim=Nostim/ntrials.lck.OT$ntrials[(ntrials.lck.OT$Condition=="Nostim")]
Stim=Stim/ntrials.lck.OT$ntrials[(ntrials.lck.OT$Condition=="Stim")]
Shortstim=Shortstim/ntrials.lck.OT$ntrials[(ntrials.lck.OT$Condition=="shortstim")]
  
#make a data frame for plotting
AC.lck.histo<- data.frame(cbind(Nostim, Stim, Shortstim))
AC.lck.histo$time<-histseq2[2:length(histseq2)]
AC.lck.histo2<-melt(AC.lck.histo,id="time")
  
ggplot(data=AC.lck.histo2, aes(x=time, y= value, colour=variable)) + 
    geom_line(size=1)+
    xlab("Onset time (s)") +
    ylab("number of ROIs per trial")+
    ggtitle("Astrocytes- onset time- ROIs per trial- lck data") + 
    scale_colour_manual(values=cbbPalette)+
  max.theme

###
# neurons and astrocytes together by density to account for different number of ROIs

ggplot(all.lck.OT[(all.lck.OT$Condition!="Nostim"),], aes(x=OnsetTime, y=..density.., colour = Channel, linetype=Condition)) +
  geom_freqpoly(binwidth = 0.5, lwd=1)+
  ggtitle("stim- neurons vs astrocytes-lck data") + 
  xlab("Onset Time (s)") + 
  xlim(-0.5,30)+
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

ggplot(all.lck.OT[(all.lck.OT$Condition!="Nostim"),], aes(x=OnsetTime, colour = Channel, linetype=Condition)) + 
  stat_ecdf() +
  ggtitle("cdf- neurons vs astrocytes-lck data- stim types") + 
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

ggplot(all.lck.OT[(all.lck.OT$Condition!="Nostim"),], aes(x=OnsetTime, colour = Channel)) + 
  stat_ecdf(lwd=1) +
  ggtitle("cdf- neurons vs astrocytes-lck data") + 
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme


# zoomed in plot of onset times- identify "FAST" astrocytes
stim.lck.OT<-subset(all.lck.OT,Condition!="Nostim")
stim.lck.OT$Group<-0
stim.lck.OT$Group[stim.lck.OT$OnsetTime<1]<-"fast"
stim.lck.OT$Group[stim.lck.OT$OnsetTime>=1]<-"delayed"

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
# stats for onset times- neurons vs astrocytes
OT.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial), stim.lck.OT,REML=FALSE)
OT.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.lck.OT,REML=FALSE)
OT.model2 = lmer(OnsetTime ~ Condition + (1|Animal) + (1|Spot) + (1|Spot_trial) , stim.lck.OT,REML=FALSE)
OT.model3 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.lck.OT,REML=FALSE)
OT.model4 = lmer(OnsetTime ~ Group_Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.lck.OT,REML=FALSE)
OT.model5 = lmer(OnsetTime ~ Condition_Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), stim.lck.OT,REML=FALSE)
OT.anova <- anova(OT.null, OT.model1,OT.model2,OT.model3,OT.model4,OT.model5)
print(OT.anova)

OT.Group_channel<- glht(OT.model4, mcp(Group_Channel= "Tukey"))
summary(OT.Group_channel)

###### 
#cyto data
# neuronal onset histogram
#histogram bins
histseq= seq(0,stimwindow,1)
# counts for each condition in the histogram
Nostim=hist(all.cyto.OT$OnsetTime[(all.cyto.OT$Channel=="RCaMP" & all.cyto.OT$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim=hist(all.cyto.OT$OnsetTime[(all.cyto.OT$Channel=="RCaMP" & all.cyto.OT$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts
Shortstim=hist(all.cyto.OT$OnsetTime[(all.cyto.OT$Channel=="RCaMP" & all.cyto.OT$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim=Nostim/ntrials.cyto.OT$ntrials[(ntrials.cyto.OT$Condition=="Nostim")]
Stim=Stim/ntrials.cyto.OT$ntrials[(ntrials.cyto.OT$Condition=="Stim")]
Shortstim=Shortstim/ntrials.cyto.OT$ntrials[(ntrials.cyto.OT$Condition=="shortstim")]

#make a data frame for plotting
Neuron.cyto.histo <- data.frame(cbind(Nostim, Stim, Shortstim))
Neuron.cyto.histo$time<-histseq[2:length(histseq)]
Neuron.cyto.histo<-melt(Neuron.cyto.histo,id="time")

ggplot(data=Neuron.cyto.histo, aes(x=time, y= value, colour=variable)) + 
  geom_line(size=1)+
  xlab("Onset time (s)") +
  ylab("number of ROIs per trial")+
  ggtitle("Neurons- onset time- ROIs per trial- cyto data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

# astrocyte cyto onset histogram
#histogram bins
histseq= seq(0,stimwindow,1)

# counts for each condition in the histogram
Nostim=hist(all.cyto.OT$OnsetTime[(all.cyto.OT$Channel=="GCaMP" & all.cyto.OT$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim=hist(all.cyto.OT$OnsetTime[(all.cyto.OT$Channel=="GCaMP" & all.cyto.OT$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts
Shortstim=hist(all.cyto.OT$OnsetTime[(all.cyto.OT$Channel=="GCaMP" & all.cyto.OT$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
Nostim=Nostim/ntrials.cyto.OT$ntrials[(ntrials.cyto.OT$Condition=="Nostim")]
Stim=Stim/ntrials.cyto.OT$ntrials[(ntrials.cyto.OT$Condition=="Stim")]
Shortstim=Shortstim/ntrials.cyto.OT$ntrials[(ntrials.cyto.OT$Condition=="shortstim")]

#make a data frame for plotting
AC.cyto.histo<- data.frame(cbind(Nostim, Stim, Shortstim))
AC.cyto.histo$time<-histseq[2:length(histseq)]
AC.cyto.histo<-melt(AC.cyto.histo,id="time")

ggplot(data=AC.cyto.histo, aes(x=time, y= value, colour=variable)) + 
  geom_line(size=1.0)+
  xlab("Onset time (s)") +
  ylab("number of ROIs per trial")+
  ggtitle("Astrocytes- onset time- ROIs per trial- cyto data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme



###
# neurons and astrocytes together by density to account for different number of ROIs
ggplot(all.cyto.OT[(all.cyto.OT$Condition!="Nostim"),], aes(x=OnsetTime, y=..density.., colour = Channel, linetype=Condition)) +
  geom_freqpoly(binwidth = 0.5, lwd=1)+
  ggtitle("stim- neurons vs astrocytes-cyto data") + 
  xlab("Onset Time (s)") + 
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

ggplot(all.cyto.OT[(all.cyto.OT$Condition!="Nostim"),], aes(x=OnsetTime, colour = Channel, linetype=Condition)) + 
  stat_ecdf() +
  ggtitle("cdf- neurons vs astrocytes-cyto data- stim types") + 
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

ggplot(all.cyto.OT[(all.cyto.OT$Condition!="Nostim"),], aes(x=OnsetTime, colour = Channel)) + 
  stat_ecdf(lwd=1) +
  ggtitle("cdf- neurons vs astrocytes-cyto data") + 
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme


# zoomed in plot of onset times- identify "FAST" astrocytes
stim.cyto.OT<-subset(all.cyto.OT,Condition!="Nostim")
stim.cyto.OT$Group<-0
stim.cyto.OT$Group[stim.cyto.OT$OnsetTime<1]<-"fast"
stim.cyto.OT$Group[stim.cyto.OT$OnsetTime>=1]<-"delayed"

ggplot(stim.cyto.OT[stim.cyto.OT$Group=="fast",], aes(x=OnsetTime, y=..density.., colour = interaction(Channel,Group), linetype=Condition)) +
  geom_freqpoly(binwidth = 0.0845, lwd=1)+
  ggtitle("stimwindow- neurons vs astrocytes-cyto data") + 
  xlab("Onset Time (s)") + 
  xlim(0,1)+
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

ggplot(stim.cyto.OT[stim.cyto.OT$Group=="fast",], aes(x=OnsetTime, y=..density.., colour = interaction(Channel,Group))) +
  geom_freqpoly(binwidth = 0.0845, lwd=1)+
  ggtitle("stimwindow- neurons vs astrocytes-cyto data") + 
  xlab("Onset Time (s)") + 
  xlim(0,1)+
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

ggplot(stim.cyto.OT[stim.cyto.OT$Group=="fast",], aes(x=OnsetTime, colour = interaction(Channel,Group), linetype=Condition)) + 
  stat_ecdf(lwd=1) +
  ggtitle("cdf- neurons vs astrocytes-fast cyto ROIs") + 
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

# mean onset times
df.OT2<- summarySE(stim.cyto.OT, measurevar = "OnsetTime", groupvars = c("Channel", "Group"))

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
ntrials.cyto.peaks<- ddply(all.lck.peaks, c("Condition"), summarise, ntrials=length(unique(trials)))
ntrials.cyto.DSP4.peaks<- ddply(all.lck.peaks, c("Condition"), summarise, ntrials=length(unique(trials)))


###
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

# mean onset times
df.pT1<- summarySE(all.lck.peaks.group, measurevar = "peakTime", groupvars = c("Channel", "Group","Condition"))


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

#####
# cyto DSP4 peak time data

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
#df.OT3$treatment <- factor(df.OT3$treatment, levels = c("fast","delayed"))

#df.OT1 = df.OT1[!(df.OT1$Channel=="RCaMP"&df.OT1$Group=="delayed"),]

ggplot(df.peaks3, aes(x=Channel,y=peakTime, fill= interaction(Channel,treatment))) +
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

#Results: not significantly different between control and DSP4 for astrocyte peak times!


#########
# CALCULATE PEAK RISE TIME:  Peak time- onset time

#all.lck.peaks, all.lck.OT

# some ROIs may not have detected onset times (if the peak started before the stimulus), these peaks will have a rise time of zero
#all.lck.peaks$OnsetTime=all.lck.peaks$peakTime

#OnsetROInames=unique(all.lck.OT$ROIs_trial)

 #for (iii in 1:length(OnsetROInames))
 #{
  # CurrentOT=all.lck.OT$OnsetTime[all.lck.OT$ROIs_trial==(OnsetROInames[iii])]
  # all.lck.peaks$OnsetTime[grepl(all.lck.peaks$ROIs_trial, OnsetROInames[iii])]=CurrentOT
 #}

#all.lck.peaks$riseTime=all.lck.peaks$peakTime-all.lck.peaks$OnsetTime

# if (nrow(CurrentPeaks)>1)
#  { 
# FirstIdx= min(CurrentPeaks$peakTime)
#  FirstPeak=CurrentPeaks[FirstIdx,]
#}
#  else
# {
#  FirstPeak=CurrentPeaks
#}
#all.lck.firstpeaks<-rbind(all.lck.firstpeaks,FirstPeak)
#}



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


# bar graph
df.dur.lck1<-summarySE(data=all.lck.peaks[all.lck.peaks$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Condition"))

df.dur.lck2<-summarySE(data=all.lck.peaks[all.lck.peaks$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Condition","ROIType"))

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

all.lck.peaks$Condition_type<-interaction(all.lck.peaks$Condition,all.lck.peaks$ROIType)
# stats for astrocyte durations
dur.lck.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],REML=FALSE)
dur.lck.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],REML=FALSE)
dur.lck.model2= lmer(Duration ~ Condition_type + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],REML=FALSE)
dur.lck.anova <- anova(dur.lck.null, dur.lck.model1,dur.lck.model2)
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
dur.lck.group.model2= lmer(Duration ~ Condition_group + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",],REML=FALSE)
dur.lck.group.anova <- anova(dur.lck.group.null, dur.lck.group.model1,dur.lck.group.model2)
print(dur.lck.group.anova)

dur.lck.group.Condition<- glht(dur.lck.group.model1, mcp(Condition= "Tukey"))
summary(dur.lck.group.Condition)

dur.lck.group.Condition_group<- glht(dur.lck.group.model2, mcp(Condition_group= "Tukey"))
summary(dur.lck.group.Condition_group)

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
df.dur.cyto3<- summarySE(control.vs.DSP4.peaks, measurevar = "Duration", groupvars = c("Channel", "Condition","treatment"))

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

######
#  Lck ROI area

lck.ROIarea<-subset(all.lck.peaks, (ROIType=="Dendrite" |ROIType=="Process"))
lck.ROIarea.group<-subset(all.lck.peaks.group, ROIType=="Dendrite" |ROIType=="Process")

lck.ROIarea$Condition<- factor(lck.ROIarea$Condition, levels=c("Nostim","shortstim","Stim"))

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

ggplot(data=df.area3[df.area3$ROIType!="Dendrite",], aes(x=ROIType, y= area, fill=Group)) + 
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


df.cyto.numPeaks1<-summarySE(control.vs.DSP4.trials, measurevar="nPeaks", groupvars=c("Channel","treatment"))
df.cyto.numPeaks2<-summarySE(control.vs.DSP4.trials, measurevar="nPeaks", groupvars=c("Condition","Channel","treatment"))
df.cyto.numPeaks3<-summarySE(control.vs.DSP4.trials.type, measurevar="nPeaks", groupvars=c("ROIType","Condition","Channel","treatment"))


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
cyto.nPeak.anova <- anova(cyto.nPeak.null, cyto.nPeak.model1, cyto.nPeak.model2,cyto.nPeak.model3,cyto.nPeak.model4)
print(cyto.nPeak.anova)

cyto.nPeak.Condition_channel<- glht(cyto.nPeak.model2B, mcp(Condition_channel= "Tukey"))
summary(cyto.nPeak.Condition_channel)

cyto.nPeak.treatment_channel<- glht(cyto.nPeak.model2B, mcp(treatment_channel= "Tukey"))
summary(cyto.nPeak.treatment_channel)

cyto.nPeak.Condition_channel_t<- glht(cyto.nPeak.model3, mcp(Condition_channel_treatment= "Tukey"))
summary(cyto.nPeak.Condition_channel_t)

cyto.nPeak.Condition_channel_t_t<- glht(cyto.nPeak.model4, mcp(Condition_channel_treatment_type= "Tukey"))
summary(cyto.nPeak.Condition_channel_t_t)



# fraction of responding ROIs per trial? or # of peaks per trial (for no stim)



########
#



#Area under the curve for each ROI with an onset time near the neuron?
# only area in second before or second after?
# based on astrocyte onset time



###### 



#####
# short stim

short$ROIType= 0
shortA<- subset(short, Channel=="GCaMP")
shortB<- subset(short, Channel=="RCaMP")

# ROITypes
shortA$ROIType[grepl("r",shortA$ROIname)]="Process"
shortA$ROIType[grepl("E",shortA$ROIname)]="Endfoot"
shortB$ROIType[grepl("r",shortB$ROIname)]="Dendrite"
shortB$ROIType[grepl("D",shortB$ROIname)]="Dendrite"
shortB$ROIType[grepl("N",shortB$ROIname)]="Neuron"

short<-rbind(shortA, shortB)
short$ROIType<- as.factor(short$ROIType)

#unique ROI names
short$ROIs_trial<-paste(short$Animal, short$Spot, short$Trial,short$ROIname, sep= "_")

short$trials<-paste(short$Animal, short$Spot, short$Trial, sep= "_")


# remove matching astrocyte process and soma ROIs
Overlap= short$overlap!=0
short2<-short[!Overlap,]

# stim onset at 5 sec
short2$peakTime<-short2$peakTime-5
short2$peakStart<-short2$peakStart-5
short2$peakStartHalf<-short2$peakStartHalf-5

#duration
short2$Duration<-short2$halfWidth*2

#######
# histogram of peak times for stim
ggplot(long2, aes(x=peakTime, fill=Channel)) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("long stim")

ggplot(short2, aes(x=peakTime, fill=Channel)) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("short stim")

nostim<-subset(stim.all2, Condition=="Nostim")
ggplot(nostim, aes(x=peakTime, fill=Channel)) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("No stim")

RCaMP<- subset(stim.all2, Channel=="RCaMP")
ggplot(RCaMP, aes(x=peakTime, fill=Condition)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("RCaMP")

GCaMP<- subset(stim.all2, Channel=="GCaMP")
ggplot(GCaMP, aes(x=peakTime, fill=Condition)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("GCaMP")


########
# consider only the trials where the neurons responded to stimulation

# what is a neuronal response to stimulation??
# I defined it as a peak within stimulus onset to 1 sec after stimulus stop
# PLUS- peak duration must be close to this window
# long stim= peak between 0 and 9 sec, duration < 11 s
# short stim= peak between 0 and 2 sec, duration < 3 s

# find responding neurons
responding.neurons_long<- subset(longstim, peakTime>0 & peakTime<9 & Duration<11 & Channel=="RCaMP")
responding.neurons_short<- subset(shortstim, peakTime>0 & peakTime<2 & Duration<3 & Channel=="RCaMP")

responding.trials_long<-unique(responding.neurons_long$trials)  # 313 trials of 423- 73.9% of trials
responding.trials_short<-unique(responding.neurons_short$trials) # 96 trials of 149- 64.4% of trials

longstim.responding<-subset(longstim, trials %in% responding.trials_long)
shortstim.responding<-subset(shortstim, trials %in% responding.trials_short)

# distribution of peaktimes
ggplot(longstim.responding, aes(x=peakTime, fill=interaction(Channel,treatment))) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("long stim responding")

ggplot(shortstim.responding, aes(x=peakTime, fill=interaction(Channel,treatment))) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("short stim responding")

#####
#neuronal population sort into high, mid, low groups

# consider all peaks with a time near 10 s
neurons_longstim.mean<- ddply(responding.neurons_long, c("Animal", "Spot", "treatment", "ROIs"), summarise, 
                              PA_mean = mean(peakAUC), nEvents = length(peakAUC),
                              Dur_mean = mean(Duration), Prom_mean = mean(prominence),
                              amp_mean = mean(amplitude), HalfDur = mean(halfWidth),
                              peakT_mean = mean(peakTime),peakHalf_mean= mean(peakStartHalf))

ggplot(neurons_longstim.mean, aes(x=Prom_mean)) + geom_histogram(binwidth=0.05, position="dodge") +
  ggtitle("all neurons during long stim (10 s window)")


Prominence_percentiles<-quantile(neurons_longstim.mean$Prom_mean, prob = seq(0, 1, length = 21), type = 5)

highresponding<-subset(neurons_longstim.mean, Prom_mean>Prominence_percentiles[20])
midresponding<-subset(neurons_longstim.mean, Prom_mean<=Prominence_percentiles[20]&Prom_mean>=Prominence_percentiles[11])
lowresponding<-subset(neurons_longstim.mean, Prom_mean<Prominence_percentiles[11])

longstim.responding$NeuronGroup <- 0
long.highresponders=unique(highresponding$ROIs)
long.midresponders=unique(midresponding$ROIs)
long.lowresponders=unique(lowresponding$ROIs)

longstim.responding$NeuronGroup[longstim.responding$ROIs %in% long.highresponders]<-"high"
longstim.responding$NeuronGroup[longstim.responding$ROIs %in% long.midresponders]<-"mid"
longstim.responding$NeuronGroup[longstim.responding$ROIs %in% long.lowresponders]<-"low"


# consider all peaks with a time near 10 s
neurons_shortstim.mean<- ddply(responding.neurons_short, c("Animal", "Spot", "treatment", "ROIs"), summarise, 
                              PA_mean = mean(peakAUC), nEvents = length(peakAUC),
                              Dur_mean = mean(Duration), Prom_mean = mean(prominence),
                              amp_mean = mean(amplitude), HalfDur = mean(halfWidth),
                              peakT_mean = mean(peakTime),peakHalf_mean= mean(peakStartHalf))

ggplot(neurons_shortstim.mean, aes(x=Prom_mean)) + geom_histogram(binwidth=0.05, position="dodge") +
  ggtitle("all neurons during short stim (10 s window)")


Prominence_percentiles<-quantile(neurons_shortstim.mean$Prom_mean, prob = seq(0, 1, length = 21), type = 5)

highresponding<-subset(neurons_shortstim.mean, Prom_mean>Prominence_percentiles[20])
midresponding<-subset(neurons_shortstim.mean, Prom_mean<=Prominence_percentiles[20]&Prom_mean>=Prominence_percentiles[11])
lowresponding<-subset(neurons_shortstim.mean, Prom_mean<Prominence_percentiles[11])

shortstim.responding$NeuronGroup <- 0
short.highresponders=unique(highresponding$ROIs)
short.midresponders=unique(midresponding$ROIs)
short.lowresponders=unique(lowresponding$ROIs)

shortstim.responding$NeuronGroup[shortstim.responding$ROIs %in% short.highresponders]<-"high"
shortstim.responding$NeuronGroup[shortstim.responding$ROIs %in% short.midresponders]<-"mid"
shortstim.responding$NeuronGroup[shortstim.responding$ROIs %in% short.lowresponders]<-"low"

# are high responders the same during long stim or short stim?
overlaping_high<-intersect(long.highresponders, short.highresponders)

#########
# LONG STIM (90Hz, 8sec)

# identify active ROIs
longstim.responding$ActivePeak <- 0

#responding neurons
farpeaks1 <- longstim.responding$peakTime>0 & longstim.responding$peakTime<9 & longstim.responding$Duration<11 & longstim.responding$Channel=="RCaMP"
longstim.responding$ActivePeak[farpeaks1] <- 1 

#responding astrocytes
farpeaks2 <- longstim.responding$peakTime>0 & longstim.responding$peakTime<12 & longstim.responding$Channel!="RCaMP"
longstim.responding$ActivePeak[farpeaks2] <- 1 

# pull out only the peaks that occur around the stimulation
longstim.stimwindow<- subset(longstim.responding, peakTime>=0 & peakTime<=20)
longstim.active<- subset(longstim.responding, ActivePeak==1)



ggplot(longstim.stimwindow, aes(x=peakTime, fill=Channel)) + geom_histogram(binwidth=0.5, position="dodge") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

ggplot(longstim.active, aes(x=peakTime, fill=Channel)) + geom_histogram(binwidth=0.5, position="dodge") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

longstim.after<- subset(longstim.responding, peakTime>20 & peakTime<80)

ggplot(longstim.after, aes(x=peakTime, fill=interaction(Channel,treatment))) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("long stim after")

#####
library(xlsx)
#write.xlsx(longstim.active, "E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/respondingROIs_longstim.xlsx")
#write.xlsx(longstim.active, "D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/respondingROIs_longstim.xlsx")

#########
# proportion of ROIs that respond per trial
activeROIs1<- ddply(longstim.responding, c("Animal", "Spot", "Channel","trials","ROIType","treatment", "ROIs_trial"), summarise, 
                    nEvents = sum(ActivePeak))

activeROIs1$ROIActive<-0
active1 <- activeROIs1$nEvents>0
activeROIs1$ROIActive[active1] <- 1 

longstim.propActive<- ddply(activeROIs1, c("Animal", "Spot", "Channel","trials","ROIType","treatment"), summarise, 
                            nSignals = sum(nEvents), nROIs= length(unique(ROIs_trial)),
                            nActiveROIs = sum(ROIActive))
longstim.propActive$propActive<-longstim.propActive$nActiveROIs/longstim.propActive$nROIs

df2A1<-summarySE(longstim.propActive, measurevar="propActive", groupvars=c("ROIType","Channel"))
df2A2<-summarySE(longstim.propActive, measurevar="propActive", groupvars=c("Channel","treatment"))

ggplot(longstim.propActive, aes(x=propActive, fill=Channel)) + geom_histogram(binwidth=0.05, position="dodge") +
  ggtitle("long stim prop active ")

df2A1$ROIType <- factor(df2A1$ROIType , levels = c("Endfoot","Process"))
ggplot(data=df2A1, aes(x=ROIType, y=propActive, fill=Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=propActive-se, ymax=propActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("ROIType") +
  ylab("Proportion of Responding ROIs") +
  scale_fill_manual(
    values=c("green", "red")) + 
  max.theme


ggplot(data=df2A2, aes(x=Channel, y=propActive, fill=Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=propActive-se, ymax=propActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Channel") +
  ylab("Proportion of Responding ROIs") +
  ggtitle("Responding ROIs for long stim") 





#########
# SHORT STIM (90Hz, 1sec)
shortstim.responding$ActivePeak <- 0
farpeaks1 <- shortstim.responding$peakTime>0 & shortstim.responding$peakTime<2 & shortstim.responding$Duration<3 & shortstim.responding$Channel=="RCaMP"
shortstim.responding$ActivePeak[farpeaks1] <- 1 

farpeaks2 <- shortstim.responding$peakTime>0 & shortstim.responding$peakTime<5 & shortstim.responding$Channel!="RCaMP"
shortstim.responding$ActivePeak[farpeaks2] <- 1 

# pull out only the peaks that occur around the stimulation
shortstim.stimwindow<- subset(shortstim.responding, peakTime>=0 & peakTime<=10)

ggplot(shortstim.stimwindow, aes(x=peakTime, fill=Channel)) + geom_histogram(binwidth=0.5, position="dodge") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

shortstim.active<- subset(shortstim.responding, ActivePeak==1)

ggplot(shortstim.active, aes(x=peakTime, fill=Channel)) + geom_histogram(binwidth=0.5, position="dodge") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

shortstim.after<- subset(shortstim.responding, peakTime>10 & peakTime<80)

ggplot(shortstim.after, aes(x=peakTime, fill=interaction(Channel,treatment))) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("short stim after")

write.xlsx(shortstim.active, "E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/respondingROIs_shortstim.xlsx")


#####
activeROIs2<- ddply(shortstim.responding, c("Animal", "Spot", "Channel","trials","ROIType","treatment", "ROIs_trial"), summarise, 
                    nEvents = sum(ActivePeak))

activeROIs2$ROIActive<-0
active2 <- activeROIs2$nEvents>0
activeROIs2$ROIActive[active2] <- 1 

shortstim.propActive<- ddply(activeROIs2, c("Animal", "Spot", "Channel","trials","ROIType","treatment"), summarise, 
                             nSignals = sum(nEvents), nROIs= length(unique(ROIs_trial)),
                             nActiveROIs = sum(ROIActive))
shortstim.propActive$propActive<-shortstim.propActive$nActiveROIs/shortstim.propActive$nROIs


df2B1<-summarySE(shortstim.propActive, measurevar="propActive", groupvars=c("ROIType","treatment"))
df2B2<-summarySE(shortstim.propActive, measurevar="propActive", groupvars=c("Channel","treatment"))

ggplot(shortstim.propActive, aes(x=propActive, fill=Channel)) + geom_histogram(binwidth=0.05, position="dodge") +
  ggtitle("short stim prop active ")


df2B1$ROIType <- factor(df2B1$ROIType , levels = c("Neuron","Dendrite","Neuropil","Endfoot","Process"))
ggplot(data=df2B1, aes(x=ROIType, y=propActive, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=propActive-se, ymax=propActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("ROIType") +
  ylab("Proportion of Responding ROIs") +
  ggtitle("Responding ROIs for short stim") 


ggplot(data=df2B2, aes(x=Channel, y=propActive, fill=Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=propActive-se, ymax=propActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Channel") +
  ylab("Proportion of Responding ROIs") +
  ggtitle("Responding ROIs for short stim") 


##########
# considering peaks in stim window
# when do AC peaks come after neuronal peak?
# bin the number of AC peaks per time point after stim?

df9A <- summarySE(longstim.stimwindow, measurevar="peakTime", groupvars=c("ROIType","Channel"))
df9B <- summarySE(shortstim.stimwindow, measurevar="peakTime", groupvars=c("ROIType","treatment"))
df9C <- summarySE(respondingNeurons_Astrocytes, measurevar="peakTime", groupvars=c("ROIType","treatment"))


df10A <- summarySE(longstim.stimwindow, measurevar="peakStart", groupvars=c("ROIType","treatment"))
df10B <- summarySE(shortstim.stimwindow, measurevar="peakStart", groupvars=c("ROIType","treatment"))

df11A <- summarySE(longstim.stimwindow, measurevar="peakStartHalf", groupvars=c("ROIType","Channel"))
df11B <- summarySE(shortstim.stimwindow, measurevar="peakStartHalf", groupvars=c("ROIType","treatment"))

df9A$ROIType <- factor(df9A$ROIType , levels = c("Neuron","Endfoot","Process"))
ggplot(data=df9A, aes(x=ROIType, y=peakTime, fill=Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("Mean Peak Time (s)") +
  scale_fill_manual(
    values=c("green", "red")) + 
  max.theme 

df9B$ROIType <- factor(df9B$ROIType , levels = c("Neuron","Dendrite","Neuropil","Endfoot","Process"))
ggplot(data=df9B, aes(x=ROIType, y=peakTime, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("Mean Peak Time (s)") +
  ggtitle("max peak time for short stim") 

ggplot(data=df10A, aes(x=ROIType, y=peakStart, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStart-se, ymax=peakStart+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("peakStart Time") +
  ggtitle("peakStart times for long stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df10B, aes(x=ROIType, y=peakStart, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStart-se, ymax=peakStart+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("peakStart Time") +
  ggtitle("peakStart times for short stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

df11A$ROIType <- factor(df11A$ROIType , levels = c("Neuron","Endfoot","Process"))
ggplot(data=df11A, aes(x=ROIType, y=peakStartHalf, fill=Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStartHalf-se, ymax=peakStartHalf+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("peakStartHalf Time (t1/2)") +
  ggtitle("peakStartHalf times for long stim") +
  scale_fill_manual(
    values=c("green", "red")) + 
  max.theme 

ggplot(data=df11B, aes(x=ROIType, y=peakStartHalf, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStartHalf-se, ymax=peakStartHalf+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("peakStartHalf Time (t1/2)") +
  ggtitle("peakStartHalf times for short stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

#############
# Likelihood-ratio test 

# long stim, only stim window
ROIType_treatment=interaction(longstim.stimwindow$ROIType,longstim.stimwindow$treatment)
pT.long.null = lmer(peakTime ~ (1|Animal) + (1|Spot), longstim.stimwindow,REML=FALSE)
pT.long.model1 = lmer(peakTime ~ ROIType + (1|Animal) + (1|Spot), longstim.stimwindow,REML=FALSE)
pT.long.model2 = lmer(peakTime ~ treatment + (1|Animal) + (1|Spot), longstim.stimwindow,REML=FALSE)
pT.long.model3 = lmer(peakTime ~ ROIType+treatment + (1|Animal) + (1|Spot), longstim.stimwindow,REML=FALSE)
pT.long.model4 = lmer(peakTime ~ ROIType_treatment + (1|Animal) + (1|Spot), longstim.stimwindow,REML=FALSE)

pT.long.anova <- anova(pT.long.null, pT.long.model1,pT.long.model2,pT.long.model3,pT.long.model4)
print(pT.long.anova)

pT.long.anova <- anova(pT.long.null, pT.long.model1)
print(pT.long.anova)

# p values
pT.pv.longstim <- glht(pT.long.model4, mcp(ROIType_treatment= "Tukey"))
summary(pT.pv.longstim)

pT.pv.longstim <- glht(pT.long.model1, mcp(ROIType= "Tukey"))
summary(pT.pv.longstim)

# short stim, only stim window
ROIType_treatment=interaction(shortstim.stimwindow$ROIType,shortstim.stimwindow$treatment)
pT.short.null = lmer(peakTime ~ (1|Animal) + (1|Spot), shortstim.stimwindow,REML=FALSE)
pT.short.model1 = lmer(peakTime ~ ROIType + (1|Animal) + (1|Spot), shortstim.stimwindow,REML=FALSE)
pT.short.model2 = lmer(peakTime ~ treatment + (1|Animal) + (1|Spot), shortstim.stimwindow,REML=FALSE)
pT.short.model3 = lmer(peakTime ~ ROIType+treatment + (1|Animal) + (1|Spot), shortstim.stimwindow,REML=FALSE)
pT.short.model4 = lmer(peakTime ~ ROIType_treatment + (1|Animal) + (1|Spot), shortstim.stimwindow,REML=FALSE)

pT.short.anova <- anova(pT.short.null, pT.short.model1,pT.short.model2,pT.short.model3,pT.short.model4)
print(pT.short.anova)

pT.short.anova <- anova(pT.short.null, pT.short.model1)
print(pT.short.anova)


# p values
pT.pv.shortstim <- glht(pT.short.model4, mcp(ROIType_treatment= "Tukey"))
summary(pT.pv.shortstim)

pT.pv.shortstim <- glht(pT.short.model1, mcp(ROIType= "Tukey"))
summary(pT.pv.shortstim)

# half maximum time
ROIType_treatment=interaction(longstim.stimwindow$ROIType,longstim.stimwindow$treatment)
HM.null = lmer(peakStartHalf ~ (1|Animal) + (1|Spot), longstim.stimwindow,REML=FALSE)
HM.model1 = lmer(peakStartHalf ~ ROIType+ (1|Animal) + (1|Spot), longstim.stimwindow,REML=FALSE)
HM.anova <- anova(HM.null, HM.model1)
print(HM.anova)
# p values
HM.pv.ROIType <- glht(HM.model1, mcp(ROIType= "Tukey"))
summary(HM.pv.ROIType)

############
# consider peaks that are close together in time (determined in Matlab)
#early<-read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/earlyGC_byAUC.csv", header=TRUE, sep = ",")
early<-read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/earlyGC_byAUC.csv", header=TRUE, sep = ",")


#proportion of astrocyte ROIs that are early/late in each treatment group, in each trial
longstim.gcamp=subset(longstim.stimwindow, Channel=="GCaMP")
longstim.gcamp$TimeGroup<-0
longstim.gcamp$TimeGroup[longstim.gcamp$ROIs_trial %in% early$names]<-"early"
longstim.gcamp$TimeGroup[!(longstim.gcamp$ROIs_trial %in% early$names)]<-"late"

longstim.rcamp=subset(longstim.stimwindow, Channel=="RCaMP")
longstim.rcamp$TimeGroup<-"neurons"
longstim.both<-rbind(longstim.gcamp, longstim.rcamp)

# Simple Pie Chart
Enum=length(early$names)
Lnum=length(unique(longstim.gcamp$ROIs_trial))-Enum
slices <- c(Enum, Lnum)
lbls <- c("Fast", "Slow")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=rainbow(length(lbls)))


#what are the mean characteristics of early astrocyte ROIs
#amplitudes, mean peak times, etc.

# amplitude
df3A <- summarySE(longstim.gcamp, measurevar="amplitude", groupvars=c("TimeGroup"),na.rm=TRUE)
df3B <- summarySE(longstim.gcamp, measurevar="amplitude", groupvars=c("TimeGroup","ROIType"),na.rm=TRUE)


ggplot(data=df3A, aes(x=TimeGroup, y=amplitude, fill=TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("TimeGroup") +
  ylab("amplitude") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

ggplot(data=df3B, aes(x=ROIType, y=amplitude, fill=TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("TimeGroup") +
  ylab("amplitude") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

#stats
amp.null = lmer(amplitude ~ (1|Animal) + (1|Spot), longstim.gcamp,REML=FALSE)
amp.model1 = lmer(amplitude ~ TimeGroup+ (1|Animal) + (1|Spot), longstim.gcamp,REML=FALSE)
amp.anova <- anova(amp.null, amp.model1)
print(amp.anova)
# p values
amp.pv.ROIType <- glht(amp.model1, mcp(TimeGroup= "Tukey"))
summary(amp.pv.ROIType)

# peakTimes
df4A <- summarySE(longstim.gcamp, measurevar="peakTime", groupvars=c("TimeGroup"),na.rm=TRUE)
df4B <- summarySE(longstim.gcamp, measurevar="peakStartHalf", groupvars=c("TimeGroup"),na.rm=TRUE)

df4C <- summarySE(longstim.both, measurevar="peakTime", groupvars=c("TimeGroup"),na.rm=TRUE)
df4D <- summarySE(longstim.both, measurevar="peakStartHalf", groupvars=c("TimeGroup"),na.rm=TRUE)
df4E <- summarySE(longstim.both, measurevar="peakStart", groupvars=c("TimeGroup"),na.rm=TRUE)

df4C$TimeGroup <- factor(df4C$TimeGroup, levels = c("neurons","early","late"))
df4D$TimeGroup <- factor(df4D$TimeGroup, levels = c("neurons","early","late"))
df4E$TimeGroup <- factor(df4E$TimeGroup, levels = c("neurons","early","late"))

ggplot(data=df4A, aes(x=TimeGroup, y=peakTime, fill=TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("TimeGroup") +
  ylab("peakTime") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

ggplot(data=df4B, aes(x=TimeGroup, y=peakStartHalf, fill=TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStartHalf-se, ymax=peakStartHalf+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("TimeGroup") +
  ylab("peakStartHalf") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

ggplot(data=df4C, aes(x=TimeGroup, y=peakTime, fill=TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("TimeGroup") +
  ylab("peakTime") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

ggplot(data=df4D, aes(x=TimeGroup, y=peakStartHalf, fill=TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStartHalf-se, ymax=peakStartHalf+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("TimeGroup") +
  ylab("peakStartHalf") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

ggplot(data=df4E, aes(x=TimeGroup, y=peakStart, fill=TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStart-se, ymax=peakStart+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("TimeGroup") +
  ylab("peakStart") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme


#stats
pSH.null = lmer(peakStartHalf ~ (1|Animal) + (1|Spot), longstim.both,REML=FALSE)
pSH.model1 = lmer(peakStartHalf ~ TimeGroup+ (1|Animal) + (1|Spot), longstim.both,REML=FALSE)
pSH.anova <- anova(pSH.null, pSH.model1)
print(pSH.anova)
# p values
pSH.pv.ROIType <- glht(pSH.model1, mcp(TimeGroup= "Tukey"))
summary(pSH.pv.ROIType)

#stats
pS.null = lmer(peakStart ~ (1|Animal) + (1|Spot), longstim.both,REML=FALSE)
pS.model1 = lmer(peakStart ~ TimeGroup+ (1|Animal) + (1|Spot), longstim.both,REML=FALSE)
pS.anova <- anova(pS.null, pS.model1)
print(pS.anova)
# p values
pS.pv.ROIType <- glht(pS.model1, mcp(TimeGroup= "Tukey"))
summary(pS.pv.ROIType)


df5A <- summarySE(longstim.gcamp, measurevar="Duration", groupvars=c("TimeGroup"),na.rm=TRUE)
ggplot(data=df5A, aes(x=TimeGroup, y=Duration, fill=TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("TimeGroup") +
  ylab("Duration") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

#stats
dur.null = lmer(Duration ~ (1|Animal) + (1|Spot), longstim.gcamp,REML=FALSE)
dur.model1 = lmer(Duration ~ TimeGroup+ (1|Animal) + (1|Spot), longstim.gcamp,REML=FALSE)
dur.anova <- anova(dur.null, dur.model1)
print(dur.anova)
# p values
dur.pv.ROIType <- glht(dur.model1, mcp(TimeGroup= "Tukey"))
summary(dur.pv.ROIType)


#what are the correlations of the ROIs with these time differences?
#spatially similar?

#####
# early short stim
early<-read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/earlyGC_byAUCshort.csv", header=TRUE, sep = ",")

#proportion of astrocyte ROIs that are early/late in each treatment group, in each trial
shortstim.gcamp=subset(shortstim.stimwindow, Channel=="GCaMP")
shortstim.gcamp$TimeGroup<-0
shortstim.gcamp$TimeGroup[shortstim.gcamp$ROIs_trial %in% early$names]<-"early"
shortstim.gcamp$TimeGroup[!(shortstim.gcamp$ROIs_trial %in% early$names)]<-"late"

shortstim.rcamp=subset(shortstim.stimwindow, Channel=="RCaMP")
shortstim.rcamp$TimeGroup<-"neurons"
shortstim.both<-rbind(shortstim.gcamp, shortstim.rcamp)

# Simple Pie Chart
Enum=length(early$names)
Lnum=length(unique(shortstim.gcamp$ROIs_trial))-Enum
slices <- c(Enum, Lnum)
lbls <- c("Fast", "Slow")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=rainbow(length(lbls)))


#what are the mean characteristics of early astrocyte ROIs
#amplitudes, mean peak times, etc.

# amplitude
df3A <- summarySE(shortstim.gcamp, measurevar="amplitude", groupvars=c("TimeGroup"),na.rm=TRUE)
df3B <- summarySE(shortstim.gcamp, measurevar="amplitude", groupvars=c("TimeGroup","ROIType"),na.rm=TRUE)


ggplot(data=df3A, aes(x=TimeGroup, y=amplitude, fill=TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("TimeGroup") +
  ylab("amplitude") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

ggplot(data=df3B, aes(x=ROIType, y=amplitude, fill=TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("TimeGroup") +
  ylab("amplitude") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

#stats
amp.null = lmer(amplitude ~ (1|Animal) + (1|Spot), shortstim.gcamp,REML=FALSE)
amp.model1 = lmer(amplitude ~ TimeGroup+ (1|Animal) + (1|Spot), shortstim.gcamp,REML=FALSE)
amp.anova <- anova(amp.null, amp.model1)
print(amp.anova)
# p values
amp.pv.ROIType <- glht(amp.model1, mcp(TimeGroup= "Tukey"))
summary(amp.pv.ROIType)

# peakTimes
df4A <- summarySE(shortstim.gcamp, measurevar="peakTime", groupvars=c("TimeGroup"),na.rm=TRUE)
df4B <- summarySE(shortstim.gcamp, measurevar="peakStartHalf", groupvars=c("TimeGroup"),na.rm=TRUE)

df4C <- summarySE(shortstim.both, measurevar="peakTime", groupvars=c("TimeGroup"),na.rm=TRUE)
df4D <- summarySE(shortstim.both, measurevar="peakStartHalf", groupvars=c("TimeGroup"),na.rm=TRUE)
df4E <- summarySE(shortstim.both, measurevar="peakStart", groupvars=c("TimeGroup"),na.rm=TRUE)

df4C$TimeGroup <- factor(df4C$TimeGroup, levels = c("neurons","early","late"))
df4D$TimeGroup <- factor(df4D$TimeGroup, levels = c("neurons","early","late"))
df4E$TimeGroup <- factor(df4E$TimeGroup, levels = c("neurons","early","late"))

ggplot(data=df4A, aes(x=TimeGroup, y=peakTime, fill=TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("TimeGroup") +
  ylab("peakTime") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

ggplot(data=df4B, aes(x=TimeGroup, y=peakStartHalf, fill=TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStartHalf-se, ymax=peakStartHalf+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("TimeGroup") +
  ylab("peakStartHalf") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

ggplot(data=df4C, aes(x=TimeGroup, y=peakTime, fill=TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("TimeGroup") +
  ylab("peakTime") +
  scale_fill_manual(
    values=c("red", "blue","green")) + 
  max.theme

ggplot(data=df4D, aes(x=TimeGroup, y=peakStartHalf, fill=TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStartHalf-se, ymax=peakStartHalf+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("TimeGroup") +
  ylab("peakStartHalf") +
  scale_fill_manual(
    values=c("red", "blue","green")) + 
  max.theme

ggplot(data=df4E, aes(x=TimeGroup, y=peakStart, fill=TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStart-se, ymax=peakStart+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("TimeGroup") +
  ylab("peakStart") +
  scale_fill_manual(
    values=c("red", "blue","green")) + 
  max.theme


#stats
pSH.null = lmer(peakStartHalf ~ (1|Animal) + (1|Spot), shortstim.both,REML=FALSE)
pSH.model1 = lmer(peakStartHalf ~ TimeGroup+ (1|Animal) + (1|Spot), shortstim.both,REML=FALSE)
pSH.anova <- anova(pSH.null, pSH.model1)
print(pSH.anova)
# p values
pSH.pv.ROIType <- glht(pSH.model1, mcp(TimeGroup= "Tukey"))
summary(pSH.pv.ROIType)

#stats
pS.null = lmer(peakStart ~ (1|Animal) + (1|Spot), shortstim.both,REML=FALSE)
pS.model1 = lmer(peakStart ~ TimeGroup+ (1|Animal) + (1|Spot), shortstim.both,REML=FALSE)
pS.anova <- anova(pS.null, pS.model1)
print(pS.anova)
# p values
pS.pv.ROIType <- glht(pS.model1, mcp(TimeGroup= "Tukey"))
summary(pS.pv.ROIType)


df5A <- summarySE(shortstim.gcamp, measurevar="Duration", groupvars=c("TimeGroup"),na.rm=TRUE)
ggplot(data=df5A, aes(x=TimeGroup, y=Duration, fill=TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("TimeGroup") +
  ylab("Duration") +
  scale_fill_manual(
    values=cbbPalette) + 
  max.theme

#stats
dur.null = lmer(Duration ~ (1|Animal) + (1|Spot), shortstim.gcamp,REML=FALSE)
dur.model1 = lmer(Duration ~ TimeGroup+ (1|Animal) + (1|Spot), shortstim.gcamp,REML=FALSE)
dur.anova <- anova(dur.null, dur.model1)
print(dur.anova)
# p values
dur.pv.ROIType <- glht(dur.model1, mcp(TimeGroup= "Tukey"))
summary(dur.pv.ROIType)



#######
# no stim peaks


