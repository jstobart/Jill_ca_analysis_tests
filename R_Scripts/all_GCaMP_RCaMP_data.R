
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
# data sets from nostim, shortstim (90Hz,1s), longstim (90Hz,8s)
# with hand selected somata and FLIKA activities for neurons AND astrocytes

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
all.lck.peaks<-all.lck.peaks[!(all.lck.peaks$ROIname=="none"),]

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


# neurons and astrocytes together by density to account for different number of ROIs
<<<<<<< HEAD


ggplot(all.lck.OT[(all.lck.OT$Condition!="Nostim"),], aes(x=OnsetTime, y=..density.., colour = Channel, linetype=Condition)) +
  geom_freqpoly(binwidth = 0.5, lwd=1)+
  ggtitle("stim- neurons vs astrocytes-lck data") + 
  xlab("Onset Time (s)") + 
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme
=======
histseq2= seq(0,12, 0.5)

GC.stim<-hist(all.lck.OT$OnsetTime[(all.lck.OT$Channel=="GCaMP" & all.lck.OT$Condition=="Stim")], breaks=histseq2, plot=FALSE)$density
RC.stim<-hist(all.lck.OT$OnsetTime[(all.lck.OT$Channel=="RCaMP" & all.lck.OT$Condition=="Stim")], breaks=histseq2, plot=FALSE)$density
>>>>>>> origin/master

lck.longstim.histo<-data.frame(cbind(GC.stim, RC.stim))
lck.longstim.histo$time<-histseq2[2:length(histseq2)]
lck.longstim.histo2<-melt(lck.longstim.histo,id="time")

<<<<<<< HEAD

# zoomed in plot of onset times
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

# mean onset times
df.OT1<- summarySE(stim.lck.OT, measurevar = "OnsetTime", groupvars = c("Channel", "Group"))

#df.OT1$OnsetTime<-df.OT1$OnsetTime*1000  # convert to ms
#df.OT1$se<-df.OT1$se*1000 # convert to ms

df.OT1$Channel <- factor(df.OT1$Channel, levels = c("RCaMP","GCaMP"))
df.OT1$Condition <- factor(df.OT1$Condition, levels = c("shortstim","Stim"))
df.OT1$Group <- factor(df.OT1$Group, levels = c("fast","delayed"))

df.OT1 = df.OT1[!(df.OT1$Channel=="RCaMP"&df.OT1$Group=="delayed"),]

ggplot(df.OT1, aes(x=interaction(Channel,Group),y=OnsetTime, fill= interaction(Channel,Group))) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(data=df2A1, aes(x=ROIType, y=propActive, fill=Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=propActive-se, ymax=propActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("ROIType") +
  ylab("Proportion of Responding ROIs") +
  scale_fill_manual(
    values=c("green", "red")) + 
  max.theme

stim.lck.OT$Group<-as.factor(stim.lck.OT$Group)
Condition_Channel=interaction(stim.lck.OT$Condition,stim.lck.OT$Channel)
Group_Channel=interaction(stim.lck.OT$Group,stim.lck.OT$Channel)
# stats for onset times- neurons vs astrocytes
OT.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial) + (1|ROIs_trial), stim.lck.OT,REML=FALSE)
OT.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Spot_trial) + (1|ROIs_trial), stim.lck.OT,REML=FALSE)
OT.model2 = lmer(OnsetTime ~ Condition + (1|Animal) + (1|Spot) + (1|Spot_trial) + (1|ROIs_trial), stim.lck.OT,REML=FALSE)
OT.model3 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|Spot_trial) + (1|ROIs_trial), stim.lck.OT,REML=FALSE)
OT.model4 = lmer(OnsetTime ~ Group_Channel + (1|Animal) + (1|Spot) + (1|Spot_trial) + (1|ROIs_trial), stim.lck.OT,REML=FALSE)
OT.model5 = lmer(OnsetTime ~ Condition_Channel + (1|Animal) + (1|Spot) + (1|Spot_trial) + (1|ROIs_trial), stim.lck.OT,REML=FALSE)
OT.anova <- anova(OT.null, OT.model1,OT.model2,OT.model3,OT.model4,OT.model5)
print(OT.anova)

OT.Group_channel<- glht(OT.model4, mcp(Group_Channel= "Tukey"))
summary(OT.Group_channel)

=======
ggplot(data=lck.longstim.histo2, aes(x=time, y= value, colour=variable)) + 
  geom_line(size=1)+
  xlab("Onset time (s)") +
  ylab("Density")+
  ggtitle("Astrocytes vs neurons- long onset time- lck data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme


ggplot(all.lck.OT[(all.lck.OT$Condition=="Stim"),], aes(x=OnsetTime, fill=Channel)) + 
  geom_density(aes(y=..scaled.., Colour=Channel, fill=Channel), alpha=0.7, adjust=1/5,size=1) +
  scale_fill_manual(values=cbbPalette[2:3])+
  max.theme

library("scales")
ggplot(all.lck.OT[(all.lck.OT$Condition=="Stim"),], aes(x=OnsetTime, fill=Channel)) + 
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  ## scale_y_continuous(labels = percent_format()) #version 3.0.9
  #scale_y_continuous(labels = percent_format())+
  ggtitle("astrocytes vs neurons- onset time- long stim")+ 
  max.theme  
>>>>>>> origin/master

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

# astrocyte lck onset histogram
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

# neurons vs. astrocytes
ggplot(all.cyto.OT[(all.cyto.OT$Condition!="Nostim"),], aes(x=OnsetTime, y=..density.., colour = Channel, linetype=Condition)) +
  geom_freqpoly(binwidth = 0.5, lwd=1)+
  ggtitle("stim- neurons vs astrocytes-cyto data") + 
  xlab("Onset Time (s)") + 
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme



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

ggplot(all.cyto.DSP4.OT[(all.cyto.DSP4.OT$Condition!="Nostim"),], aes(x=OnsetTime, y=..density.., colour = Channel, linetype=Condition)) +
  geom_freqpoly(binwidth = 0.5, lwd=1)+
  ggtitle("stim- neurons vs astrocytes-cyto.DSP4 data") + 
  xlab("Onset Time (s)") + 
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
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

######

# consider ONLY the first peak for each ROI after the start of stimulation
#ROInames=unique(all.lck.peaks$ROIs_trial)
#all.lck.firstpeaks=data.frame()
#for (ii in 1:length(ROInames))
#{
 # CurrentROI=ROInames[ii]
#  CurrentPeaks=subset(all.lck.peaks, ROIs_trial %in% CurrentROI)
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


######
# Duration histograms

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

<<<<<<< HEAD
# boxplot

all.lck.peaks$Condition <- factor(all.lck.peaks$Condition, levels = c("Nostim","shortstim","Stim"))

=======
levels(all.lck.peaks$Condition) <- c("Nostim","shortstim","Stim")
>>>>>>> origin/master
ggplot(data=all.lck.peaks[all.lck.peaks$Channel=="GCaMP",], aes(x=Condition, y= Duration, colour=Condition)) + 
  geom_boxplot(size=1)+
  ylab("Duration (s)")+
  ggtitle("Astrocytes- duration - lck data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

df.dur<-summarySE(data=all.lck.peaks[all.lck.peaks$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Condition"))

ggplot(data=df.dur, aes(x=Condition, y= Duration, colour=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
    max.theme

all.lck.peaks$Condition_type<-interaction(all.lck.peaks$Condition,all.lck.peaks$ROIType)
# stats for astrocyte durations
dur.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],REML=FALSE)
dur.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],REML=FALSE)
dur.model2= lmer(Duration ~ Condition_type + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],REML=FALSE)
dur.anova <- anova(dur.null, dur.model1,dur.model2)
print(dur.anova)

dur.Condition<- glht(dur.model1, mcp(Condition= "Tukey"))
summary(dur.Condition)

dur.Condition_type<- glht(dur.model2, mcp(Condition_type= "Tukey"))
summary(dur.Condition_type)

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

<<<<<<< HEAD
all.cyto.peaks$Condition <- factor(all.cyto.peaks$Condition, levels = c("Nostim","shortstim","Stim"))
=======
>>>>>>> origin/master
ggplot(data=all.cyto.peaks[all.cyto.peaks$Channel=="GCaMP",], aes(x=Condition, y= Duration, colour=Condition)) + 
  geom_boxplot(size=1)+
  ylab("Duration (s)")+
  ggtitle("Astrocytes- duration - cyto data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

# boxplot

df.dur2<-summarySE(data=all.cyto.peaks[all.cyto.peaks$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Condition"))

ggplot(data=df.dur2, aes(x=Condition, y= Duration, colour=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

all.cyto.peaks$Condition_type<-interaction(all.cyto.peaks$Condition,all.cyto.peaks$ROIType)
# stats for astrocyte durations
dur.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.cyto.peaks[all.cyto.peaks$Channel=="GCaMP",],REML=FALSE)
dur.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.cyto.peaks[all.cyto.peaks$Channel=="GCaMP",],REML=FALSE)
dur.model2= lmer(Duration ~ Condition_type + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.cyto.peaks[all.cyto.peaks$Channel=="GCaMP",],REML=FALSE)
dur.anova <- anova(dur.null, dur.model1,dur.model2)
print(dur.anova)

dur.Condition<- glht(dur.model1, mcp(Condition= "Tukey"))
summary(dur.Condition)

dur.Condition_type<- glht(dur.model2, mcp(Condition_type= "Tukey"))
summary(dur.Condition_type)

######
# astrocyte cyto DSP4 duration

#histogram bins
histseq= seq(0,80, 1)

# counts for each condition in the histogram
Nostim=hist(all.cyto.DSP4.peaks$Duration[(all.cyto.DSP4.peaks$Channel=="GCaMP" & all.cyto.DSP4.peaks$Condition=="Nostim")], breaks=histseq, plot=FALSE)$counts
Stim=hist(all.cyto.DSP4.peaks$Duration[(all.cyto.DSP4.peaks$Channel=="GCaMP" & all.cyto.DSP4.peaks$Condition=="Stim")], breaks=histseq, plot=FALSE)$counts
Shortstim=hist(all.cyto.DSP4.peaks$Duration[(all.cyto.DSP4.peaks$Channel=="GCaMP" & all.cyto.DSP4.peaks$Condition=="shortstim")], breaks=histseq, plot=FALSE)$counts

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
  xlab("Duration (s)") +
  ylab("number of peaks per trial")+
  ggtitle("Astrocytes- duration- ROIs per trial- cyto DSP4 data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

<<<<<<< HEAD
all.cyto.DSP4.peaks$Condition <- factor(all.cyto.DSP4.peaks$Condition, levels = c("Nostim","shortstim","Stim"))
=======
>>>>>>> origin/master
ggplot(data=all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="GCaMP",], aes(x=Condition, y= Duration, colour=Condition)) + 
  geom_boxplot(size=1)+
  ylab("Duration (s)")+
  ggtitle("Astrocytes- duration - cyto data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

# boxplot

df.dur2<-summarySE(data=all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Condition"))

ggplot(data=df.dur2, aes(x=Condition, y= Duration, colour=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

all.cyto.DSP4.peaks$Condition_type<-interaction(all.cyto.DSP4.peaks$Condition,all.cyto.DSP4.peaks$ROIType)
# stats for astrocyte durations
dur.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="GCaMP",],REML=FALSE)
dur.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="GCaMP",],REML=FALSE)
dur.model2= lmer(Duration ~ Condition_type + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.cyto.DSP4.peaks[all.cyto.DSP4.peaks$Channel=="GCaMP",],REML=FALSE)
dur.anova <- anova(dur.null, dur.model1,dur.model2)
print(dur.anova)

dur.Condition<- glht(dur.model1, mcp(Condition= "Tukey"))
summary(dur.Condition)

dur.Condition_type<- glht(dur.model2, mcp(Condition_type= "Tukey"))
summary(dur.Condition_type)

######

# amplitude histograms

# lck data

ggplot(data=all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],x=amplitude) + 
  geom_histogram(binwidth=0.1)
  
# astrocyte lck duration histogram
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

<<<<<<< HEAD
ggplot(all.lck.peaks[(all.lck.peaks$Channel=="GCaMP"&all.lck.peaks$Condition!="Nostim"),], aes(x=amplitude, y=..density.., colour = ROIType)) +
  geom_freqpoly(binwidth = 0.1, lwd=1)+
  ggtitle("stim- astrocytes-lck") + 
  xlab("amplitude") + 
  scale_colour_manual(values=c(cbbPalette[3],cbbPalette[2]))+
  max.theme

=======
levels(all.lck.peaks$Condition) <- c("Nostim","shortstim","Stim")
>>>>>>> origin/master
ggplot(data=all.lck.peaks[all.lck.peaks$Channel=="GCaMP",], aes(x=Condition, y= amplitude, colour=Condition)) + 
  geom_boxplot(size=1)+
  ylab("amplitude dF/F")+
  ggtitle("Astrocytes- amplitude- - lck data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme

<<<<<<< HEAD
df.amp<-summarySE(data=all.lck.peaks[all.lck.peaks$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("Condition"))
=======
######
# ROI area
levels(all.lck.peaks$Condition) <- c("Nostim","shortstim","Stim")
ggplot(data=all.lck.peaks, aes(x=Condition, y= area, colour=Condition)) + 
  geom_boxplot(size=1)+
  ylab("area (sq m)")+
  ggtitle("Astrocytes- roi area- - lck data") + 
  scale_colour_manual(values=cbbPalette)+
  max.theme


df1A<- summarySE(all.lck.peaks, measurevar="area", groupvars=c("Condition","Channel"))
df1B<- summarySE(all.lck.peaks, measurevar="area", groupvars=c("Condition","Channel","ROIType"))
>>>>>>> origin/master

ggplot(data=df.amp, aes(x=Condition, y= amplitude, colour=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

all.lck.peaks$Condition_type<-interaction(all.lck.peaks$Condition,all.lck.peaks$ROIType)
# stats for astrocyte durations
amplitude.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],REML=FALSE)
amplitude.model1 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],REML=FALSE)
amplitude.model2= lmer(amplitude~ Condition_type + (1|Animal) + (1|Spot) + (1|trials) + (1|ROIs_trial), all.lck.peaks[all.lck.peaks$Channel=="GCaMP",],REML=FALSE)
amplitude.anova <- anova(amplitude.null, amplitude.model1,amplitude.model2)
print(amplitude.anova)

amplitude.Condition<- glht(amplitude.model1, mcp(Condition= "Tukey"))
summary(amplitude.Condition)

amplitude.Condition_type<- glht(amplitude.model2, mcp(Condition_type= "Tukey"))
summary(amplitude.Condition_type)

######





########
# onset time comparisons- neurons vs. astrocytes
longstim.OT.comp$compType<-paste(longstim.OT$N_ROIType, longstim.OT$A_ROIType, sep= "_")
shortstim.OT.comp$compType<-paste(shortstim.OT$N_ROIType, shortstim.OT$A_ROIType, sep= "_")
nostim.OT.comp$compType<-paste(nostim.OT$N_ROIType, nostim.OT$A_ROIType, sep= "_")


# get rid of onsets from the last few seconds of the trial (probably not a complete peak and we can't measure it)

nostim.OT=nostim.OT[nostim.OT$N_Onset<43,]
nostim.OT=nostim.OT[nostim.OT$A_Onset<43,]

# subset data to 3 sec on either side (so AC peaks 3 sec before or after neuronal peaks)
nostim.OT.small<-subset(nostim.OT, TimeDiff<3 & TimeDiff>-3)

nostim.OT.close<-subset(nostim.OT.small, distance<5)


######
# histograms of neuron-astrocyte time differences
ggplot(nostim.OT, aes(x=TimeDiff)) + geom_histogram(binwidth=0.0845, position="dodge") +
  ggtitle("onset time differences- all peaks") +
  max.theme

ggplot(nostim.OT, aes(x=TimeDiff, fill=A_ROIType)) + geom_histogram(binwidth=0.0845, position="dodge") +
  ggtitle("onset time differences- all peaks") +
  max.theme

# histograms of time differences
ggplot(nostim.OT.small, aes(x=TimeDiff)) + geom_histogram(binwidth=0.0845, position="dodge") +
  ggtitle("onset time differences- peaks close to zero time difference") + 
  max.theme

library('scales')
ggplot(nostim.OT.small, aes(x=TimeDiff, fill=A_ROIType)) + 
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  ## scale_y_continuous(labels = percent_format()) #version 3.0.9
  scale_y_continuous(labels = percent_format())+
  ggtitle("percent distribution- onset time differences")+ 
  max.theme  
  

# density of time differences
ggplot(nostim.OT.small, aes(x=TimeDiff)) + geom_density(aes(group=A_ROIType, colour=A_ROIType, fill=A_ROIType), alpha=0.3, adjust=1/3,size=1) +
  max.theme

ggplot(nostim.OT, aes(x=TimeDiff)) + geom_density(aes(group=A_ROIType, colour=A_ROIType, fill=A_ROIType), alpha=0.3, adjust=1/5,size=1) +
  max.theme

ggplot(nostim.OT.small, aes(x=TimeDiff)) + geom_density(aes(group=compType, colour=compType, fill=compType), alpha=0.3, adjust=1/5,size=1) +
  max.theme

ggplot(nostim.OT.small, aes(x=TimeDiff)) + geom_histogram(binwidth = 0.0845,aes( y=..density..,fill=A_ROIType)) +
  #geom_density(aes(group=A_ROIType), size=1) +
  ggtitle("onset time differences- peaks close to zero time difference") + 
  max.theme

ggplot(nostim.OT.small, aes(x=TimeDiff, fill=A_ROIType)) + 
  geom_histogram(aes(y = ..density..),binwidth=0.0845, alpha = 0.7) + 
  geom_density(aes(fill = A_ROIType), alpha = 0.5) 


ggplot(nostim.OT.small, aes(x=TimeDiff, y=..density..,colour=A_ROIType)) + stat_bin(geom="step", binwidth=(0.0845*2))+max.theme


#aes(y=..count../sum(..count..))

# histograms of ROI distances
ggplot(nostim.OT.small, aes(x=distance)) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("distance between ROIs- close peaks") + max.theme

ggplot(nostim.OT.small, aes(x=distance, fill=A_ROIType)) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("distance between ROIs- close peaks") + max.theme


#Distance between ROIs that are compared

ggplot(nostim.OT.small, aes(x=distance, y=TimeDiff, colour=compType)) +
  geom_point()+ max.theme

ggplot(nostim.OT.small, aes(x=distance, y=TimeDiff)) +
  geom_point(alpha=1/20)+ max.theme

ggplot(nostim.OT.small, aes(x=distance, y=TimeDiff)) +
  geom_count()+ max.theme

ggplot(nostim.OT.small, aes(x=distance, y=TimeDiff)) +
  geom_hex(bins=20)+ max.theme


# mean time diffs and mean distances
df1A1<-summarySE(nostim.OT.small, measurevar="TimeDiff", groupvars=c("compType"))
df1A2<-summarySE(nostim.OT.small, measurevar="TimeDiff", groupvars=c("A_ROIType"))

df2A1<-summarySE(nostim.OT.small, measurevar="distance", groupvars=c("compType"))
df2A2<-summarySE(nostim.OT.small, measurevar="distance", groupvars=c("A_ROIType"))

######
Before<- subset(nostim.OT.small, TimeDiff<=0)
After <- subset(nostim.OT.small, TimeDiff>=0)
Zero <- subset(nostim.OT.small, TimeDiff==0)


df1A3<-summarySE(Before, measurevar="TimeDiff", groupvars=c("A_ROIType"))
df1A4<-summarySE(After, measurevar="TimeDiff", groupvars=c("A_ROIType"))

df2A3<-summarySE(Before, measurevar="distance", groupvars=c("A_ROIType"))
df2A4<-summarySE(After, measurevar="distance", groupvars=c("A_ROIType"))

ggplot(After, aes(x=TimeDiff)) + geom_histogram(aes(group=A_ROIType, colour=A_ROIType, fill=A_ROIType), binwidth=0.0845)+
  max.theme

ggplot(After, aes(x=TimeDiff)) + geom_density(aes(group=A_ROIType, colour=A_ROIType, fill=A_ROIType), alpha=0.3, adjust=1/2,size=1) +
  max.theme

# stats for after group
# endfeet delayed compared to processes
OT.null = lmer(TimeDiff ~ (1|Animal) + (1|Spot) + (1|trials) + (1|N_ROI) + (1|A_ROI), After,REML=FALSE)
OT.model1 = lmer(TimeDiff ~ A_ROIType + (1|Animal) + (1|Spot) + (1|trials) +(1|N_ROI) + (1|A_ROI), After,REML=FALSE)
OT.anova <- anova(OT.null, OT.model1)
print(OT.anova)

OT.after.nostim<- glht(OT.model1, mcp(A_ROIType= "Tukey"))
summary(OT.after.nostim)



######

#percent astrocyte ROIs per field of view with at least one comparison

# find total number of ROIs per trial
# the find total number of ROIs after neurons, or before neurons
gcamp.nostim<- subset(nostim, Channel=="GCaMP")
nostim.OT.small$trials<- paste(nostim.OT.small$Animal, nostim.OT.small$Spot, nostim.OT.small$Trial, sep="_")
nostim.OT.small$ROIType<-nostim.OT.small$A_ROIType
# astrocyte ROIs per trial- total
ACROI_trial<- ddply(gcamp.nostim, c("Animal","Spot","trials"), summarise, AC_ROInum= length(unique(ROIs_trial)))

# astrocyte ROIs per trial- ROI type
ACROI_type_trial<- ddply(gcamp.nostim, c("Animal","Spot","trials","ROIType"), summarise, AC_ROInum= length(unique(ROIs_trial)))

# astrocyte ROIs with time differences around zero- total
ACROI_trial.OT.small<-ddply(nostim.OT.small, c("Animal","Spot","trials"), summarise, AC_ROInum= length(unique(A_ROI)))

# astrocyte ROIs with time differences around zero- ROI type
ACROI_type_trial.OT.small<- ddply(nostim.OT.small, c("Animal","Spot","trials","ROIType"), summarise, AC_ROInum= length(unique(A_ROI)))

# group data together for proportion calculation
numBefore<-ddply(Before, c("Animal","Spot","trials"), summarise, AC_ROInum= length(unique(A_ROI)))
numBefore.type<-ddply(Before, c("Animal","Spot","trials", "ROIType"), summarise, AC_ROInum= length(unique(A_ROI)))

numAfter<-ddply(After, c("Animal","Spot","trials"), summarise, AC_ROInum= length(unique(A_ROI)))
numAfter.type<-ddply(After, c("Animal","Spot","trials", "ROIType"), summarise, AC_ROInum= length(unique(A_ROI)))

numZero<-ddply(Zero, c("Animal","Spot","trials"), summarise, AC_ROInum= length(unique(A_ROI)))
numZero.type<-ddply(Zero, c("Animal","Spot","trials", "ROIType"), summarise, AC_ROInum= length(unique(A_ROI)))


#number of ROIs per trial 

Trialnames <-as.character(unique(ACROI_trial$trials))

#create dataframes that will be made
tot.ACnum.small <-data.frame()
bef.ACnum.small <-data.frame()
aft.ACnum.small <-data.frame()
zero.ACnum.small <-data.frame()

for (ii in 1:length(Trialnames))
{
  name =Trialnames[ii]
  subset1 = subset(ACROI_trial, trials == name) # total number of ROIs for this trial
  subset3 = subset(ACROI_trial.OT.small, trials == name) # number of ROIs with small time differences
  subset5 = subset(numBefore, trials == name) # number of ROIs with time before neurons- total
  subset7 = subset(numAfter, trials == name) # number of ROIs with time after neurons
  subset9 = subset(numZero, trials == name) # number of ROIs with zero

  # count total number of time diff ROIs per trial
  # use row data from trial if there are no comparisons
  if ((nrow(subset3) ==0)==TRUE)
  {subset3<- head(subset1,1)
  subset3$AC_ROInum = 0
  }
  subset3$numTimeDiffROIs<-0
  subset3$numTimeDiffROIs<-subset3$AC_ROInum/subset1$AC_ROInum
  tot.ACnum.small <-rbind(tot.ACnum.small, subset3)
  
  
  if ((nrow(subset5) ==0)==TRUE)
  {subset5<- head(subset1,1)
  subset5$AC_ROInum = 0
  }
  subset5$numTimeDiffROIs<-0
  subset5$numTimeDiffROIs<-subset5$AC_ROInum/subset1$AC_ROInum
  bef.ACnum.small  <-rbind(bef.ACnum.small, subset5)
  
  if ((nrow(subset7) ==0)==TRUE)
  {subset7<- head(subset1,1)
  subset7$AC_ROInum = 0
  }
  subset7$numTimeDiffROIs<-0
  subset7$numTimeDiffROIs<-subset7$AC_ROInum/subset1$AC_ROInum
  aft.ACnum.small  <-rbind(aft.ACnum.small, subset7)
  
  if ((nrow(subset9) ==0)==TRUE)
  {subset9<- head(subset1,1)
  subset9$AC_ROInum = 0
  }
  subset9$numTimeDiffROIs<-0
  subset9$numTimeDiffROIs<-subset9$AC_ROInum/subset1$AC_ROInum
  zero.ACnum.small  <-rbind(zero.ACnum.small, subset9)
  
}
 
df.numTotal<-summarySE(data=tot.ACnum.small, measurevar = "numTimeDiffROIs")
df.numBef <- summarySE(data=bef.ACnum.small, measurevar = "numTimeDiffROIs") 
df.numAft<- summarySE(data=aft.ACnum.small, measurevar = "numTimeDiffROIs") 
df.numZero <- summarySE(data=zero.ACnum.small, measurevar = "numTimeDiffROIs") 








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


