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

##############################

long.cyto.DSP4 <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cytoGC&RC_2D_longstim_28_04_2017.csv", header=TRUE, sep = ",")
short.cyto.DSP4<- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cytoGC&RC_2D_shortstim_28_04_2017.csv", header=TRUE, sep = ",")
nostim.cyto.DSP4<- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cytoGC&RC_2D_nostim_28_04_2017.csv", header=TRUE, sep = ",")


longstim.cyto.DSP4.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cyto_longstim_firstonset&AUC.csv", header=TRUE, sep = ",")
shortstim.cyto.DSP4.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cyto_shortstim_firstonset&AUC.csv", header=TRUE, sep = ",")
nostim.cyto.DSP4.OT <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cyto_nostim_firstonset&AUC.csv", header=TRUE, sep = ",")

#########################
#Home files

long.cyto.DSP4 <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cytoGC&RC_2D_longstim_28_04_2017.csv", header=TRUE, sep = ",")
short.cyto.DSP4<- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cytoGC&RC_2D_shortstim_28_04_2017.csv", header=TRUE, sep = ",")
nostim.cyto.DSP4<- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cytoGC&RC_2D_nostim_28_04_2017.csv", header=TRUE, sep = ",")


longstim.cyto.DSP4.OT <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cyto_longstim_firstonset&AUC.csv", header=TRUE, sep = ",")
shortstim.cyto.DSP4.OT <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cyto_shortstim_firstonset&AUC.csv", header=TRUE, sep = ",")
nostim.cyto.DSP4.OT <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cyto_nostim_firstonset&AUC.csv", header=TRUE, sep = ",")

########
# cytoGCaMP6s DSP4 data

all.cyto.DSP4.peaks<-rbind(nostim.cyto.DSP4,long.cyto.DSP4,short.cyto.DSP4)
all.cyto.DSP4.OT<-rbind(nostim.cyto.DSP4.OT,longstim.cyto.DSP4.OT,shortstim.cyto.DSP4.OT)

# remove the data frames that are combined
rm(nostim.cyto.DSP4, long.cyto.DSP4, short.cyto.DSP4, nostim.cyto.DSP4.OT, long.cyto.DSP4.OT, short.cyto.DSP4.OT)


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

## NOTE: define a stimulation window
# for distributions: histograms

stimwindow=15

stim.cyto.DSP4.OT.dist<-subset(all.cyto.DSP4.OT,Condition!="shortstim" & OnsetTime<stimwindow)
short.cyto.DSP4.OT.dist<-subset(all.cyto.DSP4.OT,Condition!="Stim" & OnsetTime<stimwindow)


# for median and mean calculations

# no stim vs 8 s stim- neuronal window=9s, AC window= 15 s for peak time, neuronal window=2s, AC window=12 s for onset
# no stim vs 1 s stim- neuronal window=2s, AC window= 10 s for peak time, neuronal window=2s, AC window=8 s for onset

LongN_PTwind2=9
LongAC_PTwind2=15

LongN_OTwind2=2
LongAC_OTwind2=12

Long_PTwind=12

Long_OTwind=10


ShortN_PTwind=2
ShortAC_PTwind=10

ShortN_OTwind=2
ShortAC_OTwind=8


#cyto DSP4
stim.cyto.DSP4.OT<-subset(all.cyto.DSP4.OT,Condition!="shortstim")
stim.cyto.DSP4.OT.R<-subset(stim.cyto.DSP4.OT, Channel=="RCaMP" & OnsetTime<=LongN_OTwind)
stim.cyto.DSP4.OT.G<-subset(stim.cyto.DSP4.OT, Channel=="GCaMP" & OnsetTime<=LongAC_OTwind)

stim.cyto.DSP4.OT.window<-rbind(stim.cyto.DSP4.OT.R, stim.cyto.DSP4.OT.G)

short.cyto.DSP4.OT<-subset(all.cyto.DSP4.OT,Condition!="Stim")
short.cyto.DSP4.OT.R<-subset(short.cyto.DSP4.OT, Channel=="RCaMP" & OnsetTime<=ShortN_OTwind)
short.cyto.DSP4.OT.G<-subset(short.cyto.DSP4.OT, Channel=="GCaMP" & OnsetTime<=ShortAC_OTwind)

short.cyto.DSP4.OT.window<-rbind(short.cyto.DSP4.OT.R, short.cyto.DSP4.OT.G)

#cyto DSP4
stim.cyto.DSP4.PT<-subset(all.cyto.DSP4.peaks,Condition!="shortstim")
stim.cyto.DSP4.PT.R<-subset(stim.cyto.DSP4.PT, Channel=="RCaMP" & peakTime<=LongN_PTwind & peakTime>=0 & Duration<50)
stim.cyto.DSP4.PT.G<-subset(stim.cyto.DSP4.PT, Channel=="GCaMP" & peakTime<=LongAC_PTwind & peakTime>=0 & Duration<50)

stim.cyto.DSP4.peaks.window<-rbind(stim.cyto.DSP4.PT.R, stim.cyto.DSP4.PT.G)

short.cyto.DSP4.PT<-subset(all.cyto.DSP4.peaks,Condition!="Stim")
short.cyto.DSP4.PT.R<-subset(short.cyto.DSP4.PT, Channel=="RCaMP" & peakTime<=ShortN_PTwind & peakTime>=0 & Duration<50)
short.cyto.DSP4.PT.G<-subset(short.cyto.DSP4.PT, Channel=="GCaMP" & peakTime<=ShortAC_PTwind & peakTime>=0 & Duration<50)

short.cyto.DSP4.peaks.window<-rbind(short.cyto.DSP4.PT.R, short.cyto.DSP4.PT.G)


#######

histseq= seq(0,15,0.5)
#cyto data

ntrials.cyto.DSP4.OT.long.R<- ddply(stim.cyto.DSP4.OT.dist[stim.cyto.DSP4.OT.dist$Channel=="RCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))
ntrials.cyto.DSP4.OT.long.G<- ddply(stim.cyto.DSP4.OT.dist[stim.cyto.DSP4.OT.dist$Channel=="GCaMP",], c("Condition"), summarise, ntrials=length(unique(Spot_trial)))


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
