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

# peak data
long.peaks.74<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/74_peaks_longtrials.csv",  header=TRUE, sep = ",")
long.peaks.92<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/92_peaks_longtrials.csv",  header=TRUE, sep = ",")
long.peaks.94<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/94_peaks_longtrials.csv",  header=TRUE, sep = ",")
long.peaks.95<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/95_peaks_longtrials.csv",  header=TRUE, sep = ",")
long.peaks.96<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/96_peaks_longtrials.csv",  header=TRUE, sep = ",")
long.peaks.Alice<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/Alice_peaks_longtrials.csv",  header=TRUE, sep = ",")
#long.peaks.crazy8<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/crazy8_peaks_longtrials.csv",  header=TRUE, sep = ",")

# OT data
long.OT.74<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/74_onset_time_longtrials.csv",  header=TRUE, sep = ",")
long.OT.92<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/92_onset_time_longtrials.csv",  header=TRUE, sep = ",")
long.OT.94<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/94_onset_time_longtrials.csv",  header=TRUE, sep = ",")
long.OT.95<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/95_onset_time_longtrials.csv",  header=TRUE, sep = ",")
long.OT.96<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/96_onset_time_longtrials.csv",  header=TRUE, sep = ",")
long.OT.Alice<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/Alice_onset_time_longtrials.csv",  header=TRUE, sep = ",")
#long.OT.crazy8<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/crazy8_onset_time_longtrials.csv",  header=TRUE, sep = ",")

#field data
long.field.74<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/74_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
long.field.92<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/92_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
long.field.94<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/94_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
long.field.95<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/95_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
long.field.96<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/96_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
long.field.Alice<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/Alice_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
#long.field.crazy8<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/crazy8_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
##########

lsm.options(pbkrtest.limit = 100000)

#  merge data
all.OT<-rbind(long.OT.74, long.OT.92, long.OT.94,long.OT.95, long.OT.96, long.OT.Alice)#, long.OT.crazy8)
all.peaks<-rbind(long.peaks.74, long.peaks.92, long.peaks.94,long.peaks.95, long.peaks.96, long.peaks.Alice)#, long.peaks.crazy8)
all.field<-rbind(long.field.74, long.field.92, long.field.94,long.field.95, long.field.96, long.field.Alice)#, long.field.crazy8)

all.peaks$Duration<- all.peaks$halfWidth*2

#define shRNA for onset times
# 1- Control including Alice and Crazy8, 2- Alice and Crazy 8 are separate from 74 and 92
all.OT$shRNA1="Control"
all.OT$shRNA2="Control"

all.OT$shRNA1[grepl("94",all.OT$Animal)]="KD"
all.OT$shRNA1[grepl("95",all.OT$Animal)]="KD"
all.OT$shRNA1[grepl("96",all.OT$Animal)]="KD"
all.OT$shRNA2[grepl("94",all.OT$Animal)]="KD"
all.OT$shRNA2[grepl("95",all.OT$Animal)]="KD"
all.OT$shRNA2[grepl("96",all.OT$Animal)]="KD"
all.OT$shRNA2[grepl("74",all.OT$Animal)]="NS"
all.OT$shRNA2[grepl("92",all.OT$Animal)]="NS"

#define shRNA for peaks
# 1- Control including Alice and Crazy8, 2- Alice and Crazy 8 are separate from 74 and 92
all.peaks$shRNA1="Control"
all.peaks$shRNA2="Control"

all.peaks$shRNA1[grepl("94",all.peaks$Animal)]="KD"
all.peaks$shRNA1[grepl("95",all.peaks$Animal)]="KD"
all.peaks$shRNA1[grepl("96",all.peaks$Animal)]="KD"
all.peaks$shRNA2[grepl("94",all.peaks$Animal)]="KD"
all.peaks$shRNA2[grepl("95",all.peaks$Animal)]="KD"
all.peaks$shRNA2[grepl("96",all.peaks$Animal)]="KD"
all.peaks$shRNA2[grepl("74",all.peaks$Animal)]="NS"
all.peaks$shRNA2[grepl("92",all.peaks$Animal)]="NS"

#define shRNA for onset times
# 1- Control including Alice and Crazy8, 2- Alice and Crazy 8 are separate from 74 and 92
all.field$shRNA1="Control"
all.field$shRNA2="Control"

all.field$shRNA1[grepl("94",all.field$Animal)]="KD"
all.field$shRNA1[grepl("95",all.field$Animal)]="KD"
all.field$shRNA1[grepl("96",all.field$Animal)]="KD"
all.field$shRNA2[grepl("94",all.field$Animal)]="KD"
all.field$shRNA2[grepl("95",all.field$Animal)]="KD"
all.field$shRNA2[grepl("96",all.field$Animal)]="KD"
all.field$shRNA2[grepl("74",all.field$Animal)]="NS"
all.field$shRNA2[grepl("92",all.field$Animal)]="NS"

# set these new variables as factors so we can do stats on them
all.peaks$shRNA1<-as.factor(all.peaks$shRNA1)
all.peaks$shRNA2<-as.factor(all.peaks$shRNA2)
all.OT$shRNA1<-as.factor(all.OT$shRNA1)
all.OT$shRNA2<-as.factor(all.OT$shRNA2)
all.field$shRNA1<-as.factor(all.field$shRNA1)
all.field$shRNA2<-as.factor(all.field$shRNA2)

# set the order of groups for plots
all.OT$shRNA1<-factor(all.OT$shRNA1,levels=c("Control","KD"))
all.OT$shRNA2<-factor(all.OT$shRNA2,levels=c("Control","NS", "KD"))
all.peaks$shRNA1<-factor(all.peaks$shRNA1,levels=c("Control","KD"))
all.peaks$shRNA2<-factor(all.peaks$shRNA2,levels=c("Control","NS", "KD"))
all.field$shRNA1<-factor(all.field$shRNA1,levels=c("Control","KD"))
all.field$shRNA2<-factor(all.field$shRNA2,levels=c("Control","NS", "KD"))

#unique ROI names
all.OT$ROIs_trial<-paste(all.OT$Animal, all.OT$Spot, all.OT$Trial,all.OT$ROI, sep= "_")
all.OT$ROIs_trial_Cond<-paste(all.OT$ROIs_trial, all.OT$Condition, sep= "_")

all.peaks$ROIs_trial<-paste(all.peaks$Animal, all.peaks$Spot, all.peaks$Trial,all.peaks$roiName, sep= "_")
all.peaks$ROIs_trial_Cond<-paste(all.peaks$ROIs_trial, all.peaks$Condition, sep= "_")

all.field$Spot_trial<-paste(all.field$animalname, all.field$Spot, all.field$Trial, sep= "_")
all.field$Spot_trial_Cond<-paste(all.field$Spot_trial, all.field$Cond, sep= "_")


# count number of trials per spot
Spot.ntrials<-ddply(all.peaks, c("Animal","shRNA1","Spot","Condition"), summarise, nTrials=length(unique(Trial)))
Spot.ntrials$Ani_Spot_Cond<-paste(Spot.ntrials$Animal, Spot.ntrials$Spot, Spot.ntrials$Condition, sep="_")


########################
# ONSET TIME
#remove entries with no onset time

all.OT<-subset(all.OT, OnsetTime!="NaN")
#all.OT.astrocyte<- subset(all.OT, Channel=="GCaMP")
#all.OT.neuron<- subset(all.OT, Channel=="RCaMP")

###
# neuronal responses to stimulation
NeuronalStim<-subset(all.OT, Channel=="RCaMP" & Condition=="stim" & OnsetTime<8)

# should have an onset time in 8 s stimulus
ggplot(NeuronalStim,aes(x=OnsetTime,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP onset times between 0 and 8 s from stim trials")+
  max.theme

# should have an onset time in 8 s stimulus
ggplot(NeuronalStim,aes(x=OnsetTime,y=..density..,fill=shRNA2)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP onset times between 0 and 8 s from stim trials")+
  max.theme


# look at the percentiles for neuronal responses
Neuron95Onset.KD<-quantile(NeuronalStim$OnsetTime[(NeuronalStim$shRNA1=="KD")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
NeuronPT50.KD<-Neuron95Onset.KD[[11]]
print(NeuronPT50.KD)

#All Controls
Neuron95Onset.All_Con<-quantile(NeuronalStim$OnsetTime[(NeuronalStim$shRNA1=="Control")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
NeuronPT50.All_Con<-Neuron95Onset.All_Con[[11]]
print(NeuronPT50.All_Con)

# alice and crazy 8
Neuron95Onset.Control<-quantile(NeuronalStim$OnsetTime[(NeuronalStim$shRNA2=="Control")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
NeuronPT50.Control<-Neuron95Onset.Control[[11]]
print(NeuronPT50.Control)

# non-silencing
Neuron95Onset.NS<-quantile(NeuronalStim$OnsetTime[(NeuronalStim$shRNA2=="NS")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
NeuronPT50.NS<-Neuron95Onset.NS[[11]]
print(NeuronPT50.NS)


# all responding neurons
Neuron95Onset<-quantile(NeuronalStim$OnsetTime, prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
NeuronPT50<-Neuron95Onset[[11]]
print(NeuronPT50)

# time thresold to consider an astrocyte to be fast:
fastTh<-NeuronPT50

#######
#plot more distributions

# GCaMP onset times
ggplot(all.OT[(all.OT$Channel=="GCaMP"& all.OT$Condition=="stim"),],aes(x=OnsetTime,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP onset times")+
  max.theme

ggplot(all.OT[(all.OT$Channel=="GCaMP"& all.OT$Condition=="stim"),],aes(x=OnsetTime,y=..density..,fill=shRNA2)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP onset times")+
  max.theme


ggplot(all.peaks[(all.peaks$Channel=="GCaMP"& all.peaks$Condition=="stim"),],aes(x=amplitude,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP amplitudes")+
  max.theme

ggplot(all.peaks[(all.peaks$Channel=="RCaMP"& all.peaks$Condition=="stim"),],aes(x=amplitude,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP amplitudes")+
  max.theme

ggplot(all.peaks[(all.peaks$Channel=="RCaMP"& all.peaks$Condition=="stim"),],aes(x=amplitude,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP duration")+
  max.theme


#########
# find fast and delayed astrocyte microdomains


# remove ROIs with no peaks
all.peaks<-all.peaks[!(all.peaks$peakType=="NoPeak"),]

# add onset time information to the peak data table
all.peaks<-merge(all.peaks, all.OT[, c("ROIs_trial_Cond", "OnsetTime")], by="ROIs_trial_Cond", all.x=TRUE)


# remove matching dendrite and neuronal soma ROIs
Overlap= all.peaks$overlap!=0
all.peaks<-all.peaks[!Overlap,]


#####
# things to analyze
# field: percentage of active pixels, total number of ROIs per field per trial, total number of peaks per trial per field
# timing of signals- proportion of fast microdomains, onset time for fast MDs
# peaks- amplitde, duration of signals,
# is there a change in spontaneous activity in neurons or astrocytes in the no stim trials

# calculate the frequency (# of signals per trial)

########
#amplitude
df.amp1<-summarySE(all.peaks, measurevar = "amplitude", groupvars = c("Channel","Condition", "shRNA1"))
df.amp2<-summarySE(lck.peaks.window, measurevar = "amplitude", groupvars = c("Channel", "ROIType","Condition", "Genotype"))



ggplot(df.amp1, aes(x=interaction(Channel,shRNA1),y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  max.theme


