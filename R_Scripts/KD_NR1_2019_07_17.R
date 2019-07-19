library("lme4")
library("lmerTest")
library("lattice")
library("plyr")
library("dplyr")
library("ggplot2")
library("emmeans")
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

# Jill's work files

# peak data
long.peaks.74<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/74_peaks_longtrials_bis.csv",  header=TRUE, sep = ",")
long.peaks.92<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/92_peaks_longtrials_bis.csv",  header=TRUE, sep = ",")
long.peaks.94<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/94_peaks_longtrials_bis.csv",  header=TRUE, sep = ",")
long.peaks.95<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/95_peaks_longtrials_bis.csv",  header=TRUE, sep = ",")
long.peaks.96<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/96_peaks_longtrials_bis.csv",  header=TRUE, sep = ",")
long.peaks.Alice<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/Alice_peaks_longtrials_bis.csv",  header=TRUE, sep = ",")
long.peaks.crazy8<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/crazy8_peaks_longtrials_bis.csv",  header=TRUE, sep = ",")

# OT data
long.OT.74<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/74_onset_time_longtrials_bis.csv",  header=TRUE, sep = ",")
long.OT.92<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/92_onset_time_longtrials_bis.csv",  header=TRUE, sep = ",")
long.OT.94<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/94_onset_time_longtrials_bis.csv",  header=TRUE, sep = ",")
long.OT.95<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/95_onset_time_longtrials_bis.csv",  header=TRUE, sep = ",")
long.OT.96<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/96_onset_time_longtrials_bis.csv",  header=TRUE, sep = ",")
long.OT.Alice<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/Alice_onset_time_longtrials_bis.csv",  header=TRUE, sep = ",")
long.OT.crazy8<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/crazy8_onset_time_longtrials_bis.csv",  header=TRUE, sep = ",")

#field data
long.field.74<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/74_Lck_field_longtrials_bis.csv",  header=TRUE, sep = ",")
long.field.92<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/92_Lck_field_longtrials_bis.csv",  header=TRUE, sep = ",")
long.field.94<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/94_Lck_field_longtrials_bis.csv",  header=TRUE, sep = ",")
long.field.95<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/95_Lck_field_longtrials_bis.csv",  header=TRUE, sep = ",")
long.field.96<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/96_Lck_field_longtrials_bis.csv",  header=TRUE, sep = ",")
long.field.Alice<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/Alice_Lck_field_longtrials_bis.csv",  header=TRUE, sep = ",")
long.field.crazy8<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/crazy8_Lck_field_longtrials_bis.csv",  header=TRUE, sep = ",")

##########

# Jill's home files

# peak data
long.peaks.74<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/74_peaks_longtrials_bis.csv",  header=TRUE, sep = ",")
long.peaks.92<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/92_peaks_longtrials_bis.csv",  header=TRUE, sep = ",")
long.peaks.94<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/94_peaks_longtrials_bis.csv",  header=TRUE, sep = ",")
long.peaks.95<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/95_peaks_longtrials_bis.csv",  header=TRUE, sep = ",")
long.peaks.96<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/96_peaks_longtrials_bis.csv",  header=TRUE, sep = ",")
long.peaks.Alice<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/Alice_peaks_longtrials_bis.csv",  header=TRUE, sep = ",")
long.peaks.crazy8<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/crazy8_peaks_longtrials_bis.csv",  header=TRUE, sep = ",")

# OT data
long.OT.74<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/74_onset_time_longtrials_bis.csv",  header=TRUE, sep = ",")
long.OT.92<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/92_onset_time_longtrials_bis.csv",  header=TRUE, sep = ",")
long.OT.94<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/94_onset_time_longtrials_bis.csv",  header=TRUE, sep = ",")
long.OT.95<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/95_onset_time_longtrials_bis.csv",  header=TRUE, sep = ",")
long.OT.96<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/96_onset_time_longtrials_bis.csv",  header=TRUE, sep = ",")
long.OT.Alice<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/Alice_onset_time_longtrials_bis.csv",  header=TRUE, sep = ",")
long.OT.crazy8<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/crazy8_onset_time_longtrials_bis.csv",  header=TRUE, sep = ",")

#field data
long.field.74<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/74_Lck_field_longtrials_bis.csv",  header=TRUE, sep = ",")
long.field.92<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/92_Lck_field_longtrials_bis.csv",  header=TRUE, sep = ",")
long.field.94<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/94_Lck_field_longtrials_bis.csv",  header=TRUE, sep = ",")
long.field.95<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/95_Lck_field_longtrials_bis.csv",  header=TRUE, sep = ",")
long.field.96<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/96_Lck_field_longtrials_bis.csv",  header=TRUE, sep = ",")
long.field.Alice<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/Alice_Lck_field_longtrials_bis.csv",  header=TRUE, sep = ",")
long.field.crazy8<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/complete_data/crazy8_Lck_field_longtrials_bis.csv",  header=TRUE, sep = ",")
##########

#lsm.options(pbkrtest.limit = 100000)
framerate=12.7917

#  merge data
all.OT<-rbind(long.OT.74, long.OT.92, long.OT.94,long.OT.95, long.OT.96, long.OT.Alice, long.OT.crazy8)
all.peaks<-rbind(long.peaks.74, long.peaks.92, long.peaks.94,long.peaks.95, long.peaks.96, long.peaks.Alice, long.peaks.crazy8)
all.field<-rbind(long.field.74, long.field.92, long.field.94,long.field.95, long.field.96, long.field.Alice, long.field.crazy8)

all.peaks$Duration<- all.peaks$halfWidth*2
#add baseline time to peaks table
all.peaks$BL_time<-5

# adjust peak time and duration
all.peaks$peakTime<- all.peaks$peakTime-all.peaks$BL_time  # now time= 0s is the start of stimulation
all.peaks$peakStart<- all.peaks$peakStart-all.peaks$BL_time
all.peaks$peakStartHalf<- all.peaks$peakStartHalf-all.peaks$BL_time


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

all.field$shRNA1[grepl("94",all.field$animalname)]="KD"
all.field$shRNA1[grepl("95",all.field$animalname)]="KD"
all.field$shRNA1[grepl("96",all.field$animalname)]="KD"
all.field$shRNA2[grepl("94",all.field$animalname)]="KD"
all.field$shRNA2[grepl("95",all.field$animalname)]="KD"
all.field$shRNA2[grepl("96",all.field$animalname)]="KD"
all.field$shRNA2[grepl("74",all.field$animalname)]="NS"
all.field$shRNA2[grepl("92",all.field$animalname)]="NS"

# set these new variables as factors so we can do stats on them
all.peaks$shRNA1<-as.factor(all.peaks$shRNA1)
all.peaks$shRNA2<-as.factor(all.peaks$shRNA2)
all.OT$shRNA1<-as.factor(all.OT$shRNA1)
all.OT$shRNA2<-as.factor(all.OT$shRNA2)
all.field$shRNA1<-as.factor(all.field$shRNA1)
all.field$shRNA2<-as.factor(all.field$shRNA2)

######
# days post injection
all.OT$DaysPostInjection="47"

all.OT$DaysPostInjection[grepl("06_15",all.OT$Spot)]="28"
all.OT$DaysPostInjection[grepl("06_21",all.OT$Spot)]="35"
all.OT$DaysPostInjection[grepl("06_22",all.OT$Spot)]="35"
all.OT$DaysPostInjection[grepl("06_28",all.OT$Spot)]="42"
all.OT$DaysPostInjection[grepl("07_02",all.OT$Spot)]="46"

#peaks
all.peaks$DaysPostInjection="47"

all.peaks$DaysPostInjection[grepl("06_15",all.peaks$Spot)]="28"
all.peaks$DaysPostInjection[grepl("06_21",all.peaks$Spot)]="35"
all.peaks$DaysPostInjection[grepl("06_22",all.peaks$Spot)]="35"
all.peaks$DaysPostInjection[grepl("06_28",all.peaks$Spot)]="42"
all.peaks$DaysPostInjection[grepl("07_02",all.peaks$Spot)]="46"

#field
all.field$DaysPostInjection="47"

all.field$DaysPostInjection[grepl("06_15",all.field$Spot)]="28"
all.field$DaysPostInjection[grepl("06_21",all.field$Spot)]="35"
all.field$DaysPostInjection[grepl("06_22",all.field$Spot)]="35"
all.field$DaysPostInjection[grepl("06_28",all.field$Spot)]="42"
all.field$DaysPostInjection[grepl("07_02",all.field$Spot)]="46"

# set these new variables as factors so we can do stats on them
all.peaks$DaysPostInjection<-as.factor(all.peaks$DaysPostInjection)
all.OT$DaysPostInjection<-as.factor(all.OT$DaysPostInjection)
all.field$DaysPostInjection<-as.factor(all.field$DaysPostInjection)


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


all.field$Spot_trial<-paste(all.field$animalname, all.field$Spot, all.field$trialname, sep= "_")
all.field$Spot_trial_Cond<-paste(all.field$Spot_trial, all.field$Cond, sep= "_")


# count number of trials per spot
Spot.ntrials<-ddply(all.peaks, c("Animal","shRNA1","Spot","Condition"), summarise, nTrials=length(unique(Trial)))
Spot.ntrials$Ani_Spot_Cond<-paste(Spot.ntrials$Animal, Spot.ntrials$Spot, Spot.ntrials$Condition, sep="_")


# remove matching dendrite and neuronal soma ROIs
Overlap= all.peaks$overlap!=0
all.peaks<-all.peaks[!Overlap,]

#consider onsly positive amplitude and accurate time scale
all.peaks.window<-subset(all.peaks, peakTime>0 & peakTime<12 & amplitude>0)

##########################
# FIELD DATA

######
# fraction of active pixels from all the astrocyte pixels
# with peaks near stimulus
all.field$nFluoPix[all.field$nFluoPix==0]<-128*128  #adjust spots where no pixels were above the threshold

all.field$FracActive=all.field$nActivePix/all.field$nFluoPix

all.field.spot<-ddply(all.field, c("Spot","animalname", "Cond", "shRNA1", "shRNA2", "DaysPostInjection", "Response_Score"), summarise, 
                      meanFracActive=mean(FracActive))

df.FracActive<- summarySE(all.field.spot, measurevar = "meanFracActive", groupvars = c("Cond","shRNA1"))
df.FracActive.day<-summarySE(all.field.spot, measurevar = "meanFracActive", groupvars = c("Cond","shRNA1", "DaysPostInjection" ))

df.ResponseScore<- summarySE(all.field.spot, measurevar = "Response_Score", groupvars = c("Cond","shRNA1"))
df.ResponseScore.day<- summarySE(all.field.spot, measurevar = "Response_Score", groupvars = c("Cond","shRNA1","DaysPostInjection"))

ggplot(data=df.FracActive, aes(x=shRNA1, y= meanFracActive, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanFracActive-se, ymax=meanFracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.FracActive.day, aes(x=interaction(DaysPostInjection, shRNA1), y= meanFracActive, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanFracActive-se, ymax=meanFracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.ResponseScore, aes(x=shRNA1, y= Response_Score, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Response_Score-se, ymax=Response_Score+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ResponseScore") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.ResponseScore.day, aes(x=interaction(DaysPostInjection,shRNA1), y= Response_Score, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Response_Score-se, ymax=Response_Score+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ResponseScore") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


#################################
# PEAK DATA

#####

# consider all peaks from the whole trial, regardless of stimulation time

#remove entries with no peaks
#all.peaks.window<-subset(all.peaks.window, peakTime!="NaN")


all.peaks.RC<- subset(all.peaks.window, Channel=="RCaMP")

all.peaks.GC<- subset(all.peaks.window, Channel=="GCaMP")

# amplitude histograms
ggplot(all.peaks.RC[(all.peaks.RC$Condition=="stim"),],aes(x=amplitude,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("RCaMP amplitudes")+
  max.theme

ggplot(all.peaks.GC[(all.peaks.GC$Condition=="stim"),],aes(x=amplitude,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.2,position="dodge") +
  ggtitle("GCaMP amplitudes")+
  max.theme

#duration histograms
ggplot(all.peaks.RC[(all.peaks.RC$Condition=="stim"),],aes(x=Duration,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=1, position="dodge") +
  ggtitle("RCaMP signal duration")+
  max.theme

ggplot(all.peaks.GC[(all.peaks.GC$Condition=="stim"),],aes(x=Duration,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=1,position="dodge") +
  ggtitle("GCaMP signal duration")+
  max.theme

#mean amplitude of signals for each ROI
df.amp1.RC<- summarySE(all.peaks.RC, measurevar = "amplitude", groupvars = c("Condition","shRNA1"))
df.amp2.RC<- summarySE(all.peaks.RC, measurevar = "amplitude", groupvars = c("Condition","shRNA2"))
df.amp.day.RC<-summarySE(all.peaks.RC, measurevar = "amplitude", groupvars = c("Condition","shRNA1", "DaysPostInjection" ))

df.amp1.GC<- summarySE(all.peaks.GC, measurevar = "amplitude", groupvars = c("Condition","shRNA1"))
df.amp2.GC<- summarySE(all.peaks.GC, measurevar = "amplitude", groupvars = c("Condition","shRNA2"))
df.amp.day.GC<-summarySE(all.peaks.GC, measurevar = "amplitude", groupvars = c("Condition","shRNA1", "DaysPostInjection" ))

ggplot(data=df.amp1.RC, aes(x=shRNA1, y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp1.GC, aes(x=shRNA1, y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp2.RC, aes(x=shRNA2, y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp2.GC, aes(x=shRNA2, y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp.day.RC, aes(x=interaction(shRNA1, DaysPostInjection), y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette) + 
  ggtitle("RCaMP") +
  max.theme

ggplot(data=df.amp.day.GC, aes(x=interaction(shRNA1, DaysPostInjection), y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette) + 
  ggtitle("GCaMP")+
  max.theme

## STATS

shRNA1_peaks.RC= interaction(all.peaks.RC$Condition,all.peaks.RC$shRNA1)
shRNA1_day_peaks.RC= interaction(all.peaks.RC$Condition,all.peaks.RC$shRNA1, 
                                 all.peaks.RC$DaysPostInjection)

shRNA1_peaks.GC= interaction(all.peaks.GC$Condition,all.peaks.GC$shRNA1)
shRNA1_day_peaks.GC= interaction(all.peaks.GC$Condition,all.peaks.GC$shRNA1, 
                                 all.peaks.GC$DaysPostInjection)

#RCaMP
amp.RC.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.RC,REML=FALSE)
amp.RC.model1 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.RC,REML=FALSE)
amp.RC.model2 = lmer(amplitude ~ shRNA1_peaks.RC + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.RC,REML=FALSE)
amp.RC.model3 = lmer(amplitude ~ shRNA1_day_peaks.RC+ (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.RC,REML=FALSE)
amp.RC.anova <- anova(amp.RC.null, amp.RC.model1, amp.RC.model2, amp.RC.model3)
print(amp.RC.anova)

amp.RC.shRNA1<- glht(amp.RC.model2, mcp(shRNA1_peaks.RC= "Tukey"))
summary(amp.RC.shRNA1)

amp.RC.shRNA1_day<- glht(amp.RC.model3, mcp(shRNA1_day_peaks.RC= "Tukey"))
summary(amp.RC.shRNA1_day)

#GCaMP
amp.GC.null = lmer(amplitude ~ (1|Animal) + (1|Spot), all.peaks.GC,REML=FALSE)
amp.GC.model1 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot), all.peaks.GC,REML=FALSE)
amp.GC.model2 = lmer(amplitude ~ shRNA1_peaks.GC + (1|Animal) + (1|Spot), all.peaks.GC,REML=FALSE)
amp.GC.model3 = lmer(amplitude ~ shRNA1_day_peaks.GC+ (1|Animal) + (1|Spot), all.peaks.GC,REML=FALSE)
amp.GC.anova <- anova(amp.GC.null, amp.GC.model1, amp.GC.model2, amp.GC.model3)
print(amp.GC.anova)

amp.GC.shRNA1<- glht(amp.GC.model2, mcp(shRNA1_peaks.GC= "Tukey"))
summary(amp.GC.shRNA1)

amp.GC.shRNA1_day<- glht(amp.GC.model3, mcp(shRNA1_day_peaks.GC= "Tukey"))
summary(amp.GC.shRNA1_day)

#outlier <- boxplot.stats(all.peaks.RC$amplitude)
#ROInum.8strial<-subset(ROInum.8strial, !(ROIsPerTrial %in% outlier))

#mean Duration
df.dur1.RC<- summarySE(all.peaks.RC, measurevar = "Duration", groupvars = c("Condition","shRNA1"))
df.dur2.RC<- summarySE(all.peaks.RC, measurevar = "Duration", groupvars = c("Condition","shRNA2"))
df.dur.day.RC<-summarySE(all.peaks.RC, measurevar = "Duration", groupvars = c("Condition","shRNA1", "DaysPostInjection" ))

df.dur1.GC<- summarySE(all.peaks.GC, measurevar = "Duration", groupvars = c("Condition","shRNA1"))
df.dur2.GC<- summarySE(all.peaks.GC, measurevar = "Duration", groupvars = c("Condition","shRNA2"))
df.dur.day.GC<-summarySE(all.peaks.GC, measurevar = "Duration", groupvars = c("Condition","shRNA1", "DaysPostInjection" ))

ggplot(data=df.dur1.RC, aes(x=shRNA1, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("RCdur") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

gplot(data=df.dur1.GC, aes(x=shRNA1, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("GCdur") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.dur2.RC, aes(x=shRNA2, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("RCdur") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.dur2.GC, aes(x=shRNA2, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("GCdur") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.dur.day.RC, aes(x=interaction(shRNA1, DaysPostInjection), y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
  ggtitle("RCdur") +
  max.theme

ggplot(data=df.dur.day.GC, aes(x=interaction(shRNA1, DaysPostInjection), y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
  ggtitle("GCdur")+
  max.theme

## STATS


#RCaMP
dur.RC.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.RC,REML=FALSE)
dur.RC.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.RC,REML=FALSE)
dur.RC.model2 = lmer(Duration ~ shRNA1_peaks.RC + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.RC,REML=FALSE)
dur.RC.model3 = lmer(Duration ~ shRNA1_day_peaks.RC+ (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.RC,REML=FALSE)
dur.RC.anova <- anova(dur.RC.null, dur.RC.model1, dur.RC.model2, dur.RC.model3)
print(dur.RC.anova)

dur.RC.shRNA1<- glht(dur.RC.model2, mcp(shRNA1_peaks.RC= "Tukey"))
summary(dur.RC.shRNA1)

dur.RC.shRNA1_day<- glht(dur.RC.model3, mcp(shRNA1_day_peaks.RC= "Tukey"))
summary(dur.RC.shRNA1_day)

#GCaMP
dur.GC.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.GC,REML=FALSE)
dur.GC.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.GC,REML=FALSE)
dur.GC.model2 = lmer(Duration ~ shRNA1_peaks.GC + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.GC,REML=FALSE)
dur.GC.model3 = lmer(Duration ~ shRNA1_day_peaks.GC+ (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.GC,REML=FALSE)
dur.GC.anova <- anova(dur.GC.null, dur.GC.model1, dur.GC.model2, dur.GC.model3)
print(dur.GC.anova)

dur.GC.shRNA1<- glht(dur.GC.model2, mcp(shRNA1_peaks.GC= "Tukey"))
summary(dur.GC.shRNA1)

dur.GC.shRNA1_day<- glht(dur.GC.model3, mcp(shRNA1_day_peaks.GC= "Tukey"))
summary(dur.GC.shRNA1_day)


##############
# which fields have neurons that respond to stimulation??


########################
# ONSET TIME
#remove entries with no onset time

all.OT<-subset(all.OT, OnsetTime!="NaN")
all.peaks<-subset(all.peaks, peaktime>0 & peaktime<12 & amplitude>0)

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




ggplot(all.peaks[(all.peaks$Channel=="RCaMP"& all.peaks$Condition=="stim"),],aes(x=Duration,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP duration")+
  max.theme

ggplot(all.peaks[(all.peaks$Channel=="GCaMP"& all.peaks$Condition=="stim"),],aes(x=Duration,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP duration")+
  max.theme

ggplot(all.peaks[(all.peaks$Channel=="GCaMP"& all.peaks$Condition=="stim"),],aes(x=peakTime,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP peak times")+
  max.theme

ggplot(all.peaks[(all.peaks$Channel=="RCaMP"& all.peaks$Condition=="stim"),],aes(x=peakTime,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP peak times")+
  max.theme


#####
# things to analyze
# field: percentage of active pixels, total number of ROIs per field per trial, total number of peaks per trial per field
# timing of signals- proportion of fast microdomains, onset time for fast MDs
# peaks- amplitde, duration of signals,
# is there a change in spontaneous activity in neurons or astrocytes in the no stim trials

# calculate the frequency (# of signals per trial)


#########
# find fast and delayed astrocyte microdomains

# isolate stim trials because we only care about the onset time after stimulation
stim.OT<- all.OT[all.OT$Condition=="stim",]

# identify "FAST" astrocytes
stim.OT$Group<-0
stim.OT$Group[stim.OT$OnsetTime<fastTh]<-"fast"
stim.OT$Group[stim.OT$OnsetTime>=fastTh]<-"delayed"

# remove ROIs with no peaks
all.peaks<-all.peaks[!(all.peaks$peakType=="NoPeak"),]

stim.peaks<-subset(all.peaks, Condition=="stim")

# add onset time information to the peak data table
stim.peaks<-merge(stim.peaks, stim.OT[, c("ROIs_trial_Cond", "OnsetTime", "Group")], by="ROIs_trial_Cond", all.x=TRUE)



###################

# number of active ROIs per field of view


all.peaks$Animal_Spot<- paste(all.peaks$Animal, all.peaks$Spot, sep="_")

# number of ROIs in each trial for each field of view (across the whole trial) with a peak during no stim and stim
# only consider during the stimulus (and the same time window in no stim trials)
all.peaks.8s<-all.peaks[(all.peaks$peakTime>0 & all.peaks$peakTime<15),]

all.peaks.8s$Channel <- factor(all.peaks.8s$Channel, levels = c("RCaMP","GCaMP"))

ROInum.8strial<-ddply(all.peaks.8s, c("Animal","Spot","shRNA1","shRNA2", "Condition","Channel","Animal_Spot","DaysPostInjection"), summarise, nROIs=length(unique(ROIs_trial_Cond)))

stim.peaks.8s<-stim.peaks[(stim.peaks$peakTime>0 & stim.peaks$peakTime<15),]
ROInum.8strial.group<-ddply(stim.peaks.8s, c("Animal","Spot","shRNA1","shRNA2","Condition","Channel","DaysPostInjection","Group"), summarise, nROIs=length(unique(ROIs_trial_Cond)))


# add in number of trials
ROInum.8strial$Ani_Spot_Cond<-paste(ROInum.8strial$Animal_Spot, ROInum.8strial$Condition, sep="_")
ROInum.8strial.group$Ani_Spot_Cond<-paste(ROInum.8strial.group$Animal, ROInum.8strial.group$Spot, ROInum.8strial.group$Condition, sep="_")

ROInum.8strial<-merge(ROInum.8strial, Spot.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)
ROInum.8strial$ROIsPerTrial<-ROInum.8strial$nROIs/ROInum.8strial$nTrials

ROInum.8strial.group<-merge(ROInum.8strial.group, Spot.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)
ROInum.8strial.group$ROIsPerTrial<-ROInum.8strial.group$nROIs/ROInum.8strial.group$nTrials

#remove NaNs
ROInum.8strial<-subset(ROInum.8strial, Spot!="NaN")
ROInum.8strial.group<-subset(ROInum.8strial.group, Group!="NaN")

# mean
df.ROInum.8strial<-summarySE(ROInum.8strial, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition","shRNA1"))
df.ROInum.8strial.day<-summarySE(ROInum.8strial, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition","shRNA1","DaysPostInjection"))
df.ROInum.8strial.group<-summarySE(ROInum.8strial.group, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition","shRNA1","Group"))
df.ROInum.8strial.group.day<-summarySE(ROInum.8strial.group, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition","shRNA1","Group","DaysPostInjection"))

# plots
# RCaMP
ggplot(df.ROInum.8strial[df.ROInum.8strial$Channel=="RCaMP",], aes(x=shRNA1,y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("RCaMP ROIs")+
  max.theme

# GCaMP
ggplot(df.ROInum.8strial[df.ROInum.8strial$Channel=="GCaMP",], aes(x=shRNA1,y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("GCaMP ROIs")+
  max.theme


# GCaMP fast vs delayed
ggplot(df.ROInum.8strial.group[df.ROInum.8strial.group$Channel=="GCaMP",], aes(x=Group,y=ROIsPerTrial, fill= shRNA1)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("GCaMP ROIs")+
  max.theme

# paired line plots
df.ROInum.8strial$ROIsPerTrialMean<-df.ROInum.8strial$ROIsPerTrial
df.ROInum.8strial$Chan_Cond_shRNA<-interaction(df.ROInum.8strial$Channel, df.ROInum.8strial$Condition,df.ROInum.8strial$shRNA1)
ROInum.8strial$Chan_Cond_shRNA<-interaction(ROInum.8strial$Channel, ROInum.8strial$Condition,ROInum.8strial$shRNA1)

ROInum.8strial<-merge(ROInum.8strial, df.ROInum.8strial[, c("Chan_Cond_shRNA", "ROIsPerTrialMean","se")], by="Chan_Cond_shRNA", all.x=TRUE)

df.ROInum.8strial.group$ROIsPerTrialMean<-df.ROInum.8strial.group$ROIsPerTrial
df.ROInum.8strial.group$Chan_Cond_shRNA_Group<-interaction(df.ROInum.8strial.group$Channel, 
                                                              df.ROInum.8strial.group$Condition,
                                                              df.ROInum.8strial.group$shRNA1,
                                                              df.ROInum.8strial.group$Group)
ROInum.8strial.group$Chan_Cond_shRNA_Group<-interaction(ROInum.8strial.group$Channel, 
                                                           ROInum.8strial.group$Condition,
                                                           ROInum.8strial.group$shRNA1,
                                                           ROInum.8strial.group$Group)

ROInum.8strial.group<-merge(ROInum.8strial.group, df.ROInum.8strial.group[, c("Chan_Cond_shRNA_Group", "ROIsPerTrialMean","se")], by="Chan_Cond_shRNA_Group", all.x=TRUE)

# plots
ggplot(df.ROInum.8strial[df.ROInum.8strial$Channel=="RCaMP",],aes(x=shRNA1,y=ROIsPerTrial, group=Condition)) +
  geom_point(aes(x=shRNA1,y=ROIsPerTrial),stat="identity", size=3)+
  geom_line(aes(x=shRNA1, y=ROIsPerTrial, colour=Condition))+
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("RCaMP ROIs")+
  max.theme

# paired plot control
ggplot(ROInum.8strial[(ROInum.8strial$Channel=="RCaMP" & ROInum.8strial$shRNA1=="Control"),], aes(x=Condition, y = ROIsPerTrial)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(x=Condition, y=ROIsPerTrial,group=Animal_Spot), colour="#b5b5b5")+
  ggtitle("RCaMP control")+
  geom_point(aes(x=Condition, y=ROIsPerTrialMean), size = 5, colour="#7b3294")+
  geom_line(aes(x=Condition, y=ROIsPerTrialMean,group=Animal_Spot), size=1.5, colour="#7b3294")+
  geom_errorbar(aes(x=Condition,ymin=ROIsPerTrialMean-se, ymax=ROIsPerTrialMean+se), colour="#7b3294", width=0.2,  size=1.5,position=position_dodge(.9)) +
  max.theme

ggplot(ROInum.8strial[(ROInum.8strial$Channel=="RCaMP" & ROInum.8strial$shRNA1=="KD"),], aes(x=Condition, y = ROIsPerTrial)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(x=Condition, y=ROIsPerTrial,group=Animal_Spot), colour="#b5b5b5")+
  ggtitle("RCaMP knockdown")+
  geom_point(aes(x=Condition, y=ROIsPerTrialMean), size = 5, colour="#7b3294")+
  geom_line(aes(x=Condition, y=ROIsPerTrialMean,group=Animal_Spot), size=1.5, colour="#7b3294")+
  geom_errorbar(aes(x=Condition,ymin=ROIsPerTrialMean-se, ymax=ROIsPerTrialMean+se), colour="#7b3294", width=0.2,  size=1.5,position=position_dodge(.9)) +
  max.theme

#GCaMP
ggplot(df.ROInum.8strial[df.ROInum.8strial$Channel=="GCaMP",],aes(x=shRNA1,y=ROIsPerTrial, group=Condition)) +
  geom_point(aes(x=shRNA1,y=ROIsPerTrial),stat="identity", size=3)+
  geom_line(aes(x=shRNA1, y=ROIsPerTrial, colour=Condition))+
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("GCaMP ROIs")+
  max.theme

# paired plot
ggplot(ROInum.8strial[(ROInum.8strial$Channel=="GCaMP" & ROInum.8strial$shRNA1=="Control"),], aes(x=Condition, y = ROIsPerTrial)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(x=Condition, y=ROIsPerTrial,group=Animal_Spot), colour="#b5b5b5")+
  ggtitle("GCaMP control")+
  geom_point(aes(x=Condition, y=ROIsPerTrialMean), size = 5, colour="#008837")+
  geom_line(aes(x=Condition, y=ROIsPerTrialMean,group=Animal_Spot), size=1.5, colour="#008837")+
  geom_errorbar(aes(x=Condition,ymin=ROIsPerTrialMean-se, ymax=ROIsPerTrialMean+se), colour="#008837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  max.theme

ggplot(ROInum.8strial[(ROInum.8strial$Channel=="GCaMP" & ROInum.8strial$shRNA1=="KD"),], aes(x=Condition, y = ROIsPerTrial)) +
  geom_point(shape = 21,size = 3, colour="#b5b5b5") +
  geom_line(aes(x=Condition, y=ROIsPerTrial,group=Animal_Spot), colour="#b5b5b5")+
  ggtitle("GCaMP knockdown")+
  geom_point(aes(x=Condition, y=ROIsPerTrialMean), size = 5, colour="#008837")+
  geom_line(aes(x=Condition, y=ROIsPerTrialMean,group=Animal_Spot), size=1.5, colour="#008837")+
  geom_errorbar(aes(x=Condition,ymin=ROIsPerTrialMean-se, ymax=ROIsPerTrialMean+se), colour="#008837", width=0.2,  size=1.5,position=position_dodge(.9)) +
  max.theme


# only fast GcaMP

fastROInum<- subset(ROInum.8strial.group, Channel=="GCaMP" & Condition=="stim" & Group=="fast")

ggplot(fastROInum, aes(x=shRNA1,y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("Lck- fast ROIs per FOV") +
  max.theme

ggplot(df.ROInum.8strial.group[df.ROInum.8strial.group$Channel=="GCaMP" & 
                                 df.ROInum.8strial.group$Group=="fast" &
                                 df.ROInum.8strial.group$Condition=="stim",],
       aes(x=shRNA1,y=ROIsPerTrial, colour=shRNA1)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("fast GCaMP ROIs") +
  max.theme

# only delayed GcaMP

delayedROInum<- subset(ROInum.8strial.group, Channel=="GCaMP" & Condition=="stim" & Group=="delayed")

ggplot(delayedROInum, aes(x=shRNA1,y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("Lck- delayed ROIs per FOV") +
  max.theme

ggplot(df.ROInum.8strial.group[df.ROInum.8strial.group$Channel=="GCaMP" & 
                                 df.ROInum.8strial.group$Group=="delayed" &
                                 df.ROInum.8strial.group$Condition=="stim",],
       aes(x=shRNA1,y=ROIsPerTrial, colour=shRNA1)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("delayed GCaMP ROIs") +
  max.theme

# days post injection
# RCaMP
ggplot(df.ROInum.8strial.day[df.ROInum.8strial.day$Channel=="RCaMP",], aes(x=interaction(DaysPostInjection,shRNA1),y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("RCaMP ROIs")+
  max.theme

# GCaMP
ggplot(df.ROInum.8strial.day[df.ROInum.8strial.day$Channel=="GCaMP",], aes(x=interaction(DaysPostInjection,shRNA1),y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("GCaMP ROIs")+
  max.theme


# GCaMP fast vs delayed
ggplot(df.ROInum.8strial.group.day[df.ROInum.8strial.group.day$Channel=="GCaMP",], aes(x=interaction(DaysPostInjection,Group),y=ROIsPerTrial, fill= shRNA1)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("GCaMP ROIs")+
  max.theme


#line plots
ggplot(df.ROInum.8strial.day[(df.ROInum.8strial.day$Channel=="RCaMP" & df.ROInum.8strial.day$Condition=="stim"),],aes(x=DaysPostInjection,y=ROIsPerTrial, group=shRNA1)) +
  geom_point(aes(x=DaysPostInjection,y=ROIsPerTrial),stat="identity", size=3)+
  geom_line(aes(x=DaysPostInjection, y=ROIsPerTrial, colour=shRNA1))+
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("RCaMP ROIs- Stim")+
  max.theme

ggplot(df.ROInum.8strial.day[(df.ROInum.8strial.day$Channel=="RCaMP" & df.ROInum.8strial.day$Condition=="nostim"),],aes(x=DaysPostInjection,y=ROIsPerTrial, group=shRNA1)) +
  geom_point(aes(x=DaysPostInjection,y=ROIsPerTrial),stat="identity", size=3)+
  geom_line(aes(x=DaysPostInjection, y=ROIsPerTrial, colour=shRNA1))+
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("RCaMP ROIs- No Stim")+
  max.theme

ggplot(df.ROInum.8strial.day[(df.ROInum.8strial.day$Channel=="GCaMP" & df.ROInum.8strial.day$Condition=="stim"),],aes(x=DaysPostInjection,y=ROIsPerTrial, group=shRNA1)) +
  geom_point(aes(x=DaysPostInjection,y=ROIsPerTrial),stat="identity", size=3)+
  geom_line(aes(x=DaysPostInjection, y=ROIsPerTrial, colour=shRNA1))+
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("GCaMP ROIs- Stim")+
  max.theme

ggplot(df.ROInum.8strial.day[(df.ROInum.8strial.day$Channel=="GCaMP" & df.ROInum.8strial.day$Condition=="nostim"),],aes(x=DaysPostInjection,y=ROIsPerTrial, group=shRNA1)) +
  geom_point(aes(x=DaysPostInjection,y=ROIsPerTrial),stat="identity", size=3)+
  geom_line(aes(x=DaysPostInjection, y=ROIsPerTrial, colour=shRNA1))+
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("GCaMP ROIs- No Stim")+
  max.theme

###############
# stats for active ROI number per trials per FOV

# RCaMP
Condition_shRNA1_RC= interaction(ROInum.8strial$Condition[ROInum.8strial$Channel=="RCaMP"],ROInum.8strial$shRNA1[ROInum.8strial$Channel=="RCaMP"])
nROI.RC.stim.null = lmer(ROIsPerTrial ~ (1|Animal) + (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="RCaMP",],REML=FALSE)
nROI.RC.stim.model1 = lmer(ROIsPerTrial ~ Condition + (1|Animal)+ (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="RCaMP",],REML=FALSE)
nROI.RC.stim.model2 = lmer(ROIsPerTrial ~ Condition_shRNA1_RC + (1|Animal)+ (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="RCaMP",],REML=FALSE)
nROI.RC.stim.anova <- anova(nROI.RC.stim.null, nROI.RC.stim.model1,nROI.RC.stim.model2)
print(nROI.RC.stim.anova)

nROI.RC.stim.Cond_shRNA1<- glht(nROI.RC.stim.model2, mcp(Condition_shRNA1_RC= "Tukey"))
summary(nROI.RC.stim.Cond_shRNA1)


# GCaMP KD vs controls, no stim vs stim
Condition_shRNA1_GC= interaction(ROInum.8strial$Condition[ROInum.8strial$Channel=="GCaMP"],
                                   ROInum.8strial$shRNA1[ROInum.8strial$Channel=="GCaMP"])
Condition_shRNA1_day_GC= interaction(ROInum.8strial$Condition[ROInum.8strial$Channel=="GCaMP"],
                                 ROInum.8strial$shRNA1[ROInum.8strial$Channel=="GCaMP"],
                                 ROInum.8strial$DaysPostInjection[ROInum.8strial$Channel=="GCaMP"])

nROI.GC.stim.null = lmer(ROIsPerTrial ~ (1|Animal) + (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="GCaMP",],REML=FALSE)
nROI.GC.stim.model1 = lmer(ROIsPerTrial ~ Condition + (1|Animal)+ (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="GCaMP",],REML=FALSE)
nROI.GC.stim.model2 = lmer(ROIsPerTrial ~ shRNA1 + (1|Animal)+ (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="GCaMP",],REML=FALSE)
nROI.GC.stim.model3 = lmer(ROIsPerTrial ~ Condition_shRNA1_GC + (1|Animal)+ (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="GCaMP",],REML=FALSE)
nROI.GC.stim.model4 = lmer(ROIsPerTrial ~ Condition_shRNA1_day_GC + (1|Animal)+ (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="GCaMP",],REML=FALSE)
nROI.GC.stim.anova <- anova(nROI.GC.stim.null, nROI.GC.stim.model1,
                            nROI.GC.stim.model2,nROI.GC.stim.model3,
                            nROI.GC.stim.model4)
print(nROI.GC.stim.anova)

nROI.GC.shRNA1<- glht(nROI.GC.stim.model2, mcp(shRNA1= "Tukey"))
summary(nROI.GC.shRNA1)

nROI.GC.stim.Cond_shRNA1<- glht(nROI.GC.stim.model3, mcp(Condition_shRNA1_GC= "Tukey"))
summary(nROI.GC.stim.Cond_shRNA1)


# groups
GC.stim<-subset(ROInum.8strial.group, Channel=="GCaMP" & Condition=="stim")

shRNA1_Group_GC= interaction(GC.stim$shRNA1, GC.stim$Group)
shRNA1_Group_day_GC= interaction(GC.stim$shRNA1, GC.stim$Group, GC.stim$DaysPostInjection)

nROI.GC.group.null = lmer(ROIsPerTrial ~ (1|Animal) + (1|Spot), GC.stim,REML=FALSE)
nROI.GC.group.model1 = lmer(ROIsPerTrial ~ Group + (1|Animal)+ (1|Spot), GC.stim,REML=FALSE)
nROI.GC.group.model2 = lmer(ROIsPerTrial ~ shRNA1 + (1|Animal)+ (1|Spot), GC.stim,REML=FALSE)
nROI.GC.group.model3 = lmer(ROIsPerTrial ~ shRNA1_Group_GC + (1|Animal)+ (1|Spot), GC.stim,REML=FALSE)
nROI.GC.group.model4 = lmer(ROIsPerTrial ~ shRNA1_Group_day_GC + (1|Animal)+ (1|Spot), GC.stim,REML=FALSE)
nROI.GC.group.anova <- anova(nROI.GC.group.null, nROI.GC.group.model1,
                             nROI.GC.group.model2,nROI.GC.group.model3,
                             nROI.GC.group.model4)
print(nROI.GC.group.anova)

nROI.GC.shRNA1_Group<- glht(nROI.GC.group.model3, mcp(shRNA1_Group_GC= "Tukey"))
summary(nROI.GC.shRNA1_Group)

nROI.GC.shRNA1_Group.day<- glht(nROI.GC.group.model4, mcp(shRNA1_Group_day_GC= "Tukey"))
summary(nROI.GC.shRNA1_Group.day)


# fast GCaMP only
GC.stim.fast<-subset(GC.stim, Group=="fast")

shRNA1_day= interaction(GC.stim.fast$shRNA1, GC.stim.fast$DaysPostInjection)

nROI.fast.day.null = lmer(ROIsPerTrial ~ (1|Animal), GC.stim.fast,REML=FALSE)
nROI.fast.day.model1 = lmer(ROIsPerTrial ~ shRNA1 + (1|Animal), GC.stim.fast,REML=FALSE)
nROI.fast.day.model2 = lmer(ROIsPerTrial ~ shRNA1_day + (1|Animal), GC.stim.fast,REML=FALSE)
nROI.fast.day.anova <- anova(nROI.fast.day.null, nROI.fast.day.model1,
                             nROI.fast.day.model2)
print(nROI.fast.day.anova)

nROI.GC.shRNA1_day<- glht(nROI.fast.day.model2, mcp(shRNA1_day= "Tukey"))
summary(nROI.GC.shRNA1_day)




