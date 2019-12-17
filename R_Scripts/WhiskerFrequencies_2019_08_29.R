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
#cbbPalette <- c("#D55E00","#009E73","#E69F00","#56B4E9","#CC79A7","#F0E442")
cbbPalette <- c("#D55E00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#E69F00", "#CC79A7")

########################
# load data

# Jill's home files
# peak data
long.peaks.Alice<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/WhiskerFrequencies/Results/FilesforR/Alice_peaks_frequencies.csv",  header=TRUE, sep = ",")
long.peaks.crazy8<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/WhiskerFrequencies/Results/FilesforR/Crazy8_peaks_frequencies.csv",  header=TRUE, sep = ",")

# OT data
long.OT.Alice<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/WhiskerFrequencies/Results/FilesforR/Alice_onset_time_frequencies.csv",  header=TRUE, sep = ",")
long.OT.crazy8<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/WhiskerFrequencies/Results/FilesforR/Crazy8_onset_time_frequencies.csv",  header=TRUE, sep = ",")

#field data
long.field.Alice<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/WhiskerFrequencies/Results/FilesforR/Alice_Lck_field_frequencies.csv",  header=TRUE, sep = ",")
long.field.crazy8<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/WhiskerFrequencies/Results/FilesforR/Crazy8_Lck_field_frequencies.csv",  header=TRUE, sep = ",")


###############
# Jill's work files
# peak data
long.peaks.Alice<-read.table("G:/Data/GCaMP_RCaMP/WhiskerFrequencies/Results/FilesforR/Alice_peaks_frequencies.csv",  header=TRUE, sep = ",")
long.peaks.crazy8<-read.table("G:/Data/GCaMP_RCaMP/WhiskerFrequencies/Results/FilesforR/Crazy8_peaks_frequencies.csv",  header=TRUE, sep = ",")

# OT data
long.OT.Alice<-read.table("G:/Data/GCaMP_RCaMP/WhiskerFrequencies/Results/FilesforR/Alice_onset_time_frequencies.csv",  header=TRUE, sep = ",")
long.OT.crazy8<-read.table("G:/Data/GCaMP_RCaMP/WhiskerFrequencies/Results/FilesforR/Crazy8_onset_time_frequencies.csv",  header=TRUE, sep = ",")

#field data
long.field.Alice<-read.table("G:/Data/GCaMP_RCaMP/WhiskerFrequencies/Results/FilesforR/Alice_Lck_field_frequencies.csv",  header=TRUE, sep = ",")
long.field.crazy8<-read.table("G:/Data/GCaMP_RCaMP/WhiskerFrequencies/Results/FilesforR/Crazy8_Lck_field_frequencies.csv",  header=TRUE, sep = ",")



##########

#lsm.options(pbkrtest.limit = 100000)
framerate=12.7917

#  merge data
all.OT<-rbind(long.OT.Alice, long.OT.crazy8)
all.peaks<-rbind(long.peaks.Alice, long.peaks.crazy8)
all.field<-rbind(long.field.Alice, long.field.crazy8)

all.peaks$Duration<- all.peaks$halfWidth*2
#add baseline time to peaks table
all.peaks$BL_time<-2

# adjust peak time and duration
all.peaks$peakTime<- all.peaks$peakTime-all.peaks$BL_time  # now time= 0s is the start of stimulation
all.peaks$peakStart<- all.peaks$peakStart-all.peaks$BL_time
all.peaks$peakStartHalf<- all.peaks$peakStartHalf-all.peaks$BL_time




#unique ROI names
all.OT$ROIs_trial<-paste(all.OT$Animal, all.OT$Spot, all.OT$Trial,all.OT$ROI, sep= "_")
all.OT$ROIs_trial_Cond<-paste(all.OT$ROIs_trial, all.OT$Condition, sep= "_")

all.peaks$ROIs_trial<-paste(all.peaks$Animal, all.peaks$Spot, all.peaks$Trial,all.peaks$roiName, sep= "_")
all.peaks$ROIs_trial_Cond<-paste(all.peaks$ROIs_trial, all.peaks$Condition, sep= "_")


all.field$Spot_trial<-paste(all.field$animalname, all.field$Spot, all.field$trialname, sep= "_")
all.field$Spot_trial_Cond<-paste(all.field$Spot_trial, all.field$Cond, sep= "_")
all.field$Spot_animals<-paste(all.field$Spot, all.field$animalname, sep= "_")

# count number of trials per spot
Spot.ntrials<-ddply(all.peaks, c("Animal","Spot","Condition"), summarise, nTrials=length(unique(Trial)))
Spot.ntrials$Ani_Spot_Cond<-paste(Spot.ntrials$Animal, Spot.ntrials$Spot, Spot.ntrials$Condition, sep="_")


# remove matching dendrite and neuronal soma ROIs
Overlap= all.peaks$overlap!=0
all.peaks<-all.peaks[!Overlap,]

#consider onsly positive amplitude and accurate time scale
all.peaks.window<-subset(all.peaks, peakTime>0 & peakTime<7 & amplitude>0 & Duration<10)

##########################
# FIELD DATA

# remove random trials of NaNs
all.field<-all.field[!is.na(all.field$pixelsize),]

######
# fraction of active pixels from all the astrocyte pixels
# with peaks near stimulus
all.field$nFluoPix[all.field$nFluoPix==0]<-all.field$nTotalPix[all.field$nFluoPix==0]  #adjust spots where no pixels were above the threshold

# 

all.field$FracActive=all.field$nActivePix/all.field$nFluoPix
all.field$FracFluo=all.field$nFluoPix/all.field$nTotalPix

# account for the fact that this was less time
all.field$FracActive<-all.field$FracActive*2  #fraction per minute!!!!!
all.field$FracFluo<-all.field$FracFluo*2  #fraction per minute!!!!!

all.field.spot<-ddply(all.field, c("Spot","animalname", "Cond", "Response_Score"), summarise, 
                      meanFracActive=mean(FracActive), meanFluoPix = mean(nFluoPix), meanFracFluo = mean(FracFluo))

df.FracActive<- summarySE(all.field.spot, measurevar = "meanFracActive", groupvars = c("Cond"))

df.ResponseScore<- summarySE(all.field.spot, measurevar = "Response_Score", groupvars = c("Cond"))

df.nFluoPix<-summarySE(all.field.spot, measurevar = "meanFluoPix", groupvars = c("Cond"))



# graphs fraction of active pixels
ggplot(data=df.FracActive, aes(x=Cond, y= meanFracActive, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanFracActive-se, ymax=meanFracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



# graphs response probabliity (potential for same region to appear in mutliple trials)
ggplot(data=df.ResponseScore, aes(x=Cond, y= Response_Score, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Response_Score-se, ymax=Response_Score+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ResponseScore") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


######
#stats

field.FA.null = lmer(FracActive ~ (1|animalname) + (1|Spot), all.field,REML=FALSE)
field.FA.model1 = lmer(FracActive ~ Cond + (1|animalname) + (1|Spot), all.field,REML=FALSE)
field.FA.anova <- anova(field.FA.null, field.FA.model1)
print(field.FA.anova)

field.FA.pv<- glht(field.FA.model1, mcp(Cond= "Tukey"))
summary(field.FA.pv)


field.RR.null = lmer(Response_Score ~ (1|animalname) + (1|Spot), all.field,REML=FALSE)
field.RR.model1 = lmer(Response_Score ~ Cond + (1|animalname) + (1|Spot), all.field,REML=FALSE)
field.RR.anova <- anova(field.RR.null, field.RR.model1)
print(field.RR.anova)

field.RR.pv<- glht(field.RR.model1, mcp(Cond= "Tukey"))
summary(field.RR.pv)

#################################
# PEAK DATA 

#for peaks that make sense (duration less than trial length, etc.)

all.peaks<-subset(all.peaks, amplitude>0 & Duration<10 & peakTime<5)


# RCaMP peaks (nostim and stim) near stimulation time
all.peaks.RC<- subset(all.peaks.window, Channel=="RCaMP")


# GCaMP peaks (nostim and stim) near stimulation time
all.peaks.GC<- subset(all.peaks.window, Channel=="GCaMP")

######
# amplitude histograms
ggplot(all.peaks.RC,aes(x=amplitude,y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("RCaMP stim amplitudes")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(all.peaks.GC,aes(x=amplitude,y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.2,position="dodge") +
  ggtitle("GCaMP stim amplitudes")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme


#duration histograms
ggplot(all.peaks.RC,aes(x=Duration,y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("RCaMP signal duration")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(all.peaks.GC,aes(x=Duration,y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.2,position="dodge") +
  ggtitle("GCaMP signal duration")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# peak time
ggplot(all.peaks.GC,aes(x=peakTime,y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP peak times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(all.peaks.RC,aes(x=peakTime,y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP peak times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme


#############

#mean amplitude of signals for each ROI

#RCaMP
df.amp1.RC<- summarySE(all.peaks.RC, measurevar = "amplitude", groupvars = c("Condition"))

# RCaMP graphs
ggplot(data=df.amp1.RC, aes(x=Condition, y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


# GCaMP amplitude
df.amp1.GC<- summarySE(all.peaks.GC, measurevar = "amplitude", groupvars = c("Condition"))



ggplot(data=df.amp1.GC, aes(x=Condition, y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



#####################
## STATS

#RCaMP 
amp.RC.S.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.RC,REML=FALSE)
amp.RC.S.model1 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.RC,REML=FALSE)
amp.RC.S.anova <- anova(amp.RC.S.null, amp.RC.S.model1)
print(amp.RC.S.anova)

amp.RC.S.pv<- glht(amp.RC.S.model1, mcp(Condition= "Tukey"))
summary(amp.RC.S.pv)


#GCaMP stim
amp.GC.S.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.GC,REML=FALSE)
amp.GC.S.model1 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.GC,REML=FALSE)
amp.GC.S.anova <- anova(amp.GC.S.null, amp.GC.S.model1)
print(amp.GC.S.anova)

amp.GC.S.pv<- glht(amp.GC.S.model1, mcp(Condition= "Tukey"))
summary(amp.GC.S.pv)


#########
#mean Duration of signals for each ROI

#RCaMP
df.dur1.RC<- summarySE(all.peaks.RC, measurevar = "Duration", groupvars = c("Condition"))


# RCaMP graphs
ggplot(data=df.dur1.RC, aes(x=Condition, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


# GCaMP Duration
df.dur1.GC<- summarySE(all.peaks.GC, measurevar = "Duration", groupvars = c("Condition"))


ggplot(data=df.dur1.GC, aes(x=Condition, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



#####################
## STATS

#RCaMP nostim
dur.RC.S.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.RC,REML=FALSE)
dur.RC.S.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.RC,REML=FALSE)
dur.RC.S.anova <- anova(dur.RC.S.null, dur.RC.S.model1)
print(dur.RC.S.anova)

dur.RC.S.pv<- glht(dur.RC.S.model1, mcp(Condition= "Tukey"))
summary(dur.RC.S.pv)



#GCaMP stim
dur.GC.S.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.GC,REML=FALSE)
dur.GC.S.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.GC,REML=FALSE)
dur.GC.S.anova <- anova(dur.GC.S.null, dur.GC.S.model1)
print(dur.GC.S.anova)

dur.GC.S.pv<- glht(dur.GC.S.model1, mcp(Condition= "Tukey"))
summary(dur.GC.S.pv)

##########

# for peaks shortly after the start of stimulation

stim.peaks.RC<- subset(all.peaks.RC, Condition!="0Hz")
stim.peaks.GC<- subset(all.peaks.GC, Condition!="0Hz")

# GCaMP peaks (nostim and stim) near stimulation time
all.peaks.GC<- subset(all.peaks.window, Channel=="GCaMP")

# peak time
ggplot(all.peaks.GC, aes(x=Condition,y=peakTime, fill= Condition)) +
  geom_boxplot()+
  ylab("Latency to Peak Max (s)") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


ggplot(all.peaks.RC, aes(x=Condition,y=peakTime, fill= Condition)) +
  geom_boxplot()+
  ylab("Latency to Peak Max (s)") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



# peak times
#RCaMP stim peak time
pT.RC.S.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
pT.RC.S.model1 = lmer(peakTime ~ Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
pT.RC.S.anova <- anova(pT.RC.S.null, pT.RC.S.model1)
print(pT.RC.S.anova)

pT.RC.S.Condition<- glht(pT.RC.S.model1, mcp(Condition= "Tukey"))
summary(pT.RC.S.Condition)


#GCaMP stim peak time
pT.GC.S.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
pT.GC.S.model1 = lmer(peakTime ~ Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
pT.GC.S.anova <- anova(pT.GC.S.null, pT.GC.S.model1)
print(pT.GC.S.anova)

pT.GC.S.Condition<- glht(pT.GC.S.model1, mcp(Condition= "Tukey"))
summary(pT.GC.S.Condition)

########################
# ONSET TIME
#remove entries with no onset time

all.OT<-subset(all.OT, OnsetTime!="NaN")

stim.OT.GC<-subset(all.OT, Channel=="GCaMP" & Condition !="0Hz")
stim.OT.RC<-subset(all.OT, Channel=="RCaMP" & Condition !="0Hz")

###
# neuronal responses to stimulation
NeuronalStim<-subset(all.OT, Channel=="RCaMP" & Condition!="0Hz" & OnsetTime<3)

# should have an onset time in 8 s stimulus
ggplot(NeuronalStim,aes(x=OnsetTime,y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP onset times between 0 and 8 s from stim trials")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme


# all responding neurons
Neuron95Onset<-quantile(NeuronalStim$OnsetTime, prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
NeuronPT50<-Neuron95Onset[[11]]
print(NeuronPT50)

# time thresold to consider an astrocyte to be fast:
fastTh<-NeuronPT50

###
#plot more distributions

# GCaMP onset times
ggplot(stim.OT.GC,aes(x=OnsetTime,y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP onset times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme



##
# onset time boxplots

stim.OT.GC.window<-subset(stim.OT.GC, OnsetTime<3)
stim.OT.RC.window<-subset(stim.OT.RC, OnsetTime<3)

ggplot(stim.OT.GC.window, aes(x=Condition,y=OnsetTime, fill= Condition)) +
  geom_boxplot()+
  ylab("Onset Latency (s)") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.OT.RC.window, aes(x=Condition,y=OnsetTime, fill= Condition)) +
  geom_boxplot()+
  ylab("Onset Latency (s)") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



######
## STATS

#RCaMP onset time
OT.RC.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot), stim.OT.RC.window,REML=FALSE)
OT.RC.model1 = lmer(OnsetTime ~ Condition + (1|Animal) + (1|Spot), stim.OT.RC.window,REML=FALSE)
OT.RC.anova <- anova(OT.RC.null, OT.RC.model1)
print(OT.RC.anova)

OT.RC.pv<- glht(OT.RC.model1, mcp(Condition= "Tukey"))
summary(OT.RC.pv)


#GCaMP onset time
OT.GC.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot), stim.OT.GC.window,REML=FALSE)
OT.GC.model1 = lmer(OnsetTime ~ Condition + (1|Animal) + (1|Spot), stim.OT.GC.window,REML=FALSE)
OT.GC.anova <- anova(OT.GC.null, OT.GC.model1)
print(OT.GC.anova)

OT.GC.pv<- glht(OT.GC.model1, mcp(Condition= "Tukey"))
summary(OT.GC.pv)



#########
# find fast and delayed astrocyte microdomains


# identify "FAST" astrocytes
stim.OT.GC$Group<-0
stim.OT.GC$Group[stim.OT.GC$OnsetTime<fastTh]<-"fast"
stim.OT.GC$Group[stim.OT.GC$OnsetTime>=fastTh]<-"delayed"


# add onset time information to the peak data table
stim.peaks.GC<-merge(stim.peaks.GC, stim.OT.GC[, c("ROIs_trial_Cond", "OnsetTime", "Group")], by="ROIs_trial_Cond", all.x=TRUE)


# remove peaks that don't have a corresponding onset time

stim.peaks.GC2<-stim.peaks.GC[!is.na(stim.peaks.GC$OnsetTime),]

###################

# number of active ROIs per field of view


stim.peaks.GC2$Animal_Spot<- paste(stim.peaks.GC2$Animal, stim.peaks.GC2$Spot, sep="_")
all.peaks.window$Animal_Spot<- paste(all.peaks.window$Animal, all.peaks.window$Spot, sep="_")

# number of ROIs in each trial for each field of view (across the whole trial) with a peak during no stim and stim
# only consider during the stimulus (and the same time window in no stim trials)

ROInum.8strial<-ddply(all.peaks.window, c("Animal","Spot", "Condition","Channel","Animal_Spot"), summarise, nROIs=length(unique(ROIs_trial_Cond)))

# only GCaMP stim peaks with fast or delayed
ROInum.8strial.group<-ddply(stim.peaks.GC2, c("Animal","Spot","Condition","Channel","Group"), summarise, nROIs=length(unique(ROIs_trial_Cond)))


# add in number of trials
ROInum.8strial$Ani_Spot_Cond<-paste(ROInum.8strial$Animal_Spot, ROInum.8strial$Condition, sep="_")
ROInum.8strial.group$Ani_Spot_Cond<-paste(ROInum.8strial.group$Animal, ROInum.8strial.group$Spot, ROInum.8strial.group$Condition, sep="_")

ROInum.8strial<-merge(ROInum.8strial, Spot.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)
ROInum.8strial$ROIsPerTrial<-ROInum.8strial$nROIs/ROInum.8strial$nTrials

ROInum.8strial.group<-merge(ROInum.8strial.group, Spot.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)
ROInum.8strial.group$ROIsPerTrial<-ROInum.8strial.group$nROIs/ROInum.8strial.group$nTrials

ROInum.8strial.stim<-subset(ROInum.8strial, Condition!="0Hz")
ROInum.8strial.nostim<-subset(ROInum.8strial, Condition=="0Hz")

# means
df.ROInum.8strial.RC1<-summarySE(ROInum.8strial[ROInum.8strial$Channel=="RCaMP",], measurevar = "ROIsPerTrial", groupvars = c("Condition"))

df.ROInum.8strial.GC1<-summarySE(ROInum.8strial[ROInum.8strial$Channel=="GCaMP",], measurevar = "ROIsPerTrial", groupvars = c("Condition"))
df.ROInum.8strial.GC2<-summarySE(ROInum.8strial.stim[ROInum.8strial.stim$Channel=="GCaMP",], measurevar = "ROIsPerTrial", groupvars = c("Condition"))

# plots
# RCaMP

ggplot(ROInum.8strial[ROInum.8strial$Channel=="RCaMP",], aes(x=Condition,y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


# GCaMP
ggplot(ROInum.8strial[ROInum.8strial$Channel=="GCaMP",], aes(x=Condition,y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(ROInum.8strial.stim[ROInum.8strial.stim$Channel=="GCaMP",], aes(x=Condition,y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("GCaMP stim") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme




###############
# stats for active ROI number per trials per FOV

# RCaMP stim
ROInum.RC.S.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.8strial[ROInum.8strial$Channel=="RCaMP",],REML=FALSE)
ROInum.RC.S.model1 = lmer(ROIsPerTrial ~ Condition + (1|Animal), ROInum.8strial[ROInum.8strial$Channel=="RCaMP",],REML=FALSE)
ROInum.RC.S.anova <- anova(ROInum.RC.S.null, ROInum.RC.S.model1)
print(ROInum.RC.S.anova)

ROInum.RC.S.Condition<- glht(ROInum.RC.S.model1, mcp(Condition= "Tukey"))
summary(ROInum.RC.S.Condition)

# GCaMP stim
ROInum.GC.S.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.8strial[ROInum.8strial$Channel=="GCaMP",],REML=FALSE)
ROInum.GC.S.model1 = lmer(ROIsPerTrial ~ Condition + (1|Animal), ROInum.8strial[ROInum.8strial$Channel=="GCaMP",],REML=FALSE)
ROInum.GC.S.anova <- anova(ROInum.GC.S.null, ROInum.GC.S.model1)
print(ROInum.GC.S.anova)

ROInum.GC.S.Condition<- glht(ROInum.GC.S.model1, mcp(Condition= "Tukey"))
summary(ROInum.GC.S.Condition)





###############
# GCaMP fast vs delayed

df.ROInum.8strial.group1<-summarySE(ROInum.8strial.group, measurevar = "ROIsPerTrial", groupvars = c("Condition","Group"))

# plots


ggplot(df.ROInum.8strial.group1, aes(x=Group,y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("GCaMP ROIs")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme


#boxplots
ggplot(ROInum.8strial.group, aes(x=Group,y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


# only fast GcaMP

fastROInum<- subset(ROInum.8strial.group, Group=="fast")

ggplot(fastROInum, aes(x=Condition,y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("Lck- fast ROIs per FOV") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(df.ROInum.8strial.group1[df.ROInum.8strial.group1$Group=="fast",], aes(x=Condition,y=ROIsPerTrial, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("fast GCaMP ROIs") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# only delayed GcaMP

delayedROInum<- subset(ROInum.8strial.group, Group=="delayed")

ggplot(delayedROInum, aes(x=Condition,y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("Lck- delayed ROIs per FOV") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(df.ROInum.8strial.group1[df.ROInum.8strial.group1$Group=="delayed",], aes(x=Condition,y=ROIsPerTrial, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("delayed GCaMP ROIs") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

######
# stats

# fast GCaMP only

nROI.fast.day.null = lmer(ROIsPerTrial ~ (1|Animal), fastROInum,REML=FALSE)
nROI.fast.day.model1 = lmer(ROIsPerTrial ~ Condition + (1|Animal), fastROInum,REML=FALSE)
nROI.fast.day.anova <- anova(nROI.fast.day.null, nROI.fast.day.model1)
print(nROI.fast.day.anova)

nROI.GC.Condition.fast<- glht(nROI.fast.day.model1, mcp(Condition= "Tukey"))
summary(nROI.GC.Condition.fast)



# delayed GCaMP only

nROI.delayed.day.null = lmer(ROIsPerTrial ~ (1|Animal), delayedROInum,REML=FALSE)
nROI.delayed.day.model1 = lmer(ROIsPerTrial ~ Condition + (1|Animal), delayedROInum,REML=FALSE)
nROI.delayed.day.anova <- anova(nROI.delayed.day.null, nROI.delayed.day.model1)
print(nROI.delayed.day.anova)

nROI.GC.Condition.delayed<- glht(nROI.delayed.day.model1, mcp(Condition= "Tukey"))
summary(nROI.GC.Condition.delayed)



# wilcox test for median compairsons

wilcox.test(fastROInum$ROIsPerTrial[fastROInum$Condition=="10Hz"], 
            fastROInum$ROIsPerTrial[fastROInum$Condition=="90Hz"])


######
# compare amplitude and duration of fast signals

#########
#mean Duration of signals for each ROI

#GCaMP
df.dur2.GC<- summarySE(stim.peaks.GC2, measurevar = "Duration", groupvars = c("Condition","Group"))
df.dur3.GC<- summarySE(stim.peaks.GC2, measurevar = "Duration", groupvars = c("Condition"))


# graphs
ggplot(data=df.dur2.GC, aes(x=Group, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# gCaMP graphs
ggplot(data=df.dur3.GC, aes(x=Condition, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("GCaMP- responding peaks") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


#amplitude
df.amp2.GC<- summarySE(stim.peaks.GC2, measurevar = "amplitude", groupvars = c("Condition","Group"))
df.amp3.GC<- summarySE(stim.peaks.GC2, measurevar = "amplitude", groupvars = c("Condition"))


# graphs
ggplot(data=df.amp2.GC, aes(x=Group, y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# gCaMP graphs
ggplot(data=df.amp3.GC, aes(x=Condition, y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("GCaMP- responding peaks") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


#####################
## STATS

Cond_group<-interaction(stim.peaks.GC2$Condition, stim.peaks.GC2$Group)

#GCaMP stim
dur2.GC.S.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC2,REML=FALSE)
dur2.GC.S.model1 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC2,REML=FALSE)
dur2.GC.S.model2 = lmer(Duration ~ Cond_group + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC2,REML=FALSE)
dur2.GC.S.anova <- anova(dur2.GC.S.null, dur2.GC.S.model1,dur2.GC.S.model2)
print(dur2.GC.S.anova)

dur2.GC.S.pv<- glht(dur2.GC.S.model1, mcp(Condition= "Tukey"))
summary(dur2.GC.S.pv)

dur2.GC.S.pv2<- glht(dur2.GC.S.model2, mcp(Cond_group= "Tukey"))
summary(dur2.GC.S.pv2)

