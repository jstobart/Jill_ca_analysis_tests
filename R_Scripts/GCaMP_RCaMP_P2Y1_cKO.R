
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

lck <- read.table("G:/ZurichData/Astrocyte_Calcium/P2Y1_Mice/Results/lckGC_P2Y1_timepoints_peaks_05_2018.csv", header=TRUE, sep = ",")
cyto <- read.table("G:/ZurichData/Astrocyte_Calcium/P2Y1_Mice/Results/cytoGC_P2Y1_timepoints_peaks_05_2018.csv", header=TRUE, sep = ",")

lck <- read.table("H:/P2Y1_Data/Results/lckGC_P2Y1_timepoints_peaks_05_2018.csv", header=TRUE, sep = ",")
cyto <- read.table("H:/P2Y1_Data/Results/cytoGC_P2Y1_timepoints_peaks_05_2018.csv", header=TRUE, sep = ",")

##########

astro.peaks<-rbind(lck, cyto)

astro.peaks$Spot_trial_TP<-paste(astro.peaks$Spot, astro.peaks$Trial, astro.peaks$TimePoint, sep="_")

astro.peaks$Spot_trial_TP_Cond<-paste(astro.peaks$Spot_trial_TP,astro.peaks$Condition, sep="_")

astro.peaks$Spot_trial_TP_Cond_ROI<-paste(astro.peaks$Spot_trial_TP_Cond,astro.peaks$roiName, sep="_")

astro.peaks$Spot_trial_TP_ROI<-paste(astro.peaks$Spot_trial_TP,astro.peaks$roiName, sep="_")

astro.peaks$Spot_ROI<-paste(astro.peaks$Spot,astro.peaks$roiName, sep="_")

# ROITypes
astro.peaks$ROIType[grepl("r",astro.peaks$roiName)]="Process"
astro.peaks$ROIType[grepl("E",astro.peaks$roiName)]="Endfoot"
astro.peaks$ROIType[grepl("S",astro.peaks$roiName)]="Soma"

astro.peaks$ROIType<- as.factor(astro.peaks$ROIType)

#Time Point Type
astro.peaks$TimeGroup[grepl("B",astro.peaks$TimePoint)]="before"
astro.peaks$TimeGroup[grepl("T",astro.peaks$TimePoint)]="after"


astro.peaks$TimeGroup<- as.factor(astro.peaks$TimeGroup)
astro.peaks$TimeGroup<- factor(astro.peaks$TimeGroup, levels = c("before","after"))

# remove matching astrocyte process and soma ROIs
Overlap= astro.peaks$overlap!=0
astro.peaks<-astro.peaks[!Overlap,]

# adjust peak time and duration
astro.peaks$peakTime<- astro.peaks$peakTime-2.11
astro.peaks$peakStart<- astro.peaks$peakStart-2.11
astro.peaks$peakStartHalf<- astro.peaks$peakStartHalf-2.11
astro.peaks$Duration<- astro.peaks$halfWidth*2

# drop peaks that occur before the start of stimulation
astro.peaks<-subset(astro.peaks,peakTime>0)

######
#distributions

# amplitude
#Lck
ggplot(astro.peaks[astro.peaks$Channel=="LckGCaMP",], aes(x=Condition,y=amplitude, fill= TimeGroup)) +
  geom_boxplot(notch=TRUE)+
  ylab("Lck amplitude") +
  ggtitle("notched")+
  max.theme

#cyto
ggplot(astro.peaks[astro.peaks$Channel=="GCaMP",], aes(x=Condition,y=amplitude, fill= TimeGroup)) +
  geom_boxplot(notch=TRUE)+
  ylab("cyto amplitude") +
  ggtitle("notched")+
  max.theme


# duration
#Lck
ggplot(astro.peaks[astro.peaks$Channel=="LckGCaMP",], aes(x=Condition,y=Duration, fill= TimeGroup)) +
  geom_boxplot(notch=TRUE)+
  ylab("Lck Duration") +
  ggtitle("notched")+
  max.theme

#cyto
ggplot(astro.peaks[astro.peaks$Channel=="GCaMP",], aes(x=Condition,y=Duration, fill= TimeGroup)) +
  geom_boxplot(notch=TRUE)+
  ylab("cyto Duration") +
  ggtitle("notched")+
  max.theme


######
#means of all peaks

# time groups
df.amp.lck<-summarySE(astro.peaks[astro.peaks$Channel=="LckGCaMP",], measurevar = "amplitude", groupvars = c("Condition", "TimeGroup"))
df.amp.cyto<-summarySE(astro.peaks[astro.peaks$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("Condition", "TimeGroup"))

ggplot(df.amp.lck, aes(x=Condition,y=amplitude, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean amplitude") +
  max.theme

ggplot(df.amp.lck[df.amp.lck$Condition=="nostim",], aes(x=Condition,y=amplitude, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean amplitude") +
  max.theme

ggplot(df.amp.cyto, aes(x=Condition,y=amplitude, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean amplitude") +
  max.theme

#stats
# Likelihood-ratio test amplitude

# only LCK, no stim
amp.GC.resp.null = lmer(amplitude ~ (1|Spot) + (1|Spot_ROI), astro.peaks[(astro.peaks$Channel=="LckGCaMP" & astro.peaks$Condition=="nostim"),],REML=FALSE)
amp.GC.resp.model2 = lmer(amplitude ~ TimeGroup + (1|Spot) + (1|Spot_ROI), astro.peaks[(astro.peaks$Channel=="LckGCaMP" & astro.peaks$Condition=="nostim"),],REML=FALSE)
amp.GC.resp.anova <- anova(amp.GC.resp.null,amp.GC.resp.model2)
print(amp.GC.resp.anova)

# p values
amp.GC.resp <- lsmeans(amp.GC.resp.model2, pairwise ~ TimeGroup, glhargs=list())
summary(amp.GC.resp)

# time points
df.amp.lck2<-summarySE(astro.peaks[astro.peaks$Channel=="LckGCaMP",], measurevar = "amplitude", groupvars = c("Condition", "TimePoint"))
df.amp.cyto2<-summarySE(astro.peaks[astro.peaks$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("Condition", "TimePoint"))

ggplot(df.amp.lck2, aes(x=Condition,y=amplitude, fill= TimePoint)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean amplitude") +
  max.theme

ggplot(df.amp.cyto2, aes(x=Condition,y=amplitude, fill= TimePoint)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean amplitude") +
  max.theme

##########
# only responding peaks

astro.peaks.resp<-subset(astro.peaks, peakTime<20 & Condition=="WPstim")

#means of stim peaks

# time groups
df.amp.lck3<-summarySE(astro.peaks.resp[astro.peaks.resp$Channel=="LckGCaMP",], measurevar = "amplitude", groupvars = c("TimeGroup"))
df.amp.cyto3<-summarySE(astro.peaks.resp[astro.peaks.resp$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("TimeGroup"))

ggplot(df.amp.lck3, aes(x=TimeGroup,y=amplitude, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean amplitude") +
  max.theme

ggplot(df.amp.cyto3, aes(x=TimeGroup,y=amplitude, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean amplitude") +
  max.theme

#####
#means duration

df.dur.lck<-summarySE(astro.peaks[astro.peaks$Channel=="LckGCaMP",], measurevar = "Duration", groupvars = c("Condition", "TimeGroup"))
df.dur.cyto<-summarySE(astro.peaks[astro.peaks$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Condition", "TimeGroup"))

ggplot(df.dur.lck, aes(x=Condition,y=Duration, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Duration") +
  max.theme

ggplot(df.amp.cyto, aes(x=Condition,y=amplitude, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean amplitude") +
  max.theme


df.dur.lck2<-summarySE(astro.peaks[astro.peaks$Channel=="LckGCaMP",], measurevar = "Duration", groupvars = c("Condition", "TimePoint"))
df.dur.cyto2<-summarySE(astro.peaks[astro.peaks$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Condition", "TimePoint"))

ggplot(df.dur.lck2, aes(x=Condition,y=Duration, fill= TimePoint)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Duration") +
  max.theme

ggplot(df.dur.cyto2, aes(x=Condition,y=Duration, fill= TimePoint)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Duration") +
  max.theme



#######
# signals per field of view

# count the number of peaks per field of view for each spot at each time point

Signals.trial.spot<- ddply(astro.peaks, c("Animal", "Spot","Channel", "TimePoint", "TimeGroup", "Condition","Trial","Spot_trial_TP"), summarise, 
                          nPeaks = length(amplitude))

Signals.spot<-ddply(Signals.trial.spot, c("Animal", "Spot", "Channel", "TimePoint", "TimeGroup", "Condition"), summarise, 
                    meanPeaks = mean(nPeaks))

# mean signals
df.lck.sig.spot<-summarySE(Signals.spot[Signals.spot$Channel=="LckGCaMP",], measurevar = "meanPeaks", groupvars = c("Condition", "TimeGroup"))

ggplot(df.lck.sig.spot, aes(x=Condition,y=meanPeaks, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanPeaks-se, ymax=meanPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peaks per FOV") +
  max.theme

######
# number of ROIs in each trial for each field of view (across the whole trial)

ROInum.trial.spot<- ddply(astro.peaks, c("Animal", "Spot","Channel", "TimePoint", "TimeGroup", "Condition","Trial","Spot_trial_TP"), summarise, 
                           nROIs = length(unique(roiName)))

ROInum.spot<-ddply(ROInum.trial.spot, c("Animal", "Spot", "Channel", "TimePoint", "TimeGroup", "Condition"), summarise, 
                    meanROIs = mean(nROIs))

# mean
df.lck.ROInum1<-summarySE(ROInum.spot[ROInum.spot$Channel=="LckGCaMP",], measurevar = "meanROIs", groupvars = c("Condition","TimeGroup"))
df.lck.ROInum2<-summarySE(ROInum.spot[ROInum.spot$Channel=="LckGCaMP",], measurevar = "meanROIs", groupvars = c("Condition","TimePoint"))

ggplot(df.lck.ROInum1, aes(x=Condition,y=meanROIs, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanROIs-se, ymax=meanROIs+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIs/FOV") +
  max.theme

ggplot(df.lck.ROInum1[df.lck.ROInum1$Condition=="WPstim",], aes(x=Condition,y=meanROIs, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanROIs-se, ymax=meanROIs+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIs/FOV") +
  max.theme

# only LCK, WP stim
ROInum.GC.null = lmer(meanROIs ~ (1|Spot), ROInum.spot[(ROInum.spot$Channel=="LckGCaMP" & ROInum.spot$Condition=="WPstim"),],REML=FALSE)
ROInum.GC.model2 = lmer(meanROIs ~ TimeGroup + (1|Spot), ROInum.spot[(ROInum.spot$Channel=="LckGCaMP" & ROInum.spot$Condition=="WPstim"),],REML=FALSE)
ROInum.GC.anova <- anova(ROInum.GC.null,ROInum.GC.model2)
print(ROInum.GC.anova)

# p values
ROInum.GC.pv <- lsmeans(ROInum.GC.model2, pairwise ~ TimeGroup, glhargs=list())
summary(ROInum.GC.pv)

###############
#RCaMP

RCaMP <- read.table("G:/ZurichData/Astrocyte_Calcium/P2Y1_Mice/Results/RC_P2Y1_timepoints_peaks_05_2018.csv", header=TRUE, sep = ",")
RCaMP <- read.table("H:/P2Y1_Data/Results/RC_P2Y1_timepoints_peaks_05_2018.csv", header=TRUE, sep = ",")


RCaMP$Spot_trial_TP<-paste(RCaMP$Spot, RCaMP$Trial, RCaMP$TimePoint, sep="_")

RCaMP$Spot_ROI<-paste(RCaMP$Spot, RCaMP$roiName, sep="_")

RCaMP$Spot_trial_TP_Cond<-paste(RCaMP$Spot_trial_TP,RCaMP$Condition, sep="_")

RCaMP$Spot_trial_TP_Cond_ROI<-paste(RCaMP$Spot_trial_TP_Cond, RCaMP$roiName, sep="_")

RCaMP$Spot_trial_TP_ROI<-paste(RCaMP$Spot_trial_TP,RCaMP$roiName, sep="_")

RCaMP$Spot_TP_ROI<-paste(RCaMP$Spot, RCaMP$TimePoint, RCaMP$roiName, sep="_")


# ROITypes
RCaMP$ROIType[grepl("N",RCaMP$roiName)]="Neuron"

RCaMP$ROIType<- as.factor(RCaMP$ROIType)

#Time Point Type
RCaMP$TimeGroup[grepl("B",RCaMP$TimePoint)]="before"
RCaMP$TimeGroup[grepl("T",RCaMP$TimePoint)]="after"


RCaMP$TimeGroup<- as.factor(RCaMP$TimeGroup)
RCaMP$TimeGroup<- factor(RCaMP$TimeGroup, levels = c("before","after"))


# adjust peak time and duration
RCaMP$peakTime<- RCaMP$peakTime-2.11
RCaMP$peakStart<- RCaMP$peakStart-2.11
RCaMP$peakStartHalf<- RCaMP$peakStartHalf-2.11
RCaMP$Duration<- RCaMP$halfWidth*2

# drop peaks that occur before the start of stimulation
#RCaMP<-subset(RCaMP,peakTime>0)

# add in zeros for NaNs  where there was no peak detected
test3 = is.na(RCaMP$prominence)
RCaMP$peakAUC[test3] = 0
RCaMP$prominence[test3] = 0
RCaMP$amplitude[test3] = 0
RCaMP$Duration[test3] = 0
RCaMP$numPeaks[test3] = 0


######
#distributions

# amplitude
ggplot(RCaMP, aes(x=Condition,y=amplitude, fill= TimeGroup)) +
  geom_boxplot(notch=TRUE)+
  ylab("RCaMP amplitude") +
  ggtitle("notched")+
  max.theme


# duration
ggplot(RCaMP, aes(x=Condition,y=Duration, fill= TimeGroup)) +
  geom_boxplot(notch=TRUE)+
  ylab("RCaMP Duration") +
  ggtitle("notched")+
  max.theme

# peak time
ggplot(RCaMP, aes(x=Condition,y=peakTime, fill= TimeGroup)) +
  geom_boxplot(notch=TRUE)+
  ylab("RCaMP peak time (s)") +
  ggtitle("notched")+
  max.theme

###########

# mean values per neuron across all trials from each time point

# sum number of ROIs per field of view
RCaMP.ROIs<-ddply(RCaMP, c("Animal","Spot","Condition","TimePoint","TimeGroup", "Spot_ROI","Spot_TP_ROI"), summarise, 
                  nPeaks=sum(numPeaks), nTrials=length(unique(Trial)),
                  meanAmp=mean(amplitude), meanDur=mean(Duration),
                  meanpeakAUC=mean(peakAUC))

RCaMP.ROIs$Freq<-RCaMP.ROIs$nPeaks/RCaMP.ROIs$nTrials

# sum number of ROIs per field of view
respRCaMP<-subset(RCaMP, peakTime<10 & Condition=="WPstim")
respNeurons<-unique(respRCaMP$Spot_TP_ROI)
respSpotROI<-unique(respRCaMP$Spot_ROI)

RCaMP.respPeaks<-subset(RCaMP, Spot_ROI %in% respSpotROI)

RCaMP.respROIs<-ddply(RCaMP.respPeaks, c("Animal","Spot","Condition","TimePoint","TimeGroup", "Spot_ROI","Spot_TP_ROI"), summarise, 
                  nPeaks=sum(numPeaks), nTrials=length(unique(Trial)),
                  meanAmp=mean(amplitude), meanDur=mean(Duration),
                  meanpeakAUC=mean(peakAUC))

RCaMP.respROIs$Freq<-RCaMP.respROIs$nPeaks/RCaMP.respROIs$nTrials

# responding ROIs

###########
#distributions

# amplitude
ggplot(RCaMP.respROIs, aes(x=Condition,y=meanAmp, fill= TimeGroup)) +
  geom_boxplot(notch=TRUE)+
  ylab("RCaMP ROIs amplitude") +
  ggtitle("notched")+
  max.theme


# duration
ggplot(RCaMP.respROIs, aes(x=Condition,y=meanDur, fill= TimeGroup)) +
  geom_boxplot(notch=TRUE)+
  ylab("RCaMP ROIs Duration") +
  ggtitle("notched")+
  max.theme

# frequency (signals/trial)
ggplot(RCaMP.respROIs, aes(x=Condition,y=Freq, fill= TimeGroup)) +
  geom_boxplot(notch=TRUE)+
  ylab("signals/trial") +
  ggtitle("notched")+
  max.theme

ggplot(RCaMP.respROIs, aes(x=Condition,y=Freq, fill= TimePoint)) +
  geom_boxplot(notch=TRUE)+
  ylab("signals/trial") +
  ggtitle("notched")+
  max.theme

###########
#means

df.amp1.RC<-summarySE(RCaMP.respROIs, measurevar = "meanAmp", groupvars = c("Condition", "TimeGroup"))
df.amp2.RC<-summarySE(RCaMP.respROIs, measurevar = "meanAmp", groupvars = c("Condition", "TimePoint"))
df.amp3.RC<-summarySE(RCaMP.ROIs, measurevar = "meanAmp", groupvars = c("Condition", "TimeGroup"))
df.amp4.RC<-summarySE(RCaMP.ROIs, measurevar = "meanAmp", groupvars = c("Condition", "TimePoint"))

ggplot(df.amp1.RC, aes(x=Condition,y=meanAmp, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanAmp-se, ymax=meanAmp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean ROI amplitude") +
  max.theme

ggplot(df.amp2.RC, aes(x=Condition,y=meanAmp, fill= TimePoint)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanAmp-se, ymax=meanAmp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean ROI amplitude") +
  max.theme

ggplot(df.amp3.RC, aes(x=Condition,y=meanAmp, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanAmp-se, ymax=meanAmp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean ROI amplitude") +
  max.theme

ggplot(df.amp3.RC[df.amp3.RC$Condition=="WPstim",], aes(x=Condition,y=meanAmp, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanAmp-se, ymax=meanAmp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean ROI amplitude") +
  max.theme

#means

df.freq1.RC<-summarySE(RCaMP.respROIs, measurevar = "Freq", groupvars = c("Condition", "TimeGroup"))
df.freq2.RC<-summarySE(RCaMP.respROIs, measurevar = "Freq", groupvars = c("Condition", "TimePoint"))
df.freq3.RC<-summarySE(RCaMP.ROIs, measurevar = "Freq", groupvars = c("Condition", "TimeGroup"))

ggplot(df.freq1.RC, aes(x=Condition,y=Freq, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Freq-se, ymax=Freq+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean signals/trial") +
  max.theme

ggplot(df.freq2.RC, aes(x=Condition,y=Freq, fill= TimePoint)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Freq-se, ymax=Freq+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean signals/trial") +
  max.theme

ggplot(df.freq3.RC, aes(x=Condition,y=Freq, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Freq-se, ymax=Freq+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean signals/trial") +
  max.theme

#############
#stats
# Likelihood-ratio test amplitude
amp.RC.resp.null = lmer(meanAmp ~ (1|Spot) + (1|Spot_ROI), RCaMP.respROIs,REML=FALSE)
amp.RC.resp.model1 = lmer(meanAmp ~ Condition + (1|Spot) + (1|Spot_ROI), RCaMP.respROIs,REML=FALSE)
amp.RC.resp.model2 = lmer(meanAmp ~ TimeGroup + (1|Spot) + (1|Spot_ROI), RCaMP.respROIs,REML=FALSE)
amp.RC.resp.model3 = lmer(meanAmp ~ Condition*TimeGroup + (1|Spot) + (1|Spot_ROI), RCaMP.respROIs,REML=FALSE)
amp.RC.resp.anova <- anova(amp.RC.resp.null, amp.RC.resp.model1, amp.RC.resp.model2, amp.RC.resp.model3)
print(amp.RC.resp.anova)

# p values
amp.RC.resp <- lsmeans(amp.RC.resp.model3, pairwise ~ Condition*TimeGroup, glhargs=list())
summary(amp.RC.resp)


# Likelihood-ratio test amplitude
freq.RC.resp.null = lmer(Freq ~ (1|Spot) + (1|Spot_ROI), RCaMP.respROIs,REML=FALSE)
freq.RC.resp.model1 = lmer(Freq ~ Condition + (1|Spot) + (1|Spot_ROI), RCaMP.respROIs,REML=FALSE)
freq.RC.resp.model2 = lmer(Freq ~ TimeGroup + (1|Spot) + (1|Spot_ROI), RCaMP.respROIs,REML=FALSE)
freq.RC.resp.model3 = lmer(Freq ~ Condition*TimeGroup + (1|Spot) + (1|Spot_ROI), RCaMP.respROIs,REML=FALSE)
freq.RC.resp.anova <- anova(freq.RC.resp.null, freq.RC.resp.model1, freq.RC.resp.model2, freq.RC.resp.model3)
print(freq.RC.resp.anova)

# p values
freq.RC.resp <- lsmeans(freq.RC.resp.model3, pairwise ~ Condition*TimeGroup, glhargs=list())
summary(freq.RC.resp)


###########
# fraction response (active trials)
RCaMP$ActivePeak <- 0
farpeaks1 <- RCaMP$peakTime>=0 & RCaMP$peakTime<=10
RCaMP$ActivePeak[farpeaks1] <- 1 


# count number of peaks per trial in the first 30 sec
fraction.response<- ddply(RCaMP, c("Animal", "Spot", "TimePoint", "TimeGroup", "Condition","roiName"), summarise, 
                      nActive = sum(ActivePeak), NTrials = length(unique(Trial)))
                      
fraction.response$frac.resp <- fraction.response$nActive/fraction.response$NTrials

#histogram of responses
ggplot(fraction.response[fraction.response$TimeGroup=="before",], aes(x=frac.resp, fill=Condition)) + geom_histogram(aes(y=..count../sum(..count..)), binwidth=0.2,position="dodge")+
  ggtitle("Fraction of active trials before")

ggplot(fraction.response[fraction.response$TimeGroup=="1month",], aes(x=frac.resp, fill=Condition)) + geom_histogram(aes(y=..count../sum(..count..)), binwidth=0.2,position="dodge")+
  ggtitle("Fraction of active trials 1 month after TAM")

ggplot(fraction.response[fraction.response$TimeGroup=="2months",], aes(x=frac.resp, fill=Condition)) + geom_histogram(aes(y=..count../sum(..count..)), binwidth=0.2,position="dodge")+
  ggtitle("Fraction of active trials 2 months after TAM")

# mean percent response- i.e. the fraction of trials where a peak is detected in the first 30 sec
df.fracResp1.RC <- summarySE(fraction.response, measurevar="frac.resp", groupvars=c("Condition", "TimeGroup"))
df.fracResp2.RC <- summarySE(fraction.response, measurevar="frac.resp", groupvars=c("Condition","TimePoint"))


ggplot(df.fracResp1.RC, aes(x=Condition,y=frac.resp, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=frac.resp-se, ymax=frac.resp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("frac.resp") +
  max.theme

ggplot(df.fracResp2.RC, aes(x=Condition,y=frac.resp, fill= TimePoint)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=frac.resp-se, ymax=frac.resp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("frac.resp") +
  max.theme

###################

# change in amplitude or frequency over time

#add frac.resp to ROI info
fraction.response$Spot_TP_ROI_Cond<-paste(fraction.response$Spot, fraction.response$TimePoint, fraction.response$roiName,fraction.response$Condition, sep="_")
RCaMP.ROIs$Spot_TP_ROI_Cond<-paste(RCaMP.ROIs$Spot_TP_ROI,RCaMP.ROIs$Condition, sep="_")

RCaMP.ROIs<-merge(RCaMP.ROIs, fraction.response[,c("Spot_TP_ROI_Cond","frac.resp")], by ="Spot_TP_ROI_Cond", all.x=FALSE)



# find each ROI

#amplitude table
amplitude<-RCaMP.ROIs[,c("Spot_ROI","TimePoint","Condition","meanAmp")]

amp_wide <- dcast(amplitude, Spot_ROI + Condition ~ TimePoint)

amp_wide$baselineMean<- rowMeans(amp_wide[c('DayBL01', 'DayBL02')], na.rm=TRUE)
amp_wide$BL01_change<-amp_wide$DayBL01/amp_wide$baselineMean 
amp_wide$BL02_change<-amp_wide$DayBL02/amp_wide$baselineMean 
amp_wide$TM10_change<-amp_wide$DayTM10/amp_wide$baselineMean 
amp_wide$TM28_change<-amp_wide$DayTM28/amp_wide$baselineMean 
amp_wide$TM42_change<-amp_wide$DayTM42/amp_wide$baselineMean 
amp_wide$TM80_change<-amp_wide$DayTM80/amp_wide$baselineMean 

# convert Inf to 10
amp_wide[mapply(is.infinite, amp_wide)] <- 10
# convert NaN to 0
#amp_wide[mapply(is.na, amp_wide)] <- 0


#plot lines for each ROI?
#mean of each time point?



frequency<-RCaMP.ROIs[,c("Spot_ROI","TimePoint","Condition","Freq")]
duration<-RCaMP.ROIs[,c("Spot_ROI","TimePoint","Condition","meanDur")]
frac.resp.df<-RCaMP.ROIs[,c("Spot_ROI","TimePoint","Condition","frac.resp")]
peakAUC<-RCaMP.ROIs[,c("Spot_ROI","TimePoint","Condition","meanpeakAUC")]



summary.freq <- summarySE(frequency, measurevar="Freq", groupvars=c("Condition","TimePoint"))


ggplot(df.fracResp1.RC, aes(x=Condition,y=frac.resp, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=frac.resp-se, ymax=frac.resp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("frac.resp") +
  max.theme

#############
#stats
# Likelihood-ratio test frequency
frac.resp.null = lmer(frac.resp ~ (1|Animal) + (1|Spot) + (1|ROI), fraction.response,REML=FALSE)
frac.resp.model1 = lmer(frac.resp ~ Condition + (1|Animal) + (1|Spot) + (1|ROI), fraction.response,REML=FALSE)
frac.resp.model2 = lmer(frac.resp~ Condition*ROIType + (1|Animal) + (1|Spot) + (1|ROI), fraction.response,REML=FALSE)
frac.resp.model3 = lmer(frac.resp~ Condition*Layer + (1|Animal) + (1|Spot) + (1|ROI), fraction.response,REML=FALSE)
frac.resp.model4 = lmer(frac.resp~ Condition*Layer*ROIType + (1|Animal) + (1|Spot) + (1|ROI), fraction.response,REML=FALSE)
frac.resp.anova <- anova(frac.resp.null, frac.resp.model1,frac.resp.model2,frac.resp.model3,frac.resp.model4)
print(frac.resp.anova)
# p values
frac.resp.pv.stim <- lsmeans(frac.resp.model1, pairwise ~ Condition, glhargs=list())
summary(frac.resp.pv.stim)
frac.resp.pv.stim_type <- lsmeans(frac.resp.model2, pairwise ~ Condition*ROIType, glhargs=list())
summary(frac.resp.pv.stim_type)
frac.resp.pv.stim_layer <- lsmeans(frac.resp.model3, pairwise ~ Condition*Layer, glhargs=list())
summary(frac.resp.pv.stim_layer)
frac.resp.pv.stim_layer_type <- lsmeans(frac.resp.model4, pairwise ~ Condition*Layer*ROIType, glhargs=list())
summary(frac.resp.pv.stim_layer_type)




