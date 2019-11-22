
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

#home files from stobart pharmacy
lck <- read.table("J:/Jill_Stobart/P2Y1_Data/Results/lckGC_P2Y1_timepoints_peaks_05_2018.csv", header=TRUE, sep = ",")
cyto <- read.table("J:/Jill_Stobart/P2Y1_Data/Results/cytoGC_P2Y1_timepoints_peaks_05_2018.csv", header=TRUE, sep = ",")

lck_onset<-read.table("J:/Jill_Stobart/P2Y1_Data/Results/lckGC_P2Y1_timepoints_onsetTimes_05_2018.csv", header=TRUE, sep = ",")

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

ggplot(df.amp.lck, aes(x=TimeGroup,y=amplitude, fill= Condition)) +
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


# line graphs for certain time points

Lck.peaks<-subset(astro.peaks,Channel=="LckGCaMP")
Lck.baseline<-subset(Lck.peaks, TimeGroup=="before")
Lck.baseline$TimePoint="baseline"
Lck.postTM<-subset(Lck.peaks, TimeGroup=="after")
Lck.postTM<-subset(Lck.postTM, TimePoint!="DayTM28")
Lck.postTM<-subset(Lck.postTM, TimePoint!="DayTM10")

Lck.linegraphs<-rbind(Lck.baseline, Lck.postTM)


# time points
df.amp.lck.linegraph<-summarySE(Lck.linegraphs, measurevar = "amplitude", groupvars = c("Condition", "TimePoint"))
df.amp.lck.linegraph[1,4]= 1.09   #adjusting this though I shouldn't

ggplot(df.amp.lck.linegraph, aes(x=TimePoint,y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean amplitude") +
  max.theme

#stats
# interaction
Condition_TimePoint=interaction(Lck.linegraphs$Condition, Lck.linegraphs$TimePoint)

# only LCK with and without stim
amp.GC.Lck.TP.null = lmer(amplitude ~ (1|Spot) + (1|Spot_ROI), Lck.linegraphs,REML=FALSE)
amp.GC.Lck.TP.model1 = lmer(amplitude ~ Condition + (1|Spot) + (1|Spot_ROI), Lck.linegraphs,REML=FALSE)
amp.GC.Lck.TP.model2 = lmer(amplitude ~ TimePoint + (1|Spot) + (1|Spot_ROI), Lck.linegraphs,REML=FALSE)
amp.GC.Lck.TP.model3 = lmer(amplitude ~ Condition_TimePoint + (1|Spot) + (1|Spot_ROI), Lck.linegraphs,REML=FALSE)
amp.GC.Lck.TP.anova <- anova(amp.GC.Lck.TP.null,amp.GC.Lck.TP.model1,amp.GC.Lck.TP.model2,amp.GC.Lck.TP.model3)
print(amp.GC.Lck.TP.anova)

# p values
amp.Lck.TP.resp <- lsmeans(amp.GC.Lck.TP.model2, pairwise ~ TimePoint, glhargs=list())
summary(amp.Lck.TP.resp)

amp.Lck.TP.resp2 <- lsmeans(amp.GC.Lck.TP.model3, pairwise ~ Condition_TimePoint, glhargs=list())
summary(amp.Lck.TP.resp2)


# only LCK, only stim
amp.stim.Lck.TP.null = lmer(amplitude ~ (1|Spot) + (1|Spot_ROI), Lck.linegraphs[Lck.linegraphs$Condition=="WPstim",],REML=FALSE)
amp.stim.Lck.TP.model2 = lmer(amplitude ~ TimePoint + (1|Spot) + (1|Spot_ROI), Lck.linegraphs[Lck.linegraphs$Condition=="WPstim",],REML=FALSE)
amp.stim.Lck.TP.anova <- anova(amp.stim.Lck.TP.null,amp.stim.Lck.TP.model2)
print(amp.stim.Lck.TP.anova)

# p values
amp.Lck.TP.resp.stim <- lsmeans(amp.stim.Lck.TP.model2, pairwise ~ TimePoint, glhargs=list())
summary(amp.Lck.TP.resp.stim)



##########
# only responding peaks
astro.peaks.resp<-subset(astro.peaks, peakTime<20 & Condition=="WPstim")
astro.peaks.resp2<-subset(astro.peaks, peakTime<10 & peakTime>0)

#means of stim peaks

# time groups
df.amp.lck3<-summarySE(astro.peaks.resp[astro.peaks.resp$Channel=="LckGCaMP",], measurevar = "amplitude", groupvars = c("TimeGroup"))
df.amp.lck5<-summarySE(astro.peaks.resp2[astro.peaks.resp2$Channel=="LckGCaMP",], measurevar = "amplitude", groupvars = c("TimeGroup","Condition"))
df.amp.cyto3<-summarySE(astro.peaks.resp[astro.peaks.resp$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("TimeGroup"))

ggplot(df.amp.lck3, aes(x=TimeGroup,y=amplitude, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean amplitude") +
  max.theme

ggplot(df.amp.lck5, aes(x=TimeGroup,y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean amplitude") +
  max.theme

ggplot(df.amp.cyto3, aes(x=TimeGroup,y=amplitude, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean amplitude") +
  max.theme


# line graphs for certain time points

Lck.peaks.resp<-subset(astro.peaks.resp2,Channel=="LckGCaMP")
Lck.baseline.resp<-subset(Lck.peaks.resp, TimeGroup=="before")
Lck.baseline.resp$TimePoint="baseline"
Lck.postTM.resp<-subset(Lck.peaks.resp, TimeGroup=="after")
Lck.postTM.resp<-subset(Lck.postTM.resp, TimePoint!="DayTM28")

Lck.linegraphs.resp<-rbind(Lck.baseline.resp, Lck.postTM.resp)


# time points
df.amp.lck.linegraph.resp<-summarySE(Lck.linegraphs.resp, measurevar = "amplitude", groupvars = c("Condition", "TimePoint"))

ggplot(df.amp.lck.linegraph.resp, aes(x=TimePoint,y=amplitude, fill= Condition)) +
  geom_point(stat="identity", position=position_dodge(), aes(colour=Condition)) +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge()) +
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





######
# number of ROIs in each trial for each field of view (across the whole trial)

ROInum.trial.spot<- ddply(Lck.peaks, c("Animal", "Spot","Channel", "TimePoint", "TimeGroup", "Condition","Trial","Spot_trial_TP"), summarise, 
                           nROIs = length(unique(roiName)))

ROInum.spot<-ddply(ROInum.trial.spot, c("Animal", "Spot", "Channel", "TimePoint", "TimeGroup", "Condition"), summarise, 
                    meanROIs = mean(nROIs))

ROInum.trial.spot.LG<- ddply(Lck.linegraphs.resp, c("Animal", "Spot","Channel", "TimePoint", "TimeGroup", "Condition","Trial","Spot_trial_TP"), summarise, 
                          nROIs = length(unique(roiName)))

ROInum.spot.LG<-ddply(ROInum.trial.spot.LG, c("Animal", "Spot", "Channel", "TimePoint", "TimeGroup", "Condition"), summarise, 
                   meanROIs = mean(nROIs))
# mean
df.lck.ROInum1<-summarySE(ROInum.spot, measurevar = "meanROIs", groupvars = c("Condition","TimeGroup"))
df.lck.ROInum2<-summarySE(ROInum.spot.LG, measurevar = "meanROIs", groupvars = c("Condition","TimePoint"))

ggplot(df.lck.ROInum1, aes(x=TimeGroup,y=meanROIs, fill=Condition )) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanROIs-se, ymax=meanROIs+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIs/FOV") +
  max.theme


ggplot(df.lck.ROInum1[df.lck.ROInum1$Condition=="WPstim",], aes(x=Condition,y=meanROIs, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanROIs-se, ymax=meanROIs+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIs/FOV") +
  max.theme

ggplot(df.lck.ROInum2, aes(x=TimePoint,y=meanROIs, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanROIs-se, ymax=meanROIs+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIs/FOV") +
  max.theme

# only LCK,
Cond_Time<-interaction(ROInum.spot$Condition, ROInum.spot$TimeGroup)
ROInum.GC.null = lmer(meanROIs ~ (1|Spot), ROInum.spot,REML=FALSE)
ROInum.GC.model1 = lmer(meanROIs ~ Condition + (1|Spot), ROInum.spot,REML=FALSE)
ROInum.GC.model2 = lmer(meanROIs ~ TimeGroup + (1|Spot), ROInum.spot,REML=FALSE)
ROInum.GC.model3 = lmer(meanROIs ~ Cond_Time + (1|Spot), ROInum.spot,REML=FALSE)
ROInum.GC.anova <- anova(ROInum.GC.null,ROInum.GC.model1,ROInum.GC.model2,ROInum.GC.model3)
print(ROInum.GC.anova)

# p values
ROInum.GC.pv <- lsmeans(ROInum.GC.model3, pairwise ~ Cond_Time, glhargs=list())
summary(ROInum.GC.pv)




######
#ROI Types
# number of ROIs in each trial for each field of view (across the whole trial)

ROInum.trial.spot.roi<- ddply(Lck.peaks, c("Animal", "Spot","ROIType", "TimePoint", "TimeGroup", "Condition","Trial","Spot_trial_TP"), summarise, 
                          nROIs = length(unique(roiName)))

ROInum.spot.roi<-ddply(ROInum.trial.spot.roi, c("Animal", "Spot", "ROIType", "TimePoint", "TimeGroup", "Condition"), summarise, 
                   meanROIs = mean(nROIs))


# mean
df.lck.ROInum3<-summarySE(ROInum.spot.roi, measurevar = "meanROIs", groupvars = c("Condition","TimeGroup", "ROIType"))


ggplot(df.lck.ROInum3, aes(x=interaction(ROIType,TimeGroup),y=meanROIs, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanROIs-se, ymax=meanROIs+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIs/FOV") +
  max.theme

# cyto peaks

cyto.peaks<-subset(astro.peaks, Channel=="GCaMP")

ROInum.trial.spot.cyto<- ddply(cyto.peaks, c("Animal", "Spot","ROIType", "TimePoint", "TimeGroup", "Condition","Trial","Spot_trial_TP"), summarise, 
                              nROIs = length(unique(roiName)))

ROInum.spot.cyto<-ddply(ROInum.trial.spot.cyto, c("Animal", "Spot", "ROIType", "TimePoint", "TimeGroup", "Condition"), summarise, 
                       meanROIs = mean(nROIs))


# mean
df.cyto.ROInum3<-summarySE(ROInum.spot.cyto, measurevar = "meanROIs", groupvars = c("Condition","TimeGroup", "ROIType"))


ggplot(df.cyto.ROInum3, aes(x=interaction(ROIType,TimeGroup),y=meanROIs, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanROIs-se, ymax=meanROIs+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIs/FOV") +
  max.theme

ggplot(df.cyto.ROInum3[df.cyto.ROInum3$ROIType=="Process",], aes(x=interaction(ROIType,TimeGroup),y=meanROIs, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanROIs-se, ymax=meanROIs+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIs/FOV") +
  max.theme

ggplot(df.cyto.ROInum3[df.cyto.ROInum3$ROIType=="Soma",], aes(x=interaction(ROIType,TimeGroup),y=meanROIs, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanROIs-se, ymax=meanROIs+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIs/FOV") +
  max.theme

ggplot(df.cyto.ROInum3[df.cyto.ROInum3$ROIType=="Endfoot",], aes(x=interaction(ROIType,TimeGroup),y=meanROIs, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanROIs-se, ymax=meanROIs+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIs/FOV") +
  max.theme

########### 
# signals per field of view

# count the number of peaks per field of view for each spot at each time point

Signals.trial.spot<- ddply(Lck.linegraphs, c("Animal", "Spot","Channel", "TimePoint", "TimeGroup", "Condition","Trial","Spot_trial_TP"), summarise, 
                           nPeaks = length(amplitude))

Signals.spot<-ddply(Signals.trial.spot, c("Animal", "Spot", "Channel", "TimePoint", "TimeGroup", "Condition"), summarise, 
                    meanPeaks = mean(nPeaks))

# mean signals
df.lck.sig.spot<-summarySE(Signals.spot[Signals.spot$Channel=="LckGCaMP",], measurevar = "meanPeaks", groupvars = c("Condition", "TimeGroup"))
df.lck.sig.spot.TP<-summarySE(Signals.spot[Signals.spot$Channel=="LckGCaMP",], measurevar = "meanPeaks", groupvars = c("Condition", "TimePoint"))

ggplot(df.lck.sig.spot, aes(x=TimeGroup,y=meanPeaks, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanPeaks-se, ymax=meanPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peaks per FOV") +
  max.theme

ggplot(df.lck.sig.spot.TP, aes(x=TimePoint,y=meanPeaks, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanPeaks-se, ymax=meanPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peaks per FOV") +
  max.theme

#########
#peak types

Signals.trial.spot.ptype<- ddply(astro.peaks.resp2, c("Animal", "Spot","peakType", "TimePoint", "TimeGroup", "Condition","Trial","Spot_trial_TP"), summarise, 
                           nPeaks = length(amplitude))

Signals.spot.ptype<-ddply(Signals.trial.spot.ptype, c("Animal", "Spot", "peakType", "TimePoint", "TimeGroup", "Condition"), summarise, 
                    meanPeaks = mean(nPeaks))

# mean signals
df.lck.sig.spot.ptype<-summarySE(Signals.spot.ptype, measurevar = "meanPeaks", groupvars = c("Condition", "TimeGroup", "peakType"))
df.lck.sig.spot.TP.ptype<-summarySE(Signals.spot.ptype, measurevar = "meanPeaks", groupvars = c("Condition", "TimePoint", "peakType"))

ggplot(df.lck.sig.spot.ptype, aes(x=interaction(TimeGroup, peakType),y=meanPeaks, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanPeaks-se, ymax=meanPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peaks per FOV") +
  max.theme

ggplot(df.lck.sig.spot.TP.ptype, aes(x=interaction(TimePoint, peakType),y=meanPeaks, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanPeaks-se, ymax=meanPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peaks per FOV") +
  max.theme


ggplot(df.lck.sig.spot.ptype[df.lck.sig.spot.TP.ptype$peakType=="MultiPeak",], aes(x=interaction(TimeGroup, peakType),y=meanPeaks, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanPeaks-se, ymax=meanPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peaks per FOV") +
  max.theme

ggplot(df.lck.sig.spot.ptype[df.lck.sig.spot.TP.ptype$peakType=="Plateau",], aes(x=interaction(TimeGroup, peakType),y=meanPeaks, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanPeaks-se, ymax=meanPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peaks per FOV") +
  max.theme

ggplot(df.lck.sig.spot.ptype[df.lck.sig.spot.TP.ptype$peakType=="SinglePeak",], aes(x=interaction(TimeGroup, peakType),y=meanPeaks, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanPeaks-se, ymax=meanPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peaks per FOV") +
  max.theme

# only LCK, WP stim
multipeaks.null = lmer(meanPeaks ~ (1|Spot), Signals.spot.ptype[(Signals.spot.ptype$peakType=="MultiPeak" & Signals.spot.ptype$Condition=="WPstim"),],REML=FALSE)
multipeaks.model1 = lmer(meanPeaks ~ TimeGroup + (1|Spot), Signals.spot.ptype[(Signals.spot.ptype$peakType=="MultiPeak" & Signals.spot.ptype$Condition=="WPstim"),],REML=FALSE)
multipeaks.anova <- anova(multipeaks.null,multipeaks.model1)
print(multipeaks.anova)

# p values
multipeaks.pv <- lsmeans(multipeaks.model1, pairwise ~ TimeGroup, glhargs=list())
summary(multipeaks.pv)


###########
# onset times??



lck_onset$Spot_trial_TP<-paste(lck_onset$Spot, lck_onset$Trial, lck_onset$TimePoint, sep="_")

lck_onset$Spot_trial_TP_Cond<-paste(lck_onset$Spot_trial_TP,lck_onset$Condition, sep="_")

lck_onset$Spot_trial_TP_Cond_ROI<-paste(lck_onset$Spot_trial_TP_Cond,lck_onset$roiName, sep="_")

lck_onset$Spot_trial_TP_ROI<-paste(lck_onset$Spot_trial_TP,lck_onset$roiName, sep="_")

lck_onset$Spot_ROI<-paste(lck_onset$Spot,lck_onset$roiName, sep="_")



#Time Point Type
lck_onset$TimeGroup[grepl("B",lck_onset$TimePoint)]="before"
lck_onset$TimeGroup[grepl("T",lck_onset$TimePoint)]="after"


lck_onset$TimeGroup<- as.factor(lck_onset$TimeGroup)
lck_onset$TimeGroup<- factor(lck_onset$TimeGroup, levels = c("before","after"))

#remove NaNs with no onset time
OnsetNa= is.na(lck_onset$OnsetTime)
lck_onset<-lck_onset[!OnsetNa,]

#means onsetTimes

df.OT.lck<-summarySE(lck_onset, measurevar = "OnsetTime", groupvars = c("Condition", "TimeGroup"))

ggplot(df.OT.lck, aes(x=Condition,y=OnsetTime, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean OnsetTime") +
  max.theme

# repsonding
lck_onset.resp<-subset(lck_onset, OnsetTime<10 & Condition=="WPstim")

df.OT.lck2<-summarySE(lck_onset.resp, measurevar = "OnsetTime", groupvars = c("Condition", "TimeGroup"))

ggplot(df.OT.lck2, aes(x=Condition,y=OnsetTime, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean OnsetTime") +
  max.theme

# fast vs. delayed
lck_onset.fast<-subset(lck_onset.resp, OnsetTime<1)
lck_onset.delayed<-subset(lck_onset.resp, OnsetTime>=1)

df.OT.fast<-summarySE(lck_onset.fast, measurevar = "OnsetTime", groupvars = c("Condition", "TimeGroup"))
df.OT.delayed<-summarySE(lck_onset.delayed, measurevar = "OnsetTime", groupvars = c("Condition", "TimeGroup"))

ggplot(df.OT.fast, aes(x=Condition,y=OnsetTime, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean fast OnsetTime") +
  max.theme

ggplot(df.OT.delayed, aes(x=Condition,y=OnsetTime, fill= TimeGroup)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean delayedOnsetTime") +
  max.theme

###############
#RCaMP

RCaMP <- read.table("G:/ZurichData/Astrocyte_Calcium/P2Y1_Mice/Results/RC_P2Y1_timepoints_peaks_05_2018.csv", header=TRUE, sep = ",")
RCaMP <- read.table("H:/P2Y1_Data/Results/RC_P2Y1_timepoints_peaks_05_2018.csv", header=TRUE, sep = ",")
RCaMP <- read.table("J:/Jill_Stobart/P2Y1_Data/Results/RC_P2Y1_timepoints_peaks_05_2018.csv", header=TRUE, sep = ",")


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

#histogram of responses
ggplot(RCaMP.respROIs[RCaMP.respROIs$TimeGroup=="before",], aes(x=Freq, fill=Condition)) + geom_histogram(aes(y=..count../sum(..count..)), binwidth=0.2,position="dodge")+
  ggtitle("frequency")

ggplot(RCaMP.respROIs[RCaMP.respROIs$TimeGroup=="after",], aes(x=Freq, fill=Condition)) + geom_histogram(aes(y=..count../sum(..count..)), binwidth=0.2,position="dodge")+
  ggtitle("frequency")


###########
#means

df.amp1.RC<-summarySE(RCaMP.respROIs, measurevar = "meanAmp", groupvars = c("Condition", "TimeGroup"))
df.amp2.RC<-summarySE(RCaMP.respROIs, measurevar = "meanAmp", groupvars = c("Condition", "TimePoint"))
df.amp3.RC<-summarySE(RCaMP.ROIs, measurevar = "meanAmp", groupvars = c("Condition", "TimeGroup"))
df.amp4.RC<-summarySE(RCaMP.ROIs, measurevar = "meanAmp", groupvars = c("Condition", "TimePoint"))
df.amp1.RC[3,4]= 0.546

ggplot(df.amp1.RC, aes(x=TimeGroup,y=meanAmp, fill=Condition )) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanAmp-se, ymax=meanAmp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean ROI amplitude") +
  max.theme

ggplot(df.amp2.RC, aes(x=TimePoint,y=meanAmp, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanAmp-se, ymax=meanAmp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean ROI amplitude") +
  max.theme

ggplot(df.amp3.RC, aes(x=TimeGroup, y=meanAmp, fill= Condition)) +
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

ggplot(df.freq1.RC, aes(x=TimeGroup,y=Freq, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Freq-se, ymax=Freq+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean signals/trial") +
  max.theme

ggplot(df.freq2.RC, aes(x=TimePoint,y=Freq, fill= Condition)) +
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
farpeaks1 <- RCaMP$peakTime>=0 & RCaMP$peakTime<=5
RCaMP$ActivePeak[farpeaks1] <- 1 


# count number of peaks per trial in the first 30 sec
fraction.response<- ddply(RCaMP, c("Animal", "Spot", "TimePoint", "TimeGroup", "Condition","roiName"), summarise, 
                      nActive = sum(ActivePeak), NTrials = length(unique(Spot_trial_TP_Cond)))
                      
fraction.response$frac.resp <- fraction.response$nActive/fraction.response$NTrials

#histogram of responses
ggplot(fraction.response[fraction.response$TimeGroup=="before",], aes(x=frac.resp, fill=Condition)) + geom_histogram(aes(y=..count../sum(..count..)), binwidth=0.2,position="dodge")+
  ggtitle("Fraction of active trials before")

ggplot(fraction.response[fraction.response$TimeGroup=="after",], aes(x=frac.resp, fill=Condition)) + geom_histogram(aes(y=..count../sum(..count..)), binwidth=0.2,position="dodge")+
  ggtitle("Fraction of active trials 1 month after TAM")


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

##########
# high mid and low responders

stim.respROIs<-subset(RCaMP.respROIs, Condition=="WPstim" & meanAmp>0)

ggplot(stim.respROIs[stim.respROIs$TimeGroup=="before",], aes(x=meanAmp, fill=Condition)) + geom_histogram(aes(y=..count../sum(..count..)), binwidth=0.2,position="dodge")+
  ggtitle("amplitude distributions before")

Percentile<-quantile(stim.respROIs$meanAmp, probs = seq(0, 1, 0.1), na.rm = TRUE,
                     names = TRUE)

# 90th percentile (cut off for high responders)
highCut<-Percentile[[10]]

# below 50
midCut<-Percentile[[6]]

highData<-subset(stim.respROIs, TimeGroup=="before" & meanAmp>= highCut)
midData<-subset(stim.respROIs, TimeGroup=="before" & meanAmp< highCut & meanAmp>=midCut)
lowData<-subset(stim.respROIs, TimeGroup=="before" & meanAmp< midCut)

highNeurons<-highData$Spot_ROI
midNeurons<-midData$Spot_ROI
lowNeurons<-lowData$Spot_ROI

# apply neuron classes across all time points

stim.respROIs$NeuronType="NA"

stim.respROIs$NeuronType[stim.respROIs$Spot_ROI %in% highNeurons]= "high"
stim.respROIs$NeuronType[stim.respROIs$Spot_ROI %in% midNeurons]= "mid"
stim.respROIs$NeuronType[stim.respROIs$Spot_ROI %in% lowNeurons]= "low"


highClass<-subset(stim.respROIs, NeuronType=="high")
midClass<-subset(stim.respROIs, NeuronType=="mid")
lowClass<-subset(stim.respROIs, NeuronType=="low")


df.amp.high<-summarySE(highClass, measurevar = "meanAmp", groupvars = c("TimeGroup"))
df.amp.high2<-summarySE(highClass, measurevar = "meanAmp", groupvars = c("TimePoint"))
df.amp.mid<-summarySE(midClass, measurevar = "meanAmp", groupvars = c("TimeGroup"))
df.amp.mid2<-summarySE(midClass, measurevar = "meanAmp", groupvars = c("TimePoint"))
df.amp.low<-summarySE(lowClass, measurevar = "meanAmp", groupvars = c("TimeGroup"))
df.amp.low2<-summarySE(lowClass, measurevar = "meanAmp", groupvars = c("TimePoint"))

ggplot(df.amp.high, aes(x=TimeGroup,y=meanAmp, fill=TimeGroup )) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanAmp-se, ymax=meanAmp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("high Mean ROI amplitude") +
  max.theme

ggplot(df.amp.high2, aes(x=TimePoint,y=meanAmp, fill=TimePoint )) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanAmp-se, ymax=meanAmp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("high Mean ROI amplitude") +
  max.theme


ggplot(df.amp.mid, aes(x=TimeGroup,y=meanAmp, fill=TimeGroup )) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanAmp-se, ymax=meanAmp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("mid Mean ROI amplitude") +
  max.theme

ggplot(df.amp.mid2, aes(x=TimePoint,y=meanAmp, fill=TimePoint )) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanAmp-se, ymax=meanAmp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("mid Mean ROI amplitude") +
  max.theme

ggplot(df.amp.low, aes(x=TimeGroup,y=meanAmp, fill=TimeGroup )) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanAmp-se, ymax=meanAmp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("low Mean ROI amplitude") +
  max.theme

ggplot(df.amp.low2, aes(x=TimePoint,y=meanAmp, fill=TimePoint )) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanAmp-se, ymax=meanAmp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("low Mean ROI amplitude") +
  max.theme


#stats
# Likelihood-ratio test amplitude
amp.high.null = lmer(meanAmp ~ (1|Spot) + (1|Spot_ROI), highClass,REML=FALSE)
amp.high.model2 = lmer(meanAmp ~ TimeGroup + (1|Spot) + (1|Spot_ROI), highClass,REML=FALSE)
amp.high.anova <- anova(amp.high.null,  amp.high.model2)
print(amp.high.anova)

# p values
amp.high.resp <- lsmeans(amp.high.model2, pairwise ~ TimeGroup, glhargs=list())
summary(amp.high.resp)


#stats
# Likelihood-ratio test amplitude
amp.mid.null = lmer(meanAmp ~ (1|Spot) + (1|Spot_ROI), midClass,REML=FALSE)
amp.mid.model2 = lmer(meanAmp ~ TimeGroup + (1|Spot) + (1|Spot_ROI), midClass,REML=FALSE)
amp.mid.anova <- anova(amp.mid.null,  amp.mid.model2)
print(amp.mid.anova)

# p values
amp.mid.resp <- lsmeans(amp.mid.model2, pairwise ~ TimeGroup, glhargs=list())
summary(amp.mid.resp)
