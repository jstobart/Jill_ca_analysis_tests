
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

########################

# theme for plots
max.theme <- theme_classic() + 
  theme(
    text=element_text(size=12),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14, face="bold"),
    axis.title.y=element_text(vjust=1),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14, face="bold"))


###########
# NOTES


########################
# relative "active ROI" between groups based on peak auc and frequency

# whole frame and automatic RCaMP ROI selection:
stim.all <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/S&LStim_LckGC&RC_14_02_2017.csv", header=TRUE, sep = ",")

lsm.options(pbkrtest.limit = 100000)

# treatment
stim.all$treatment<-"Control"

# stim onset at 5 sec
stim.all$peakTime<-stim.all$peakTime-5
stim.all$peakStart<-stim.all$peakStart-5
stim.all$peakStartHalf<-stim.all$peakStartHalf-5

#duration
stim.all$Duration<-stim.all$halfWidth*2

# ROITypes
stim.all$ROIType= 0
stim.all$ROIType[grepl("r",stim.all$ROIname)]="Process"
stim.all$ROIType[grepl("E",stim.all$ROIname)]="Endfoot"
stim.all$ROIType[grepl("D",stim.all$ROIname)]="Dendrite"
stim.all$ROIType[grepl("N",stim.all$ROIname)]="Neuron"
stim.all$ROIType[grepl("np",stim.all$ROIname)]="Neuropil"

stim.all$ROIType<- as.factor(stim.all$ROIType)


#unique ROI names
stim.all$ROIs_trial<-paste(stim.all$Animal, stim.all$Spot, stim.all$Trial,stim.all$ROIname, sep= "_")

stim.all$ROIs<-paste(stim.all$Animal, stim.all$Spot, stim.all$ROIname, sep= "_")

stim.all$trials<-paste(stim.all$Animal, stim.all$Spot, stim.all$Trial, sep= "_")


# exclude RG10 cause it moves too much
#stim.all<- subset(stim.all, Animal!="RG10")



# remove matching astrocyte process and soma ROIs
Overlap= stim.all$overlap!=0

OverlapROIs<-unique(stim.all$ROIs_trial[Overlap])

stim.all2<-stim.all[!Overlap,]


#######
# histogram of peak times for stim
longstim<-subset(stim.all2, Condition=="Stim")
ggplot(longstim, aes(x=peakTime, fill=Channel)) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("long stim")

shortstim<-subset(stim.all2, Condition=="shortstim")
ggplot(shortstim, aes(x=peakTime, fill=Channel)) + geom_histogram(binwidth=2, position="dodge") +
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

#considering treatment
RCaMP.longstim<- subset(RCaMP, Condition=="Stim")
ggplot(RCaMP.longstim, aes(x=peakTime, fill=treatment)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("RCaMP long stim")

RCaMP.shortstim<- subset(RCaMP, Condition=="shortstim")
ggplot(RCaMP.shortstim, aes(x=peakTime, fill=treatment)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("RCaMP short stim")

#considering treatment
GCaMP.longstim<- subset(GCaMP, Condition=="Stim")
ggplot(GCaMP.longstim, aes(x=peakTime, fill=treatment)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("GCaMP long stim")

GCaMP.shortstim<- subset(GCaMP, Condition=="shortstim")
ggplot(GCaMP.shortstim, aes(x=peakTime, fill=treatment)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("GCaMP short stim")
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

responding.trials_long<-unique(responding.neurons_long$trials)  # 215 trials of 253- 84.9% of trials
responding.trials_short<-unique(responding.neurons_short$trials) # 46 trials of 121- 38.0% of trials

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
farpeaks1 <- longstim.responding$peakTime>0 & longstim.responding$peakTime<9 & longstim.responding$Duration<11 & longstim.responding$ROIType=="Neuron"
longstim.responding$ActivePeak[farpeaks1] <- 1 

#responding astrocytes
farpeaks2 <- longstim.responding$peakTime>0 & longstim.responding$peakTime<20 & longstim.responding$ROIType!="Neuron"
longstim.responding$ActivePeak[farpeaks2] <- 1 

# pull out only the peaks that occur around the stimulation
longstim.stimwindow<- subset(longstim.responding, peakTime>=0 & peakTime<=20)



ggplot(longstim.stimwindow, aes(x=peakTime, fill=interaction(Channel,treatment))) + geom_histogram(binwidth=0.5, position="dodge") +
  ggtitle("long stim responding")

longstim.after<- subset(longstim.responding, peakTime>20 & peakTime<80)

ggplot(longstim.after, aes(x=peakTime, fill=interaction(Channel,treatment))) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("long stim after")

#####
library(xlsx)
respondingNeurons_Astrocytes=subset(longstim.responding, ActivePeak==1)
write.xlsx(respondingNeurons_Astrocytes, "E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/respondingROIs_longstim.xlsx")


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

df2A1<-summarySE(longstim.propActive, measurevar="propActive", groupvars=c("ROIType","treatment"))
df2A2<-summarySE(longstim.propActive, measurevar="propActive", groupvars=c("Channel","treatment"))

ggplot(longstim.propActive, aes(x=propActive, fill=Channel)) + geom_histogram(binwidth=0.05, position="dodge") +
  ggtitle("long stim prop active ")

df2A1$ROIType <- factor(df2A1$ROIType , levels = c("Neuron","Dendrite","Neuropil","Endfoot","Process"))
ggplot(data=df2A1, aes(x=ROIType, y=propActive, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=propActive-se, ymax=propActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("ROIType") +
  ylab("Proportion of Responding ROIs") +
  ggtitle("Responding ROIs for long stim") 


ggplot(data=df2A2, aes(x=Channel, y=propActive, fill=Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=propActive-se, ymax=propActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Channel") +
  ylab("Proportion of Responding ROIs") +
  ggtitle("Responding ROIs for long stim") 





#########
# SHORT STIM (90Hz, 1sec)
shortstim.responding$ActivePeak <- 0
farpeaks1 <- shortstim.responding$peakTime>0 & shortstim.responding$peakTime<2 & shortstim.responding$Duration<3 & shortstim.responding$ROIType=="Neuron"
shortstim.responding$ActivePeak[farpeaks1] <- 1 

farpeaks2 <- shortstim.responding$peakTime>0 & shortstim.responding$peakTime<20 & shortstim.responding$ROIType!="Neuron"
shortstim.responding$ActivePeak[farpeaks2] <- 1 

# pull out only the peaks that occur around the stimulation
shortstim.stimwindow<- subset(shortstim.responding, peakTime>=0 & peakTime<=10)

ggplot(shortstim.stimwindow, aes(x=peakTime, fill=interaction(Channel,treatment))) + geom_histogram(binwidth=0.5, position="dodge") +
  ggtitle("short stim responding")

ggplot(shortstim.stimwindow, aes(x=peakTime, y=prominence, color=interaction(Channel,treatment))) + geom_point() +
  ggtitle("short stim responding")

shortstim.after<- subset(shortstim.responding, peakTime>10 & peakTime<80)

ggplot(shortstim.after, aes(x=peakTime, fill=interaction(Channel,treatment))) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("short stim after")

#respondingNeurons_AstrocytesShort=subset(shortstim.responding, ActivePeak==1)
#write.xlsx(respondingNeurons_Astrocytes, "E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/respondingROIs_shortstim.xlsx")



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


df2A1$ROIType <- factor(df2A1$ROIType , levels = c("Neuron","Dendrite","Neuropil","Endfoot","Process"))
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

df9A <- summarySE(longstim.stimwindow, measurevar="peakTime", groupvars=c("ROIType","treatment"))
df9B <- summarySE(shortstim.stimwindow, measurevar="peakTime", groupvars=c("ROIType","treatment"))

df10A <- summarySE(longstim.stimwindow, measurevar="peakStart", groupvars=c("ROIType","treatment"))
df10B <- summarySE(shortstim.stimwindow, measurevar="peakStart", groupvars=c("ROIType","treatment"))

df11A <- summarySE(longstim.stimwindow, measurevar="peakStartHalf", groupvars=c("ROIType","treatment"))
df11B <- summarySE(shortstim.stimwindow, measurevar="peakStartHalf", groupvars=c("ROIType","treatment"))

df9A$ROIType <- factor(df9A$ROIType , levels = c("Neuron","Dendrite","Neuropil","Endfoot","Process"))
ggplot(data=df9A, aes(x=ROIType, y=peakTime, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("Mean Peak Time (s)") +
  ggtitle("max peak time for long stim") 

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


ggplot(data=df11A, aes(x=ROIType, y=peakStartHalf, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStartHalf-se, ymax=peakStartHalf+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("peakStartHalf Time (t1/2)") +
  ggtitle("peakStartHalf times for long stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
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





#######
# peak features

# amplitude
df2A <- summarySE(longstim.stimwindow, measurevar="amplitude", groupvars=c("ROIType","treatment"),na.rm=TRUE)


ggplot(data=df2A, aes(x=ROIType, y=amplitude, fill=treatment)) +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

# prominence

df3A <- summarySE(longstim.stimwindow, measurevar="prominence", groupvars=c("ROIType","treatment"))

ggplot(data=df3A, aes(x=treatment, y=prominence, fill=ROIType)) +
  geom_errorbar(aes(ymin=prominence-se, ymax=prominence+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("prominence") +
  max.theme

# duration
df4A <- summarySE(longstim.stimwindow, measurevar="Duration", groupvars=c("ROIType","treatment"))

ggplot(data=df4A, aes(x=treatment, y=Duration, fill=ROIType)) +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Duration") +
  max.theme



#frequency (within the same trial)

# nostim frequency
#aggregate data by trial
nostim.trials<- ddply(nostim, c("Animal", "Spot", "trials","Channel","ROIType","treatment","ROIs_trial"), summarise, 
                               PA_mean = mean(peakAUC), nEvents = length(peakAUC),
                               Dur_mean = mean(Duration), Prom_mean = mean(prominence),
                               amp_mean = mean(amplitude), HalfDur = mean(halfWidth),
                              freq = sum(numPeaks), area_mean= mean(area))
nostim.trials$peaks_min=nostim.trials$freq/1.5
nostim.trials$signals_min=nostim.trials$nEvents/1.5

df5A <- summarySE(nostim.trials, measurevar="peaks_min", groupvars=c("Channel"), na.rm = T)
df5B <- summarySE(nostim.trials, measurevar="signals_min", groupvars=c("Channel"),na.rm = T)

ggplot(data=df4A, aes(x=treatment, y=Duration, fill=Condition)) +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Duration") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme
######################
