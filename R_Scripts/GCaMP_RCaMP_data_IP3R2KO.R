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
# TO DO:
# onset time distribution
#fract active pixels
#onset times, peak times of fast vs. delayed

#amp, duration RC and GC


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

all.lck.peaks <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/Peaks_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
all.lck.OT<-read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/OnsetTimes_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")

##### 
#home files

all.lck.peaks <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/Peaks_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
all.lck.OT<-read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/OnsetTimes_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")


##########

lsm.options(pbkrtest.limit = 100000)

# only consider IP3R2KO
all.lck.peaks<-all.lck.peaks[grepl("IP",all.lck.peaks$Animal),]
all.lck.OT<-all.lck.OT[grepl("IP",all.lck.OT$Animal),]

all.lck.OT$Spot_trial_Cond<-paste(all.lck.OT$Spot_trial, all.lck.OT$Condition, sep="_")
all.lck.OT$ROIs_Cond<-paste(all.lck.OT$ROIs_trial, all.lck.OT$Condition, sep="_")


# REMOVE duplicate entries from onset time and peak time data frames 
# only the first entry will be used
all.lck.OT2=all.lck.OT[order(all.lck.OT$OnsetTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.lck.OT<-distinct(all.lck.OT2, ROIs_Cond, .keep_all = TRUE)

#adjust onset time for the data with 2 s before stimulation included:

rm(all.lck.OT2)

# Genotype
all.lck.peaks$Genotype="IP3R2_WT"
all.lck.peaks$Genotype[grepl("IPRG1",all.lck.peaks$Animal)]="IP3R2_KO"
all.lck.peaks$Genotype[grepl("IPRG4",all.lck.peaks$Animal)]="IP3R2_KO"

all.lck.OT$Genotype="IP3R2_WT"
all.lck.OT$Genotype[grepl("IPRG1",all.lck.OT$Animal)]="IP3R2_KO"
all.lck.OT$Genotype[grepl("IPRG4",all.lck.OT$Animal)]="IP3R2_KO"


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





# neuronal responses to stimulation
NeuronalStim<-subset(all.lck.OT, Channel=="RCaMP" & Condition=="Stim")

# should have an onset time in 8 s stimulus
ggplot(NeuronalStim[NeuronalStim$OnsetTime<10,],aes(x=OnsetTime,y=..density..,fill=Genotype)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP onset times between 0 and 10 s from stim trials")+
  max.theme


Neuron95Onset<-quantile(NeuronalStim$OnsetTime[NeuronalStim$OnsetTime<8], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
print(Neuron95Onset)


AstroStim<-subset(all.lck.OT, Channel=="GCaMP" & Condition=="Stim")

# should have an onset time in 8 s stimulus
ggplot(AstroStim[AstroStim$OnsetTime<15,],aes(x=OnsetTime,fill=Genotype)) +
  geom_histogram(binwidth=(0.084*5), position="dodge") +
  ggtitle("Lck-GCaMP onset times between 0 and 15 s from stim trials")+
  max.theme


######
# for median and mean calculations

# no stim vs 8 s stim- 
#neuronal window=9s, AC window= 15 s for peak time, 
#neuronal window=2s, AC window=12 s for onset


LongN_PTwind2=9
LongAC_PTwind2=15

LongN_OTwind2=2
LongAC_OTwind2=12

# remove data that is outside the above windows
# lck
stim.lck.OT.R<-subset(all.lck.OT, Channel=="RCaMP" & OnsetTime<=LongN_OTwind2)
stim.lck.OT.G<-subset(all.lck.OT, Channel=="GCaMP" & OnsetTime<=LongAC_OTwind2)

stim.lck.OT.window<-rbind(stim.lck.OT.R, stim.lck.OT.G)


# peak times
# lck
stim.lck.PT.R<-subset(all.lck.peaks, Channel=="RCaMP" & peakTime<=LongN_PTwind2 & peakTime>=0 & Duration<45)
stim.lck.PT.G<-subset(all.lck.peaks, Channel=="GCaMP" & peakTime<=LongAC_PTwind2 & peakTime>=0 & Duration<45)

stim.lck.peaks.window<-rbind(stim.lck.PT.R, stim.lck.PT.G)


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
#stim.both.alldata$Channel <- factor(stim.both.alldata$Channel, levels = c("cyto_RCaMP","cyto_GCaMP","lck_RCaMP","lck_GCaMP"))


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
# active pixels from any ROI with a signal (not just those during stimulation)

# number of ROIs in each trial for each field of view (across the time window (2 s for neurons, 15 s for AC))

stim.lck.alldata$Channel <- factor(stim.lck.alldata$Channel, levels = c("RCaMP","GCaMP"))

spot.lck.stim<-ddply(stim.lck.alldata, c("Animal","Spot","Genotype","Condition","Channel"), summarise, nROIs=length(unique(ROIs_Cond)), nFluoPix=nFluoPix,
                      nActivePix=nActivePix, PixelSize=pixelsize)

# add in number of trials
spot.lck.stim$Ani_Spot_Cond<-paste(spot.lck.stim$Animal, spot.lck.stim$Spot, spot.lck.stim$Condition, sep="_")
Spot.lck.ntrials$Ani_Spot_Cond<-paste(Spot.lck.ntrials$Animal, Spot.lck.ntrials$Spot, Spot.lck.ntrials$Condition, sep="_")

spot.lck.stim<-merge(spot.lck.stim, Spot.lck.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)

spot.lck.stim$ROIsPerTrial<-spot.lck.stim$nROIs/spot.lck.stim$nTrials

spot.lck.stim$FracActive=spot.lck.stim$nActivePix/spot.lck.stim$nFluoPix
spot.lck.stim$FracActivePerROI=spot.lck.stim$FracActive/spot.lck.stim$nROIs

# mean fraction of active pixels
df.FracActive<- summarySE(spot.lck.stim, measurevar = "FracActive", groupvars = c("Condition"))

ggplot(data=df.FracActive, aes(x=Condition, y= FracActive, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=FracActive-se, ymax=FracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# mean number of total ROIs per trial
df.lck.ROInum.mean<-summarySE(spot.lck.stim, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Genotype","Condition"))


ggplot(df.lck.ROInum.mean, aes(x=interaction(Genotype,Channel),y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view") +
  max.theme

Condition_Channel2= interaction(spot.lck.stim$Condition,spot.lck.stim$Channel)
Condition_Channel_Genotype= interaction(spot.lck.stim$Condition,spot.lck.stim$Channel,spot.lck.stim$Genotype)
nROI.lck.stim.null = lmer(ROIsPerTrial ~ (1|Animal), spot.lck.stim,REML=FALSE)
nROI.lck.stim.model1 = lmer(ROIsPerTrial~ Channel + (1|Animal), spot.lck.stim,REML=FALSE)
nROI.lck.stim.model2 = lmer(ROIsPerTrial ~ Condition + (1|Animal), spot.lck.stim,REML=FALSE)
nROI.lck.stim.model3 = lmer(ROIsPerTrial ~ Genotype + (1|Animal), spot.lck.stim,REML=FALSE)
nROI.lck.stim.model4 = lmer(ROIsPerTrial ~ Condition_Channel2 + (1|Animal), spot.lck.stim,REML=FALSE)
nROI.lck.stim.model5 = lmer(ROIsPerTrial ~ Condition_Channel2 + Genotype + (1|Animal), spot.lck.stim,REML=FALSE)
nROI.lck.stim.model6 = lmer(ROIsPerTrial ~ Condition_Channel_Genotype + (1|Animal), spot.lck.stim,REML=FALSE)
nROI.lck.stim.anova <- anova(nROI.lck.stim.null, nROI.lck.stim.model1,nROI.lck.stim.model2,nROI.lck.stim.model3,
                             nROI.lck.stim.model4,nROI.lck.stim.model5,nROI.lck.stim.model6)
print(nROI.lck.stim.anova)

nROI.lck.stim.Cond_Channel_Genotype<- glht(nROI.lck.stim.model6, mcp(Condition_Channel_Genotype= "Tukey"))
summary(nROI.lck.stim.Cond_Channel_Genotype)



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
df.OT2<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata$Channel=="GCaMP",], measurevar = "OnsetTime", groupvars = c("Group", "Genotype"))
df.OT3<- summarySE(stim.lck.compdata.STIM, measurevar = "OnsetTime", groupvars = c("ROIType","Channel_Group", "Genotype"))

ggplot(df.OT1, aes(x=interaction(Genotype,Channel_Group),y=OnsetTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.OT2, aes(x=Group,y=OnsetTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.OT3, aes(x=interaction(Channel_Group,ROIType),y=OnsetTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  scale_fill_manual(values=cbbPalette)+
  max.theme



# stats
Group_Channel_Type_Gen=interaction(stim.lck.alldata$Group,stim.lck.alldata$Channel,stim.lck.alldata$ROIType, stim.lck.alldata$Genotype)
Group_Channel_Gen=interaction(stim.lck.alldata$Group,stim.lck.alldata$Channel,stim.lck.alldata$Genotype)

# stats for onset times- neurons vs astrocytes
OT.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model2 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model3 = lmer(OnsetTime ~ Group_Channel_Gen + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model4 = lmer(OnsetTime ~ Group_Channel_Type_Gen + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.anova <- anova(OT.null, OT.model1,OT.model2,OT.model3,OT.model4)
print(OT.anova)

OT.Group_channel_Gen<- glht(OT.model3, mcp(Group_Channel_Gen= "Tukey"))
summary(OT.Group_channel_Gen)

OT.Group_Channel_Type_Gen<- glht(OT.model4, mcp(Group_Channel_Type_Gen= "Tukey"))
summary(OT.Group_Channel_Type_Gen)


Group_Channel_Type=interaction(stim.lck.compdata.STIM$Group,stim.lck.compdata.STIM$Channel,stim.lck.compdata.STIM$ROIType)
Group_Channel=interaction(stim.lck.compdata.STIM$Channel, stim.lck.compdata.STIM$Group)

# stats for onset times- neurons vs astrocytes
OT.stim.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.model3 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.model4 = lmer(OnsetTime ~ Group_Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.model6 = lmer(OnsetTime ~ Group_Channel_Type + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.anova <- anova(OT.stim.null, OT.stim.model1,OT.stim.model3,OT.stim.model4,OT.stim.model6)
print(OT.stim.anova)

OT.stim.Group_channel<- glht(OT.stim.model4, mcp(Group_Channel= "Tukey"))
summary(OT.stim.Group_channel)

OT.stim.Group_channel_type<- glht(OT.stim.model6, mcp(Group_Channel_Type= "Tukey"))
summary(OT.stim.Group_channel_type)

summary(OT.stim.model4)

# check residuals for linearity
plot(fitted(OT.stim.model4), residuals(OT.stim.model4),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(OT.stim.model4), residuals(OT.stim.model4)), col=46, lwd=2.5)
















########
#amplitude
df.amp1A<-summarySE(stim.both.alldata[stim.both.alldata$Channel=="cyto_GCaMP",], measurevar = "amplitude", groupvars = c("Channel", "Condition"))
df.amp1B<-summarySE(stim.both.alldata[stim.both.alldata$Channel=="lck_GCaMP",], measurevar = "amplitude", groupvars = c("Channel", "Condition"))

df.amp2<-summarySE(stim.both.alldata, measurevar = "amplitude", groupvars = c("Channel", "Group","Condition"))
df.amp3<- summarySE(stim.lck.compdata, measurevar = "amplitude", groupvars = c("Channel_Group","Condition"))
df.amp4<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim",], measurevar = "amplitude", groupvars = c("Channel_Group"))
df.amp5<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim"& stim.lck.compdata$Channel=="lck_GCaMP",], measurevar = "amplitude", groupvars = c("ROIType","Channel_Group"))
df.amp6<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="lck_RCaMP",], measurevar = "amplitude", groupvars = c("Group"))

ggplot(df.amp1A, aes(x=Condition,y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp1B, aes(x=Condition,y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp2, aes(x=interaction(Channel,Group),y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp3, aes(x=Channel_Group,y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp4, aes(x=Channel_Group,y=amplitude, fill= Channel_Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp5, aes(x=Channel_Group,y=amplitude, fill= ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.amp6, aes(x=Group,y=amplitude, fill= Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(stim.lck.compdata, aes(x=Channel_Group,y=amplitude, fill=Channel_Group)) +
  geom_boxplot(notch=TRUE)+
  ylab("amplitude") +
  ggtitle("lck data- fast vs delayed")+
  max.theme

ggplot(stim.lck.compdata, aes(x=amplitude, y=..density..,fill=Channel_Group)) +
  geom_histogram(binwidth=0.5, position="dodge")+
  ylab("density") +
  xlim(-2,15)+
  ggtitle("lck data- fast vs delayed")+
  max.theme

ggplot(stim.lck.alldata[stim.lck.alldata$Channel_Group=="lck_GCaMP.fast",], aes(x=amplitude, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.5, position="dodge")+
  ylab("density") +
  xlim(-2,15)+
  ggtitle("lck data- fast no stim vs stim")+
  max.theme

# for supplementary figure
ggplot(stim.both.alldata[stim.both.alldata$Channel=="cyto_GCaMP",], aes(x=amplitude, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.5, position="dodge")+
  ylab("density") +
  xlim(0,16)+
  ggtitle("cyto data- no stim vs stim")+
  max.theme

ggplot(stim.both.alldata[stim.both.alldata$Channel=="lck_GCaMP",], aes(x=amplitude, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.5, position="dodge")+
  ylab("density") +
  xlim(0,12)+
  ggtitle("lck data- no stim vs stim")+
  max.theme

Cond_Channel=interaction(stim.both.alldata$Condition,stim.both.alldata$Channel)
amp.both.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), stim.both.alldata,REML=FALSE)
amp.both.model1 = lmer(amplitude ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.both.alldata,REML=FALSE)
amp.both.model2 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.both.alldata,REML=FALSE)
amp.both.model3 = lmer(amplitude ~ Cond_Channel + (1|Animal) + (1|Spot) + (1|trials), stim.both.alldata,REML=FALSE)

amp.both.anova <- anova(amp.both.null, amp.both.model1,amp.both.model2,amp.both.model3)
print(amp.both.anova)

amp.both.Cond_channel<- glht(amp.both.model3, mcp(Cond_Channel= "Tukey"))
summary(amp.both.Cond_channel)

## lck data and group (fast, delayed)
comp_Channel_Group=interaction(stim.lck.compdata$Group,stim.lck.compdata$Channel)
comp_Channel_Group_Cond=interaction(stim.lck.compdata$Group,stim.lck.compdata$Channel,stim.lck.compdata$Condition)
comp_Channel_Group_Cond_Type=interaction(stim.lck.compdata$Group,stim.lck.compdata$Channel,
                                         stim.lck.compdata$Condition,stim.lck.compdata$ROIType)


# stats for onset times- neurons vs astrocytes
amp.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata,REML=FALSE)
amp.model1 = lmer(amplitude ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata,REML=FALSE)
amp.model2 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata,REML=FALSE)
amp.model3 = lmer(amplitude ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata,REML=FALSE)
amp.model4 = lmer(amplitude ~ comp_Channel_Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata,REML=FALSE)
amp.model5 = lmer(amplitude ~ comp_Channel_Group_Cond + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata,REML=FALSE)
amp.model6 = lmer(amplitude ~ comp_Channel_Group_Cond_Type + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata,REML=FALSE)
amp.anova <- anova(amp.null, amp.model1,amp.model2,amp.model3,amp.model4,amp.model5,amp.model6)
print(amp.anova)

amp.Group_channel<- glht(amp.model5, mcp(comp_Channel_Group_Cond= "Tukey"))
summary(amp.Group_channel)

amp.Group_channel_ty<- glht(amp.model6, mcp(comp_Channel_Group_Cond_Type= "Tukey"))
summary(amp.Group_channel_ty)

# because we used different sensors, separate RCaMP and GCaMP

#lck-GCaMP
Group_Cond_GC=interaction(stim.lck.compdata$Group[stim.lck.compdata$Channel=="lck_GCaMP"],stim.lck.compdata$Condition[stim.lck.compdata$Channel=="lck_GCaMP"])
Group_Cond_Type_GC=interaction(stim.lck.compdata$Group[stim.lck.compdata$Channel=="lck_GCaMP"],
                               stim.lck.compdata$Condition[stim.lck.compdata$Channel=="lck_GCaMP"],
                               stim.lck.compdata$ROIType[stim.lck.compdata$Channel=="lck_GCaMP"])

amp.null.GC = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="lck_GCaMP",],REML=FALSE)
amp.model2.GC  = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata[stim.lck.compdata$Channel=="lck_GCaMP",],REML=FALSE)
amp.model3.GC  = lmer(amplitude ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="lck_GCaMP",],REML=FALSE)
amp.model5.GC  = lmer(amplitude ~ Group_Cond_GC + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="lck_GCaMP",],REML=FALSE)
amp.model6.GC  = lmer(amplitude ~ Group_Cond_Type_GC + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="lck_GCaMP",],REML=FALSE)
amp.anova.GC  <- anova(amp.null.GC, amp.model2.GC,amp.model3.GC,amp.model5.GC,amp.model6.GC)
print(amp.anova.GC)

amp.Group_Cond.GC<- glht(amp.model5.GC, mcp(Group_Cond_GC= "Tukey"))
summary(amp.Group_Cond.GC)

amp.Group_channel_ty.GC<- glht(amp.model6.GC, mcp(Group_Cond_Type_GC= "Tukey"))
summary(amp.Group_channel_ty.GC)

#RCaMP
amp.null.RC = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="lck_RCaMP",],REML=FALSE)
amp.model2.RC  = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata[stim.lck.compdata$Channel=="lck_RCaMP",],REML=FALSE)
amp.anova.RC  <- anova(amp.null.RC, amp.model2.RC)
print(amp.anova.RC)

amp.Cond.RC<- glht(amp.model2.RC, mcp(Condition= "Tukey"))
summary(amp.Cond.RC)


#only consider STIM case
amp.null.stim = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="lck_RCaMP",],REML=FALSE)
amp.model1.stim  = lmer(amplitude ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="lck_RCaMP",],REML=FALSE)

amp.anova.stim  <- anova(amp.null.stim, amp.model1.stim)
print(amp.anova.stim)

amp.stim.group<- glht(amp.model1.stim, mcp(Group= "Tukey"))
summary(amp.stim.group)

########
#duration
df.Dur1<-summarySE(stim.both.alldata, measurevar = "Duration", groupvars = c("Channel", "Condition"))
df.Dur2<-summarySE(stim.both.alldata, measurevar = "Duration", groupvars = c("Channel", "Group","Condition"))
df.Dur3<- summarySE(stim.lck.compdata, measurevar = "Duration", groupvars = c("Channel_Group","Condition"))
df.Dur4<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim",], measurevar = "Duration", groupvars = c("Channel_Group"))
df.Dur5<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim"& stim.lck.compdata$Channel=="lck_GCaMP",], measurevar = "Duration", groupvars = c("ROIType","Channel_Group"))
df.Dur6<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="lck_RCaMP",], measurevar = "Duration", groupvars = c("Group"))

ggplot(df.Dur2, aes(x=interaction(Channel,Group),y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.Dur3, aes(x=Channel_Group,y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.Dur4, aes(x=Channel_Group,y=Duration, fill= Channel_Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.Dur5, aes(x=Channel_Group,y=Duration, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(df.Dur6, aes(x=Group,y=Duration, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(stim.lck.compdata, aes(x=Channel_Group,y=Duration, fill=Channel_Group)) +
  geom_boxplot(notch=TRUE)+
  ylab("Duration") +
  ggtitle("lck data- fast vs delayed")+
  max.theme

ggplot(stim.lck.compdata, aes(x=Duration, y=..density..,fill=Channel_Group)) +
  geom_histogram(binwidth=0.5, position="dodge")+
  ylab("density") +
  xlim(-2,15)+
  ggtitle("lck data- fast vs delayed")+
  max.theme

ggplot(stim.lck.alldata[stim.lck.alldata$Channel_Group=="lck_GCaMP.fast",], aes(x=Duration, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=0.5, position="dodge")+
  ylab("density") +
  xlim(-2,15)+
  ggtitle("lck data- fast no stim vs stim")+
  max.theme

# for supplementary figure
ggplot(stim.both.alldata[stim.both.alldata$Channel=="cyto_GCaMP",], aes(x=Duration, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=1, position="dodge")+
  ylab("density") +
  xlim(0,50)+
  ggtitle("cyto data- no stim vs stim")+
  max.theme

ggplot(stim.both.alldata[stim.both.alldata$Channel=="lck_GCaMP",], aes(x=Duration, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=1, position="dodge")+
  ylab("density") +
  xlim(0,50)+
  ggtitle("lck data- no stim vs stim")+
  max.theme

# stats for duration- neurons vs astrocytes
dur.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
dur.model1 = lmer(Duration ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
dur.model2 = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.alldata,REML=FALSE)
dur.model3 = lmer(Duration ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
dur.model4 = lmer(Duration ~ Channel_Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
dur.model5 = lmer(Duration ~ Group_Channel_Cond + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
dur.model6 = lmer(Duration ~ Group_Channel_Cond_Type + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
dur.anova <- anova(dur.null, dur.model1,dur.model2,dur.model3,dur.model4,dur.model5,dur.model6)
print(dur.anova)

dur.Group_channel<- glht(dur.model5, mcp(Group_Channel_Cond= "Tukey"))
summary(dur.Group_channel)

dur.Group_channel_ty<- glht(dur.model6, mcp(Group_Channel_Cond_Type= "Tukey"))
summary(dur.Group_channel_ty)



# because we used different sensors, separate RCaMP and GCaMP

#lck-GCaMP
dur.null.GC = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="lck_GCaMP",],REML=FALSE)
dur.model2.GC  = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata[stim.lck.compdata$Channel=="lck_GCaMP",],REML=FALSE)
dur.model3.GC  = lmer(Duration ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="lck_GCaMP",],REML=FALSE)
dur.model5.GC  = lmer(Duration ~ Group_Cond_GC + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="lck_GCaMP",],REML=FALSE)
dur.model6.GC  = lmer(Duration ~ Group_Cond_Type_GC + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="lck_GCaMP",],REML=FALSE)
dur.anova.GC  <- anova(dur.null.GC, dur.model2.GC,dur.model3.GC,dur.model5.GC,dur.model6.GC)
print(dur.anova.GC)

dur.Group_Cond.GC<- glht(dur.model5.GC, mcp(Group_Cond_GC= "Tukey"))
summary(dur.Group_Cond.GC)

dur.Group_channel_ty.GC<- glht(dur.model6.GC, mcp(Group_Cond_Type_GC= "Tukey"))
summary(dur.Group_channel_ty.GC)

#RCaMP
dur.null.RC = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata[stim.lck.compdata$Channel=="lck_RCaMP",],REML=FALSE)
dur.model2.RC  = lmer(Duration ~ Condition + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata[stim.lck.compdata$Channel=="lck_RCaMP",],REML=FALSE)
dur.anova.RC  <- anova(dur.null.RC, dur.model2.RC)
print(dur.anova.RC)

dur.Cond.RC<- glht(dur.model2.RC, mcp(Condition= "Tukey"))
summary(dur.Cond.RC)


#only consider STIM case
dur.null.stim = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="lck_RCaMP",],REML=FALSE)
dur.model1.stim  = lmer(Duration ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel!="lck_RCaMP",],REML=FALSE)

dur.anova.stim  <- anova(dur.null.stim, dur.model1.stim)
print(dur.anova.stim)

dur.stim.group<- glht(dur.model1.stim, mcp(Group= "Tukey"))
summary(dur.stim.group)

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





######
# distance to the nearest active neuron
# NOTE: no effect (see below)- distance is similar between groups and ROI types

respGC$minDistance=NA
for (ii in 1:nrow(respGC))
{
  ROIx=respGC$ROIs_trial[ii]
  subset1=subset(GCaMP_RCaMP.groups, ROIs_trial==ROIx)
  if (nrow(subset1)>0)
  {
    respGC$minDistance[ii]=min(subset1$MinDistance)
  }
}

dfDis<- summarySE(respGC, measurevar="minDistance", groupvars=c("Group"),na.rm = TRUE)
dfDis2<- summarySE(respGC, measurevar="minDistance", groupvars=c("ROIType"),na.rm = TRUE)
dfDis3<- summarySE(respGC, measurevar="minDistance", groupvars=c("ROIType","Group"),na.rm = TRUE)


ggplot(stim.lck.alldata[stim.lck.alldata$Channel_Group=="lck_GCaMP.fast",], aes(x=area, y=..density..,fill=Condition)) +
  geom_histogram(binwidth=5, position="dodge")+
  ylab("density") +
  xlim(-2,200)+
  ggtitle("lck data- fast no stim vs stim")+
  max.theme

Group_ROIType=interaction(respGC$Group, respGC$ROIType)
# stats for onset times- neurons vs astrocytes
minDis.lck.type.null = lmer(minDistance ~ (1|Animal) + (1|Spot) + (1|trials), respGC,REML=FALSE)
minDis.lck.type.model1 = lmer(minDistance ~ Group + (1|Animal) + (1|Spot) + (1|trials), respGC,REML=FALSE)
minDis.lck.type.model2 = lmer(minDistance ~ ROIType + (1|Animal) + (1|Spot) + (1|trials), respGC,REML=FALSE)
minDis.lck.type.model3 = lmer(minDistance ~ Group_ROIType + (1|Animal) + (1|Spot) + (1|trials),respGC,REML=FALSE)
minDis.lck.type.anova <- anova(minDis.lck.type.null, minDis.lck.type.model1, 
                               minDis.lck.type.model2,minDis.lck.type.model3)
print(minDis.lck.type.anova)





######


########
#mean peak amplitude for each ROI with peaks from the stim window
# use this to calculate high, mid and low responding neuronal somas

ROIwise<- ddply(all.lck.peaks, c("Animal", "Spot","Trial","Channel","roiName","Condition","trials", "ROIs_trial" ,"ROIType"), 
                summarise, meanAmp=mean(amplitude), meanDur= mean(Duration), meanProm=mean(prominence),
                meanPAUC=mean(peakAUC), nPeaks= length(amplitude))

ROIwise$roiNameUnique<-paste(ROIwise$Animal,ROIwise$Spot,ROIwise$roiName, sep="_")
ROIwise.Neurons<-subset(ROIwise, Condition!="Nostim" & ROIType=="Neuron")
ROIwise.Dendrites<-subset(ROIwise, Condition!="Nostim" & ROIType=="Dendrite")

# MEAN info for each ROI for trials of short OR long stim
NeuronSomas<-ddply(ROIwise.Neurons, c("Animal", "Spot", "roiNameUnique"),
                   summarise, meanAmp2=mean(meanAmp), meanProm2=mean(meanProm),
                   meanDur2=mean(meanDur), meanPAUC2=mean(meanPAUC))
NeuronDendrites<-ROIwise.Dendrites

# histograms for neurons
ggplot(NeuronSomas, aes(x=meanAmp2)) + geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("neuronal soma mean amplitude- both stim types averaged") +
  max.theme

ggplot(NeuronDendrites, aes(x=meanAmp)) + geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("neuronal dendrites mean amplitude- both stim types averaged") +
  max.theme


# calculate responders
neuronS_percentiles<-quantile(NeuronSomas$meanAmp2, prob = seq(0, 1, length = 21), type = 5)
neuronD_percentiles<-quantile(NeuronDendrites$meanAmp, prob = seq(0, 1, length = 21), type = 5)

NeuronSomas$responders=0
NeuronSomas$responders[NeuronSomas$meanAmp2>neuronS_percentiles[20]]="high"
NeuronSomas$responders[NeuronSomas$meanAmp2<neuronS_percentiles[11]]="low"
NeuronSomas$responders[NeuronSomas$meanAmp2<=neuronS_percentiles[20] & NeuronSomas$meanAmp2>=neuronS_percentiles[11]]="mid"

NeuronDendrites$responders=0
NeuronDendrites$responders[NeuronDendrites$meanAmp>neuronD_percentiles[20]]="high"
NeuronDendrites$responders[NeuronDendrites$meanAmp<neuronD_percentiles[11]]="low"
NeuronDendrites$responders[NeuronDendrites$meanAmp<=neuronD_percentiles[20] & NeuronDendrites$meanAmp>=neuronD_percentiles[11]]="mid"

#dendrites vs somas

df.amp.somas<- summarySE(data=NeuronSomas, measurevar = "meanAmp2", groupvars = c("responders"))
df.amp.dend<- summarySE(data=NeuronDendrites, measurevar = "meanAmp", groupvars = c("responders"))


#Add the responder type info to the peak table

# enter responding ROI information- high, mid or low responding neurons
all.lck.peaks.group$Nresponders=0
all.lck.peaks.Nresp=data.frame()
for (ii in 1:nrow(NeuronDendrites))
{
  Dend=NeuronDendrites$ROIs_trial[ii]
  RespType=NeuronDendrites$responders[ii]
  subset1=subset(all.lck.peaks.group, ROIs_trial==Dend)
  if (nrow(subset1)>0)
  {
    subset1$Nresponders=RespType
    all.lck.peaks.Nresp<-rbind(all.lck.peaks.Nresp, subset1)
  }
}

all.lck.peaks.group$roiNameUnique<-paste(all.lck.peaks.group$Animal,all.lck.peaks.group$Spot,all.lck.peaks.group$roiName, sep="_")
#GCaMP_RCaMP.groups$Nresponders=0
for (ii in 1:nrow(NeuronSomas))
{
  Soma=NeuronSomas$roiNameUnique[ii]
  RespType=NeuronSomas$responders[ii]
  subset1=subset(all.lck.peaks.group, roiNameUnique==Soma)
  if (nrow(subset1)>0)
  {
    subset1$Nresponders=RespType
    subset1$roiNameUnique=NULL
    all.lck.peaks.Nresp<-rbind(all.lck.peaks.Nresp, subset1)
  }
}

all.lck.peaks.group$roiNameUnique=NULL
all.lck.peaks.Nresp<-rbind(all.lck.peaks.Nresp, all.lck.peaks.group[all.lck.peaks.group$Channel=="GCaMP",])

###### 

#####

GCaMP$Condition <- factor(GCaMP$Condition, levels = c("shortstim","Stim"))

df.fast.amp1<- summarySE(GCaMP, measurevar = "amplitude", groupvars = c("Group"))
df.fast.amp2<- summarySE(GCaMP, measurevar = "amplitude", groupvars = c("Group","Condition"))
df.fast.amp3<- summarySE(GCaMP, measurevar = "amplitude", groupvars = c("Group","ROIType"))
df.fast.amp4<- summarySE(GCaMP, measurevar = "amplitude", groupvars = c("Group","Condition","ROIType"))


ggplot(data=df.fast.amp1, aes(x=Group, y= amplitude, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.fast.amp2, aes(x=Condition, y= amplitude, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.fast.amp4, aes(x=interaction(ROIType,Condition), y= amplitude, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


ggplot(data=df.fast.amp3, aes(x=ROIType, y= amplitude, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


GCaMP$Group_type<-interaction(GCaMP$Group, GCaMP$ROIType)
GCaMP$Group_type_Cond<-interaction(GCaMP$Group, GCaMP$ROIType, GCaMP$Condition)

group.amp.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP,REML=FALSE)
group.amp.model1 = lmer(amplitude  ~ Group + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP,REML=FALSE)
group.amp.model2= lmer(amplitude ~ Condition_group + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP, REML=FALSE)
group.amp.model3= lmer(amplitude ~ Group_type + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP, REML=FALSE)
group.amp.model4= lmer(amplitude ~ Group_type_Cond + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP, REML=FALSE)
group.amp.anova <- anova(group.amp.null, group.amp.model1, 
                         group.amp.model2,group.amp.model3,
                         group.amp.model4)
print(group.amp.anova)

group.amp.pv<- glht(group.amp.model1, mcp(Group= "Tukey"))
summary(group.amp.pv)

group.amp.cond.pv<- glht(group.amp.model2, mcp(Condition_group= "Tukey"))
summary(group.amp.cond.pv)

group.amp.Group_type.pv<- glht(group.amp.model3, mcp(Group_type= "Tukey"))
summary(group.amp.Group_type.pv)



#####
#duration

df.fast.dur1<- summarySE(GCaMP, measurevar = "Duration", groupvars = c("Group"))
df.fast.dur2<- summarySE(GCaMP, measurevar = "Duration", groupvars = c("Group","Condition"))
df.fast.dur3<- summarySE(GCaMP, measurevar = "Duration", groupvars = c("Group","ROIType"))
df.fast.dur4<- summarySE(GCaMP, measurevar = "Duration", groupvars = c("Group","Condition","ROIType"))


ggplot(data=df.fast.dur1, aes(x=Group, y= Duration, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.fast.dur2, aes(x=Condition, y= Duration, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.fast.dur4, aes(x=interaction(ROIType,Condition), y= Duration, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


ggplot(data=df.fast.dur3, aes(x=ROIType, y= Duration, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


group.dur.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP,REML=FALSE)
group.dur.model1 = lmer(Duration  ~ Group + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP,REML=FALSE)
group.dur.model2= lmer(Duration ~ Condition_group + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP, REML=FALSE)
group.dur.model3= lmer(Duration ~ Group_type + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP, REML=FALSE)
group.dur.model4= lmer(Duration ~ Group_type_Cond + (1|Animal) + (1|Spot) + (1|ROIs_trial), GCaMP, REML=FALSE)
group.dur.anova <- anova(group.dur.null, group.dur.model1, 
                         group.dur.model2,group.dur.model3,
                         group.dur.model4)
print(group.dur.anova)

group.dur.pv<- glht(group.dur.model1, mcp(Group= "Tukey"))
summary(group.dur.pv)

group.dur.cond.pv<- glht(group.dur.model2, mcp(Condition_group= "Tukey"))
summary(group.dur.cond.pv)

group.dur.Group_type.pv<- glht(group.dur.model3, mcp(Group_type= "Tukey"))
summary(group.dur.Group_type.pv)





