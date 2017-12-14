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
#Pharmacology
# Trazodone, Prazosin, Atropine, Meterogoline

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

all.lck.peaks <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/Peaks_pharmacology_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
all.lck.OT<-read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/OnsetTimes_pharmacology_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")

##### 
#home files

pharm.lck.peaks <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/Peaks_pharmacology_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
pharm.lck.OT<-read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/OnsetTimes_pharmacology_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")

control.lck.peaks <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/Peaks_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
control.lck.OT<-read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/OnsetTimes_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")

##########

lsm.options(pbkrtest.limit = 100000)


control.lck.peaks$Drug="Control"
control.lck.peaks<-subset(control.lck.peaks, Animal!="IPRG1")
control.lck.peaks<-subset(control.lck.peaks, Animal!="IPRG4")

control.lck.OT$ROIs_trial=NULL
control.lck.OT$Spot_trial=NULL
control.lck.OT$ROIType= NULL

#unique ROI names
all.lck.OT$ROIs_trial<-paste(all.lck.OT$Animal, all.lck.OT$Spot, all.lck.OT$Trial,all.lck.OT$ROI, sep= "_")
all.lck.OT$Spot_trial<-paste(all.lck.OT$Animal, all.lck.OT$Spot, all.lck.OT$Trial, sep= "_")
all.lck.OT$Spot_trial_Cond<-paste(all.lck.OT$Spot_trial, all.lck.OT$Condition, sep="_")
all.lck.OT$ROIs_Cond<-paste(all.lck.OT$ROIs_trial, all.lck.OT$Condition, sep="_")


# REMOVE duplicate entries from onset time and peak time data frames 
# only the first entry will be used
all.lck.OT2=all.lck.OT[order(all.lck.OT$OnsetTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.lck.OT<-distinct(all.lck.OT2, ROIs_Cond, .keep_all = TRUE)

#adjust onset time for the data with 2 s before stimulation included:

rm(all.lck.OT2)


all.lck.peaks<-rbind(control.lck.peaks, pharm.lck.peaks)

all.lck.peaks$ROIType= "none"
all.lck.peaksA<- subset(all.lck.peaks, Channel=="GCaMP")
all.lck.peaksB<- subset(all.lck.peaks, Channel=="RCaMP")

# ROITypes
all.lck.peaksA$ROIType[grepl("r",all.lck.peaksA$roiName)]="Process"
all.lck.peaksA$ROIType[grepl("E",all.lck.peaksA$roiName)]="Endfoot"
all.lck.peaksB$ROIType[grepl("r",all.lck.peaksB$roiName)]="Dendrite"
all.lck.peaksB$ROIType[grepl("D",all.lck.peaksB$roiName)]="Dendrite"
all.lck.peaksB$ROIType[grepl("N",all.lck.peaksB$roiName)]="Neuron"

all.lck.peaksB<-subset(all.lck.peaksB, ROIType!="none")
#all.lck.peaksB$ROIType[grepl("S",all.lck.peaksB$roiName)]="SmoothMuscle"

all.lck.peaks<-rbind(all.lck.peaksA, all.lck.peaksB)
rm(all.lck.peaksA, all.lck.peaksB)
all.lck.peaks$ROIType<- as.factor(all.lck.peaks$ROIType)

#unique ROI names
all.lck.peaks$ROIs_trial<-paste(all.lck.peaks$Animal, all.lck.peaks$Spot, all.lck.peaks$Trial,all.lck.peaks$roiName, sep= "_")
all.lck.peaks$trials<-paste(all.lck.peaks$Animal, all.lck.peaks$Spot, all.lck.peaks$Trial, sep= "_")
all.lck.peaks$trials_Cond<-paste(all.lck.peaks$trials, all.lck.peaks$Condition, sep= "_")

all.lck.peaks$ROIs_Cond<-paste(all.lck.peaks$ROIs_trial, all.lck.peaks$Condition, sep= "_")


# count number of trials per spot
Spot.lck.ntrials<-ddply(all.lck.peaks, c("Animal","Drug","Spot","Condition"), summarise, nTrials=length(unique(Trial)))


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

all.lck.OT$Drug<- factor(all.lck.OT$Drug, levels=c("Control","Atropine","Prazosin","Trazodone"))
all.lck.peaks$Drug<- factor(all.lck.peaks$Drug, levels=c("Control","Atropine","Prazosin","Trazodone"))



# neuronal responses to stimulation
NeuronalStim<-subset(all.lck.OT, Channel=="RCaMP" & Condition=="Stim")

# should have an onset time in 8 s stimulus
ggplot(NeuronalStim[NeuronalStim$OnsetTime<10,],aes(x=OnsetTime,y=..density..,fill=Drug)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP onset times between 0 and 10 s from stim trials")+
  max.theme


Neuron95Onset<-quantile(NeuronalStim$OnsetTime[NeuronalStim$OnsetTime<8], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
print(Neuron95Onset)


AstroStim<-subset(all.lck.OT, Channel=="GCaMP" & Condition=="Stim")

# should have an onset time in 8 s stimulus
ggplot(AstroStim[AstroStim$OnsetTime<15,],aes(x=OnsetTime,fill=Drug)) +
  geom_histogram(binwidth=(0.084*5), position="dodge") +
  ggtitle("Lck-GCaMP onset times between 0 and 15 s from stim trials")+
  max.theme

# 
# Onset time histograms- normalized to the number of trials

stimwindow=15

GC.lck.OT.dist<-subset(all.lck.OT,Condition=="Stim" & OnsetTime<stimwindow & Channel=="GCaMP")

# KO vs WT stim
ntrials.GC.KOvsWT<- ddply(GC.lck.OT.dist, c("Drug"), summarise, ntrials=length(unique(Spot_trial)))


#histogram bins
histseq= seq(0,15,0.5)
KO.A=0
WT.A=0
zeroRow<-data.frame(cbind(KO.A, WT.A))

# neuronal lck onset histogram
# counts for each condition in the histogram
KO.A=hist(GC.lck.OT.dist$OnsetTime[GC.lck.OT.dist$Drug=="IP3R2_KO"], breaks=histseq, plot=FALSE)$counts
WT.A=hist(GC.lck.OT.dist$OnsetTime[GC.lck.OT.dist$Drug=="IP3R2_WT"], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
KO.A=KO.A/ntrials.GC.KOvsWT$ntrials[(ntrials.GC.KOvsWT$Drug=="IP3R2_KO")]
WT.A=WT.A/ntrials.GC.KOvsWT$ntrials[(ntrials.GC.KOvsWT$Drug=="IP3R2_WT")]

#make a data frame for plotting
lck.histo <- data.frame(cbind(KO.A, WT.A))
lck.histo2<-rbind(zeroRow,lck.histo)
lck.histo2$time<-histseq

ggplot(NULL, aes(x=time))+
  geom_line(data=lck.histo2, aes(y=KO.A, color="KO.A")) +
  geom_line(data=lck.histo2, aes(y=WT.A, color="WT.A")) +
  ggtitle("IP3R2 KO vs WT astrocytes-lck data") + 
  xlab("Onset Time (s)") + 
  ylab("Astrocytes peaks/trial") + 
  xlim(-0.5,15) +
  max.theme




#compare distributions (what is plotted in the figure)

# ks test- astrocytes
OT.lck.GC.KOvsWT.kstest<- ks.test(KO.A,WT.A)
print(OT.lck.GC.KOvsWT.kstest)




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

rm(stim.lck.OT.G,stim.lck.OT.R,stim.lck.PT.G, stim.lck.PT.R)
#######
# peak times only for ROIs with onset times?

stim.lck.alldata<-merge(stim.lck.peaks.window, stim.lck.OT.window[, c("ROIs_Cond", "OnsetTime","TraceAUC1","TraceAUC10")], by="ROIs_Cond")


# onset times
ggplot(stim.lck.alldata, aes(x=interaction(Channel, Drug),y=OnsetTime, fill= Condition)) +
  geom_boxplot(notch=TRUE)+
  ylab("OnsetTimes Time (s)") +
  ggtitle("matched data= peak times and onset times")+
  max.theme

#peak times
ggplot(stim.lck.alldata, aes(x=interaction(Channel, Drug),y=peakTime, fill= Condition)) +
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

GCaMP.KO<-subset(stim.lck.alldata,Channel=="GCaMP" & Drug=="IP3R2_KO" & Condition=="Stim")
GCaMP.WT<-subset(stim.lck.alldata,Channel=="GCaMP" & Drug=="IP3R2_WT" & Condition=="Stim")
RCaMP.KO<-subset(stim.lck.alldata,Channel=="RCaMP" & Drug=="IP3R2_KO" & Condition=="Stim")
RCaMP.WT<-subset(stim.lck.alldata,Channel=="RCaMP" & Drug=="IP3R2_WT" & Condition=="Stim")

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

Spot.lck.stim<-ddply(stim.lck.alldata, c("Animal","Spot","Drug","Condition","Channel","nFluoPix","nActivePix","pixelsize"), summarise,
                     nROIs=length(unique(ROIs_Cond)))

# add in number of trials
Spot.lck.stim$Ani_Spot_Cond<-paste(Spot.lck.stim$Animal, Spot.lck.stim$Spot, Spot.lck.stim$Condition, sep="_")
Spot.lck.ntrials$Ani_Spot_Cond<-paste(Spot.lck.ntrials$Animal, Spot.lck.ntrials$Spot, Spot.lck.ntrials$Condition, sep="_")

Spot.lck.stim<-merge(Spot.lck.stim, Spot.lck.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)

Spot.lck.stim$ROIsPerTrial<-Spot.lck.stim$nROIs/Spot.lck.stim$nTrials

Spot.lck.stim$FracActive=Spot.lck.stim$nActivePix/Spot.lck.stim$nFluoPix
Spot.lck.stim$FracActivePerROI=Spot.lck.stim$FracActive/Spot.lck.stim$nROIs

# mean fraction of active pixels
df.FracActive<- summarySE(Spot.lck.stim[Spot.lck.stim$Channel=="GCaMP",], measurevar = "FracActive", groupvars = c("Condition", "Drug"))
df.FracActivePerROI<- summarySE(Spot.lck.stim[Spot.lck.stim$Channel=="GCaMP",], measurevar = "FracActivePerROI", groupvars = c("Condition", "Drug"))

ggplot(data=df.FracActive, aes(x=Condition, y= FracActive, fill=Drug)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=FracActive-se, ymax=FracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.FracActivePerROI, aes(x=Condition, y= FracActivePerROI, fill=Drug)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=FracActivePerROI-se, ymax=FracActivePerROI+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive Per ROI") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

#only astrocytes
Cond_Drug= interaction(Spot.lck.stim$Condition[Spot.lck.stim$Channel=="GCaMP"],Spot.lck.stim$Drug[Spot.lck.stim$Channel=="GCaMP"])
FracActive.lck.stim.null = lmer(FracActive ~ (1|Animal), Spot.lck.stim[Spot.lck.stim$Channel=="GCaMP",],REML=FALSE)
FracActive.lck.stim.model1 = lmer(FracActive~ Condition+ (1|Animal), Spot.lck.stim[Spot.lck.stim$Channel=="GCaMP",],REML=FALSE)
FracActive.lck.stim.model2 = lmer(FracActive ~ Drug + (1|Animal), Spot.lck.stim[Spot.lck.stim$Channel=="GCaMP",],REML=FALSE)
FracActive.lck.stim.model3 = lmer(FracActive ~ condition + Drug + (1|Animal), Spot.lck.stim[Spot.lck.stim$Channel=="GCaMP",],REML=FALSE)
FracActive.lck.stim.model4 = lmer(FracActive ~ Condition_Drug + (1|Animal), Spot.lck.stim[Spot.lck.stim$Channel=="GCaMP",],REML=FALSE)
FracActive.lck.stim.anova <- anova(FracActive.lck.stim.null, FracActive.lck.stim.model1,FracActive.lck.stim.model2,
                                 FracActive.lck.stim.model3,FracActive.lck.stim.model4)
print(FracActive.lck.stim.anova)

FracActive.lck.stim.Cond_Drug<- glht(FracActive.lck.stim.model4, mcp(Cond_Drug= "Tukey"))
summary(FracActive.lck.stim.Cond_Drug)

####
# mean number of total ROIs per trial
df.lck.ROInum.mean<-summarySE(Spot.lck.stim, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Drug","Condition"))


ggplot(df.lck.ROInum.mean, aes(x=interaction(Drug,Channel),y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view") +
  max.theme

Condition_Channel2= interaction(Spot.lck.stim$Condition,Spot.lck.stim$Channel)
Condition_Channel_Drug= interaction(Spot.lck.stim$Condition,Spot.lck.stim$Channel,Spot.lck.stim$Drug)
nROI.lck.stim.null = lmer(ROIsPerTrial ~ (1|Animal), Spot.lck.stim,REML=FALSE)
nROI.lck.stim.model1 = lmer(ROIsPerTrial~ Channel + (1|Animal), Spot.lck.stim,REML=FALSE)
nROI.lck.stim.model2 = lmer(ROIsPerTrial ~ Condition + (1|Animal), Spot.lck.stim,REML=FALSE)
nROI.lck.stim.model3 = lmer(ROIsPerTrial ~ Drug + (1|Animal), Spot.lck.stim,REML=FALSE)
nROI.lck.stim.model4 = lmer(ROIsPerTrial ~ Condition_Channel2 + (1|Animal), Spot.lck.stim,REML=FALSE)
nROI.lck.stim.model5 = lmer(ROIsPerTrial ~ Condition_Channel2 + Drug + (1|Animal), Spot.lck.stim,REML=FALSE)
nROI.lck.stim.model6 = lmer(ROIsPerTrial ~ Condition_Channel_Drug + (1|Animal), Spot.lck.stim,REML=FALSE)
nROI.lck.stim.anova <- anova(nROI.lck.stim.null, nROI.lck.stim.model1,nROI.lck.stim.model2,nROI.lck.stim.model3,
                             nROI.lck.stim.model4,nROI.lck.stim.model5,nROI.lck.stim.model6)
print(nROI.lck.stim.anova)

nROI.lck.stim.Cond_Channel_Drug<- glht(nROI.lck.stim.model6, mcp(Condition_Channel_Drug= "Tukey"))
summary(nROI.lck.stim.Cond_Channel_Drug)



# number of fastROIs per trial
FastROInum.lck.stim<-ddply(stim.lck.alldata[stim.lck.alldata$Condition=="Stim",], c("Animal","Spot","Drug","Condition","Channel","Group","Channel_Group"), summarise, nROIs=length(OnsetTime))

# add in number of trials
FastROInum.lck.stim$Ani_Spot_Cond<-paste(FastROInum.lck.stim$Animal, FastROInum.lck.stim$Spot, FastROInum.lck.stim$Condition, sep="_")

FastROInum.lck.stim<-merge(FastROInum.lck.stim, Spot.lck.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)

FastROInum.lck.stim$ROIsPerTrial<-FastROInum.lck.stim$nROIs/FastROInum.lck.stim$nTrials

# mean number of total ROIs per trial
df.lck.FastROInum.mean<-summarySE(FastROInum.lck.stim, measurevar = "ROIsPerTrial", groupvars = c("Channel_Group", "Drug"))
df.lck.FastROInum.GC<-summarySE(FastROInum.lck.stim[FastROInum.lck.stim$Channel=="GCaMP",], measurevar = "ROIsPerTrial", groupvars = c("Group", "Drug"))


ggplot(df.lck.FastROInum.mean, aes(x=Channel_Group,y=ROIsPerTrial, fill= Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during stimulation") +
  max.theme

ggplot(df.lck.FastROInum.GC, aes(x=Group,y=ROIsPerTrial, fill= Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during stimulation") +
  max.theme

#only astrocytes
Group_Drug= interaction(FastROInum.lck.stim$Group[FastROInum.lck.stim$Channel=="GCaMP"],FastROInum.lck.stim$Drug[FastROInum.lck.stim$Channel=="GCaMP"])
FastnROI.lck.stim.null = lmer(ROIsPerTrial ~ (1|Animal), FastROInum.lck.stim[FastROInum.lck.stim$Channel=="GCaMP",],REML=FALSE)
FastnROI.lck.stim.model1 = lmer(ROIsPerTrial~ Group + (1|Animal), FastROInum.lck.stim[FastROInum.lck.stim$Channel=="GCaMP",],REML=FALSE)
FastnROI.lck.stim.model2 = lmer(ROIsPerTrial ~ Drug + (1|Animal), FastROInum.lck.stim[FastROInum.lck.stim$Channel=="GCaMP",],REML=FALSE)
FastnROI.lck.stim.model3 = lmer(ROIsPerTrial ~ Group + Drug + (1|Animal), FastROInum.lck.stim[FastROInum.lck.stim$Channel=="GCaMP",],REML=FALSE)
FastnROI.lck.stim.model4 = lmer(ROIsPerTrial ~ Group_Drug + (1|Animal), FastROInum.lck.stim[FastROInum.lck.stim$Channel=="GCaMP",],REML=FALSE)
FastnROI.lck.stim.anova <- anova(FastnROI.lck.stim.null, FastnROI.lck.stim.model1,FastnROI.lck.stim.model2,
                                 FastnROI.lck.stim.model3,FastnROI.lck.stim.model4)
print(FastnROI.lck.stim.anova)

FastnROI.lck.stim.Group_Drug<- glht(FastnROI.lck.stim.model4, mcp(Group_Drug= "Tukey"))
summary(FastnROI.lck.stim.Group_Drug)


#########
# mean onset times
df.OT1<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim",], measurevar = "OnsetTime", groupvars = c("Channel_Group", "Drug"))
df.OT2<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",], measurevar = "OnsetTime", groupvars = c("Group", "Drug"))
df.OT3<- summarySE(stim.lck.compdata.STIM, measurevar = "OnsetTime", groupvars = c("ROIType","Channel_Group", "Drug"))

ggplot(df.OT1, aes(x=interaction(Drug,Channel_Group),y=OnsetTime, fill=Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme

ggplot(df.OT2, aes(x=Group,y=OnsetTime, fill=Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme

ggplot(df.OT3, aes(x=interaction(Channel_Group,ROIType),y=OnsetTime, fill=Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme



# stats
Group_Channel_Type_Cond_Drug=interaction(stim.lck.alldata$Group,stim.lck.alldata$Channel,stim.lck.alldata$ROIType, 
                                        stim.lck.alldata$Drug,stim.lck.alldata$Condition)
Group_Channel_Cond_Drug=interaction(stim.lck.alldata$Group,stim.lck.alldata$Channel,stim.lck.alldata$Drug,
                                   stim.lck.alldata$Condition)

# stats for onset times- neurons vs astrocytes
OT.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model2 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model3 = lmer(OnsetTime ~ Group_Channel_Cond_Drug + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.model4 = lmer(OnsetTime ~ Group_Channel_Type_Cond_Drug + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
OT.anova <- anova(OT.null, OT.model1,OT.model2,OT.model3,OT.model4)
print(OT.anova)

OT.Group_channel_Drug<- glht(OT.model3, mcp(Group_Channel_Cond_Drug= "Tukey"))
summary(OT.Group_channel_Drug)

OT.Group_Channel_Type_Drug<- glht(OT.model4, mcp(Group_Channel_Type_Cond_Drug= "Tukey"))
summary(OT.Group_Channel_Type_Drug)


Group_Channel_Type_Drug=interaction(stim.lck.compdata.STIM$Group,stim.lck.compdata.STIM$Channel,
                                   stim.lck.compdata.STIM$ROIType, stim.lck.compdata.STIM$Drug)
Group_Channel_Drug=interaction(stim.lck.compdata.STIM$Channel, stim.lck.compdata.STIM$Group, stim.lck.compdata.STIM$Drug)

# stats for onset times- neurons vs astrocytes
OT.stim.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.model3 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.model4 = lmer(OnsetTime ~ Group_Channel_Drug + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.model6 = lmer(OnsetTime ~ Group_Channel_Type_Drug + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.anova <- anova(OT.stim.null, OT.stim.model1,OT.stim.model3,OT.stim.model4,OT.stim.model6)
print(OT.stim.anova)

OT.stim.Group_channel_Drug<- glht(OT.stim.model4, mcp(Group_Channel_Drug= "Tukey"))
summary(OT.stim.Group_channel_Drug)

OT.stim.Group_channel_type_Drug<- glht(OT.stim.model6, mcp(Group_Channel_Type_Drug= "Tukey"))
summary(OT.stim.Group_channel_type_Drug)

summary(OT.stim.model4)

# check residuals for linearity
plot(fitted(OT.stim.model4), residuals(OT.stim.model4),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(OT.stim.model4), residuals(OT.stim.model4)), col=46, lwd=2.5)





#############
# mean peak time
# mean onset times
df.PT1<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim",], measurevar = "peakTime", groupvars = c("Channel_Group", "Drug"))
df.PT2<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",], measurevar = "peakTime", groupvars = c("Group", "Drug"))
df.PT3<- summarySE(stim.lck.compdata.STIM, measurevar = "peakTime", groupvars = c("ROIType","Channel_Group", "Drug"))

ggplot(df.PT1, aes(x=interaction(Drug,Channel_Group),y=peakTime, fill=Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Peak Max Time (s)") +
  max.theme

ggplot(df.PT2, aes(x=Group,y=peakTime, fill=Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Peak Max Time (s)") +
  max.theme

ggplot(df.PT3, aes(x=interaction(Channel_Group,ROIType),y=peakTime, fill=Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Peak Max Time (s)") +
  max.theme



# stats
# stats for onset times- neurons vs astrocytes
PT.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
PT.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
PT.model2 = lmer(peakTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
PT.model3 = lmer(peakTime ~ Group_Channel_Cond_Drug + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
PT.model4 = lmer(peakTime ~ Group_Channel_Type_Cond_Drug + (1|Animal) + (1|Spot) + (1|trials), stim.lck.alldata,REML=FALSE)
PT.anova <- anova(PT.null, PT.model1,PT.model2,PT.model3,PT.model4)
print(PT.anova)

PT.Group_channel_Drug<- glht(PT.model3, mcp(Group_Channel_Drug= "Tukey"))
summary(PT.Group_channel_Drug)

PT.Group_Channel_Type_Drug<- glht(PT.model4, mcp(Group_Channel_Type_Drug= "Tukey"))
summary(PT.Group_Channel_Type_Drug)


Group_Channel_Type_Drug=interaction(stim.lck.compdata.STIM$Group,stim.lck.compdata.STIM$Channel,
                                   stim.lck.compdata.STIM$ROIType, stim.lck.compdata.STIM$Drug)
Group_Channel_Drug=interaction(stim.lck.compdata.STIM$Channel, stim.lck.compdata.STIM$Group, stim.lck.compdata.STIM$Drug)

# stats for onset times- neurons vs astrocytes
PT.stim.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
PT.stim.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
PT.stim.model3 = lmer(peakTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
PT.stim.model4 = lmer(peakTime ~ Group_Channel_Drug + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
PT.stim.model6 = lmer(peakTime ~ Group_Channel_Type_Drug + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM,REML=FALSE)
PT.stim.anova <- anova(PT.stim.null, PT.stim.model1,PT.stim.model3,PT.stim.model4,PT.stim.model6)
print(PT.stim.anova)

PT.stim.Group_channel_Drug<- glht(PT.stim.model4, mcp(Group_Channel_Drug= "Tukey"))
summary(PT.stim.Group_channel_Drug)

PT.stim.Group_channel_type_Drug<- glht(PT.stim.model6, mcp(Group_Channel_Type_Drug= "Tukey"))
summary(PT.stim.Group_channel_type_Drug)


########
#amplitude
df.amp1<-summarySE(stim.lck.alldata, measurevar = "amplitude", groupvars = c("Channel","Condition", "Drug"))
df.amp2<-summarySE(stim.lck.alldata, measurevar = "amplitude", groupvars = c("Channel", "ROIType","Condition", "Drug"))

df.amp3A<- summarySE(stim.lck.compdata, measurevar = "amplitude", groupvars = c("Channel_Group","Drug","Condition"))
df.amp3B<- summarySE(stim.lck.compdata.STIM, measurevar = "amplitude", groupvars = c("Channel_Group","Drug"))

df.amp4<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("Group","Drug"))
df.amp5<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",], measurevar = "amplitude", groupvars = c("Group","ROIType","Drug"))


ggplot(df.amp1, aes(x=interaction(Channel,Drug),y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  max.theme

ggplot(df.amp2, aes(x=interaction(Channel,interaction(Drug, ROIType)),y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  max.theme

ggplot(df.amp3A, aes(x=interaction(Channel_Group,Drug),y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  max.theme

ggplot(df.amp3B, aes(x=Channel_Group,y=amplitude, fill= Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude during stim trials") +
  max.theme

ggplot(df.amp4, aes(x=Group,y=amplitude, fill= Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude during stim trials") +
  max.theme

ggplot(df.amp5, aes(x=interaction(Group, ROIType),y=amplitude, fill= Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude during stim trials") +
  max.theme





#lck-GCaMP ONLY
Group_Drug=interaction(stim.lck.compdata.STIM$Group[stim.lck.compdata.STIM$Channel=="GCaMP"],stim.lck.compdata.STIM$Drug[stim.lck.compdata.STIM$Channel=="GCaMP"])
Group_Type_Drug=interaction(stim.lck.compdata.STIM$Group[stim.lck.compdata.STIM$Channel=="GCaMP"],
                               stim.lck.compdata.STIM$Drug[stim.lck.compdata.STIM$Channel=="GCaMP"],
                               stim.lck.compdata.STIM$ROIType[stim.lck.compdata.STIM$Channel=="GCaMP"])

amp.null.GC = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
amp.model2.GC  = lmer(amplitude ~ Drug + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
amp.model3.GC  = lmer(amplitude ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
amp.model5.GC  = lmer(amplitude ~ Group_Drug + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
amp.model6.GC  = lmer(amplitude ~ Group_Type_Drug + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
amp.anova.GC  <- anova(amp.null.GC, amp.model2.GC,amp.model3.GC,amp.model5.GC,amp.model6.GC)
print(amp.anova.GC)

amp.Group_Drug.GC<- glht(amp.model5.GC, mcp(Group_Drug= "Tukey"))
summary(amp.Group_Drug.GC)

amp.Group_Drug_ty.GC<- glht(amp.model6.GC, mcp(Group_Type_Drug= "Tukey"))
summary(amp.Group_Drug_ty.GC)


########
#duration
df.Dur1<-summarySE(stim.lck.alldata, measurevar = "Duration", groupvars = c("Channel","Condition", "Drug"))
df.Dur2<-summarySE(stim.lck.alldata, measurevar = "Duration", groupvars = c("Channel", "ROIType","Condition", "Drug"))

df.Dur3A<- summarySE(stim.lck.compdata, measurevar = "Duration", groupvars = c("Channel_Group","Drug","Condition"))
df.Dur3B<- summarySE(stim.lck.compdata.STIM, measurevar = "Duration", groupvars = c("Channel_Group","Drug"))

df.Dur4<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Group","Drug"))
df.Dur5<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",], measurevar = "Duration", groupvars = c("Group","ROIType","Drug"))


ggplot(df.Dur1, aes(x=interaction(Channel,Drug),y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  max.theme

ggplot(df.Dur2, aes(x=interaction(Channel,interaction(Drug, ROIType)),y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  max.theme

ggplot(df.Dur3A, aes(x=interaction(Channel_Group,Drug),y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  max.theme

ggplot(df.Dur3B, aes(x=Channel_Group,y=Duration, fill= Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration during stim trials") +
  max.theme

ggplot(df.Dur4, aes(x=Group,y=Duration, fill= Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration during stim trials") +
  max.theme

ggplot(df.Dur5, aes(x=interaction(Group, ROIType),y=Duration, fill= Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration during stim trials") +
  max.theme





#lck-GCaMP ONLY
Group_Drug=interaction(stim.lck.compdata.STIM$Group[stim.lck.compdata.STIM$Channel=="GCaMP"],stim.lck.compdata.STIM$Drug[stim.lck.compdata.STIM$Channel=="GCaMP"])
Group_Type_Drug=interaction(stim.lck.compdata.STIM$Group[stim.lck.compdata.STIM$Channel=="GCaMP"],
                                stim.lck.compdata.STIM$Drug[stim.lck.compdata.STIM$Channel=="GCaMP"],
                                stim.lck.compdata.STIM$ROIType[stim.lck.compdata.STIM$Channel=="GCaMP"])

Dur.null.GC = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
Dur.model2.GC  = lmer(Duration ~ Drug + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
Dur.model3.GC  = lmer(Duration ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
Dur.model5.GC  = lmer(Duration ~ Group_Drug + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
Dur.model6.GC  = lmer(Duration ~ Group_Type_Drug + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
Dur.anova.GC  <- anova(Dur.null.GC, Dur.model2.GC,Dur.model3.GC,Dur.model5.GC,Dur.model6.GC)
print(Dur.anova.GC)

Dur.Group_Drug.GC<- glht(Dur.model5.GC, mcp(Group_Drug= "Tukey"))
summary(Dur.Group_Drug.GC)

Dur.Group_Drug_ty.GC<- glht(Dur.model6.GC, mcp(Group_Type_Drug= "Tukey"))
summary(Dur.Group_Drug_ty.GC)

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



