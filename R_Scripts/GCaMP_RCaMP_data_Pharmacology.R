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
# Trazodone, Prazosin, Atropine, Metergoline, DSP4

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

pharm.lck.peaks <- read.table("D:/Data/GCaMP_RCaMP/Revision/Lck_GCaMP/FilesforR/Pharmacology/Peaks_pharmacology_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
pharm.lck.OT<-read.table("D:/Data/GCaMP_RCaMP/Revision/Lck_GCaMP/FilesforR/Pharmacology/OnsetTimes_pharmacology_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")

lck.peaks1 <- read.table("D:/Data/GCaMP_RCaMP/Revision/Lck_GCaMP/FilesforR/Control_untreated/Peaks_1stCohort_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
lck.OT1<-read.table("D:/Data/GCaMP_RCaMP/Revision/Lck_GCaMP/FilesforR/Control_untreated/OnsetTimes_1stCohort_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
lck.peaks2 <- read.table("D:/Data/GCaMP_RCaMP/Revision/Lck_GCaMP/FilesforR/Control_untreated/Peaks_2ndCohort_Lck_nostim_vs_longstim_01_2018.csv", header=TRUE, sep = ",")
lck.OT2<-read.table("D:/Data/GCaMP_RCaMP/Revision/Lck_GCaMP/FilesforR/Control_untreated/OnsetTimes_2ndCohort_Lck_nostim_vs_longstim_01_2018.csv", header=TRUE, sep = ",")

##########
lsm.options(pbkrtest.limit = 100000)

#  ONLY CONSIDER DATA FROM THE SAME SPOTS IN CONTROL

# remove the DMSO_control data for now
pharm.lck.OT<-subset(pharm.lck.OT, Drug!="DMSO_Control")
pharm.lck.peaks<-subset(pharm.lck.peaks, Drug!="DMSO_Control")

lck.peaks1$Drug="Control"
lck.peaks2$Drug="Control"

# reshape data tables (remove data we don't need)
lck.peaks2$Genotype<-NULL
lck.OT2$Genotype<-NULL
lck.OT1$Spot_trial<-NULL
lck.OT1$ROIs_trial<-NULL
lck.OT1$ROIType<-NULL
names(lck.OT1)[names(lck.OT1)=="Overlap"] <- "overlap"
lck.OT1$Drug="Control"
lck.OT2$Drug="Control"

control.OT<-rbind(lck.OT1, lck.OT2)
control.peaks<-rbind(lck.peaks1, lck.peaks2)

rm(lck.OT1, lck.OT2, lck.peaks1, lck.peaks2)

# get a list of spots that are found in the pharamcology table
pharm.lck.OT$Animal_Spot<-paste(pharm.lck.OT$Animal, pharm.lck.OT$Spot, sep="_")
pharm.lck.peaks$Animal_Spot<-paste(pharm.lck.peaks$Animal, pharm.lck.peaks$Spot, sep="_")
control.OT$Animal_Spot<-paste(control.OT$Animal, control.OT$Spot, sep="_")
control.peaks$Animal_Spot<-paste(control.peaks$Animal, control.peaks$Spot, sep="_")

SpotList= unique(pharm.lck.OT$Animal_Spot)

control.OT<-subset(control.OT, Animal_Spot %in% SpotList)
control.peaks<-subset(control.peaks, Animal_Spot %in% SpotList)

all.lck.OT<-rbind(control.OT, pharm.lck.OT)
all.lck.peaks<-rbind(control.peaks, pharm.lck.peaks)

rm(control.OT, pharm.lck.OT, control.peaks, pharm.lck.peaks)

#unique ROI names
all.lck.OT$ROIs_trial<-paste(all.lck.OT$Animal, all.lck.OT$Spot, all.lck.OT$Trial,all.lck.OT$ROI, sep= "_")
all.lck.OT$Spot_trial_Drug<-paste(all.lck.OT$Animal, all.lck.OT$Spot, all.lck.OT$Trial,all.lck.OT$Drug, sep= "_")
all.lck.OT$Spot_trial_Drug_Cond<-paste(all.lck.OT$Spot_trial_Drug, all.lck.OT$Condition, sep="_")
all.lck.OT$ROIs_Cond_Drug<-paste(all.lck.OT$ROIs_trial, all.lck.OT$Condition, all.lck.OT$Drug, sep="_")


# REMOVE duplicate entries from onset time and peak time data frames 
# only the first entry will be used
all.lck.OT2=all.lck.OT[order(all.lck.OT$OnsetTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.lck.OT<-distinct(all.lck.OT2, ROIs_Cond_Drug, .keep_all = TRUE)

#remove entries with no onset time
all.lck.OT<-subset(all.lck.OT, OnsetTime!="NaN")

all.lck.OT$ROIType= "none"
all.lck.OTA<- subset(all.lck.OT, Channel=="GCaMP")
all.lck.OTB<- subset(all.lck.OT, Channel=="RCaMP")

# ROITypes
all.lck.OTA$ROIType[grepl("r",all.lck.OTA$ROI)]="Process"
all.lck.OTA$ROIType[grepl("E",all.lck.OTA$ROI)]="Endfoot"
all.lck.OTB$ROIType[grepl("r",all.lck.OTB$ROI)]="Dendrite"
all.lck.OTB$ROIType[grepl("D",all.lck.OTB$ROI)]="Dendrite"
all.lck.OTB$ROIType[grepl("N",all.lck.OTB$ROI)]="Neuron"

all.lck.OTB<-subset(all.lck.OTB, ROIType!="none")

all.lck.OT<-rbind(all.lck.OTA, all.lck.OTB)
all.lck.OT$ROIType<- as.factor(all.lck.OT$ROIType)

rm(all.lck.OT2, all.lck.OTA, all.lck.OTB)


# peak data
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

all.lck.peaks<-rbind(all.lck.peaksA, all.lck.peaksB)
rm(all.lck.peaksA, all.lck.peaksB)
all.lck.peaks$ROIType<- as.factor(all.lck.peaks$ROIType)

#unique ROI names
all.lck.peaks$ROIs_trial<-paste(all.lck.peaks$Animal, all.lck.peaks$Spot, all.lck.peaks$Trial,all.lck.peaks$roiName, sep= "_")
all.lck.peaks$Spot_trial_Drug<-paste(all.lck.peaks$Animal, all.lck.peaks$Spot, all.lck.peaks$Trial,all.lck.peaks$Drug, sep= "_")
all.lck.peaks$Spot_trial_Drug_Cond<-paste(all.lck.peaks$Spot_trial_Drug, all.lck.peaks$Condition, sep="_")
all.lck.peaks$ROIs_Cond_Drug<-paste(all.lck.peaks$ROIs_trial, all.lck.peaks$Condition, all.lck.peaks$Drug, sep="_")

# count number of trials per spot
Spot.lck.ntrials<-ddply(all.lck.peaks, c("Animal","Drug","Spot","Condition"), summarise, nTrials=length(unique(Trial)))


# remove ROIs with no peaks
all.lck.peaks<-all.lck.peaks[!(all.lck.peaks$peakType=="NoPeak"),]


# remove matching astrocyte process and soma ROIs
Overlap= all.lck.peaks$overlap!=0
all.lck.peaks<-all.lck.peaks[!Overlap,]
#OverlapROIs<-unique(nostim$ROIs_trial[Overlap])


## COMBINE ONSET TIME AND PEAK TIME DATA TABLES
#add baseline time, onset time, trace auc to peaks table
all.lck.OT$BL_time<-all.lck.OT$baseline/all.lck.OT$FrameRate
all.lck.peaks<-merge(all.lck.peaks, all.lck.OT[, c("ROIs_Cond_Drug", "BL_time", "OnsetTime", "TraceAUC1", "TraceAUC10")], by="ROIs_Cond_Drug", all.x=TRUE)

#only consider peaks with a valid onset time
all.lck.peaks<-subset(all.lck.peaks, complete.cases(all.lck.peaks$BL_time))

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
all.lck.peaks<-distinct(all.lck.peaks3, ROIs_Cond_Drug,.keep_all = TRUE)


all.lck.peaks$Drug<- as.factor(all.lck.peaks$Drug)
all.lck.peaks$Drug<- factor(all.lck.peaks$Drug, levels=c("Control","Atropine","Metergoline","Trazodone","Prazosin","DSP4"))


rm(all.lck.peaks2, all.lck.peaks3, all.lck.OT)



## DSP4 data

DSP4<- subset(all.lck.peaks, Drug=="DSP4")
DSP4_spots<-unique(DSP4$Animal_Spot)
control.DSP4<-subset(all.lck.peaks, (Animal_Spot %in% DSP4_spots) & Drug=="Control")

DSP4.data<-rbind(DSP4, control.DSP4)

all.lck.peaks<-subset(all.lck.peaks, Drug!="DSP4")
all.lck.peaks<-subset(all.lck.peaks, Animal!="IPRG6")

rm(DSP4,control.DSP4)


######
# does neuronal spontaneous activity change?- at any point in nostim trials
NeuronalNostim<-subset(all.lck.peaks, Channel=="RCaMP" & Condition=="Nostim")

Neur.NS.Spot<-ddply(NeuronalNostim, c("Animal","Spot","Drug","Condition", "Animal_Spot"), summarise, nROIs=length(OnsetTime),
                 meanAmp=mean(amplitude))
Neur.NS.Spot$Ani_Spot_Cond_Drug<-paste(Neur.NS.Spot$Animal_Spot, Neur.NS.Spot$Condition, Neur.NS.Spot$Drug, sep="_")
Spot.lck.ntrials$Ani_Spot_Cond_Drug<-paste(Spot.lck.ntrials$Animal, Spot.lck.ntrials$Spot,Spot.lck.ntrials$Condition,Spot.lck.ntrials$Drug, sep="_")
Neur.NS.Spot<-merge(Neur.NS.Spot, Spot.lck.ntrials[, c("Ani_Spot_Cond_Drug", "nTrials")], by="Ani_Spot_Cond_Drug", all.x=TRUE)
Neur.NS.Spot$ROIsPerTrial<-Neur.NS.Spot$nROIs/Neur.NS.Spot$nTrials

#DSP4
NeuronalNostim.DSP4<-subset(DSP4.data, Channel=="RCaMP" & Condition=="Nostim")

Neur.NS.Spot.DSP4<-ddply(NeuronalNostim.DSP4, c("Animal","Spot","Drug","Condition", "Animal_Spot"), summarise, nROIs=length(OnsetTime),
                    meanAmp=mean(amplitude))
Neur.NS.Spot.DSP4$Ani_Spot_Cond_Drug<-paste(Neur.NS.Spot.DSP4$Animal_Spot, Neur.NS.Spot.DSP4$Condition, Neur.NS.Spot.DSP4$Drug, sep="_")
Neur.NS.Spot.DSP4<-merge(Neur.NS.Spot.DSP4, Spot.lck.ntrials[, c("Ani_Spot_Cond_Drug", "nTrials")], by="Ani_Spot_Cond_Drug", all.x=TRUE)
Neur.NS.Spot.DSP4$ROIsPerTrial<-Neur.NS.Spot.DSP4$nROIs/Neur.NS.Spot.DSP4$nTrials


# mean
df.NeurNS.ROInum<-summarySE(Neur.NS.Spot, measurevar = "ROIsPerTrial", groupvars = c("Drug"))
df.NeurNS.amp<-summarySE(Neur.NS.Spot, measurevar = "meanAmp", groupvars = c("Drug"))

ggplot(Neur.NS.Spot, aes(x=Drug, y = ROIsPerTrial)) +
  geom_boxplot()+
  ggtitle("Spontaneous Neurons")+
  max.theme

ggplot(Neur.NS.Spot, aes(x=Drug, y = meanAmp)) +
  geom_boxplot()+
  ggtitle("Spontaneous Neurons")+ 
  max.theme

ggplot(Neur.NS.Spot.DSP4, aes(x=Drug, y = ROIsPerTrial)) +
  geom_boxplot()+
  ggtitle("Spontaneous Neurons DSP4")+ 
  max.theme

ggplot(Neur.NS.Spot.DSP4, aes(x=Drug, y = meanAmp)) +
  geom_boxplot()+
  ggtitle("Spontaneous Neurons DSP4")+ 
  max.theme

#both data sets
both.spont<-rbind(Neur.NS.Spot, Neur.NS.Spot.DSP4[Neur.NS.Spot.DSP4$Drug=="DSP4",])

ggplot(both.spont, aes(x=Drug, y = ROIsPerTrial)) +
  geom_boxplot()+
  ggtitle("Spontaneous Neurons all")+ 
  max.theme


# nROI
RC.spont.nROI.null = lmer(ROIsPerTrial ~ (1|Animal) + (1|Spot), Neur.NS.Spot,REML=FALSE)
RC.spont.nROI.model1 = lmer(ROIsPerTrial ~ Drug + (1|Animal)+ (1|Spot), Neur.NS.Spot, REML=FALSE)
RC.spont.nROI.anova <- anova(RC.spont.nROI.null, RC.spont.nROI.model1)
print(RC.spont.nROI.anova)

RC.spont.nROI<- glht(RC.spont.nROI.model1, mcp(Drug= "Tukey"))
summary(RC.spont.nROI)

# amp
RC.spont.amp.null = lmer(meanAmp ~ (1|Animal) + (1|Spot), Neur.NS.Spot,REML=FALSE)
RC.spont.amp.model1 = lmer(meanAmp ~ Drug + (1|Animal)+ (1|Spot), Neur.NS.Spot, REML=FALSE)
RC.spont.amp.anova <- anova(RC.spont.amp.null, RC.spont.amp.model1)
print(RC.spont.amp.anova)

RC.spont.amp<- glht(RC.spont.amp.model1, mcp(Drug= "Tukey"))
summary(RC.spont.amp)



# nROI .DSP4
RC.spont.DSP.nROI.null = lmer(ROIsPerTrial ~ (1|Animal) + (1|Spot), Neur.NS.Spot.DSP4,REML=FALSE)
RC.spont.DSP.nROI.model1 = lmer(ROIsPerTrial ~ Drug + (1|Animal)+ (1|Spot), Neur.NS.Spot.DSP4, REML=FALSE)
RC.spont.DSP.nROI.anova <- anova(RC.spont.DSP.nROI.null, RC.spont.DSP.nROI.model1)
print(RC.spont.DSP.nROI.anova)

RC.spont.DSP.nROI<- glht(RC.spont.DSP.nROI.model1, mcp(Drug= "Tukey"))
summary(RC.spont.DSP.nROI)

# amp .DSP4
RC.spont.DSP.amp.null = lmer(meanAmp ~ (1|Animal) + (1|Spot), Neur.NS.Spot.DSP4,REML=FALSE)
RC.spont.DSP.amp.model1 = lmer(meanAmp ~ Drug + (1|Animal)+ (1|Spot), Neur.NS.Spot.DSP4, REML=FALSE)
RC.spont.DSP.amp.anova <- anova(RC.spont.DSP.amp.null, RC.spont.DSP.amp.model1)
print(RC.spont.DSP.amp.anova)

RC.spont.DSP.amp<- glht(RC.spont.DSP.amp.model1, mcp(Drug= "Tukey"))
summary(RC.spont.DSP.amp)



############
# neuronal responses to stimulation
NeuronalStim<-subset(all.lck.peaks, Channel=="RCaMP" & Condition=="Stim")

# should have an onset time in 8 s stimulus
ggplot(NeuronalStim[NeuronalStim$OnsetTime<10,],aes(x=OnsetTime,y=..density..,fill=Drug)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP onset times between 0 and 10 s from stim trials")+
  max.theme

NeuronalStim.8s=subset(NeuronalStim, OnsetTime<8)

ggplot(NeuronalStim.8s,aes(x=Drug,y=OnsetTime,fill=Drug)) +
  geom_boxplot() +
  ggtitle("RCaMP onset times between 0 and 8 s from stim trials")+
  max.theme

# neuronal 50th percentile for each drug:
Neuron95Onset.C<-quantile(NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Control")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
#print(Neuron95Onset.C)
NeuronPT50.C1<-Neuron95Onset.C[[11]]

# atropine
Neuron95Onset.A<-quantile(NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Atropine")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
NeuronPT50.A<-Neuron95Onset.A[[11]]

# Metergoline
Neuron95Onset.M<-quantile(NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Metergoline")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
NeuronPT50.M<-Neuron95Onset.M[[11]]

# Trazodone
Neuron95Onset.T<-quantile(NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Trazodone")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
NeuronPT50.T<-Neuron95Onset.T[[11]]

# Prazosin
Neuron95Onset.P<-quantile(NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Prazosin")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
NeuronPT50.P<-Neuron95Onset.P[[11]]

print(NeuronPT50.C1)
print(NeuronPT50.A)
print(NeuronPT50.M)
print(NeuronPT50.T)
print(NeuronPT50.P)

# DSP4
NeuronalStim.8s.DSP4=subset(DSP4.data, OnsetTime<8 & Channel=="RCaMP" & Condition=="Stim")

Neuron95Onset.C2<-quantile(NeuronalStim.8s.DSP4$OnsetTime[(NeuronalStim.8s.DSP4$Drug=="Control")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
NeuronPT50.C2<-Neuron95Onset.C2[[11]]

Neuron95Onset.D<-quantile(NeuronalStim.8s.DSP4$OnsetTime[(NeuronalStim.8s.DSP4$Drug=="DSP4")], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
NeuronPT50.D<-Neuron95Onset.D[[11]]

print(NeuronPT50.C2)
print(NeuronPT50.D)


ggplot(NeuronalStim.8s.DSP4,aes(x=Drug,y=OnsetTime,fill=Drug)) +
  geom_boxplot() +
  ggtitle("RCaMP onset times between 0 and 8 s from stim trials")+
  max.theme



median.test <- function(x, y){
  z <- c(x, y)
  g <- rep(1:2, c(length(x), length(y)))
  m <- median(z)
  fisher.test(z < m, g)$p.value
}

median.test(NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Control")], NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Atropine")])

median.test(NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Control")], NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Metergoline")])

median.test(NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Control")], NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Prazosin")])

median.test(NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Control")], NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Trazodone")])


kruskal.test(list(NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Control")], 
                  NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Atropine")], 
                  NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Metergoline")],
                  NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Prazosin")],
                  NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Trazodone")]))

library('kSamples')
ad.test(NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Control")],NeuronalStim.8s$OnsetTime[(NeuronalStim.8s$Drug=="Atropine")])

#####
# Astrocytes
AstroStim<-subset(all.lck.peaks, Channel=="GCaMP" & Condition=="Stim" & OnsetTime>0)

# should have an onset time in 8 s stimulus
ggplot(AstroStim[AstroStim$OnsetTime<15,],aes(x=OnsetTime,y=..density..,fill=Drug)) +
  geom_histogram(binwidth=(0.084*5), position="dodge") +
  ggtitle("Lck-GCaMP onset times between 0 and 15 s from stim trials")+
  max.theme

ggplot(AstroStim[AstroStim$OnsetTime<15,],aes(x=Drug,y=OnsetTime,fill=Drug)) +
  geom_boxplot() +
  ggtitle("GCaMP onset times between 0 and 8 s from stim trials")+
  max.theme


AstroStim$Log_OT<-log(AstroStim$OnsetTime)

df.logOT<-summarySE(data=AstroStim, measurevar = "Log_OT", groupvars = c("Drug"),na.rm = TRUE)
df.logOT$threshold<-exp(df.logOT$Log_OT)+(2.33*exp(df.logOT$sd))


AstroNoStim<-subset(all.lck.peaks, Channel=="GCaMP" & Condition=="Nostim")

# should have an onset time in 8 s stimulus
ggplot(AstroNoStim[AstroNoStim$OnsetTime<15,],aes(x=OnsetTime,y=..density..,fill=Drug)) +
  geom_histogram(binwidth=(0.084*5), position="dodge") +
  ggtitle("Lck-GCaMP onset times between 0 and 15 s from NO stim trials")+
  max.theme

#rm(AstroNoStim, AstroStim)

####### 
# Onset time histograms- normalized to the number of trials

stimwindow=20

GC.lck.OT.dist<-subset(AstroStim, OnsetTime<stimwindow)
GC.DSP4.dist<-subset(DSP4.data, OnsetTime<stimwindow)

# number of trials
ntrials.GC.Pharmacology<- ddply(GC.lck.OT.dist, c("Drug"), summarise, ntrials=length(unique(Spot_trial_Drug)))

ntrials.GC.DSP4<- ddply(GC.DSP4.dist, c("Drug"), summarise, ntrials=length(unique(Spot_trial_Drug)))


#histogram bins
histseq= seq(0,20,0.5)
drugC=0
drugA=0
drugP=0
drugT=0
drugM=0
zeroRow<-data.frame(cbind(drugC, drugA, drugP, drugT,drugM))

# neuronal lck onset histogram
# counts for each condition in the histogram
drugC=hist(GC.lck.OT.dist$OnsetTime[GC.lck.OT.dist$Drug=="Control"], breaks=histseq, plot=FALSE)$counts
drugA=hist(GC.lck.OT.dist$OnsetTime[GC.lck.OT.dist$Drug=="Atropine"], breaks=histseq, plot=FALSE)$counts
drugP=hist(GC.lck.OT.dist$OnsetTime[GC.lck.OT.dist$Drug=="Prazosin"], breaks=histseq, plot=FALSE)$counts
drugT=hist(GC.lck.OT.dist$OnsetTime[GC.lck.OT.dist$Drug=="Trazodone"], breaks=histseq, plot=FALSE)$counts
drugM=hist(GC.lck.OT.dist$OnsetTime[GC.lck.OT.dist$Drug=="Metergoline"], breaks=histseq, plot=FALSE)$counts

# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
drugC=drugC/ntrials.GC.Pharmacology$ntrials[(ntrials.GC.Pharmacology$Drug=="Control")]
drugA=drugA/ntrials.GC.Pharmacology$ntrials[(ntrials.GC.Pharmacology$Drug=="Atropine")]
drugP=drugP/ntrials.GC.Pharmacology$ntrials[(ntrials.GC.Pharmacology$Drug=="Prazosin")]
drugT=drugT/ntrials.GC.Pharmacology$ntrials[(ntrials.GC.Pharmacology$Drug=="Trazodone")]
drugM=drugM/ntrials.GC.Pharmacology$ntrials[(ntrials.GC.Pharmacology$Drug=="Metergoline")]
drugD=drugD/ntrials.GC.Pharmacology$ntrials[(ntrials.GC.Pharmacology$Drug=="DSP4")]

#make a data frame for plotting
lck.histo <- data.frame(cbind(drugC, drugA, drugP, drugT,drugM,drugD))
lck.histo2<-rbind(zeroRow,lck.histo)
lck.histo2$time<-histseq

ggplot(NULL, aes(x=time))+
  geom_line(data=lck.histo2, aes(y=drugC, color="Control")) +
  geom_line(data=lck.histo2, aes(y=drugA, color="Atropine")) +
  geom_line(data=lck.histo2, aes(y=drugP, color="Prazosin")) +
  geom_line(data=lck.histo2, aes(y=drugT, color="Trazodone")) +
  geom_line(data=lck.histo2, aes(y=drugM, color="Metergoline")) +
  ggtitle("Pharmacology astrocytes-lck data") + 
  xlab("Onset Time (s)") + 
  ylab("Astrocytes peaks/trial") + 
  xlim(-0.5,20) +
  max.theme



# DSP4 distributions
#histogram bins
histseq= seq(0,20,0.5)
drugC=0
drugD=0
zeroRow<-data.frame(cbind(drugC, drugD))

# neuronal lck onset histogram
# counts for each condition in the histogram
drugC=hist(GC.DSP4.dist$OnsetTime[GC.DSP4.dist$Drug=="Control"], breaks=histseq, plot=FALSE)$counts
drugD=hist(GC.DSP4.dist$OnsetTime[GC.DSP4.dist$Drug=="DSP4"], breaks=histseq, plot=FALSE)$counts


# normalized: divide each bin (number of ROIs) by the total number of trials for this condition
drugC=drugC/ntrials.GC.DSP4$ntrials[(ntrials.GC.DSP4$Drug=="Control")]
drugD=drugD/ntrials.GC.DSP4$ntrials[(ntrials.GC.DSP4$Drug=="DSP4")]

#make a data frame for plotting
lck.histo <- data.frame(cbind(drugC, drugD))
lck.histo2<-rbind(zeroRow,lck.histo)
lck.histo2$time<-histseq

ggplot(NULL, aes(x=time))+
  geom_line(data=lck.histo2, aes(y=drugC, color="Control")) +
  geom_line(data=lck.histo2, aes(y=drugD, color="DSP4")) +
  ggtitle("Pharmacology astrocytes- DSP data") + 
  xlab("Onset Time (s)") + 
  ylab("Astrocytes peaks/trial") + 
  xlim(-0.5,20) +
  max.theme


#compare distributions (what is plotted in the figure)

# ks test- astrocytes
OT.lck.GC.pharm.kstest<- ks.test(drugC,drugA, drugP, drugT,drugM,drugD)
print(OT.lck.GC.pharma.kstest)

#library("SuppDist")
library("kSamples")

OT.lck.GC.pharm.adtest<- ad.test(drugC,drugA, drugP, drugT,drugM,drugD)
print(OT.lck.GC.pharm.adtest)


######
# for median and mean calculations

# no stim vs 8 s stim- 
#neuronal window=9s, AC window= 15 s for peak time, 
#neuronal window=2s, AC window=12 s for onset


LongN_OTwind=8
LongAC_OTwind=12.76
fastTh<-1.09

# remove data that is outside the above windows
# lck
stim.lck.OT.R<-subset(all.lck.peaks, Channel=="RCaMP" & OnsetTime<=LongN_OTwind & peakTime>=0 & Duration<45)
stim.lck.OT.G<-subset(all.lck.peaks, Channel=="GCaMP" & OnsetTime<=LongAC_OTwind & peakTime>=0 & Duration<45)

lck.peaks.window<-rbind(stim.lck.OT.R, stim.lck.OT.G)

rm(stim.lck.OT.G,stim.lck.OT.R)

# DSP4
stim.lck.OT.R<-subset(DSP4.data, Channel=="RCaMP" & OnsetTime<=LongN_OTwind & peakTime>=0 & Duration<45)
stim.lck.OT.G<-subset(DSP4.data, Channel=="GCaMP" & OnsetTime<=LongAC_OTwind & peakTime>=0 & Duration<45)

DSP4.data.window<-rbind(stim.lck.OT.R, stim.lck.OT.G)

rm(stim.lck.OT.G,stim.lck.OT.R)


# onset times
ggplot(lck.peaks.window[(lck.peaks.window$Channel=="RCaMP" & lck.peaks.window$Condition=="Stim"),], aes(x=interaction(Channel, Drug),y=OnsetTime, fill= Condition)) +
  geom_boxplot()+
  ylab("Onset Time (s)") +
  ggtitle("time window 0-12.76s, rcamp, stim")+
  max.theme

ggplot(lck.peaks.window[(lck.peaks.window$Channel=="GCaMP" & lck.peaks.window$Condition=="Stim"),], aes(x=interaction(Channel, Drug),y=OnsetTime, fill= Condition)) +
  geom_boxplot()+
  ylab("Onset Time (s)") +
  ggtitle("time window 0-12.76s, gcamp, stim")+
  max.theme


#peak times
ggplot(lck.peaks.window[(lck.peaks.window$Channel=="RCaMP" & lck.peaks.window$Condition=="Stim"),], aes(x=interaction(Channel, Drug),y=peakTime, fill= Condition)) +
  geom_boxplot()+
  ylab("peak Time (s)") +
  ggtitle("time window 0-12.76s, rcamp, stim")+
  max.theme

ggplot(lck.peaks.window[(lck.peaks.window$Channel=="RCaMP" & lck.peaks.window$Condition=="Stim"),], aes(x=interaction(Channel, Drug),y=peakTime, fill= Condition)) +
  geom_boxplot()+
  ylab("peak Time (s)") +
  ggtitle("time window 0-12.76s, rcamp, stim")+
  max.theme

ggplot(lck.peaks.window[(lck.peaks.window$Channel=="GCaMP" & lck.peaks.window$Condition=="Stim"),], aes(x=interaction(Channel, Drug),y=peakTime, fill= Condition)) +
  geom_boxplot()+
  ylab("peak Time (s)") +
  ggtitle("time window 0-12.76s, gcamp, stim")+
  max.theme


# DSP4 data
# onset times
ggplot(DSP4.data.window, aes(x=interaction(Channel, Drug),y=OnsetTime, fill= Condition)) +
  geom_boxplot(notch=TRUE)+
  ylab("OnsetTimes Time (s)") +
  ggtitle("matched data= peak times and onset times")+
  max.theme

#peak times
ggplot(DSP4.data.window, aes(x=interaction(Channel, Drug),y=peakTime, fill= Condition)) +
  geom_boxplot(notch=TRUE)+
  ylab("peak Time (s)") +
  ggtitle("matched data= peak times and onset times")+
  max.theme



######

# identify "FAST" astrocytes
lck.peaks.window$Group<-0
lck.peaks.window$Group[lck.peaks.window$OnsetTime<fastTh]<-"fast_MDs"
lck.peaks.window$Group[lck.peaks.window$OnsetTime>=fastTh]<-"delayed_MDs"

lck.peaks.window$Group <- factor(lck.peaks.window$Group, levels = c("fast_MDs","delayed_MDs"))
lck.peaks.window$Channel <- factor(lck.peaks.window$Channel, levels = c("RCaMP","GCaMP"))


# fast vs delayed
lck.peaks.window$Channel_Group<-interaction(lck.peaks.window$Channel, lck.peaks.window$Group)
lck.peaks.window$Channel_Group<-as.factor(lck.peaks.window$Channel_Group)

stim.lck.compdata<-lck.peaks.window[!(lck.peaks.window$Channel=="RCaMP"& lck.peaks.window$Group=="delayed_MDs"),]

# take out the effect of Condition
# we are only interested in stim case
stim.lck.compdata.STIM<-subset(stim.lck.compdata, Condition=="Stim")


# DSP4 data
# identify "FAST" astrocytes
DSP4.data.window$Group<-0
DSP4.data.window$Group[DSP4.data.window$OnsetTime<fastTh]<-"fast_MDs"
DSP4.data.window$Group[DSP4.data.window$OnsetTime>=fastTh]<-"delayed_MDs"

DSP4.data.window$Group <- factor(DSP4.data.window$Group, levels = c("fast_MDs","delayed_MDs"))
DSP4.data.window$Channel <- factor(DSP4.data.window$Channel, levels = c("RCaMP","GCaMP"))


# fast vs delayed
DSP4.data.window$Channel_Group<-interaction(DSP4.data.window$Channel, DSP4.data.window$Group)
DSP4.data.window$Channel_Group<-as.factor(DSP4.data.window$Channel_Group)

DSP4.data.compdata<-DSP4.data.window[!(DSP4.data.window$Channel=="RCaMP"& DSP4.data.window$Group=="delayed_MDs"),]

# take out the effect of Condition
# we are only interested in stim case
DSP4.data.compdata.STIM<-subset(DSP4.data.compdata, Condition=="Stim")



######
                      
                      


# number of ROIs in each trial for each field of view (across the whole trial)
lck.8strial<-lck.peaks.window#[lck.peaks.window$OnsetTime<8,]
DSP4.8strial<-DSP4.data.window#[DSP4.data.window$OnsetTime<8,]

lck.8strial$Channel <- factor(lck.8strial$Channel, levels = c("RCaMP","GCaMP"))
DSP4.8strial$Channel <- factor(DSP4.8strial$Channel, levels = c("RCaMP","GCaMP"))

ROInum.8strial<-ddply(lck.8strial, c("Animal","Spot","Drug","Condition","Channel","Animal_Spot"), summarise, nROIs=length(OnsetTime))
ROInum.8strial.group<-ddply(lck.8strial, c("Animal","Spot","Drug","Condition","Channel","Animal_Spot","Group"), summarise, nROIs=length(OnsetTime))

DSP4.ROInum.8strial<-ddply(DSP4.8strial, c("Animal","Spot","Drug","Condition","Channel","Animal_Spot"), summarise, nROIs=length(OnsetTime))
DSP4.ROInum.8strial.group<-ddply(DSP4.8strial, c("Animal","Spot","Drug","Condition","Channel","Animal_Spot","Group"), summarise, nROIs=length(OnsetTime))

# add in number of trials
ROInum.8strial$Ani_Spot_Cond_Drug<-paste(ROInum.8strial$Animal_Spot, ROInum.8strial$Condition,ROInum.8strial$Drug, sep="_")
ROInum.8strial.group$Ani_Spot_Cond_Drug<-paste(ROInum.8strial.group$Animal_Spot, ROInum.8strial.group$Condition,ROInum.8strial.group$Drug, sep="_")
DSP4.ROInum.8strial.group$Ani_Spot_Cond_Drug<-paste(DSP4.ROInum.8strial.group$Animal_Spot, DSP4.ROInum.8strial.group$Condition,DSP4.ROInum.8strial.group$Drug, sep="_")
DSP4.ROInum.8strial$Ani_Spot_Cond_Drug<-paste(DSP4.ROInum.8strial$Animal_Spot, DSP4.ROInum.8strial$Condition,DSP4.ROInum.8strial$Drug, sep="_")

Spot.lck.ntrials$Ani_Spot_Cond_Drug<-paste(Spot.lck.ntrials$Animal, Spot.lck.ntrials$Spot,Spot.lck.ntrials$Condition,Spot.lck.ntrials$Drug, sep="_")
ROInum.8strial<-merge(ROInum.8strial, Spot.lck.ntrials[, c("Ani_Spot_Cond_Drug", "nTrials")], by="Ani_Spot_Cond_Drug", all.x=TRUE)
ROInum.8strial$ROIsPerTrial<-ROInum.8strial$nROIs/ROInum.8strial$nTrials
ROInum.8strial.group<-merge(ROInum.8strial.group, Spot.lck.ntrials[, c("Ani_Spot_Cond_Drug", "nTrials")], by="Ani_Spot_Cond_Drug", all.x=TRUE)
ROInum.8strial.group$ROIsPerTrial<-ROInum.8strial.group$nROIs/ROInum.8strial.group$nTrials

DSP4.ROInum.8strial<-merge(DSP4.ROInum.8strial, Spot.lck.ntrials[, c("Ani_Spot_Cond_Drug", "nTrials")], by="Ani_Spot_Cond_Drug", all.x=TRUE)
DSP4.ROInum.8strial$ROIsPerTrial<-DSP4.ROInum.8strial$nROIs/DSP4.ROInum.8strial$nTrials
DSP4.ROInum.8strial.group<-merge(DSP4.ROInum.8strial.group, Spot.lck.ntrials[, c("Ani_Spot_Cond_Drug", "nTrials")], by="Ani_Spot_Cond_Drug", all.x=TRUE)
DSP4.ROInum.8strial.group$ROIsPerTrial<-DSP4.ROInum.8strial.group$nROIs/DSP4.ROInum.8strial.group$nTrials

# mean
df.ROInum.8strial<-summarySE(ROInum.8strial, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition","Drug"))
df.ROInum.8strial.group<-summarySE(ROInum.8strial.group, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition","Drug","Group"))

df.ROInum.8strial<-summarySE(ROInum.8strial, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition","Drug"))
df.ROInum.8strial.group<-summarySE(ROInum.8strial.group, measurevar = "ROIsPerTrial", groupvars = c("Channel", "Condition","Drug","Group"))


# plots
#boxplot
ggplot(ROInum.8strial, aes(x=interaction(Drug, Channel),y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("Lck- ROIs per FOV") +
  max.theme

ggplot(data=df.ROInum.8strial, aes(x=interaction(Drug, Channel), y= ROIsPerTrial, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIsPerTrial") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(df.ROInum.8strial[df.ROInum.8strial$Channel=="RCaMP",], aes(x=Drug,y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  max.theme

ggplot(df.ROInum.8strial[df.ROInum.8strial$Channel=="RCaMP",],aes(x=Drug,y=ROIsPerTrial, group=Condition)) +
  geom_point(aes(x=Drug,y=ROIsPerTrial),stat="identity", size=3)+
  geom_line(aes(x=Drug, y=ROIsPerTrial, colour=Condition))+
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("RCaMP ROIs")+
  max.theme

ggplot(ROInum.8strial[ROInum.8strial$Channel=="RCaMP",], aes(x=interaction(Drug, Channel),y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("Lck- ROIs per FOV") +
  max.theme

# gcamp
ggplot(df.ROInum.8strial[df.ROInum.8strial$Channel=="GCaMP",], aes(x=Drug,y=ROIsPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  max.theme

ggplot(df.ROInum.8strial[df.ROInum.8strial$Channel=="GCaMP",],aes(x=Drug,y=ROIsPerTrial, group=Condition)) +
  geom_point(aes(x=Drug,y=ROIsPerTrial),stat="identity", size=3)+
  geom_line(aes(x=Drug, y=ROIsPerTrial, colour=Condition))+
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("GCaMP ROIs") +
  max.theme

ggplot(ROInum.8strial[ROInum.8strial$Channel=="GCaMP",], aes(x=interaction(Drug, Channel),y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("Lck- ROIs per FOV") +
  max.theme

# only fast GcaMP
ggplot(ROInum.8strial.group[(ROInum.8strial.group$Channel=="GCaMP" & ROInum.8strial.group$Condition=="Stim" & ROInum.8strial.group$Group=="fast_MDs"),], aes(x=interaction(Drug, Channel),y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("Lck- fast ROIs per FOV") +
  max.theme

#only delayed GCaMP
ggplot(ROInum.8strial.group[(ROInum.8strial.group$Channel=="GCaMP" & ROInum.8strial.group$Condition=="Stim" & ROInum.8strial.group$Group=="delayed_MDs"),], aes(x=interaction(Drug, Channel),y=ROIsPerTrial, fill= Condition)) +
  geom_boxplot()+
  ylab("Lck- delayed ROIs per FOV") +
  max.theme

# RCaMP
Condition_Drug_RC= interaction(ROInum.8strial$Condition[ROInum.8strial$Channel=="RCaMP"],ROInum.8strial$Drug[ROInum.8strial$Channel=="RCaMP"])
nROI.RC.stim.null = lmer(ROIsPerTrial ~ (1|Animal) + (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="RCaMP",],REML=FALSE)
nROI.RC.stim.model1 = lmer(ROIsPerTrial ~ Condition + (1|Animal)+ (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="RCaMP",],REML=FALSE)
nROI.RC.stim.model2 = lmer(ROIsPerTrial ~ Condition_Drug_RC + (1|Animal)+ (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="RCaMP",],REML=FALSE)
nROI.RC.stim.anova <- anova(nROI.RC.stim.null, nROI.RC.stim.model1,nROI.RC.stim.model2)
print(nROI.RC.stim.anova)

nROI.RC.stim.Cond_Drug<- glht(nROI.RC.stim.model2, mcp(Condition_Drug_RC= "Tukey"))
summary(nROI.RC.stim.Cond_Drug)

# ALL DRUGS ARE P<0.01 less than control for stim case
# atropine and metergoline do not have sig difference between no stim and stim

# GCaMP
Condition_Drug_GC= interaction(ROInum.8strial$Condition[ROInum.8strial$Channel=="GCaMP"],ROInum.8strial$Drug[ROInum.8strial$Channel=="GCaMP"])
nROI.GC.stim.null = lmer(ROIsPerTrial ~ (1|Animal) + (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="GCaMP",],REML=FALSE)
nROI.GC.stim.model1 = lmer(ROIsPerTrial ~ Condition + (1|Animal)+ (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="GCaMP",],REML=FALSE)
nROI.GC.stim.model2 = lmer(ROIsPerTrial ~ Condition_Drug_GC + (1|Animal)+ (1|Spot), ROInum.8strial[ROInum.8strial$Channel=="GCaMP",],REML=FALSE)
nROI.GC.stim.anova <- anova(nROI.GC.stim.null, nROI.GC.stim.model1,nROI.GC.stim.model2)
print(nROI.GC.stim.anova)

nROI.GC.stim.Cond_Drug<- glht(nROI.GC.stim.model2, mcp(Condition_Drug_GC= "Tukey"))
summary(nROI.GC.stim.Cond_Drug)

# metergoline, trazodone, prazosin- no sig difference between stim and control

# trazodone, and prazosin- sig fewer ROIs compared to control


# fast GC
Condition_Drug_GC_fast= interaction(ROInum.8strial.group$Condition[(ROInum.8strial.group$Channel=="GCaMP" & ROInum.8strial.group$Group=="fast_MDs")],
                               ROInum.8strial.group$Drug[(ROInum.8strial.group$Channel=="GCaMP"& ROInum.8strial.group$Group=="fast_MDs")])
nROI.GC.fast.stim.null = lmer(ROIsPerTrial ~ (1|Animal) + (1|Spot), 
                              ROInum.8strial.group[(ROInum.8strial.group$Channel=="GCaMP" & ROInum.8strial.group$Group=="fast_MDs"),],
                              REML=FALSE)
nROI.GC.fast.stim.model1 = lmer(ROIsPerTrial ~ Condition + (1|Animal) + (1|Spot), 
                              ROInum.8strial.group[(ROInum.8strial.group$Channel=="GCaMP" & ROInum.8strial.group$Group=="fast_MDs"),],
                              REML=FALSE)
nROI.GC.fast.stim.model2 = lmer(ROIsPerTrial ~ Drug + (1|Animal) + (1|Spot), 
                              ROInum.8strial.group[(ROInum.8strial.group$Channel=="GCaMP" & ROInum.8strial.group$Group=="fast_MDs"),],
                              REML=FALSE)
nROI.GC.fast.stim.model3 = lmer(ROIsPerTrial ~ Condition_Drug_GC_fast + (1|Animal) + (1|Spot), 
                                ROInum.8strial.group[(ROInum.8strial.group$Channel=="GCaMP" & ROInum.8strial.group$Group=="fast_MDs"),],
                                REML=FALSE)
nROI.GC.fast.stim.anova <- anova(nROI.GC.fast.stim.null, nROI.GC.fast.stim.model1,
                                 nROI.GC.fast.stim.model2, nROI.GC.fast.stim.model3)
print(nROI.GC.fast.stim.anova)

nROI.GC.fast.stim.Cond_Drug<- glht(nROI.GC.fast.stim.model3, mcp(Condition_Drug_GC_fast= "Tukey"))
summary(nROI.GC.fast.stim.Cond_Drug)

# no difference in ROI number across fast AC types


# delayed GC
Condition_Drug_GC_delayed= interaction(ROInum.8strial.group$Condition[(ROInum.8strial.group$Channel=="GCaMP" & ROInum.8strial.group$Group=="delayed_MDs")],
                                    ROInum.8strial.group$Drug[(ROInum.8strial.group$Channel=="GCaMP"& ROInum.8strial.group$Group=="delayed_MDs")])
nROI.GC.delayed.stim.null = lmer(ROIsPerTrial ~ (1|Animal) + (1|Spot), 
                              ROInum.8strial.group[(ROInum.8strial.group$Channel=="GCaMP" & ROInum.8strial.group$Group=="delayed_MDs"),],
                              REML=FALSE)
nROI.GC.delayed.stim.model1 = lmer(ROIsPerTrial ~ Condition + (1|Animal) + (1|Spot), 
                                ROInum.8strial.group[(ROInum.8strial.group$Channel=="GCaMP" & ROInum.8strial.group$Group=="delayed_MDs"),],
                                REML=FALSE)
nROI.GC.delayed.stim.model2 = lmer(ROIsPerTrial ~ Drug + (1|Animal) + (1|Spot), 
                                ROInum.8strial.group[(ROInum.8strial.group$Channel=="GCaMP" & ROInum.8strial.group$Group=="delayed_MDs"),],
                                REML=FALSE)
nROI.GC.delayed.stim.model3 = lmer(ROIsPerTrial ~ Condition_Drug_GC_delayed + (1|Animal) + (1|Spot), 
                                ROInum.8strial.group[(ROInum.8strial.group$Channel=="GCaMP" & ROInum.8strial.group$Group=="delayed_MDs"),],
                                REML=FALSE)
nROI.GC.delayed.stim.anova <- anova(nROI.GC.delayed.stim.null, nROI.GC.delayed.stim.model1,
                                 nROI.GC.delayed.stim.model2, nROI.GC.delayed.stim.model3)
print(nROI.GC.delayed.stim.anova)

nROI.GC.delayed.stim.Cond_Drug<- glht(nROI.GC.delayed.stim.model3, mcp(Condition_Drug_GC_delayed= "Tukey"))
summary(nROI.GC.delayed.stim.Cond_Drug)







######
# fraction of active pixels from all the astrocyte pixels
# active pixels from any ROI with a signal (not just those during stimulation)

# number of ROIs in each trial for each field of view (across the time window (2 s for neurons, 15 s for AC))
lck.peaks.window$Channel <- factor(lck.peaks.window$Channel, levels = c("RCaMP","GCaMP"))

Spot.lck.stim<-ddply(lck.peaks.window, c("Animal","Spot","Drug","Trial","Condition","Channel","nFluoPix","nActivePix","pixelsize"), summarise,
                     nROIs=length(unique(ROIs_Cond_Drug)))

# add in number of trials
Spot.lck.stim$Ani_Spot_Cond<-paste(Spot.lck.stim$Animal, Spot.lck.stim$Spot, Spot.lck.stim$Condition,Spot.lck.stim$Drug, sep="_")
Spot.lck.ntrials$Ani_Spot_Cond<-paste(Spot.lck.ntrials$Animal, Spot.lck.ntrials$Spot, Spot.lck.ntrials$Condition,Spot.lck.ntrials$Drug, sep="_")

Spot.lck.stim<-merge(Spot.lck.stim, Spot.lck.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)

Spot.lck.stim$ROIsPerTrial<-Spot.lck.stim$nROIs/Spot.lck.stim$nTrials

Spot.lck.stim$FracActive=Spot.lck.stim$nActivePix/Spot.lck.stim$nFluoPix
Spot.lck.stim$FracActivePerROI=Spot.lck.stim$FracActive/Spot.lck.stim$nROIs

# mean fraction of active pixels
df.FracActive<- summarySE(Spot.lck.stim[Spot.lck.stim$Channel=="GCaMP",], measurevar = "FracActive", groupvars = c("Condition", "Drug"))
df.FracActivePerROI<- summarySE(Spot.lck.stim[Spot.lck.stim$Channel=="GCaMP",], measurevar = "FracActivePerROI", groupvars = c("Condition", "Drug"))

ggplot(data=df.FracActive, aes(x=Drug, y= FracActive, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=FracActive-se, ymax=FracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.FracActivePerROI, aes(x=Drug, y= FracActivePerROI, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=FracActivePerROI-se, ymax=FracActivePerROI+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive Per ROI") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

#only astrocytes
Cond_Drug= interaction(Spot.lck.stim$Condition[Spot.lck.stim$Channel=="GCaMP"],Spot.lck.stim$Drug[Spot.lck.stim$Channel=="GCaMP"])
FracActive.lck.stim.null = lmer(FracActive ~ (1|Animal)+ (1|Spot), Spot.lck.stim[Spot.lck.stim$Channel=="GCaMP",],REML=FALSE)
FracActive.lck.stim.model1 = lmer(FracActive~ Condition+ (1|Animal)+ (1|Spot), Spot.lck.stim[Spot.lck.stim$Channel=="GCaMP",],REML=FALSE)
FracActive.lck.stim.model2 = lmer(FracActive ~ Drug + (1|Animal)+ (1|Spot), Spot.lck.stim[Spot.lck.stim$Channel=="GCaMP",],REML=FALSE)
FracActive.lck.stim.model3 = lmer(FracActive ~ Condition + Drug + (1|Animal)+ (1|Spot), Spot.lck.stim[Spot.lck.stim$Channel=="GCaMP",],REML=FALSE)
FracActive.lck.stim.model4 = lmer(FracActive ~ Cond_Drug + (1|Animal)+ (1|Spot), Spot.lck.stim[Spot.lck.stim$Channel=="GCaMP",],REML=FALSE)
FracActive.lck.stim.anova <- anova(FracActive.lck.stim.null, FracActive.lck.stim.model1,FracActive.lck.stim.model2,
                                 FracActive.lck.stim.model3,FracActive.lck.stim.model4)
print(FracActive.lck.stim.anova)

FracActive.lck.stim.Cond_Drug<- glht(FracActive.lck.stim.model4, mcp(Cond_Drug= "Tukey"))
summary(FracActive.lck.stim.Cond_Drug)



######

# calculate the change in amplitude, onset time, AUC, etc. with drug since it is the same cells in each FOV

#########
# mean onset times
df.OT1<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim",], measurevar = "OnsetTime", groupvars = c("Channel_Group", "Drug"))
df.OT2<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",], measurevar = "OnsetTime", groupvars = c("Group", "Drug"))
df.OT3<- summarySE(stim.lck.compdata.STIM, measurevar = "OnsetTime", groupvars = c("ROIType","Channel_Group", "Drug"))
df.OT4<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="RCaMP",], measurevar = "OnsetTime", groupvars = c("Group", "Drug"))

ggplot(df.OT1, aes(x=Channel_Group,y=OnsetTime, fill=Drug)) +
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

ggplot(df.OT4, aes(x=Drug,y=OnsetTime, fill=Drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  ggtitl("RCaMP")
  max.theme

# stats
Group_Channel_Type_Cond_Drug=interaction(lck.peaks.window$Group,lck.peaks.window$Channel,lck.peaks.window$ROIType, 
                                        lck.peaks.window$Drug,lck.peaks.window$Condition)
Group_Channel_Cond_Drug=interaction(lck.peaks.window$Group,lck.peaks.window$Channel,lck.peaks.window$Drug,
                                   lck.peaks.window$Condition)

# stats for onset times- neurons vs astrocytes
OT.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Trial), lck.peaks.window,REML=FALSE)
OT.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Trial), lck.peaks.window,REML=FALSE)
OT.model2 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|Trial), lck.peaks.window,REML=FALSE)
OT.model3 = lmer(OnsetTime ~ Group_Channel_Cond_Drug + (1|Animal) + (1|Spot) + (1|Trial), lck.peaks.window,REML=FALSE)
OT.model4 = lmer(OnsetTime ~ Group_Channel_Type_Cond_Drug + (1|Animal) + (1|Spot) + (1|Trial), lck.peaks.window,REML=FALSE)
OT.anova <- anova(OT.null, OT.model1,OT.model2,OT.model3,OT.model4)
print(OT.anova)

#OT.Group_channel_Drug<- glht(OT.model3, mcp(Group_Channel_Cond_Drug= "Tukey"))
#summary(OT.Group_channel_Drug)

#OT.Group_Channel_Type_Drug<- glht(OT.model4, mcp(Group_Channel_Type_Cond_Drug= "Tukey"))
#summary(OT.Group_Channel_Type_Drug)


Group_Channel_Type_Drug=interaction(stim.lck.compdata.STIM$Group,stim.lck.compdata.STIM$Channel,
                                   stim.lck.compdata.STIM$ROIType, stim.lck.compdata.STIM$Drug)
Group_Channel_Drug=interaction(stim.lck.compdata.STIM$Channel, stim.lck.compdata.STIM$Group, stim.lck.compdata.STIM$Drug)

# stats for onset times- neurons vs astrocytes
OT.stim.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Trial), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Trial), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.model3 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|Trial), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.model4 = lmer(OnsetTime ~ Group_Channel_Drug + (1|Animal) + (1|Spot) + (1|Trial), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.model6 = lmer(OnsetTime ~ Group_Channel_Type_Drug + (1|Animal) + (1|Spot) + (1|Trial), stim.lck.compdata.STIM,REML=FALSE)
OT.stim.anova <- anova(OT.stim.null, OT.stim.model1,OT.stim.model3,OT.stim.model4,OT.stim.model6)
print(OT.stim.anova)

OT.stim.Group_channel_Drug<- glht(OT.stim.model4, mcp(Group_Channel_Drug= "Tukey"))
summary(OT.stim.Group_channel_Drug)

#OT.stim.Group_channel_type_Drug<- glht(OT.stim.model6, mcp(Group_Channel_Type_Drug= "Tukey"))
#summary(OT.stim.Group_channel_type_Drug)

summary(OT.stim.model4)

# check residuals for linearity
plot(fitted(OT.stim.model4), residuals(OT.stim.model4),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(OT.stim.model4), residuals(OT.stim.model4)), col=46, lwd=2.5)


# GCaMP fast vs delayed only
# stats for onset times- neurons vs astrocytes
Group_Drug=interaction(stim.lck.compdata.STIM$Group[stim.lck.compdata.STIM$Channel=="GCaMP"],
                       stim.lck.compdata.STIM$Drug[stim.lck.compdata.STIM$Channel=="GCaMP"])
OT.stim2.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Trial), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
OT.stim2.model1 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|Trial), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
OT.stim2.model2 = lmer(OnsetTime ~ Group_Drug + (1|Animal) + (1|Spot) + (1|Trial), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
OT.stim2.anova <- anova(OT.stim2.null, OT.stim2.model1, OT.stim2.model2)
print(OT.stim2.anova)

OT.stim.Group_Drug<- glht(OT.stim2.model2, mcp(Group_Drug= "Tukey"))
summary(OT.stim.Group_Drug)


#fast Neurons only
# stats for onset times- neurons vs astrocytes
OT.stim3.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Trial), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="RCaMP",],REML=FALSE)
OT.stim3.model2 = lmer(OnsetTime ~ Drug + (1|Animal) + (1|Spot) + (1|Trial), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="RCaMP",],REML=FALSE)
OT.stim3.anova <- anova(OT.stim3.null, OT.stim3.model2)
print(OT.stim3.anova)

OT.stim.NeuronDrug<- glht(OT.stim3.model2, mcp(Drug= "Tukey"))
summary(OT.stim.NeuronDrug)

#############
# mean peak time
df.PT1<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim",], measurevar = "peakTime", groupvars = c("Channel_Group", "Drug"))
df.PT2<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",], measurevar = "peakTime", groupvars = c("Group", "Drug"))
df.PT3<- summarySE(stim.lck.compdata.STIM, measurevar = "peakTime", groupvars = c("ROIType","Channel_Group", "Drug"))

ggplot(df.PT1, aes(x=Channel_Group,y=peakTime, fill=Drug)) +
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
PT.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|Trial), lck.peaks.window,REML=FALSE)
PT.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|Trial), lck.peaks.window,REML=FALSE)
PT.model2 = lmer(peakTime ~ Group + (1|Animal) + (1|Spot) + (1|Trial), lck.peaks.window,REML=FALSE)
PT.model3 = lmer(peakTime ~ Group_Channel_Cond_Drug + (1|Animal) + (1|Spot) + (1|Trial), lck.peaks.window,REML=FALSE)
PT.model4 = lmer(peakTime ~ Group_Channel_Type_Cond_Drug + (1|Animal) + (1|Spot) + (1|Trial), lck.peaks.window,REML=FALSE)
PT.anova <- anova(PT.null, PT.model1,PT.model2,PT.model3,PT.model4)
print(PT.anova)

#PT.Group_channel_Drug<- glht(PT.model3, mcp(Group_Channel_Cond_Drug= "Tukey"))
#summary(PT.Group_channel_Drug)

#PT.Group_Channel_Type_Drug<- glht(PT.model4, mcp(Group_Channel_Type_Cond_Drug= "Tukey"))
#summary(PT.Group_Channel_Type_Drug)


#Group_Channel_Type_Drug=interaction(stim.lck.compdata.STIM$Group,stim.lck.compdata.STIM$Channel,
 #                                  stim.lck.compdata.STIM$ROIType, stim.lck.compdata.STIM$Drug)
#Group_Channel_Drug=interaction(stim.lck.compdata.STIM$Channel, stim.lck.compdata.STIM$Group, stim.lck.compdata.STIM$Drug)

# stats for onset times- neurons vs astrocytes
PT.stim.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|Trial), stim.lck.compdata.STIM,REML=FALSE)
PT.stim.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|Trial), stim.lck.compdata.STIM,REML=FALSE)
PT.stim.model3 = lmer(peakTime ~ Group + (1|Animal) + (1|Spot) + (1|Trial), stim.lck.compdata.STIM,REML=FALSE)
PT.stim.model4 = lmer(peakTime ~ Group_Channel_Drug + (1|Animal) + (1|Spot) + (1|Trial), stim.lck.compdata.STIM,REML=FALSE)
PT.stim.model6 = lmer(peakTime ~ Group_Channel_Type_Drug + (1|Animal) + (1|Spot) + (1|Trial), stim.lck.compdata.STIM,REML=FALSE)
PT.stim.anova <- anova(PT.stim.null, PT.stim.model1,PT.stim.model3,PT.stim.model4,PT.stim.model6)
print(PT.stim.anova)

PT.stim.Group_channel_Drug<- glht(PT.stim.model4, mcp(Group_Channel_Drug= "Tukey"))
summary(PT.stim.Group_channel_Drug)

#PT.stim.Group_channel_type_Drug<- glht(PT.stim.model6, mcp(Group_Channel_Type_Drug= "Tukey"))
#summary(PT.stim.Group_channel_type_Drug)


########
#amplitude
df.amp1<-summarySE(lck.peaks.window, measurevar = "amplitude", groupvars = c("Channel","Condition", "Drug"))
df.amp2<-summarySE(lck.peaks.window, measurevar = "amplitude", groupvars = c("Channel", "ROIType","Condition", "Drug"))

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
df.Dur1<-summarySE(lck.peaks.window, measurevar = "Duration", groupvars = c("Channel","Condition", "Drug"))
df.Dur2<-summarySE(lck.peaks.window, measurevar = "Duration", groupvars = c("Channel", "ROIType","Condition", "Drug"))

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
# proportion of ROIs that are fast or delayed

#what about proportion of each group based on ROI type??

GCaMP.C<-subset(lck.peaks.window,Channel=="GCaMP" & Drug=="Control" & Condition=="Stim")
GCaMP.A<-subset(lck.peaks.window,Channel=="GCaMP" & Drug=="Atropine" & Condition=="Stim")
GCaMP.T<-subset(lck.peaks.window,Channel=="GCaMP" & Drug=="Trazodone" & Condition=="Stim")
GCaMP.P<-subset(lck.peaks.window,Channel=="GCaMP" & Drug=="Prazosin" & Condition=="Stim")
GCaMP.M<-subset(lck.peaks.window,Channel=="GCaMP" & Drug=="Metergoline" & Condition=="Stim")
GCaMP.D<-subset(lck.peaks.window,Channel=="GCaMP" & Drug=="DSP4" & Condition=="Stim")


RCaMP.C<-subset(lck.peaks.window,Channel=="RCaMP" & Drug=="Control" & Condition=="Stim")
RCaMP.A<-subset(lck.peaks.window,Channel=="RCaMP" & Drug=="Atropine" & Condition=="Stim")
RCaMP.T<-subset(lck.peaks.window,Channel=="RCaMP" & Drug=="Trazodone" & Condition=="Stim")
RCaMP.P<-subset(lck.peaks.window,Channel=="RCaMP" & Drug=="Prazosin" & Condition=="Stim")
RCaMP.M<-subset(lck.peaks.window,Channel=="RCaMP" & Drug=="Metergoline" & Condition=="Stim")
RCaMP.D<-subset(lck.peaks.window,Channel=="RCaMP" & Drug=="DSP4" & Condition=="Stim")

#control gcamp
allROIs.G.C<-length(unique(GCaMP.C$ROIs_trial))
fastROIs.G.C<-length(unique(GCaMP.C$ROIs_trial[GCaMP.C$Group=="fast_MDs"]))
delayedROIs.G.C<-length(unique(GCaMP.C$ROIs_trial[GCaMP.C$Group=="delayed_MDs"]))

propfast.G.C=fastROIs.G.C/allROIs.G.C
propdelayed.G.C=delayedROIs.G.C/allROIs.G.C

# atropine gcamp
allROIs.G.A<-length(unique(GCaMP.A$ROIs_trial))
fastROIs.G.A<-length(unique(GCaMP.A$ROIs_trial[GCaMP.A$Group=="fast_MDs"]))
delayedROIs.G.A<-length(unique(GCaMP.A$ROIs_trial[GCaMP.A$Group=="delayed_MDs"]))

propfast.G.A=fastROIs.G.A/allROIs.G.A
propdelayed.G.A=delayedROIs.G.A/allROIs.G.A

# prazosin gcamp
allROIs.G.P<-length(unique(GCaMP.P$ROIs_trial))
fastROIs.G.P<-length(unique(GCaMP.P$ROIs_trial[GCaMP.P$Group=="fast_MDs"]))
delayedROIs.G.P<-length(unique(GCaMP.P$ROIs_trial[GCaMP.P$Group=="delayed_MDs"]))

propfast.G.P=fastROIs.G.P/allROIs.G.P
propdelayed.G.P=delayedROIs.G.P/allROIs.G.P

# trazodone gcamp
allROIs.G.T<-length(unique(GCaMP.T$ROIs_trial))
fastROIs.G.T<-length(unique(GCaMP.T$ROIs_trial[GCaMP.T$Group=="fast_MDs"]))
delayedROIs.G.T<-length(unique(GCaMP.T$ROIs_trial[GCaMP.T$Group=="delayed_MDs"]))

propfast.G.T=fastROIs.G.T/allROIs.G.T
propdelayed.G.T=delayedROIs.G.T/allROIs.G.T

# metergoline gcamp
allROIs.G.M<-length(unique(GCaMP.M$ROIs_trial))
fastROIs.G.M<-length(unique(GCaMP.M$ROIs_trial[GCaMP.M$Group=="fast_MDs"]))
delayedROIs.G.M<-length(unique(GCaMP.M$ROIs_trial[GCaMP.M$Group=="delayed_MDs"]))

propfast.G.M=fastROIs.G.M/allROIs.G.M
propdelayed.G.M=delayedROIs.G.M/allROIs.G.M

# DSP4  gcamp
allROIs.G.D<-length(unique(GCaMP.D$ROIs_trial))
fastROIs.G.D<-length(unique(GCaMP.D$ROIs_trial[GCaMP.D$Group=="fast_MDs"]))
delayedROIs.G.D<-length(unique(GCaMP.D$ROIs_trial[GCaMP.D$Group=="delayed_MDs"]))

propfast.G.D=fastROIs.G.D/allROIs.G.D
propdelayed.G.D=delayedROIs.G.D/allROIs.G.D


# aggregate the spot info for all the peaks within the time window used to define fast/delayed

# Supplementary Fig.
# ROIwise for neurons- mean onset time, peak time, amplitude, across all trials
# calculate the difference between stim and each drug

#Spot wise- for neurons and astrocytes- stim and no stim
# mean num of ROIs, onset time, peak time, amplitude, AUC
#calculate the difference between stim and each drug

Spot.pharma.GC<-ddply(lck.peaks.window, 
                      c("Animal","Spot","Drug","Condition","Channel", "Animal_Spot"), 
                      summarise, nROIs=length(OnsetTime), meanOnset=mean(OnsetTime),
                      meanPeakT=mean(peakTime), meanAmp=mean(amplitude),
                      meanDur=mean(Duration), meanAUC1=mean(TraceAUC1),
                      meanAUC10=mean(TraceAUC10))

Spot.pharma.GC.group<-ddply(lck.peaks.window, 
                            c("Animal","Spot","Drug","Condition","Channel", "Animal_Spot","Group"), 
                            summarise, nROIs=length(OnsetTime), meanOnset=mean(OnsetTime),
                            meanPeakT=mean(peakTime), meanAmp=mean(amplitude),
                            meanDur=mean(Duration), meanAUC1=mean(TraceAUC1),
                            meanAUC10=mean(TraceAUC10))


Spot.pharma.GC.diff.OT2 <- dcast(Spot.pharma.GC, Animal + Spot + Animal_Spot + Channel + Condition ~ Drug, value.var="meanOnset")
Spot.pharma.GC.diff.OT2$norm_Control<- (Spot.pharma.GC.diff.OT2$Control/Spot.pharma.GC.diff.OT2$Control)
Spot.pharma.GC.diff.OT2$norm_Atropine<- (Spot.pharma.GC.diff.OT2$Atropine/Spot.pharma.GC.diff.OT2$Control)
Spot.pharma.GC.diff.OT2$norm_Metergoline<- (Spot.pharma.GC.diff.OT2$Metergoline/Spot.pharma.GC.diff.OT2$Control)
Spot.pharma.GC.diff.OT2$norm_Trazodone<- (Spot.pharma.GC.diff.OT2$Trazodone/Spot.pharma.GC.diff.OT2$Control)
Spot.pharma.GC.diff.OT2$norm_Prazosin<- (Spot.pharma.GC.diff.OT2$Prazosin/Spot.pharma.GC.diff.OT2$Control)


diff.data.short.OT <- Spot.pharma.GC.diff.OT2 [c("Animal","Spot","Animal_Spot","Channel","Condition","norm_Control",
                                                 "norm_Atropine", "norm_Metergoline", "norm_Trazodone", "norm_Prazosin")]
Spot.diff.OT<- melt(diff.data.short.OT)


ggplot(Spot.diff.OT, aes(x=interaction(variable, Condition), y=value, fill=Channel))+
  geom_boxplot()+
  ggtitle("onset time differences")+
  max.theme

df.OT.diff<-summarySE(Spot.diff.OT, measurevar = "value", groupvars = c("Channel", "Condition","variable"))
