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
#look at actaRCaMP vs endfeet

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

control.lck.peaks <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/Peaks_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
control.lck.OT<-read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/OnsetTimes_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")

pharmacology.lck.peaks <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/Peaks_pharmacology_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
pharmacology.lck.OT<-read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/OnsetTimes_pharmacology_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")

VSMC.peaks <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/Peaks_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")


##### 
#home files

control.lck.peaks <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/Peaks_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")
control.lck.OT<-read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/OnsetTimes_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")

##########

lsm.options(pbkrtest.limit = 100000)


control.lck.peaks<-subset(control.lck.peaks, Animal=="ARG2")
control.lck.OT<- subset(control.lck.OT, Animal=="ARG2")

pharmacology.lck.peaks<-subset(pharmacology.lck.peaks, Animal=="ARG2")
pharmacology.lck.OT<- subset(pharmacology.lck.OT, Animal=="ARG2")

control.lck.OT$ROIs_trial=NULL
control.lck.OT$Spot_trial=NULL
control.lck.OT$ROIType= NULL
control.lck.OT$overlap=control.lck.OT$Overlap
control.lck.OT$Overlap=NULL
control.lck.OT$Drug="Control"

all.lck.OT<-rbind(control.lck.OT, pharmacology.lck.OT)

#unique ROI names
all.lck.OT$ROIs_trial<-paste(all.lck.OT$Animal, all.lck.OT$Spot, all.lck.OT$Trial,all.lck.OT$ROI, sep= "_")
all.lck.OT$Spot_trial<-paste(all.lck.OT$Animal, all.lck.OT$Spot, all.lck.OT$Trial, sep= "_")
all.lck.OT$Spot_trial_Cond<-paste(all.lck.OT$Spot_trial, all.lck.OT$Condition, sep="_")
all.lck.OT$ROIs_Cond<-paste(all.lck.OT$ROIs_trial, all.lck.OT$Condition, sep="_")
all.lck.OT$ROIs_Cond_Drug<-paste(all.lck.OT$ROIs_trial, all.lck.OT$Condition,all.lck.OT$Drug, sep="_")


# REMOVE duplicate entries from onset time and peak time data frames 
# only the first entry will be used
#all.lck.OT2=all.lck.OT[order(all.lck.OT$OnsetTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
#all.lck.OT<-distinct(all.lck.OT2, ROIs_Cond_Drug, .keep_all = TRUE)

#remove entries with no onset time
all.lck.OT<-subset(all.lck.OT, OnsetTime!="NaN")

rm(all.lck.OT2)

control.lck.peaks$Drug="Control"

all.lck.peaks<-rbind(control.lck.peaks, pharmacology.lck.peaks)

all.lck.peaks$ROIType= "none"
all.lck.peaksA<- subset(all.lck.peaks, Channel=="GCaMP")
all.lck.peaksB<- subset(all.lck.peaks, Channel=="RCaMP")

# ROITypes
all.lck.peaksA$ROIType[grepl("r",all.lck.peaksA$roiName)]="Process"
all.lck.peaksA$ROIType[grepl("E",all.lck.peaksA$roiName)]="Endfoot"
all.lck.peaksB$ROIType[grepl("r",all.lck.peaksB$roiName)]="Dendrite"
all.lck.peaksB$ROIType[grepl("D",all.lck.peaksB$roiName)]="Dendrite"
all.lck.peaksB$ROIType[grepl("N",all.lck.peaksB$roiName)]="Neuron"
all.lck.peaksB$ROIType[grepl("S",all.lck.peaksB$roiName)]="SmoothMuscle"

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
# since there is the same ROI name for different drugs:
all.lck.peaks$ROIs_Cond_Drug<-paste(all.lck.peaks$ROIs_Cond, all.lck.peaks$Drug, sep="_")
all.lck.peaks$FullROIs<-paste(all.lck.peaks$Animal, all.lck.peaks$Spot, all.lck.peaks$roiName, sep="_")

# count number of trials per spot
Spot.lck.ntrials<-ddply(all.lck.peaks, c("Animal","Drug","Spot","Condition"), summarise, nTrials=length(unique(Trial)))


# remove ROIs with no peaks
all.lck.peaks<-all.lck.peaks[!(all.lck.peaks$peakType=="NoPeak"),]


# remove matching astrocyte process and soma ROIs
Overlap= all.lck.peaks$overlap!=0
all.lck.peaks<-all.lck.peaks[!Overlap,]
#OverlapROIs<-unique(nostim$ROIs_trial[Overlap])


#add baseline time to peaks table
all.lck.peaks$BL_time=10
all.lck.OT$BL_time=10


# adjust peak time and duration
all.lck.peaks$peakTime<- all.lck.peaks$peakTime-all.lck.peaks$BL_time
all.lck.peaks$peakStart<- all.lck.peaks$peakStart-all.lck.peaks$BL_time
all.lck.peaks$peakStartHalf<- all.lck.peaks$peakStartHalf-all.lck.peaks$BL_time
all.lck.peaks$Duration<- all.lck.peaks$halfWidth*2



# drop peaks that occur before the start of stimulation
#all.lck.peaks2<-subset(all.lck.peaks,peakTime>0)

# only the first entry will be used
#all.lck.peaks3<-all.lck.peaks2[order(all.lck.peaks2$peakTime),] # sort by ascending peak time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
#all.lck.peaks<-distinct(all.lck.peaks3, ROIs_Cond_Drug,.keep_all = TRUE)

rm(control.lck.OT, control.lck.peaks, pharmacology.lck.OT, pharmacology.lck.peaks, all.lck.peaks2, all.lck.peaks3)


# remove all data except astrocyte endfeet and smooth muscle
all.lck.peaks<-subset(all.lck.peaks, ROIType=="Endfoot" | ROIType=="SmoothMuscle")

# remove trazodone data
all.lck.peaks<-subset(all.lck.peaks, Drug!="Trazodone")
all.lck.OT<-subset(all.lck.OT, Drug!="Trazodone")

all.lck.OT$Drug<- factor(all.lck.OT$Drug, levels=c("Control","Atropine","Prazosin"))
all.lck.peaks$Drug<- factor(all.lck.peaks$Drug, levels=c("Control","Atropine","Prazosin"))


SpotsToInclude= c("17_11_14_spot2", "17_11_16_spot4", "17_11_20_spot5",
                  "17_11_21_spot6", "17_11_22_spot7")

all.lck.peaks<- subset(all.lck.peaks, Spot %in% SpotsToInclude)

####
# mean number of peaks per ROI

ROIwise<-ddply(all.lck.peaks, c("Animal","Spot","Drug","Condition","Channel","ROIs_Cond_Drug","FullROIs"), summarise, nPeaks=length(peakTime))


df.lck.peaknum<-summarySE(ROIwise, measurevar = "nPeaks", groupvars = c("Channel", "Drug","Condition"))


ggplot(df.lck.peaknum, aes(x=interaction(Drug,Channel),y=nPeaks, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=nPeaks-se, ymax=nPeaks+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("nPeaksPerTriall per field of view") +
  max.theme

# only RCaMP (smooth muscle)
RCaMP.ROIs<- subset(ROIwise, Channel=="RCaMP")
Condition_Drug= interaction(RCaMP.ROIs$Condition, RCaMP.ROIs$Drug)
nPeaks.lck.stim.null = lmer(nPeaks ~ (1|Spot) + (1|FullROIs), RCaMP.ROIs,REML=FALSE)
nPeaks.lck.stim.model2 = lmer(nPeaks ~ Condition + (1|Spot)+ (1|FullROIs), RCaMP.ROIs,REML=FALSE)
nPeaks.lck.stim.model3 = lmer(nPeaks ~ Drug +  (1|Spot)+ (1|FullROIs), RCaMP.ROIs,REML=FALSE)
nPeaks.lck.stim.model4 = lmer(nPeaks ~ Condition_Drug + (1|Spot)+ (1|FullROIs), RCaMP.ROIs,REML=FALSE)
nPeaks.lck.stim.anova <- anova(nPeaks.lck.stim.null, nPeaks.lck.stim.model2,nPeaks.lck.stim.model3,
                               nPeaks.lck.stim.model4)
print(nPeaks.lck.stim.anova)

nPeaks.lck.stim.Cond_Channel_Drug<- glht(nPeaks.lck.stim.model4, mcp(Condition_Drug= "Tukey"))
summary(nPeaks.lck.stim.Cond_Channel_Drug)



########
#amplitude
df.amp1<-summarySE(all.lck.peaks, measurevar = "amplitude", groupvars = c("Channel","Condition", "Drug"))
df.amp2<-summarySE(stim.lck.alldata, measurevar = "amplitude", groupvars = c("Channel", "ROIType","Condition", "Drug"))



ggplot(df.amp1, aes(x=interaction(Channel,Drug),y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
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

# calculate the change in amplitude, onset time etc. with drug since it is the same cells in each FOV



