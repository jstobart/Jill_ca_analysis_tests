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
# linescans with genotypes and Pharmacology
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

all.lck.OT<- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/LinescanOnsets_allMice_Lck_nostim_vs_longstim_01_2018.csv", header=TRUE, sep = ",")
#linescans2<- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/LinescanOnsets_2ndCohort_Lck_nostim_vs_longstim_01_2018.csv", header=TRUE, sep = ",")

##### 
#home files

all.lck.OT<-read.table("D:/Data/GCaMP_RCaMP/Revision/Lck_GCaMP/FilesforR/Linescans/LinescanOnsets_allMice_Lck_nostim_vs_longstim_01_2018.csv", header=TRUE, sep = ",")


##########

lsm.options(pbkrtest.limit = 100000)

#all.lck.OT<-rbind(linescans1, linescans2)

#unique ROI names
all.lck.OT$ROIs_trial<-paste(all.lck.OT$Animal, all.lck.OT$Spot, all.lck.OT$Trial,all.lck.OT$ROI, sep= "_")
all.lck.OT$Spot_trial<-paste(all.lck.OT$Animal, all.lck.OT$Spot, all.lck.OT$Trial, sep= "_")
all.lck.OT$Spot_trial_Cond<-paste(all.lck.OT$Spot_trial, all.lck.OT$Condition, sep="_")
all.lck.OT$ROIs_Cond<-paste(all.lck.OT$ROIs_trial, all.lck.OT$Condition, sep="_")


# REMOVE ROIs with no onset time detected
all.lck.OT<-subset(all.lck.OT, !(is.na(OnsetTime)))

# REMOVE duplicate entries from the same ROI (i.e. multiple peaks with different onset times) 
all.lck.OT2=all.lck.OT[order(all.lck.OT$OnsetTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain)
all.lck.OT<-distinct(all.lck.OT2, ROIs_Cond, .keep_all = TRUE)

rm(all.lck.OT2, linescans1, linescans2)


# count number of trials per spot
Spot.lck.ntrials<-ddply(all.lck.OT, c("Animal","drug","Genotype", "Spot","Condition"), summarise, nTrials=length(unique(Trial)))
Spot.lck.ntrials$Ani_Spot_Cond<-paste(Spot.lck.ntrials$Animal, Spot.lck.ntrials$Spot, Spot.lck.ntrials$Condition, sep="_")

#add baseline time to peaks table
all.lck.OT$BL_time<-10
all.lck.OT$Channel<-factor(all.lck.OT$Channel, levels= c("RCaMP","GCaMP"))

######

# identify "FAST" astrocytes
all.lck.OT$Group<-0
all.lck.OT$Group[all.lck.OT$OnsetTime<1]<-"fast_MDs"
all.lck.OT$Group[all.lck.OT$OnsetTime>=1]<-"delayed_MDs"

all.lck.OT$Group <- factor(all.lck.OT$Group, levels = c("fast_MDs","delayed_MDs"))
all.lck.OT$drug<- factor(all.lck.OT$drug, levels=c("Control","Atropine","Prazosin","Trazodone"))

# fast vs delayed
all.lck.OT$Channel_Group<-interaction(all.lck.OT$Channel, all.lck.OT$Group)
all.lck.OT$Channel_Group<-as.factor(all.lck.OT$Channel_Group)

#percentile distribution on neurons & astrocytes

NeuronalStim<-subset(all.lck.OT, Channel=="RCaMP" & Condition=="Stim")
Neuron95Onset<-quantile(NeuronalStim$OnsetTime[NeuronalStim$OnsetTime<8], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
print(Neuron95Onset)

AstrocyteStim<-subset(all.lck.OT, Channel=="GCaMP" & Condition=="Stim")
Astrocyte95Onset<-quantile(AstrocyteStim$OnsetTime[AstrocyteStim$OnsetTime<8], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
print(Astrocyte95Onset)

######
# for median and mean calculations

# no stim vs 8 s stim- 
#neuronal window=2s, AC window=12 s for onset

LongN_OTwind2=2
LongAC_OTwind2=12

# remove data that is outside the above windows
# lck
stim.lck.OT.R<-subset(all.lck.OT, Channel=="RCaMP" & OnsetTime<=LongN_OTwind2)
stim.lck.OT.G<-subset(all.lck.OT, Channel=="GCaMP" & OnsetTime<=LongAC_OTwind2)

stim.lck.OT.window<-rbind(stim.lck.OT.R, stim.lck.OT.G)


rm(stim.lck.OT.G,stim.lck.OT.R)


# removed delayed neurons
stim.lck.OT.window.compdata<-stim.lck.OT.window[!(stim.lck.OT.window$Channel=="RCaMP"& stim.lck.OT.window$Group=="delayed_MDs"),]




#######################
# astrocytes vs. neurons onset times only from stim trials

# data from wildtype and IP3R2_WT
all.lck.OT$Genotype[all.lck.OT$Animal=="IPRG7"]="IP3R2_KO"
all.lck.OT$Genotype[all.lck.OT$Animal=="IPRG6"]="IP3R2_WT"

control.all<-subset(all.lck.OT, drug=="Control" & Genotype!="IP3R2_KO")
control.all.stim<-subset(all.lck.OT, drug=="Control" & Genotype!="IP3R2_KO" & Condition=="Stim")
#control.window.stim.compdata<-subset(stim.lck.OT.window.compdata, drug=="Control" & Genotype!="IP3R2_KO" & Condition=="Stim")

#control.window.stim<-subset(stim.lck.OT.window, drug=="Control" & Genotype!="IP3R2_KO" & Condition=="Stim")

#onset time histogram
ggplot(control.all.stim[control.all.stim$OnsetTime<20,],aes(x=OnsetTime,y=..density..,fill=Channel)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP vs RCaMP onset times")+
  max.theme

ggplot(control.window.stim.compdata[control.window.stim.compdata$Group!="delayed_MDs",], 
       aes(x=OnsetTime,fill=Channel)) +
  geom_histogram(position="dodge") +
  ggtitle("GCaMP vs RCaMP onset times")+
  max.theme

######
# number of peaks detected in nostim and stim
control.OT.20strial<-control.all[control.all$OnsetTime<20,]
control.OT.20strial.stim<-control.all.stim[control.all.stim$OnsetTime<20,]

control.OT.20strial$Channel <- factor(control.OT.20strial$Channel, levels = c("RCaMP","GCaMP"))

peakNum.OT.20strial<-ddply(control.OT.20strial, c("Animal","Spot","Condition","Channel","Group"), summarise, nPeaks=length(OnsetTime))

# add in number of trials
peakNum.OT.20strial$Ani_Spot_Cond<-paste(peakNum.OT.20strial$Animal, peakNum.OT.20strial$Spot, peakNum.OT.20strial$Condition, sep="_")
peakNum.OT.20strial<-merge(peakNum.OT.20strial, Spot.lck.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)
peakNum.OT.20strial$PeaksPerTrial<-peakNum.OT.20strial$nPeaks/peakNum.OT.20strial$nTrials

peakNum.OT.20strial.compdata<-peakNum.OT.20strial[!(peakNum.OT.20strial$Channel=="RCaMP"& peakNum.OT.20strial$Group=="delayed_MDs"),]

# mean
df.lck.peakNum.20strial<-summarySE(peakNum.OT.20strial, measurevar = "PeaksPerTrial", groupvars = c("Channel", "Condition"))
df.lck.peakNum.20strial.Group<-summarySE(peakNum.OT.20strial.compdata, measurevar = "PeaksPerTrial", groupvars = c("Channel", "Condition","Group"))

ggplot(df.lck.peakNum.20strial, aes(x=Channel,y=PeaksPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=PeaksPerTrial-se, ymax=PeaksPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num Peaks/trial per field of view") +
  max.theme

ggplot(df.lck.peakNum.20strial.Group, aes(x=interaction(Channel,Group),y=PeaksPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=PeaksPerTrial-se, ymax=PeaksPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num Peaks/trial per field of view, Group") +
  max.theme


Condition_Channel2= interaction(peakNum.OT.20strial.compdata$Condition,peakNum.OT.20strial.compdata$Channel)
Condition_Channel_Group= interaction(peakNum.OT.20strial.compdata$Condition,peakNum.OT.20strial.compdata$Channel,peakNum.OT.20strial.compdata$Group)

nPeaks.lck.stim.null = lmer(PeaksPerTrial ~ (1|Animal), peakNum.OT.20strial.compdata,REML=FALSE)
nPeaks.lck.stim.model1 = lmer(PeaksPerTrial~ Channel + (1|Animal), peakNum.OT.20strial.compdata,REML=FALSE)
nPeaks.lck.stim.model2 = lmer(PeaksPerTrial ~ Condition + (1|Animal), peakNum.OT.20strial.compdata,REML=FALSE)
nPeaks.lck.stim.model3 = lmer(PeaksPerTrial ~ Condition_Channel2 + (1|Animal), peakNum.OT.20strial.compdata,REML=FALSE)
nPeaks.lck.stim.model4 = lmer(PeaksPerTrial ~ Condition_Channel_Group + (1|Animal), peakNum.OT.20strial.compdata,REML=FALSE)
nPeaks.lck.stim.anova <- anova(nPeaks.lck.stim.null, nPeaks.lck.stim.model1,nPeaks.lck.stim.model2,
                               nPeaks.lck.stim.model3,nPeaks.lck.stim.model4)
print(nPeaks.lck.stim.anova)

nPeaks.lck.stim.Cond_Channel_Group<- glht(nPeaks.lck.stim.model4, mcp(Condition_Channel_Group= "Tukey"))
summary(nPeaks.lck.stim.Cond_Channel_Group)


# not significant between stim and no stim for any group

#########
# mean onset times- control

control.OT.20strial.compdata<-control.OT.20strial.stim[!(control.OT.20strial.stim$Channel=="RCaMP"
                                                    & control.OT.20strial.stim$Group=="delayed_MDs"),]


df.OT1<- summarySE(control.OT.20strial.stim, measurevar = "OnsetTime", groupvars = c("Channel","Condition"))
df.OT2<- summarySE(control.OT.20strial.compdata, measurevar = "OnsetTime", groupvars = c("Channel_Group", "Channel"))
df.OT2.fast<-summarySE(control.OT.20strial.compdata[control.OT.20strial.compdata$Group!="delayed_MDs",], 
                       measurevar = "OnsetTime", groupvars = c("Channel_Group", "Channel"))

ggplot(df.OT1, aes(x=Channel,y=OnsetTime, fill=Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)- neurons vs astros in stim window") +
  ggtitle("Mean Onset Time (s)- neurons vs astros in stim window") +
  max.theme

ggplot(df.OT2, aes(x=Channel_Group,y=OnsetTime, fill=Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme

ggplot(df.OT2.fast, aes(x=Channel_Group,y=OnsetTime, fill=Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme



# stats for onset times- neurons vs astrocytes
OT.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial), control.OT.20strial.compdata,REML=FALSE)
OT.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), control.OT.20strial.compdata,REML=FALSE)
OT.model2 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|Spot_trial), control.OT.20strial.compdata,REML=FALSE)
OT.model3 = lmer(OnsetTime ~ Channel_Group + (1|Animal) + (1|Spot) + (1|Spot_trial), control.OT.20strial.compdata,REML=FALSE)
OT.anova <- anova(OT.null, OT.model1,OT.model2,OT.model3)
print(OT.anova)

OT.Channel_Group<- glht(OT.model3, mcp(Channel_Group= "Tukey"))
summary(OT.Channel_Group)


summary(OT.model3)

# check residuals for linearity
plot(fitted(OT.model3), residuals(OT.model3),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(OT.model3), residuals(OT.model3)), col=46, lwd=2.5)








#######################
# IP3R2 knockouts

IP3.all<-subset(all.lck.OT, drug=="Control" & Genotype!="NaN")
IP3.all.stim<-subset(all.lck.OT, drug=="Control" & Genotype!="NaN" & Condition=="Stim")
#IP3.window.stim<-subset(stim.lck.OT.window, drug=="Control" & Genotype!="NaN" & Condition=="Stim")
#IP3.window.stim.compdata<-subset(stim.lck.OT.window.compdata, drug=="Control" & Genotype!="NaN" & Condition=="Stim")

#onset time histogram
ggplot(IP3.all.stim[(IP3.all.stim$OnsetTime<15),],aes(x=OnsetTime,y=..density..,fill=interaction(Channel,Genotype))) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP vs RCaMP onset times- IP3 WTs")+
  max.theme

ggplot(IP3.all.stim[(IP3.all.stim$OnsetTime<15 & IP3.all.stim$Genotype=="IP3R2_KO"),],aes(x=OnsetTime,y=..density..,fill=interaction(Channel,Genotype))) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP vs RCaMP onset times- IP3 KOs")+
  max.theme


#########
# number of peaks detected in nostim and stim
IP3.OT.20strial<-IP3.all[IP3.all$OnsetTime<20,]

IP3.OT.20strial$Channel <- factor(IP3.OT.20strial$Channel, levels = c("RCaMP","GCaMP"))

peakNum.IP3.20strial<-ddply(IP3.OT.20strial, c("Animal","Spot","Condition","Channel","Group","Genotype"), summarise, nPeaks=length(OnsetTime))

# add in number of trials
peakNum.IP3.20strial$Ani_Spot_Cond<-paste(peakNum.IP3.20strial$Animal, peakNum.IP3.20strial$Spot, peakNum.IP3.20strial$Condition, sep="_")
peakNum.IP3.20strial<-merge(peakNum.IP3.20strial, Spot.lck.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)
peakNum.IP3.20strial$PeaksPerTrial<-peakNum.IP3.20strial$nPeaks/peakNum.IP3.20strial$nTrials

peakNum.IP3.20strial.compdata<-peakNum.IP3.20strial[!(peakNum.IP3.20strial$Channel=="RCaMP"& peakNum.IP3.20strial$Group=="delayed_MDs"),]

# mean
df.IP3.peakNum.20strial<-summarySE(peakNum.IP3.20strial, measurevar = "PeaksPerTrial", groupvars = c("Channel", "Condition","Genotype"))
df.IP3.peakNum.20strial.Group<-summarySE(peakNum.IP3.20strial.compdata, measurevar = "PeaksPerTrial", groupvars = c("Channel", "Condition","Group","Genotype"))

ggplot(df.IP3.peakNum.20strial, aes(x=interaction(Channel, Genotype),y=PeaksPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=PeaksPerTrial-se, ymax=PeaksPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num Peaks/trial per field of view") +
  max.theme

ggplot(df.IP3.peakNum.20strial.Group, aes(x=interaction(Channel,Group, Genotype),y=PeaksPerTrial, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=PeaksPerTrial-se, ymax=PeaksPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num Peaks/trial per field of view, Group") +
  max.theme

ggplot(df.IP3.peakNum.20strial.Group[df.IP3.peakNum.20strial.Group$Condition=="Stim",], aes(x=interaction(Channel,Group),y=PeaksPerTrial, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=PeaksPerTrial-se, ymax=PeaksPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num Peaks/trial per field of view, Group") +
  max.theme

Condition_Channel_Gen= interaction(peakNum.IP3.20strial.compdata$Condition,
                                peakNum.IP3.20strial.compdata$Channel,
                                peakNum.IP3.20strial.compdata$Genotype)
Condition_Channel_Group_Gen= interaction(peakNum.IP3.20strial.compdata$Condition,
                                     peakNum.IP3.20strial.compdata$Channel,
                                     peakNum.IP3.20strial.compdata$Group,
                                     peakNum.IP3.20strial.compdata$Genotype)

nPeaks.IP3.stim.null = lmer(PeaksPerTrial ~ (1|Animal), peakNum.IP3.20strial.compdata,REML=FALSE)
nPeaks.IP3.stim.model1 = lmer(PeaksPerTrial~ Channel + (1|Animal), peakNum.IP3.20strial.compdata,REML=FALSE)
nPeaks.IP3.stim.model2A= lmer(PeaksPerTrial ~ Condition + (1|Animal), peakNum.IP3.20strial.compdata,REML=FALSE)
nPeaks.IP3.stim.model2B = lmer(PeaksPerTrial ~ Genotype + (1|Animal), peakNum.IP3.20strial.compdata,REML=FALSE)
nPeaks.IP3.stim.model3 = lmer(PeaksPerTrial ~ Condition_Channel_Gen + (1|Animal), peakNum.IP3.20strial.compdata,REML=FALSE)
nPeaks.IP3.stim.model4 = lmer(PeaksPerTrial ~ Condition_Channel_Group_Gen + (1|Animal), peakNum.IP3.20strial.compdata,REML=FALSE)
nPeaks.IP3.stim.anova <- anova(nPeaks.IP3.stim.null, nPeaks.IP3.stim.model1,nPeaks.IP3.stim.model2A,
                               nPeaks.IP3.stim.model2B, nPeaks.IP3.stim.model3,nPeaks.IP3.stim.model4)
print(nPeaks.IP3.stim.anova)

nPeaks.IP3.stim.Cond_Channel_Group<- glht(nPeaks.IP3.stim.model4, mcp(Condition_Channel_Group_Gen= "Tukey"))
summary(nPeaks.IP3.stim.Cond_Channel_Group)


# not significant between stim and no stim for any group


##########
# mean onset times- IP3
IP3.OT.20strial.stim<-IP3.OT.20strial[IP3.OT.20strial$Condition=="Stim",]
IP3.OT.20strial.compdata<-IP3.OT.20strial.stim[!(IP3.OT.20strial.stim$Channel=="RCaMP"
                                                    & IP3.OT.20strial.stim$Group=="delayed_MDs"),]

df.IP3.OT1<- summarySE(IP3.OT.20strial.compdata, measurevar = "OnsetTime", groupvars = c("Channel", "Genotype"))
df.IP3.OT2<- summarySE(IP3.OT.20strial.compdata, measurevar = "OnsetTime", groupvars = c("Channel_Group", "Channel", "Genotype"))
df.IP3.OT2.fast<-summarySE(IP3.OT.20strial.compdata[IP3.OT.20strial.compdata$Group!="delayed_MDs",], 
                       measurevar = "OnsetTime", groupvars = c("Channel_Group", "Channel", "Genotype"))

ggplot(df.IP3.OT1, aes(x=Channel,y=OnsetTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  ggtitle("Mean Onset Time (s)- neurons vs astros in stim window") +
  max.theme

ggplot(df.IP3.OT2, aes(x=Channel_Group,y=OnsetTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme

ggplot(df.IP3.OT2.fast, aes(x=Channel_Group,y=OnsetTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme


# stats for onset times- neurons vs astrocytes
Channel_Group_Genotype= interaction(IP3.OT.20strial.compdata$Channel, IP3.OT.20strial.compdata$Group, IP3.OT.20strial.compdata$Genotype)
OT.IP3.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial), IP3.OT.20strial.compdata,REML=FALSE)
OT.IP3.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), IP3.OT.20strial.compdata,REML=FALSE)
OT.IP3.model2 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|Spot_trial), IP3.OT.20strial.compdata,REML=FALSE)
OT.IP3.model3 = lmer(OnsetTime ~ Genotype + (1|Animal) + (1|Spot) + (1|Spot_trial), IP3.OT.20strial.compdata,REML=FALSE)
OT.IP3.model4 = lmer(OnsetTime ~ Channel_Group_Genotype + (1|Animal) + (1|Spot) + (1|Spot_trial), IP3.OT.20strial.compdata,REML=FALSE)
OT.IP3.anova <- anova(OT.IP3.null, OT.IP3.model1, OT.IP3.model2, OT.IP3.model3, OT.IP3.model4)
print(OT.IP3.anova)

OT.IP3.Channel_Group<- glht(OT.IP3.model4, mcp(Channel_Group_Genotype= "Tukey"))
summary(OT.IP3.Channel_Group)

# check residuals for linearity
plot(fitted(OT.IP3.model4), residuals(OT.model3),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(OT.IP3.model4), residuals(OT.IP3.model4)), col=46, lwd=2.5)




#######################
# Pharmacology

Pharm.all.stim<-subset(all.lck.OT, Genotype!="IP3R2_KO" & Condition=="Stim")
Pharm.window.stim<-subset(stim.lck.OT.window, Genotype!="IP3R2_KO" & Condition=="Stim")
Pharm.window.stim.compdata<-subset(stim.lck.OT.window.compdata, Genotype!="IP3R2_KO" & Condition=="Stim")

#onset time histogram
ggplot(Pharm.all.stim[(Pharm.all.stim$OnsetTime<15 & Pharm.all.stim$drug=="IP3"),],aes(x=OnsetTime,y=..density..,fill=interaction(Channel,drug))) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP vs RCaMP onset times- IP3")+
  max.theme

ggplot(Pharm.all.stim[(Pharm.all.stim$OnsetTime<15 & Pharm.all.stim$drug=="Atropine"),],aes(x=OnsetTime,y=..density..,fill=interaction(Channel,drug))) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP vs RCaMP onset times- Atropine")+
  max.theme

ggplot(Pharm.all.stim[(Pharm.all.stim$OnsetTime<15 & Pharm.all.stim$drug=="Prazosin"),],aes(x=OnsetTime,y=..density..,fill=interaction(Channel,drug))) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP vs RCaMP onset times- Prazosin")+
  max.theme

ggplot(Pharm.all.stim[(Pharm.all.stim$OnsetTime<15 & Pharm.all.stim$drug=="Trazodone"),],aes(x=OnsetTime,y=..density..,fill=interaction(Channel,drug))) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP vs RCaMP onset times- Trazadone")+
  max.theme

ggplot(Pharm.all.stim[(Pharm.all.stim$OnsetTime<15 & Pharm.all.stim$drug=="Meterogoline"),],aes(x=OnsetTime,y=..density..,fill=interaction(Channel,drug))) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP vs RCaMP onset times- Meterogoline")+
  max.theme

#########
# mean onset times- Pharm
df.OT1<- summarySE(Pharm.window.stim, measurevar = "OnsetTime", groupvars = c("Channel", "drug"))
df.OT2<- summarySE(Pharm.window.stim.compdata, measurevar = "OnsetTime", groupvars = c("Channel_Group", "Channel", "drug"))
df.OT2.fast<-summarySE(Pharm.window.stim.compdata[Pharm.window.stim.compdata$Group!="delayed_MDs",], 
                       measurevar = "OnsetTime", groupvars = c("Channel_Group", "Channel", "drug"))

ggplot(df.OT1, aes(x=Channel,y=OnsetTime, fill=drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  ggtitle("Mean Onset Time (s)- neurons vs astros in stim window") +
  max.theme

ggplot(df.OT2, aes(x=Channel_Group,y=OnsetTime, fill=drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme

ggplot(df.OT2.fast, aes(x=Channel_Group,y=OnsetTime, fill=drug)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme


Channel_Group_drug= interaction(Pharm.window.stim.compdata$Channel, Pharm.window.stim.compdata$Group, Pharm.window.stim.compdata$drug)
# stats for onset times- neurons vs astrocytes
OT.Pharm.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial), Pharm.window.stim.compdata,REML=FALSE)
OT.Pharm.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), Pharm.window.stim.compdata,REML=FALSE)
OT.Pharm.model2 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|Spot_trial), Pharm.window.stim.compdata,REML=FALSE)
OT.Pharm.model3 = lmer(OnsetTime ~ drug + (1|Animal) + (1|Spot) + (1|Spot_trial), Pharm.window.stim.compdata,REML=FALSE)
OT.Pharm.model4 = lmer(OnsetTime ~ Channel_Group + drug + (1|Animal) + (1|Spot) + (1|Spot_trial), Pharm.window.stim.compdata,REML=FALSE)
OT.Pharm.model5 = lmer(OnsetTime ~ Channel_Group_drug + (1|Animal) + (1|Spot) + (1|Spot_trial), Pharm.window.stim.compdata,REML=FALSE)
OT.Pharm.anova <- anova(OT.Pharm.null, OT.Pharm.model1, OT.Pharm.model2, OT.Pharm.model3,
                      OT.Pharm.model4, OT.Pharm.model5)
print(OT.Pharm.anova)

OT.Channel_Group_drug<- glht(OT.Pharm.model5, mcp(Channel_Group_drug= "Tukey"))
summary(OT.Channel_Group_drug)


summary(OT.Pharm.model5)

# check residuals for linearity
plot(fitted(OT.Pharm.model5), residuals(OT.Pharm.model5),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(OT.Pharm.model5), residuals(OT.Pharm.model5)), col=46, lwd=2.5)

