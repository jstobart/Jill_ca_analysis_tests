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

all.lck.OT<- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/LinescanOnsets_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")

##### 
#home files

all.lck.OT<-read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/LinescanOnsets_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")


##########

lsm.options(pbkrtest.limit = 100000)

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

#rm(all.lck.OT2)


# count number of trials per spot
Spot.lck.ntrials<-ddply(all.lck.OT, c("Animal","drug","Genotype", "Spot","Condition"), summarise, nTrials=length(unique(Trial)))

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

#percentile distribution on neuronsl

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

control.all.stim<-subset(all.lck.OT, drug=="Control" & Genotype!="IP3R2_KO" & Condition=="Stim")
control.window.stim.compdata<-subset(stim.lck.OT.window.compdata, drug=="Control" & Genotype!="IP3R2_KO" & Condition=="Stim")

control.window.stim<-subset(stim.lck.OT.window, drug=="Control" & Genotype!="IP3R2_KO" & Condition=="Stim")

#onset time histogram
ggplot(control.all.stim[control.all.stim$OnsetTime<15,],aes(x=OnsetTime,y=..density..,fill=Channel)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP vs RCaMP onset times")+
  max.theme

ggplot(control.window.stim.compdata[control.window.stim.compdata$Group!="delayed_MDs",], 
       aes(x=OnsetTime,fill=Channel)) +
  geom_histogram(position="dodge") +
  ggtitle("GCaMP vs RCaMP onset times")+
  max.theme

#########
# mean onset times- control
df.OT1<- summarySE(control.window.stim, measurevar = "OnsetTime", groupvars = c("Channel"))
df.OT2<- summarySE(control.window.stim.compdata, measurevar = "OnsetTime", groupvars = c("Channel_Group", "Channel"))
df.OT2.fast<-summarySE(control.window.stim.compdata[control.window.stim.compdata$Group!="delayed_MDs",], 
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
OT.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial), control.window.stim.compdata,REML=FALSE)
OT.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), control.window.stim.compdata,REML=FALSE)
OT.model2 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|Spot_trial), control.window.stim.compdata,REML=FALSE)
OT.model3 = lmer(OnsetTime ~ Channel_Group + (1|Animal) + (1|Spot) + (1|Spot_trial), control.window.stim.compdata,REML=FALSE)
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

IP3.all.stim<-subset(all.lck.OT, drug=="Control" & Genotype!="NaN" & Condition=="Stim")
IP3.window.stim<-subset(stim.lck.OT.window, drug=="Control" & Genotype!="NaN" & Condition=="Stim")
IP3.window.stim.compdata<-subset(stim.lck.OT.window.compdata, drug=="Control" & Genotype!="NaN" & Condition=="Stim")

#onset time histogram
ggplot(IP3.all.stim[(IP3.all.stim$OnsetTime<15 & IP3.all.stim$Genotype=="IP3R2_WT"),],aes(x=OnsetTime,y=..density..,fill=interaction(Channel,Genotype))) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP vs RCaMP onset times- IP3 WTs")+
  max.theme

ggplot(IP3.all.stim[(IP3.all.stim$OnsetTime<15 & IP3.all.stim$Genotype=="IP3R2_KO"),],aes(x=OnsetTime,y=..density..,fill=interaction(Channel,Genotype))) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP vs RCaMP onset times- IP3 KOs")+
  max.theme


#########
# mean onset times- IP3
df.OT1<- summarySE(IP3.window.stim, measurevar = "OnsetTime", groupvars = c("Channel", "Genotype"))
df.OT2<- summarySE(IP3.window.stim.compdata, measurevar = "OnsetTime", groupvars = c("Channel_Group", "Channel", "Genotype"))
df.OT2.fast<-summarySE(IP3.window.stim.compdata[IP3.window.stim.compdata$Group!="delayed_MDs",], 
                       measurevar = "OnsetTime", groupvars = c("Channel_Group", "Channel", "Genotype"))

ggplot(df.OT1, aes(x=Channel,y=OnsetTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  ggtitle("Mean Onset Time (s)- neurons vs astros in stim window") +
  max.theme

ggplot(df.OT2, aes(x=Channel_Group,y=OnsetTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme

ggplot(df.OT2.fast, aes(x=Channel_Group,y=OnsetTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Onset Time (s)") +
  max.theme


Channel_Group_Genotype= interaction(IP3.window.stim.compdata$Channel, IP3.window.stim.compdata$Group, IP3.window.stim.compdata$Genotype)
# stats for onset times- neurons vs astrocytes
OT.IP3.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot) + (1|Spot_trial), IP3.window.stim.compdata,REML=FALSE)
OT.IP3.model1 = lmer(OnsetTime ~ Channel + (1|Animal) + (1|Spot) + (1|Spot_trial), IP3.window.stim.compdata,REML=FALSE)
OT.IP3.model2 = lmer(OnsetTime ~ Group + (1|Animal) + (1|Spot) + (1|Spot_trial), IP3.window.stim.compdata,REML=FALSE)
OT.IP3.model3 = lmer(OnsetTime ~ Genotype + (1|Animal) + (1|Spot) + (1|Spot_trial), IP3.window.stim.compdata,REML=FALSE)
OT.IP3.model4 = lmer(OnsetTime ~ Channel_Group + Genotype + (1|Animal) + (1|Spot) + (1|Spot_trial), IP3.window.stim.compdata,REML=FALSE)
OT.IP3.model5 = lmer(OnsetTime ~ Channel_Group_Genotype + (1|Animal) + (1|Spot) + (1|Spot_trial), IP3.window.stim.compdata,REML=FALSE)
OT.IP3.anova <- anova(OT.IP3.null, OT.IP3.model1, OT.IP3.model2, OT.IP3.model3,
                      OT.IP3.model4, OT.IP3.model5)
print(OT.IP3.anova)

OT.Channel_Group_Genotype<- glht(OT.IP3.model5, mcp(Channel_Group_Genotype= "Tukey"))
summary(OT.Channel_Group_Genotype)


summary(OT.IP3.model5)

# check residuals for linearity
plot(fitted(OT.IP3.model5), residuals(OT.IP3.model5),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(OT.IP3.model5), residuals(OT.IP3.model5)), col=46, lwd=2.5)





#######################
# Pharmacology

Pharm.all.stim<-subset(all.lck.OT, Genotype!="IP3R2_KO" & Condition=="Stim")
Pharm.window.stim<-subset(stim.lck.OT.window, Genotype!="IP3R2_KO" & Condition=="Stim")
Pharm.window.stim.compdata<-subset(stim.lck.OT.window.compdata, Genotype!="IP3R2_KO" & Condition=="Stim")

#onset time histogram
ggplot(Pharm.all.stim[(Pharm.all.stim$OnsetTime<15 & Pharm.all.stim$drug=="Control"),],aes(x=OnsetTime,y=..density..,fill=interaction(Channel,drug))) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP vs RCaMP onset times- Control")+
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

