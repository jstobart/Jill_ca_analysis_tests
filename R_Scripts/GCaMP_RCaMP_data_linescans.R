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

# fast vs delayed
all.lck.OT$Channel_Group<-interaction(all.lck.OT$Channel, all.lck.OT$Group)
all.lck.OT$Channel_Group<-as.factor(all.lck.OT$Channel_Group)

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

ggplot(control.stim[control.stim$OnsetTime<15,],aes(x=OnsetTime,y=..density..,fill=Channel)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP vs RCaMP onset times")+
  max.theme
#########
# mean onset times
df.OT1<- summarySE(control.window.stim, measurevar = "OnsetTime", groupvars = c("Channel"))
df.OT2<- summarySE(control.all.stim.compdata, measurevar = "OnsetTime", groupvars = c("Channel_Group", "Channel"))

summarySE(stim.lck.compdata.STIM, measurevar = "OnsetTime", groupvars = c("ROIType","Channel_Group", "Drug"))

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



