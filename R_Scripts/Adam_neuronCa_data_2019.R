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
#IP3R2KO

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

neuron1<- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/FilesforR/Peaks_allMice_Lck_nostim_vs_longstim_12_2017.csv", header=TRUE, sep = ",")


# add animal name, FOV name, trial name, genotype name

##########

lsm.options(pbkrtest.limit = 100000)

# join data sets
all.lck.peaks<-rbind(lck.peaks1,lck.peaks2)


# remove the data frames that are combined
rm(lck.peaks1,lck.peaks2,lck.OT1, lck.OT2)

# only consider IP3R2KO
all.lck.peaks<-all.lck.peaks[grepl("IP",all.lck.peaks$Animal),]

all.lck.OTB$ROIType[grepl("r",all.lck.OTB$ROI)]="Dendrite"

rm(all.lck.OT2,all.lck.OTA, all.lck.OTB)


#unique ROI names
all.lck.peaks$ROIs_trial<-paste(all.lck.peaks$Animal, all.lck.peaks$Spot, all.lck.peaks$Trial,all.lck.peaks$roiName, sep= "_")


# count number of trials per spot
Spot.lck.ntrials<-ddply(all.lck.peaks, c("Animal","Genotype","Spot","Condition"), summarise, nTrials=length(unique(Trial)))
Spot.lck.ntrials$Ani_Spot_Cond<-paste(Spot.lck.ntrials$Animal, Spot.lck.ntrials$Spot, Spot.lck.ntrials$Condition, sep="_")


# remove ROIs with no peaks
all.lck.peaks<-all.lck.peaks[!(all.lck.peaks$peakType=="NoPeak"),]


# adjust peak time and duration
all.lck.peaks$peakTime<- all.lck.peaks$peakTime-all.lck.peaks$BL_time
all.lck.peaks$peakStart<- all.lck.peaks$peakStart-all.lck.peaks$BL_time
all.lck.peaks$peakStartHalf<- all.lck.peaks$peakStartHalf-all.lck.peaks$BL_time
all.lck.peaks$Duration<- all.lck.peaks$halfWidth*2



# drop peaks that occur before the start of stimulation
all.lck.peaks2<-subset(all.lck.peaks,peakTime>0)


# plot some distributions

#peak times
ggplot(lck.peaks.window,y=peakTime, fill= Condition) +
  geom_boxplot()+
  ylab("peak Time (s)") +
  ggtitle("time window 0-12.76s, rcamp, stim")+
  max.theme


ggplot(lck.peaks.window,y=peakAUC, fill= Condition) +
  geom_boxplot()+
  ylab("peak Time (s)") +
  ggtitle("time window 0-12.76s, rcamp, stim")+
  max.theme

ggplot(lck.peaks.window,y=amplitude, fill= Condition) +
  geom_boxplot()+
  ylab("peak Time (s)") +
  ggtitle("time window 0-12.76s, rcamp, stim")+
  max.theme



#############
# mean peak time
# mean onset times
df.PT1<- summarySE(stim.lck.compdata[stim.lck.compdata$Condition=="Stim",], measurevar = "peakTime", groupvars = c("Channel_Group", "Genotype"))
df.PT2<- summarySE(stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",], measurevar = "peakTime", groupvars = c("Group", "Genotype"))
df.PT3<- summarySE(stim.lck.compdata.STIM, measurevar = "peakTime", groupvars = c("ROIType","Channel_Group", "Genotype"))

ggplot(df.PT1, aes(x=interaction(Genotype,Channel_Group),y=peakTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Peak Max Time (s)") +
  max.theme

ggplot(df.PT2, aes(x=Group,y=peakTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Peak Max Time (s)") +
  max.theme

ggplot(df.PT3, aes(x=interaction(Channel_Group,ROIType),y=peakTime, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Mean Peak Max Time (s)") +
  max.theme



# stats
# stats for onset times- neurons vs astrocytes
PT.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials), lck.peaks.window,REML=FALSE)
PT.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot) + (1|trials), lck.peaks.window,REML=FALSE)
PT.model2 = lmer(peakTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), lck.peaks.window,REML=FALSE)
PT.model3 = lmer(peakTime ~ Group_Channel_Cond_Gen + (1|Animal) + (1|Spot) + (1|trials), lck.peaks.window,REML=FALSE)
PT.model4 = lmer(peakTime ~ Group_Channel_Type_Cond_Gen + (1|Animal) + (1|Spot) + (1|trials), lck.peaks.window,REML=FALSE)
PT.anova <- anova(PT.null, PT.model1,PT.model2,PT.model3,PT.model4)
print(PT.anova)

PT.Group_channel_Gen<- glht(PT.model3, mcp(Group_Channel_Cond_Gen= "Tukey"))
summary(PT.Group_channel_Gen)

#PT.Group_Channel_Type_Gen<- glht(PT.model4, mcp(Group_Channel_Type_Cond_Gen= "Tukey"))
#summary(PT.Group_Channel_Type_Gen)


# stats for onset times- only ASTROCYTES
PT.stim.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
PT.stim.model1 = lmer(peakTime ~ Genotype + (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
PT.stim.model3 = lmer(peakTime ~ Group + (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
PT.stim.model4 = lmer(peakTime ~ Group_Gen + (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
PT.stim.model6 = lmer(peakTime ~ Group_Type_Gen + (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
PT.stim.anova <- anova(PT.stim.null, PT.stim.model1,PT.stim.model3,PT.stim.model4,PT.stim.model6)
print(PT.stim.anova)

PT.stim.Group_Gen<- glht(PT.stim.model4, mcp(Group_Gen= "Tukey"))
summary(PT.stim.Group_Gen)

#PT.stim.Group_type_Gen<- glht(PT.stim.model6, mcp(Group_Type_Gen= "Tukey"))
#summary(PT.stim.Group_type_Gen)


########
#amplitude
df.amp1<-summarySE(lck.peaks.window, measurevar = "amplitude", groupvars = c("Channel","Condition", "Genotype"))
df.amp2<-summarySE(lck.peaks.window, measurevar = "amplitude", groupvars = c("Channel", "ROIType","Condition", "Genotype"))

df.amp3A<- summarySE(stim.lck.compdata, measurevar = "amplitude", groupvars = c("Channel_Group","Genotype","Condition"))
df.amp3B<- summarySE(stim.lck.compdata.STIM, measurevar = "amplitude", groupvars = c("Channel_Group","Genotype"))

df.amp4<- summarySE(astro.compdata.stim, measurevar = "amplitude", groupvars = c("Group","Genotype"))
df.amp5<- summarySE(astro.compdata.stim, measurevar = "amplitude", groupvars = c("Group","ROIType","Genotype"))


ggplot(df.amp1, aes(x=interaction(Channel,Genotype),y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  max.theme

ggplot(df.amp2, aes(x=interaction(Channel,interaction(Genotype, ROIType)),y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  max.theme

ggplot(df.amp3A, aes(x=interaction(Channel_Group,Genotype),y=amplitude, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  max.theme

ggplot(df.amp3B, aes(x=Channel_Group,y=amplitude, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude during stim trials") +
  max.theme

ggplot(df.amp4, aes(x=Group,y=amplitude, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude during stim trials") +
  max.theme


ggplot(df.amp5, aes(x=interaction(Group, ROIType),y=amplitude, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude during stim trials") +
  max.theme





#lck-GCaMP ONLY
amp.null.GC = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
amp.model2.GC  = lmer(amplitude ~ Genotype + (1|Animal) + (1|Spot) + (1|trials) , astro.compdata.stim,REML=FALSE)
amp.model3.GC  = lmer(amplitude ~ Group + (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
amp.model5.GC  = lmer(amplitude ~ Group_Gen + (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
amp.model6.GC  = lmer(amplitude ~ Group_Type_Gen + (1|Animal) + (1|Spot) + (1|trials), astro.compdata.stim,REML=FALSE)
amp.anova.GC  <- anova(amp.null.GC, amp.model2.GC,amp.model3.GC,amp.model5.GC,amp.model6.GC)
print(amp.anova.GC)

amp.Group_Gen.GC<- glht(amp.model5.GC, mcp(Group_Gen= "Tukey"))
summary(amp.Group_Gen.GC)

amp.Group_gen_ty.GC<- glht(amp.model6.GC, mcp(Group_Type_Gen= "Tukey"))
summary(amp.Group_gen_ty.GC)


# KNOCKOUTS have significantly lower amplitude that WT!
# still a significant difference between fast and delayed WTs



########
#duration
df.Dur1<-summarySE(lck.peaks.window, measurevar = "Duration", groupvars = c("Channel","Condition", "Genotype"))
df.Dur2<-summarySE(lck.peaks.window, measurevar = "Duration", groupvars = c("Channel", "ROIType","Condition", "Genotype"))

df.Dur3A<- summarySE(stim.lck.compdata, measurevar = "Duration", groupvars = c("Channel_Group","Genotype","Condition"))
df.Dur3B<- summarySE(stim.lck.compdata.STIM, measurevar = "Duration", groupvars = c("Channel_Group","Genotype"))

df.Dur4<- summarySE(astro.compdata.stim, measurevar = "Duration", groupvars = c("Group","Genotype"))
df.Dur5<- summarySE(astro.compdata.stim, measurevar = "Duration", groupvars = c("Group","ROIType","Genotype"))


ggplot(df.Dur1, aes(x=interaction(Channel,Genotype),y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  max.theme

ggplot(df.Dur2, aes(x=interaction(Channel,interaction(Genotype, ROIType)),y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  max.theme

ggplot(df.Dur3A, aes(x=interaction(Channel_Group,Genotype),y=Duration, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  max.theme

ggplot(df.Dur3B, aes(x=Channel_Group,y=Duration, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration during stim trials") +
  max.theme

ggplot(df.Dur4, aes(x=Group,y=Duration, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration during stim trials") +
  max.theme

ggplot(df.Dur5, aes(x=interaction(Group, ROIType),y=Duration, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration during stim trials") +
  max.theme





#lck-GCaMP ONLY
Group_Genotype=interaction(stim.lck.compdata.STIM$Group[stim.lck.compdata.STIM$Channel=="GCaMP"],stim.lck.compdata.STIM$Genotype[stim.lck.compdata.STIM$Channel=="GCaMP"])
Group_Type_Genotype=interaction(stim.lck.compdata.STIM$Group[stim.lck.compdata.STIM$Channel=="GCaMP"],
                                stim.lck.compdata.STIM$Genotype[stim.lck.compdata.STIM$Channel=="GCaMP"],
                                stim.lck.compdata.STIM$ROIType[stim.lck.compdata.STIM$Channel=="GCaMP"])

Dur.null.GC = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
Dur.model2.GC  = lmer(Duration ~ Genotype + (1|Animal) + (1|Spot) + (1|trials) , stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
Dur.model3.GC  = lmer(Duration ~ Group + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
Dur.model5.GC  = lmer(Duration ~ Group_Genotype + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
Dur.model6.GC  = lmer(Duration ~ Group_Type_Genotype + (1|Animal) + (1|Spot) + (1|trials), stim.lck.compdata.STIM[stim.lck.compdata.STIM$Channel=="GCaMP",],REML=FALSE)
Dur.anova.GC  <- anova(Dur.null.GC, Dur.model2.GC,Dur.model3.GC,Dur.model5.GC,Dur.model6.GC)
print(Dur.anova.GC)

Dur.Group_Gen.GC<- glht(Dur.model5.GC, mcp(Group_Genotype= "Tukey"))
summary(Dur.Group_Gen.GC)

Dur.Group_gen_ty.GC<- glht(Dur.model6.GC, mcp(Group_Type_Genotype= "Tukey"))
summary(Dur.Group_gen_ty.GC)


#######
#peakAUC
df.peakAUC1<-summarySE(lck.peaks.window, measurevar = "peakAUC", groupvars = c("Channel","Condition", "Genotype"), na.rm=TRUE)
df.peakAUC2<-summarySE(lck.peaks.window, measurevar = "peakAUC", groupvars = c("Channel", "ROIType","Condition", "Genotype"),na.rm=TRUE)

df.peakAUC3A<- summarySE(stim.lck.compdata, measurevar = "peakAUC", groupvars = c("Channel_Group","Genotype","Condition"),na.rm=TRUE)
df.peakAUC3B<- summarySE(stim.lck.compdata.STIM, measurevar = "peakAUC", groupvars = c("Channel_Group","Genotype"),na.rm=TRUE)

df.peakAUC4<- summarySE(astro.compdata.stim, measurevar = "peakAUC", groupvars = c("Group","Genotype"),na.rm=TRUE)
df.peakAUC5<- summarySE(astro.compdata.stim, measurevar = "peakAUC", groupvars = c("Group","ROIType","Genotype"),na.rm=TRUE)


ggplot(df.peakAUC1, aes(x=interaction(Channel,Genotype),y=peakAUC, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakAUC-se, ymax=peakAUC+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakAUC") +
  max.theme

ggplot(df.peakAUC2, aes(x=interaction(Channel,interaction(Genotype, ROIType)),y=peakAUC, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakAUC-se, ymax=peakAUC+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakAUC") +
  max.theme

ggplot(df.peakAUC3A, aes(x=interaction(Channel_Group,Genotype),y=peakAUC, fill= Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakAUC-se, ymax=peakAUC+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakAUC") +
  max.theme

ggplot(df.peakAUC3B, aes(x=Channel_Group,y=peakAUC, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakAUC-se, ymax=peakAUC+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakAUC peakAUCing stim trials") +
  max.theme

ggplot(df.peakAUC4, aes(x=Group,y=peakAUC, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakAUC-se, ymax=peakAUC+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakAUC during stim trials") +
  max.theme

ggplot(df.peakAUC5, aes(x=interaction(Group, ROIType),y=peakAUC, fill= Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakAUC-se, ymax=peakAUC+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakAUC during stim trials") +
  max.theme


