
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

########################
#Correlation of RCaMP and GCaMP Signals during long stim
#CorrData <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/LongStim_Correlations.csv", header=TRUE, sep = ",")
CorrData <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LongStim_Correlations.csv", header=TRUE, sep = ",")

CorrData$ROI_Y_type<-as.character(CorrData$ROI_Y_type)
CorrData$ROI_Y_type[grepl("D",CorrData$ROI_Y)]="Dendrite"
CorrData$ROI_Y_type<-as.factor(CorrData$ROI_Y_type)
CorrData$CompChannel<-paste(CorrData$ChannelX, CorrData$ChannelY, sep= "_")
CorrData$CompType<-paste(CorrData$ROI_X_type, CorrData$ROI_Y_type, sep= "_")


ggplot(CorrData, aes(x=Long_Corr, fill=CompType)) + geom_histogram(binwidth=0.05, position="dodge") +
  ggtitle("Distribution of Correlation- All comparisons")

ggplot(CorrData, aes(x=Long_Corr, fill=CompChannel)) + geom_histogram(binwidth=0.05, position="dodge") +
  ggtitle("Distribution of Correlation- All comparisons")

ggplot(CorrData, aes(x=Short_Corr, fill=CompChannel)) + geom_histogram(binwidth=0.05, position="dodge") +
  ggtitle("Distribution of stim window Correlation- All comparisons")

ggplot(CorrData, aes(x=Lag, fill=CompChannel)) + geom_histogram(binwidth=0.5, position="dodge") +
  ggtitle("Distribution of lags- All comparisons")

ggplot(CorrData, aes(x=xCorr, fill=CompChannel)) + geom_histogram(binwidth=0.1, position="dodge") +
  ggtitle("Distribution of CrossCorrelations- All comparisons")

fisher.r2z <- function(r) { 0.5 * (log(1+r) - log(1-r)) }
round(fisher.r2z(seq(0,1,0.01)), 4)

CorrData$Fisher=fisher.r2z(CorrData$Short_Corr)

df1A<- summarySE(CorrData, measurevar="Short_Corr", groupvars=c("CompChannel"))
df1B<- summarySE(CorrData, measurevar="xCorr", groupvars=c("CompChannel"))
df1C<- summarySE(CorrData, measurevar="Lag", groupvars=c("CompChannel"))


#####
# consider ONLY RCaMP and GCaMP Comparisons

GCaMP_RCaMP<-subset(CorrData, CompChannel=="GCaMP_RCaMP")
GCaMP_RCaMP$ROIs_trial<-paste(GCaMP_RCaMP$Animal, GCaMP_RCaMP$Spot, GCaMP_RCaMP$Trial,GCaMP_RCaMP$ROI_X, sep= "_")
#remove random spots with gcamp neurons
ntypes=c("Neuron_Neuron","Neuron_Neuropil","Neuropil_Neuron","Neuropil_Neuropil")
GCaMP_RCaMP<-subset(GCaMP_RCaMP, !(CompType %in% ntypes))

ggplot(GCaMP_RCaMP, aes(x=MinDistance, y=xCorr,fill=CompType)) + geom_point(aes(colour=CompType))

df1A<- summarySE(GCaMP_RCaMP, measurevar="Short_Corr", groupvars=c("CompType"))
df1B<- summarySE(GCaMP_RCaMP, measurevar="xCorr", groupvars=c("CompType"))

ggplot(data=df1A, aes(x=CompType, y=Short_Corr, fill=CompType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("CompType") +
  ylab("Short_Corr") +
  ggtitle("Short_Corr for gcamp rcamp in all trials") +
max.theme

ggplot(data=df1B, aes(x=CompType, y=xCorr, fill=CompType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("CompType") +
  ylab("xCorr") +
  ggtitle("xCorr for gcamp rcamp in all trials") +
  max.theme

#####
#higher correlated data during stimulation window

#sort by cross correlation
highcorr<- subset(GCaMP_RCaMP, xCorr>0.5)

df2A<- summarySE(highcorr, measurevar="Short_Corr", groupvars=c("CompType"))
df2B<- summarySE(highcorr, measurevar="Lag", groupvars=c("CompType"))

ggplot(highcorr, aes(x=xCorr, fill=CompType)) + geom_histogram(binwidth=0.1, position="dodge") +
  ggtitle("Distribution of CrossCorrelations- highly correlated AC and N")

ggplot(highcorr, aes(x=MinDistance, fill=CompType)) + geom_histogram(binwidth=5, position="dodge") +
  ggtitle("Distribution of distance- highly correlated AC and N")

ggplot(highcorr, aes(x=Lag, fill=CompType)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of lags- highly correlated AC and N")

ggplot(highcorr, aes(x=MinDistance, y=xCorr,fill=CompType)) + geom_point(aes(colour=CompType))

##################
# sort by lag

zerolag<-subset(GCaMP_RCaMP, Lag==0)
zerolagROIs<-unique(zerolag$ROIs_trial)

ggplot(zerolag, aes(x=MinDistance, fill=CompType)) + geom_histogram(binwidth=5, position="dodge") +
  ggtitle("Distribution of distance- zerolag rois")

ggplot(zerolag, aes(x=MinDistance, y=xCorr,fill=CompType)) + geom_point(aes(colour=CompType))
ggplot(zerolag, aes(x=MinDistance, y=Short_Corr,fill=CompType)) + geom_point(aes(colour=CompType))


#####
# extract data for "early" astrocytes identified from the traces
# compare to "late" astrocytes (all other responding ROIs)

responders<-longstim.stimwindow

#earlyAUC<-read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/EarlyGC_byAUC.csv", header=TRUE, sep = ",")
earlyAUC<-read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/EarlyGC_byAUC.csv", header=TRUE, sep = ",")


respondingGC<-subset(responders, Channel=="GCaMP")
respondingRC<-subset(responders, Channel=="RCaMP")

GCaMP_RCaMP$RCaMP_ROIs<-paste(GCaMP_RCaMP$Animal, GCaMP_RCaMP$Spot, GCaMP_RCaMP$Trial,GCaMP_RCaMP$ROI_Y, sep= "_")

respGC.corr<-subset(GCaMP_RCaMP, ROIs_trial %in% unique(respondingGC$ROIs_trial))
respGC.corr<-subset(respGC.corr, RCaMP_ROIs %in% respondingRC$ROIs_trial)

Corr_EarlyAUC<-subset(GCaMP_RCaMP,  ROIs_trial %in% earlyAUC$names)

Corr_LateAUC<-subset(respGC.corr,  !(ROIs_trial %in% earlyAUC$names))

Corr_EarlyAUC$Group<-"Early"
Corr_LateAUC$Group<-"Late"


GCaMP.CorrAUC<-rbind(Corr_EarlyAUC,Corr_LateAUC)

GCaMP.CorrAUC$Group<-as.factor(GCaMP.CorrAUC$Group)


df3A<- summarySE(GCaMP.CorrAUC, measurevar="Short_Corr", groupvars=c("Group"))
df3B<- summarySE(GCaMP.CorrAUC, measurevar="Short_Corr", groupvars=c("Group","CompType"))
df3C<- summarySE(GCaMP.CorrAUC, measurevar="Fisher", groupvars=c("Group"))

ggplot(data=df3A, aes(x=Group, y=Short_Corr, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Group") +
  ylab("Short_Corr") +
  ggtitle("Short_Corr for responding GCaMP") +
  max.theme

ggplot(data=df3B , aes(x=CompType, y=Short_Corr, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Group") +
  ylab("Short_Corr") +
  ggtitle("Short_Corr for responding GCaMP") +
  max.theme

ggplot(data=df3C, aes(x=Group, y=Fisher, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Fisher-se, ymax=Fisher+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Group") +
  ylab("Fisher") +
  ggtitle("Fisher for responding GCaMP") +
  max.theme


df4A<- summarySE(GCaMP.CorrAUC, measurevar="xCorr", groupvars=c("Group"))
df4B<- summarySE(GCaMP.CorrAUC, measurevar="xCorr", groupvars=c("Group","CompType"))


ggplot(data=df4A, aes(x=Group, y=xCorr, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Group") +
  ylab("xCorr") +
  ggtitle("xCorr for responding GCaMP") +
  max.theme

ggplot(data=df4B , aes(x=CompType, y=xCorr, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Group") +
  ylab("xCorr") +
  ggtitle("xCorr for responding GCaMP") +
  max.theme




df5A<- summarySE(GCaMP.CorrAUC, measurevar="Lag", groupvars=c("Group"))
df5B<- summarySE(GCaMP.CorrAUC, measurevar="Lag", groupvars=c("Group","CompType"))

ggplot(data=df5A, aes(x=Group, y=Lag, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Group") +
  ylab("Lag") +
  ggtitle("Lag for responding GCaMP") +
  max.theme

ggplot(data=df5B , aes(x=CompType, y=Lag, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Group") +
  ylab("Lag") +
  ggtitle("Lag for responding GCaMP") +
  max.theme


#########
#stats
# Likelihood-ratio test 

# ROI by AUC
Group_CompType=interaction(GCaMP.CorrAUC$Group, GCaMP.CorrAUC$CompType)

#correlation of only the stim window traces (short corr)
shortCorr.null = lmer(Short_Corr ~ (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)
shortCorr.model1 = lmer(Short_Corr ~ Group + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)
shortCorr.model2 = lmer(Short_Corr ~ CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)
shortCorr.model3 = lmer(Short_Corr ~ CompType+Group + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)
shortCorr.model4 = lmer(Short_Corr ~ Group_CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)

shortCorr.anova <- anova(shortCorr.null,shortCorr.model1,shortCorr.model2,shortCorr.model3,shortCorr.model4)
print(shortCorr.anova)

# p values
shortCorr.Group <- glht(shortCorr.model1, mcp(Group= "Tukey"))
summary(shortCorr.Group)

shortCorr.Group_types <- glht(shortCorr.model4, mcp(Group_CompType= "Tukey"))
summary(shortCorr.Group_types)


#cross correlation
xCorr.null = lmer(xCorr ~ (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)
xCorr.model1 = lmer(xCorr ~ Group + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)
xCorr.model2 = lmer(xCorr ~ CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)
xCorr.model3 = lmer(xCorr ~ CompType+Group + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)
xCorr.model4 = lmer(xCorr ~ Group_CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)

xCorr.anova <- anova(xCorr.null,xCorr.model1,xCorr.model2,xCorr.model3,xCorr.model4)
print(xCorr.anova)

# p values
xCorr.Group <- glht(xCorr.model1, mcp(Group= "Tukey"))
summary(xCorr.Group)

xCorr.Group_types <- glht(xCorr.model4, mcp(Group_CompType= "Tukey"))
summary(xCorr.Group_types)

#lag
Lag.null = lmer(Lag ~ (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)
Lag.model1 = lmer(Lag ~ Group + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)
Lag.model2 = lmer(Lag ~ CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)
Lag.model3 = lmer(Lag ~ CompType+Group + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)
Lag.model4 = lmer(Lag ~ Group_CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)

Lag.anova <- anova(Lag.null,Lag.model1,Lag.model2,Lag.model3,Lag.model4)
print(Lag.anova)

# p values
Lag.Group <- glht(Lag.model1, mcp(Group= "Tukey"))
summary(Lag.Group)

Lag.Group_types <- glht(Lag.model4, mcp(Group_CompType= "Tukey"))
summary(Lag.Group_types)


#####
#distance

df6A<- summarySE(GCaMP.CorrAUC, measurevar="MinDistance", groupvars=c("Group"))
df6B<- summarySE(GCaMP.CorrAUC, measurevar="MinDistance", groupvars=c("Group","CompType"))

ggplot(data=df6A, aes(x=Group, y=MinDistance, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=MinDistance-se, ymax=MinDistance+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Group") +
  ylab("MinDistance") +
  ggtitle("MinDistance for responding GCaMP") +
  max.theme

ggplot(data=df6B , aes(x=CompType, y=MinDistance, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=MinDistance-se, ymax=MinDistance+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Group") +
  ylab("MinDistance") +
  ggtitle("MinDistance for responding GCaMP") +
  max.theme



#MinDistance
Dis.null = lmer(MinDistance ~ (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)
Dis.model1 = lmer(MinDistance ~ Group + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)
Dis.model2 = lmer(MinDistance ~ CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)
Dis.model3 = lmer(MinDistance ~ CompType+Group + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)
Dis.model4 = lmer(MinDistance ~ Group_CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP.CorrAUC,REML=FALSE)

Dis.anova <- anova(Dis.null,Dis.model1,Dis.model2,Dis.model3,Dis.model4)
print(Dis.anova)

# p values
Dis.Group <- glht(Dis.model1, mcp(Group= "Tukey"))
summary(Dis.Group)

Dis.Group_types <- glht(Dis.model4, mcp(Group_CompType= "Tukey"))
summary(Dis.Group_types)

#####
ggplot(data=Corr_EarlyAUC, aes(x=Lag, fill=CompType)) +
  geom_density(alpha=.2) 
# box plot
ggplot(data=Corr_EarlyAUC, aes(x=CompType, y=Lag, fill=CompType)) +
  geom_boxplot() 

# cumulative fraction plot
ggplot(Corr_EarlyAUC, aes(x=Lag, colour = CompType)) + stat_ecdf()

######

