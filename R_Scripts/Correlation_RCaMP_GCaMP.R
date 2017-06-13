
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

########################
#Correlation of RCaMP and GCaMP Signals during long stim
longstim.corr <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LongStim_Correlations_fixedDis.csv", header=TRUE, sep = ",")
shortstim.corr <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/shortstim_Correlations_fixedDis.csv", header=TRUE, sep = ",")

longstim.corr <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LongStim_Correlations_fixedDis.csv", header=TRUE, sep = ",")
shortstim.corr  <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/shortstim_Correlations_fixedDis.csv", header=TRUE, sep = ",")

CorrData<- rbind(longstim.corr, shortstim.corr)
#CorrData<- shortstim.corr



CorrData$ROI_Y_type<-as.character(CorrData$ROI_Y_type)
CorrData$ROI_Y_type[grepl("D",CorrData$ROI_Y)]="Dendrite"
CorrData$ROI_Y_type<-as.factor(CorrData$ROI_Y_type)
CorrData$CompChannel<-paste(CorrData$ChannelX, CorrData$ChannelY, sep= "_")
CorrData$CompType<-paste(CorrData$ROI_X_type, CorrData$ROI_Y_type, sep= "_")


ggplot(CorrData, aes(x=Long_Corr, fill=CompType)) + geom_density(alpha=0.5)+
  ggtitle("Distribution of Correlation- All comparisons") + max.theme

ggplot(CorrData, aes(x=Short_Corr, fill=CompType)) + geom_density(alpha=0.5)+
  ggtitle("Distribution of Correlation- All comparisons")  + max.theme

ggplot(CorrData, aes(x=Long_Corr, fill=CompChannel)) + geom_histogram(binwidth=0.05, position="dodge") +
  ggtitle("Distribution of Correlation- All comparisons") + max.theme

ggplot(CorrData, aes(x=Short_Corr, fill=CompChannel)) + geom_density(alpha=0.5) +
  ggtitle("Distribution of stim window Correlation- All comparisons") + max.theme

ggplot(CorrData, aes(x=xCorr, fill=CompChannel)) + geom_density(alpha=0.5) +
  ggtitle("Distribution of stim window Correlation- All comparisons") + max.theme

ggplot(CorrData, aes(x=Lag, fill=CompChannel)) + geom_density(alpha=0.5) +
  ggtitle("Distribution of stim window Correlation- All comparisons") + max.theme




# to fix non-normally distributed data
fisher.r2z <- function(r) { 0.5 * (log(1+r) - log(1-r)) }
round(fisher.r2z(seq(0,1,0.01)), 4)

CorrData$FisherShort=fisher.r2z(CorrData$Short_Corr)
CorrData$FisherxCorr=fisher.r2z(CorrData$xCorr)

ggplot(CorrData, aes(x=FisherxCorr, fill=CompChannel)) + geom_density(alpha=0.2) +
  ggtitle("Distribution of Fisher xCorr values") + max.theme



df1A<- summarySE(CorrData, measurevar="Short_Corr", groupvars=c("CompChannel"))
df1B<- summarySE(CorrData, measurevar="FisherShort", groupvars=c("CompChannel"))
df1C<- summarySE(CorrData, measurevar="xCorr", groupvars=c("CompChannel"))
df1D<- summarySE(CorrData, measurevar="FisherxCorr", groupvars=c("CompChannel"))
df1E<- summarySE(CorrData, measurevar="Lag", groupvars=c("CompChannel"))


#####
# consider ONLY RCaMP and GCaMP Comparisons

#histograms
ggplot(CorrData[CorrData$CompChannel=="GCaMP_RCaMP",], aes(x=Short_Corr, y=..density..,fill=CompChannel)) + geom_histogram(binwidth=0.03, position="dodge") +
  ggtitle("short corr-gcamp vs rcamp") + max.theme

ggplot(CorrData[CorrData$CompChannel=="GCaMP_RCaMP",], aes(x=Short_Corr, fill=CompChannel)) + geom_density()+
  ggtitle("density-") + max.theme

ggplot(CorrData[CorrData$CompChannel=="GCaMP_RCaMP",], aes(x=xCorr, y=..density..,fill=CompChannel)) + geom_histogram(binwidth=0.03, position="dodge") +
  ggtitle("xcorr-gcamp vs rcamp") + max.theme

ggplot(CorrData[CorrData$CompChannel=="GCaMP_RCaMP",], aes(x=Lag, y=..density..,fill=CompChannel)) + geom_histogram(binwidth=0.5, position="dodge") +
  ggtitle("xcorr-gcamp vs rcamp") + max.theme



# consider only GCaMP and RCaMP ROIs that respond to stim in the stim window
GCaMP_RCaMP<-subset(CorrData, CompChannel=="GCaMP_RCaMP")
GCaMP_RCaMP$ROIs_trial<-paste(GCaMP_RCaMP$Animal, GCaMP_RCaMP$Spot, GCaMP_RCaMP$Trial,GCaMP_RCaMP$ROI_X, sep= "_")

#remove random spots with gcamp neurons
ntypes=c("Neuron_Neuron","Neuron_Neuropil","Neuropil_Neuron","Neuropil_Neuropil")
GCaMP_RCaMP<-subset(GCaMP_RCaMP, !(CompType %in% ntypes))



#####
#higher correlated data during stimulation window

#sort by cross correlation
highcorr<- subset(GCaMP_RCaMP, Short_Corr>0.3)

dfhighcorrA<- summarySE(highcorr, measurevar="Short_Corr", groupvars=c("CompType"))
dfhighcorrB<- summarySE(highcorr, measurevar="Lag", groupvars=c("CompType"))

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

responding<-subset(all.lck.peaks.group, Group=="fast"| Group=="delayed")
#respondingRC<-subset(all.lck.peaks.group, (Group=="fast"| Group=="delayed") & Channel=="RCaMP")

GCaMP_RCaMP$RCaMP_ROIs<-paste(GCaMP_RCaMP$Animal, GCaMP_RCaMP$Spot, GCaMP_RCaMP$Trial,GCaMP_RCaMP$ROI_Y, sep= "_")

#make a list of gcamp ROIs (fast vs delayed)

#unique(df[c("yad", "per")])

respGC<-subset(responding, Channel=="GCaMP") 
respRC<-subset(responding, Channel=="RCaMP") 

# list of responding astrocytes and their corresponding group
respGC2<-unique(respGC[c("ROIs_trial","Group")])
respRC2<-unique(respRC[c("ROIs_trial","Group")])

# only correlations from responding astrocytes and neurons
GCaMP_RCaMP<-subset(GCaMP_RCaMP, ROIs_trial %in% respGC2$ROIs_trial)
GCaMP_RCaMP<-subset(GCaMP_RCaMP, RCaMP_ROIs %in% unique(respRC$ROIs_trial))

# put group information into correlation table ("fast or delayed")
GCaMP_RCaMP.group=data.frame()
for (ii in 1:nrow(respGC2))
{
  ROIx=respGC2$ROIs_trial[ii]
  Groupx=respGC2$Group[ii]
  Groupx=paste(Groupx,"A")
  subset1=subset(GCaMP_RCaMP, ROIs_trial==ROIx)
  if (nrow(subset1)>0)
  {
    subset1$GroupX=Groupx
    GCaMP_RCaMP.group<-rbind(GCaMP_RCaMP.group, subset1)
  }
}

# RCaMP group info
GCaMP_RCaMP.groups=data.frame()
for (ii in 1:nrow(respRC2))
{
  ROIy=respRC2$ROIs_trial[ii]
  Groupy=respRC2$Group[ii]
  Groupy=paste(Groupy,"N")
  subset1=subset(GCaMP_RCaMP.group, RCaMP_ROIs==ROIy)
  if (nrow(subset1)>0)
  {
    subset1$GroupY=Groupy
    GCaMP_RCaMP.groups<-rbind(GCaMP_RCaMP.groups, subset1)
  }
}



# enter responding ROI information- high, mid or low responding neurons
GCaMP_RCaMP.groups$Nresponders=0
GCaMP_RCaMP.Nresp=data.frame()
for (ii in 1:nrow(NeuronDendrites))
{
  Dend=NeuronDendrites$ROIs_trial[ii]
  RespType=NeuronDendrites$responders[ii]
  subset1=subset(GCaMP_RCaMP.groups, RCaMP_ROIs==Dend)
  if (nrow(subset1)>0)
  {
    subset1$Nresponders=RespType
    GCaMP_RCaMP.Nresp<-rbind(GCaMP_RCaMP.Nresp, subset1)
  }
}

GCaMP_RCaMP.groups$ROInameUnique<-paste(GCaMP_RCaMP.groups$Animal,GCaMP_RCaMP.groups$Spot,GCaMP_RCaMP.groups$ROI_Y, sep="_")
GCaMP_RCaMP.Nresp2=data.frame()
for (ii in 1:nrow(NeuronSomas))
{
  Soma=NeuronSomas$ROInameUnique[ii]
  RespType=NeuronSomas$responders[ii]
  subset1=subset(GCaMP_RCaMP.groups, ROInameUnique==Soma)
  if (nrow(subset1)>0)
  {
    subset1$Nresponders=RespType
    subset1$ROInameUnique=NULL
    GCaMP_RCaMP.Nresp2<-rbind(GCaMP_RCaMP.Nresp2, subset1)
  }
}

GCaMP_RCaMP.Nresp$ROInameUnique=NULL

GCaMP_RCaMP.Nresp<-rbind(GCaMP_RCaMP.Nresp,GCaMP_RCaMP.Nresp2)

#####
# histograms


ggplot(GCaMP_RCaMP.Nresp, aes(x=Short_Corr, fill=interaction(Nresponders,GroupX))) + geom_density(alpha=0.3)+
  ggtitle("density-") + max.theme

ggplot(GCaMP_RCaMP.Nresp[GCaMP_RCaMP.Nresp$GroupX=="fast A",], aes(x=Short_Corr, fill=Condition)) + geom_density(alpha=0.3)+
  ggtitle("density-") + max.theme

#change order for plots
GCaMP_RCaMP.Nresp$GroupX<- as.factor(GCaMP_RCaMP.Nresp$GroupX)
GCaMP_RCaMP.Nresp$GroupY<- as.factor(GCaMP_RCaMP.Nresp$GroupY)
GCaMP_RCaMP.Nresp$Nresponders<- as.factor(GCaMP_RCaMP.Nresp$Nresponders)

GCaMP_RCaMP.Nresp$CompType <- factor(GCaMP_RCaMP.Nresp$CompType, levels = c("Endfeet_Neuron","Endfeet_Dendrite", "Process_Neuron","Process_Dendrite"))
GCaMP_RCaMP.Nresp$GroupX <- factor(GCaMP_RCaMP.Nresp$GroupX, levels = c("fast A","delayed A"))
GCaMP_RCaMP.Nresp$GroupY <- factor(GCaMP_RCaMP.Nresp$GroupY, levels = c("fast N","delayed N"))
GCaMP_RCaMP.Nresp$Nresponders <- factor(GCaMP_RCaMP.Nresp$Nresponders, levels = c("low","mid", "high"))


# bar graphs
df4A<- summarySE(GCaMP_RCaMP.Nresp, measurevar="Short_Corr", groupvars=c("CompType"))
df4B<- summarySE(GCaMP_RCaMP.Nresp, measurevar="Short_Corr", groupvars=c("GroupX"))
df4C<- summarySE(GCaMP_RCaMP.Nresp, measurevar="Short_Corr", groupvars=c("GroupX","Condition"))
df4E<- summarySE(GCaMP_RCaMP.Nresp, measurevar="Short_Corr", groupvars=c("Nresponders"))
df4F<- summarySE(GCaMP_RCaMP.Nresp, measurevar="Short_Corr", groupvars=c("GroupX","Nresponders"))
df4G<- summarySE(GCaMP_RCaMP.Nresp, measurevar="Short_Corr", groupvars=c("GroupX","GroupY"))
df4G2<- summarySE(GCaMP_RCaMP.Nresp, measurevar="Short_Corr", groupvars=c("GroupX","GroupY","CompType"))


ggplot(data=df4B, aes(x=GroupX, y=Short_Corr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Short_Corr for responding GCaMP") +
  max.theme

ggplot(data=df4C, aes(x=Condition, y=Short_Corr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Short_Corr for responding GCaMP") +
  max.theme

ggplot(data=df4A, aes(x=CompType, y=Short_Corr, fill=CompType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Astrocyte group vs neuron group for responding GCaMP") +
  max.theme

ggplot(data=df4E, aes(x=Nresponders, y=Short_Corr, fill=Nresponders)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Astrocyte group vs neuron group for responding GCaMP") +
  max.theme

ggplot(data=df4F, aes(x=Nresponders, y=Short_Corr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Astrocyte group vs neuron group for responding GCaMP") +
  max.theme

ggplot(data=df4G, aes(x=GroupY, y=Short_Corr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Astrocyte group vs neuron group for responding GCaMP") +
  max.theme

ggplot(data=df4G2, aes(x=interaction(GroupY,CompType), y=Short_Corr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Astrocyte group vs neuron group for responding GCaMP") +
  max.theme

ggplot(data=df4G[df4G$GroupY=="fast N",], aes(x=GroupX, y=Short_Corr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Astrocyte group vs neuron group for responding GCaMP") +
  max.theme

#stats
# interactions and factors
GroupX_CompType=interaction(GCaMP_RCaMP.Nresp$GroupX, GCaMP_RCaMP.Nresp$CompType)
GroupX_Condition=interaction(GCaMP_RCaMP.Nresp$GroupX, GCaMP_RCaMP.Nresp$Condition)
GroupX_Nresp=interaction(GCaMP_RCaMP.Nresp$GroupX, GCaMP_RCaMP.Nresp$Nresponders)
Condition_Nresp=interaction(GCaMP_RCaMP.Nresp$Condition, GCaMP_RCaMP.Nresp$Nresponders)
GroupX_GroupY=interaction(GCaMP_RCaMP.Nresp$GroupX, GCaMP_RCaMP.Nresp$GroupY)
#+ (1|ROIs_trial) + (1|RCaMP_ROIs)
#correlation of only the stim window traces (short corr)
shortCorr.null = lmer(Short_Corr ~ (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
shortCorr.model1 = lmer(Short_Corr ~ GroupX + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
shortCorr.model1B = lmer(Short_Corr ~ GroupX_Condition + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
shortCorr.model2 = lmer(Short_Corr ~ CompType + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
shortCorr.model3 = lmer(Short_Corr ~ CompType+GroupX + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
shortCorr.model4 = lmer(Short_Corr ~ GroupX_CompType + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
shortCorr.model5 = lmer(Short_Corr ~ Nresponders+ (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
shortCorr.model5B = lmer(Short_Corr ~ Condition_Nresp+ (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
shortCorr.model6 = lmer(Short_Corr ~ GroupX_Nresp + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
shortCorr.model7 = lmer(Short_Corr ~ GroupX_GroupY + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)

shortCorr.anova <- anova(shortCorr.null,shortCorr.model1,shortCorr.model1B,shortCorr.model2,shortCorr.model3,
                         shortCorr.model4,shortCorr.model5,shortCorr.model5B,shortCorr.model6,shortCorr.model7)
print(shortCorr.anova)

# p values
shortCorr.GroupX <- glht(shortCorr.model1, mcp(GroupX= "Tukey"))
summary(shortCorr.GroupX)

shortCorr.GroupX_types <- glht(shortCorr.model4, mcp(GroupX_CompType= "Tukey"))
summary(shortCorr.GroupX_types)

shortCorr.CompType.pvalue <- glht(shortCorr.model2, mcp(CompType= "Tukey"))
summary(shortCorr.CompType.pvalue)

shortCorr.Nresp.pvalue <- glht(shortCorr.model5, mcp(Nresponders= "Tukey"))
summary(shortCorr.Nresp.pvalue)

shortCorr.NrespType.pvalue <- glht(shortCorr.model6, mcp(GroupX_Nresp= "Tukey"))
summary(shortCorr.NrespType.pvalue)

shortCorr.groupxgroupy.pvalue <- glht(shortCorr.model7, mcp(GroupX_GroupY= "Tukey"))
summary(shortCorr.groupxgroupy.pvalue)

######
df5A<- summarySE(GCaMP_RCaMP.Nresp, measurevar="xCorr", groupvars=c("GroupX"))
df5B1<- summarySE(GCaMP_RCaMP.Nresp, measurevar="xCorr", groupvars=c("CompType"))
df5B2<- summarySE(GCaMP_RCaMP.Nresp, measurevar="xCorr", groupvars=c("GroupX","CompType"))
df5C<- summarySE(GCaMP_RCaMP.Nresp, measurevar="xCorr", groupvars=c("Nresponders"))
df5D<- summarySE(GCaMP_RCaMP.Nresp, measurevar="xCorr", groupvars=c("GroupX","Nresponders"))
df5G<- summarySE(GCaMP_RCaMP.Nresp, measurevar="xCorr", groupvars=c("GroupX","GroupY"))
df5G2<- summarySE(GCaMP_RCaMP.Nresp, measurevar="xCorr", groupvars=c("GroupX","GroupY","CompType"))

ggplot(GCaMP_RCaMP.Nresp, aes(x=xCorr, fill=GroupX)) + geom_density(alpha=0.3)+
  ggtitle("density-") + max.theme


ggplot(data=df5A, aes(x=GroupX, y=xCorr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte groups") +
  ylab("xCorr") +
  max.theme

ggplot(data=df5B1 , aes(x=CompType, y=xCorr, fill=CompType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("xCorr") +
  ggtitle("xCorr for responding GCaMP") +
  max.theme

ggplot(data=df5B2 , aes(x=CompType, y=xCorr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("xCorr") +
  ggtitle("xCorr for responding GCaMP") +
  max.theme

ggplot(data=df5C, aes(x=Nresponders, y=xCorr, fill=Nresponders)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("xCorr") +
  ggtitle("Astrocyte group vs neuron group for responding GCaMP") +
  max.theme

ggplot(data=df5D, aes(x=Nresponders, y=xCorr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("xCorr") +
  ggtitle("Astrocyte group vs neuron group for responding GCaMP") +
  max.theme


ggplot(data=df5G, aes(x=GroupY, y=xCorr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("xCorr") +
  ggtitle("Astrocyte group vs neuron group for responding GCaMP") +
  max.theme

ggplot(data=df5G2, aes(x=interaction(GroupY,CompType), y=xCorr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("xCorr") +
  ggtitle("Astrocyte group vs neuron group for responding GCaMP") +
  max.theme


ggplot(data=df5G[df5G$GroupY=="fast N",], aes(x=GroupX, y=xCorr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("xCorr") +
  ggtitle("Astrocyte group vs neuron group for responding GCaMP") +
  max.theme
#stats
# Likelihood-ratio test 

#cross correlation
xCorr.null = lmer(xCorr ~ (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
xCorr.model2 = lmer(xCorr ~ GroupX + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
xCorr.model3 = lmer(xCorr ~ CompType + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
xCorr.model4 = lmer(xCorr ~ CompType+GroupX + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
xCorr.model5 = lmer(xCorr ~ GroupX_CompType + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
xCorr.model6 = lmer(xCorr ~ Nresponders + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
xCorr.model7 = lmer(xCorr ~ GroupX_Nresp + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
xCorr.model8 = lmer(xCorr ~ GroupX_GroupY + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)


xCorr.anova <- anova(xCorr.null,xCorr.model2,xCorr.model3,
                     xCorr.model4,xCorr.model5,xCorr.model6,
                     xCorr.model7,xCorr.model8)
print(xCorr.anova)

# p values

xCorr.CompType.pvalue <- glht(xCorr.model3, mcp(CompType= "Tukey"))
summary(xCorr.CompType.pvalue)

xCorr.GroupX <- glht(xCorr.model2, mcp(GroupX= "Tukey"))
summary(xCorr.GroupX)

xCorr.GroupX_types <- glht(xCorr.model5, mcp(GroupX_CompType= "Tukey"))
summary(xCorr.GroupX_types)

xCorr.Nresponders<- glht(xCorr.model6, mcp(Nresponders= "Tukey"))
summary(xCorr.Nresponders)

xCorr.Groupx_Nresp<- glht(xCorr.model7, mcp(GroupX_Nresp= "Tukey"))
summary(xCorr.Groupx_Nresp)

xCorr.Groupx_groupy<- glht(xCorr.model8, mcp(GroupX_GroupY= "Tukey"))
summary(xCorr.Groupx_groupy)

######
#Lag
df6A<- summarySE(GCaMP_RCaMP.Nresp, measurevar="Lag", groupvars=c("GroupX"))
df6B<- summarySE(GCaMP_RCaMP.Nresp, measurevar="Lag", groupvars=c("GroupX","CompType"))
df6C<- summarySE(GCaMP_RCaMP.Nresp, measurevar="Lag", groupvars=c("CompType"))

df6D<- summarySE(GCaMP_RCaMP.Nresp, measurevar="Lag", groupvars=c("Nresponders"))
df6E<- summarySE(GCaMP_RCaMP.Nresp, measurevar="Lag", groupvars=c("GroupX","Nresponders"))
df6F<- summarySE(GCaMP_RCaMP.Nresp, measurevar="Lag", groupvars=c("GroupX","GroupY"))


ggplot(data=df6A, aes(x=GroupX, y=Lag, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Groups") +
  ylab("Lag") +
  ggtitle("Lag for responding GCaMP") +
  max.theme

ggplot(data=df6B , aes(x=CompType, y=Lag, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Group") +
  ylab("Lag") +
  ggtitle("Lag for responding GCaMP") +
  max.theme


ggplot(data=df6C, aes(x=CompType, y=Lag, fill=CompType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("CompType") +
  ylab("Lag") +
  ggtitle("Lag for responding GCaMP") +
  max.theme

ggplot(data=df6D, aes(x=Nresponders, y=Lag, fill=Nresponders)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Nresponders") +
  ylab("Lag") +
  ggtitle("Lag for responding GCaMP") +
  max.theme

ggplot(data=df6E, aes(x=GroupX, y=Lag, fill=Nresponders)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Nresponders") +
  ylab("Lag") +
  ggtitle("Lag for responding GCaMP") +
  max.theme

ggplot(data=df6F, aes(x=GroupY, y=Lag, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Nresponders") +
  ylab("Lag") +
  ggtitle("Lag for responding GCaMP") +
  max.theme

ggplot(data=df6F[df6F$GroupY=="fast N",], aes(x=GroupX, y=Lag, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Nresponders") +
  ylab("Lag") +
  ggtitle("Lag for responding GCaMP") +
  max.theme

#Stats
#lag
Lag.null = lmer(Lag ~ (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
Lag.model1 = lmer(Lag ~ GroupX + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
Lag.model2 = lmer(Lag ~ CompType + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
Lag.model3 = lmer(Lag ~ CompType+GroupX + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
Lag.model4 = lmer(Lag ~ GroupX_CompType + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
Lag.model5 = lmer(Lag ~ Nresponders + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
Lag.model6 = lmer(Lag ~ GroupX+Nresponders + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
Lag.model7 = lmer(Lag ~ GroupX_Nresp + (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)
Lag.model8 = lmer(Lag ~ GroupX_GroupY+ (1|Animal) + (1|Spot), GCaMP_RCaMP.Nresp,REML=FALSE)

Lag.anova <- anova(Lag.null,Lag.model1,Lag.model2,Lag.model3,
                   Lag.model4,Lag.model5,Lag.model6,Lag.model7,Lag.model8)
print(Lag.anova)

# p values
Lag.GroupX <- glht(Lag.model1, mcp(GroupX= "Tukey"))
summary(Lag.GroupX)

lag.Nresp_pvalue <- glht(Lag.model5, mcp(Nresponders= "Tukey"))
summary(lag.Nresp_pvalue)

lag.CompType.pvalue <- glht(Lag.model2, mcp(CompType= "Tukey"))
summary(lag.CompType.pvalue)

Lag.GroupX_types <- glht(Lag.model4, mcp(GroupX_CompType= "Tukey"))
summary(Lag.GroupX_types)

Lag.GroupX_Nresp <- glht(Lag.model7, mcp(GroupX_Nresp= "Tukey"))
summary(Lag.GroupX_Nresp)

Lag.GroupX_GroupY <- glht(Lag.model8, mcp(GroupX_GroupY= "Tukey"))
summary(Lag.GroupX_GroupY)

#####
#distance

df7A<- summarySE(GCaMP_RCaMP.Nresp, measurevar="MinDistance", groupvars=c("GroupX"))
df7B<- summarySE(GCaMP_RCaMP.Nresp, measurevar="MinDistance", groupvars=c("Nresponders"))
df7C<- summarySE(GCaMP_RCaMP.Nresp, measurevar="MinDistance", groupvars=c("CompType"))
df7D<- summarySE(GCaMP_RCaMP.Nresp, measurevar="MinDistance", groupvars=c("GroupX","CompType"))
df7E<- summarySE(GCaMP_RCaMP.Nresp, measurevar="MinDistance", groupvars=c("GroupX","Nresponders"))
df7F<- summarySE(GCaMP_RCaMP.Nresp, measurevar="MinDistance", groupvars=c("GroupX","GroupY"))

ggplot(data=df7A, aes(x=GroupX, y=MinDistance, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=MinDistance-se, ymax=MinDistance+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Group") +
  ylab("MinDistance") +
  ggtitle("MinDistance for responding GCaMP") +
  max.theme

ggplot(data=df7D , aes(x=CompType, y=MinDistance, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=MinDistance-se, ymax=MinDistance+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Group") +
  ylab("MinDistance") +
  ggtitle("MinDistance for responding GCaMP") +
  max.theme

ggplot(data=df7E , aes(x=Nresponders, y=MinDistance, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=MinDistance-se, ymax=MinDistance+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Group") +
  ylab("MinDistance") +
  ggtitle("MinDistance for responding GCaMP") +
  max.theme

ggplot(data=df7B, aes(x=Nresponders, y=MinDistance, fill=Nresponders)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=MinDistance-se, ymax=MinDistance+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("ROI_Y_type") +
  ylab("MinDistance") +
  ggtitle("MinDistance for responding GCaMP") +
  max.theme

ggplot(data=df7C, aes(x=CompType, y=MinDistance, fill=CompType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=MinDistance-se, ymax=MinDistance+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("CompType") +
  ylab("MinDistance") +
  ggtitle("MinDistance for responding GCaMP") +
  max.theme

ggplot(data=df7F, aes(x=GroupY, y=MinDistance, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=MinDistance-se, ymax=MinDistance+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("CompType") +
  ylab("MinDistance") +
  ggtitle("MinDistance for responding GCaMP") +
  max.theme


#MinDistance
Dis.null = lmer(MinDistance ~ (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.Nresp,REML=FALSE)
Dis.model1 = lmer(MinDistance ~ GroupX + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.Nresp,REML=FALSE)
Dis.model2 = lmer(MinDistance ~ CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.Nresp,REML=FALSE)
Dis.model3 = lmer(MinDistance ~ CompType+GroupX + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.Nresp,REML=FALSE)
Dis.model4 = lmer(MinDistance ~ GroupX_CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.Nresp,REML=FALSE)
Dis.model5 = lmer(MinDistance ~ Nresponders + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.Nresp,REML=FALSE)
Dis.model6 = lmer(MinDistance ~ GroupX+Nresponders + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.Nresp,REML=FALSE)
Dis.model7 = lmer(MinDistance ~ GroupX_Nresp + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.Nresp,REML=FALSE)

Dis.anova <- anova(Dis.null,Dis.model1,Dis.model2,Dis.model3,Dis.model4,
                   Dis.model5,Dis.model6,Dis.model7)
print(Dis.anova)

# p values
Dis.Group <- glht(Dis.model1, mcp(GroupX= "Tukey"))
summary(Dis.Group)

Dis.comptypes <- glht(Dis.model2, mcp(CompType= "Tukey"))
summary(Dis.comptypes)

Dis.Group_types <- glht(Dis.model4, mcp(GroupX_CompType= "Tukey"))
summary(Dis.Group_types)

Dis.Nresp <- glht(Dis.model5, mcp(Nresponders= "Tukey"))
summary(Dis.Nresp)

Dis.group_Nresp <- glht(Dis.model7, mcp(GroupX_Nresp= "Tukey"))
summary(Dis.group_Nresp)


#####
ggplot(data=Corr_EarlyAUC, aes(x=Lag, fill=CompType)) +
  geom_density(alpha=.2) 
# box plot
ggplot(data=Corr_EarlyAUC, aes(x=CompType, y=Lag, fill=CompType)) +
  geom_boxplot() 

# cumulative fraction plot
ggplot(Corr_EarlyAUC, aes(x=Lag, colour = CompType)) + stat_ecdf()

######

