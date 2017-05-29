
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
longstim.corr <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LongStim_Correlations_fixedDis.csv", header=TRUE, sep = ",")
shortstim.corr <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/shortstim_Correlations_fixedDis.csv", header=TRUE, sep = ",")

#CorrData <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LongStim_Correlations.csv", header=TRUE, sep = ",")
#CorrData <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LongStim_Correlations.csv", header=TRUE, sep = ",")

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




GCaMP_RCaMP<-subset(CorrData, CompChannel=="GCaMP_RCaMP")
GCaMP_RCaMP$ROIs_trial<-paste(GCaMP_RCaMP$Animal, GCaMP_RCaMP$Spot, GCaMP_RCaMP$Trial,GCaMP_RCaMP$ROI_X, sep= "_")

#remove random spots with gcamp neurons
ntypes=c("Neuron_Neuron","Neuron_Neuropil","Neuropil_Neuron","Neuropil_Neuropil")
GCaMP_RCaMP<-subset(GCaMP_RCaMP, !(CompType %in% ntypes))

#ggplot(GCaMP_RCaMP, aes(x=MinDistance, y=xCorr,fill=CompType)) + geom_point(aes(colour=CompType))
#ggplot(GCaMP_RCaMP, aes(x=MinDistance, y=Short_Corr,fill=CompType)) + geom_point(aes(colour=CompType))

df2A1<- summarySE(GCaMP_RCaMP.groups, measurevar="Short_Corr", groupvars=c("CompType"))
df2A2<- summarySE(GCaMP_RCaMP, measurevar="Short_Corr", groupvars=c("ROI_Y_type"))

df2B1<- summarySE(GCaMP_RCaMP, measurevar="xCorr", groupvars=c("CompType"))
df2B2<- summarySE(GCaMP_RCaMP, measurevar="xCorr", groupvars=c("ROI_Y_type"))

df2C1<- summarySE(GCaMP_RCaMP, measurevar="Lag", groupvars=c("CompType"))
df2C2<- summarySE(GCaMP_RCaMP, measurevar="Lag", groupvars=c("ROI_Y_type"))

ggplot(data=df2A1, aes(x=CompType, y=Short_Corr, fill=CompType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("CompType") +
  ylab("Short_Corr") +
  ggtitle("Short_Corr for gcamp rcamp in all trials") +
  scale_fill_manual(values=cbbPalette)+
max.theme

ggplot(data=df2A2, aes(x=ROI_Y_type, y=Short_Corr, fill=ROI_Y_type)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("NeuronType") +
  ylab("Short_Corr") +
  ggtitle("Short_Corr for gcamp rcamp in all trials") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(data=df2B1, aes(x=CompType, y=xCorr, fill=CompType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("CompType") +
  ylab("xCorr") +
  ggtitle("xCorr for gcamp rcamp in all trials") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(data=df2B2, aes(x=ROI_Y_type, y=xCorr, fill=ROI_Y_type)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("ROI_Y_type") +
  ylab("xCorr") +
  ggtitle("xCorr for gcamp rcamp in all trials") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(data=df2C1, aes(x=CompType, y=Lag, fill=CompType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Lag") +
  ggtitle("Lag for gcamp rcamp in all trials") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

ggplot(data=df2C2, aes(x=ROI_Y_type, y=Lag, fill=ROI_Y_type)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("ROI_Y_type") +
  ylab("Lag") +
  ggtitle("Lag for gcamp rcamp in all trials") +
  scale_fill_manual(values=cbbPalette)+
  max.theme

Group_CompType=interaction(GCaMP.CorrAUC$Group, GCaMP.CorrAUC$CompType)

#correlation of only the stim window traces (short corr)
# is there a difference between astrocytes and neuronal processes or dendrites
shortCorr.type.null = lmer(Short_Corr ~ (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP,REML=FALSE)
shortCorr.type.model1 = lmer(Short_Corr ~ ROI_Y_type + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP,REML=FALSE)
shortCorr.type.model2 = lmer(Short_Corr ~ CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP,REML=FALSE)

shortCorr.type.anova <- anova(shortCorr.type.null,shortCorr.type.model1,shortCorr.type.model2)
print(shortCorr.type.anova)

# p values
shortCorr.ROI_y_pvalue <- glht(shortCorr.type.model1, mcp(ROI_Y_type= "Tukey"))
summary(shortCorr.ROI_y_pvalue)

shortCorr.CompType.pvalue <- glht(shortCorr.type.model2, mcp(CompType= "Tukey"))
summary(shortCorr.CompType.pvalue)


#xcorr
xCorr.type.null = lmer(xCorr ~ (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP,REML=FALSE)
xCorr.type.model1 = lmer(xCorr ~ ROI_Y_type + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP,REML=FALSE)
xCorr.type.model2 = lmer(xCorr ~ CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP,REML=FALSE)

xCorr.type.anova <- anova(xCorr.type.null,xCorr.type.model1,xCorr.type.model2)
print(xCorr.type.anova)

# p values
xCorr.ROI_y_pvalue <- glht(xCorr.type.model1, mcp(ROI_Y_type= "Tukey"))
summary(xCorr.ROI_y_pvalue)

xCorr.CompType.pvalue <- glht(xCorr.type.model2, mcp(CompType= "Tukey"))
summary(xCorr.CompType.pvalue)


#Lag
lag.type.null = lmer(Lag ~ (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP,REML=FALSE)
lag.type.model1 = lmer(Lag ~ ROI_Y_type + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP,REML=FALSE)
lag.type.model2 = lmer(Lag ~ CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP,REML=FALSE)

lag.type.anova <- anova(lag.type.null,lag.type.model1,lag.type.model2)
print(lag.type.anova)

# p values
lag.ROI_y_pvalue <- glht(lag.type.model1, mcp(ROI_Y_type= "Tukey"))
summary(lag.ROI_y_pvalue)

lag.CompType.pvalue <- glht(lag.type.model2, mcp(CompType= "Tukey"))
summary(lag.CompType.pvalue)



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
  subset1=subset(GCaMP_RCaMP.group, RCaMP_ROIs==ROIy)
  if (nrow(subset1)>0)
  {
    subset1$GroupY=Groupy
    GCaMP_RCaMP.groups<-rbind(GCaMP_RCaMP.groups, subset1)
  }
}



# enter responding ROI information- high, mid or low responding neurons

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

for (ii in 1:nrow(NeuronSomas))
{
  Soma=NeuronSomas$ROInameUnique[ii]
  RespType=NeuronSomas$responders[ii]
  subset1=subset(GCaMP_RCaMP.groups, ROInameUnique==Soma)
  if (nrow(subset1)>0)
  {
    subset1$Nresponders=RespType
    GCaMP_RCaMP.Nresp<-rbind(GCaMP_RCaMP.Nresp, subset1)
  }
}


df4A<- summarySE(GCaMP_RCaMP.groups, measurevar="Short_Corr", groupvars=c("GroupX"))
df4B<- summarySE(GCaMP_RCaMP.groups, measurevar="Short_Corr", groupvars=c("GroupX","CompType"))
df4D<- summarySE(GCaMP_RCaMP.groups, measurevar="Short_Corr", groupvars=c("GroupX","GroupY"))
df4C<- summarySE(GCaMP_RCaMP.groups, measurevar="FisherShort", groupvars=c("GroupX"))

ggplot(data=df4A, aes(x=GroupX, y=Short_Corr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Short_Corr for responding GCaMP") +
  max.theme

ggplot(data=df4B , aes(x=CompType, y=Short_Corr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Short_Corr for responding GCaMP") +
  max.theme

ggplot(data=df4C, aes(x=GroupX, y=FisherShort, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=FisherShort-se, ymax=FisherShort+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Group") +
  ylab("Fisher") +
  ggtitle("Fisher for responding GCaMP") +
  max.theme

ggplot(data=df4D, aes(x=GroupX, y=Short_Corr, fill=GroupY)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Short_Corr-se, ymax=Short_Corr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("Short_Corr") +
  ggtitle("Astrocyte group vs neuron group for responding GCaMP") +
  max.theme


df5A<- summarySE(GCaMP_RCaMP.groups, measurevar="xCorr", groupvars=c("GroupX"))
df5B<- summarySE(GCaMP_RCaMP.groups, measurevar="xCorr", groupvars=c("GroupX","CompType"))
df5C<- summarySE(GCaMP_RCaMP.groups, measurevar="FisherxCorr", groupvars=c("GroupX"))
df5D<- summarySE(GCaMP_RCaMP.groups, measurevar="xCorr", groupvars=c("GroupX","GroupY"))


ggplot(data=df5A, aes(x=GroupX, y=xCorr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte groups") +
  ylab("xCorr") +
  ggtitle("xCorr for responding GCaMP") +
  max.theme

ggplot(data=df5B , aes(x=CompType, y=xCorr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("xCorr") +
  ggtitle("xCorr for responding GCaMP") +
  max.theme

ggplot(data=df5C, aes(x=GroupX, y=FisherxCorr, fill=GroupX)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=FisherxCorr-se, ymax=FisherxCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("xCorr") +
  ggtitle("xCorr for responding GCaMP") +
  max.theme

ggplot(data=df5D, aes(x=GroupX, y=xCorr, fill=GroupY)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=xCorr-se, ymax=xCorr+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Astrocyte Group") +
  ylab("xCorr") +
  ggtitle("xCorr for responding GCaMP") +
  max.theme


df6A<- summarySE(GCaMP_RCaMP.groups, measurevar="Lag", groupvars=c("GroupX"))
df6B<- summarySE(GCaMP_RCaMP.groups, measurevar="Lag", groupvars=c("GroupX","CompType"))
df6C<- summarySE(GCaMP_RCaMP.groups, measurevar="Lag", groupvars=c("GroupX","GroupY"))

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


ggplot(data=df6C, aes(x=GroupX, y=Lag, fill=GroupY)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag-se, ymax=Lag+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Ästorycte groups") +
  ylab("Lag") +
  ggtitle("Lag for responding GCaMP") +
  max.theme


#########
#stats
# Likelihood-ratio test 

# ROI by AUC
GroupX_CompType=interaction(GCaMP_RCaMP.groups$GroupX, GCaMP_RCaMP.groups$CompType)
GroupX_GroupY=interaction(GCaMP_RCaMP.groups$GroupX, GCaMP_RCaMP.groups$GroupY)

#correlation of only the stim window traces (short corr)
shortCorr.null = lmer(Short_Corr ~ (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
shortCorr.model1 = lmer(Short_Corr ~ GroupX + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
shortCorr.model2 = lmer(Short_Corr ~ CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
shortCorr.model3 = lmer(Short_Corr ~ CompType+GroupX + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
shortCorr.model4 = lmer(Short_Corr ~ GroupX_CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
shortCorr.model5 = lmer(Short_Corr ~ GroupX+GroupY + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
shortCorr.model6 = lmer(Short_Corr ~ GroupX_GroupY + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)

shortCorr.anova <- anova(shortCorr.null,shortCorr.model1,shortCorr.model2,shortCorr.model3,
                         shortCorr.model4,shortCorr.model5,shortCorr.model6)
print(shortCorr.anova)

# p values
shortCorr.GroupX <- glht(shortCorr.model1, mcp(GroupX= "Tukey"))
summary(shortCorr.GroupX)

shortCorr.GroupX_types <- glht(shortCorr.model4, mcp(GroupX_CompType= "Tukey"))
summary(shortCorr.GroupX_types)

shortCorr.GroupX_GroupY<- glht(shortCorr.model6, mcp(GroupX_GroupY= "Tukey"))
summary(shortCorr.GroupX_GroupY)


#cross correlation
xCorr.null = lmer(xCorr ~ (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
xCorr.model1 = lmer(xCorr ~ GroupX + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
xCorr.model2 = lmer(xCorr ~ CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
xCorr.model3 = lmer(xCorr ~ CompType+GroupX + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
xCorr.model4 = lmer(xCorr ~ GroupX_CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
xCorr.model5 = lmer(xCorr ~ GroupX + GroupY + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
xCorr.model6 = lmer(xCorr ~ GroupX_GroupY + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)

xCorr.anova <- anova(xCorr.null,xCorr.model1,xCorr.model2,xCorr.model3,
                     xCorr.model4,xCorr.model5,xCorr.model6)
print(xCorr.anova)

# p values
xCorr.GroupX <- glht(xCorr.model1, mcp(GroupX= "Tukey"))
summary(xCorr.GroupX)

xCorr.GroupX_types <- glht(xCorr.model4, mcp(GroupX_CompType= "Tukey"))
summary(xCorr.GroupX_types)

xCorr.GroupX_GroupY <- glht(xCorr.model6, mcp(GroupX_GroupY= "Tukey"))
summary(xCorr.GroupX_GroupY)

#lag
Lag.null = lmer(Lag ~ (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
Lag.model1 = lmer(Lag ~ GroupX + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
Lag.model2 = lmer(Lag ~ CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
Lag.model3 = lmer(Lag ~ CompType+GroupX + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
Lag.model4 = lmer(Lag ~ GroupX_CompType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
Lag.model5 = lmer(Lag ~ GroupX+GroupY + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)
Lag.model6 = lmer(Lag ~ GroupX_GroupY + (1|Animal) + (1|Spot)+ (1|ROIs_trial), GCaMP_RCaMP.groups,REML=FALSE)

Lag.anova <- anova(Lag.null,Lag.model1,Lag.model2,Lag.model3,
                   Lag.model4,Lag.model5,Lag.model6)
print(Lag.anova)

# p values
Lag.GroupX <- glht(Lag.model1, mcp(GroupX= "Tukey"))
summary(Lag.GroupX)

Lag.GroupX_types <- glht(Lag.model4, mcp(GroupX_CompType= "Tukey"))
summary(Lag.GroupX_types)

Lag.GroupX_GroupY <- glht(Lag.model6, mcp(GroupX_GroupY= "Tukey"))
summary(Lag.GroupX_GroupY)

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

