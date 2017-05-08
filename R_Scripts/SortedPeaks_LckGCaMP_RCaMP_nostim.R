
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
library('hexbin')

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
# data sets from nostim, shortstim (90Hz,1s), longstim (90Hz,8s)
# with hand selected somata and FLIKA activities for neurons AND astrocytes

########################
# load data

#nostim <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_nostim_05_04_2017.csv", header=TRUE, sep = ",")
#short <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_shortstim_05_04_2017.csv", header=TRUE, sep = ",")

nostim <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_nostim_28_04_2017.csv", header=TRUE, sep = ",")
short <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_shortstim_28_04_2017.csv", header=TRUE, sep = ",")
long <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/LckGC&RC_2D_longstim_28_04_2017.csv", header=TRUE, sep = ",")

lsm.options(pbkrtest.limit = 100000)

# onset time comparisons for nostim data
nostim.OT <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/nostim_onset_comparisons.csv", header=TRUE, sep = ",")


# exclude the neuropil ROIs, because they were hand selected and not necessary
nostim<-nostim[!(nostim$ROIname=="np"),]
nostim<-nostim[!(nostim$ROIname=="none"),]

short<-short[!(short$ROIname=="np"),]
long<-long[!(long$ROIname=="np"),]

# no stim peak data

nostim$ROIType= 0
nostimA<- subset(nostim, Channel=="GCaMP")
nostimB<- subset(nostim, Channel=="RCaMP")

# ROITypes
nostimA$ROIType[grepl("r",nostimA$ROIname)]="Process"
nostimA$ROIType[grepl("E",nostimA$ROIname)]="Endfoot"
nostimB$ROIType[grepl("r",nostimB$ROIname)]="Dendrite"
nostimB$ROIType[grepl("D",nostimB$ROIname)]="Dendrite"
nostimB$ROIType[grepl("N",nostimB$ROIname)]="Neuron"

nostim<-rbind(nostimA, nostimB)
nostim$ROIType<- as.factor(nostim$ROIType)

#unique ROI names
nostim$ROIs_trial<-paste(nostim$Animal, nostim$Spot, nostim$Trial,nostim$ROIname, sep= "_")

nostim$trials<-paste(nostim$Animal, nostim$Spot, nostim$Trial, sep= "_")


# remove matching astrocyte process and soma ROIs
Overlap= nostim$overlap!=0
nostim<-nostim[!Overlap,]
OverlapROIs<-unique(nostim$ROIs_trial[Overlap])


######
# No stim onset time comparisons
nostim.OT$compType<-paste(nostim.OT$N_ROIType, nostim.OT$A_ROIType, sep= "_")

# get rid of overlaping ROIs (likely the same region?)
nostim.OT<-subset(nostim.OT, !(N_ROI %in% OverlapROIs))
nostim.OT<-subset(nostim.OT, !(A_ROI %in% OverlapROIs))

# get rid of onsets from the last few seconds of the trial (probably not a complete peak and we can't measure it)

nostim.OT=nostim.OT[nostim.OT$N_Onset<43,]
nostim.OT=nostim.OT[nostim.OT$A_Onset<43,]

# subset data to 3 sec on either side (so AC peaks 3 sec before or after neuronal peaks)
nostim.OT.small<-subset(nostim.OT, TimeDiff<3 & TimeDiff>-3)

nostim.OT.close<-subset(nostim.OT.small, distance<5)


######
# histograms of time differences
ggplot(nostim.OT, aes(x=TimeDiff)) + geom_histogram(binwidth=0.0845, position="dodge") +
  ggtitle("onset time differences- all peaks") +
  max.theme

ggplot(nostim.OT, aes(x=TimeDiff, fill=A_ROIType)) + geom_histogram(binwidth=0.0845, position="dodge") +
  ggtitle("onset time differences- all peaks") +
  max.theme

# histograms of time differences
ggplot(nostim.OT.small, aes(x=TimeDiff)) + geom_histogram(binwidth=0.0845, position="dodge") +
  ggtitle("onset time differences- peaks close to zero time difference") + 
  max.theme

library('scales')
ggplot(nostim.OT.small, aes(x=TimeDiff, fill=A_ROIType)) + 
  geom_bar(aes(y = (..count..)/sum(..count..))) + 
  ## scale_y_continuous(labels = percent_format()) #version 3.0.9
  scale_y_continuous(labels = percent_format())+
  ggtitle("percent distribution- onset time differences")+ 
  max.theme  
  

# density of time differences
ggplot(nostim.OT.small, aes(x=TimeDiff)) + geom_density(aes(group=A_ROIType, colour=A_ROIType, fill=A_ROIType), alpha=0.3, adjust=1/3,size=1) +
  max.theme

ggplot(nostim.OT, aes(x=TimeDiff)) + geom_density(aes(group=A_ROIType, colour=A_ROIType, fill=A_ROIType), alpha=0.3, adjust=1/5,size=1) +
  max.theme

ggplot(nostim.OT.small, aes(x=TimeDiff)) + geom_density(aes(group=compType, colour=compType, fill=compType), alpha=0.3, adjust=1/5,size=1) +
  max.theme

ggplot(nostim.OT.small, aes(x=TimeDiff)) + geom_histogram(binwidth = 0.0845,aes( y=..density..,fill=A_ROIType)) +
  #geom_density(aes(group=A_ROIType), size=1) +
  ggtitle("onset time differences- peaks close to zero time difference") + 
  max.theme

ggplot(nostim.OT.small, aes(x=TimeDiff, fill=A_ROIType)) + 
  geom_histogram(aes(y = ..density..),binwidth=0.0845, alpha = 0.7) + 
  geom_density(aes(fill = A_ROIType), alpha = 0.5) 


ggplot(nostim.OT.small, aes(x=TimeDiff, y=..density..,colour=A_ROIType)) + stat_bin(geom="step", binwidth=(0.0845*2))+max.theme


#aes(y=..count../sum(..count..))

# histograms of ROI distances
ggplot(nostim.OT.small, aes(x=distance)) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("distance between ROIs- close peaks") + max.theme

ggplot(nostim.OT.small, aes(x=distance, fill=A_ROIType)) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("distance between ROIs- close peaks") + max.theme


#Distance between ROIs that are compared

ggplot(nostim.OT.small, aes(x=distance, y=TimeDiff, colour=compType)) +
  geom_point()+ max.theme

ggplot(nostim.OT.small, aes(x=distance, y=TimeDiff)) +
  geom_point(alpha=1/20)+ max.theme

ggplot(nostim.OT.small, aes(x=distance, y=TimeDiff)) +
  geom_count()+ max.theme

ggplot(nostim.OT.small, aes(x=distance, y=TimeDiff)) +
  geom_hex(bins=20)+ max.theme


# mean time diffs and mean distances
df1A1<-summarySE(nostim.OT.small, measurevar="TimeDiff", groupvars=c("compType"))
df1A2<-summarySE(nostim.OT.small, measurevar="TimeDiff", groupvars=c("A_ROIType"))

df2A1<-summarySE(nostim.OT.small, measurevar="distance", groupvars=c("compType"))
df2A2<-summarySE(nostim.OT.small, measurevar="distance", groupvars=c("A_ROIType"))

######
Before<- subset(nostim.OT.small, TimeDiff<=0)
After <- subset(nostim.OT.small, TimeDiff>=0)
Zero <- subset(nostim.OT.small, TimeDiff==0)


df1A3<-summarySE(Before, measurevar="TimeDiff", groupvars=c("A_ROIType"))
df1A4<-summarySE(After, measurevar="TimeDiff", groupvars=c("A_ROIType"))

df2A3<-summarySE(Before, measurevar="distance", groupvars=c("A_ROIType"))
df2A4<-summarySE(After, measurevar="distance", groupvars=c("A_ROIType"))

ggplot(After, aes(x=TimeDiff)) + geom_histogram(aes(group=A_ROIType, colour=A_ROIType, fill=A_ROIType), binwidth=0.0845)+
  max.theme

ggplot(After, aes(x=TimeDiff)) + geom_density(aes(group=A_ROIType, colour=A_ROIType, fill=A_ROIType), alpha=0.3, adjust=1/2,size=1) +
  max.theme

# stats for after group
# endfeet delayed compared to processes
OT.null = lmer(TimeDiff ~ (1|Animal) + (1|Spot) + (1|trials) + (1|N_ROI) + (1|A_ROI), After,REML=FALSE)
OT.model1 = lmer(TimeDiff ~ A_ROIType + (1|Animal) + (1|Spot) + (1|trials) +(1|N_ROI) + (1|A_ROI), After,REML=FALSE)
OT.anova <- anova(OT.null, OT.model1)
print(OT.anova)

OT.after.nostim<- glht(OT.model1, mcp(A_ROIType= "Tukey"))
summary(OT.after.nostim)



######

#percent astrocyte ROIs per field of view with at least one comparison

# find total number of ROIs per trial
# the find total number of ROIs after neurons, or before neurons
gcamp.nostim<- subset(nostim, Channel=="GCaMP")
nostim.OT.small$trials<- paste(nostim.OT.small$Animal, nostim.OT.small$Spot, nostim.OT.small$Trial, sep="_")
nostim.OT.small$ROIType<-nostim.OT.small$A_ROIType
# astrocyte ROIs per trial- total
ACROI_trial<- ddply(gcamp.nostim, c("Animal","Spot","trials"), summarise, AC_ROInum= length(unique(ROIs_trial)))

# astrocyte ROIs per trial- ROI type
ACROI_type_trial<- ddply(gcamp.nostim, c("Animal","Spot","trials","ROIType"), summarise, AC_ROInum= length(unique(ROIs_trial)))

# astrocyte ROIs with time differences around zero- total
ACROI_trial.OT.small<-ddply(nostim.OT.small, c("Animal","Spot","trials"), summarise, AC_ROInum= length(unique(A_ROI)))

# astrocyte ROIs with time differences around zero- ROI type
ACROI_type_trial.OT.small<- ddply(nostim.OT.small, c("Animal","Spot","trials","ROIType"), summarise, AC_ROInum= length(unique(A_ROI)))

# group data together for proportion calculation
numBefore<-ddply(Before, c("Animal","Spot","trials"), summarise, AC_ROInum= length(unique(A_ROI)))
numBefore.type<-ddply(Before, c("Animal","Spot","trials", "ROIType"), summarise, AC_ROInum= length(unique(A_ROI)))

numAfter<-ddply(After, c("Animal","Spot","trials"), summarise, AC_ROInum= length(unique(A_ROI)))
numAfter.type<-ddply(After, c("Animal","Spot","trials", "ROIType"), summarise, AC_ROInum= length(unique(A_ROI)))

numZero<-ddply(Zero, c("Animal","Spot","trials"), summarise, AC_ROInum= length(unique(A_ROI)))
numZero.type<-ddply(Zero, c("Animal","Spot","trials", "ROIType"), summarise, AC_ROInum= length(unique(A_ROI)))


#number of ROIs per trial 

Trialnames <-as.character(unique(ACROI_trial$trials))

#create dataframes that will be made
tot.ACnum.small <-data.frame()
bef.ACnum.small <-data.frame()
aft.ACnum.small <-data.frame()
zero.ACnum.small <-data.frame()

for (ii in 1:length(Trialnames))
{
  name =Trialnames[ii]
  subset1 = subset(ACROI_trial, trials == name) # total number of ROIs for this trial
  subset3 = subset(ACROI_trial.OT.small, trials == name) # number of ROIs with small time differences
  subset5 = subset(numBefore, trials == name) # number of ROIs with time before neurons- total
  subset7 = subset(numAfter, trials == name) # number of ROIs with time after neurons
  subset9 = subset(numZero, trials == name) # number of ROIs with zero

  # count total number of time diff ROIs per trial
  # use row data from trial if there are no comparisons
  if ((nrow(subset3) ==0)==TRUE)
  {subset3<- head(subset1,1)
  subset3$AC_ROInum = 0
  }
  subset3$numTimeDiffROIs<-0
  subset3$numTimeDiffROIs<-subset3$AC_ROInum/subset1$AC_ROInum
  tot.ACnum.small <-rbind(tot.ACnum.small, subset3)
  
  
  if ((nrow(subset5) ==0)==TRUE)
  {subset5<- head(subset1,1)
  subset5$AC_ROInum = 0
  }
  subset5$numTimeDiffROIs<-0
  subset5$numTimeDiffROIs<-subset5$AC_ROInum/subset1$AC_ROInum
  bef.ACnum.small  <-rbind(bef.ACnum.small, subset5)
  
  if ((nrow(subset7) ==0)==TRUE)
  {subset7<- head(subset1,1)
  subset7$AC_ROInum = 0
  }
  subset7$numTimeDiffROIs<-0
  subset7$numTimeDiffROIs<-subset7$AC_ROInum/subset1$AC_ROInum
  aft.ACnum.small  <-rbind(aft.ACnum.small, subset7)
  
  if ((nrow(subset9) ==0)==TRUE)
  {subset9<- head(subset1,1)
  subset9$AC_ROInum = 0
  }
  subset9$numTimeDiffROIs<-0
  subset9$numTimeDiffROIs<-subset9$AC_ROInum/subset1$AC_ROInum
  zero.ACnum.small  <-rbind(zero.ACnum.small, subset9)
  
}
 
df.numTotal<-summarySE(data=tot.ACnum.small, measurevar = "numTimeDiffROIs")
df.numBef <- summarySE(data=bef.ACnum.small, measurevar = "numTimeDiffROIs") 
df.numAft<- summarySE(data=aft.ACnum.small, measurevar = "numTimeDiffROIs") 
df.numZero <- summarySE(data=zero.ACnum.small, measurevar = "numTimeDiffROIs") 





#####
#ROI area


###### 
# ROI data
AfterROIs.AC<-unique(After$A_ROI)
AfterROIs.N<-unique(After$N_ROI)

AfterPeaks.AC<-subset(nostim2, ROIs_trial %in% AfterROIs.AC)# %in% ROIs_trial)
AfterPeaks.N<-subset(nostim2, ROIs_trial %in% AfterROIs.N)# %in% ROIs_trial)

AfterPeaks<-rbind(AfterPeaks.AC,AfterPeaks.N)
