
#onset time comparisons- neurons vs. astrocytes

longstim.OT.comp <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/longstim_firstonset_comparisons.csv", header=TRUE, sep = ",")
shortstim.OT.comp <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/shortstim_firstonset_comparisons.csv", header=TRUE, sep = ",")
nostim.OT.comp <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/nostim_firstonset_comparisons.csv", header=TRUE, sep = ",")

longstim.OT.comp <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/longstim_firstonset_comparisons_fixedDis.csv", header=TRUE, sep = ",")
shortstim.OT.comp <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/shortstim_firstonset_comparisons_fixedDis.csv", header=TRUE, sep = ",")
nostim.OT.comp <- read.table("D:/Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/nostim_firstonset_comparisons_fixedDis.csv", header=TRUE, sep = ",")


longstim.OT.comp$compType<-paste(longstim.OT.comp$N_ROIType, longstim.OT.comp$A_ROIType, sep= "_")
shortstim.OT.comp$compType<-paste(shortstim.OT.comp$N_ROIType, shortstim.OT.comp$A_ROIType, sep= "_")
nostim.OT.comp$compType<-paste(nostim.OT.comp$N_ROIType, nostim.OT.comp$A_ROIType, sep= "_")


# get rid of onsets from the last few seconds of the trial (probably not a complete peak and we can't measure it)

#nostim.OT=nostim.OT[nostim.OT$N_Onset<43,]
#nostim.OT=nostim.OT[nostim.OT$A_Onset<43,]

# subset data to 3 sec on either side (so AC peaks 3 sec before or after neuronal peaks)
#nostim.OT.small<-subset(nostim.OT, TimeDiff<3 & TimeDiff>-3)
#nostim.OT.close<-subset(nostim.OT.small, distance<5)


######
# histograms of neuron-astrocyte time differences
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




