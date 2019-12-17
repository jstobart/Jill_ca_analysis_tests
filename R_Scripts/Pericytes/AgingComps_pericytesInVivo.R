library("lme4")
library("lmerTest")
library("lattice")
library("plyr")
library("dplyr")
library("ggplot2")
library("emmeans")
library("Rmisc")
#library("MASS")
library("multcomp")
library("reshape2")
library("tidyr")
#library("data.table")
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


MJ_baseline1<- read.table("J:/Megan Rodriguez/Analysis in Vivo/MJ/Calcium - Cell Scan/times/Results/MJ_peaks_baseline.csv", header=TRUE, sep = ",")
#pull out the spot name
pos= regexpr('SPOT', MJ_baseline1$Spot)
MJ_baseline1$SpotName<-substr(MJ_baseline1 $Spot,pos, pos+5)

MJ_baseline2<- read.table("J:/Yewande Anozie/Fall 2019/MJ/MJ Calcium Results/MJ_peaks_Nov7.csv", header=TRUE, sep = ",")
#pull out the spot name
pos= regexpr('Spot', MJ_baseline2$Spot)
MJ_baseline2$SpotName<-substr(MJ_baseline2$Spot,pos, pos+5)

Spot1<-subset(MJ_baseline1, SpotName=="SPOT1")
spot1names<-unique(Spot1$Spot)
spot1names1<-spot1names[1]
Spot1$time<-"6 months"
Spot1$time[Spot1$Spot==spot1names1]<-"5 months"
Spot1_10<-subset(MJ_baseline2, SpotName=="Spot1")
Spot1_10$SpotName="SPOT1"
Spot1_10$time<-"10 months"
Spot1<-rbind(Spot1, Spot1_10)

Spot3<-subset(MJ_baseline1, SpotName=="SPOT3")
spot3names<-unique(Spot3$Spot)
spot3names1<-spot3names[1]
Spot3$time<-"6 months"
Spot3$time[Spot3$Spot==spot3names1]<-"5 months"
Spot3_10<-subset(MJ_baseline2, SpotName=="Spot3")
Spot3_10$SpotName="SPOT3"
Spot3_10$time<-"10 months"
Spot3<-rbind(Spot3, Spot3_10)

Spot6<-subset(MJ_baseline1, SpotName=="SPOT6")
spot6names<-unique(Spot6$Spot)
spot6names1<-spot6names[1]
Spot6$time<-"6 months"
Spot6$time[Spot6$Spot==spot6names1]<-"5 months"
Spot6_10<-subset(MJ_baseline2, SpotName=="Spot6")
Spot6_10$SpotName="SPOT6"
Spot6_10$time<-"10 months"
Spot6<-rbind(Spot6, Spot6_10)

Spot4<-subset(MJ_baseline1, SpotName=="SPOT4")
spot4names<-unique(Spot4$Spot)
spot4names1<-spot4names[1]
Spot4$time<-"6 months"
Spot4$time[Spot4$Spot==spot4names1]<-"5 months"
Spot4_10<-subset(MJ_baseline2, SpotName=="Spot4")
Spot4_10$SpotName="SPOT4"
Spot4_10$time<-"10 months"
Spot4<-rbind(Spot4, Spot4_10)

Spot5<-subset(MJ_baseline1, SpotName=="SPOT5")
Spot5$time<-"6 months"
Spot5_10<-subset(MJ_baseline2, SpotName=="Spot5")
Spot5_10$SpotName="SPOT5"
Spot5_10$time<-"10 months"
Spot5<-rbind(Spot5, Spot5_10)

calcium<-rbind(Spot1, Spot3, Spot4, Spot5, Spot6)

calcium$SpotROI<-paste(calcium$SpotName, calcium$roiName)

# set factors
calcium$time <-as.factor(calcium$time) 
calcium$time<- factor(calcium$time, levels = c("5 months","6 months", "10 months"))

# remove 5 months just for compairson sake

calcium<-subset(calcium, time!="5 months")

# average per ROI.... but the ROI names are not the same!

calcium.ROI<-ddply(calcium, c("SpotROI", "time"), summarise, 
                      meanAmp=mean(amplitude, na.rm = TRUE), meanDur = mean(halfWidth, na.rm = TRUE), meanAUC= mean(peakAUC, na.rm = TRUE), numEvents = length(peakAUC))

# average per field of view

calcium.Spot<-ddply(calcium, c("SpotName", "time"), summarise, 
                   meanAmp=mean(amplitude, na.rm = TRUE), meanDur = mean(halfWidth, na.rm = TRUE), meanAUC= mean(peakAUC, na.rm = TRUE), numEvents = length(peakAUC))

calcium.Spot$meanAmp[is.na(calcium.Spot$meanAmp)]=0
calcium.Spot$meanDur[is.na(calcium.Spot$meanDur)]=0
calcium.Spot$meanAUC[is.na(calcium.Spot$meanAUC)]=0



#duration
df1A <- summarySE(calcium, measurevar="halfWidth", groupvars=c("time"), na.rm = TRUE)
df1B <- summarySE(calcium.ROI, measurevar="meanDur", groupvars=c("time"),na.rm = TRUE)
df1C <- summarySE(calcium.Spot, measurevar="meanDur", groupvars=c("time"))

ggplot(data=df1A, aes(x=time, y=halfWidth, fill=time)) +
  geom_errorbar(aes(ymin=halfWidth-se, ymax=halfWidth+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("half width (s)") +
  scale_fill_manual(
    values=c("black", "red", "blue")) +
  max.theme

ggplot(data=df1B, aes(x=time, y=meanDur, fill=time)) +
  geom_errorbar(aes(ymin=meanDur-se, ymax=meanDur+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("half width (s)") +
  scale_fill_manual(
    values=c("black", "red","blue")) +
  max.theme 

ggplot(data=df1C, aes(x=time, y=meanDur, fill=time)) +
  geom_errorbar(aes(ymin=meanDur-se, ymax=meanDur+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("mean half width (s) per spot") +
  scale_fill_manual(
    values=c("black", "red","blue")) +
  max.theme 

#amplitude
df2A <- summarySE(calcium, measurevar="amplitude", groupvars=c("time"), na.rm = TRUE)
df2B <- summarySE(calcium.ROI, measurevar="meanAmp", groupvars=c("time"),na.rm = TRUE)
df2C <- summarySE(calcium.Spot, measurevar="meanAmp", groupvars=c("time"))


ggplot(data=df2A, aes(x=time, y=amplitude, fill=time)) +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude") +
  scale_fill_manual(
    values=c("black", "red", "blue")) +
  max.theme

ggplot(data=df2B, aes(x=time, y=meanAmp, fill=time)) +
  geom_errorbar(aes(ymin=meanAmp-se, ymax=meanAmp+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amp") +
  scale_fill_manual(
    values=c("black", "red", "blue")) +
  max.theme 

ggplot(data=df2C, aes(x=time, y=meanAmp, fill=time)) +
  geom_errorbar(aes(ymin=meanAmp-se, ymax=meanAmp+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("mean amp per Spot") +
  scale_fill_manual(
    values=c("black", "red", "blue")) +
  max.theme 

#peak AUC
df3A <- summarySE(calcium, measurevar="peakAUC", groupvars=c("time"),na.rm = TRUE)
df3B <- summarySE(calcium.ROI, measurevar="meanAUC", groupvars=c("time"),na.rm = TRUE)
df3C <- summarySE(calcium.Spot, measurevar="meanAUC", groupvars=c("time"),na.rm = TRUE)

ggplot(data=df3A, aes(x=time, y=peakAUC, fill=time)) +
  geom_errorbar(aes(ymin=peakAUC-se, ymax=peakAUC+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("peakAUC") +
  scale_fill_manual(
    values=c("black", "red")) +
  max.theme

ggplot(data=df3B, aes(x=time, y=meanAUC, fill=time)) +
  geom_errorbar(aes(ymin=meanAUC-se, ymax=meanAUC+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("peakAUC") +
  scale_fill_manual(
    values=c("black", "red")) +
  max.theme

ggplot(data=df3C, aes(x=time, y=meanAUC, fill=time)) +
  geom_errorbar(aes(ymin=meanAUC-se, ymax=meanAUC+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("peakAUC") +
  scale_fill_manual(
    values=c("black", "red")) +
  max.theme

# number of events
df4B <- summarySE(calcium.ROI, measurevar="numEvents", groupvars=c("time"))
df4C <- summarySE(calcium.Spot, measurevar="numEvents", groupvars=c("time"))

ggplot(data=df4B, aes(x=time, y=numEvents, fill=time)) +
  geom_errorbar(aes(ymin=numEvents-se, ymax=numEvents+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("numEvents per ROI") +
  scale_fill_manual(
    values=c("black", "red")) +
  max.theme

ggplot(data=df4C, aes(x=time, y=numEvents, fill=time)) +
  geom_errorbar(aes(ymin=numEvents-se, ymax=numEvents+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("numEvents per Spot") +
  scale_fill_manual(
    values=c("black", "red")) +
  max.theme

#####
# stats
dur.null = lmer(meanDur ~ (1|SpotName), calcium.Spot,REML=FALSE)
dur.model1 = lmer(meanDur~ time + (1|SpotName), calcium.Spot,REML=FALSE)
dur.anova <- anova(dur.null, dur.model1)
print(dur.anova)

amp.null = lmer(meanAmp~ (1|SpotName), calcium.Spot,REML=FALSE)
amp.model1 = lmer(meanAmp~ time + (1|SpotName), calcium.Spot, REML=FALSE)
amp.anova <- anova(amp.null, amp.model1)
print(amp.anova)

AUC.null = lmer(meanAUC~ (1|SpotName), calcium.Spot,REML=FALSE)
AUC.model1 = lmer(meanAUC~ time + (1|SpotName), calcium.Spot, REML=FALSE)
AUC.anova <- anova(AUC.null, AUC.model1)
print(AUC.anova)

#num events
NE.null = lmer(numEvents ~ (1|SpotName), calcium.Spot,REML=FALSE)
NE.model1 = lmer(numEvents ~ time + (1|SpotName), calcium.Spot,REML=FALSE)
NE.anova <- anova(NE.null, NE.model1)
print(NE.anova)


######

# diameter

MJ_baseline1 <- read.table("J:/Megan Rodriguez/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_19/2019_08_19_SPOT1_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline2 <- read.table("J:/Megan Rodriguez/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_19/2019_08_19_SPOT1_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline3 <- read.table("J:/Megan Rodriguez/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_19/2019_08_19_SPOT1_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline4 <- read.table("J:/Megan Rodriguez/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_19/2019_08_19_SPOT6_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline5 <- read.table("J:/Megan Rodriguez//Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_19/2019_08_19_SPOT3_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline6 <- read.table("J:/Megan Rodriguez//Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_22/2019_07_22_SPOT4_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline7 <- read.table("J:/Megan Rodriguez//Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_22/2019_07_22_SPOT4_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline8 <- read.table("J:/Megan Rodriguez//Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_22/2019_07_22_SPOT5_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")

MJ_baseline9 <- read.table("J:/Yewande Anozie/Fall 2019/MJ/MJ Blood Flow Results/2019_11_07_SPOT1_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline10 <- read.table("J:/Yewande Anozie/Fall 2019/MJ/MJ Blood Flow Results/2019_11_07_SPOT1_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline11 <- read.table("J:/Yewande Anozie/Fall 2019/MJ/MJ Blood Flow Results/2019_11_07_SPOT1_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline12 <- read.table("J:/Yewande Anozie/Fall 2019/MJ/MJ Blood Flow Results/2019_11_07_SPOT6_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline13 <- read.table("J:/Yewande Anozie/Fall 2019/MJ/MJ Blood Flow Results/2019_11_07_SPOT3_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline14 <- read.table("J:/Yewande Anozie/Fall 2019/MJ/MJ Blood Flow Results/2019_11_07_SPOT4_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline15 <- read.table("J:/Yewande Anozie/Fall 2019/MJ/MJ Blood Flow Results/2019_11_07_SPOT4_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline16 <- read.table("J:/Yewande Anozie/Fall 2019/MJ/MJ Blood Flow Results/2019_11_07_SPOT5_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")

MJ_BL1<-summarySE(MJ_baseline1, measurevar="diameter", na.rm = TRUE)
MJ_BL1$Animal<-"MJ"
MJ_BL1$Vessel<-"SPOT1_v1"
MJ_BL1$time<-"6 months"

MJ_BL2<-summarySE(MJ_baseline2, measurevar="diameter", na.rm = TRUE)
MJ_BL2$Animal<-"MJ"
MJ_BL2$Vessel<-"SPOT1_v2"
MJ_BL2$time<-"6 months"

MJ_BL3<-summarySE(MJ_baseline3, measurevar="diameter", na.rm = TRUE)
MJ_BL3$Animal<-"MJ"
MJ_BL3$Vessel<-"SPOT1_v3"
MJ_BL3$time<-"6 months"

MJ_BL4<-summarySE(MJ_baseline4, measurevar="diameter", na.rm = TRUE)
MJ_BL4$Animal<-"MJ"
MJ_BL4$Vessel<-"SPOT6_v1"
MJ_BL4$time<-"6 months"

MJ_BL5<-summarySE(MJ_baseline5, measurevar="diameter", na.rm = TRUE)
MJ_BL5$Animal<-"MJ"
MJ_BL5$Vessel<-"SPOT3_v1"
MJ_BL5$time<-"6 months"

MJ_BL6<-summarySE(MJ_baseline6, measurevar="diameter", na.rm = TRUE)
MJ_BL6$Animal<-"MJ"
MJ_BL6$Vessel<-"SPOT4_v1"
MJ_BL6$time<-"6 months"

MJ_BL7<-summarySE(MJ_baseline7, measurevar="diameter", na.rm = TRUE)
MJ_BL7$Animal<-"MJ"
MJ_BL7$Vessel<-"SPOT4_v2"
MJ_BL7$time<-"6 months"

MJ_BL8<-summarySE(MJ_baseline8, measurevar="diameter", na.rm = TRUE)
MJ_BL8$Animal<-"MJ"
MJ_BL8$Vessel<-"SPOT5_v3"
MJ_BL8$time<-"6 months"


MJ_BL9<-summarySE(MJ_baseline9, measurevar="diameter", na.rm = TRUE)
MJ_BL9$Animal<-"MJ"
MJ_BL9$Vessel<-"SPOT1_v1"
MJ_BL9$time<-"10 months"

MJ_BL10<-summarySE(MJ_baseline10, measurevar="diameter", na.rm = TRUE)
MJ_BL10$Animal<-"MJ"
MJ_BL10$Vessel<-"SPOT1_v2"
MJ_BL10$time<-"10 months"

MJ_BL11<-summarySE(MJ_baseline11, measurevar="diameter", na.rm = TRUE)
MJ_BL11$Animal<-"MJ"
MJ_BL11$Vessel<-"SPOT1_v3"
MJ_BL11$time<-"10 months"

MJ_BL12<-summarySE(MJ_baseline12, measurevar="diameter", na.rm = TRUE)
MJ_BL12$Animal<-"MJ"
MJ_BL12$Vessel<-"SPOT6_v1"
MJ_BL12$time<-"10 months"

MJ_BL13<-summarySE(MJ_baseline13, measurevar="diameter", na.rm = TRUE)
MJ_BL13$Animal<-"MJ"
MJ_BL13$Vessel<-"SPOT3_v1"
MJ_BL13$time<-"10 months"

MJ_BL14<-summarySE(MJ_baseline14, measurevar="diameter", na.rm = TRUE)
MJ_BL14$Animal<-"MJ"
MJ_BL14$Vessel<-"SPOT4_v1"
MJ_BL14$time<-"10 months"

MJ_BL15<-summarySE(MJ_baseline15, measurevar="diameter", na.rm = TRUE)
MJ_BL15$Animal<-"MJ"
MJ_BL15$Vessel<-"SPOT4_v2"
MJ_BL15$time<-"10 months"

MJ_BL16<-summarySE(MJ_baseline16, measurevar="diameter", na.rm = TRUE)
MJ_BL16$Animal<-"MJ"
MJ_BL16$Vessel<-"SPOT5_v3"
MJ_BL16$time<-"10 months"

diameter<-rbind(MJ_BL1,MJ_BL2,MJ_BL3,MJ_BL4,MJ_BL5,MJ_BL6, MJ_BL7, MJ_BL8, MJ_BL9, 
                MJ_BL10, MJ_BL11, MJ_BL12, MJ_BL13, MJ_BL14, MJ_BL15, MJ_BL16)

diameter$time<-as.factor(diameter$time)
diameter$time<-factor(diameter$time, levels = c("6 months", "10 months"))

df5 <- summarySE(diameter, measurevar="diameter", groupvars=c("time"))

ggplot(data=df5, aes(x=time, y=diameter, fill=time)) +
  geom_errorbar(aes(ymin=diameter-se, ymax=diameter+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("diameter") +
  scale_fill_manual(
    values=c("black", "red")) +
  max.theme

#diameter stats
diam.null = lmer(diameter ~ (1|Vessel), diameter,REML=FALSE)
diam.model1 = lmer(diameter ~ time + (1|Vessel), diameter,REML=FALSE)
diam.anova <- anova(diam.null, diam.model1)
print(diam.anova)


# diameter oscillations

MJ_baseline_oscill <- read.table("J:/Megan Rodriguez/Analysis In Vivo/MJ/Line Scan - Diameter/Oscillations/Results/DiameterOscillationsAnalysis_baseline.csv", header=TRUE, sep = ",")

T1<-MJ_baseline_oscill[5:9 ,]
T1$time<-"T1"

T2<-MJ_baseline_oscill[18:22 ,]
T2$time<-"T2"

MJ_oscill<-rbind(T1, T2)


df5 <- summarySE(MJ_oscill, measurevar="peaksPerMin", groupvars=c("time"))

ggplot(data=df5, aes(x=time, y=peaksPerMin, fill=time)) +
  geom_errorbar(aes(ymin=peaksPerMin-se, ymax=peaksPerMin+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("peaksPerMin") +
  scale_fill_manual(
    values=c("black", "red")) +
  max.theme

#diameter stats
diam.null = lmer(diameter ~ (1|Vessel), diameter,REML=FALSE)
diam.model1 = lmer(diameter ~ time + (1|Vessel), diameter,REML=FALSE)
diam.anova <- anova(diam.null, diam.model1)
print(diam.anova)
