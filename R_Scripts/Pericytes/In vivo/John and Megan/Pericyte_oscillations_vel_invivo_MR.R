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

########################

# load files
#baseline
RR_baseline1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/Velocity/Results/VelocityAnalysis_baseline.csv", header=TRUE, sep = ",")

MJ_baseline1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/Velocity/Results/VelocityAnalysis_baseline.csv", header=TRUE, sep = ",")

#nimodipine
RR_nimodipine1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/Velocity/Results/VelocityAnalysis_nimodipine.csv", header=TRUE, sep = ",")

MJ_nimodipine1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/Velocity/Results/VelocityAnalysis_nimodipine.csv", header=TRUE, sep = ",")

#pyr3
RR_pyr3_1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/Velocity/Results/VelocityAnalysis_pyr3.csv", header=TRUE, sep = ",")

MJ_pyr3_1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/Velocity/Results/VelocityAnalysis_pyr3.csv", header=TRUE, sep = ",")

#combine data together for each treatment

AllData<-rbind(RR_baseline1,RR_nimodipine1,RR_pyr3_1,MJ_baseline1,MJ_nimodipine1,MJ_pyr3_1)
AllData$Drug<-as.factor(AllData$Drug)

#make a unqiue ROI name
AllData$SpotName <-paste(AllData$Animal, AllData$Vessel,sep= "_")

AllData$peaksPerMin[AllData$peaksPerMin=="NaN"]=0
AllData$peakProminence[AllData$peakProminence=="NaN"]=0
AllData$peakWidth[AllData$peakWidth=="NaN"]=0
AllData$CycleTime[AllData$CycleTime=="NaN"]=0

###############################
#histograms
ggplot(AllData, aes(x=peaksPerMin, fill=Drug))+ geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of Number of Peaks")

ggplot(AllData, aes(x=peakProminence, fill=Drug))+ geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of Peak Prominence")

ggplot(AllData, aes(x=peakWidth, fill=Drug)) + geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("Distribution of Peak Width") 

ggplot(AllData, aes(x=CycleTime, fill=Drug)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of Cycle Time") 

################
#peaksPerMin
df1A <- summarySE(AllData, measurevar="peaksPerMin", groupvars=c("Drug"))
df1B <- summarySE(AllData, measurevar="peaksPerMin", groupvars=c("Drug","Animal"))

ggplot(data=df1A, aes(x=Drug, y=peaksPerMin, fill=Drug)) +
  geom_errorbar(aes(ymin=peaksPerMin-se, ymax=peaksPerMin+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("peaks/min") +
  scale_fill_manual(
    values=c("black", "red", "blue")) +
  max.theme

ggplot(data=df1B, aes(x=Animal, y=peaksPerMin, fill=Drug)) +
  geom_errorbar(aes(ymin=peaksPerMin-se, ymax=peaksPerMin+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min") +
  scale_fill_manual(
    values=c("black", "red", "blue")) +
  max.theme

#duration
df2A <- summarySE(AllData, measurevar="peakWidth", groupvars=c("Drug"))

ggplot(data=df2A, aes(x=Drug, y=peakWidth, fill=Drug)) +
  geom_errorbar(aes(ymin=peakWidth-se, ymax=peakWidth+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("peak duration (s)") +
  scale_fill_manual(
    values=c("black", "red", "blue")) +
  max.theme

#amplitude
df3A <- summarySE(AllData, measurevar="peakProminence", groupvars=c("Drug"))

ggplot(data=df3A, aes(x=Drug, y=peakProminence, fill=Drug)) +
  geom_errorbar(aes(ymin=peakProminence-se, ymax=peakProminence+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude (dF/F)") +
  scale_fill_manual(
    values=c("black", "red", "blue")) +
  max.theme

######
#stats

#peaksPerMinuency
peaksPerMin.null = lmer(peaksPerMin ~ (1|SpotName), AllData,REML=FALSE)
peaksPerMin.model1 = lmer(peaksPerMin~ Drug + (1|SpotName), AllData,REML=FALSE)
peaksPerMin.anova <- anova(peaksPerMin.null, peaksPerMin.model1)
print(peaksPerMin.anova)

#amplitude
peakProminence.null = lmer(peakProminence ~ (1|SpotName), AllData,REML=FALSE)
peakProminence.model1 = lmer(peakProminence~ Drug + (1|SpotName), AllData,REML=FALSE)
peakProminence.anova <- anova(peakProminence.null, peakProminence.model1)
print(peakProminence.anova)

# look at residuals
#plot(peakProminence.model1)

# p values
peakProminence.drug <- lsmeans(peakProminence.model1, pairwise ~ Drug, glhargs=list())
summary(peakProminence.drug)

