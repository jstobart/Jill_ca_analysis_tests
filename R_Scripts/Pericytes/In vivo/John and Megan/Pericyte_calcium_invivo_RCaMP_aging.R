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

#MJ_baseline1 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Calcium Results/RCaMP 59721 R MJ/Conditions/Results/MJ_peaks_baseline.csv", header=TRUE, sep = ",")

MJ_9month <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Calcium Results/RCaMP 59721 R MJ/2019_11_07/RCaMP_59721_peaks_9months_02_2020.csv", header=TRUE, sep = ",")

MJ_11month <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Calcium Results/RCaMP 59721 R MJ/2020_01_06/RCaMP_59721_peaks_11months_02_2020.csv", header=TRUE, sep = ",")

MJ_7month <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Calcium Results/RCaMP 59721 R MJ/2019_08_19/Results/RCaMP_59721_peaks_7months_02_2020.csv", header=TRUE, sep = ",")


MJ_baseline1$Condition="baseline"
MJ_9month$Condition="month9"
MJ_7month$Condition="month7"
MJ_11month$Condition=NULL
MJ_11month$Condition="month11"
MJ_11month$Animal="MJ"

# adjust all spot names to be the same
MJ_baseline1$SpotName<-paste("spot",str_sub(MJ_baseline1$Spot, -1, -1), sep="")
MJ_7month$SpotName<-paste("spot",str_sub(MJ_7month$Spot, -1, -1), sep="")
MJ_9month$SpotName<-paste("spot",str_sub(MJ_9month$Spot, -1, -1), sep="")


#combine all data
AllData<-rbind(MJ_7month,MJ_9month, MJ_11month)
AllData$Condition<-as.factor(AllData$Condition)

AllData$Condition<-factor(AllData$Condition, levels=c("month7", "month9", "month11"))

#make a unqiue ROI name
AllData$ROI_Spot <-paste(AllData$Animal, AllData$SpotName, AllData$roiName, sep= "_")
AllData$Spot_trial <-paste(AllData$SpotName, AllData$Trial, sep= "_")

AllData$Duration<-AllData$halfWidth*2


# separate out Condition data from spots that only have baseline
#ConditionROIs<-unique(AllData$ROI_Spot[AllData$Condition=="pyr3"])
#ConditionData<-subset(AllData, ROI_Spot %in% ConditionROIs)

#Baseline Data
#BL_Data<-subset(AllData, Condition=="baseline")

###############################
#histograms
ggplot(AllData, aes(x=Duration, fill=Condition))+ geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of ROI duration")

ggplot(AllData, aes(x=halfWidth, fill=Condition))+ geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of ROI duration")

ggplot(AllData, aes(x=amplitude, fill=Condition)) + geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("Distribution of ROI amplitudes") 

ggplot(AllData, aes(x=peakAUC, fill=Condition)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of peak AUC") 

#########


# pool data for each cell (mean amplitude, total number of signals, etc.)
ROI.means<- ddply(AllData, c("Animal", "SpotName","Condition","ROI_Spot"), 
                  summarise, AUC_mean = mean(peakAUC, na.rm=TRUE), dur_mean = mean(halfWidth,na.rm=TRUE), 
                  amp_mean = mean(amplitude,na.rm=TRUE), nEvents = length(amplitude), nTrials = length(unique(Spot_trial)))

ROI.means$dur_mean[ROI.means$AUC_mean=="NaN"]=0
ROI.means$amp_mean[ROI.means$AUC_mean=="NaN"]=0
ROI.means$nEvents[ROI.means$AUC_mean=="NaN"]=0
ROI.means$AUC_mean[ROI.means$AUC_mean=="NaN"]=0

ROI.means$freq<-ROI.means$nEvents/(ROI.means$nTrials*1) # number of signals/min

Spot.means<- ddply(ROI.means, c("Animal", "SpotName","Condition"), 
                  summarise, AUC_mean2 = mean(AUC_mean, na.rm=TRUE), dur_mean2 = mean(dur_mean,na.rm=TRUE), 
                  amp_mean2 = mean(amp_mean,na.rm=TRUE), freq_mean= mean(freq, na.rm=TRUE))

#only spots that are in all conditions
Conditions<-subset(Spot.means, SpotName %in% MJ_11month$SpotName)

################
#frequency
df1A <- summarySE(ROI.means, measurevar="freq", groupvars=c("Condition"))
df1B <- summarySE(Spot.means, measurevar="freq_mean", groupvars=c("Condition"))
df1C <- summarySE(Conditions, measurevar="freq_mean", groupvars=c("Condition"))

ggplot(data=df1A, aes(x=Condition, y=freq, fill=Condition)) +
  geom_errorbar(aes(ymin=freq-se, ymax=freq+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red", "blue","green")) +
  max.theme

ggplot(data=df1B, aes(x=Condition, y=freq_mean, fill=Condition)) +
  geom_errorbar(aes(ymin=freq_mean-se, ymax=freq_mean+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red", "blue","green")) +
  max.theme

ggplot(data=df1C, aes(x=Condition, y=freq_mean, fill=Condition)) +
  geom_errorbar(aes(ymin=freq_mean-se, ymax=freq_mean+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red", "blue","green")) +
  max.theme

#duration
df2A <- summarySE(ROI.means, measurevar="dur_mean", groupvars=c("Condition"))

ggplot(data=df2A, aes(x=Condition, y=dur_mean, fill=Condition)) +
  geom_errorbar(aes(ymin=dur_mean-se, ymax=dur_mean+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("half width (s)") +
  scale_fill_manual(
    values=c("black", "red", "blue","green")) +
  max.theme

#amplitude
df3A <- summarySE(ROI.means, measurevar="amp_mean", groupvars=c("Condition"))
df3C <- summarySE(Conditions, measurevar="amp_mean2", groupvars=c("Condition"))

ggplot(data=df3A, aes(x=Condition, y=amp_mean, fill=Condition)) +
  geom_errorbar(aes(ymin=amp_mean-se, ymax=amp_mean+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude (dF/F)") +
  scale_fill_manual(
    values=c("black", "red", "blue", "green")) +
  max.theme

ggplot(data=df3C, aes(x=Condition, y=amp_mean2, fill=Condition)) +
  geom_errorbar(aes(ymin=amp_mean2-se, ymax=amp_mean2+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude (dF/F)") +
  scale_fill_manual(
    values=c("black", "red", "blue", "green")) +
  max.theme


######
#stats

#frequency
freq.null = lmer(freq ~ (1|SpotName), ROI.means,REML=FALSE)
freq.model1 = lmer(freq~ Condition + (1|SpotName), ROI.means,REML=FALSE)
freq.anova <- anova(freq.null, freq.model1)
print(freq.anova)

freq2.null = lmer(freq_mean ~ (1|SpotName), Conditions,REML=FALSE)
freq2.model1 = lmer(freq_mean~ Condition + (1|SpotName), Conditions,REML=FALSE)
freq2.anova <- anova(freq2.null, freq2.model1)
print(freq2.anova)

# p values
freq2.Condition <- lsmeans(freq2.model1, pairwise ~ Condition, glhargs=list())
summary(freq2.Condition)

#amplitude
amp.null = lmer(amp_mean ~ (1|SpotName), ROI.means,REML=FALSE)
amp.model1 = lmer(amp_mean~ Condition + (1|SpotName), ROI.means,REML=FALSE)
amp.anova <- anova(amp.null, amp.model1)
print(amp.anova)

# look at residuals
#plot(amp.model1)

# p values
amp.Condition <- lsmeans(amp.model1, pairwise ~ Condition, glhargs=list())
summary(amp.Condition)


#####
# baseline data


#histograms
ggplot(BL_Data, aes(x=Duration, fill=Channel))+ geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of ROI duration")

ggplot(BL_Data, aes(x=halfWidth, fill=Channel))+ geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of ROI duration")

ggplot(BL_Data, aes(x=amplitude, fill=Channel)) + geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("Distribution of ROI amplitudes") 

ggplot(BL_Data, aes(x=peakAUC, fill=Channel)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of peak AUC") 
