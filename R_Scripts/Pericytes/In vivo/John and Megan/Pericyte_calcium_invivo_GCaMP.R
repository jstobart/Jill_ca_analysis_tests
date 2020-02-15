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

month9_65 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Calcium Results/GCaMP 65626 L/GCaMP_65626_peaks_9months_02_2020.csv", header=TRUE, sep = ",")

month11_65  <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Calcium Results/GCaMP 65626 L/GCaMP_65626_peaks_11months_02_2020.csv", header=TRUE, sep = ",")

baseline_82 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Calcium Results/GCaMP 82093/GCaMP_82093_peaks_baseline_02_2020.csv", header=TRUE, sep = ",")
nimodipine_82 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Calcium Results/GCaMP 82093/GCaMP_82093_peaks_nimodipine_02_2020.csv", header=TRUE, sep = ",")
pyr3_82 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Calcium Results/GCaMP 82093/GCaMP_82093_peaks_Pyr3_02_2020.csv", header=TRUE, sep = ",")

baseline_79 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Calcium Results/GCaMP 79335/GCaMP_79335_peaks_baseline_02_2020.csv", header=TRUE, sep = ",")
nimodipine_79 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Calcium Results/GCaMP 79335/GCaMP_79335_peaks_nimodipine_02_2020.csv", header=TRUE, sep = ",")
pyr3_79 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Calcium Results/GCaMP 79335/GCaMP_79335_peaks_Pyr3_02_2020.csv", header=TRUE, sep = ",")


# adjust all spot names to be the same
#MJ_baseline1$SpotName<-paste("spot",str_sub(MJ_baseline1$Spot, -1, -1), sep="")


#combine all data
DrugData<-rbind(baseline_79,nimodipine_79,pyr3_79)
AgeData<-rbind(month9_65, month11_65)

#remove spot 2 because it is a problem
DrugData<-subset(DrugData, !(Animal=="79335" & Spot=="spot2"))

#make a unqiue ROI name
DrugData$ROI_Spot <-paste(DrugData$Animal, DrugData$Spot, DrugData$roiName, sep= "_")
DrugData$Spot_trial <-paste(DrugData$Spot, DrugData$Trial, sep= "_")

DrugData$Duration<-DrugData$halfWidth*2

AgeData$ROI_Spot <-paste(AgeData$Animal, AgeData$Spot, AgeData$roiName, sep= "_")
AgeData$Spot_trial <-paste(AgeData$Spot, AgeData$Trial, sep= "_")

AgeData$Duration<-AgeData$halfWidth*2


#Baseline Data
#BL_Data<-subset(AllData, Condition=="baseline")

###############################
#histograms
ggplot(DrugData, aes(x=halfWidth, fill=Condition))+ geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of ROI duration")

ggplot(DrugData, aes(x=amplitude, fill=Condition)) + geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("Distribution of ROI amplitudes") 

ggplot(DrugData, aes(x=peakAUC, fill=Condition)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of peak AUC") 


ggplot(AgeData, aes(x=halfWidth, fill=Condition))+ geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of ROI duration")

ggplot(AgeData, aes(x=amplitude, fill=Condition)) + geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("Distribution of ROI amplitudes") 

ggplot(AgeData, aes(x=peakAUC, fill=Condition)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of peak AUC") 



#########
# drugs

# pool data for each cell (mean amplitude, total number of signals, etc.)
ROI.means<- ddply(DrugData, c("Animal", "Spot","Condition","ROI_Spot"), 
                  summarise, AUC_mean = mean(peakAUC, na.rm=TRUE), dur_mean = mean(halfWidth,na.rm=TRUE), 
                  amp_mean = mean(amplitude,na.rm=TRUE), nEvents = length(amplitude), nTrials = length(unique(Spot_trial)))

ROI.means$dur_mean[ROI.means$AUC_mean=="NaN"]=0
ROI.means$amp_mean[ROI.means$AUC_mean=="NaN"]=0
ROI.means$nEvents[ROI.means$AUC_mean=="NaN"]=0
ROI.means$AUC_mean[ROI.means$AUC_mean=="NaN"]=0

ROI.means$freq<-ROI.means$nEvents/(ROI.means$nTrials*1) # number of signals/min

Spot.means<- ddply(ROI.means, c("Animal", "Spot","Condition"), 
                  summarise, AUC_mean2 = mean(AUC_mean, na.rm=TRUE), dur_mean2 = mean(dur_mean,na.rm=TRUE), 
                  amp_mean2 = mean(amp_mean,na.rm=TRUE), freq_mean= mean(freq, na.rm=TRUE))


# separate out Condition data from spots that only have baseline
ConditionROIs<-unique(ROI.means$ROI_Spot[ROI.means$Condition=="Pyr3"])
ConditionData<-subset(ROI.means, ROI_Spot %in% ConditionROIs)

################
#frequency
df1A <- summarySE(ROI.means, measurevar="freq", groupvars=c("Condition"))
df1B <- summarySE(Spot.means, measurevar="freq_mean", groupvars=c("Condition"))
df1C <- summarySE(ConditionData, measurevar="freq", groupvars=c("Condition"))

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

ggplot(data=df1C, aes(x=Condition, y=freq, fill=Condition)) +
  geom_errorbar(aes(ymin=freq-se, ymax=freq+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
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
df3C <- summarySE(ConditionData, measurevar="amp_mean", groupvars=c("Condition"))

ggplot(data=df3A, aes(x=Condition, y=amp_mean, fill=Condition)) +
  geom_errorbar(aes(ymin=amp_mean-se, ymax=amp_mean+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude (dF/F)") +
  scale_fill_manual(
    values=c("black", "red", "blue", "green")) +
  max.theme

ggplot(data=df3C, aes(x=Condition, y=amp_mean, fill=Condition)) +
  geom_errorbar(aes(ymin=amp_mean-se, ymax=amp_mean+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude (dF/F)") +
  scale_fill_manual(
    values=c("black", "red", "blue", "green")) +
  max.theme


######
#stats

#frequency
freq.null = lmer(freq ~ (1|Animal) + (1|Spot), ROI.means,REML=FALSE)
freq.model1 = lmer(freq~ Condition + (1|Animal) + (1|Spot), ROI.means,REML=FALSE)
freq.anova <- anova(freq.null, freq.model1)
print(freq.anova)

freq.null = lmer(freq ~ (1|Spot), ROI.means,REML=FALSE)
freq.model1 = lmer(freq~ Condition + (1|Spot), ROI.means,REML=FALSE)
freq.anova <- anova(freq.null, freq.model1)
print(freq.anova)

freq2.null = lmer(freq~ (1|Animal) + (1|Spot), ConditionData,REML=FALSE)
freq2.model1 = lmer(freq~ Condition + (1|Animal) + (1|Spot), ConditionData,REML=FALSE)
freq2.anova <- anova(freq2.null, freq2.model1)
print(freq2.anova)

# p values
freq2.Condition <- lsmeans(freq2.model1, pairwise ~ Condition, glhargs=list())
summary(freq2.Condition)

#amplitude
amp.null = lmer(amp_mean ~ (1|Animal) + (1|Spot), ROI.means,REML=FALSE)
amp.model1 = lmer(amp_mean~ Condition + (1|Animal) + (1|Spot), ROI.means,REML=FALSE)
amp.anova <- anova(amp.null, amp.model1)
print(amp.anova)

amp.null = lmer(amp_mean ~ (1|Spot), ROI.means,REML=FALSE)
amp.model1 = lmer(amp_mean~ Condition + (1|Spot), ROI.means,REML=FALSE)
amp.anova <- anova(amp.null, amp.model1)
print(amp.anova)

# look at residuals
#plot(amp.model1)

# p values
amp.Condition <- lsmeans(amp.model1, pairwise ~ Condition, glhargs=list())
summary(amp.Condition)


#####
# aging

# pool data for each cell (mean amplitude, total number of signals, etc.)
ROI.means.age<- ddply(AgeData, c("Animal", "Spot","Condition","ROI_Spot"), 
                  summarise, AUC_mean = mean(peakAUC, na.rm=TRUE), dur_mean = mean(halfWidth,na.rm=TRUE), 
                  amp_mean = mean(amplitude,na.rm=TRUE), nEvents = length(amplitude), nTrials = length(unique(Spot_trial)))

ROI.means.age$dur_mean[ROI.means.age$AUC_mean=="NaN"]=0
ROI.means.age$amp_mean[ROI.means.age$AUC_mean=="NaN"]=0
ROI.means.age$nEvents[ROI.means.age$AUC_mean=="NaN"]=0
ROI.means.age$AUC_mean[ROI.means.age$AUC_mean=="NaN"]=0

ROI.means.age$freq<-ROI.means.age$nEvents/(ROI.means.age$nTrials*1) # number of signals/min

Spot.means.age<- ddply(ROI.means.age, c("Animal", "Spot","Condition"), 
                   summarise, AUC_mean2 = mean(AUC_mean, na.rm=TRUE), dur_mean2 = mean(dur_mean,na.rm=TRUE), 
                   amp_mean2 = mean(amp_mean,na.rm=TRUE), freq_mean= mean(freq, na.rm=TRUE))


# separate out Condition data from spots that only have baseline
ConditionROIs.age<-unique(ROI.means.age$ROI_Spot[ROI.means.age$Condition=="11months"])
ConditionData.age<-subset(ROI.means.age, ROI_Spot %in% ConditionROIs.age)

################
#frequency
df1A.age <- summarySE(ROI.means.age, measurevar="freq", groupvars=c("Condition"))
df1B.age <- summarySE(Spot.means.age, measurevar="freq_mean", groupvars=c("Condition"))
df1C.age <- summarySE(ConditionData.age, measurevar="freq", groupvars=c("Condition"))

ggplot(data=df1A.age, aes(x=Condition, y=freq, fill=Condition)) +
  geom_errorbar(aes(ymin=freq-se, ymax=freq+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red", "blue","green")) +
  max.theme

ggplot(data=df1B.age, aes(x=Condition, y=freq_mean, fill=Condition)) +
  geom_errorbar(aes(ymin=freq_mean-se, ymax=freq_mean+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red", "blue","green")) +
  max.theme

ggplot(data=df1C.age, aes(x=Condition, y=freq, fill=Condition)) +
  geom_errorbar(aes(ymin=freq-se, ymax=freq+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red", "blue","green")) +
  max.theme

#duration
df2A.age <- summarySE(ROI.means.age, measurevar="dur_mean", groupvars=c("Condition"))

ggplot(data=df2A.age, aes(x=Condition, y=dur_mean, fill=Condition)) +
  geom_errorbar(aes(ymin=dur_mean-se, ymax=dur_mean+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("half width (s)") +
  scale_fill_manual(
    values=c("black", "red", "blue","green")) +
  max.theme

#amplitude
df3A.age <- summarySE(ROI.means.age, measurevar="amp_mean", groupvars=c("Condition"))
df3C.age <- summarySE(ConditionData.age, measurevar="amp_mean", groupvars=c("Condition"))

ggplot(data=df3A.age, aes(x=Condition, y=amp_mean, fill=Condition)) +
  geom_errorbar(aes(ymin=amp_mean-se, ymax=amp_mean+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude (dF/F)") +
  scale_fill_manual(
    values=c("black", "red", "blue", "green")) +
  max.theme

ggplot(data=df3C.age, aes(x=Condition, y=amp_mean, fill=Condition)) +
  geom_errorbar(aes(ymin=amp_mean-se, ymax=amp_mean+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude (dF/F)") +
  scale_fill_manual(
    values=c("black", "red", "blue", "green")) +
  max.theme

######
#stats

#frequency
freq.age.null = lmer(freq ~ (1|Spot), ROI.means.age,REML=FALSE)
freq.age.model1 = lmer(freq~ Condition + (1|Spot), ROI.means.age,REML=FALSE)
freq.age.anova <- anova(freq.age.null, freq.age.model1)
print(freq.age.anova)

freq2.age.null = lmer(freq~ (1|Spot), ConditionData.age,REML=FALSE)
freq2.age.model1 = lmer(freq~ Condition + (1|Spot), ConditionData.age,REML=FALSE)
freq2.age.anova <- anova(freq2.age.null, freq2.age.model1)
print(freq2.age.anova)

# p values
freq2.age.Condition <- lsmeans(freq2.age.model1, pairwise ~ Condition, glhargs=list())
summary(freq2.age.Condition)

#amplitude
amp.age.null = lmer(amp_mean ~ (1|Spot), ROI.means.age,REML=FALSE)
amp.age.model1 = lmer(amp_mean~ Condition + (1|Spot), ROI.means.age,REML=FALSE)
amp.age.anova <- anova(amp.age.null, amp.age.model1)
print(amp.anova)

# look at residuals
#plot(amp.model1)

# p values
amp.age.Condition <- lsmeans(amp.age.model1, pairwise ~ Condition, glhargs=list())
summary(amp.age.Condition)
