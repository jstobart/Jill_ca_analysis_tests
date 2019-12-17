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
RCaMPSlice_baseline_2019_06_21 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/RCaMP Baseline PYR3 files/RcampSlice_peaks_06_21_2019_baseline.csv", header=TRUE, sep = ",")
RCaMPSlice_baseline_2019_06_24 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/RCaMP Baseline PYR3 files/RCaMPSlice_peaks_06_24_2019_Baseline.csv", header=TRUE, sep = ",")
RCaMPSlice_baseline_2019_06_25 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/RCaMP Baseline PYR3 files/RCaMPSlice_peaks_06_25_2019_Baseline.csv", header=TRUE, sep = ",")
RCaMPSlice_baseline_2019_06_26 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/RCaMP Baseline PYR3 files/RCaMPSlice_peaks_06_26_2019_Baseline.csv", header=TRUE, sep = ",")
RCaMPSlice_baseline_2019_07_09 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/RCaMP Baseline PYR3 files/RCaMPSlice_peaks_07_09_2019_Baseline.csv", header=TRUE, sep = ",")

#pyr3
RCaMPSlice_pyr3_2019_06_21 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/RCaMP Baseline PYR3 files/RcampSlice_peaks_06_21_2019_pyr3.csv", header=TRUE, sep = ",")
RCaMPSlice_pyr3_2019_06_24 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/RCaMP Baseline PYR3 files/RCaMPSlice_peaks_06_24_2019_PYR3.csv", header=TRUE, sep = ",")
RCaMPSlice_pyr3_2019_06_25 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/RCaMP Baseline PYR3 files/RCaMPSlice_peaks_06_25_2019_PYR3.csv", header=TRUE, sep = ",")
RCaMPSlice_pyr3_2019_06_26 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/RCaMP Baseline PYR3 files/RCaMPSlice_peaks_06_26_2019_PYR3.csv", header=TRUE, sep = ",")
RCaMPSlice_pyr3_2019_07_09 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/RCaMP Baseline PYR3 files/RCaMPSlice_peaks_07_09_2019_PYR3.csv", header=TRUE, sep = ",")



#pyr3
#RR_pyr3_1<-

#combine data together for each treatment
baseline<-rbind(RCaMPSlice_baseline_2019_06_26, RCaMPSlice_baseline_2019_06_25, RCaMPSlice_baseline_2019_06_24, RCaMPSlice_baseline_2019_06_21, RCaMPSlice_baseline_2019_07_09)
pyr3<-rbind(RCaMPSlice_pyr3_2019_06_21, RCaMPSlice_pyr3_2019_06_24, RCaMPSlice_pyr3_2019_06_25, RCaMPSlice_pyr3_2019_06_26, RCaMPSlice_pyr3_2019_07_09)



baseline$Drug="baseline"
pyr3$Drug="pyr3"



AllData<-rbind(baseline,pyr3)
AllData$Drug<-as.factor(AllData$Drug)

#make a unqiue ROI name
AllData$ROI_Spot <-paste(AllData$Spot, AllData$roiName, sep= "_")
AllData$Spot_trial <-paste(AllData$Spot, AllData$Trial, sep= "_")

AllData$Duration<-AllData$halfWidth*2

###############################
#histograms
ggplot(AllData, aes(x=Duration, fill=Drug))+ geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of ROI duration")

ggplot(AllData, aes(x=halfWidth, fill=Drug))+ geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of ROI duration")

ggplot(AllData, aes(x=amplitude, fill=Drug)) + geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("Distribution of ROI amplitudes") 

ggplot(AllData, aes(x=peakAUC, fill=Drug)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of peak AUC") 

#########


# pool data for each cell (mean amplitude, total number of signals, etc.)
ROI.means<- ddply(AllData, c("Animal", "Spot","Drug","ROI_Spot"), 
                   summarise, AUC_mean = mean(peakAUC, na.rm=TRUE), dur_mean = mean(halfWidth,na.rm=TRUE), 
                   amp_mean = mean(amplitude,na.rm=TRUE), nEvents = length(amplitude), nTrials = length(unique(Spot_trial)))

ROI.means$freq<-ROI.means$nEvents/((ROI.means$nTrials*120)/60) # number of signals/min

ROI.means$dur_mean[ROI.means$AUC_mean=="NaN"]=0
ROI.means$amp_mean[ROI.means$AUC_mean=="NaN"]=0
ROI.means$nEvents[ROI.means$AUC_mean=="NaN"]=0
ROI.means$AUC_mean[ROI.means$AUC_mean=="NaN"]=0

################
#frequency
df1A <- summarySE(ROI.means, measurevar="freq", groupvars=c("Drug"))

ggplot(data=df1A, aes(x=Drug, y=freq, fill=Drug)) +
  geom_errorbar(aes(ymin=freq-se, ymax=freq+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red")) +
  max.theme

#duration
df2A <- summarySE(ROI.means, measurevar="dur_mean", groupvars=c("Drug"))

ggplot(data=df2A, aes(x=Drug, y=dur_mean, fill=Drug)) +
  geom_errorbar(aes(ymin=dur_mean-se, ymax=dur_mean+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("half width (s)") +
  scale_fill_manual(
    values=c("black", "red")) +
  max.theme

#amplitude
df3A <- summarySE(ROI.means, measurevar="amp_mean", groupvars=c("Drug"))

ggplot(data=df3A, aes(x=Drug, y=amp_mean, fill=Drug)) +
  geom_errorbar(aes(ymin=amp_mean-se, ymax=amp_mean+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude (dF/F)") +
  scale_fill_manual(
    values=c("black", "red")) +
  max.theme

######
#stats

#amplitude
amp.null = lmer(amp_mean ~ (1|Spot), ROI.means,REML=FALSE)
amp.model1 = lmer(amp_mean~ Drug + (1|Spot), ROI.means,REML=FALSE)
amp.anova <- anova(amp.null, amp.model1)
print(amp.anova)

amp.drug <- lsmeans(amp.model1, pairwise ~ Drug, glhargs=list())
summary(amp.drug)

#frequency
freq.null = lmer(freq ~ (1|Spot), ROI.means,REML=FALSE)
freq.model1 = lmer(freq~ Drug + (1|Spot), ROI.means,REML=FALSE)
freq.anova <- anova(freq.null, freq.model1)
print(freq.anova)

freq.drug <- lsmeans(freq.model1, pairwise ~ Drug, glhargs=list())
summary(freq.drug)

#duration
dur.null = lmer(dur_mean ~ (1|Spot), ROI.means,REML=FALSE)
dur.model1 = lmer(dur_mean~ Drug + (1|Spot), ROI.means,REML=FALSE)
dur.anova <- anova(dur.null, dur.model1)
print(dur.anova)

dur.drug <- lsmeans(dur.model1, pairwise ~ Drug, glhargs=list())
summary(dur.drug)


# look at residuals
#plot(amp.model1)
plot(freq.model1)
plot(dur.model1)
# p values
amp.drug <- lsmeans(amp.model1, pairwise ~ Drug, glhargs=list())
summary(amp.drug)


