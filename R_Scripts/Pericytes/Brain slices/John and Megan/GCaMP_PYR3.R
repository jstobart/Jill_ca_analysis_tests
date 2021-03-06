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
GCaMPSlice_baseline_2019_07_11 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_11_2019_Baseline.csv", header=TRUE, sep = ",")
GCaMPSlice_baseline_2019_07_12 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_12_2019_Baseline.csv", header=TRUE, sep = ",")
GCaMPSlice_baseline_2019_07_15 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_15_2019_Baseline.csv", header=TRUE, sep = ",")
GCaMPSlice_baseline_2019_07_16 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_16_2019_Baseline.csv", header=TRUE, sep = ",")
GCaMPSlice_baseline_2019_07_17 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_17_2019_baseline.csv", header=TRUE, sep = ",")
#GCaMPSlice_TBX_2019_07_17 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_17_2019_TBX.csv", header=TRUE, sep = ",")

#pyr3
GCaMPSlice_pyr3_2019_07_11 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_11_2019_PYR3.csv", header=TRUE, sep = ",")
GCaMPSlice_pyr3_2019_07_12 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_12_2019_PYR3.csv", header=TRUE, sep = ",")
GCaMPSlice_pyr3_2019_07_15 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_15_2019_PYR3.csv", header=TRUE, sep = ",")
GCaMPSlice_pyr3_2019_07_16 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_16_2019_pyr3.csv", header=TRUE, sep = ",")
GCaMPSlice_pyr3_2019_07_17 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_17_2019_pyr3.csv", header=TRUE, sep = ",")
#GCaMPSlice_pyr3_2019_07_17 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_17_2019_PTBX.csv", header=TRUE, sep = ",")

#pyr3tbx
GCaMPSlice_pyr3tbx_2019_07_11 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_11_2019_PTBX.csv", header=TRUE, sep = ",")
GCaMPSlice_pyr3tbx_2019_07_12 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_12_2019_PTBX.csv", header=TRUE, sep = ",")
GCaMPSlice_pyr3tbx_2019_07_15 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_15_2019_PTBX.csv", header=TRUE, sep = ",")
GCaMPSlice_pyr3tbx_2019_07_16 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_16_2019_PTBX.csv", header=TRUE, sep = ",")
GCaMPSlice_pyr3tbx_2019_07_17 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_17_2019_PTBX.csv", header=TRUE, sep = ",")

#tbx
GCaMPSlice_tbx_2019_07_11 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_11_2019_TBX.csv", header=TRUE, sep = ",")
GCaMPSlice_tbx_2019_07_12 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_12_2019_TBX.csv", header=TRUE, sep = ",")
GCaMPSlice_tbx_2019_07_15 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_15_2019_TBX.csv", header=TRUE, sep = ",")
GCaMPSlice_tbx_2019_07_16 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_16_2019_TBX.csv", header=TRUE, sep = ",")
GCaMPSlice_tbx_2019_07_17 <- read.table("C:/Users/John PC/Documents/Summer Research 2019/R-analysis/GCaMP PYR3/GCaMPSlice_peaks_07_17_2019_TBX.csv", header=TRUE, sep = ",")


#pyr3
#RR_pyr3_1<-

#combine data together for each treatment
baseline<-rbind(GCaMPSlice_baseline_2019_07_11, GCaMPSlice_baseline_2019_07_12, GCaMPSlice_baseline_2019_07_15, GCaMPSlice_baseline_2019_07_16, GCaMPSlice_baseline_2019_07_17)
pyr3<-rbind(GCaMPSlice_pyr3_2019_07_11, GCaMPSlice_pyr3_2019_07_12, GCaMPSlice_pyr3_2019_07_15, GCaMPSlice_pyr3_2019_07_16, GCaMPSlice_pyr3_2019_07_17)
pyr3tbx<-rbind(GCaMPSlice_pyr3tbx_2019_07_11, GCaMPSlice_pyr3tbx_2019_07_12, GCaMPSlice_pyr3tbx_2019_07_15, GCaMPSlice_pyr3tbx_2019_07_16, GCaMPSlice_pyr3tbx_2019_07_17)
tbx<-rbind(GCaMPSlice_tbx_2019_07_11, GCaMPSlice_tbx_2019_07_12, GCaMPSlice_tbx_2019_07_15, GCaMPSlice_tbx_2019_07_15, GCaMPSlice_tbx_2019_07_16, GCaMPSlice_tbx_2019_07_17)


baseline$Drug="baseline"
pyr3$Drug="pyr3"
pyr3tbx$Drug="pyr3 + tbx"
tbx$Drug="tbx"


AllData<-rbind(baseline,pyr3,pyr3tbx,tbx)
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
    values=c("black", "red", "blue", "cyan")) +
  max.theme

#duration
df2A <- summarySE(ROI.means, measurevar="dur_mean", groupvars=c("Drug"))

ggplot(data=df2A, aes(x=Drug, y=dur_mean, fill=Drug)) +
  geom_errorbar(aes(ymin=dur_mean-se, ymax=dur_mean+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("half width (s)") +
  scale_fill_manual(
    values=c("black", "red", "blue", "cyan")) +
  max.theme

#amplitude
df3A <- summarySE(ROI.means, measurevar="amp_mean", groupvars=c("Drug"))

ggplot(data=df3A, aes(x=Drug, y=amp_mean, fill=Drug)) +
  geom_errorbar(aes(ymin=amp_mean-se, ymax=amp_mean+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude (dF/F)") +
  scale_fill_manual(
    values=c("black", "red", "blue", "cyan")) +
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


