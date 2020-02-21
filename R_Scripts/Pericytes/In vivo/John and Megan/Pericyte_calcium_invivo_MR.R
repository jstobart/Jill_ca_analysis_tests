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
RR_baseline1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Cell Scan - Calcium/Drug/Results/RR_peaks_baseline.csv", header=TRUE, sep = ",")

MJ_baseline1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Calcium - Cell Scan/Drugs/Results/MJ_peaks_baseline.csv", header=TRUE, sep = ",")

L_baseline1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/L/Cell Scan - Calcium/2019_05_27/Results/L_peaks_05_27_2019_045_046.csv", header=TRUE, sep = ",")
L_baseline2 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/L/Cell Scan - Calcium/2019_05_28/Results/L_peaks_05_28_2019_051.csv", header=TRUE, sep = ",")
L_baseline3 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/L/Cell Scan - Calcium/2019_06_04/Results/L_peaks_06_04_2019_spot1_2_3.csv", header=TRUE, sep = ",")
L_baseline4 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/L/Cell Scan - Calcium/2019_06_14/Results/L_peaks_06_14_2019_spot1_2.csv", header=TRUE, sep = ",")
L_baseline5 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/L/Cell Scan - Calcium/2019_07_08/Results/L_peaks_07_08_2019_spot1_2.csv", header=TRUE, sep = ",")

JJ_baseline1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/JJ/Cell Scan - Calcium/2019_06_07/Results/JJ_peaks_06_07_2019_008_009.csv", header=TRUE, sep = ",")
JJ_baseline2 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/JJ/Cell Scan - Calcium/2019_06_11/Results/JJ_peaks_06_11_2019_spot2.csv", header=TRUE, sep = ",")

OM_baseline1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/Old Man/Cell Scan - Calcium/2019_07_31/Results/OldMan_peaks_07_31_2019_Spot1.csv", header=TRUE, sep = ",")
OM_baseline2 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/Old Man/Cell Scan - Calcium/2019_08_13/Results/OM_peaks_08_13_2019_SPOT1_2_3.csv", header=TRUE, sep = ",")

NL_baseline1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/New L/Cell Scan - Calcium/2019_07_31/Results/NewL_peaks_07_31_2019_SPOT1_2.csv", header=TRUE, sep = ",")

JJ_baseline1$Channel<-"GCaMP"
JJ_baseline2$Channel<-"GCaMP"
OM_baseline1$Channel<-"GCaMP"
OM_baseline2$Channel<-"GCaMP"
NL_baseline1$Channel<-"GCaMP"

#nimodipine
RR_nimodipine1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Cell Scan - Calcium/Drug/Results/RR_peaks_nimodipine.csv", header=TRUE, sep = ",")

MJ_nimodipine1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Calcium - Cell Scan/Drugs/Results/MJ_peaks_nimodipine.csv", header=TRUE, sep = ",")

L_nimodipine1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/L/Cell Scan - Calcium/2019_07_11 Nimodipine/Results/L_peaks_07_11_2019_spot1.csv", header=TRUE, sep = ",")

#pyr3
RR_pyr3_1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Cell Scan - Calcium/Drug/Results/RR_peaks_pyr3.csv", header=TRUE, sep = ",")

MJ_pyr3_1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Calcium - Cell Scan/Drugs/Results/MJ_peaks_pyr3.csv", header=TRUE, sep = ",")

#combine data together for each treatment
baseline<-rbind(RR_baseline1,MJ_baseline1,L_baseline1,L_baseline2,L_baseline3,L_baseline4,L_baseline5,JJ_baseline1,JJ_baseline2,OM_baseline1,OM_baseline2,NL_baseline1)
nimodipine<-rbind(RR_nimodipine1,MJ_nimodipine1,L_nimodipine1)
pyr3<-rbind(RR_pyr3_1,MJ_pyr3_1)


baseline$Drug="baseline"
nimodipine$Drug="nimodipine"
pyr3$Drug="pyr3"


AllData<-rbind(baseline,nimodipine,pyr3)
AllData$Drug<-as.factor(AllData$Drug)

#pull out the spot name
pos= regexpr('SPOT', AllData$Spot)
AllData$SpotName<-substr(AllData$Spot,pos, pos+5)

#make a unqiue ROI name
AllData$ROI_Spot <-paste(AllData$Animal, AllData$SpotName, AllData$roiName, sep= "_")
AllData$Spot_trial <-paste(AllData$SpotName, AllData$Trial, sep= "_")

AllData$Duration<-AllData$halfWidth*2

AllData$Animal[AllData$Animal=="OldMan"]="OM"

# separate out drug data from spots that only have baseline
DrugROIs<-unique(AllData$ROI_Spot[AllData$Drug=="pyr3"])
DrugData<-subset(AllData, ROI_Spot %in% DrugROIs)

#Baseline Data
BL_Data<-subset(AllData, Drug=="baseline")

###############################
#histograms
ggplot(DrugData, aes(x=Duration, fill=Drug))+ geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of ROI duration")

ggplot(DrugData, aes(x=halfWidth, fill=Drug))+ geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of ROI duration")

ggplot(DrugData, aes(x=amplitude, fill=Drug)) + geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("Distribution of ROI amplitudes") 

ggplot(DrugData, aes(x=peakAUC, fill=Drug)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of peak AUC") 

#########


# pool data for each cell (mean amplitude, total number of signals, etc.)
ROI.means<- ddply(DrugData, c("Animal", "Spot","Drug","ROI_Spot"), 
                  summarise, AUC_mean = mean(peakAUC, na.rm=TRUE), dur_mean = mean(halfWidth,na.rm=TRUE), 
                  amp_mean = mean(amplitude,na.rm=TRUE), nEvents = length(amplitude), nTrials = length(unique(Spot_trial)))

ROI.means$dur_mean[ROI.means$AUC_mean=="NaN"]=0
ROI.means$amp_mean[ROI.means$AUC_mean=="NaN"]=0
ROI.means$nEvents[ROI.means$AUC_mean=="NaN"]=0
ROI.means$AUC_mean[ROI.means$AUC_mean=="NaN"]=0

ROI.means$freq<-ROI.means$nEvents/((ROI.means$nTrials*120)/60) # number of signals/min



################
#frequency
df1A <- summarySE(ROI.means, measurevar="freq", groupvars=c("Drug"))
df1B <- summarySE(ROI.means, measurevar="freq", groupvars=c("Drug","Animal"))

ggplot(data=df1A, aes(x=Drug, y=freq, fill=Drug)) +
  geom_errorbar(aes(ymin=freq-se, ymax=freq+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red", "blue")) +
  max.theme

ggplot(data=df1B, aes(x=Animal, y=freq, fill=Drug)) +
  geom_errorbar(aes(ymin=freq-se, ymax=freq+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red", "blue")) +
  max.theme

#duration
df2A <- summarySE(ROI.means, measurevar="dur_mean", groupvars=c("Drug"))

ggplot(data=df2A, aes(x=Drug, y=dur_mean, fill=Drug)) +
  geom_errorbar(aes(ymin=dur_mean-se, ymax=dur_mean+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("half width (s)") +
  scale_fill_manual(
    values=c("black", "red", "blue")) +
  max.theme

#amplitude
df3A <- summarySE(ROI.means, measurevar="amp_mean", groupvars=c("Drug"))

ggplot(data=df3A, aes(x=Drug, y=amp_mean, fill=Drug)) +
  geom_errorbar(aes(ymin=amp_mean-se, ymax=amp_mean+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude (dF/F)") +
  scale_fill_manual(
    values=c("black", "red", "blue")) +
  max.theme

######
#stats

#frequency
freq.null = lmer(freq ~ (1|Spot), ROI.means,REML=FALSE)
freq.model1 = lmer(freq~ Drug + (1|Spot), ROI.means,REML=FALSE)
freq.anova <- anova(freq.null, freq.model1)
print(freq.anova)

#amplitude
amp.null = lmer(amp_mean ~ (1|Spot), ROI.means,REML=FALSE)
amp.model1 = lmer(amp_mean~ Drug + (1|Spot), ROI.means,REML=FALSE)
amp.anova <- anova(amp.null, amp.model1)
print(amp.anova)

# look at residuals
#plot(amp.model1)

# p values
amp.drug <- lsmeans(amp.model1, pairwise ~ Drug, glhargs=list())
summary(amp.drug)


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
