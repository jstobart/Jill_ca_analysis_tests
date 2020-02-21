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
RR_baseline1 <- read.table("C:/Data/Results/Pericytes/RR/Results/Ca/2019_07_22/Results/RR_peaks_07_22_2019_SPOT4_5.csv", header=TRUE, sep = ",")

#nimodipine
RR_nimodipine1 <- read.table("C:/Data/Results/Pericytes/RR/Results/Ca/2019_08_02 Nimodipine/Results/RR_peaks_08_02_2019_SPOT4_5.csv", header=TRUE, sep = ",")

#pyr3
RR_pyr3_1<-

#combine data together for each treatment
baseline<-rbind(RR_baseline1)
nimodipine<-rbind(RR_nimodipine1)



baseline$Drug="baseline"
nimodipine$Drug="nimodipine"

# combine all the data together
AllData<-rbind(baseline,nimodipine)
AllData$Drug<-is.factor(AllData$Drug)

#pull out the spot name
pos= regexpr('SPOT', AllData$Spot)
AllData$SpotName<-substr(AllData$Spot,pos, pos+5)

#make a unqiue ROI name
AllData$ROI_Spot <-paste(AllData$Animal, AllData$SpotName, AllData$roiName, sep= "_")

AllData$Duration<-AllData$halfWidth*2

###############################
#histograms
ggplot(AllData, aes(x=Duration, fill=Drug))+ geom_histogram(binwidth=0.5, position="dodge") +
  ggtitle("Distribution of ROI duration")

ggplot(AllData, aes(x=amplitude, fill=Drug)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of ROI amplitudes") + xlim(0, 40)


#########
# pool data for each ROI (mean amplitude, total number of signals, etc.)
AllData.means<- ddply(AllData, c("Animal", "Spot","Drug","FOVname","celltype", "ROIType"), 
                   summarise, AUC_mean = mean(auc, na.rm=TRUE), dur_mean = mean(duration,na.rm=TRUE), 
                   amp_mean = mean(Max_amplitude,na.rm=TRUE),
                   nEvents = length(Max_amplitude), volume_mean=mean(volume), area_mean=mean(area), 
                   propDist=mean(centroidDis_Traveled), propRate=mean(centroidProp_Rate))

baseline.means$freq<-(baseline.means$nEvents/baseline.summary$img_duration[1])*60 # number of signals/min

################
#frequency
df1A <- summarySE(baseline.means, measurevar="freq", groupvars=c("ROIType"))
df1B <- summarySE(baseline.means, measurevar="freq", groupvars=c("celltype","ROIType"))
df1C <- summarySE(baseline.means, measurevar="freq", groupvars=c("Drug"))
df1D <- summarySE(baseline.means, measurevar="freq", groupvars=c("Drug","ROIType"))

ggplot(data=df1A, aes(x=ROIType, y=freq, fill=ROIType)) +
  geom_errorbar(aes(ymin=freq-se, ymax=freq+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red")) +
  max.theme

ggplot(data=df1B, aes(x=celltype, y=freq, fill=ROIType)) +
  geom_errorbar(aes(ymin=freq-se, ymax=freq+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red")) + 
  max.theme

ggplot(data=df1C, aes(x=Drug, y=freq, fill=Drug)) +
  geom_errorbar(aes(ymin=freq-se, ymax=freq+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red")) + 
  max.theme

ggplot(data=df1D, aes(x=Drug, y=freq, fill=ROIType)) +
  geom_errorbar(aes(ymin=freq-se, ymax=freq+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red")) + 
  max.theme

#frequency
# mean of each cell
freq.null = lmer(freq ~ (1|animalname) + (1|Spot), baseline.means,REML=FALSE)
freq.model1 = lmer(freq~ ROIType + (1|animalname) + (1|Spot), baseline.means,REML=FALSE)
freq.model2 = lmer(freq~ Drug + (1|animalname) + (1|Spot), baseline.means,REML=FALSE)
freq.model3 = lmer(freq~ ROIType + Drug + (1|animalname) + (1|Spot), baseline.means,REML=FALSE)
freq.model4 = lmer(freq~ ROIType * Drug + (1|animalname) + (1|Spot), baseline.means,REML=FALSE)
freq.anova <- anova(freq.null, freq.model1,freq.model2,freq.model3,freq.model4)
print(freq.anova)

# look at residuals
plot(freq.model1)
plot(freq.model2)
plot(freq.model3)
plot(freq.model4)

# p values
freq.ROIType <- lsmeans(freq.model1, pairwise ~ ROIType, glhargs=list())
summary(freq.ROIType)

freq.ROIType_PC <- lsmeans(freq.model4, pairwise ~ ROIType*Drug, glhargs=list())
summary(freq.ROIType_PC)


########
#amplitude
df2A <- summarySE(baseline, measurevar="Max_amplitude", groupvars=c("ROIType"),na.rm=TRUE)
df2B <- summarySE(baseline, measurevar="Max_amplitude", groupvars=c("celltype","ROIType"),na.rm=TRUE)
df2C <- summarySE(baseline, measurevar="Max_amplitude", groupvars=c("Drug"),na.rm=TRUE)
df2D <- summarySE(baseline, measurevar="Max_amplitude", groupvars=c("Drug","ROIType"),na.rm=TRUE)


ggplot(data=df2A, aes(x=ROIType, y=Max_amplitude, fill=ROIType)) +
  geom_errorbar(aes(ymin=Max_amplitude-se, ymax=Max_amplitude+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude") +
  scale_fill_manual(
    values=c("black", "red")) +
  max.theme

ggplot(data=df2B, aes(x=celltype, y=Max_amplitude, fill=ROIType)) +
  geom_errorbar(aes(ymin=Max_amplitude-se, ymax=Max_amplitude+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude") +
  scale_fill_manual(
    values=c("black", "red")) + 
  max.theme

ggplot(data=df2C, aes(x=Drug, y=Max_amplitude, fill=Drug)) +
  geom_errorbar(aes(ymin=Max_amplitude-se, ymax=Max_amplitude+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude") +
  scale_fill_manual(
    values=c("black", "red")) + 
  max.theme

ggplot(data=df2D, aes(x=Drug, y=Max_amplitude, fill=ROIType)) +
  geom_errorbar(aes(ymin=Max_amplitude-se, ymax=Max_amplitude+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude") +
  scale_fill_manual(
    values=c("black", "red")) + 
  max.theme


ggplot(data=baseline, aes(x=Drug, y= Max_amplitude, fill=ROIType)) +
  geom_boxplot() +
  scale_fill_manual(
    values=c("black", "red")) + 
  max.theme


#amplitude stats
# all ROIs
amplitude.null = lmer(Max_amplitude ~ (1|animalname) + (1|Spot), baseline,REML=FALSE)
amplitude.model1 = lmer(Max_amplitude~ ROIType + (1|animalname) + (1|Spot), baseline,REML=FALSE)
amplitude.model2 = lmer(Max_amplitude~ Drug + (1|animalname) + (1|Spot), baseline,REML=FALSE)
amplitude.model3 = lmer(Max_amplitude~ ROIType + Drug + (1|animalname) + (1|Spot), baseline,REML=FALSE)
amplitude.model4 = lmer(Max_amplitude~ ROIType * Drug + (1|animalname) + (1|Spot), baseline,REML=FALSE)
amplitude.anova <- anova(amplitude.null, amplitude.model1,amplitude.model2,amplitude.model3,amplitude.model4)
print(amplitude.anova)

# look at residuals
plot(amplitude.model1)
plot(amplitude.model2)
plot(amplitude.model3)
plot(amplitude.model4)


dotplot(ranef(amplitude.null, condVar = TRUE))$animalname
dotplot(ranef(amplitude.model1, condVar = TRUE))$animalname
dotplot(ranef(amplitude.model2, condVar = TRUE))$animalname
dotplot(ranef(amplitude.model3, condVar = TRUE))$animalname
dotplot(ranef(amplitude.model4, condVar = TRUE))$animalname

dotplot(ranef(amplitude.model4, condVar = TRUE))$Spot

# p values
amplitude.ROIType <- lsmeans(amplitude.model1, pairwise ~ ROIType, glhargs=list())
summary(amplitude.ROIType)

amplitude.ROIType_PC <- lsmeans(amplitude.model4, pairwise ~ ROIType*Drug, glhargs=list())
summary(amplitude.ROIType_PC)

