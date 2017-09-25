
library("lme4")
library("lmerTest")
library("lattice")
library("plyr")
library("ggplot2")
library("lsmeans")
library("Rmisc")
library("MASS")
library("multcomp")
library("reshape2")
library("data.table")
library("Hmisc")
library("stringr")
library("GGally")

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
# peak data

# load files
baseline1 <- read.delim("E:/Data/Pericyte_project/Two-photon-data/Calcium/Results/ROI_table_18_08_2017.csv", header=TRUE, sep = "\t")
baseline.summary1 <- read.delim("E:/Data/Pericyte_project/Two-photon-data/Calcium/Results/Summary_table_18_08_2017.csv", header=TRUE, sep = "\t")

baseline2 <- read.delim("E:/Data/Pericyte_project/Two-photon-data/Calcium/Results/ROI_table_18_08_2017.csv", header=TRUE, sep = "\t")
baseline.summary2 <- read.delim("E:/Data/Pericyte_project/Two-photon-data/Calcium/Results/Summary_table_18_08_2017.csv", header=TRUE, sep = "\t")

baseline1$Drug="Drug1"
baseline2$Drug="Drug2"

baseline<-rbind(baseline1,baseline2)

#########
# unique animal and spot name
baseline$FOVname <-paste(baseline$animalname, baseline$Spot,sep= "_")
baseline.summary $FOVname <-paste(baseline.summary $animalname, baseline.summary $Spot,sep= "_")

# unique ROI name
baseline$ROIname <-paste(baseline$animalname, baseline$Spot, baseline$ROI, sep= "_")

# ROIType
baseline$ROIType<-"Process"
baseline$ROIType[baseline$is_soma==1]="Soma"
baseline$ROIType<- as.factor(baseline$ROIType)

baseline.summary$ROIType<-"Process"
baseline.summary$ROIType[baseline.summary$is_soma==1]="Soma"
baseline.summary$ROIType<- as.factor(baseline.summary$ROIType)

baseline$ROIType<- factor(baseline$ROIType,levels = c("Soma", "Process"))
#baseline$celltype<- factor(baseline$celltype,levels = c("string"))

###############################
#histograms
ggplot(baseline, aes(x=volume, fill=interaction(ROIType,Drug))) + geom_histogram(binwidth=100, position="dodge") +
  ggtitle("Distribution of ROI volumes")

ggplot(baseline, aes(x=duration, fill=interaction(ROIType,Drug)))+ geom_histogram(binwidth=0.5, position="dodge") +
  ggtitle("Distribution of ROI duration")

ggplot(baseline, aes(x=Max_amplitude, fill=interaction(ROIType,Drug))) + geom_histogram(binwidth=0.5, position="dodge") +
  ggtitle("Distribution of ROI amplitudes")


#########
# pool data for each cell (mean amplitude, total number of signals, etc.)
baseline.means<- ddply(baseline, c("animalname", "Spot","Drug","FOVname","celltype", "ROIType"), 
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
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red")) + 
  max.theme

ggplot(data=df2C, aes(x=Drug, y=Max_amplitude, fill=Drug)) +
  geom_errorbar(aes(ymin=Max_amplitude-se, ymax=Max_amplitude+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red")) + 
  max.theme

ggplot(data=df2D, aes(x=Drug, y=Max_amplitude, fill=ROIType)) +
  geom_errorbar(aes(ymin=Max_amplitude-se, ymax=Max_amplitude+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
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
# p values
amplitude.ROIType <- lsmeans(amplitude.model1, pairwise ~ ROIType, glhargs=list())
summary(amplitude.ROIType)

amplitude.ROIType_PC <- lsmeans(amplitude.model4, pairwise ~ ROIType*Drug, glhargs=list())
summary(amplitude.ROIType_PC)

