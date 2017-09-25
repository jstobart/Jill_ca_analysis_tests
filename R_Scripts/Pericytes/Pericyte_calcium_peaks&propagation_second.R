
library("lme4")
library("lmerTest")
library("lattice")
library("plyr")
library("ggplot2")
library("gplots")
library("lsmeans")
library("Rmisc")
library("MASS")
library("multcomp")
library("reshape2")
library("data.table")
library("Hmisc")
library("stringr")
library("spatstat")

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
peaks <- read.table("D:/cre pdgfrb_GCaMP/Data ordered for automated analysis/Results/peaks_05_04_2017.csv", header=TRUE, sep = ",")

# unique ROI name
peaks$ROIs <-paste(peaks$Animal, peaks$Spot, peaks$ROIname, sep= "_")

# ROIType
peaks$ROIType<-"Soma"
process=grepl("P*",peaks$ROIname)
peaks$ROIType[process]="Process"
peaks$ROIType<- as.factor(peaks$ROIType)

peaks$ROIType<- factor(peaks$ROIType,levels = c("Soma", "Process"))
peaks$CellType<- factor(peaks$CellType,levels = c("SMC","Mesh", "Strand","Bifurcation","Transitional"))

#duration
peaks$Duration = peaks$halfWidth*2


###############################

# aggregate trials of all ROIs
all.trials<- ddply(peaks, c("Animal", "Spot", "ROIs","trial","CellType", "ROIType"), 
                   summarise, PA_mean = mean(peakAUC),Dur_mean = mean(Duration), Prom_mean = mean(prominence),
                   amp_mean = mean(amplitude), nEvents = length(peakAUC))

all.trials$freq<-all.trials$nEvents/1.5 # number of signals/min

# aggregate ROIs
all.ROIs<- ddply(all.trials, c("Animal", "Spot", "ROIs","CellType","ROIType"), summarise, 
                 PA_mean2 = mean(PA_mean),freq_mean2 = mean(freq),
                 Dur_mean2 = mean(Dur_mean),Prom_mean2 = mean(Prom_mean), 
                 amp_mean2 = mean(amp_mean), TotEvents = sum(nEvents),
                 MeanEvents = mean(nEvents), N= length(PA_mean))



################
#histograms
ggplot(all.ROIs, aes(x=freq_mean2, fill=ROIType)) + geom_histogram(binwidth=0.5, position="dodge") +
  ggtitle("Distribution of frequencies")
ggplot(all.ROIs, aes(x=freq_mean2, fill=CellType)) + geom_histogram(binwidth=0.5, position="dodge") +
  ggtitle("Distribution of frequencies")
ggplot(all.ROIs, aes(x=freq_mean2, fill=interaction(CellType,ROIType))) + geom_histogram(binwidth=0.5, position="dodge") +
  ggtitle("Distribution of frequencies")


############################
#frequency
df1A <- summarySE(all.ROIs, measurevar="freq_mean2", groupvars=c("ROIType"))
df1B <- summarySE(all.ROIs, measurevar="freq_mean2", groupvars=c("CellType"))
df1C <- summarySE(all.ROIs, measurevar="freq_mean2", groupvars=c("CellType","ROIType"))


ggplot(data=df1A, aes(x=ROIType, y=freq_mean2, fill=ROIType)) +
  geom_errorbar(aes(ymin=freq_mean2-se, ymax=freq_mean2+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red")) +
  max.theme

ggplot(data=df1B, aes(x=CellType, y=freq_mean2, fill=CellType)) +
  geom_errorbar(aes(ymin=freq_mean2-se, ymax=freq_mean2+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red","blue", "green", "yellow")) + 
  max.theme

ggplot(data=df1C, aes(x=CellType, y=freq_mean2, fill=ROIType)) +
  geom_errorbar(aes(ymin=freq_mean2-se, ymax=freq_mean2+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red")) + 
  max.theme

#frequency
# all ROIs
freq.null = lmer(freq_mean2 ~ (1|Animal) + (1|Spot), all.ROIs,REML=FALSE)
freq.model1 = lmer(freq_mean2~ ROIType + (1|Animal) + (1|Spot), all.ROIs,REML=FALSE)
freq.model2 = lmer(freq_mean2~ CellType + (1|Animal) + (1|Spot), all.ROIs,REML=FALSE)
freq.model3 = lmer(freq_mean2~ ROIType + CellType + (1|Animal) + (1|Spot), all.ROIs,REML=FALSE)
freq.model4 = lmer(freq_mean2~ ROIType * CellType + (1|Animal) + (1|Spot), all.ROIs,REML=FALSE)
freq.anova <- anova(freq.null, freq.model1,freq.model2,freq.model3,freq.model4)
print(freq.anova)
# p values
freq.ROIType <- lsmeans(freq.model1, pairwise ~ ROIType, glhargs=list())
summary(freq.ROIType)

freq.CellType <- lsmeans(freq.model2, pairwise ~ CellType, glhargs=list())
summary(freq.CellType)

freq.ROIType_PC <- lsmeans(freq.model4, pairwise ~ ROIType*CellType, glhargs=list())
summary(freq.ROIType_PC)



#amplitude
df2A <- summarySE(all.ROIs, measurevar="amp_mean2", groupvars=c("ROIType"))
df2B <- summarySE(all.ROIs, measurevar="amp_mean2", groupvars=c("CellType"))
df2C <- summarySE(all.ROIs, measurevar="amp_mean2", groupvars=c("CellType","ROIType"))


ggplot(data=df2A, aes(x=ROIType, y=amp_mean2, fill=ROIType)) +
  geom_errorbar(aes(ymin=amp_mean2-se, ymax=amp_mean2+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude") +
  scale_fill_manual(
    values=c("black", "red")) +
  max.theme

ggplot(data=df2B, aes(x=CellType, y=amp_mean2, fill=CellType)) +
  geom_errorbar(aes(ymin=amp_mean2-se, ymax=amp_mean2+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude") +
  scale_fill_manual(
    values=c("black", "red","blue")) +
  max.theme

ggplot(data=df2C, aes(x=CellType, y=amp_mean2, fill=ROIType)) +
  geom_errorbar(aes(ymin=amp_mean2-se, ymax=amp_mean2+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude") +
  scale_fill_manual(
    values=c("black", "red")) +
  max.theme

#amplitude
# all ROIs
amp.null = lmer(amp_mean2 ~ (1|Animal) + (1|Spot), all.ROIs,REML=FALSE)
amp.model1 = lmer(amp_mean2~ ROIType + (1|Animal) + (1|Spot), all.ROIs,REML=FALSE)
amp.model2 = lmer(amp_mean2~ ROIType + (1|Animal) + (1|Spot), all.ROIs,REML=FALSE)
amp.model3 = lmer(amp_mean2~ ROIType + CellType + (1|Animal) + (1|Spot), all.ROIs,REML=FALSE)
amp.model4 = lmer(amp_mean2~ ROIType * CellType + (1|Animal) + (1|Spot), all.ROIs,REML=FALSE)
amp.anova <- anova(amp.null, amp.model1,amp.model2,amp.model3,amp.model4)
print(amp.anova)
# p values
amp.ROIType <- lsmeans(amp.model1, pairwise ~ ROIType, glhargs=list())
summary(amp.ROIType)

amp.CellType <- lsmeans(amp.model2, pairwise ~ CellType, glhargs=list())
summary(amp.CellType)

amp.ROIType_PC <- lsmeans(amp.model4, pairwise ~ ROIType*CellType, glhargs=list())
summary(amp.ROIType_PC)


#duration
df3A <- summarySE(all.ROIs, measurevar="Dur_mean2", groupvars=c("ROIType"))
df3B <- summarySE(all.ROIs, measurevar="Dur_mean2", groupvars=c("CellType"))
df3C <- summarySE(all.ROIs, measurevar="Dur_mean2", groupvars=c("CellType","ROIType"))

ggplot(data=df3A, aes(x=ROIType, y=Dur_mean2, fill=ROIType)) +
  geom_errorbar(aes(ymin=Dur_mean2-se, ymax=Dur_mean2+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("duration (s)") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme

ggplot(data=df3B, aes(x=CellType, y=Dur_mean2, fill=CellType)) +
  geom_errorbar(aes(ymin=Dur_mean2-se, ymax=Dur_mean2+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("duration (s)") +
  scale_fill_manual(
    values=c("black", "red","blue", "green", "yellow"))+
  max.theme

ggplot(data=df3C, aes(x=CellType, y=Dur_mean2, fill=ROIType)) +
  geom_errorbar(aes(ymin=Dur_mean2-se, ymax=Dur_mean2+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("duration (s)") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme

#duration
# all ROIs
Dur.null = lmer(Dur_mean2 ~ (1|Animal) + (1|Spot), all.ROIs,REML=FALSE)
Dur.model1 = lmer(Dur_mean2~ ROIType + (1|Animal) + (1|Spot), all.ROIs,REML=FALSE)
Dur.model2 = lmer(Dur_mean2~ CellType + CellType + (1|Animal) + (1|Spot), all.ROIs,REML=FALSE)
Dur.model3 = lmer(Dur_mean2~ ROIType + CellType + (1|Animal) + (1|Spot), all.ROIs,REML=FALSE)
Dur.model4 = lmer(Dur_mean2~ ROIType * CellType + (1|Animal) + (1|Spot), all.ROIs,REML=FALSE)
Dur.anova <- anova(Dur.null, Dur.model1,Dur.model2,Dur.model3,Dur.model4)
print(Dur.anova)
# p values
Dur.ROIType <- lsmeans(Dur.model1, pairwise ~ ROIType, glhargs=list())
summary(Dur.ROIType)

Dur.CellType <- lsmeans(Dur.model2, pairwise ~ CellType, glhargs=list())
summary(Dur.CellType)

Dur.ROIType_PC <- lsmeans(Dur.model4, pairwise ~ ROIType*CellType, glhargs=list())
summary(Dur.ROIType_PC)


#########################
# PROPAGATION

# load files
propagation <- read.table("D:/cre pdgfrb_GCaMP/Data ordered for automated analysis/Results/propagation_05_04_2017.csv", header=TRUE, sep = ",")

# unique spot name
propagation$Spot2 <-paste(propagation$Animal, propagation$Spot, sep= "_")

# Distance of propagation (um)
df4 <- summarySE(propagation, measurevar="distance", groupvars=c("CellType"))

ggplot(data=df4, aes(x=CellType, y=distance, fill=CellType)) +
  geom_errorbar(aes(ymin=distance-se, ymax=distance+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("propagation distance") +
  scale_fill_manual(
    values=c("black", "red","blue", "green", "yellow")) +
  max.theme

# Duration
df5 <- summarySE(propagation, measurevar="duration", groupvars=c("CellType"))

ggplot(data=df5, aes(x=CellType, y=duration, fill=CellType)) +
  geom_errorbar(aes(ymin=duration-se, ymax=duration+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("propagation duration") +
  scale_fill_manual(
    values=c("black", "red","blue", "green", "yellow")) +
  max.theme

# Rate of propagation- um/s
df6 <- summarySE(propagation, measurevar="propRate", groupvars=c("CellType"))

ggplot(data=df6, aes(x=CellType, y=propRate, fill=CellType)) +
  geom_errorbar(aes(ymin=propRate-se, ymax=propRate+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("propagation propRate") +
  scale_fill_manual(
    values=c("black", "red","blue", "green", "yellow")) +
  max.theme

# volume of propagation
df7 <- summarySE(propagation, measurevar="volume", groupvars=c("CellType"))

ggplot(data=df7, aes(x=CellType, y=volume, fill=CellType)) +
  geom_errorbar(aes(ymin=volume-se, ymax=volume+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("propagation volume") +
  scale_fill_manual(
    values=c("black", "red","blue", "green", "yellow")) +
  max.theme

# area of propagation
df8 <- summarySE(propagation, measurevar="area", groupvars=c("CellType"))

ggplot(data=df8, aes(x=CellType, y=area, fill=CellType)) +
  geom_errorbar(aes(ymin=area-se, ymax=area+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("propagation area") +
  scale_fill_manual(
    values=c("black", "red","blue", "green", "yellow")) +
  max.theme


#####
# propagation stats

# distance
dis.null = lmer(distance ~ (1|Animal) + (1|Spot2), propagation,REML=FALSE)
dis.model1 = lmer(distance~ CellType + (1|Animal) + (1|Spot2), propagation,REML=FALSE)
dis.anova <- anova(dis.null, dis.model1)
print(dis.anova)
# p values
dis.CellType <- lsmeans(dis.model1, pairwise ~ CellType, glhargs=list())
summary(dis.CellType)

# duration
dur.null = lmer(duration ~ (1|Animal) + (1|Spot2), propagation,REML=FALSE)
dur.model1 = lmer(duration ~ CellType + (1|Animal) + (1|Spot2), propagation,REML=FALSE)
dur.anova <- anova(dur.null, dur.model1)
print(dur.anova)
# p values
dur.CellType <- lsmeans(dur.model1, pairwise ~ CellType, glhargs=list())
summary(dur.CellType)

# propagation rate
PR.null = lmer(propRate ~ (1|Animal) + (1|Spot2), propagation,REML=FALSE)
PR.model1 = lmer(propRate ~ CellType + (1|Animal) + (1|Spot2), propagation,REML=FALSE)
PR.anova <- anova(PR.null, PR.model1)
print(PR.anova)
# p values
PR.CellType <- lsmeans(PR.model1, pairwise ~ CellType, glhargs=list())
summary(PR.CellType)

# volume
vol.null = lmer(volume ~ (1|Animal) + (1|Spot2), propagation,REML=FALSE)
vol.model1 = lmer(volume ~ CellType + (1|Animal) + (1|Spot2), propagation,REML=FALSE)
vol.anova <- anova(vol.null, vol.model1)
print(vol.anova)
# p values
vol.CellType <- lsmeans(vol.model1, pairwise ~ CellType, glhargs=list())
summary(vol.CellType)

# area
area.null = lmer(area ~ (1|Animal) + (1|Spot2), propagation,REML=FALSE)
area.model1 = lmer(area ~ CellType + (1|Animal) + (1|Spot2), propagation,REML=FALSE)
area.anova <- anova(area.null, area.model1)
print(area.anova)
# p values
area.CellType <- lsmeans(area.model1, pairwise ~ CellType, glhargs=list())
summary(area.CellType)


######
# Correlation of branches and soma

# load files
#corr<- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/Lck_GCaMP6f/Results/correlation.csv", header=TRUE, sep = ",")



#LinCorr
#xCorr
#Lag
#MinDistance
