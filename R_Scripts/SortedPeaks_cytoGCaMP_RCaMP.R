
library("lme4")
library("lmerTest")
library("lattice")
library("plyr")
library("ggplot2")
library("gplots")
library("lsmeans")
library("bear")
library("MASS")
library("multcomp")
library("reshape2")
library("data.table")
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
    legend.title=element_text(size=14, face="bold"))


########################
# relative "active ROI" between groups based on peak auc and frequency

# whole frame and automatic RCaMP ROI selection:
stim.all <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/Control_Peaks_3Conds.csv", header=TRUE, sep = ",")
auc.all <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/Control_TraceAUC_10sWindow_3Conds.csv", header=TRUE, sep = ",")

# stim onset at 5 sec
stim.all$peakTime<-stim.all$peakTime-5

# ROITypes
stim.all$ROIType= 0
stim.all$ROIType[grepl("r",stim.all$ROIname)]="Process"
stim.all$ROIType[grepl("E",stim.all$ROIname)]="Endfoot"
stim.all$ROIType[grepl("S",stim.all$ROIname)]="Soma"
stim.all$ROIType[grepl("N",stim.all$ROIname)]="Neuron"

stim.all$ROIType<- as.factor(stim.all$ROIType)


auc.all$ROIType= 0
auc.all$ROIType[grepl("r",auc.all$ROI)]="Process"
auc.all$ROIType[grepl("E",auc.all$ROI)]="Endfoot"
auc.all$ROIType[grepl("S",auc.all$ROI)]="Soma"
auc.all$ROIType[grepl("N",auc.all$ROI)]="Neuron"

auc.all$ROIType<- as.factor(auc.all$ROIType)


#unique ROI names

stim.all$ROIs_trial<-paste(stim.all$Animal, stim.all$Spot, stim.all$Trial,stim.all$ROIname, sep= "_")
auc.all$ROIs_trial<-paste(auc.all$Animal, auc.all$Spot, auc.all$Trial, auc.all$ROIname, sep= "_")

stim.all$ROIs<-paste(stim.all$Animal, stim.all$Spot, stim.all$ROIname, sep= "_")
auc.all$ROIs<-paste(auc.all$Animal, auc.all$Spot, auc.all$ROIname, sep= "_")

# histogram of peak times for stim

ggplot(stim.all, aes(x=peakTime, fill=Channel)) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("Distribution of peak times GCaMP and RCaMP")

#####
# adjust peak start time for stim onset
stim.all$peakStart2<-stim.all$peakStart-5

ggplot(stim.all, aes(x=peakStart2, fill=Channel)) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("Distribution of peak start times GCaMP and RCaMP")

#######
# remove matching process and soma ROIs
Animals = as.character(unique(stim.all$Animal))
stim.all2 <- data.frame()
for (iAnimal in 1:length(Animals))
{
  CurrentAnimal = Animals[iAnimal]
  AnimalSubset = subset(stim.all,Animal ==CurrentAnimal)
  Spots2 = as.character(unique(AnimalSubset$Spot))
  for (ii in 1:length(Spots2))
  {
    name = Spots2[ii]
    spotsubset = subset(AnimalSubset, Spot == name)
    # find ROIs with matches
    matchROIs = grepl("ROI*",spotsubset$Comparison)
    # find ROI number from matching ROIs
    comparisons <-as.character(spotsubset$Comparison[matchROIs])
    allROImatch<- unique(str_sub(comparisons,start=-2))
    if (length(allROImatch) ==0)
    { stim.all2<- rbind(stim.all2,spotsubset)
    } else {      
      for (iii in 1:length(allROImatch))
      {ROInum =as.numeric(allROImatch[iii])  
       ROInum =sprintf("%02d", ROInum) 
       ROItag = paste("ROI",ROInum,sep= "")
       # find the process ROIs
       ROIind= grepl(ROItag,spotsubset$ROI)
       spotsubset <- spotsubset[!ROIind,]     
      }
      stim.all2<- rbind(stim.all2,spotsubset)
    }
  }
}

#########
# trace AUC for first 10 s
df1A <- summarySE(auc.all, measurevar="AUC10s", groupvars=c("Condition","Channel"))
df1B <- summarySE(auc.all, measurevar="AUC10s", groupvars=c("Condition","ROIType"))

ggplot(data=df1A, aes(x=Channel, y=AUC10s, fill=Condition)) +
  geom_errorbar(aes(ymin=AUC10s-se, ymax=AUC10s+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("AUC for 10s of trace") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df1B, aes(x=ROIType, y=AUC10s, fill=Condition)) +
  geom_errorbar(aes(ymin=AUC10s-se, ymax=AUC10s+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("AUC for 10s of trace") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

####################
# aggregate peaks from soma, endfoot, neuron ROIs
struct.peaks<- subset(stim.all, ROIType!="Process")
struc.trials<- ddply(struct.peaks, c("Animal", "Spot", "Condition","Channel","ROI","Trial","ROIType"), summarise, 
                   PA_mean = mean(peak_auc), freq_mean = sum(numPeaks)/95*60,
                   Dur_mean = mean(Duration), Prom_mean = mean(prominence),
                   amp_mean = mean(amplitude), nEvents = sum(numPeaks), 
                   numActive = sum(ActivePeak))

autoRCaMP.trials$ActiveTrial<-0
peaks = autoRCaMP.trials$numActive>0
autoRCaMP.trials$ActiveTrial[peaks] <- 1

# aggregate ROI (no timepoints)
autoRCaMP.ROI<- ddply(autoRCaMP.trials, c("Animal", "Spot", "Condition","Channel","ROI","Layer","xLoc","yLoc","ROIArea"), summarise, 
                 PA_mean2 = mean(PA_mean), PA_SD = sd(PA_mean),
                 freq_mean2 = mean(freq_mean), freq_SD = sd(freq_mean),
                 Dur_mean2 = mean(Dur_mean), Dur_SD = sd(Dur_mean),
                 Prom_mean2 = mean(Prom_mean), Prom_SD = sd(Prom_mean),
                 amp_mean2 = mean(amp_mean), amp_SD = sd(amp_mean),TotEvents = sum(nEvents),
                 TotActive = sum(ActiveTrial), Ntrials= length(amp_mean))

autoRCaMP.ROI$frac.resp <- autoRCaMP.ROI$TotActive/autoRCaMP.ROI$Ntrials


# aggregate trials of all ROI (with timepoints)
wholeframe.trials<- ddply(wholeframe, c("Animal", "Spot", "Condition","Channel","ROI","Trial","Layer","xLoc","yLoc","ROIArea"), summarise, 
                         PA_mean = mean(peak_auc), freq_mean = sum(numPeaks)/95*60,
                         Dur_mean = mean(Duration), Prom_mean = mean(prominence),
                         amp_mean = mean(amplitude), nEvents = sum(numPeaks), 
                         numActive = sum(ActivePeak))

wholeframe.trials$ActiveTrial<-0
peaks = wholeframe.trials$numActive>0
wholeframe.trials$ActiveTrial[peaks] <- 1

# aggregate ROI (no timepoints)
wholeframe.ROI<- ddply(wholeframe.trials, c("Animal", "Spot", "Condition","Channel","ROI","Layer","xLoc","yLoc","ROIArea"), summarise, 
                      PA_mean2 = mean(PA_mean), PA_SD = sd(PA_mean),
                      freq_mean2 = mean(freq_mean), freq_SD = sd(freq_mean),
                      Dur_mean2 = mean(Dur_mean), Dur_SD = sd(Dur_mean),
                      Prom_mean2 = mean(Prom_mean), Prom_SD = sd(Prom_mean),
                      amp_mean2 = mean(amp_mean), amp_SD = sd(amp_mean),TotEvents = sum(nEvents),
                      TotActive = sum(ActiveTrial), Ntrials= length(amp_mean))

wholeframe.ROI$frac.resp <- wholeframe.ROI$TotActive/wholeframe.ROI$Ntrials

wholeframe.ROI$ROIType= 'wholeframe'
autoRCaMP.ROI$ROIType= 'autoRCaMP'

############
library(xlsx)
write.xlsx(wholeframe.ROI, "E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/wholeframeData.xlsx")
write.xlsx(autoRCaMP.ROI, "E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/autoRCaMPData.xlsx")

#################


#################
# remove unactive neurons
activeNeurons<- subset(HCRCaMP.ROI, Channel=="RCaMP" & Condition=="Stim" & PA_mean2>0)

activeNeurons2 <-subset(HCRCaMP.ROI, ROI %in% activeNeurons$ROI)

astrocytes <-subset(HCRCaMP.ROI, Channel=="cytoGCaMP")

HCRCaMP.ROI2<- rbind(activeNeurons2, astrocytes)

# drop shortstim
#short=HCRCaMP.ROI2$Condition=="shortstim"
#HCRCaMP.ROI2=HCRCaMP.ROI2[!short,]

#short=wholeframe.ROI$Condition=="shortstim"
#wholeframe.ROI=wholeframe.ROI[!short,]
#######
#all.ROI2 <- rbind(wholeframe.ROI, autoRCaMP.ROI, HCRCaMP.ROI2)

# mean percent response- i.e. the fraction of trials where a peak is detected in the first 30 sec
df1A <- summarySE(wholeframe.ROI, measurevar="frac.resp", groupvars=c("Condition","Channel"))
df1B <- summarySE(autoRCaMP.ROI, measurevar="frac.resp", groupvars=c("Condition","Channel"))
df1C <- summarySE(HCRCaMP.ROI2, measurevar="frac.resp", groupvars=c("Condition","Channel","ROIType"))

#####
ggplot(data=df1A, aes(x=ROIType, y=freq_mean2, fill=ROIType)) +
  geom_errorbar(aes(ymin=freq_mean2-se, ymax=freq_mean2+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Signals/min/ROI") +
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(data=df1A, aes(x=Channel, y=frac.resp, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=frac.resp-se, ymax=frac.resp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Condition") +
  ylab("frac.resp") +
  ggtitle("frac.resp for WF ROIs") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

ggplot(data=df1B, aes(x=Channel, y=frac.resp, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=frac.resp-se, ymax=frac.resp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Condition") +
  ylab("frac.resp") +
  ggtitle("frac.resp for autoRCaMP ROIs") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

ggplot(data=df1C, aes(x=interaction(Channel,ROIType), y=frac.resp, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=frac.resp-se, ymax=frac.resp+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Condition") +
  ylab("frac.resp") +
  ggtitle("frac.resp for HandClick ROIs") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

#############
#fraction of response
# Likelihood-ratio test 

#Whole frame
frac.resp.null = lmer(frac.resp ~ (1|Animal) + (1|Spot) + (1|ROI), wholeframe.ROI,REML=FALSE)
frac.resp.model1 = lmer(frac.resp ~ Condition + (1|Animal) + (1|Spot)+ (1|ROI), wholeframe.ROI,REML=FALSE)
frac.resp.model3A = lmer(frac.resp~ Condition+Channel + (1|Animal) + (1|Spot)+ (1|ROI), wholeframe.ROI,REML=FALSE)
frac.resp.model3B = lmer(frac.resp~ Condition*Channel + (1|Animal) + (1|Spot)+ (1|ROI), wholeframe.ROI,REML=FALSE)
frac.resp.anova <- anova(frac.resp.null, frac.resp.model1,frac.resp.model3A,frac.resp.model3B)
print(frac.resp.anova)
# p values
frac.resp.pv.stim <- lsmeans(frac.resp.model1, pairwise ~ Condition, glhargs=list())
summary(frac.resp.pv.stim)
frac.resp.pv.stim_type <- lsmeans(frac.resp.model3A, pairwise ~ Condition+Channel, glhargs=list())
summary(frac.resp.pv.stim_type)

# hand circled neurons and FLIKA
frac.resp.null = lmer(frac.resp ~ (1|Animal) + (1|Spot) + (1|ROI), HCRCaMP.ROI2,REML=FALSE)
frac.resp.model1 = lmer(frac.resp ~ Condition + (1|Animal) + (1|Spot)+ (1|ROI), HCRCaMP.ROI2,REML=FALSE)
frac.resp.model3A = lmer(frac.resp~ Condition+Channel + (1|Animal) + (1|Spot)+ (1|ROI), HCRCaMP.ROI2,REML=FALSE)
frac.resp.model3B = lmer(frac.resp~ Condition*Channel + (1|Animal) + (1|Spot)+ (1|ROI), HCRCaMP.ROI2,REML=FALSE)
frac.resp.anova <- anova(frac.resp.null, frac.resp.model1,frac.resp.model3A,frac.resp.model3B)
print(frac.resp.anova)
# p values
frac.resp.pv.stim <- lsmeans(frac.resp.model1, pairwise ~ Condition, glhargs=list())
summary(frac.resp.pv.stim)
frac.resp.pv.stim_type <- lsmeans(frac.resp.model3A, pairwise ~ Condition+Channel, glhargs=list())
summary(frac.resp.pv.stim_type)

######
df2A <- summarySE(wholeframe.ROI, measurevar="freq_mean2", groupvars=c("Condition","Channel"))
df2B <- summarySE(HCRCaMP.ROI2, measurevar="freq_mean2", groupvars=c("Condition","Channel"))

#####
ggplot(data=df2A, aes(x=Channel, y=freq_mean2, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=freq_mean2-se, ymax=freq_mean2+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Condition") +
  ylab("Events per min") +
  ggtitle("Frequency for WF ROIs") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

ggplot(data=df2B, aes(x=Channel, y=freq_mean2, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=freq_mean2-se, ymax=freq_mean2+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Condition") +
  ylab("Events per min") +
  ggtitle("Frequency for handclick ROIs") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

#####
#frequency
#Whole frame
freq.null = lmer(freq_mean2 ~ (1|Animal) + (1|Spot) + (1|ROI), wholeframe.ROI,REML=FALSE)
freq.model1 = lmer(freq_mean2 ~ Condition + (1|Animal) + (1|Spot)+ (1|ROI), wholeframe.ROI,REML=FALSE)
freq.model3A = lmer(freq_mean2 ~ Condition+Channel + (1|Animal) + (1|Spot)+ (1|ROI), wholeframe.ROI,REML=FALSE)
freq.model3B = lmer(freq_mean2 ~ Condition*Channel + (1|Animal) + (1|Spot)+ (1|ROI), wholeframe.ROI,REML=FALSE)
freq.anova <- anova(freq.null, freq.model1,freq.model3A,freq.model3B)
print(freq.anova)
# p values
freq.pv.stim <- lsmeans(freq.model1, pairwise ~ Condition, glhargs=list())
summary(freq.pv.stim)
freq.pv.stim_type <- lsmeans(freq.model3A, pairwise ~ Condition+Channel, glhargs=list())
summary(freq.pv.stim_type)

# hand circled neurons and FLIKA
freq.null = lmer(freq_mean2 ~ (1|Animal) + (1|Spot) + (1|ROI), HCRCaMP.ROI2,REML=FALSE)
freq.model1 = lmer(freq_mean2 ~ Condition + (1|Animal) + (1|Spot)+ (1|ROI), HCRCaMP.ROI2,REML=FALSE)
freq.model3A = lmer(freq_mean2 ~ Condition+Channel + (1|Animal) + (1|Spot)+ (1|ROI), HCRCaMP.ROI2,REML=FALSE)
freq.model3B = lmer(freq_mean2 ~ Condition*Channel + (1|Animal) + (1|Spot)+ (1|ROI), HCRCaMP.ROI2,REML=FALSE)
freq.anova <- anova(freq.null, freq.model1,freq.model3A,freq.model3B)
print(freq.anova)
# p values
freq.pv.stim <- lsmeans(freq.model1, pairwise ~ Condition, glhargs=list())
summary(freq.pv.stim)
freq.pv.stim_type <- lsmeans(freq.model3B, pairwise ~ Condition*Channel, glhargs=list())
summary(freq.pv.stim_type)


#############
df4A<- summarySE(wholeframe.ROI, measurevar="amp_mean2", groupvars=c("Condition","Channel"))
df4B<- summarySE(HCRCaMP.ROI2, measurevar="amp_mean2", groupvars=c("Condition","Channel"))

##########
ggplot(data=df4A, aes(x=Channel, y=amp_mean2, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amp_mean2-se, ymax=amp_mean2+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Channel") +
  ylab("amp") +
  ggtitle("amp for WF ROIs") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

ggplot(data=df4B, aes(x=Channel, y=amp_mean2, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amp_mean2-se, ymax=amp_mean2+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Channel") +
  ylab("amp") +
  ggtitle("amp for handclicked ROIs") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme


#############
#amplitude
#Whole frame
amp.null = lmer(amp_mean2 ~ (1|Animal) + (1|Spot) + (1|ROI), wholeframe.ROI,REML=FALSE)
amp.model1 = lmer(amp_mean2 ~ Condition + (1|Animal) + (1|Spot)+ (1|ROI), wholeframe.ROI,REML=FALSE)
amp.model3A = lmer(amp_mean2 ~ Condition+Channel + (1|Animal) + (1|Spot)+ (1|ROI), wholeframe.ROI,REML=FALSE)
amp.model3B = lmer(amp_mean2 ~ Condition*Channel + (1|Animal) + (1|Spot)+ (1|ROI), wholeframe.ROI,REML=FALSE)
amp.anova <- anova(amp.null, amp.model1,amp.model3A,amp.model3B)
print(amp.anova)
# p values
amp.pv.stim <- lsmeans(amp.model1, pairwise ~ Condition, glhargs=list())
summary(amp.pv.stim)
amp.pv.stim_type <- lsmeans(amp.model3A, pairwise ~ Condition+Channel, glhargs=list())
summary(amp.pv.stim_type)

# hand circled neurons and FLIKA
amp.null = lmer(amp_mean2 ~ (1|Animal) + (1|Spot) + (1|ROI), HCRCaMP.ROI2,REML=FALSE)
amp.model1 = lmer(amp_mean2 ~ Condition + (1|Animal) + (1|Spot)+ (1|ROI), HCRCaMP.ROI2,REML=FALSE)
amp.model3A = lmer(amp_mean2 ~ Condition+Channel + (1|Animal) + (1|Spot)+ (1|ROI), HCRCaMP.ROI2,REML=FALSE)
amp.model3B = lmer(amp_mean2 ~ Condition*Channel + (1|Animal) + (1|Spot)+ (1|ROI), HCRCaMP.ROI2,REML=FALSE)
amp.anova <- anova(amp.null, amp.model1,amp.model3A,amp.model3B)
print(amp.anova)
# p values
amp.pv.stim <- lsmeans(amp.model1, pairwise ~ Condition, glhargs=list())
summary(amp.pv.stim)
amp.pv.stim_type <- lsmeans(amp.model3B, pairwise ~ Condition*Channel, glhargs=list())
summary(amp.pv.stim_type)

############
df3A <- summarySE(wholeframe.ROI, measurevar="Dur_mean2", groupvars=c("Condition","Channel"))
df3B <- summarySE(HCRCaMP.ROI2, measurevar="Dur_mean2", groupvars=c("Condition","Channel"))

###############
ggplot(data=df3A, aes(x=Channel, y=Dur_mean2, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Dur_mean2-se, ymax=Dur_mean2+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Channel") +
  ylab("Duration") +
  ggtitle("Duration for all ROIs") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

ggplot(data=df3B, aes(x=Channel, y=Dur_mean2, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Dur_mean2-se, ymax=Dur_mean2+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Channel") +
  ylab("Duration") +
  ggtitle("Duration for all ROIs") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

#amplitude
#Whole frame
dur.null = lmer(Dur_mean2 ~ (1|Animal) + (1|Spot) + (1|ROI), wholeframe.ROI,REML=FALSE)
dur.model1 = lmer(Dur_mean2 ~ Condition + (1|Animal) + (1|Spot)+ (1|ROI), wholeframe.ROI,REML=FALSE)
dur.model3A = lmer(Dur_mean2 ~ Condition+Channel + (1|Animal) + (1|Spot)+ (1|ROI), wholeframe.ROI,REML=FALSE)
dur.model3B = lmer(Dur_mean2 ~ Condition*Channel + (1|Animal) + (1|Spot)+ (1|ROI), wholeframe.ROI,REML=FALSE)
dur.anova <- anova(dur.null, dur.model1,dur.model3A,dur.model3B)
print(dur.anova)
# p values
dur.pv.stim <- lsmeans(dur.model1, pairwise ~ Condition, glhargs=list())
summary(dur.pv.stim)
dur.pv.stim_type <- lsmeans(dur.model3A, pairwise ~ Condition+Channel, glhargs=list())
summary(dur.pv.stim_type)

# hand circled neurons and FLIKA
dur.null = lmer(Dur_mean2 ~ (1|Animal) + (1|Spot) + (1|ROI), HCRCaMP.ROI2,REML=FALSE)
dur.model1 = lmer(Dur_mean2~ Condition + (1|Animal) + (1|Spot)+ (1|ROI), HCRCaMP.ROI2,REML=FALSE)
dur.model3A = lmer(Dur_mean2~ Condition+Channel + (1|Animal) + (1|Spot)+ (1|ROI), HCRCaMP.ROI2,REML=FALSE)
dur.model3B = lmer(Dur_mean2~ Condition*Channel + (1|Animal) + (1|Spot)+ (1|ROI), HCRCaMP.ROI2,REML=FALSE)
dur.anova <- anova(dur.null, dur.model1,dur.model3A,dur.model3B)
print(dur.anova)
# p values
dur.pv.stim <- lsmeans(dur.model1, pairwise ~ Condition, glhargs=list())
summary(dur.pv.stim)
dur.pv.stim_type <- lsmeans(dur.model3B, pairwise ~ Condition*Channel, glhargs=list())
summary(dur.pv.stim_type)

######################

# considering peaks
# when do AC peaks come after neuronal peak?
# bin the number of AC peaks per time point after stim?

peaks<-subset(stim.all2, peakTime>-5)
# drop shortstim
#short=peaks$Condition=="shortstim"
#peaks=peaks[!short,]

ggplot(peaks, aes(x=peakTime, fill=Condition)) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("peakTimes- All comparisons")


peaks.stim<-subset(peaks, Condition=="Stim")

ggplot(peaks.stim, aes(x=peakTime, fill=Channel)) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("peakTimes- All comparisons")


# short stim & long stim
N_short <- subset(stim.all2, peakTime>0 & Channel=="RCaMP" & Condition=="shortstim")
N_long <- subset(stim.all2, peakTime>0 & Channel=="RCaMP" & Condition=="Stim")
N_peaks <- rbind(N_short, N_long)
N_peaks$comb<-paste(N_peaks$Animal, N_peaks$Spot, N_peaks$Trial, sep= "_")
N_peaks$comb2<-paste(N_peaks$Animal, N_peaks$Spot, N_peaks$Trial, N_peaks$Condition, sep= "_")
N_peaks$onset <- N_peaks$peakTime- N_peaks$halfWidth
N_stim_peaks <- subset(N_peaks, peakTime<10)

AC_short <- subset(stim.all2, Channel=="cytoGCaMP" & Condition =="shortstim")
AC_long <- subset(stim.all2, Channel=="cytoGCaMP" & Condition =="Stim")
AC_peaks <- rbind(AC_short, AC_long)
AC_peaks$comb<-paste(AC_peaks$Animal, AC_peaks$Spot, AC_peaks$Trial, sep= "_")
AC_peaks$comb2<-paste(AC_peaks$Animal, AC_peaks$Spot, AC_peaks$Trial,AC_peaks$Condition, sep= "_")
AC_peaks$onset <- AC_peaks$peakTime- AC_peaks$halfWidth
AC_stim_peaks <- subset(AC_peaks, peakTime<15)
AC_stim_peaks <- subset(AC_stim_peaks , peakTime>0)

# mean peak time in the
stim_peaks <- rbind(N_stim_peaks, AC_stim_peaks)

df9A <- summarySE(stim_peaks, measurevar="peakTime", groupvars=c("Condition","Channel"))
df10A <- summarySE(stim_peaks, measurevar="onset", groupvars=c("Condition","Channel"))

ggplot(data=df9A, aes(x=interaction(Channel,Condition), y=peakTime, fill=Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("Channel") +
  ylab("Mean Peak Time (s)") +
  ggtitle("timepeaks for stim") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

ggplot(data=df10A, aes(x=interaction(Channel,Condition), y=onset, fill=Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=onset-se, ymax=onset+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("Channel") +
  ylab("onset Time (t1/2)") +
  ggtitle("timepeaks for stim") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

# hand circled neurons and FLIKA
pT.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|comb), stim_peaks,REML=FALSE)
pT.model1 = lmer(peakTime ~ Channel + (1|Animal) + (1|Spot)+ (1|comb), stim_peaks,REML=FALSE)
pT.model2 = lmer(peakTime ~ Channel+Condition + (1|Animal) + (1|Spot)+ (1|comb), stim_peaks,REML=FALSE)
pT.model3 = lmer(peakTime ~ Channel*Condition + (1|Animal) + (1|Spot)+ (1|comb), stim_peaks,REML=FALSE)
pT.anova <- anova(pT.null, pT.model1,pT.model2,pT.model3)
print(pT.anova)
# p values
pT.pv.stim1 <- lsmeans(pT.model1, pairwise ~ Channel, glhargs=list())
summary(pT.pv.stim1)
pT.pv.stim2 <- lsmeans(pT.model3, pairwise ~ Channel*Condition, glhargs=list())
summary(pT.pv.stim2)

# half maximum time
HM.null = lmer(onset ~ (1|Animal) + (1|Spot) + (1|comb), stim_peaks,REML=FALSE)
HM.model1 = lmer(onset ~ Channel + (1|Animal) + (1|Spot)+ (1|comb), stim_peaks,REML=FALSE)
HM.model2 = lmer(onset ~ Channel+Condition + (1|Animal) + (1|Spot)+ (1|comb), stim_peaks,REML=FALSE)
HM.model3 = lmer(onset ~ Channel*Condition + (1|Animal) + (1|Spot)+ (1|comb), stim_peaks,REML=FALSE)
HM.anova <- anova(HM.null, HM.model1,HM.model2,HM.model3)
print(HM.anova)
# p values
HM.pv.stim1 <- lsmeans(HM.model1, pairwise ~ Channel, glhargs=list())
summary(HM.pv.stim1)
HM.pv.stim2 <- lsmeans(HM.model3, pairwise ~ Channel*Condition, glhargs=list())
summary(HM.pv.stim2)

###########
# BUT! We should correlate specific AC peaks with specific neuron peaks from the SAME TRIAL

# count the number of AC peaks that occur in the same trials and neuronal peaks
# bin by timepoint

ActiveNeuronTrials= unique(N_stim_peaks$comb2)

NeuronTrials = data.frame()
Stim_pooled = data.frame()
Dis= data.frame()
for (ipeaks in 1:length(ActiveNeuronTrials))
{
  CurrentTrial= ActiveNeuronTrials[ipeaks]
  ACsubset = subset(AC_peaks, comb2 == CurrentTrial)
  ACsubset2 = subset(AC_stim_peaks, comb2 == CurrentTrial)
  Neuronsub = subset(N_stim_peaks, comb2==CurrentTrial)
  # subset AC data based on peak time
  AC0 = subset(ACsubset, peakTime<0)
  AC10 = subset(ACsubset, peakTime<10 & peakTime>=0)
  AC20 = subset(ACsubset, peakTime<20 & peakTime>=10)
  AC40 = subset(ACsubset, peakTime<40 & peakTime>=20)
  AC60 = subset(ACsubset, peakTime<60 & peakTime>=40)
  
  # count AC peaks in each bin
  Data= Neuronsub[1,]
  Data$AC0 = sum(length(AC0$peakTime))
  Data$AC10 = sum(length(AC10$peakTime))
  Data$AC20 = sum(length(AC20$peakTime))
  Data$AC40 = sum(length(AC40$peakTime))
  Data$AC60 = sum(length(AC60$peakTime))
  NeuronTrials = rbind(NeuronTrials, Data)
  
  # calculate mean peakTime and onset for neuronal signals and astrocyte signals in each trial
  if (length(ACsubset2[,1])>0)
  {
    Data_pooled1<- ddply(Neuronsub, c("Animal", "Spot", "Condition","comb2"), summarise, 
                         Channel="Neuron", PT_mean = mean(peakTime), PT_SD = sd(peakTime),
                         onset_mean = mean(onset), onset_SD = sd(onset),
                         amp_mean = mean(amplitude), amp_SD = sd(amplitude),
                         Peaknum = sum(length(amplitude)))
    Data_pooled2<- ddply(ACsubset2, c("Animal", "Spot", "Condition","comb2"), summarise, 
                         Channel="AC", PT_mean = mean(peakTime), PT_SD = sd(peakTime),
                         onset_mean = mean(onset), onset_SD = sd(onset),
                         amp_mean = mean(amplitude), amp_SD = sd(amplitude),
                         Peaknum = sum(length(amplitude)))
    #Data_pooled1$response = Neuron
    #Data_pooled2$response = AC
    Data_pooled<- rbind(Data_pooled1,Data_pooled2)
    Stim_pooled<- rbind(Data_pooled, Stim_pooled)
  }
}

df11A <- summarySE(Stim_pooled, measurevar="PT_mean", groupvars=c("Condition","Channel"))
df12A <- summarySE(Stim_pooled, measurevar="onset", groupvars=c("Condition","Channel"))

#Interaction plots
library(nlme)
RCaMP = subset(Stim_pooled, Channel=="Neuron")
GCaMP = subset(Stim_pooled, Channel=="AC")
for (iSpot in 1:length(RCaMP[,1]))
  {
  Data1=RCaMP[iSpot,]
  Data2=GCaMP[iSpot,]
  Data= rbind(Data1,Data2)
  interaction.plot(Data$Channel, Data$Condition, Data$PT_mean, ylab = "Peak Time", xlab = "Channel")
}



#### ?????????? Do in matlab?

# calculate the distance between each neuron with a signal and all responding AC ROIs
Neuron_unique = length(unique(Neuronsub$ROI))
AC_unique = length(unique(ACsubset2$ROI))
for (iROI in 1:Neuron_unique)
{
  for (xROI in 1:AC_unique)
  {
    Dis
    Dis$comb2[iROI]=Neuronsub$comb2[iROI]
    Dis$ROI_1 = Neuronsub$ROI[iROI]
  }
}
