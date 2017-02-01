
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
    legend.title=element_text(size=14, face="bold"))


###########
# NOTES

# Did DSP4 injection work?  Have a look at the raw data
# Can we see differences in the response to stim or spontaneous activity?
# Check immediately following stim but also in the period after

# remember- DSP4 images were from the same spots on different days


# sort out the trials where there is a neuronal response and only consider those
# number of trials with a response with DSP4 and not

# some process ROIs are very large


#is it a problem that I have mutliple measurements for soma, EF, and neurons from different trials, but only one number for processes?

# Trials with arousal? Before and after changes?

# Can I find astrocyte signals that correlate with spontaneous or stimulus evoked neuronal signals?



# time correlation of neuronal and astrocyte signals?
# consider each trial individually


# exclude mouse RG10 because it twitched too much?


########################
# relative "active ROI" between groups based on peak auc and frequency

# whole frame and automatic RCaMP ROI selection:
peaks.control <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/S&LStim_cGC&RC_01_30_2017.csv", header=TRUE, sep = ",")
auc.control <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/awake_cGC&RC_AUC_01_2017.csv", header=TRUE, sep = ",")
peaks.DSP4 <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cGC&RC_01_2017.csv", header=TRUE, sep = ",")
auc.DSP4 <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cGC&RC_AUC_01_2017.csv", header=TRUE, sep = ",")

#peaks.control <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/awake_cGC&RC_01_2017.csv", header=TRUE, sep = ",")
#auc.control <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/awake_cGC&RC_AUC_01_2017.csv", header=TRUE, sep = ",")
#peaks.DSP4 <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cGC&RC_01_2017.csv", header=TRUE, sep = ",")
#auc.DSP4 <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/DSP4_cGC&RC_AUC_01_2017.csv", header=TRUE, sep = ",")

lsm.options(pbkrtest.limit = 100000)

# treatment
peaks.control$treatment<-"Control"
auc.control$treatment<-"Control"
peaks.DSP4$treatment<-"DSP4"
auc.DSP4$treatment<-"DSP4"

stim.all<-rbind(peaks.control, peaks.DSP4)
auc.all<-rbind(auc.control, auc.DSP4)


######
# exclude DSP4 data FOR NOW

stim.all<-subset(stim.all, treatment=="Control")
auc.all<-subset(auc.all, treatment=="Control")

stim.all$treatment<-as.factor(stim.all$treatment)
auc.all$treatment<-as.factor(auc.all$treatment)

# stim onset at 5 sec
stim.all$peakTime<-stim.all$peakTime-5
stim.all$peakStart<-stim.all$peakStart-5
stim.all$peakStartHalf<-stim.all$peakStartHalf-5

#duration
stim.all$Duration<-stim.all$halfWidth*2

# ROITypes
stim.all$ROIType= 0
stim.all$ROIType[grepl("r",stim.all$ROIname)]="Process"
stim.all$ROIType[grepl("E",stim.all$ROIname)]="Endfoot"
stim.all$ROIType[grepl("S",stim.all$ROIname)]="Soma"
stim.all$ROIType[grepl("N",stim.all$ROIname)]="Neuron"
stim.all$ROIType[grepl("np",stim.all$ROIname)]="Neuropil"

stim.all$ROIType<- as.factor(stim.all$ROIType)


auc.all$ROIType= 0
auc.all$ROIType[grepl("r",auc.all$ROI)]="Process"
auc.all$ROIType[grepl("E",auc.all$ROI)]="Endfoot"
auc.all$ROIType[grepl("S",auc.all$ROI)]="Soma"
auc.all$ROIType[grepl("N",auc.all$ROI)]="Neuron"
auc.all$ROIType[grepl("np",auc.all$ROIname)]="Neuropil"

auc.all$ROIType<- as.factor(auc.all$ROIType)


#unique ROI names
stim.all$ROIs_trial<-paste(stim.all$Animal, stim.all$Spot, stim.all$Trial,stim.all$ROIname, sep= "_")
auc.all$ROIs_trial<-paste(auc.all$Animal, auc.all$Spot, auc.all$Trial, auc.all$ROI, sep= "_")

stim.all$ROIs<-paste(stim.all$Animal, stim.all$Spot, stim.all$ROIname, sep= "_")
auc.all$ROIs<-paste(auc.all$Animal, auc.all$Spot, auc.all$ROI, sep= "_")

stim.all$trials<-paste(stim.all$Animal, stim.all$Spot, stim.all$Trial, sep= "_")
auc.all$trials<-paste(auc.all$Animal, auc.all$Spot, auc.all$Trial, sep= "_")


# exclude RG10 cause it moves too much
#stim.all<- subset(stim.all, Animal!="RG10")



# remove matching astrocyte process and soma ROIs
Overlap= stim.all$overlap!=0

OverlapROIs<-unique(stim.all$ROIs_trial[Overlap])

stim.all2<-stim.all[!Overlap,]

##########
# remove ROIs from AUC data frame
#AUC_overlap<- subset(auc.all, ROIs_trial !%in% OverlapROIs)

# make rows with no peak zero instead of NaN
#test3 = is.na(stim.all2$amplitude)
#stim.all2$peakAUC[test3] = 0
#stim.all2$prominence[test3] = 0
#stim.all2$amplitude[test3] = 0
#stim.all2$peakTime[test3] = -5
#stim.all2$halfWidth[test3] = 0
#stim.all2$fullWidth[test3] = 0
#stim.all2$Duration[test3] = 0
#stim.all2$numPeaks[test3] = 0



#######
# histogram of peak times for stim
longstim<-subset(stim.all2, Condition=="Stim")
ggplot(longstim, aes(x=peakTime, fill=Channel)) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("long stim")

shortstim<-subset(stim.all2, Condition=="shortstim")
ggplot(shortstim, aes(x=peakTime, fill=Channel)) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("short stim")

nostim<-subset(stim.all2, Condition=="Nostim")
ggplot(nostim, aes(x=peakTime, fill=Channel)) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("No stim")

RCaMP<- subset(stim.all2, Channel=="RCaMP")
ggplot(RCaMP, aes(x=peakTime, fill=Condition)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("RCaMP")

GCaMP<- subset(stim.all2, Channel=="GCaMP")
ggplot(GCaMP, aes(x=peakTime, fill=Condition)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("GCaMP")

#considering treatment
RCaMP.longstim<- subset(RCaMP, Condition=="Stim")
ggplot(RCaMP.longstim, aes(x=peakTime, fill=treatment)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("RCaMP long stim")

RCaMP.shortstim<- subset(RCaMP, Condition=="shortstim")
ggplot(RCaMP.shortstim, aes(x=peakTime, fill=treatment)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("RCaMP short stim")

#considering treatment
GCaMP.longstim<- subset(GCaMP, Condition=="Stim")
ggplot(GCaMP.longstim, aes(x=peakTime, fill=treatment)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("GCaMP long stim")

GCaMP.shortstim<- subset(GCaMP, Condition=="shortstim")
ggplot(GCaMP.shortstim, aes(x=peakTime, fill=treatment)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("GCaMP short stim")
###########

# remove data with really large prominences
stim.all3<- subset(stim.all2, prominence<15)

ggplot(stim.all3, aes(x=prominence, fill=treatment)) + geom_histogram(binwidth=1, position="dodge")

# outliers in each GCaMP or RCaMP group

outlierKD <- function(dt, var) {
  var_name <- eval(substitute(var),eval(dt))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow=c(2, 2), oma=c(0,0,3,0))
  boxplot(var_name, main="With outliers")
  hist(var_name, main="With outliers", xlab=NA, ylab=NA)
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main="Without outliers")
  hist(var_name, main="Without outliers", xlab=NA, ylab=NA)
  title("Outlier Check", outer=TRUE)
  na2 <- sum(is.na(var_name))
  cat("Outliers identified:", na2 - na1, "n")
  cat("Propotion (%) of outliers:", round((na2 - na1) / sum(!is.na(var_name))*100, 1), "n")
  cat("Mean of the outliers:", round(mo, 2), "n")
  m2 <- mean(var_name, na.rm = T)
  cat("Mean without removing outliers:", round(m1, 2), "n")
  cat("Mean if we remove outliers:", round(m2, 2), "n")
  response <- readline(prompt="Do you want to remove outliers and to replace with NA? [yes/no]: ")
  if(response == "y" | response == "yes"){
    dt[as.character(substitute(var))] <- invisible(var_name)
    assign(as.character(as.list(match.call())$dt), dt, envir = .GlobalEnv)
    cat("Outliers successfully removed", "n")
    return(invisible(dt))
  } else{
    cat("Nothing changed", "n")
    return(invisible(var_name))
  }
}

#source("http://goo.gl/UUyEzD")
outlierKD(GCaMP, prominence)
yes

outlierKD(RCaMP, prominence)
yes

stim.all3<-rbind(GCaMP, RCaMP)

#remove ROIs with NaNs
#stim.all3 = stim.all3[complete.cases(stim.all3$prominence),]

#remove ROIs with no name
#stim.all3$ROIname <- as.character(stim.all3$ROIname)
#stim.all3 = stim.all3[(stim.all3$ROIname!=""),]
#stim.all3$ROIname <- as.factor(stim.all3$ROIname)




########
# consider only the trials where the neurons responded to stimulation

# what is a neuronal response to stimulation??
# I defined it as a peak within stimulus onset to 1 sec after stimulus stop
# PLUS- peak duration must be close to this window
# long stim= peak between 0 and 9 sec, duration < 11 s
# short stim= peak between 0 and 2 sec, duration < 3 s

# find responding neurons
responding.neurons_long<- subset(longstim, peakTime>0 & peakTime<8 & Duration<11 & ROIType=="Neuron")
responding.neurons_short<- subset(shortstim, peakTime>0 & peakTime<2 & Duration<3 & ROIType=="Neuron")

responding.trials_long<-unique(responding.neurons_long$trials)  # 225 trials of 265- 85% of trials
responding.trials_short<-unique(responding.neurons_short$trials) # 102 trials of 272- 38% of trials

longstim.responding<-subset(longstim, trials %in% responding.trials_long)
shortstim.responding<-subset(shortstim, trials %in% responding.trials_short)

# distribution of peaktimes
ggplot(longstim.responding, aes(x=peakTime, fill=interaction(Channel,treatment))) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("long stim responding")

ggplot(shortstim.responding, aes(x=peakTime, fill=interaction(Channel,treatment))) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("short stim responding")

#####
#neuronal population sort into high, mid, low groups

# consider all peaks with a time near 10 s
neurons_longstim.mean<- ddply(responding.neurons_long, c("Animal", "Spot", "treatment", "ROIs"), summarise, 
                              PA_mean = mean(peakAUC), nEvents = length(peakAUC),
                              Dur_mean = mean(Duration), Prom_mean = mean(prominence),
                              amp_mean = mean(amplitude), HalfDur = mean(halfWidth),
                              peakT_mean = mean(peakTime),peakHalf_mean= mean(peakStartHalf))

ggplot(neurons_longstim.mean, aes(x=Prom_mean)) + geom_histogram(binwidth=0.05, position="dodge") +
  ggtitle("all neurons during long stim (10 s window)")


Prominence_percentiles<-quantile(neurons_longstim.mean$Prom_mean, prob = seq(0, 1, length = 21), type = 5)

highresponding<-subset(neurons_longstim.mean, Prom_mean>Prominence_percentiles[20])
midresponding<-subset(neurons_longstim.mean, Prom_mean<=Prominence_percentiles[20]&Prom_mean>=Prominence_percentiles[11])
lowresponding<-subset(neurons_longstim.mean, Prom_mean<Prominence_percentiles[11])

longstim.responding$NeuronGroup <- 0
highresponders=unique(highresponding$ROIs)
midresponders=unique(midresponding$ROIs)
lowresponders=unique(lowresponding$ROIs)

longstim.responding$NeuronGroup[longstim.responding$ROIs %in% highresponders]<-"high"
longstim.responding$NeuronGroup[longstim.responding$ROIs %in% midresponders]<-"mid"
longstim.responding$NeuronGroup[longstim.responding$ROIs %in% lowresponders]<-"low"


# consider all peaks with a time near 10 s
neurons_shortstim.mean<- ddply(responding.neurons_short, c("Animal", "Spot", "treatment", "ROIs"), summarise, 
                              PA_mean = mean(peakAUC), nEvents = length(peakAUC),
                              Dur_mean = mean(Duration), Prom_mean = mean(prominence),
                              amp_mean = mean(amplitude), HalfDur = mean(halfWidth),
                              peakT_mean = mean(peakTime),peakHalf_mean= mean(peakStartHalf))

ggplot(neurons_shortstim.mean, aes(x=Prom_mean)) + geom_histogram(binwidth=0.05, position="dodge") +
  ggtitle("all neurons during short stim (10 s window)")


Prominence_percentiles<-quantile(neurons_shortstim.mean$Prom_mean, prob = seq(0, 1, length = 21), type = 5)

highresponding<-subset(neurons_shortstim.mean, Prom_mean>Prominence_percentiles[20])
midresponding<-subset(neurons_shortstim.mean, Prom_mean<=Prominence_percentiles[20]&Prom_mean>=Prominence_percentiles[11])
lowresponding<-subset(neurons_shortstim.mean, Prom_mean<Prominence_percentiles[11])

shortstim.responding$NeuronGroup <- 0
highresponders=unique(highresponding$ROIs)
midresponders=unique(midresponding$ROIs)
lowresponders=unique(lowresponding$ROIs)

shortstim.responding$NeuronGroup[shortstim.responding$ROIs %in% highresponders]<-"high"
shortstim.responding$NeuronGroup[shortstim.responding$ROIs %in% midresponders]<-"mid"
shortstim.responding$NeuronGroup[shortstim.responding$ROIs %in% lowresponders]<-"low"

#########
# LONG STIM (90Hz, 8sec)

# identify active ROIs
longstim.responding$ActivePeak <- 0

#responding neurons
farpeaks1 <- longstim.responding$peakTime>0 & longstim.responding$peakTime<8 & longstim.responding$Duration<11 & longstim.responding$ROIType=="Neuron"
longstim.responding$ActivePeak[farpeaks1] <- 1 

#responding astrocytes
farpeaks2 <- longstim.responding$peakTime>0 & longstim.responding$peakTime<20 & longstim.responding$ROIType!="Neuron"
longstim.responding$ActivePeak[farpeaks2] <- 1 

# pull out only the peaks that occur around the stimulation
longstim.stimwindow<- subset(longstim.responding, peakTime>=0 & peakTime<=20)



ggplot(longstim.stimwindow, aes(x=peakTime, fill=interaction(Channel,treatment))) + geom_histogram(binwidth=0.5, position="dodge") +
  ggtitle("long stim responding")

longstim.after<- subset(longstim.responding, peakTime>20 & peakTime<80)

ggplot(longstim.after, aes(x=peakTime, fill=interaction(Channel,treatment))) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("long stim after")

#####
library(xlsx)
respondingNeurons_Astrocytes=subset(longstim.responding, ActivePeak==1)
write.xlsx(respondingNeurons_Astrocytes, "E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/respondingROIs_longstim.xlsx")


#########
# proportion of ROIs that respond per trial
activeROIs1<- ddply(longstim.responding, c("Animal", "Spot", "Channel","trials","ROIType","treatment", "ROIs_trial"), summarise, 
                    nEvents = sum(ActivePeak))

activeROIs1$ROIActive<-0
active1 <- activeROIs1$nEvents>0
activeROIs1$ROIActive[active1] <- 1 

longstim.propActive<- ddply(activeROIs1, c("Animal", "Spot", "Channel","trials","ROIType","treatment"), summarise, 
                            nSignals = sum(nEvents), nROIs= length(unique(ROIs_trial)),
                            nActiveROIs = sum(ROIActive))
longstim.propActive$propActive<-longstim.propActive$nActiveROIs/longstim.propActive$nROIs

df2A1<-summarySE(longstim.propActive, measurevar="propActive", groupvars=c("ROIType","treatment"))
df2A2<-summarySE(longstim.propActive, measurevar="propActive", groupvars=c("Channel","treatment"))

ggplot(longstim.propActive, aes(x=propActive, fill=Channel)) + geom_histogram(binwidth=0.05, position="dodge") +
  ggtitle("long stim prop active ")

df2A1$ROIType <- factor(df2A1$ROIType , levels = c("Neuron","Endfoot","Soma","Process"))
ggplot(data=df2A1, aes(x=ROIType, y=propActive, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=propActive-se, ymax=propActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("ROIType") +
  ylab("Proportion of Responding ROIs") +
  ggtitle("Responding ROIs for long stim") 


ggplot(data=df2A2, aes(x=Channel, y=propActive, fill=Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=propActive-se, ymax=propActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Channel") +
  ylab("Proportion of Responding ROIs") +
  ggtitle("Responding ROIs for long stim") 





#########
# SHORT STIM (90Hz, 1sec)
shortstim.responding$ActivePeak <- 0
farpeaks1 <- shortstim.responding$peakTime>0 & shortstim.responding$peakTime<2 & shortstim.responding$Duration<3 & shortstim.responding$ROIType=="Neuron"
shortstim.responding$ActivePeak[farpeaks1] <- 1 

farpeaks2 <- shortstim.responding$peakTime>0 & shortstim.responding$peakTime<20 & shortstim.responding$ROIType!="Neuron"
shortstim.responding$ActivePeak[farpeaks2] <- 1 

# pull out only the peaks that occur around the stimulation
shortstim.stimwindow<- subset(shortstim.responding, peakTime>=0 & peakTime<=10)

ggplot(shortstim.stimwindow, aes(x=peakTime, fill=interaction(Channel,treatment))) + geom_histogram(binwidth=0.5, position="dodge") +
  ggtitle("short stim responding")

ggplot(shortstim.stimwindow, aes(x=peakTime, y=prominence, color=interaction(Channel,treatment))) + geom_point() +
  ggtitle("short stim responding")

shortstim.after<- subset(shortstim.responding, peakTime>10 & peakTime<80)

ggplot(shortstim.after, aes(x=peakTime, fill=interaction(Channel,treatment))) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("short stim after")

respondingNeurons_AstrocytesShort=subset(shortstim.responding, ActivePeak==1)
write.xlsx(respondingNeurons_Astrocytes, "E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/respondingROIs_shortstim.xlsx")



activeROIs2<- ddply(shortstim.responding, c("Animal", "Spot", "Channel","trials","ROIType","treatment", "ROIs_trial"), summarise, 
                    nEvents = sum(ActivePeak))

activeROIs2$ROIActive<-0
active2 <- activeROIs2$nEvents>0
activeROIs2$ROIActive[active2] <- 1 

shortstim.propActive<- ddply(activeROIs2, c("Animal", "Spot", "Channel","trials","ROIType","treatment"), summarise, 
                             nSignals = sum(nEvents), nROIs= length(unique(ROIs_trial)),
                             nActiveROIs = sum(ROIActive))
shortstim.propActive$propActive<-shortstim.propActive$nActiveROIs/shortstim.propActive$nROIs


df2B1<-summarySE(shortstim.propActive, measurevar="propActive", groupvars=c("ROIType","treatment"))
df2B2<-summarySE(shortstim.propActive, measurevar="propActive", groupvars=c("Channel","treatment"))

ggplot(shortstim.propActive, aes(x=propActive, fill=Channel)) + geom_histogram(binwidth=0.05, position="dodge") +
  ggtitle("short stim prop active ")


df2B1$ROIType <- factor(df2B1$ROIType , levels = c("Neuron","Endfoot","Soma","Process"))
ggplot(data=df2B1, aes(x=ROIType, y=propActive, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=propActive-se, ymax=propActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("ROIType") +
  ylab("Proportion of Responding ROIs") +
  ggtitle("Responding ROIs for short stim") 


ggplot(data=df2B2, aes(x=Channel, y=propActive, fill=Channel)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=propActive-se, ymax=propActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Channel") +
  ylab("Proportion of Responding ROIs") +
  ggtitle("Responding ROIs for short stim") 



##########

## Signal characteristics for all peaks together
respondingAUC.long<- subset(auc.all, trials %in% responding.trials_long & Condition=="Stim")
respondingAUC.short<- subset(auc.all, trials %in% responding.trials_short & Condition=="shortstim")

# trace AUC for first 10 s
df1A <- summarySE(auc.all, measurevar="AUC10s", groupvars=c("Condition","ROIType"),na.rm=TRUE)
df1B <- summarySE(auc.all, measurevar="AUC10s", groupvars=c("treatment","Condition","ROIType"),na.rm=TRUE)
df1C <- summarySE(respondingAUC.long, measurevar="AUC10s", groupvars=c("treatment","ROIType"),na.rm=TRUE)
df1D <- summarySE(respondingAUC.short, measurevar="AUC10s", groupvars=c("treatment","ROIType"),na.rm=TRUE)

ggplot(data=df1A, aes(x=ROIType, y=AUC10s, fill=Condition)) +
  geom_errorbar(aes(ymin=AUC10s-se, ymax=AUC10s+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("AUC for 10s of trace") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df1B, aes(x=interaction(ROIType,treatment), y=AUC10s, fill=Condition)) +
  geom_errorbar(aes(ymin=AUC10s-se, ymax=AUC10s+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("AUC for 10s of trace") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df1C, aes(x=ROIType, y=AUC10s, fill=treatment)) +
  geom_errorbar(aes(ymin=AUC10s-se, ymax=AUC10s+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("AUC for 10s of trace") +
  ggtitle("responding long stim")+
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df1D, aes(x=ROIType, y=AUC10s, fill=treatment)) +
  geom_errorbar(aes(ymin=AUC10s-se, ymax=AUC10s+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("AUC for 10s of trace") +
  ggtitle("responding short stim")+
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

# hand circled neurons and FLIK
Cond_ROI_treat=interaction(auc.all$Condition, auc.all$treatment, auc.all$ROIType)
auc.null = lmer(AUC10s ~ (1|Animal) + (1|Spot) , auc.all,REML=FALSE)
auc.model1 = lmer(AUC10s ~ Condition + (1|Animal) + (1|Spot), auc.all,REML=FALSE)
auc.model2A = lmer(AUC10s ~ Condition+ROIType + (1|Animal) + (1|Spot), auc.all,REML=FALSE)
auc.model2B = lmer(AUC10s ~ Condition*ROIType + (1|Animal) + (1|Spot), auc.all,REML=FALSE)
auc.model3A = lmer(AUC10s ~ Condition+ROIType+treatment + (1|Animal) + (1|Spot), auc.all,REML=FALSE)
auc.model3B = lmer(AUC10s ~ Cond_ROI_treat + (1|Animal) + (1|Spot), auc.all,REML=FALSE)
auc.anova <- anova(auc.null, auc.model1,auc.model2A,auc.model2B,auc.model3A,auc.model3B)
print(auc.anova)
# p values
auc.pv.stim2 <- glht(auc.model3B, mcp(Cond_ROI_treat= "Tukey"))
summary(auc.pv.stim2)

##########
# considering peaks in stim window
# when do AC peaks come after neuronal peak?
# bin the number of AC peaks per time point after stim?

df9A <- summarySE(longstim.stimwindow, measurevar="peakTime", groupvars=c("ROIType","treatment"))
df9B <- summarySE(shortstim.stimwindow, measurevar="peakTime", groupvars=c("ROIType","treatment"))

df10A <- summarySE(longstim.stimwindow, measurevar="peakStart", groupvars=c("ROIType","treatment"))
df10B <- summarySE(shortstim.stimwindow, measurevar="peakStart", groupvars=c("ROIType","treatment"))

df11A <- summarySE(longstim.stimwindow, measurevar="peakStartHalf", groupvars=c("ROIType","treatment"))
df11B <- summarySE(shortstim.stimwindow, measurevar="peakStartHalf", groupvars=c("ROIType","treatment"))

df9A$ROIType <- factor(df9A$ROIType , levels = c("Neuron","Neuropil","Endfoot","Soma","Process"))
ggplot(data=df9A, aes(x=ROIType, y=peakTime, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("Mean Peak Time (s)") +
  ggtitle("max peak time for long stim") 

df9B$ROIType <- factor(df9B$ROIType , levels = c("Neuron","Neuropil","Endfoot","Soma","Process"))
ggplot(data=df9B, aes(x=ROIType, y=peakTime, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("Mean Peak Time (s)") +
  ggtitle("max peak time for short stim") 

ggplot(data=df10A, aes(x=ROIType, y=peakStart, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStart-se, ymax=peakStart+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("peakStart Time") +
  ggtitle("peakStart times for long stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df10B, aes(x=ROIType, y=peakStart, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStart-se, ymax=peakStart+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("peakStart Time") +
  ggtitle("peakStart times for short stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme


ggplot(data=df11A, aes(x=ROIType, y=peakStartHalf, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStartHalf-se, ymax=peakStartHalf+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("peakStartHalf Time (t1/2)") +
  ggtitle("peakStartHalf times for long stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df11B, aes(x=ROIType, y=peakStartHalf, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStartHalf-se, ymax=peakStartHalf+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("peakStartHalf Time (t1/2)") +
  ggtitle("peakStartHalf times for short stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

#############
# Likelihood-ratio test 

# long stim, only stim window
ROIType_treatment=interaction(longstim.stimwindow$ROIType,longstim.stimwindow$treatment)
pT.long.null = lmer(peakTime ~ (1|Animal) + (1|Spot), longstim.stimwindow,REML=FALSE)
pT.long.model1 = lmer(peakTime ~ ROIType + (1|Animal) + (1|Spot), longstim.stimwindow,REML=FALSE)
pT.long.model2 = lmer(peakTime ~ treatment + (1|Animal) + (1|Spot), longstim.stimwindow,REML=FALSE)
pT.long.model3 = lmer(peakTime ~ ROIType+treatment + (1|Animal) + (1|Spot), longstim.stimwindow,REML=FALSE)
pT.long.model4 = lmer(peakTime ~ ROIType_treatment + (1|Animal) + (1|Spot), longstim.stimwindow,REML=FALSE)

pT.long.anova <- anova(pT.long.null, pT.long.model1,pT.long.model2,pT.long.model3,pT.long.model4)
print(pT.long.anova)

pT.long.anova <- anova(pT.long.null, pT.long.model1)
print(pT.long.anova)

# p values
pT.pv.longstim <- glht(pT.long.model4, mcp(ROIType_treatment= "Tukey"))
summary(pT.pv.longstim)

pT.pv.longstim <- glht(pT.long.model1, mcp(ROIType= "Tukey"))
summary(pT.pv.longstim)

# short stim, only stim window
ROIType_treatment=interaction(shortstim.stimwindow$ROIType,shortstim.stimwindow$treatment)
pT.short.null = lmer(peakTime ~ (1|Animal) + (1|Spot), shortstim.stimwindow,REML=FALSE)
pT.short.model1 = lmer(peakTime ~ ROIType + (1|Animal) + (1|Spot), shortstim.stimwindow,REML=FALSE)
pT.short.model2 = lmer(peakTime ~ treatment + (1|Animal) + (1|Spot), shortstim.stimwindow,REML=FALSE)
pT.short.model3 = lmer(peakTime ~ ROIType+treatment + (1|Animal) + (1|Spot), shortstim.stimwindow,REML=FALSE)
pT.short.model4 = lmer(peakTime ~ ROIType_treatment + (1|Animal) + (1|Spot), shortstim.stimwindow,REML=FALSE)

pT.short.anova <- anova(pT.short.null, pT.short.model1,pT.short.model2,pT.short.model3,pT.short.model4)
print(pT.short.anova)

pT.short.anova <- anova(pT.short.null, pT.short.model1)
print(pT.short.anova)


# p values
pT.pv.shortstim <- glht(pT.short.model4, mcp(ROIType_treatment= "Tukey"))
summary(pT.pv.shortstim)

pT.pv.shortstim <- glht(pT.short.model1, mcp(ROIType= "Tukey"))
summary(pT.pv.shortstim)

# half maximum time
ROIType_treatment=interaction(longstim.stimwindow$ROIType,longstim.stimwindow$treatment)
HM.null = lmer(peakStartHalf ~ (1|Animal) + (1|Spot), longstim.stimwindow,REML=FALSE)
HM.model1 = lmer(peakStartHalf ~ ROIType+ (1|Animal) + (1|Spot), longstim.stimwindow,REML=FALSE)
HM.anova <- anova(HM.null, HM.model1)
print(HM.anova)
# p values
HM.pv.ROIType <- glht(HM.model1, mcp(ROIType= "Tukey"))
summary(HM.pv.ROIType)





#######
# peak features

# amplitude
df2A <- summarySE(longstim.stimwindow, measurevar="amplitude", groupvars=c("ROIType","treatment"),na.rm=TRUE)
df2B <- summarySE(stim.all3, measurevar="amplitude", groupvars=c("Condition","treatment"),na.rm=TRUE)

df2B <- summarySE(longstim.peaks, measurevar="amplitude", groupvars=c("Condition","ROIType","treatment"))
df2B <- summarySE(stim.all3, measurevar="amplitude", groupvars=c("Condition","ROIType","treatment"))

df2C <- summarySE(longstim.peaks_after, measurevar="amplitude", groupvars=c("Condition","ROIType","treatment"))

ggplot(data=df2A, aes(x=ROIType, y=amplitude, fill=treatment)) +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df2B, aes(x=interaction(Condition,ROIType), y=amplitude, fill=treatment)) +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude") +
  scale_fill_manual(
    values=c("black", "red", "blue","green")) + 
  max.theme

ggplot(data=df2C, aes(x=interaction(Condition,ROIType), y=amplitude, fill=treatment)) +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude") +
  scale_fill_manual(
    values=c("black", "red", "blue","green")) + 
  max.theme


# prominence

df3A <- summarySE(longstim.peaks, measurevar="prominence", groupvars=c("Condition","treatment"))
df3B <- summarySE(longstim.peaks, measurevar="prominence", groupvars=c("Condition","ROIType","treatment"))

ggplot(data=df3A, aes(x=treatment, y=prominence, fill=Condition)) +
  geom_errorbar(aes(ymin=prominence-se, ymax=prominence+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("prominence") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df3B, aes(x=interaction(Condition,ROIType), y=prominence, fill=treatment)) +
  geom_errorbar(aes(ymin=prominence-se, ymax=prominence+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("prominence") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme


# duration
df4A <- summarySE(longstim.peaks, measurevar="Duration", groupvars=c("Condition","treatment"))
df4B <- summarySE(longstim.peaks, measurevar="Duration", groupvars=c("Condition","ROIType","treatment"))

ggplot(data=df4A, aes(x=treatment, y=Duration, fill=Condition)) +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Duration") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df4B, aes(x=interaction(Condition,ROIType), y=Duration, fill=treatment)) +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Duration") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

#frequency (within the same trial)

# nostim frequency
#aggregate data by trial
nostim.trials<- ddply(nostim, c("Animal", "Spot", "trials","Channel","ROIType","treatment","ROIs_trial"), summarise, 
                               PA_mean = mean(peakAUC), nEvents = length(peakAUC),
                               Dur_mean = mean(Duration), Prom_mean = mean(prominence),
                               amp_mean = mean(amplitude), HalfDur = mean(halfWidth),
                              freq = sum(numPeaks), area_mean= mean(area))
nostim.trials$peaks_min=nostim.trials$freq/1.5
nostim.trials$signals_min=nostim.trials$nEvents/1.5

df5A <- summarySE(nostim.trials, measurevar="peaks_min", groupvars=c("Channel"), na.rm = T)
df5B <- summarySE(nostim.trials, measurevar="signals_min", groupvars=c("Channel"),na.rm = T)

ggplot(data=df4A, aes(x=treatment, y=Duration, fill=Condition)) +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Duration") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme
######################

# considering peaks and triggering astrocyte response to neuronal peak onset

# find the mean neuronal peak onset (peak time at half max) for each trial and set this to zero
# subtract this time from all astrocyte peak maximums

# aggregate neurons by trial
longstim.window.trials<- ddply(longstim.stimwindow, c("Animal", "Spot", "trials","ROIType","treatment"), summarise, 
                               PA_mean = mean(peakAUC), nEvents = length(peakAUC), nROIs= length(unique(ROIname)),
                               Dur_mean = mean(Duration), Prom_mean = mean(prominence),
                               amp_mean = mean(amplitude), HalfDur = mean(halfWidth),
                               freq_mean = sum(numPeaks)/length(unique(ROIname)), peakT_mean = mean(peakTime),
                               peakHalf_mean= mean(peakStartHalf), area_mean= mean(area),Active=sum(ActivePeak))

neuronal_meanPeakTimes<-subset(longstim.window.trials, ROIType=="Neuron")
  
neuropil_meanPeakTimes<-subset(longstim.window.trials, ROIType=="Neuropil")


# set neuronal peak max time to zero
longstim.NeuronalPeakTime<-data.frame()

for (iTrials in 1:length(neuronal_meanPeakTimes$trials))
{
  CurrentTrial= neuronal_meanPeakTimes$trials[iTrials]
  CurrentTime= neuronal_meanPeakTimes$peakT_mean[iTrials]
  
  Datasubset = subset(longstim.stimwindow, trials == CurrentTrial)
  Datasubset$peakTime=Datasubset$peakTime-CurrentTime
  Datasubset$peakStart=Datasubset$peakStart-CurrentTime
  Datasubset$peakStartHalf=Datasubset$peakStartHalf-CurrentTime

  longstim.NeuronalPeakTime<- rbind(longstim.NeuronalPeakTime,Datasubset)
  }


# set neuropil peak max time to zero
longstim.NeuropilPeakTime<-data.frame()

for (iTrials in 1:length(neuronal_meanPeakTimes$trials))
{
  CurrentTrial= neuropil_meanPeakTimes$trials[iTrials]
  CurrentTime= neuropil_meanPeakTimes$peakT_mean[iTrials]
  
  Datasubset = subset(longstim.stimwindow, trials == CurrentTrial)
  Datasubset$peakTime=Datasubset$peakTime-CurrentTime
  Datasubset$peakStart=Datasubset$peakStart-CurrentTime
  Datasubset$peakStartHalf=Datasubset$peakStartHalf-CurrentTime
  
  longstim.NeuropilPeakTime<- rbind(longstim.NeuropilPeakTime,Datasubset)
}


# set neuronal half max time to zero
longstim.NeuronalHalfTime<-data.frame()

for (iTrials in 1:length(neuronal_meanPeakTimes$trials))
{
  CurrentTrial= neuronal_meanPeakTimes$trials[iTrials]
  CurrentTime= neuronal_meanPeakTimes$peakHalf_mean[iTrials]
  
  Datasubset = subset(longstim.stimwindow, trials == CurrentTrial)
  Datasubset$peakTime=Datasubset$peakTime-CurrentTime
  Datasubset$peakStart=Datasubset$peakStart-CurrentTime
  Datasubset$peakStartHalf=Datasubset$peakStartHalf-CurrentTime
  
  longstim.NeuronalHalfTime<- rbind(longstim.NeuronalHalfTime,Datasubset)
}


# set neuropil half max time to zero
longstim.NeuropilHalfTime<-data.frame()

for (iTrials in 1:length(neuronal_meanPeakTimes$trials))
{
  CurrentTrial= neuropil_meanPeakTimes$trials[iTrials]
  CurrentTime= neuropil_meanPeakTimes$peakHalf_mean[iTrials]
  
  Datasubset = subset(longstim.stimwindow, trials == CurrentTrial)
  Datasubset$peakTime=Datasubset$peakTime-CurrentTime
  Datasubset$peakStart=Datasubset$peakStart-CurrentTime
  Datasubset$peakStartHalf=Datasubset$peakStartHalf-CurrentTime
  
  longstim.NeuropilHalfTime<- rbind(longstim.NeuropilHalfTime,Datasubset)
}

######
# Plots
df9C <- summarySE(longstim.NeuronalPeakTime, measurevar="peakTime", groupvars=c("ROIType","treatment"))
df9D <- summarySE(longstim.NeuropilPeakTime, measurevar="peakTime", groupvars=c("ROIType","treatment"))
df9E <- summarySE(longstim.NeuronalHalfTime, measurevar="peakTime", groupvars=c("ROIType","treatment"))
df9F <- summarySE(longstim.NeuropilHalfTime, measurevar="peakTime", groupvars=c("ROIType","treatment"))

df9C$ZeroPoint<-"NeuronalPeakTime"
df9D$ZeroPoint<-"NeuropilPeakTime"

df9G<-rbind(df9C,df9D)

df9C$ROIType <- factor(df9C$ROIType , levels = c("Neuron","Neuropil","Endfoot","Soma","Process"))
ggplot(data=df9C, aes(x=ROIType, y=peakTime, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("Mean Peak Time (s)") +
  ggtitle("max peak time set to neuronal peak time for long stim") 

df9D$ROIType <- factor(df9D$ROIType , levels = c("Neuron","Neuropil","Endfoot","Soma","Process"))
ggplot(data=df9D, aes(x=ROIType, y=peakTime, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("Mean Peak Time (s)") +
  ggtitle("max peak time set to neuropil peak time for long stim") 


df9E$ROIType <- factor(df9E$ROIType , levels = c("Neuron","Neuropil","Endfoot","Soma","Process"))
ggplot(data=df9E, aes(x=ROIType, y=peakTime, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("Mean Peak Time (s)") +
  ggtitle("max peak time set to neuronal half onset time for long stim") 

df9F$ROIType <- factor(df9F$ROIType , levels = c("Neuron","Neuropil","Endfoot","Soma","Process"))
ggplot(data=df9F, aes(x=ROIType, y=peakTime, fill=ROIType)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("Mean Peak Time (s)") +
  ggtitle("max peak time set to neuropil half onset time for long stim") 


df9G$ROIType <- factor(df9G$ROIType , levels = c("Neuron","Neuropil","Endfoot","Soma","Process"))
ggplot(data=df9G, aes(x=ROIType, y=peakTime, fill=ZeroPoint)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("Mean Peak Time (s)") +
  ggtitle("max peak time set to neuropil half onset time for long stim") 


longstim.NeuronalPeakTime$ZeroPoint<-"Neuron_PeakTime"
longstim.NeuropilPeakTime$ZeroPoint<-"Neuropil_PeakTime"

Adjusted.peaktimes<-rbind(longstim.NeuronalPeakTime,longstim.NeuropilPeakTime)

Adjusted.peaktimes$ZeroPoint<-as.factor(Adjusted.peaktimes$ZeroPoint)
Adjusted.peaktimes$ZeroPoint_ROIType<-interaction(Adjusted.peaktimes$ROIType,Adjusted.peaktimes$ZeroPoint)

# hand circled neurons and FLIKA
pT.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), Adjusted.peaktimes,REML=FALSE)
pT.model1 = lmer(peakTime ~ ROIType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), Adjusted.peaktimes,REML=FALSE)
pT.model2A = lmer(peakTime ~ ZeroPoint + (1|Animal) + (1|Spot)+ (1|ROIs_trial), Adjusted.peaktimes,REML=FALSE)
pT.model2B = lmer(peakTime ~ ZeroPoint+ROIType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), Adjusted.peaktimes,REML=FALSE)
pT.model3B = lmer(peakTime ~ ZeroPoint_ROIType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), Adjusted.peaktimes,REML=FALSE)
pT.anova <- anova(pT.null, pT.model1,pT.model2A,pT.model2B,pT.model3B)
print(pT.anova)
# p values
pT.pv.longstim <- glht(pT.model1, mcp(ROIType= "Tukey"))
summary(pT.pv.longstim)

pT.pv.longstim2 <- glht(pT.model3B, mcp(ZeroPoint_ROIType= "Tukey"))
summary(pT.pv.longstim2)


########
ggplot(data=df10A, aes(x=ROIType, y=peakStart, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStart-se, ymax=peakStart+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("peakStart Time") +
  ggtitle("peakStart times for long stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df10B, aes(x=ROIType, y=peakStart, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStart-se, ymax=peakStart+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("peakStart Time") +
  ggtitle("peakStart times for short stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme


ggplot(data=df11A, aes(x=ROIType, y=peakStartHalf, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStartHalf-se, ymax=peakStartHalf+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("peakStartHalf Time (t1/2)") +
  ggtitle("peakStartHalf times for long stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df11B, aes(x=ROIType, y=peakStartHalf, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStartHalf-se, ymax=peakStartHalf+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("ROIType") +
  ylab("peakStartHalf Time (t1/2)") +
  ggtitle("peakStartHalf times for short stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme




##########
shortstim.window.trials<- ddply(shortstim.stimwindow, c("Animal", "Spot", "trials","ROIType","treatment"), summarise, 
                                PA_mean = mean(peakAUC), nEvents = length(peakAUC), nROIs= length(unique(ROIname)),
                                Dur_mean = mean(Duration), Prom_mean = mean(prominence),
                                amp_mean = mean(amplitude), HalfDur = mean(halfWidth),
                                freq_mean = sum(numPeaks)/length(unique(ROIname)), peakT_mean = mean(peakTime),
                                peakHalf_mean= mean(peakStartHalf), area_mean= mean(area))



  
  
  
  
df9A <- summarySE(longstim.peaks, measurevar="peakTime", groupvars=c("Condition","ROIType","treatment"))
df10A <- summarySE(longstim.peaks, measurevar="peakStart", groupvars=c("Condition","ROIType","treatment"))
df11A <- summarySE(longstim.peaks, measurevar="peakStartHalf", groupvars=c("Condition","ROIType","treatment"))

ggplot(data=df9A, aes(x=interaction(ROIType,Condition), y=peakTime, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("Channel") +
  ylab("Mean Peak Time (s)") +
  ggtitle("timepeaks for stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df10A, aes(x=interaction(ROIType,Condition), y=peakStart, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStart-se, ymax=peakStart+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("Channel") +
  ylab("peakStart Time") +
  ggtitle("peakStart times for stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df11A, aes(x=interaction(ROIType,Condition), y=peakStartHalf, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStartHalf-se, ymax=peakStartHalf+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("Channel") +
  ylab("peakStartHalf Time (t1/2)") +
  ggtitle("peakStartHalf times for stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

#############
# Likelihood-ratio test 

# hand circled neurons and FLIKA
pT.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), longstim.peaks,REML=FALSE)
pT.model1 = lmer(peakTime ~ Condition + (1|Animal) + (1|Spot)+ (1|ROIs_trial), longstim.peaks,REML=FALSE)
pT.model2A = lmer(peakTime ~ Condition+ROIType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), longstim.peaks,REML=FALSE)
pT.model2B = lmer(peakTime ~ Condition*ROIType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), longstim.peaks,REML=FALSE)
pT.model3A = lmer(peakTime ~ Condition+ROIType+treatment + (1|Animal) + (1|Spot)+ (1|ROIs_trial), longstim.peaks,REML=FALSE)
pT.model3B = lmer(peakTime ~ Condition*ROIType*treatment + (1|Animal) + (1|Spot)+ (1|ROIs_trial), longstim.peaks,REML=FALSE)
pT.anova <- anova(pT.null, pT.model1,pT.model2A,pT.model2B,pT.model3A,pT.model3B)
print(pT.anova)
# p values
pT.pv.stim2 <- lsmeans(pT.model3B, pairwise ~ Condition*ROIType*treatment, glhargs=list())
summary(pT.pv.stim2)

# half maximum time
HM.null = lmer(onset ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), longstim.peaks,REML=FALSE)
HM.model1 = lmer(onset ~ Condition + (1|Animal) + (1|Spot)+ (1|ROIs_trial), longstim.peaks,REML=FALSE)
HM.model2A = lmer(onset ~ Condition+ROIType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), longstim.peaks,REML=FALSE)
HM.model2B = lmer(onset ~ Condition*ROIType + (1|Animal) + (1|Spot)+ (1|ROIs_trial), longstim.peaks,REML=FALSE)
HM.model3A = lmer(onset ~ Condition+ROIType+treatment + (1|Animal) + (1|Spot)+ (1|ROIs_trial), longstim.peaks,REML=FALSE)
HM.model3B = lmer(onset ~ Condition*ROIType*treatment + (1|Animal) + (1|Spot)+ (1|ROIs_trial), longstim.peaks,REML=FALSE)
HM.anova <- anova(HM.null, HM.model1,HM.model2A,HM.model2B,HM.model3A,HM.model3B)
print(HM.anova)
# p values
HM.pv.stim2 <- lsmeans(HM.model3B, pairwise ~ Condition*ROIType*treatment, glhargs=list())
summary(HM.pv.stim2)


#amplitude
amp.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), longstim.peaks,REML=FALSE)
amp.model1 = lmer(amplitude ~ Condition + (1|Animal) + (1|Spot) + (1|ROIs_trial), longstim.peaks,REML=FALSE)
amp.model2A = lmer(amplitude~ Condition+ROIType + (1|Animal) + (1|Spot) + (1|ROIs_trial), longstim.peaks,REML=FALSE)
amp.model2B = lmer(amplitude~ Condition*ROIType + (1|Animal) + (1|Spot) + (1|ROIs_trial), longstim.peaks,REML=FALSE)
amp.model3A = lmer(amplitude~ Condition+ROIType+treatment + (1|Animal) + (1|Spot) + (1|ROIs_trial), longstim.peaks,REML=FALSE)
amp.model3B = lmer(amplitude~ Condition*ROIType*treatment + (1|Animal) + (1|Spot) + (1|ROIs_trial), longstim.peaks,REML=FALSE)
amp.anova1 <- anova(amp.null, amp.model1,amp.model2A,amp.model2B,amp.model3A,amp.model3B)
print(amp.anova1)
# p values
amp.pv.stim_type <- lsmeans(amp.model3B, pairwise ~ Condition*ROIType*treatment, glhargs=list())
summary(amp.pv.stim_type)

##########
# SHORT STIM (90Hz, 1 sec)

shortstim.all<-subset(stim.all3, Condition!="Stim")
shortstim.peaks<- subset(shortstim.all, peakTime>=0 & peakTime<=5)
#shortstim.peaks<- subset(stim.all3, peakTime>5 & peakTime<80)


# amplitude
df2A <- summarySE(shortstim.peaks, measurevar="amplitude", groupvars=c("Condition","treatment"))
df2B <- summarySE(shortstim.peaks, measurevar="amplitude", groupvars=c("Condition","ROIType","treatment"))

ggplot(data=df2A, aes(x=treatment, y=amplitude, fill=Condition)) +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df2B, aes(x=interaction(Condition,ROIType), y=amplitude, fill=treatment)) +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("amplitude") +
  scale_fill_manual(
    values=c("black", "red", "blue","green")) + 
  max.theme


# prominence

df3A <- summarySE(shortstim.peaks, measurevar="prominence", groupvars=c("Condition","treatment"))
df3B <- summarySE(shortstim.peaks, measurevar="prominence", groupvars=c("Condition","ROIType","treatment"))

ggplot(data=df3A, aes(x=treatment, y=prominence, fill=Condition)) +
  geom_errorbar(aes(ymin=prominence-se, ymax=prominence+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("prominence") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df3B, aes(x=interaction(Condition,ROIType), y=prominence, fill=treatment)) +
  geom_errorbar(aes(ymin=prominence-se, ymax=prominence+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("prominence") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme


# duration
df4A <- summarySE(shortstim.peaks, measurevar="Duration", groupvars=c("Condition","treatment"))
df4B <- summarySE(shortstim.peaks, measurevar="Duration", groupvars=c("Condition","ROIType","treatment"))

ggplot(data=df4A, aes(x=treatment, y=Duration, fill=Condition)) +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Duration") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df4B, aes(x=interaction(Condition,ROIType), y=Duration, fill=treatment)) +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("Duration") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

#frequency (within the same trial)


######################

# considering peaks
# when do AC peaks come after neuronal peak?
# bin the number of AC peaks per time point after stim?

df9A <- summarySE(shortstim.peaks, measurevar="peakTime", groupvars=c("Condition","ROIType","treatment"))
df10A <- summarySE(shortstim.peaks, measurevar="peakStart", groupvars=c("Condition","ROIType","treatment"))
df11A <- summarySE(shortstim.peaks, measurevar="peakStartHalf", groupvars=c("Condition","ROIType","treatment"))

ggplot(data=df9A, aes(x=interaction(ROIType,Condition), y=peakTime, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("Channel") +
  ylab("Mean Peak Time (s)") +
  ggtitle("timepeaks for stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df10A, aes(x=interaction(ROIType,Condition), y=peakStart, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStart-se, ymax=peakStart+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("Channel") +
  ylab("onset Time start") +
  ggtitle("onset times for stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme

ggplot(data=df11A, aes(x=interaction(ROIType,Condition), y=peakStartHalf, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakStartHalf-se, ymax=peakStartHalf+se), colour="black", width=.1,  position=position_dodge(.9)) +
  coord_flip() +
  xlab("Channel") +
  ylab("peakStartHalfTime start") +
  ggtitle("peakStartHalf times for short stim") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme


####################
# aggregate peaks from soma, endfoot, neuron ROIs
struct.peaks<- subset(stim.all3, ROIType!="Process")

# fraction response (active trials)
struct.peaks$ActivePeak <- 0
farpeaks1 <- struct.peaks$peakTime>=0 & struct.peaks$peakTime<=20
struct.peaks$ActivePeak[farpeaks1] <- 1 

struct.trials<- ddply(struct.peaks, c("Animal", "Spot", "Condition","Channel","ROIs","Trial","ROIType","treatment"), summarise, 
                   PA_mean = mean(peakAUC,na.rm=TRUE), nEvents = sum(numPeaks,na.rm=TRUE),
                   Dur_mean = mean(fullWidth,na.rm=TRUE), Prom_mean = mean(prominence,na.rm=TRUE),
                   amp_mean = mean(amplitude,na.rm=TRUE), HalfDur = mean(halfWidth,na.rm=TRUE),
                   numActive = sum(ActivePeak))

struct.trials$ActiveTrial<-0
peaks = struct.trials$numActive>0
struct.trials$ActiveTrial[peaks] <- 1
struct.trials$freq<- (struct.trials$nEvents/85)*60

# aggregate ROI (no timepoints)
struct.ROI<- ddply(struct.trials, c("Animal", "Spot", "Condition","Channel","ROIs","ROIType","treatment"), summarise, 
                 PA_mean2 = mean(PA_mean,na.rm=TRUE), PA_SD = sd(PA_mean,na.rm=TRUE),
                 freq_mean2 = mean(freq,na.rm=TRUE), freq_SD = sd(freq,na.rm=TRUE),
                 Dur_mean2 = mean(Dur_mean,na.rm=TRUE), Dur_SD = sd(Dur_mean,na.rm=TRUE),
                 Prom_mean2 = mean(Prom_mean,na.rm=TRUE), Prom_SD = sd(Prom_mean,na.rm=TRUE),
                 amp_mean2 = mean(amp_mean,na.rm=TRUE), amp_SD = sd(amp_mean,na.rm=TRUE),
                 HalfDur_mean2 = mean(HalfDur,na.rm=TRUE), HalfDur_SD = sd(HalfDur,na.rm=TRUE),
                 TotEvents = sum(nEvents),TotActive = sum(ActiveTrial), Ntrials= length(amp_mean))

struct.ROI$frac.resp <- struct.ROI$TotActive/struct.ROI$Ntrials

#######
#all.ROI2 <- rbind(wholeframe.ROI, autoRCaMP.ROI, HCRCaMP.ROI2)

# mean percent response- i.e. the fraction of trials where a peak is detected in the first 20 sec
df5A <- summarySE(struct.ROI, measurevar="frac.resp", groupvars=c("Condition","ROIType","treatment"))


ggplot(data=df5A, aes(x=interaction(ROIType,Condition), y=frac.resp, fill=treatment)) +
  geom_errorbar(aes(ymin=frac.resp-se, ymax=frac.resp+se), colour="black", width=.5, size= 1, position=position_dodge(1.0)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("frac.resp") +
  scale_fill_manual(
    values=c("black", "red", "blue")) + 
  max.theme


#fraction of response
# Likelihood-ratio test 

frac.resp.null = lmer(frac.resp ~ (1|Animal) + (1|Spot) + (1|ROI), struct.ROI,REML=FALSE)
frac.resp.model1 = lmer(frac.resp ~ Condition + (1|Animal) + (1|Spot)+ (1|ROI), struct.ROI,REML=FALSE)
frac.resp.model3A = lmer(frac.resp~ Condition+ROIType + (1|Animal) + (1|Spot)+ (1|ROI), struct.ROI,REML=FALSE)
frac.resp.model3B = lmer(frac.resp~ Condition*ROIType + (1|Animal) + (1|Spot)+ (1|ROI), struct.ROI,REML=FALSE)
frac.resp.anova <- anova(frac.resp.null, frac.resp.model1,frac.resp.model3A,frac.resp.model3B)
print(frac.resp.anova)
# p values
frac.resp.pv.stim_type <- lsmeans(frac.resp.model3B, pairwise ~ Condition*ROIType, glhargs=list())
summary(frac.resp.pv.stim_type)



######
df6A <- summarySE(struct.ROI, measurevar="freq_mean2", groupvars=c("Condition","ROIType","treatment"))


ggplot(data=df6A, aes(x=interaction(ROIType,Condition), y=freq_mean2, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=freq_mean2-se, ymax=freq_mean2+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Condition") +
  ylab("Events per min") +
  ggtitle("Frequency for struc ROIs") +
  scale_fill_manual(
    values=c("black", "red", "blue")) +
  max.theme


#frequency
freq.null = lmer(freq_mean2 ~ (1|Animal) + (1|Spot) + (1|ROI), struct.ROI,REML=FALSE)
freq.model1 = lmer(freq_mean2 ~ Condition + (1|Animal) + (1|Spot)+ (1|ROI), struct.ROI,REML=FALSE)
freq.model3A = lmer(freq_mean2 ~ Condition+ROIType + (1|Animal) + (1|Spot)+ (1|ROI), struct.ROI,REML=FALSE)
freq.model3B = lmer(freq_mean2 ~ Condition*ROIType + (1|Animal) + (1|Spot)+ (1|ROI), struct.ROI,REML=FALSE)
freq.anova <- anova(freq.null, freq.model1,freq.model3A,freq.model3B)
print(freq.anova)
# p values
freq.pv.stim_type <- lsmeans(freq.model3B, pairwise ~ Condition*ROIType, glhargs=list())
summary(freq.pv.stim_type)


#############
df7A<- summarySE(struct.ROI, measurevar="amp_mean2", groupvars=c("Condition","ROIType","treatment"))


ggplot(data=df7A, aes(x=interaction(ROIType,Condition), y=amp_mean2, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amp_mean2-se, ymax=amp_mean2+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Channel") +
  ylab("amp") +
  ggtitle("amp for WF ROIs") +
  scale_fill_manual(
    values=c("black", "red", "blue")) +
  max.theme


#amplitude
amp.null = lmer(amp_mean2 ~ (1|Animal) + (1|Spot) + (1|ROI), struct.ROI,REML=FALSE)
amp.model1 = lmer(amp_mean2 ~ Condition + (1|Animal) + (1|Spot)+ (1|ROI), struct.ROI,REML=FALSE)
amp.model3A = lmer(amp_mean2 ~ Condition+ROIType + (1|Animal) + (1|Spot)+ (1|ROI), struct.ROI,REML=FALSE)
amp.model3B = lmer(amp_mean2 ~ Condition*ROIType + (1|Animal) + (1|Spot)+ (1|ROI), struct.ROI,REML=FALSE)
amp.anova <- anova(amp.null, amp.model1,amp.model3A,amp.model3B)
print(amp.anova)
# p values
amp.pv.stim_type <- lsmeans(amp.model3B, pairwise ~ Condition*ROIType, glhargs=list())
summary(amp.pv.stim_type)

############
df8A <- summarySE(struct.ROI, measurevar="Dur_mean2", groupvars=c("Condition","ROIType","treatment"))


ggplot(data=df8A, aes(x=interaction(ROIType,Condition), y=Dur_mean2, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Dur_mean2-se, ymax=Dur_mean2+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Channel") +
  ylab("Duration") +
  ggtitle("Duration for all ROIs") +
  scale_fill_manual(
    values=c("black", "red", "blue")) +
  max.theme

#duration
dur.null = lmer(Dur_mean2 ~ (1|Animal) + (1|Spot) + (1|ROI), struct.ROI,REML=FALSE)
dur.model1 = lmer(Dur_mean2 ~ Condition + (1|Animal) + (1|Spot)+ (1|ROI), struct.ROI,REML=FALSE)
dur.model3A = lmer(Dur_mean2 ~ Condition+ROIType + (1|Animal) + (1|Spot)+ (1|ROI), struct.ROI,REML=FALSE)
dur.model3B = lmer(Dur_mean2 ~ Condition*ROIType + (1|Animal) + (1|Spot)+ (1|ROI), struct.ROI,REML=FALSE)
dur.anova <- anova(dur.null, dur.model1,dur.model3A,dur.model3B)
print(dur.anova)
# p values
dur.pv.stim_type <- lsmeans(dur.model3A, pairwise ~ Condition*ROIType, glhargs=list())
summary(dur.pv.stim_type)




###########
# BUT! We should correlate specific AC peaks with specific neuron peaks from the SAME TRIAL

# count the number of AC peaks that occur in the same trials and neuronal peaks
# bin by timepoint

NeuronStim= subset(longstim.peaks, Condition=="Stim" & Channel=="RCaMP")
ACStim= subset(longstim.peaks, Condition=="Stim" & Channel =="GCaMP")

ActiveNeuronTrials= unique(NeuronStim$trials)

NeuronTrials = data.frame()
Stim_pooled = data.frame()
Dis= data.frame()
for (ipeaks in 1:length(ActiveNeuronTrials))
{
  CurrentTrial= ActiveNeuronTrials[ipeaks]
  ACsubset = subset(ACStim, trials == CurrentTrial)
  #ACsubset2 = subset(AC_stim_peaks, trial == CurrentTrial)
  Neuronsub = subset(NeuronStim, trials==CurrentTrial)
  # subset AC data based on peak time
  AC0 = subset(ACsubset, peakTime<0)
  AC5 = subset(ACsubset, peakTime<5 & peakTime>=0)
  AC10 = subset(ACsubset, peakTime<10 & peakTime>=5)
  AC15 = subset(ACsubset, peakTime<15 & peakTime>=10)
  AC20 = subset(ACsubset, peakTime<20 & peakTime>=15)
  
  # count AC peaks in each bin
  Data= Neuronsub[1,]
  Data$AC0 = sum(length(AC0$peakTime))
  Data$AC5 = sum(length(AC5$peakTime))
  Data$AC10 = sum(length(AC10$peakTime))
  Data$AC15 = sum(length(AC15$peakTime))
  Data$AC20 = sum(length(AC20$peakTime))
  NeuronTrials = rbind(NeuronTrials, Data)
  
  # calculate mean peakTime and onset for neuronal signals and astrocyte signals in each trial
  if (length(ACsubset[,1])>0)
  {
    Data_pooled1<- ddply(Neuronsub, c("Animal", "Spot", "Condition","trials","ROIType","treatment"), summarise, 
                         Channel="Neuron", PT_mean = mean(peakTime), PT_SD = sd(peakTime),
                         onset_mean = mean(onset), onset_SD = sd(onset),
                         amp_mean = mean(amplitude), amp_SD = sd(amplitude),
                         Peaknum = sum(length(amplitude)))
    Data_pooled2<- ddply(ACsubset, c("Animal", "Spot", "Condition","trials","ROIType","treatment"), summarise, 
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

df11A <- summarySE(Stim_pooled, measurevar="PT_mean", groupvars=c("Condition","ROIType","treatment"))
df12A <- summarySE(Stim_pooled, measurevar="onset", groupvars=c("Condition","ROIType","treatment"))

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
