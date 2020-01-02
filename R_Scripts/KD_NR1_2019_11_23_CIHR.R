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

###########
# NOTES
# COLOUR BLIND FRIENDLY PALETTE FOR PLOTS
# The palette with black:
#cbbPalette <- c("#D55E00","#009E73","#E69F00","#56B4E9","#CC79A7","#F0E442")
cbbPalette <- c("#D55E00", "#56B4E9", "#009E73","#F0E442", "#0072B2", "#E69F00", "#CC79A7")

########################
# load data

#Jill's home files
# peak data
control.peaks.12<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/12_control_peaks_longtrials.csv",  header=TRUE, sep = ",")
control.peaks.14<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/14_control_peaks_longtrials.csv",  header=TRUE, sep = ",")
control.peaks.15<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/15_control_peaks_longtrials.csv",  header=TRUE, sep = ",")
evening.peaks.12<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/12_control_evening_peaks_longtrials.csv",  header=TRUE, sep = ",")
evening.peaks.14<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/14_control_evening_peaks_longtrials.csv",  header=TRUE, sep = ",")
evening.peaks.15<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/15_control_evening_peaks_longtrials.csv",  header=TRUE, sep = ",")
control.peaks.crazy8<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/Crazy8_control_peaks_longtrials.csv",  header=TRUE, sep = ",")
evening.peaks.crazy8<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/Crazy8_control_evening_peaks_longtrials.csv",  header=TRUE, sep = ",")
KD.peaks.12<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/12_KD_peaks_longtrials.csv",  header=TRUE, sep = ",")
KD.peaks.14<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/14_KD_peaks_longtrials.csv",  header=TRUE, sep = ",")
KD.peaks.15<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/15_KD_peaks_longtrials.csv",  header=TRUE, sep = ",")

# OT
control.OT.12<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/12_control_onset_time_longtrials.csv",  header=TRUE, sep = ",")
control.OT.14<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/14_control_onset_time_longtrials.csv",  header=TRUE, sep = ",")
control.OT.15<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/15_control_onset_time_longtrials.csv",  header=TRUE, sep = ",")
evening.OT.12<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/12_control_evening_onset_time_longtrials.csv",  header=TRUE, sep = ",")
evening.OT.14<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/14_control_evening_onset_time_longtrials.csv",  header=TRUE, sep = ",")
evening.OT.15<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/15_control_evening_onset_time_longtrials.csv",  header=TRUE, sep = ",")
control.OT.crazy8<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/Crazy8_control_onset_time_longtrials.csv",  header=TRUE, sep = ",")
evening.OT.crazy8<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/Crazy8_control_evening_onset_time_longtrials.csv",  header=TRUE, sep = ",")
KD.OT.12<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/12_KD_onset_time_longtrials.csv",  header=TRUE, sep = ",")
KD.OT.14<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/14_KD_onset_time_longtrials.csv",  header=TRUE, sep = ",")
KD.OT.15<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/15_KD_onset_time_longtrials.csv",  header=TRUE, sep = ",")

#field data
control.field.12<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/12_control_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
control.field.14<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/14_control_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
control.field.15<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/15_control_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
evening.field.12<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/12_control_evening_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
evening.field.14<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/14_control_evening_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
evening.field.15<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/15_control_evening_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
control.field.crazy8<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/Crazy8_control_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
evening.field.crazy8<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/Crazy8_control_evening_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
KD.field.12<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/12_KD_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
KD.field.14<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/14_KD_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
KD.field.15<-read.table("J:/Noushin Ahmadpour/2P Data/Results/FilesforR/15_KD_Lck_field_longtrials.csv",  header=TRUE, sep = ",")


##########

#lsm.options(pbkrtest.limit = 100000)

#change animal name
control.OT.14$Animal<-"14"
control.OT.15$Animal<-"15"
control.OT.crazy8$Animal<-"Crazy8"
evening.OT.14$Animal<-"14"
evening.OT.15$Animal<-"15"
evening.OT.crazy8$Animal<-"Crazy8"
KD.OT.14$Animal<-"14"
KD.OT.15$Animal<-"15"

control.peaks.14$Animal<-"14"
control.peaks.15$Animal<-"15"
control.peaks.crazy8$Animal<-"Crazy8"
evening.peaks.14$Animal<-"14"
evening.peaks.15$Animal<-"15"
evening.peaks.crazy8$Animal<-"Crazy8"
KD.peaks.14$Animal<-"14"
KD.peaks.15$Animal<-"15"

control.field.14$animalname<-"14"
control.field.15$animalname<-"15"
control.field.crazy8$animalname<-"Crazy8"
evening.field.14$animalname<-"14"
evening.field.15$animalname<-"15"
evening.field.crazy8$animalname<-"Crazy8"
KD.field.14$animalname<-"14"
KD.field.15$animalname<-"15"

#  merge data
control.OT<-rbind(control.OT.12, control.OT.14, control.OT.15, control.OT.crazy8)
evening.OT<-rbind(evening.OT.12, evening.OT.14, evening.OT.15, evening.OT.crazy8)
KD.OT<-rbind(KD.OT.12, KD.OT.14, KD.OT.15)

control.OT$Type<-"Control"
control.OT$Spot<-paste(substr(control.OT$Spot,2, 6),"control",sep= "_")
evening.OT$Type<-"Evening"
evening.OT$Spot<-paste(substr(evening.OT$Spot,2, 6),"evening",sep= "_")
KD.OT$Type<-"KD"
KD.OT$Spot<-paste(substr(KD.OT$Spot,2, 6),"KD",sep= "_")


all.OT<-rbind(control.OT, evening.OT, KD.OT)

# peaks
control.peaks<-rbind(control.peaks.12, control.peaks.14, control.peaks.15, control.peaks.crazy8)
evening.peaks<-rbind(evening.peaks.12, evening.peaks.14, evening.peaks.15, evening.peaks.crazy8)
KD.peaks<-rbind(KD.peaks.12, KD.peaks.14, KD.peaks.15)

control.peaks$Type<-"Control"
control.peaks$Spot<-paste(substr(control.peaks$Spot,2, 6),"control",sep= "_")
evening.peaks$Type<-"Evening"
evening.peaks$Spot<-paste(substr(evening.peaks$Spot,2, 6),"evening",sep= "_")
KD.peaks$Type<-"KD"
KD.peaks$Spot<-paste(substr(KD.peaks$Spot,2, 6),"KD",sep= "_")

all.peaks<-rbind(control.peaks, evening.peaks, KD.peaks)


# field
control.field<-rbind(control.field.12, control.field.14, control.field.15, control.field.crazy8)
evening.field<-rbind(evening.field.12, evening.field.14, evening.field.15, evening.field.crazy8)
KD.field<-rbind(KD.field.12, KD.field.14, KD.field.15)


control.field$Type<-"Control"
control.field$Spot<-paste(substr(control.field$Spot,2, 6),"control",sep= "_")
evening.field$Type<-"Evening"
evening.field$Spot<-paste(substr(evening.field$Spot,2, 6),"evening",sep= "_")
KD.field$Type<-"KD"
KD.field$Spot<-paste(substr(KD.field$Spot,2, 6),"KD",sep= "_")

all.field<-rbind(control.field, evening.field, KD.field)



all.peaks$Duration<- all.peaks$halfWidth*2

#add baseline time to peaks table
all.peaks$BL_time<-5.032

# adjust peak time and duration
all.peaks$peakTime<- all.peaks$peakTime-all.peaks$BL_time  # now time= 0s is the start of stimulation
all.peaks$peakStart<- all.peaks$peakStart-all.peaks$BL_time
all.peaks$peakStartHalf<- all.peaks$peakStartHalf-all.peaks$BL_time



# set these new variables as factors so we can do stats on them
all.peaks$Type<-as.factor(all.peaks$Type)
all.OT$Type<-as.factor(all.OT$Type)
all.field$Type<-as.factor(all.field$Type)

######
# days post injection
#all.OT$DaysPostInjection="47"

#all.OT$DaysPostInjection[grepl("06_15",all.OT$Spot)]="28"
#all.OT$DaysPostInjection[grepl("06_21",all.OT$Spot)]="35"
#all.OT$DaysPostInjection[grepl("06_22",all.OT$Spot)]="35"
#all.OT$DaysPostInjection[grepl("06_28",all.OT$Spot)]="42"
#all.OT$DaysPostInjection[grepl("07_02",all.OT$Spot)]="46"

#all.peaks$Timepoint[grepl("06_21",all.peaks$Spot)]="early"
#all.peaks$Timepoint[grepl("06_22",all.peaks$Spot)]="early"
#all.peaks$Timepoint[grepl("06_28",all.peaks$Spot)]="mid"


# set these new variables as factors so we can do stats on them
#all.peaks$DaysPostInjection<-as.factor(all.peaks$DaysPostInjection)
#all.OT$DaysPostInjection<-as.factor(all.OT$DaysPostInjection)
#all.field$DaysPostInjection<-as.factor(all.field$DaysPostInjection)

#all.peaks$Timepoint<-as.factor(all.peaks$Timepoint)
#all.OT$Timepoint<-as.factor(all.OT$Timepoint)
#all.field$Timepoint<-as.factor(all.field$Timepoint)

# set the order of groups for plots
all.OT$Type<-factor(all.OT$Type,levels=c("Control","KD", "Evening"))
all.peaks$Type<-factor(all.peaks$Type,levels=c("Control","KD", "Evening"))
all.field$Type<-factor(all.field$Type,levels=c("Control","KD", "Evening"))

####
#unique ROI names
all.OT$ROIs_trial<-paste(all.OT$Animal, all.OT$Spot, all.OT$Trial,all.OT$ROI, sep= "_")
all.OT$ROIs_trial_Cond<-paste(all.OT$ROIs_trial, all.OT$Condition, sep= "_")
all.OT$Spot_type<-paste(all.OT$Spot, all.OT$Type, sep= "_")

all.peaks$ROIs_trial<-paste(all.peaks$Animal, all.peaks$Spot, all.peaks$Trial,all.peaks$roiName, sep= "_")
all.peaks$ROIs_trial_Cond<-paste(all.peaks$ROIs_trial, all.peaks$Condition, sep= "_")
all.peaks$Spot_type<-paste(all.peaks$Spot, all.peaks$Type, sep= "_")

all.field$Spot_trial<-paste(all.field$animalname, all.field$Spot, all.field$trialname, sep= "_")
all.field$Spot_trial_Cond<-paste(all.field$Spot_trial, all.field$Cond, sep= "_")
all.field$Spot_type<-paste(all.field$Spot, all.field$Type, sep= "_")


# count number of trials per spot
Spot.ntrials<-ddply(all.peaks, c("Animal","Type","Spot_type","Condition"), summarise, nTrials=length(unique(Trial)))
Spot.ntrials$Ani_Spot_Cond<-paste(Spot.ntrials$Animal, Spot.ntrials$Spot, Spot.ntrials$Condition, sep="_")


# remove matching dendrite and neuronal soma ROIs
Overlap= all.peaks$overlap!=0
all.peaks<-all.peaks[!Overlap,]



##########################
# FIELD DATA

# remove random trials of NaNs
all.field<-all.field[!is.na(all.field$pixelsize),]
all.field$Spot_type<-paste(all.field$Type, all.field$Spot)

######
# fraction of active pixels from all the astrocyte pixels
# with peaks near stimulus
#all.field$nFluoPix[all.field$nFluoPix==0]<-128*128  #adjust spots where no pixels were above the threshold

all.field$FracActive=all.field$nActivePix/all.field$nFluoPix
#all.field$FracFluo=all.field$nFluoPix/all.field$nTotalPix

all.field.spot<-ddply(all.field, c("Spot","animalname", "Cond", "Type", "Spot_type","Response_Score"), summarise, 
                      meanFracActive=mean(FracActive), meanFluoPix = mean(nFluoPix), nTrial=length(FracActive))
nostim.field.spot<-subset(all.field.spot, Cond=="Nostim")
stim.field.spot<-subset(all.field.spot, Cond=="Stim")

df.FracActive<- summarySE(all.field.spot, measurevar = "meanFracActive", groupvars = c("Cond","Type"))
df.FracActive.nostim1<-summarySE(nostim.field.spot, measurevar = "meanFracActive", groupvars = c("Type"))
df.FracActive.stim1<-summarySE(stim.field.spot, measurevar = "meanFracActive", groupvars = c("Type"))


df.ResponseScore<- summarySE(all.field.spot, measurevar = "Response_Score", groupvars = c("Cond","Type"))

df.nFluoPix.nostim1<-summarySE(nostim.field.spot, measurevar = "meanFluoPix", groupvars = c("Type"))
df.nFluoPix.stim1<-summarySE(stim.field.spot, measurevar = "meanFluoPix", groupvars = c("Type"))


# graphs fraction of active pixels
ggplot(data=df.FracActive, aes(x=Type, y= meanFracActive, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanFracActive-se, ymax=meanFracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


# graphs response probabliity (potential for same region to appear in mutliple trials)
ggplot(data=df.ResponseScore, aes(x=Type, y= Response_Score, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Response_Score-se, ymax=Response_Score+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ResponseScore") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# average fluorescent pixels
ggplot(data=df.nFluoPix.nostim1, aes(x=Type, y= meanFluoPix, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanFluoPix-se, ymax=meanFluoPix+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("# of fluorescent pixels") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

## STATS

cont_vs_kd_fracActive<-subset(all.field.spot, Type!="Evening")
Cond_Type=interaction(cont_vs_kd_fracActive$Cond, cont_vs_kd_fracActive$Type)

#FracActive
fracAct.null = lmer(meanFracActive ~ (1|animalname) + (1|Spot), cont_vs_kd_fracActive,REML=FALSE)
fracAct.model1 = lmer(meanFracActive ~ Cond + (1|animalname) + (1|Spot), cont_vs_kd_fracActive,REML=FALSE)
fracAct.model2 = lmer(meanFracActive ~ Type + (1|animalname) + (1|Spot), cont_vs_kd_fracActive,REML=FALSE)
fracAct.model3 = lmer(meanFracActive ~ Cond_Type + (1|animalname) + (1|Spot), cont_vs_kd_fracActive,REML=FALSE)

fracAct.anova <- anova(fracAct.null, fracAct.model1, fracAct.model2, fracAct.model3)
print(fracAct.anova)

fracAct.Type<- glht(fracAct.model3, mcp(Cond_Type= "Tukey"))
summary(fracAct.Type)

#################################
# PEAK DATA


#consider onsly positive amplitude and accurate time scale
all.peaks.window<-subset(all.peaks, peakTime>0 & peakTime<12 & amplitude>0 & Duration<15)

# RCaMP peaks (nostim and stim) near stimulation time
all.peaks.RC<- subset(all.peaks.window, Channel=="RCaMP")

#spontaneous RCaMP peaks from whole trial
nostim.peaks.RC<-subset(all.peaks.window, Condition=="Nostim" & Channel=="RCaMP")
#nostim.peaks.RC<-subset(nostim.peaks.RC, amplitude>0 & Duration<15)

# stim RCaMP peaks only from time near stimulation
stim.peaks.RC<-subset(all.peaks.RC, Condition=="Stim")

# GCaMP peaks (nostim and stim) near stimulation time
all.peaks.GC<- subset(all.peaks.window, Channel=="GCaMP")

#spontaneous RCaMP peaks from whole trial
nostim.peaks.GC<-subset(all.peaks.window, Condition=="Nostim" & Channel=="GCaMP")
#nostim.peaks.GC<-subset(nostim.peaks.GC, amplitude>0 & Duration<15)

# stim RCaMP peaks only from time near stimulation
stim.peaks.GC<-subset(all.peaks.GC, Condition=="Stim")


######

# peakTime histograms
# all the peaks in the trial
ggplot(all.peaks[all.peaks$Condition=="Stim" & all.peaks$Channel=="RCaMP",],aes(x=peakTime,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("all RCaMP stim peak times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(all.peaks[all.peaks$Condition=="Nostim" & all.peaks$Channel=="RCaMP",],aes(x=peakTime,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("all RCaMP no stim peak times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(all.peaks[all.peaks$Condition=="Stim" & all.peaks$Channel=="GCaMP",],aes(x=peakTime,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("all GCaMP stim peak times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(all.peaks[all.peaks$Condition=="Nostim" & all.peaks$Channel=="GCaMP",],aes(x=peakTime,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("all GCaMP no stim peak times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.RC,aes(x=peakTime,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("RCaMP stim peak times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.GC,aes(x=peakTime,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.2,position="dodge") +
  ggtitle("GCaMP signal peak times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme


# amplitude histograms
ggplot(stim.peaks.RC,aes(x=amplitude,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("RCaMP stim amplitudes")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.RC[stim.peaks.RC$Type!="Evening",],aes(x=amplitude,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("RCaMP stim amplitudes")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.RC[(stim.peaks.RC$Type!="Evening" & stim.peaks.RC$amplitude>2.5),],aes(x=amplitude,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("RCaMP stim amplitudes")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(nostim.peaks.RC,aes(x=amplitude,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("RCaMP nostim amplitudes")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(nostim.peaks.RC[nostim.peaks.RC$Type!="Evening",],aes(x=amplitude,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("RCaMP nostim amplitudes")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(all.peaks.GC[(all.peaks.GC$Condition=="Stim"),],aes(x=amplitude,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.2,position="dodge") +
  ggtitle("GCaMP signal amplitude")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme



#duration histograms
ggplot(all.peaks.RC[(all.peaks.RC$Condition=="Stim"),],aes(x=Duration,y=..density..,fill=Type)) +
  geom_histogram(binwidth=1, position="dodge") +
  ggtitle("RCaMP signal duration")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(all.peaks.GC[(all.peaks.GC$Condition=="Stim"),],aes(x=Duration,y=..density..,fill=Type)) +
  geom_histogram(binwidth=1,position="dodge") +
  ggtitle("GCaMP signal duration")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

#############

#mean amplitude of signals for each ROI for all peaks

#RCaMP
df.amp1.RC<- summarySE(all.peaks.RC, measurevar = "amplitude", groupvars = c("Condition","Type"))
df.amp3.RC<-summarySE(all.peaks.RC, measurevar = "amplitude", groupvars = c("Type"))

# RCaMP spontaneous only
df.amp1.RC.nostim<- summarySE(nostim.peaks.RC, measurevar = "amplitude", groupvars = c("Type"))

# RCaMP stim only
df.amp1.RC.stim<- summarySE(stim.peaks.RC, measurevar = "amplitude", groupvars = c("Type"))

# RCaMP graphs
ggplot(data=df.amp1.RC, aes(x=Type, y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp3.RC, aes(x=Type, y= amplitude, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("RCaMP- all peaks- nostim and stim together") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


# no stim
ggplot(data=df.amp1.RC.nostim, aes(x=Type, y= amplitude, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("no stim amplitude") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(nostim.peaks.RC, aes(x=Type,y=amplitude, fill= Type)) +
  geom_boxplot()+
  ylab("no stim RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


# stim
ggplot(data=df.amp1.RC.stim, aes(x=Type, y= amplitude, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("stim amplitude") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.RC, aes(x=Type,y=amplitude, fill= Type)) +
  geom_boxplot()+
  ggtitle("stim RCaMP") +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



# GCaMP amplitude
df.amp1.GC<- summarySE(all.peaks.GC, measurevar = "amplitude", groupvars = c("Condition","Type"))
df.amp3.GC<-summarySE(all.peaks.GC, measurevar = "amplitude", groupvars = c("Type"))

# GCaMP spontaneous only
df.amp1.GC.nostim<- summarySE(nostim.peaks.GC, measurevar = "amplitude", groupvars = c("Type"))

# GCaMP stim only
df.amp1.GC.stim<- summarySE(stim.peaks.GC, measurevar = "amplitude", groupvars = c("Type"))

ggplot(data=df.amp1.GC, aes(x=Type, y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp3.GC, aes(x=Type, y= amplitude, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("GCaMP- all peaks (stim and nostim together") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


ggplot(nostim.peaks.GC, aes(x=Type,y=amplitude, fill= Type)) +
  geom_boxplot()+
  ylab("amplitude") +
  ggtitle("no stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp1.GC.nostim, aes(x=Type, y= amplitude, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("no stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# stim

ggplot(data=df.amp1.GC.stim, aes(x=Type, y= amplitude, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


ggplot(stim.peaks.GC, aes(x=Type,y=amplitude, fill= Type)) +
  geom_boxplot()+
  ylab("stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


#####################
## STATS


#RCaMP nostim
amp.RC.NS.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.RC,REML=FALSE)
amp.RC.NS.model1 = lmer(amplitude ~ Type + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.RC,REML=FALSE)
amp.RC.NS.anova <- anova(amp.RC.NS.null, amp.RC.NS.model1)
print(amp.RC.NS.anova)

amp.RC.NS.Type<- glht(amp.RC.NS.model1, mcp(Type= "Tukey"))
summary(amp.RC.NS.Type)



#RCaMP stim
amp.RC.S.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
amp.RC.S.model1 = lmer(amplitude ~ Type + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
amp.RC.S.anova <- anova(amp.RC.S.null, amp.RC.S.model1)
print(amp.RC.S.anova)

amp.RC.S.Type<- glht(amp.RC.S.model1, mcp(Type= "Tukey"))
summary(amp.RC.S.Type)


# GCaMP all peaks together
amp.GC.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.GC,REML=FALSE)
amp.GC.model1 = lmer(amplitude ~ Type + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.GC,REML=FALSE)
amp.GC.anova <- anova(amp.GC.null, amp.GC.model1)
print(amp.GC.anova)

amp.GC.Type<- glht(amp.GC.model1, mcp(Type= "Tukey"))
summary(amp.GC.Type)

#GCaMP nostim
amp.GC.NS.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.GC,REML=FALSE)
amp.GC.NS.model1 = lmer(amplitude ~ Type + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.GC,REML=FALSE)
amp.GC.NS.anova <- anova(amp.GC.NS.null, amp.GC.NS.model1)
print(amp.GC.NS.anova)

amp.GC.NS.Type<- glht(amp.GC.NS.model1, mcp(Type= "Tukey"))
summary(amp.GC.NS.Type)


#GCaMP stim
amp.GC.S.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
amp.GC.S.model1 = lmer(amplitude ~ Type + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
amp.GC.S.anova <- anova(amp.GC.S.null, amp.GC.S.model1)
print(amp.GC.S.anova)

amp.GC.S.Type<- glht(amp.GC.S.model1, mcp(Type= "Tukey"))
summary(amp.GC.S.Type)



#########
#mean Duration of signals for each ROI

#RCaMP
df.dur1.RC<- summarySE(all.peaks.RC, measurevar = "Duration", groupvars = c("Condition","Type"))

# RCdur spontaneous only
df.dur1.RC.nostim<- summarySE(nostim.peaks.RC, measurevar = "Duration", groupvars = c("Type"))

# RCaMP stim only
df.dur1.RC.stim<- summarySE(stim.peaks.RC, measurevar = "Duration", groupvars = c("Type"))

# RCaMP graphs
ggplot(data=df.dur1.RC, aes(x=Type, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



# no stim
ggplot(data=df.dur1.RC.nostim, aes(x=Type, y= Duration, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("no stim Duration") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme




# stim
ggplot(data=df.dur1.RC.stim, aes(x=Type, y= Duration, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("stim Duration") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.RC, aes(x=Type,y=Duration, fill= Type)) +
  geom_boxplot()+
  ylab("stim RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


# GCaMP Duration
df.dur1.GC<- summarySE(all.peaks.GC, measurevar = "Duration", groupvars = c("Condition","Type"))

# GCaMP spontaneous only
df.dur1.GC.nostim<- summarySE(nostim.peaks.GC, measurevar = "Duration", groupvars = c("Type"))

# GCaMP stim only
df.dur1.GC.stim<- summarySE(stim.peaks.GC, measurevar = "Duration", groupvars = c("Type"))


ggplot(data=df.dur1.GC, aes(x=Type, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



ggplot(nostim.peaks.GC, aes(x=Type,y=Duration, fill= Type)) +
  geom_boxplot()+
  ylab("no stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.dur1.GC.nostim, aes(x=Type, y= Duration, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("no stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


ggplot(data=df.dur1.GC.stim, aes(x=Type, y= Duration, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



ggplot(stim.peaks.GC, aes(x=Type,y=Duration, fill= Type)) +
  geom_boxplot()+
  ylab("stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



#####################
## STATS

#RCaMP nostim
dur.RC.NS.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.RC,REML=FALSE)
dur.RC.NS.model1 = lmer(Duration ~ Type + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.RC,REML=FALSE)
dur.RC.NS.anova <- anova(dur.RC.NS.null, dur.RC.NS.model1)
print(dur.RC.NS.anova)

dur.RC.NS.Type<- glht(dur.RC.NS.model1, mcp(Type= "Tukey"))
summary(dur.RC.NS.Type)


#RCaMP stim
dur.RC.S.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
dur.RC.S.model1 = lmer(Duration ~ Type + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
dur.RC.S.anova <- anova(dur.RC.S.null, dur.RC.S.model1)
print(dur.RC.S.anova)

dur.RC.S.Type<- glht(dur.RC.S.model1, mcp(Type= "Tukey"))
summary(dur.RC.S.Type)



#GCaMP nostim
dur.GC.NS.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.GC,REML=FALSE)
dur.GC.NS.model1 = lmer(Duration ~ Type + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.GC,REML=FALSE)
dur.GC.NS.anova <- anova(dur.GC.NS.null, dur.GC.NS.model1)
print(dur.GC.NS.anova)

dur.GC.NS.Type<- glht(dur.GC.NS.model1, mcp(Type= "Tukey"))
summary(dur.GC.NS.Type)

#GCaMP stim
dur.GC.S.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
dur.GC.S.model1 = lmer(Duration ~ Type + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
dur.GC.S.anova <- anova(dur.GC.S.null, dur.GC.S.model1)
print(dur.GC.S.anova)

dur.GC.S.Type<- glht(dur.GC.S.model1, mcp(Type= "Tukey"))
summary(dur.GC.S.Type)


#################
# PEAK TIME

#RCaMP
ggplot(all.peaks.RC[all.peaks.RC$Condition=="Stim",],aes(x=peakTime,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("all RCaMP stim peak times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

#
all.peaks.RC.PT<-subset(all.peaks.RC, peakTime>0 & peakTime<4)
df.PT1.RC<- summarySE(all.peaks.RC.PT, measurevar = "peakTime", groupvars = c("Condition","Type"))

# RCdur spontaneous only
df.PT1.RC.nostim<- summarySE(all.peaks.RC.PT[all.peaks.RC.PT$Condition=="Nostim",], measurevar = "peakTime", groupvars = c("Type"))

# RCaMP stim only
df.PT1.RC.stim<- summarySE(all.peaks.RC.PT[all.peaks.RC.PT$Condition=="Stim",], measurevar = "peakTime", groupvars = c("Type"))

# RCaMP graphs
ggplot(data=df.PT1.RC, aes(x=Type, y= peakTime, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakTime") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



# no stim
ggplot(data=df.PT1.RC.nostim, aes(x=Type, y= peakTime, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("no stim peakTime") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


# stim
ggplot(data=df.PT1.RC.stim, aes(x=Type, y= peakTime, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("stim peakTime") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.RC, aes(x=Type,y=peakTime, fill= Type)) +
  geom_boxplot()+
  ylab("stim RCaMP peakTime") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


# GCaMP peakTime
df.PT1.GC<- summarySE(all.peaks.GC, measurevar = "peakTime", groupvars = c("Condition","Type"))

# GCaMP spontaneous only
df.PT1.GC.nostim<- summarySE(nostim.peaks.GC, measurevar = "peakTime", groupvars = c("Type"))

# GCaMP stim only
df.PT1.GC.stim<- summarySE(stim.peaks.GC, measurevar = "peakTime", groupvars = c("Type"))


ggplot(data=df.PT1.GC, aes(x=Type, y= peakTime, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakTime") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



ggplot(nostim.peaks.GC, aes(x=Type,y=peakTime, fill= Type)) +
  geom_boxplot()+
  ylab("no stim GCaMP peakTime") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.PT1.GC.nostim, aes(x=Type, y= peakTime, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakTime") +
  ggtitle("no stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


ggplot(data=df.PT1.GC.stim, aes(x=Type, y= peakTime, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=peakTime-se, ymax=peakTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("peakTime") +
  ggtitle("stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



ggplot(stim.peaks.GC, aes(x=Type,y=peakTime, fill= Type)) +
  geom_boxplot()+
  ylab("stim GCaMP peakTime") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



#####################
## STATS

#RCaMP stim
PT.RC.S.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.RC.PT[all.peaks.RC.PT$Condition=="Stim",],REML=FALSE)
PT.RC.S.model1 = lmer(peakTime ~ Type + (1|Animal) + (1|Spot) + (1|ROIs_trial), all.peaks.RC.PT[all.peaks.RC.PT$Condition=="Stim",],REML=FALSE)
PT.RC.S.anova <- anova(PT.RC.S.null, PT.RC.S.model1)
print(PT.RC.S.anova)

PT.RC.S.Type<- glht(PT.RC.S.model1, mcp(Type= "Tukey"))
summary(PT.RC.S.Type)



#GCaMP nostim
PT.GC.NS.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.GC,REML=FALSE)
PT.GC.NS.model1 = lmer(peakTime ~ Type + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.GC,REML=FALSE)
PT.GC.NS.anova <- anova(PT.GC.NS.null, PT.GC.NS.model1)
print(PT.GC.NS.anova)

PT.GC.NS.Type<- glht(PT.GC.NS.model1, mcp(Type= "Tukey"))
summary(PT.GC.NS.Type)

#GCaMP stim
PT.GC.S.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
PT.GC.S.model1 = lmer(peakTime ~ Type + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
PT.GC.S.anova <- anova(PT.GC.S.null, PT.GC.S.model1)
print(PT.GC.S.anova)

PT.GC.S.Type<- glht(PT.GC.S.model1, mcp(Type= "Tukey"))
summary(PT.GC.S.Type)



########################
# ONSET TIME
#remove entries with no onset time

all.OT<-subset(all.OT, OnsetTime!="NaN")

stim.OT.GC<-subset(all.OT, Channel=="GCaMP" & Condition =="Stim")
stim.OT.RC<-subset(all.OT, Channel=="RCaMP" & Condition =="Stim")

######
# neuronal responses to stimulation
NeuronalStim<-subset(all.OT, Channel=="RCaMP" & Condition=="Stim" & OnsetTime<8)

# should have an onset time in 8 s stimulus
ggplot(NeuronalStim[NeuronalStim$Type!="Evening",],aes(x=OnsetTime,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.07, position="dodge") +
  ggtitle("RCaMP onset times between 0 and 8 s from stim trials")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(NeuronalStim[NeuronalStim$Type=="KD",],aes(x=OnsetTime,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP onset times between 0 and 8 s from stim trials")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(NeuronalStim[NeuronalStim$Type=="Control",],aes(x=OnsetTime,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP onset times between 0 and 8 s from stim trials")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# all responding neurons
Neuron95Onset<-quantile(NeuronalStim$OnsetTime, prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
NeuronPT50<-Neuron95Onset[[11]]
print(NeuronPT50)

Neuron95Onset.control<-quantile(NeuronalStim$OnsetTime[NeuronalStim$Type=="Control"], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
NeuronPT50.c<-Neuron95Onset.control[[11]]
print(NeuronPT50.c)

Neuron95Onset.KD<-quantile(NeuronalStim$OnsetTime[NeuronalStim$Type=="KD"], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
NeuronPT50.KD<-Neuron95Onset.KD[[11]]
print(NeuronPT50.KD)

# time thresold to consider an astrocyte to be fast:
#fastTh.c<-NeuronPT50.c
#fastTh.KD<-NeuronPT50.KD

fastTh.c<-1.5


#######
#plot more distributions

# GCaMP onset times
ggplot(stim.OT.GC,aes(x=OnsetTime,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP onset times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme



#####
# onset time boxplots

stim.OT.GC.window<-subset(stim.OT.GC, OnsetTime<12)
stim.OT.RC.window<-subset(stim.OT.RC, OnsetTime<8)

ggplot(stim.OT.GC.window, aes(x=Type,y=OnsetTime, fill= Type)) +
  geom_boxplot()+
  ylab("Onset Latency (s)") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.OT.RC.window[stim.OT.RC.window$Type!="Evening",], aes(x=Type,y=OnsetTime, fill= Type)) +
  geom_boxplot()+
  ylab("Onset Latency (s)") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

wilcox.test(stim.OT.RC.window$OnsetTime[stim.OT.RC.window$Type=="Control"], 
            stim.OT.RC.window$OnsetTime[stim.OT.RC.window$Type=="KD"])

wilcox.test(stim.OT.GC.window$OnsetTime[stim.OT.GC.window$Type=="Control"], 
            stim.OT.GC.window$OnsetTime[stim.OT.GC.window$Type=="KD"])

df.RC.OT<-summarySE(stim.OT.RC.window[stim.OT.RC.window$Type!="Evening",], measurevar = "OnsetTime", groupvars = c("Type"))
df.GC.OT<-summarySE(stim.OT.GC.window[stim.OT.GC.window$Type!="Evening",], measurevar = "OnsetTime", groupvars = c("Type"))


# plots

ggplot(df.RC.OT, aes(x=Type,y=OnsetTime, fill= Type)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Onset Latency (s)") +
  ggtitle("RCaMP ROIs")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(df.GC.OT, aes(x=Type,y=OnsetTime, fill= Type)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=OnsetTime-se, ymax=OnsetTime+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Onset Latency (s)") +
  ggtitle("GCaMP ROIs")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

######
## STATS

#RCaMP onset time
OT.RC.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot), stim.OT.RC.window,REML=FALSE)
OT.RC.model1 = lmer(OnsetTime ~ Type + (1|Animal) + (1|Spot), stim.OT.RC.window,REML=FALSE)
OT.RC.anova <- anova(OT.RC.null, OT.RC.model1)
print(OT.RC.anova)

OT.RC.Type<- glht(OT.RC.model1, mcp(Type= "Tukey"))
summary(OT.RC.Type)

#GCaMP onset time
OT.GC.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot), stim.OT.GC.window,REML=FALSE)
OT.GC.model1 = lmer(OnsetTime ~ Type + (1|Animal) + (1|Spot), stim.OT.GC.window,REML=FALSE)
OT.GC.anova <- anova(OT.GC.null, OT.GC.model1)
print(OT.GC.anova)

OT.GC.Type<- glht(OT.GC.model1, mcp(Type= "Tukey"))
summary(OT.GC.Type)

#########
# find fast and delayed astrocyte microdomains


# identify "FAST" astrocytes
#stim.OT.GC.cont<-subset(stim.OT.GC, Type== "Control")
#stim.OT.GC.cont$Group<-0
#stim.OT.GC.cont$Group[stim.OT.GC.cont$OnsetTime<fastTh.c]<-"fast"
#stim.OT.GC.cont$Group[stim.OT.GC.cont$OnsetTime>=fastTh.c]<-"delayed"

stim.OT.GC.window$Group<-0
stim.OT.GC.window$Group[stim.OT.GC.window$OnsetTime<fastTh.c]<-"fast"
stim.OT.GC.window$Group[stim.OT.GC.window$OnsetTime>=fastTh.c]<-"delayed"


# proportion of fast ROIs

control.fast<-subset(stim.OT.GC.window,Group=="fast" & Type=="Control")
control.del<-subset(stim.OT.GC.window,Group=="delayed" & Type=="Control")

KD.fast<-subset(stim.OT.GC.window,Group=="fast" & Type=="KD")
KD.del<-subset(stim.OT.GC.window,Group=="delayed" & Type=="KD")

stim.peaks.GC.match<-stim.peaks.GC
# add onset time information to the peak data table
stim.peaks.GC.match<-merge(stim.peaks.GC.match, stim.OT.GC.window[, c("ROIs_trial_Cond", "OnsetTime", "Group")], by="ROIs_trial_Cond", all.x=TRUE)


# remove peaks that don't have a corresponding onset time

stim.peaks.GC.match<-stim.peaks.GC.match[!is.na(stim.peaks.GC.match$OnsetTime),]

###################

# number of active ROIs per field of view
stim.peaks.GC.match$Animal_Spot<- paste(stim.peaks.GC.match$Animal, stim.peaks.GC.match$Spot, sep="_")
all.peaks.window$Animal_Spot<- paste(all.peaks.window$Animal, all.peaks.window$Spot, sep="_")

# number of ROIs in each trial for each field of view (across the whole trial) with a peak during no stim and stim
# only consider during the stimulus (and the same time window in no stim trials)

ROInum.8strial<-ddply(all.peaks.window, c("Animal","Spot","Type", "Condition","Channel","Animal_Spot"), summarise, nROIs=length(unique(ROIs_trial_Cond)))

# only GCaMP stim peaks with fast or delayed
ROInum.8strial.group<-ddply(stim.peaks.GC.match, c("Animal","Spot","Type","Condition","Channel","Group"), summarise, nROIs=length(unique(ROIs_trial_Cond)))


# add in number of trials
ROInum.8strial$Ani_Spot_Cond<-paste(ROInum.8strial$Animal_Spot, ROInum.8strial$Condition, sep="_")
ROInum.8strial.group$Ani_Spot_Cond<-paste(ROInum.8strial.group$Animal, ROInum.8strial.group$Spot, ROInum.8strial.group$Condition, sep="_")

ROInum.8strial<-merge(ROInum.8strial, Spot.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)
ROInum.8strial$ROIsPerTrial<-ROInum.8strial$nROIs/ROInum.8strial$nTrials

ROInum.8strial.group<-merge(ROInum.8strial.group, Spot.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)
ROInum.8strial.group$ROIsPerTrial<-ROInum.8strial.group$nROIs/ROInum.8strial.group$nTrials

ROInum.8strial.stim<-subset(ROInum.8strial, Condition=="Stim")
ROInum.8strial.nostim<-subset(ROInum.8strial, Condition=="Nostim")

########
# means
df.ROInum.8strial.RC1<-summarySE(ROInum.8strial[ROInum.8strial$Channel=="RCaMP",], measurevar = "ROIsPerTrial", groupvars = c("Condition","Type"))
df.ROInum.8strial.RC2<-summarySE(ROInum.8strial[ROInum.8strial$Channel=="RCaMP",], measurevar = "ROIsPerTrial", groupvars = c("Type"))

df.ROInum.8strial.NS.RC1<-summarySE(ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="RCaMP",], measurevar = "ROIsPerTrial", groupvars = c("Type"))
df.ROInum.8strial.S.RC1<-summarySE(ROInum.8strial.stim[ROInum.8strial.stim$Channel=="RCaMP",], measurevar = "ROIsPerTrial", groupvars = c("Type"))


df.ROInum.8strial.GC1<-summarySE(ROInum.8strial[ROInum.8strial$Channel=="GCaMP",], measurevar = "ROIsPerTrial", groupvars = c("Condition","Type"))

df.ROInum.8strial.NS.GC1<-summarySE(ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="GCaMP",], measurevar = "ROIsPerTrial", groupvars = c("Type"))
df.ROInum.8strial.S.GC1<-summarySE(ROInum.8strial.stim[ROInum.8strial.stim$Channel=="GCaMP",], measurevar = "ROIsPerTrial", groupvars = c("Type"))

# plots

# individual plots
# RCaMP

ggplot(ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="RCaMP",], aes(x=Type,y=ROIsPerTrial, fill= Type)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("RCaMP nostim") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


ggplot(ROInum.8strial.stim[ROInum.8strial.stim$Channel=="RCaMP",], aes(x=Type,y=ROIsPerTrial, fill= Type)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("RCaMP stim") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.ROInum.8strial.RC2[df.ROInum.8strial.RC2$Type!="Evening",], aes(x=Type, y= ROIsPerTrial, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIsPerTrial") +
  ggtitle("RCaMP- all ROIs together") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.ROInum.8strial.RC1, aes(x=Type, y= ROIsPerTrial, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIsPerTrial") +
  ggtitle("RCaMP- all ROIs together") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# GCaMP
ggplot(ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="GCaMP",], aes(x=Type,y=ROIsPerTrial, fill= Type)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("GCaMP nostim") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(ROInum.8strial.stim[ROInum.8strial.stim$Channel=="GCaMP",], aes(x=Type,y=ROIsPerTrial, fill= Type)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("GCaMP stim") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.ROInum.8strial.NS.GC1, aes(x=Type, y= ROIsPerTrial, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIsPerTrial") +
  ggtitle("nostim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.ROInum.8strial.S.GC1, aes(x=Type, y= ROIsPerTrial, fill=Type)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIsPerTrial") +
  ggtitle("stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


ggplot(data=df.ROInum.8strial.GC1, aes(x=Type, y= ROIsPerTrial, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ROIsPerTrial") +
  ggtitle("all GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

###############
# stats for active ROI number per trials per FOV

# RCaMP nostim
ROInum.RC.NS.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="RCaMP",],REML=FALSE)
ROInum.RC.NS.model1 = lmer(ROIsPerTrial ~ Type + (1|Animal), ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="RCaMP",],REML=FALSE)
ROInum.RC.NS.model2 = lmer(ROIsPerTrial ~ shRNA2 + (1|Animal), ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="RCaMP",],REML=FALSE)
ROInum.RC.NS.anova <- anova(ROInum.RC.NS.null, ROInum.RC.NS.model1, ROInum.RC.NS.model2)
print(ROInum.RC.NS.anova)

ROInum.RC.NS.Type<- glht(ROInum.RC.NS.model1, mcp(Type= "Tukey"))
summary(ROInum.RC.NS.Type)

ROInum.RC.NS.shRNA2<- glht(ROInum.RC.NS.model2, mcp(shRNA2= "Tukey"))
summary(ROInum.RC.NS.shRNA2)

# RCaMP stim
ROInum.RC.S.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.8strial.nostim[ROInum.8strial.stim$Channel=="RCaMP",],REML=FALSE)
ROInum.RC.S.model1 = lmer(ROIsPerTrial ~ Type + (1|Animal), ROInum.8strial.nostim[ROInum.8strial.stim$Channel=="RCaMP",],REML=FALSE)
ROInum.RC.S.model2 = lmer(ROIsPerTrial ~ shRNA2 + (1|Animal), ROInum.8strial.nostim[ROInum.8strial.stim$Channel=="RCaMP",],REML=FALSE)
ROInum.RC.S.anova <- anova(ROInum.RC.S.null, ROInum.RC.S.model1, ROInum.RC.S.model2)
print(ROInum.RC.S.anova)

ROInum.RC.S.Type<- glht(ROInum.RC.S.model1, mcp(Type= "Tukey"))
summary(ROInum.RC.S.Type)

ROInum.RC.S.shRNA2<- glht(ROInum.RC.S.model2, mcp(shRNA2= "Tukey"))
summary(ROInum.RC.S.shRNA2)

# GCaMP nostim
ROInum.GC.NS.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="GCaMP",],REML=FALSE)
ROInum.GC.NS.model1 = lmer(ROIsPerTrial ~ Type + (1|Animal), ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="GCaMP",],REML=FALSE)
ROInum.GC.NS.model2 = lmer(ROIsPerTrial ~ shRNA2 + (1|Animal), ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="GCaMP",],REML=FALSE)
ROInum.GC.NS.anova <- anova(ROInum.GC.NS.null, ROInum.GC.NS.model1, ROInum.GC.NS.model2)
print(ROInum.GC.NS.anova)

ROInum.GC.NS.Type<- glht(ROInum.GC.NS.model1, mcp(Type= "Tukey"))
summary(ROInum.GC.NS.Type)

ROInum.GC.NS.shRNA2<- glht(ROInum.GC.NS.model2, mcp(shRNA2= "Tukey"))
summary(ROInum.GC.NS.shRNA2)

ROInum.8strial.stim$shRNA_Timepoint=interaction(ROInum.8strial.stim$Type,ROInum.8strial.stim$Timepoint)
# GCaMP stim
ROInum.GC.S.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.8strial.stim[ROInum.8strial.stim$Channel=="GCaMP",],REML=FALSE)
ROInum.GC.S.model1 = lmer(ROIsPerTrial ~ Type + (1|Animal), ROInum.8strial.stim[ROInum.8strial.stim$Channel=="GCaMP",],REML=FALSE)
ROInum.GC.S.model2 = lmer(ROIsPerTrial ~ shRNA2 + (1|Animal), ROInum.8strial.stim[ROInum.8strial.stim$Channel=="GCaMP",],REML=FALSE)
ROInum.GC.S.model3 = lmer(ROIsPerTrial ~ shRNA_Timepoint + (1|Animal), ROInum.8strial.stim[ROInum.8strial.stim$Channel=="GCaMP",],REML=FALSE)
ROInum.GC.S.anova <- anova(ROInum.GC.S.null, ROInum.GC.S.model1, ROInum.GC.S.model2,ROInum.GC.S.model3)
print(ROInum.GC.S.anova)

ROInum.GC.S.Type<- glht(ROInum.GC.S.model1, mcp(Type= "Tukey"))
summary(ROInum.GC.S.Type)

ROInum.GC.S.shRNA2<- glht(ROInum.GC.S.model2, mcp(shRNA2= "Tukey"))
summary(ROInum.GC.S.shRNA2)

ROInum.GC.S.TP<- glht(ROInum.GC.S.model3, mcp(shRNA_Timepoint= "Tukey"))
summary(ROInum.GC.S.TP)


######
# RCaMP matched peaks and onset times
# add onset time information to the peak data table
stim.peaks.RC<-merge(stim.peaks.RC, stim.OT.RC.window[, c("ROIs_trial_Cond", "OnsetTime")], by="ROIs_trial_Cond", all.x=TRUE)


# remove peaks that don't have a corresponding onset time

stim.peaks.RC<-stim.peaks.RC[!is.na(stim.peaks.RC$OnsetTime),]

ggplot(stim.peaks.RC, aes(x=OnsetTime,y=amplitude, colour= Type)) +
  geom_point(stat="identity") +
  ggtitle("RCaMP peaks onset time vs amplitude")+
  scale_colour_manual(values=cbbPalette) + 
  max.theme

#high, mid, and low responding neurons based on control population

NeuronalAmps<-quantile(stim.peaks.RC$amplitude[stim.peaks.RC$Type=="Control"], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)

# high = 95% percentile, mid = 50-95, low = <50
high<-NeuronalAmps[[20]]
low<-NeuronalAmps[[11]]

stim.peaks.RC$Neur_type<-"mid"
stim.peaks.RC$Neur_type[stim.peaks.RC$amplitude>high]<-"high"
stim.peaks.RC$Neur_type[stim.peaks.RC$amplitude<low]<-"low"

stim.peaks.RC<-subset(stim.peaks.RC, Type!="Evening")

#amplitudes
df.neur.amp<-summarySE(stim.peaks.RC, measurevar = "amplitude", groupvars = c("Neur_type", "Type"))

# onset times less than 2 s
df.neur.OT<-summarySE(stim.peaks.RC[stim.peaks.RC$OnsetTime<1.5,], measurevar = "OnsetTime", groupvars = c("Neur_type", "Type"))

df.neur.AUC<-summarySE(stim.peaks.RC, measurevar = "peakAUC", groupvars = c("Neur_type", "Type"))


ggplot(stim.peaks.RC[stim.peaks.RC$Type!="Evening" & stim.peaks.RC$OnsetTime<3,],aes(x=OnsetTime,y=..density..,fill=Type)) +
  geom_histogram(binwidth=0.06, position="dodge") +
  ggtitle("all RCaMP stim onset times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.RC[stim.peaks.RC$Type!="Evening" & stim.peaks.RC$OnsetTime<3,],aes(x=OnsetTime,y=..density..,fill=Type)) +
  geom_density(adjust = 1/3,alpha = 0.5)+
  ggtitle("all RCaMP stim onset times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

## median onset time in the 1.5s after stim start
NeuronalOT.c<-quantile(stim.peaks.RC$OnsetTime[stim.peaks.RC$Type=="Control" & stim.peaks.RC$OnsetTime<3], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
medianOT.c<-NeuronalOT.c[[11]]

NeuronalOT.KD<-quantile(stim.peaks.RC$OnsetTime[stim.peaks.RC$Type=="KD" & stim.peaks.RC$OnsetTime<3], prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
medianOT.KD<-NeuronalOT.KD[[11]]

#######
# LOOK AT TRIAL VARIATION- CHANGE IN RESPONSE WITH STIM
# look at neuronal frequencies- multiple peaks per cell
# look at pairing fast and delayed so there are both for each spot

###############
# GCaMP fast vs delayed

df.ROInum.8strial.group1<-summarySE(ROInum.8strial.group, measurevar = "ROIsPerTrial", groupvars = c("Type","Group"))



# plots


ggplot(df.ROInum.8strial.group1, aes(x=Group,y=ROIsPerTrial, fill= Type)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("GCaMP ROIs")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

#boxplots
ggplot(ROInum.8strial.group, aes(x=Group,y=ROIsPerTrial, fill= Type)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



# only fast GcaMP

fastROInum<- subset(ROInum.8strial.group, Group=="fast")

ggplot(fastROInum, aes(x=Type,y=ROIsPerTrial, fill= Type)) +
  geom_boxplot()+
  ylab("Lck- fast ROIs per FOV") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# only delayed GcaMP

delayedROInum<- subset(ROInum.8strial.group, Group=="delayed")

ggplot(delayedROInum, aes(x=Type,y=ROIsPerTrial, fill= Type)) +
  geom_boxplot()+
  ylab("Lck- delayed ROIs per FOV") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(df.ROInum.8strial.group1[df.ROInum.8strial.group1$Group=="delayed",], aes(x=Type,y=ROIsPerTrial, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("delayed GCaMP ROIs") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

######
# stats

# fast GCaMP only

nROI.fast.day.null = lmer(ROIsPerTrial ~ (1|Animal), fastROInum,REML=FALSE)
nROI.fast.day.model1 = lmer(ROIsPerTrial ~ Type + (1|Animal), fastROInum,REML=FALSE)
nROI.fast.day.anova <- anova(nROI.fast.day.null, nROI.fast.day.model1)
print(nROI.fast.day.anova)

nROI.GC.Type.fast<- glht(nROI.fast.day.model1, mcp(Type= "Tukey"))
summary(nROI.GC.Type.fast)


# delayed GCaMP only

nROI.delayed.day.null = lmer(ROIsPerTrial ~ (1|Animal), delayedROInum,REML=FALSE)
nROI.delayed.day.model1 = lmer(ROIsPerTrial ~ Type + (1|Animal), delayedROInum,REML=FALSE)
nROI.delayed.day.anova <- anova(nROI.delayed.day.null, nROI.delayed.day.model1)
print(nROI.delayed.day.anova)

nROI.GC.Type.delayed<- glht(nROI.delayed.day.model1, mcp(Type= "Tukey"))
summary(nROI.GC.Type.delayed)



# wilcox test for median compairsons

wilcox.test(fastROInum$ROIsPerTrial[fastROInum$Type=="Control"], 
            fastROInum$ROIsPerTrial[fastROInum$Type=="KD"])

wilcox.test(fastROInum$ROIsPerTrial[fastROInum$shRNA2=="Control"], 
            fastROInum$ROIsPerTrial[fastROInum$shRNA2=="KD"])


#fast vs delayed properties

df.fast.amp1<-summarySE(stim.peaks.GC, measurevar = "amplitude", groupvars = c("Type","Group"))
df.fast.dur1<-summarySE(stim.peaks.GC, measurevar = "Duration", groupvars = c("Type","Group"))


ggplot(stim.peaks.GC, aes(x=Group, y=amplitude, fill= Type)) +
  geom_boxplot()+
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(df.fast.amp1, aes(x=Group,y=amplitude, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("GcaMP signals") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.GC, aes(x=Group, y=Duration, fill= Type)) +
  geom_boxplot()+
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(df.fast.dur1, aes(x=Group,y=Duration, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("GcaMP signals") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

Type_group<-interaction(stim.peaks.GC$Group, stim.peaks.GC$Type)
#GCaMP stim
amp2.GC.S.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
amp2.GC.S.model1 = lmer(amplitude ~ Type + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
amp2.GC.S.model2 = lmer(amplitude ~ shRNA2 + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
amp2.GC.S.model3 = lmer(amplitude ~ Type_group + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
amp2.GC.S.anova <- anova(amp2.GC.S.null, amp2.GC.S.model1, amp2.GC.S.model2, amp2.GC.S.model3)
print(amp2.GC.S.anova)

amp2.GC.S.fast <- glht(amp2.GC.S.model3 , mcp(Type_group= "Tukey"))
summary(amp2.GC.S.fast)


dur2.GC.S.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
dur2.GC.S.model1 = lmer(Duration ~ Type + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
dur2.GC.S.model2 = lmer(Duration ~ shRNA2 + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
dur2.GC.S.model3 = lmer(Duration ~ Type_group + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
dur2.GC.S.anova <- anova(dur2.GC.S.null, dur2.GC.S.model1, dur2.GC.S.model2, dur2.GC.S.model3)
print(dur2.GC.S.anova)

dur2.GC.S.fast <- glht(dur2.GC.S.model3 , mcp(Type_group= "Tukey"))
summary(dur2.GC.S.fast)


