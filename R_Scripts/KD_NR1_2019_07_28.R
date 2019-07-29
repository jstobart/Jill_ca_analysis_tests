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

# Jill's home files
# peak data
long.peaks.74<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/74_peaks_longtrials_clean.csv",  header=TRUE, sep = ",")
long.peaks.92<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/92_peaks_longtrials_clean.csv",  header=TRUE, sep = ",")
long.peaks.94<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/94_peaks_longtrials_clean.csv",  header=TRUE, sep = ",")
long.peaks.95<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/95_peaks_longtrials_clean.csv",  header=TRUE, sep = ",")
long.peaks.96<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/96_peaks_longtrials_clean.csv",  header=TRUE, sep = ",")
long.peaks.Alice<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/alice_peaks_longtrials_clean.csv",  header=TRUE, sep = ",")
long.peaks.crazy8<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/crazy8_peaks_longtrials_clean.csv",  header=TRUE, sep = ",")

# OT data
long.OT.74<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/74_onset_time_longtrials_clean.csv",  header=TRUE, sep = ",")
long.OT.92<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/92_onset_time_longtrials_clean.csv",  header=TRUE, sep = ",")
long.OT.94<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/94_onset_time_longtrials_clean.csv",  header=TRUE, sep = ",")
long.OT.95<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/95_onset_time_longtrials_clean.csv",  header=TRUE, sep = ",")
long.OT.96<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/96_onset_time_longtrials_clean.csv",  header=TRUE, sep = ",")
long.OT.Alice<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/alice_onset_time_longtrials_clean.csv",  header=TRUE, sep = ",")
long.OT.crazy8<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/crazy8_onset_time_longtrials_clean.csv",  header=TRUE, sep = ",")

#field data
long.field.74<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/74_Lck_field_longtrials_clean.csv",  header=TRUE, sep = ",")
long.field.92<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/92_Lck_field_longtrials_clean.csv",  header=TRUE, sep = ",")
long.field.94<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/94_Lck_field_longtrials_clean.csv",  header=TRUE, sep = ",")
long.field.95<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/95_Lck_field_longtrials_clean.csv",  header=TRUE, sep = ",")
long.field.96<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/96_Lck_field_longtrials_clean.csv",  header=TRUE, sep = ",")
long.field.Alice<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/alice_Lck_field_longtrials_clean.csv",  header=TRUE, sep = ",")
long.field.crazy8<-read.table("D:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/crazy8_Lck_field_longtrials_clean.csv",  header=TRUE, sep = ",")

###############
# Jill's work files
# peak data
long.peaks.74<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/74_peaks_longtrials_clean.csv",  header=TRUE, sep = ",")
long.peaks.92<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/92_peaks_longtrials_clean.csv",  header=TRUE, sep = ",")
long.peaks.94<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/94_peaks_longtrials_clean.csv",  header=TRUE, sep = ",")
long.peaks.95<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/95_peaks_longtrials_clean.csv",  header=TRUE, sep = ",")
long.peaks.96<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/96_peaks_longtrials_clean.csv",  header=TRUE, sep = ",")
long.peaks.Alice<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/alice_peaks_longtrials_clean.csv",  header=TRUE, sep = ",")
long.peaks.crazy8<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/crazy8_peaks_longtrials_clean.csv",  header=TRUE, sep = ",")

# OT data
long.OT.74<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/74_onset_time_longtrials_clean.csv",  header=TRUE, sep = ",")
long.OT.92<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/92_onset_time_longtrials_clean.csv",  header=TRUE, sep = ",")
long.OT.94<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/94_onset_time_longtrials_clean.csv",  header=TRUE, sep = ",")
long.OT.95<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/95_onset_time_longtrials_clean.csv",  header=TRUE, sep = ",")
long.OT.96<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/96_onset_time_longtrials_clean.csv",  header=TRUE, sep = ",")
long.OT.Alice<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/alice_onset_time_longtrials_clean.csv",  header=TRUE, sep = ",")
long.OT.crazy8<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/crazy8_onset_time_longtrials_clean.csv",  header=TRUE, sep = ",")

#field data
long.field.74<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/74_Lck_field_longtrials_clean.csv",  header=TRUE, sep = ",")
long.field.92<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/92_Lck_field_longtrials_clean.csv",  header=TRUE, sep = ",")
long.field.94<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/94_Lck_field_longtrials_clean.csv",  header=TRUE, sep = ",")
long.field.95<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/95_Lck_field_longtrials_clean.csv",  header=TRUE, sep = ",")
long.field.96<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/96_Lck_field_longtrials_clean.csv",  header=TRUE, sep = ",")
long.field.Alice<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/alice_Lck_field_longtrials_clean.csv",  header=TRUE, sep = ",")
long.field.crazy8<-read.table("G:/Data/GCaMP_RCaMP/NR1_KD/Results/FilesforR/clean_data/crazy8_Lck_field_longtrials_clean.csv",  header=TRUE, sep = ",")


##########
#remove data from 07_02 for 74
long.peaks.74<-long.peaks.74[!(grepl("2019_07_02",long.peaks.74$Spot)),] 
long.OT.74<-long.OT.74[!(grepl("2019_07_02",long.OT.74$Spot)),] 
long.field.74<-long.field.74[!(grepl("2019_07_02",long.field.74$Spot)),] 


# other questionable data (these spots may not have had responding neurons)
long.peaks.74<-long.peaks.74[!(grepl("2019_07_03",long.peaks.74$Spot) & grepl("spot1",long.peaks.74$Spot)),] 
long.OT.74<-long.OT.74[!(grepl("2019_07_03",long.OT.74$Spot) & grepl("spot1",long.OT.74$Spot)),] 
long.field.74<-long.field.74[!(grepl("2019_07_03",long.field.74$Spot) & grepl("spot1",long.field.74$Spot)),] 


long.peaks.92<-long.peaks.92[!(grepl("2019_06_22",long.peaks.92$Spot) & grepl("spot1",long.peaks.92$Spot)),] 
long.OT.92<-long.OT.92[!(grepl("2019_06_22",long.OT.92$Spot) & grepl("spot1",long.OT.92$Spot)),] 
long.field.92<-long.field.92[!(grepl("2019_06_22",long.field.92$Spot) & grepl("spot1",long.field.92$Spot)),] 


long.peaks.95<-long.peaks.95[!(grepl("2019_06_15",long.peaks.95$Spot)),] 
long.OT.95<-long.OT.95[!(grepl("2019_06_15",long.OT.95$Spot)),] 
long.field.95<-long.field.95[!(grepl("2019_06_15",long.field.95$Spot)),] 

long.peaks.96<-long.peaks.96[!(grepl("2019_06_15",long.peaks.96$Spot)),] 
long.OT.96<-long.OT.96[!(grepl("2019_06_15",long.OT.96$Spot)),] 
long.field.96<-long.field.96[!(grepl("2019_06_15",long.field.96$Spot)),] 

long.peaks.94<-long.peaks.94[!(grepl("2019_06_15",long.peaks.94$Spot)),] 
long.OT.94<-long.OT.94[!(grepl("2019_06_15",long.OT.94$Spot)),] 
long.field.94<-long.field.94[!(grepl("2019_06_15",long.field.94$Spot)),] 

long.peaks.94<-long.peaks.94[!(grepl("2019_06_28",long.peaks.94$Spot)),] 
long.OT.94<-long.OT.94[!(grepl("2019_06_28",long.OT.94$Spot)),] 
long.field.94<-long.field.94[!(grepl("2019_06_28",long.field.94$Spot)),] 

long.peaks.94<-long.peaks.94[!(grepl("2019_06_22",long.peaks.94$Spot) & grepl("spot2",long.peaks.94$Spot)),] 
long.OT.94<-long.OT.94[!(grepl("2019_06_22",long.OT.94$Spot) & grepl("spot2",long.OT.94$Spot)),] 
long.field.94<-long.field.94[!(grepl("2019_06_22",long.field.94$Spot) & grepl("spot2",long.field.94$Spot)),] 


##########

#lsm.options(pbkrtest.limit = 100000)
framerate=12.7917

#  merge data
all.OT<-rbind(long.OT.74, long.OT.92, long.OT.94,long.OT.95, long.OT.96, long.OT.Alice, long.OT.crazy8)
all.peaks<-rbind(long.peaks.74, long.peaks.92, long.peaks.94,long.peaks.95, long.peaks.96, long.peaks.Alice, long.peaks.crazy8)
all.field<-rbind(long.field.74, long.field.92, long.field.94,long.field.95, long.field.96, long.field.Alice, long.field.crazy8)

all.peaks$Duration<- all.peaks$halfWidth*2
#add baseline time to peaks table
all.peaks$BL_time<-5

# adjust peak time and duration
all.peaks$peakTime<- all.peaks$peakTime-all.peaks$BL_time  # now time= 0s is the start of stimulation
all.peaks$peakStart<- all.peaks$peakStart-all.peaks$BL_time
all.peaks$peakStartHalf<- all.peaks$peakStartHalf-all.peaks$BL_time


#define shRNA for onset times
# 1- Control including Alice and Crazy8, 2- Alice and Crazy 8 are separate from 74 and 92
all.OT$shRNA1="Control"
all.OT$shRNA2="Control"

all.OT$shRNA1[grepl("94",all.OT$Animal)]="KD"
all.OT$shRNA1[grepl("95",all.OT$Animal)]="KD"
all.OT$shRNA1[grepl("96",all.OT$Animal)]="KD"
all.OT$shRNA2[grepl("94",all.OT$Animal)]="KD"
all.OT$shRNA2[grepl("95",all.OT$Animal)]="KD"
all.OT$shRNA2[grepl("96",all.OT$Animal)]="KD"
all.OT$shRNA2[grepl("74",all.OT$Animal)]="NS"
all.OT$shRNA2[grepl("92",all.OT$Animal)]="NS"

#define shRNA for peaks
# 1- Control including Alice and Crazy8, 2- Alice and Crazy 8 are separate from 74 and 92
all.peaks$shRNA1="Control"
all.peaks$shRNA2="Control"

all.peaks$shRNA1[grepl("94",all.peaks$Animal)]="KD"
all.peaks$shRNA1[grepl("95",all.peaks$Animal)]="KD"
all.peaks$shRNA1[grepl("96",all.peaks$Animal)]="KD"
all.peaks$shRNA2[grepl("94",all.peaks$Animal)]="KD"
all.peaks$shRNA2[grepl("95",all.peaks$Animal)]="KD"
all.peaks$shRNA2[grepl("96",all.peaks$Animal)]="KD"
all.peaks$shRNA2[grepl("74",all.peaks$Animal)]="NS"
all.peaks$shRNA2[grepl("92",all.peaks$Animal)]="NS"

#define shRNA for onset times
# 1- Control including Alice and Crazy8, 2- Alice and Crazy 8 are separate from 74 and 92
all.field$shRNA1="Control"
all.field$shRNA2="Control"

all.field$shRNA1[grepl("94",all.field$animalname)]="KD"
all.field$shRNA1[grepl("95",all.field$animalname)]="KD"
all.field$shRNA1[grepl("96",all.field$animalname)]="KD"
all.field$shRNA2[grepl("94",all.field$animalname)]="KD"
all.field$shRNA2[grepl("95",all.field$animalname)]="KD"
all.field$shRNA2[grepl("96",all.field$animalname)]="KD"
all.field$shRNA2[grepl("74",all.field$animalname)]="NS"
all.field$shRNA2[grepl("92",all.field$animalname)]="NS"

# set these new variables as factors so we can do stats on them
all.peaks$shRNA1<-as.factor(all.peaks$shRNA1)
all.peaks$shRNA2<-as.factor(all.peaks$shRNA2)
all.OT$shRNA1<-as.factor(all.OT$shRNA1)
all.OT$shRNA2<-as.factor(all.OT$shRNA2)
all.field$shRNA1<-as.factor(all.field$shRNA1)
all.field$shRNA2<-as.factor(all.field$shRNA2)

######
# days post injection
all.OT$DaysPostInjection="47"

all.OT$DaysPostInjection[grepl("06_15",all.OT$Spot)]="28"
all.OT$DaysPostInjection[grepl("06_21",all.OT$Spot)]="35"
all.OT$DaysPostInjection[grepl("06_22",all.OT$Spot)]="35"
all.OT$DaysPostInjection[grepl("06_28",all.OT$Spot)]="42"
all.OT$DaysPostInjection[grepl("07_02",all.OT$Spot)]="46"

# early (28, 35) or late (42,46,47) days post injection
all.OT$Timepoint<-"late"
all.OT$Timepoint[grepl("06_15",all.OT$Spot)]="early"
all.OT$Timepoint[grepl("06_21",all.OT$Spot)]="early"
all.OT$Timepoint[grepl("06_22",all.OT$Spot)]="early"
all.OT$Timepoint[grepl("06_28",all.OT$Spot)]="mid"

#peaks
all.peaks$DaysPostInjection="47"

all.peaks$DaysPostInjection[grepl("06_15",all.peaks$Spot)]="28"
all.peaks$DaysPostInjection[grepl("06_21",all.peaks$Spot)]="35"
all.peaks$DaysPostInjection[grepl("06_22",all.peaks$Spot)]="35"
all.peaks$DaysPostInjection[grepl("06_28",all.peaks$Spot)]="42"
all.peaks$DaysPostInjection[grepl("07_02",all.peaks$Spot)]="46"

# early (28, 35) or late (42,46,47) days post injection
all.peaks$Timepoint<-"late"
all.peaks$Timepoint[grepl("06_15",all.peaks$Spot)]="early"
all.peaks$Timepoint[grepl("06_21",all.peaks$Spot)]="early"
all.peaks$Timepoint[grepl("06_22",all.peaks$Spot)]="early"
all.peaks$Timepoint[grepl("06_28",all.peaks$Spot)]="mid"

#field
all.field$DaysPostInjection="47"

all.field$DaysPostInjection[grepl("06_15",all.field$Spot)]="28"
all.field$DaysPostInjection[grepl("06_21",all.field$Spot)]="35"
all.field$DaysPostInjection[grepl("06_22",all.field$Spot)]="35"
all.field$DaysPostInjection[grepl("06_28",all.field$Spot)]="42"
all.field$DaysPostInjection[grepl("07_02",all.field$Spot)]="46"

# early (28, 35) or late (42,46,47) days post injection
all.field$Timepoint<-"late"
all.field$Timepoint[grepl("06_15",all.field$Spot)]="early"
all.field$Timepoint[grepl("06_21",all.field$Spot)]="early"
all.field$Timepoint[grepl("06_22",all.field$Spot)]="early"
all.field$Timepoint[grepl("06_28",all.field$Spot)]="mid"

# set these new variables as factors so we can do stats on them
all.peaks$DaysPostInjection<-as.factor(all.peaks$DaysPostInjection)
all.OT$DaysPostInjection<-as.factor(all.OT$DaysPostInjection)
all.field$DaysPostInjection<-as.factor(all.field$DaysPostInjection)

all.peaks$Timepoint<-as.factor(all.peaks$Timepoint)
all.OT$Timepoint<-as.factor(all.OT$Timepoint)
all.field$Timepoint<-as.factor(all.field$Timepoint)

# set the order of groups for plots
all.OT$shRNA1<-factor(all.OT$shRNA1,levels=c("Control","KD"))
all.OT$shRNA2<-factor(all.OT$shRNA2,levels=c("Control","NS", "KD"))
all.OT$Timepoint<-factor(all.OT$Timepoint,levels=c("early","mid", "late"))

all.peaks$shRNA1<-factor(all.peaks$shRNA1,levels=c("Control","KD"))
all.peaks$shRNA2<-factor(all.peaks$shRNA2,levels=c("Control","NS", "KD"))
all.peaks$Timepoint<-factor(all.peaks$Timepoint,levels=c("early","mid", "late"))

all.field$shRNA1<-factor(all.field$shRNA1,levels=c("Control","KD"))
all.field$shRNA2<-factor(all.field$shRNA2,levels=c("Control","NS", "KD"))
all.field$Timepoint<-factor(all.field$Timepoint,levels=c("early","mid", "late"))

#unique ROI names
all.OT$ROIs_trial<-paste(all.OT$Animal, all.OT$Spot, all.OT$Trial,all.OT$ROI, sep= "_")
all.OT$ROIs_trial_Cond<-paste(all.OT$ROIs_trial, all.OT$Condition, sep= "_")

all.peaks$ROIs_trial<-paste(all.peaks$Animal, all.peaks$Spot, all.peaks$Trial,all.peaks$roiName, sep= "_")
all.peaks$ROIs_trial_Cond<-paste(all.peaks$ROIs_trial, all.peaks$Condition, sep= "_")


all.field$Spot_trial<-paste(all.field$animalname, all.field$Spot, all.field$trialname, sep= "_")
all.field$Spot_trial_Cond<-paste(all.field$Spot_trial, all.field$Cond, sep= "_")


# count number of trials per spot
Spot.ntrials<-ddply(all.peaks, c("Animal","shRNA1","Spot","Condition"), summarise, nTrials=length(unique(Trial)))
Spot.ntrials$Ani_Spot_Cond<-paste(Spot.ntrials$Animal, Spot.ntrials$Spot, Spot.ntrials$Condition, sep="_")


# remove matching dendrite and neuronal soma ROIs
Overlap= all.peaks$overlap!=0
all.peaks<-all.peaks[!Overlap,]

#consider onsly positive amplitude and accurate time scale
all.peaks.window<-subset(all.peaks, peakTime>0 & peakTime<12 & amplitude>0 & Duration<20)

##########################
# FIELD DATA

# remove random trials of NaNs
all.field<-all.field[!is.na(all.field$pixelsize),]

######
# fraction of active pixels from all the astrocyte pixels
# with peaks near stimulus
#all.field$nFluoPix[all.field$nFluoPix==0]<-128*128  #adjust spots where no pixels were above the threshold

all.field$FracActive=all.field$nActivePix/all.field$nFluoPix
all.field$FracFluo=all.field$nFluoPix/all.field$nTotalPix

all.field.spot<-ddply(all.field, c("Spot","animalname", "Cond", "shRNA1", "shRNA2", "DaysPostInjection","Timepoint", "Response_Score"), summarise, 
                      meanFracActive=mean(FracActive), meanFluoPix = mean(nFluoPix), meanFracFluo = mean(FracFluo))
nostim.field.spot<-subset(all.field.spot, Cond=="nostim")
stim.field.spot<-subset(all.field.spot, Cond=="stim")

df.FracActive<- summarySE(all.field.spot, measurevar = "meanFracActive", groupvars = c("Cond","shRNA1"))
df.FracActive2<- summarySE(all.field.spot, measurevar = "meanFracActive", groupvars = c("Cond","shRNA2"))
df.FracActive.TP2<-summarySE(all.field.spot, measurevar = "meanFracActive", groupvars = c("Cond","shRNA2", "Timepoint" ))
df.FracActive.TP<-summarySE(all.field.spot, measurevar = "meanFracActive", groupvars = c("Cond","shRNA1", "Timepoint" ))
df.FracActive.nostim1<-summarySE(nostim.field.spot, measurevar = "meanFracActive", groupvars = c("shRNA1"))
df.FracActive.nostim2<-summarySE(nostim.field.spot, measurevar = "meanFracActive", groupvars = c("shRNA2"))
df.FracActive.nostim3<-summarySE(nostim.field.spot, measurevar = "meanFracActive", groupvars = c("shRNA1","Timepoint"))
df.FracActive.stim1<-summarySE(stim.field.spot, measurevar = "meanFracActive", groupvars = c("shRNA1"))
df.FracActive.stim2<-summarySE(stim.field.spot, measurevar = "meanFracActive", groupvars = c("shRNA2"))
df.FracActive.stim3<-summarySE(stim.field.spot, measurevar = "meanFracActive", groupvars = c("shRNA1", "Timepoint"))

df.ResponseScore<- summarySE(all.field.spot, measurevar = "Response_Score", groupvars = c("Cond","shRNA1"))
df.ResponseScore2<- summarySE(all.field.spot, measurevar = "Response_Score", groupvars = c("Cond","shRNA2"))
df.ResponseScore.TP2<- summarySE(all.field.spot, measurevar = "Response_Score", groupvars = c("Cond","shRNA2","Timepoint"))
df.ResponseScore.TP<- summarySE(all.field.spot, measurevar = "Response_Score", groupvars = c("Cond","shRNA1","Timepoint"))

df.nFluoPix.nostim1<-summarySE(nostim.field.spot, measurevar = "meanFluoPix", groupvars = c("shRNA1"))
df.nFluoPix.nostim2<-summarySE(nostim.field.spot, measurevar = "meanFluoPix", groupvars = c("shRNA2"))
df.nFluoPix.nostim3<-summarySE(nostim.field.spot, measurevar = "meanFluoPix", groupvars = c("shRNA1","Timepoint"))
df.nFluoPix.stim1<-summarySE(stim.field.spot, measurevar = "meanFluoPix", groupvars = c("shRNA1"))
df.nFluoPix.stim2<-summarySE(stim.field.spot, measurevar = "meanFluoPix", groupvars = c("shRNA2"))
df.nFluoPix.stim3<-summarySE(stim.field.spot, measurevar = "meanFluoPix", groupvars = c("shRNA1","Timepoint"))

df.FracFluo.nostim2<-summarySE(nostim.field.spot, measurevar = "meanFracFluo", groupvars = c("shRNA2"))
df.FracFluo.nostim3<-summarySE(nostim.field.spot, measurevar = "meanFracFluo", groupvars = c("shRNA1","Timepoint"))
df.FracFluo.stim2<-summarySE(stim.field.spot, measurevar = "meanFracFluo", groupvars = c("shRNA2"))
df.FracFluo.stim3<-summarySE(stim.field.spot, measurevar = "meanFracFluo", groupvars = c("shRNA1", "Timepoint"))

# graphs fraction of active pixels
ggplot(data=df.FracActive, aes(x=shRNA1, y= meanFracActive, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanFracActive-se, ymax=meanFracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.FracActive2, aes(x=shRNA2, y= meanFracActive, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanFracActive-se, ymax=meanFracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.FracActive.TP2, aes(x=interaction(Timepoint, shRNA2), y= meanFracActive, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanFracActive-se, ymax=meanFracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.FracActive.TP, aes(x=interaction(Timepoint, shRNA1), y= meanFracActive, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanFracActive-se, ymax=meanFracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# graphs response probabliity (potential for same region to appear in mutliple trials)
ggplot(data=df.ResponseScore, aes(x=shRNA1, y= Response_Score, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Response_Score-se, ymax=Response_Score+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ResponseScore") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.ResponseScore2, aes(x=shRNA2, y= Response_Score, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Response_Score-se, ymax=Response_Score+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ResponseScore") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.ResponseScore.TP2, aes(x=interaction(Timepoint,shRNA2), y= Response_Score, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Response_Score-se, ymax=Response_Score+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ResponseScore") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.ResponseScore.TP, aes(x=interaction(Timepoint,shRNA1), y= Response_Score, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Response_Score-se, ymax=Response_Score+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("ResponseScore") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# average fluorescent pixels
ggplot(data=df.nFluoPix.nostim1, aes(x=shRNA1, y= meanFluoPix, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanFluoPix-se, ymax=meanFluoPix+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("# of fluorescent pixels") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.FracActive2, aes(x=shRNA2, y= meanFracActive, fill=Cond)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=meanFracActive-se, ymax=meanFracActive+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("FracActive") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme
#################################
# PEAK DATA

#####

# RCaMP peaks (nostim and stim) near stimulation time
all.peaks.RC<- subset(all.peaks.window, Channel=="RCaMP")

#spontaneous RCaMP peaks from whole trial
nostim.peaks.RC<-subset(all.peaks, Condition=="nostim" & Channel=="RCaMP")
nostim.peaks.RC<-subset(nostim.peaks.RC, amplitude>0 & Duration<20)

# stim RCaMP peaks only from time near stimulation
stim.peaks.RC<-subset(all.peaks.RC, Condition=="stim")

# GCaMP peaks (nostim and stim) near stimulation time
all.peaks.GC<- subset(all.peaks.window, Channel=="GCaMP")

#spontaneous RCaMP peaks from whole trial
nostim.peaks.GC<-subset(all.peaks, Condition=="nostim" & Channel=="GCaMP")
nostim.peaks.GC<-subset(nostim.peaks.GC, amplitude>0 & Duration<20)

# stim RCaMP peaks only from time near stimulation
stim.peaks.GC<-subset(all.peaks.GC, Condition=="stim")


######
# amplitude histograms
ggplot(stim.peaks.RC,aes(x=amplitude,y=..density..,fill=shRNA2)) +
  geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("RCaMP stim amplitudes")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(nostim.peaks.RC,aes(x=amplitude,y=..density..,fill=shRNA2)) +
  geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("RCaMP nostim amplitudes")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.GC,aes(x=amplitude,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.2,position="dodge") +
  ggtitle("GCaMP stim amplitudes")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(nostim.peaks.GC,aes(x=amplitude,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.2,position="dodge") +
  ggtitle("GCaMP nostim amplitudes")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(nostim.peaks.GC,aes(x=amplitude,y=..density..,fill=shRNA2)) +
  geom_histogram(binwidth=0.2,position="dodge") +
  ggtitle("GCaMP nostim amplitudes")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme


#duration histograms
ggplot(all.peaks.RC[(all.peaks.RC$Condition=="stim"),],aes(x=Duration,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=1, position="dodge") +
  ggtitle("RCaMP signal duration")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(all.peaks.GC[(all.peaks.GC$Condition=="stim"),],aes(x=Duration,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=1,position="dodge") +
  ggtitle("GCaMP signal duration")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

#############

#mean amplitude of signals for each ROI

#RCaMP
df.amp1.RC<- summarySE(all.peaks.RC, measurevar = "amplitude", groupvars = c("Condition","shRNA1"))
df.amp2.RC<- summarySE(all.peaks.RC, measurevar = "amplitude", groupvars = c("Condition","shRNA2"))
df.amp.TP.RC<-summarySE(all.peaks.RC, measurevar = "amplitude", groupvars = c("Condition","shRNA1", "Timepoint" ))

# RCaMP spontaneous only
df.amp1.RC.nostim<- summarySE(nostim.peaks.RC, measurevar = "amplitude", groupvars = c("shRNA1"))
df.amp2.RC.nostim<- summarySE(nostim.peaks.RC, measurevar = "amplitude", groupvars = c("shRNA2"))
df.amp.TP.RC.nostim<-summarySE(nostim.peaks.RC, measurevar = "amplitude", groupvars = c("shRNA1", "Timepoint" ))

# RCaMP stim only
df.amp1.RC.stim<- summarySE(stim.peaks.RC, measurevar = "amplitude", groupvars = c("shRNA1"))
df.amp2.RC.stim<- summarySE(stim.peaks.RC, measurevar = "amplitude", groupvars = c("shRNA2"))
df.amp.TP.RC.stim<-summarySE(stim.peaks.RC, measurevar = "amplitude", groupvars = c("shRNA1", "Timepoint" ))


# RCaMP graphs
ggplot(data=df.amp1.RC, aes(x=shRNA1, y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp2.RC, aes(x=shRNA2, y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# no stim
ggplot(data=df.amp1.RC.nostim, aes(x=shRNA1, y= amplitude, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("no stim amplitude") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(nostim.peaks.RC, aes(x=shRNA1,y=amplitude, fill= shRNA1)) +
  geom_boxplot()+
  ylab("no stim RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp2.RC.nostim, aes(x=shRNA2, y= amplitude, fill=shRNA2)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("no stim RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(nostim.peaks.RC, aes(x=shRNA2,y=amplitude, fill= shRNA2)) +
  geom_boxplot()+
  ylab("amplitude") +
  ggtitle("RCaMP no stim") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp.TP.RC.nostim, aes(x=Timepoint, y= amplitude, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("no stim amplitude") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# stim
ggplot(data=df.amp1.RC.stim, aes(x=shRNA1, y= amplitude, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("stim amplitude") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.RC, aes(x=shRNA1,y=amplitude, fill= shRNA1)) +
  geom_boxplot()+
  ggtitle("stim RCaMP") +
  ylab("amplitude") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp2.RC.stim, aes(x=shRNA2, y= amplitude, fill=shRNA2)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("stim RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.RC, aes(x=shRNA2,y=amplitude, fill= shRNA2)) +
  geom_boxplot()+
  ylab("amplitude") +
  ggtitle("stim RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp.TP.RC.stim, aes(x=Timepoint, y= amplitude, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("stim amplitude") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


# GCaMP amplitude
df.amp1.GC<- summarySE(all.peaks.GC, measurevar = "amplitude", groupvars = c("Condition","shRNA1"))
df.amp2.GC<- summarySE(all.peaks.GC, measurevar = "amplitude", groupvars = c("Condition","shRNA2"))
df.amp.TP.GC<-summarySE(all.peaks.GC, measurevar = "amplitude", groupvars = c("Condition","shRNA1", "Timepoint" ))

# GCaMP spontaneous only
df.amp1.GC.nostim<- summarySE(nostim.peaks.GC, measurevar = "amplitude", groupvars = c("shRNA1"))
df.amp2.GC.nostim<- summarySE(nostim.peaks.GC, measurevar = "amplitude", groupvars = c("shRNA2"))
df.amp.TP.GC.nostim<-summarySE(nostim.peaks.GC, measurevar = "amplitude", groupvars = c("shRNA1", "Timepoint" ))

# GCaMP stim only
df.amp1.GC.stim<- summarySE(stim.peaks.GC, measurevar = "amplitude", groupvars = c("shRNA1"))
df.amp2.GC.stim<- summarySE(stim.peaks.GC, measurevar = "amplitude", groupvars = c("shRNA2"))
df.amp.TP.GC.stim<-summarySE(stim.peaks.GC, measurevar = "amplitude", groupvars = c("shRNA1", "Timepoint" ))


ggplot(data=df.amp1.GC, aes(x=shRNA1, y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


ggplot(data=df.amp2.GC, aes(x=shRNA2, y= amplitude, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

#no stim
ggplot(data=df.amp2.GC.nostim, aes(x=shRNA2, y= amplitude, fill=shRNA2)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("no stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(nostim.peaks.GC, aes(x=shRNA2,y=amplitude, fill= shRNA2)) +
  geom_boxplot()+
  ylab("amplitude") +
  ggtitle("no stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(nostim.peaks.GC, aes(x=shRNA1,y=amplitude, fill= shRNA1)) +
  geom_boxplot()+
  ylab("amplitude") +
  ggtitle("no stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp1.GC.nostim, aes(x=shRNA1, y= amplitude, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("no stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp.TP.GC.nostim, aes(x=Timepoint, y= amplitude, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("no stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# stim
ggplot(data=df.amp2.GC.stim, aes(x=shRNA2, y= amplitude, fill=shRNA2)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp1.GC.stim, aes(x=shRNA1, y= amplitude, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.GC, aes(x=shRNA2,y=amplitude, fill= shRNA2)) +
  geom_boxplot()+
  ylab("amplitude") +
  ggtitle("stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.GC, aes(x=shRNA1,y=amplitude, fill= shRNA1)) +
  geom_boxplot()+
  ylab("stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.amp.TP.GC.stim, aes(x=Timepoint, y= amplitude, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=amplitude-se, ymax=amplitude+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("amplitude") +
  ggtitle("stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

#####################
## STATS

peaks.RC.nostim.shRNA1_TP= interaction(nostim.peaks.RC$shRNA1, nostim.peaks.RC$Timepoint)
peaks.RC.stim.shRNA1_TP= interaction(stim.peaks.RC$shRNA1, stim.peaks.RC$Timepoint)

peaks.GC.nostim.shRNA1_TP= interaction(nostim.peaks.GC$shRNA1, nostim.peaks.GC$Timepoint)
peaks.GC.stim.shRNA1_TP= interaction(stim.peaks.GC$shRNA1, stim.peaks.GC$Timepoint)



#RCaMP nostim
amp.RC.NS.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.RC,REML=FALSE)
amp.RC.NS.model1 = lmer(amplitude ~ shRNA1 + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.RC,REML=FALSE)
amp.RC.NS.model2 = lmer(amplitude ~ shRNA2 + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.RC,REML=FALSE)
amp.RC.NS.model3 = lmer(amplitude ~ peaks.RC.nostim.shRNA1_TP + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.RC,REML=FALSE)
amp.RC.NS.anova <- anova(amp.RC.NS.null, amp.RC.NS.model1, amp.RC.NS.model2, amp.RC.NS.model3)
print(amp.RC.NS.anova)

amp.RC.NS.shRNA1<- glht(amp.RC.NS.model1, mcp(shRNA1= "Tukey"))
summary(amp.RC.NS.shRNA1)

amp.RC.NS.shRNA2<- glht(amp.RC.NS.model2, mcp(shRNA2= "Tukey"))
summary(amp.RC.NS.shRNA2)

#RCaMP stim
amp.RC.S.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
amp.RC.S.model1 = lmer(amplitude ~ shRNA1 + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
amp.RC.S.model2 = lmer(amplitude ~ shRNA2 + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
amp.RC.S.model3 = lmer(amplitude ~ peaks.RC.stim.shRNA1_TP + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
amp.RC.S.anova <- anova(amp.RC.S.null, amp.RC.S.model1, amp.RC.S.model2, amp.RC.S.model3)
print(amp.RC.S.anova)

amp.RC.S.shRNA1<- glht(amp.RC.S.model1, mcp(shRNA1= "Tukey"))
summary(amp.RC.S.shRNA1)

amp.RC.S.shRNA2<- glht(amp.RC.S.model2, mcp(shRNA2= "Tukey"))
summary(amp.RC.S.shRNA2)


#GCaMP nostim
amp.GC.NS.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.GC,REML=FALSE)
amp.GC.NS.model1 = lmer(amplitude ~ shRNA1 + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.GC,REML=FALSE)
amp.GC.NS.model2 = lmer(amplitude ~ shRNA2 + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.GC,REML=FALSE)
amp.GC.NS.model3 = lmer(amplitude ~ peaks.GC.nostim.shRNA1_TP + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.GC,REML=FALSE)
amp.GC.NS.anova <- anova(amp.GC.NS.null, amp.GC.NS.model1, amp.GC.NS.model2, amp.GC.NS.model3)
print(amp.GC.NS.anova)

amp.GC.NS.shRNA1<- glht(amp.GC.NS.model1, mcp(shRNA1= "Tukey"))
summary(amp.GC.NS.shRNA1)

amp.GC.NS.shRNA2<- glht(amp.GC.NS.model2, mcp(shRNA2= "Tukey"))
summary(amp.GC.NS.shRNA2)

#GCaMP stim
amp.GC.S.null = lmer(amplitude ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
amp.GC.S.model1 = lmer(amplitude ~ shRNA1 + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
amp.GC.S.model2 = lmer(amplitude ~ shRNA2 + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
amp.GC.S.model3 = lmer(amplitude ~ peaks.GC.stim.shRNA1_TP + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
amp.GC.S.anova <- anova(amp.GC.S.null, amp.GC.S.model1, amp.GC.S.model2, amp.GC.S.model3)
print(amp.GC.S.anova)

amp.GC.S.shRNA1<- glht(amp.GC.S.model1, mcp(shRNA1= "Tukey"))
summary(amp.GC.S.shRNA1)

amp.GC.S.shRNA2<- glht(amp.GC.S.model2, mcp(shRNA2= "Tukey"))
summary(amp.GC.S.shRNA2)

amp.GC.S.shRNA1.TP<- glht(amp.GC.S.model3, mcp(peaks.GC.stim.shRNA1_TP = "Tukey"))
summary(amp.GC.S.shRNA1.TP)

#########
#mean Duration of signals for each ROI

#RCaMP
df.dur1.RC<- summarySE(all.peaks.RC, measurevar = "Duration", groupvars = c("Condition","shRNA1"))
df.dur2.RC<- summarySE(all.peaks.RC, measurevar = "Duration", groupvars = c("Condition","shRNA2"))
df.dur.TP.RC<-summarySE(all.peaks.RC, measurevar = "Duration", groupvars = c("Condition","shRNA1", "Timepoint" ))

# RCdur spontaneous only
df.dur1.RC.nostim<- summarySE(nostim.peaks.RC, measurevar = "Duration", groupvars = c("shRNA1"))
df.dur2.RC.nostim<- summarySE(nostim.peaks.RC, measurevar = "Duration", groupvars = c("shRNA2"))
df.dur.TP.RC.nostim<-summarySE(nostim.peaks.RC, measurevar = "Duration", groupvars = c("shRNA1", "Timepoint" ))

# RCaMP stim only
df.dur1.RC.stim<- summarySE(stim.peaks.RC, measurevar = "Duration", groupvars = c("shRNA1"))
df.dur2.RC.stim<- summarySE(stim.peaks.RC, measurevar = "Duration", groupvars = c("shRNA2"))
df.dur.TP.RC.stim<-summarySE(stim.peaks.RC, measurevar = "Duration", groupvars = c("shRNA1", "Timepoint" ))


# RCaMP graphs
ggplot(data=df.dur1.RC, aes(x=shRNA1, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.dur2.RC, aes(x=shRNA2, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# no stim
ggplot(data=df.dur1.RC.nostim, aes(x=shRNA1, y= Duration, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("no stim Duration") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(nostim.peaks.RC, aes(x=shRNA1,y=Duration, fill= shRNA1)) +
  geom_boxplot()+
  ylab("no stim RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.dur2.RC.nostim, aes(x=shRNA2, y= Duration, fill=shRNA2)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("no stim RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(nostim.peaks.RC, aes(x=shRNA2,y=Duration, fill= shRNA2)) +
  geom_boxplot()+
  ylab("Duration") +
  ggtitle("no stim RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.dur.TP.RC.nostim, aes(x=Timepoint, y= Duration, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("no stim Duration") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# stim
ggplot(data=df.dur1.RC.stim, aes(x=shRNA1, y= Duration, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("stim Duration") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.RC, aes(x=shRNA1,y=Duration, fill= shRNA1)) +
  geom_boxplot()+
  ylab("stim RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.dur2.RC.stim, aes(x=shRNA2, y= Duration, fill=shRNA2)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("stim RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.RC, aes(x=shRNA2,y=Duration, fill= shRNA2)) +
  geom_boxplot()+
  ylab("Duration") +
  ggtitle("stim RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.dur.TP.RC.stim, aes(x=Timepoint, y= Duration, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("stim Duration") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


# GCaMP Duration
df.dur1.GC<- summarySE(all.peaks.GC, measurevar = "Duration", groupvars = c("Condition","shRNA1"))
df.dur2.GC<- summarySE(all.peaks.GC, measurevar = "Duration", groupvars = c("Condition","shRNA2"))
df.dur.TP.GC<-summarySE(all.peaks.GC, measurevar = "Duration", groupvars = c("Condition","shRNA1", "Timepoint" ))

# GCaMP spontaneous only
df.dur1.GC.nostim<- summarySE(nostim.peaks.GC, measurevar = "Duration", groupvars = c("shRNA1"))
df.dur2.GC.nostim<- summarySE(nostim.peaks.GC, measurevar = "Duration", groupvars = c("shRNA2"))
df.dur.TP.GC.nostim<-summarySE(nostim.peaks.GC, measurevar = "Duration", groupvars = c("shRNA1", "Timepoint" ))

# GCaMP stim only
df.dur1.GC.stim<- summarySE(stim.peaks.GC, measurevar = "Duration", groupvars = c("shRNA1"))
df.dur2.GC.stim<- summarySE(stim.peaks.GC, measurevar = "Duration", groupvars = c("shRNA2"))
df.dur.TP.GC.stim<-summarySE(stim.peaks.GC, measurevar = "Duration", groupvars = c("shRNA1", "Timepoint" ))


ggplot(data=df.dur1.GC, aes(x=shRNA1, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme


ggplot(data=df.dur2.GC, aes(x=shRNA2, y= Duration, fill=Condition)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

#no stim
ggplot(data=df.dur2.GC.nostim, aes(x=shRNA2, y= Duration, fill=shRNA2)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("no stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(nostim.peaks.GC, aes(x=shRNA2,y=Duration, fill= shRNA2)) +
  geom_boxplot()+
  ggtitle("no stim GCaMP") +
  ylab("Duration") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(nostim.peaks.GC, aes(x=shRNA1,y=Duration, fill= shRNA1)) +
  geom_boxplot()+
  ylab("no stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.dur1.GC.nostim, aes(x=shRNA1, y= Duration, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("no stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.dur.TP.GC.nostim, aes(x=Timepoint, y= Duration, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("no stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# stim
ggplot(data=df.dur2.GC.stim, aes(x=shRNA2, y= Duration, fill=shRNA2)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.dur1.GC.stim, aes(x=shRNA1, y= Duration, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.GC, aes(x=shRNA2,y=Duration, fill= shRNA2)) +
  geom_boxplot()+
  ylab("Duration") +
  ggtitle("stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.GC, aes(x=shRNA1,y=Duration, fill= shRNA1)) +
  geom_boxplot()+
  ylab("stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(data=df.dur.TP.GC.stim, aes(x=Timepoint, y= Duration, fill=shRNA1)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Duration-se, ymax=Duration+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Duration") +
  ggtitle("stim GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

#####################
## STATS

#RCaMP nostim
dur.RC.NS.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.RC,REML=FALSE)
dur.RC.NS.model1 = lmer(Duration ~ shRNA1 + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.RC,REML=FALSE)
dur.RC.NS.model2 = lmer(Duration ~ shRNA2 + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.RC,REML=FALSE)
dur.RC.NS.model3 = lmer(Duration ~ peaks.RC.nostim.shRNA1_TP + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.RC,REML=FALSE)
dur.RC.NS.anova <- anova(dur.RC.NS.null, dur.RC.NS.model1, dur.RC.NS.model2, dur.RC.NS.model3)
print(dur.RC.NS.anova)

dur.RC.NS.shRNA1<- glht(dur.RC.NS.model1, mcp(shRNA1= "Tukey"))
summary(dur.RC.NS.shRNA1)

dur.RC.NS.shRNA2<- glht(dur.RC.NS.model2, mcp(shRNA2= "Tukey"))
summary(dur.RC.NS.shRNA2)

#RCaMP stim
dur.RC.S.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
dur.RC.S.model1 = lmer(Duration ~ shRNA1 + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
dur.RC.S.model2 = lmer(Duration ~ shRNA2 + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
dur.RC.S.model3 = lmer(Duration ~ peaks.RC.stim.shRNA1_TP + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
dur.RC.S.anova <- anova(dur.RC.S.null, dur.RC.S.model1, dur.RC.S.model2, dur.RC.S.model3)
print(dur.RC.S.anova)

dur.RC.S.shRNA1<- glht(dur.RC.S.model1, mcp(shRNA1= "Tukey"))
summary(dur.RC.S.shRNA1)

dur.RC.S.shRNA2<- glht(dur.RC.S.model2, mcp(shRNA2= "Tukey"))
summary(dur.RC.S.shRNA2)


#GCaMP nostim
dur.GC.NS.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.GC,REML=FALSE)
dur.GC.NS.model1 = lmer(Duration ~ shRNA1 + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.GC,REML=FALSE)
dur.GC.NS.model2 = lmer(Duration ~ shRNA2 + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.GC,REML=FALSE)
dur.GC.NS.model3 = lmer(Duration ~ peaks.GC.nostim.shRNA1_TP + (1|Animal) + (1|Spot) + (1|ROIs_trial), nostim.peaks.GC,REML=FALSE)
dur.GC.NS.anova <- anova(dur.GC.NS.null, dur.GC.NS.model1, dur.GC.NS.model2, dur.GC.NS.model3)
print(dur.GC.NS.anova)

dur.GC.NS.shRNA1<- glht(dur.GC.NS.model1, mcp(shRNA1= "Tukey"))
summary(dur.GC.NS.shRNA1)

dur.GC.NS.shRNA2<- glht(dur.GC.NS.model2, mcp(shRNA2= "Tukey"))
summary(dur.GC.NS.shRNA2)

#GCaMP stim
dur.GC.S.null = lmer(Duration ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
dur.GC.S.model1 = lmer(Duration ~ shRNA1 + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
dur.GC.S.model2 = lmer(Duration ~ shRNA2 + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
dur.GC.S.model3 = lmer(Duration ~ peaks.GC.stim.shRNA1_TP + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
dur.GC.S.anova <- anova(dur.GC.S.null, dur.GC.S.model1, dur.GC.S.model2, dur.GC.S.model3)
print(dur.GC.S.anova)

dur.GC.S.shRNA1<- glht(dur.GC.S.model1, mcp(shRNA1= "Tukey"))
summary(dur.GC.S.shRNA1)

dur.GC.S.shRNA2<- glht(dur.GC.S.model2, mcp(shRNA2= "Tukey"))
summary(dur.GC.S.shRNA2)

dur.GC.S.shRNA1.TP<- glht(dur.GC.S.model3, mcp(peaks.GC.stim.shRNA1_TP = "Tukey"))
summary(dur.GC.S.shRNA1.TP)




########################
# ONSET TIME
#remove entries with no onset time

all.OT<-subset(all.OT, OnsetTime!="NaN")

stim.OT.GC<-subset(all.OT, Channel=="GCaMP" & Condition =="stim")
stim.OT.RC<-subset(all.OT, Channel=="RCaMP" & Condition =="stim")

######
# neuronal responses to stimulation
NeuronalStim<-subset(all.OT, Channel=="RCaMP" & Condition=="stim" & OnsetTime<8)

# should have an onset time in 8 s stimulus
ggplot(NeuronalStim,aes(x=OnsetTime,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP onset times between 0 and 8 s from stim trials")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# should have an onset time in 8 s stimulus
ggplot(NeuronalStim,aes(x=OnsetTime,y=..density..,fill=shRNA2)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP onset times between 0 and 8 s from stim trials")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme


# all responding neurons
Neuron95Onset<-quantile(NeuronalStim$OnsetTime, prob = seq(0, 1, length = 21), type = 5, na.rm=TRUE)
NeuronPT50<-Neuron95Onset[[11]]
print(NeuronPT50)

# time thresold to consider an astrocyte to be fast:
fastTh<-NeuronPT50

#######
#plot more distributions

# GCaMP onset times
ggplot(stim.OT.GC,aes(x=OnsetTime,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP onset times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.OT.GC,aes(x=OnsetTime,y=..density..,fill=shRNA2)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP onset times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme



ggplot(stim.peaks.GC,aes(x=peakTime,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("GCaMP peak times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(all.peaks[(all.peaks$Channel=="RCaMP"& all.peaks$Condition=="stim"),],aes(x=peakTime,y=..density..,fill=shRNA1)) +
  geom_histogram(binwidth=0.084, position="dodge") +
  ggtitle("RCaMP peak times")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

#####
# onset time boxplots

stim.OT.GC.window<-subset(stim.OT.GC, OnsetTime<8)
stim.OT.RC.window<-subset(stim.OT.RC, OnsetTime<8)

ggplot(stim.OT.GC.window, aes(x=shRNA1,y=OnsetTime, fill= shRNA1)) +
  geom_boxplot()+
  ylab("Onset Latency (s)") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.OT.GC.window, aes(x=shRNA2,y=OnsetTime, fill= shRNA2)) +
  geom_boxplot()+
  ylab("Onset Latency (s)") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.OT.RC.window, aes(x=shRNA1,y=OnsetTime, fill= shRNA1)) +
  geom_boxplot()+
  ylab("Onset Latency (s)") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.OT.RC.window, aes(x=shRNA2,y=OnsetTime, fill= shRNA2)) +
  geom_boxplot()+
  ylab("Onset Latency (s)") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# peak time
ggplot(stim.peaks.GC, aes(x=shRNA1,y=peakTime, fill= shRNA1)) +
  geom_boxplot()+
  ylab("Latency to Peak Max (s)") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.GC, aes(x=shRNA2,y=peakTime, fill= shRNA2)) +
  geom_boxplot()+
  ylab("Latency to Peak Max (s)") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.RC, aes(x=shRNA1,y=peakTime, fill= shRNA1)) +
  geom_boxplot()+
  ylab("Latency to Peak Max (s)") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(stim.peaks.RC, aes(x=shRNA2,y=peakTime, fill= shRNA2)) +
  geom_boxplot()+
  ylab("Latency to Peak Max (s)") +
  ggtitle("RCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

######
## STATS

#RCaMP onset time
OT.RC.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot), stim.OT.RC.window,REML=FALSE)
OT.RC.model1 = lmer(OnsetTime ~ shRNA1 + (1|Animal) + (1|Spot), stim.OT.RC.window,REML=FALSE)
OT.RC.model2 = lmer(OnsetTime ~ shRNA2 + (1|Animal) + (1|Spot), stim.OT.RC.window,REML=FALSE)
OT.RC.anova <- anova(OT.RC.null, OT.RC.model1, OT.RC.model2)
print(OT.RC.anova)

OT.RC.shRNA1<- glht(OT.RC.model1, mcp(shRNA1= "Tukey"))
summary(OT.RC.shRNA1)

OT.RC.shRNA2<- glht(OT.RC.model2, mcp(shRNA2= "Tukey"))
summary(OT.RC.shRNA2)

#GCaMP onset time
OT.GC.null = lmer(OnsetTime ~ (1|Animal) + (1|Spot), stim.OT.GC.window,REML=FALSE)
OT.GC.model1 = lmer(OnsetTime ~ shRNA1 + (1|Animal) + (1|Spot), stim.OT.GC.window,REML=FALSE)
OT.GC.model2 = lmer(OnsetTime ~ shRNA2 + (1|Animal) + (1|Spot), stim.OT.GC.window,REML=FALSE)
OT.GC.anova <- anova(OT.GC.null, OT.GC.model1, OT.GC.model2)
print(OT.GC.anova)

OT.GC.shRNA1<- glht(OT.GC.model1, mcp(shRNA1= "Tukey"))
summary(OT.GC.shRNA1)

OT.GC.shRNA2<- glht(OT.GC.model2, mcp(shRNA2= "Tukey"))
summary(OT.GC.shRNA2)


#RCaMP stim peak time
pT.RC.S.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
pT.RC.S.model1 = lmer(peakTime ~ shRNA1 + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
pT.RC.S.model2 = lmer(peakTime~ shRNA2 + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.RC,REML=FALSE)
pT.RC.S.anova <- anova(pT.RC.S.null, pT.RC.S.model1, pT.RC.S.model2)
print(pT.RC.S.anova)

pT.RC.S.shRNA1<- glht(pT.RC.S.model1, mcp(shRNA1= "Tukey"))
summary(pT.RC.S.shRNA1)

pT.RC.S.shRNA2<- glht(pT.RC.S.model2, mcp(shRNA2= "Tukey"))
summary(pT.RC.S.shRNA2)


#GCaMP stim peak time
pT.GC.S.null = lmer(peakTime ~ (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
pT.GC.S.model1 = lmer(peakTime ~ shRNA1 + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
pT.GC.S.model2 = lmer(peakTime~ shRNA2 + (1|Animal) + (1|Spot) + (1|ROIs_trial), stim.peaks.GC,REML=FALSE)
pT.GC.S.anova <- anova(pT.GC.S.null, pT.GC.S.model1, pT.GC.S.model2)
print(pT.GC.S.anova)

pT.GC.S.shRNA1<- glht(pT.GC.S.model1, mcp(shRNA1= "Tukey"))
summary(pT.GC.S.shRNA1)

pT.GC.S.shRNA2<- glht(pT.GC.S.model2, mcp(shRNA2= "Tukey"))
summary(pT.GC.S.shRNA2)


#########
# find fast and delayed astrocyte microdomains


# identify "FAST" astrocytes
stim.OT.GC$Group<-0
stim.OT.GC$Group[stim.OT.GC$OnsetTime<fastTh]<-"fast"
stim.OT.GC$Group[stim.OT.GC$OnsetTime>=fastTh]<-"delayed"


# add onset time information to the peak data table
stim.peaks.GC<-merge(stim.peaks.GC, stim.OT.GC[, c("ROIs_trial_Cond", "OnsetTime", "Group")], by="ROIs_trial_Cond", all.x=TRUE)


# remove peaks that don't have a corresponding onset time

stim.peaks.GC<-stim.peaks.GC[!is.na(stim.peaks.GC$OnsetTime),]

###################

# number of active ROIs per field of view


stim.peaks.GC$Animal_Spot<- paste(stim.peaks.GC$Animal, stim.peaks.GC$Spot, sep="_")
all.peaks.window$Animal_Spot<- paste(all.peaks.window$Animal, all.peaks.window$Spot, sep="_")

# number of ROIs in each trial for each field of view (across the whole trial) with a peak during no stim and stim
# only consider during the stimulus (and the same time window in no stim trials)

ROInum.8strial<-ddply(all.peaks.window, c("Animal","Spot","shRNA1","shRNA2", "Condition","Channel","Animal_Spot","Timepoint"), summarise, nROIs=length(unique(ROIs_trial_Cond)))

# only GCaMP stim peaks with fast or delayed
ROInum.8strial.group<-ddply(stim.peaks.GC, c("Animal","Spot","shRNA1","shRNA2","Condition","Channel","Group","Timepoint"), summarise, nROIs=length(unique(ROIs_trial_Cond)))


# add in number of trials
ROInum.8strial$Ani_Spot_Cond<-paste(ROInum.8strial$Animal_Spot, ROInum.8strial$Condition, sep="_")
ROInum.8strial.group$Ani_Spot_Cond<-paste(ROInum.8strial.group$Animal, ROInum.8strial.group$Spot, ROInum.8strial.group$Condition, sep="_")

ROInum.8strial<-merge(ROInum.8strial, Spot.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)
ROInum.8strial$ROIsPerTrial<-ROInum.8strial$nROIs/ROInum.8strial$nTrials

ROInum.8strial.group<-merge(ROInum.8strial.group, Spot.ntrials[, c("Ani_Spot_Cond", "nTrials")], by="Ani_Spot_Cond", all.x=TRUE)
ROInum.8strial.group$ROIsPerTrial<-ROInum.8strial.group$nROIs/ROInum.8strial.group$nTrials

ROInum.8strial.stim<-subset(ROInum.8strial, Condition=="stim")
ROInum.8strial.nostim<-subset(ROInum.8strial, Condition=="nostim")

# means
df.ROInum.8strial.NS.RC1<-summarySE(ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="RCaMP",], measurevar = "ROIsPerTrial", groupvars = c("shRNA1"))
df.ROInum.8strial.NS.RC2<-summarySE(ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="RCaMP",], measurevar = "ROIsPerTrial", groupvars = c("shRNA2"))
df.ROInum.8strial.RC.NS.TP<-summarySE(ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="RCaMP",], measurevar = "ROIsPerTrial", groupvars = c("shRNA1","Timepoint"))
df.ROInum.8strial.S.RC1<-summarySE(ROInum.8strial.stim[ROInum.8strial.stim$Channel=="RCaMP",], measurevar = "ROIsPerTrial", groupvars = c("shRNA1"))
df.ROInum.8strial.S.RC2<-summarySE(ROInum.8strial.stim[ROInum.8strial.stim$Channel=="RCaMP",], measurevar = "ROIsPerTrial", groupvars = c("shRNA2"))
df.ROInum.8strial.RC.S.TP<-summarySE(ROInum.8strial.stim[ROInum.8strial.stim$Channel=="RCaMP",], measurevar = "ROIsPerTrial", groupvars = c("shRNA1","Timepoint"))

df.ROInum.8strial.GC1<-summarySE(ROInum.8strial[ROInum.8strial$Channel=="GCaMP",], measurevar = "ROIsPerTrial", groupvars = c("Condition","shRNA1"))
df.ROInum.8strial.GC2<-summarySE(ROInum.8strial[ROInum.8strial$Channel=="GCaMP",], measurevar = "ROIsPerTrial", groupvars = c("Condition","shRNA2"))
df.ROInum.8strial.GC.TP<-summarySE(ROInum.8strial[ROInum.8strial$Channel=="GCaMP",], measurevar = "ROIsPerTrial", groupvars = c("Condition","shRNA1","Timepoint"))

# plots
# RCaMP

ggplot(ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="RCaMP",], aes(x=shRNA1,y=ROIsPerTrial, fill= shRNA1)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("RCaMP nostim") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="RCaMP",], aes(x=shRNA2,y=ROIsPerTrial, fill= shRNA2)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("RCaMP nostim") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(ROInum.8strial.stim[ROInum.8strial.stim$Channel=="RCaMP",], aes(x=shRNA1,y=ROIsPerTrial, fill= shRNA1)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("RCaMP stim") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(ROInum.8strial.stim[ROInum.8strial.stim$Channel=="RCaMP",], aes(x=shRNA2,y=ROIsPerTrial, fill= shRNA2)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("RCaMP stim") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



# GCaMP
ggplot(ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="GCaMP",], aes(x=shRNA1,y=ROIsPerTrial, fill= shRNA1)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("GCaMP nostim") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="GCaMP",], aes(x=shRNA2,y=ROIsPerTrial, fill= shRNA2)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("GCaMP nostim") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(ROInum.8strial.stim[ROInum.8strial.stim$Channel=="GCaMP",], aes(x=shRNA1,y=ROIsPerTrial, fill= shRNA1)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("GCaMP stim") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(ROInum.8strial.stim[ROInum.8strial.stim$Channel=="GCaMP",], aes(x=shRNA2,y=ROIsPerTrial, fill= shRNA2)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("GCaMP stim") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



###############
# stats for active ROI number per trials per FOV

# RCaMP nostim
ROInum.RC.NS.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="RCaMP",],REML=FALSE)
ROInum.RC.NS.model1 = lmer(ROIsPerTrial ~ shRNA1 + (1|Animal), ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="RCaMP",],REML=FALSE)
ROInum.RC.NS.model2 = lmer(ROIsPerTrial ~ shRNA2 + (1|Animal), ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="RCaMP",],REML=FALSE)
ROInum.RC.NS.anova <- anova(ROInum.RC.NS.null, ROInum.RC.NS.model1, ROInum.RC.NS.model2)
print(ROInum.RC.NS.anova)

ROInum.RC.NS.shRNA1<- glht(ROInum.RC.NS.model1, mcp(shRNA1= "Tukey"))
summary(ROInum.RC.NS.shRNA1)

ROInum.RC.NS.shRNA2<- glht(ROInum.RC.NS.model2, mcp(shRNA2= "Tukey"))
summary(ROInum.RC.NS.shRNA2)

# RCaMP stim
ROInum.RC.S.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.8strial.nostim[ROInum.8strial.stim$Channel=="RCaMP",],REML=FALSE)
ROInum.RC.S.model1 = lmer(ROIsPerTrial ~ shRNA1 + (1|Animal), ROInum.8strial.nostim[ROInum.8strial.stim$Channel=="RCaMP",],REML=FALSE)
ROInum.RC.S.model2 = lmer(ROIsPerTrial ~ shRNA2 + (1|Animal), ROInum.8strial.nostim[ROInum.8strial.stim$Channel=="RCaMP",],REML=FALSE)
ROInum.RC.S.anova <- anova(ROInum.RC.S.null, ROInum.RC.S.model1, ROInum.RC.S.model2)
print(ROInum.RC.S.anova)

ROInum.RC.S.shRNA1<- glht(ROInum.RC.S.model1, mcp(shRNA1= "Tukey"))
summary(ROInum.RC.S.shRNA1)

ROInum.RC.S.shRNA2<- glht(ROInum.RC.S.model2, mcp(shRNA2= "Tukey"))
summary(ROInum.RC.S.shRNA2)

# GCaMP nostim
ROInum.GC.NS.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="GCaMP",],REML=FALSE)
ROInum.GC.NS.model1 = lmer(ROIsPerTrial ~ shRNA1 + (1|Animal), ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="GCaMP",],REML=FALSE)
ROInum.GC.NS.model2 = lmer(ROIsPerTrial ~ shRNA2 + (1|Animal), ROInum.8strial.nostim[ROInum.8strial.nostim$Channel=="GCaMP",],REML=FALSE)
ROInum.GC.NS.anova <- anova(ROInum.GC.NS.null, ROInum.GC.NS.model1, ROInum.GC.NS.model2)
print(ROInum.GC.NS.anova)

ROInum.GC.NS.shRNA1<- glht(ROInum.GC.NS.model1, mcp(shRNA1= "Tukey"))
summary(ROInum.GC.NS.shRNA1)

ROInum.GC.NS.shRNA2<- glht(ROInum.GC.NS.model2, mcp(shRNA2= "Tukey"))
summary(ROInum.GC.NS.shRNA2)

# GCaMP stim
ROInum.GC.S.null = lmer(ROIsPerTrial ~ (1|Animal), ROInum.8strial.nostim[ROInum.8strial.stim$Channel=="GCaMP",],REML=FALSE)
ROInum.GC.S.model1 = lmer(ROIsPerTrial ~ shRNA1 + (1|Animal), ROInum.8strial.nostim[ROInum.8strial.stim$Channel=="GCaMP",],REML=FALSE)
ROInum.GC.S.model2 = lmer(ROIsPerTrial ~ shRNA2 + (1|Animal), ROInum.8strial.nostim[ROInum.8strial.stim$Channel=="GCaMP",],REML=FALSE)
ROInum.GC.S.anova <- anova(ROInum.GC.S.null, ROInum.GC.S.model1, ROInum.GC.S.model2)
print(ROInum.GC.S.anova)

ROInum.GC.S.shRNA1<- glht(ROInum.GC.S.model1, mcp(shRNA1= "Tukey"))
summary(ROInum.GC.S.shRNA1)

ROInum.GC.S.shRNA2<- glht(ROInum.GC.S.model2, mcp(shRNA2= "Tukey"))
summary(ROInum.GC.S.shRNA2)






###############
# GCaMP fast vs delayed

df.ROInum.8strial.group1<-summarySE(ROInum.8strial.group, measurevar = "ROIsPerTrial", groupvars = c("shRNA1","Group"))
df.ROInum.8strial.group2<-summarySE(ROInum.8strial.group, measurevar = "ROIsPerTrial", groupvars = c("shRNA2","Group"))
df.ROInum.8strial.group.TP<-summarySE(ROInum.8strial.group, measurevar = "ROIsPerTrial", groupvars = c("shRNA1","Group","Timepoint"))

# plots


ggplot(df.ROInum.8strial.group1, aes(x=Group,y=ROIsPerTrial, fill= shRNA1)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("GCaMP ROIs")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(df.ROInum.8strial.group2, aes(x=Group,y=ROIsPerTrial, fill= shRNA2)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("GCaMP ROIs")+
  scale_fill_manual(values=cbbPalette) + 
  max.theme

#boxplots
ggplot(ROInum.8strial.group, aes(x=Group,y=ROIsPerTrial, fill= shRNA1)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(ROInum.8strial.group, aes(x=Group,y=ROIsPerTrial, fill= shRNA2)) +
  geom_boxplot()+
  ylab("num ROIs/trial per field of view") +
  ggtitle("GCaMP") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme



# only fast GcaMP

fastROInum<- subset(ROInum.8strial.group, Group=="fast")

ggplot(fastROInum, aes(x=shRNA1,y=ROIsPerTrial, fill= shRNA1)) +
  geom_boxplot()+
  ylab("Lck- fast ROIs per FOV") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(fastROInum, aes(x=shRNA2,y=ROIsPerTrial, fill= shRNA2)) +
  geom_boxplot()+
  ylab("Lck- fast ROIs per FOV") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(df.ROInum.8strial.group1[df.ROInum.8strial.group1$Group=="fast",], aes(x=shRNA1,y=ROIsPerTrial, fill=shRNA1)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=ROIsPerTrial-se, ymax=ROIsPerTrial+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("num ROIs/trial per field of view during 8s stim") +
  ggtitle("fast GCaMP ROIs") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

# only delayed GcaMP

delayedROInum<- subset(ROInum.8strial.group, Group=="delayed")

ggplot(delayedROInum, aes(x=shRNA1,y=ROIsPerTrial, fill= shRNA1)) +
  geom_boxplot()+
  ylab("Lck- delayed ROIs per FOV") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(delayedROInum, aes(x=shRNA2,y=ROIsPerTrial, fill= shRNA2)) +
  geom_boxplot()+
  ylab("Lck- delayed ROIs per FOV") +
  scale_fill_manual(values=cbbPalette) + 
  max.theme

ggplot(df.ROInum.8strial.group1[df.ROInum.8strial.group1$Group=="delayed",], aes(x=shRNA1,y=ROIsPerTrial, fill=shRNA1)) +
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
nROI.fast.day.model1 = lmer(ROIsPerTrial ~ shRNA1 + (1|Animal), fastROInum,REML=FALSE)
nROI.fast.day.model2 = lmer(ROIsPerTrial ~ shRNA2 + (1|Animal), fastROInum,REML=FALSE)
nROI.fast.day.anova <- anova(nROI.fast.day.null, nROI.fast.day.model1,
                             nROI.fast.day.model2)
print(nROI.fast.day.anova)

nROI.GC.shRNA1.fast<- glht(nROI.fast.day.model1, mcp(shRNA1= "Tukey"))
summary(nROI.GC.shRNA1.fast)

nROI.GC.shRNA2.fast<- glht(nROI.fast.day.model2, mcp(shRNA2= "Tukey"))
summary(nROI.GC.shRNA2.fast)


# delayed GCaMP only

nROI.delayed.day.null = lmer(ROIsPerTrial ~ (1|Animal), delayedROInum,REML=FALSE)
nROI.delayed.day.model1 = lmer(ROIsPerTrial ~ shRNA1 + (1|Animal), delayedROInum,REML=FALSE)
nROI.delayed.day.model2 = lmer(ROIsPerTrial ~ shRNA2 + (1|Animal), delayedROInum,REML=FALSE)
nROI.delayed.day.anova <- anova(nROI.delayed.day.null, nROI.delayed.day.model1,
                             nROI.delayed.day.model2)
print(nROI.delayed.day.anova)

nROI.GC.shRNA1.delayed<- glht(nROI.delayed.day.model1, mcp(shRNA1= "Tukey"))
summary(nROI.GC.shRNA1.delayed)

nROI.GC.shRNA2.delayed<- glht(nROI.delayed.day.model2, mcp(shRNA2= "Tukey"))
summary(nROI.GC.shRNA2.delayed)



# wilcox test for median compairsons

wilcox.test(fastROInum$ROIsPerTrial[fastROInum$shRNA1=="Control"], 
            fastROInum$ROIsPerTrial[fastROInum$shRNA1=="KD"])

wilcox.test(fastROInum$ROIsPerTrial[fastROInum$shRNA2=="Control"], 
            fastROInum$ROIsPerTrial[fastROInum$shRNA2=="KD"])




