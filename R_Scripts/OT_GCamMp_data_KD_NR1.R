library("lme4")
library("lmerTest")
library("lattice")
library("plyr")
library("dplyr")
library("ggplot2")
library("lsmeans")
library("Rmisc")
#library("MASS")
library("multcomp")
library("reshape2")
library("tidyr")
#library("data.table")
library("Hmisc")
library("stringr")

#################
#Pharmacology
# Trazodone, Prazosin, Atropine, Metergoline, DSP4

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
cbbPalette <- c("#000000","#D55E00","#009E73","#E69F00","#56B4E9","#CC79A7","#F0E442")

########################
# load data

# peak data
long.peaks.74<-read.table("C:\Data\Results\FilesforR\74_peaks_longtrials.csv",  header=TRUE, sep = ",")
long.peaks.92<-read.table("C:\Data\Results\FilesforR\92_peaks_longtrials.csv",  header=TRUE, sep = ",")
long.peaks.94<-read.table("C:\Data\Results\FilesforR\94_peaks_longtrials.csv",  header=TRUE, sep = ",")
long.peaks.95<-read.table("C:\Data\Results\FilesforR\95_peaks_longtrials.csv",  header=TRUE, sep = ",")
long.peaks.96<-read.table("C:\Data\Results\FilesforR\96_peaks_longtrials.csv",  header=TRUE, sep = ",")
long.peaks.Alice<-read.table("C:\Data\Results\FilesforR\Alice_peaks_longtrials.csv",  header=TRUE, sep = ",")
long.peaks.crazy8<-read.table("C:\Data\Results\FilesforR\crazy8_peaks_longtrials.csv",  header=TRUE, sep = ",")

# OT data
long.OT.74<-read.table("C:\Data\Results\FilesforR\74_onset_time_longtrials.csv",  header=TRUE, sep = ",")
long.OT.92<-read.table("C:\Data\Results\FilesforR\92_onset_time_longtrials.csv",  header=TRUE, sep = ",")
long.OT.94<-read.table("C:\Data\Results\FilesforR\94_onset_time_longtrials.csv",  header=TRUE, sep = ",")
long.OT.95<-read.table("C:\Data\Results\FilesforR\95_onset_time_longtrials.csv",  header=TRUE, sep = ",")
long.OT.96<-read.table("C:\Data\Results\FilesforR\96_onset_time_longtrials.csv",  header=TRUE, sep = ",")
long.OT.Alice<-read.table("C:\Data\Results\FilesforR\Alice_onset_time_longtrials.csv",  header=TRUE, sep = ",")
long.OT.crazy8<-read.table("C:\Data\Results\FilesforR\crazy8_onset_time_longtrials.csv",  header=TRUE, sep = ",")

#field data
long.field.74<-read.table("C:\Data\Results\FilesforR\74_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
long.field.92<-read.table("C:\Data\Results\FilesforR\92_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
long.field.94<-read.table("C:\Data\Results\FilesforR\94_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
long.field.95<-read.table("C:\Data\Results\FilesforR\95_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
long.field.96<-read.table("C:\Data\Results\FilesforR\96_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
long.field.Alice<-read.table("C:\Data\Results\FilesforR\Alice_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
long.field.crazy8<-read.table("C:\Data\Results\FilesforR\crazy8_Lck_field_longtrials.csv",  header=TRUE, sep = ",")
##########

lsm.options(pbkrtest.limit = 100000)

#  merge data
all.OT<-rbind(long.OT.74, long.OT.92, long.OT.94,long.OT.95, long.OT.96, long.OT.Alice, long.OT.crazy8)
all.peaks<-rbind(long.peaks.74, long.peaks.92, long.peaks.94,long.peaks.95, long.peaks.96, long.peaks.Alice, long.peaks.crazy8)
all.field<-rbind(long.field.74, long.field.92, long.field.94,long.field.95, long.field.96, long.field.Alice, long.field.crazy8)

# need to define shRNA

#unique ROI names #needhelp and
all.OT$ROIs_trial<-paste(all.OT$Animal, all.OT$Spot, all.OT$Trial,all.OT$ROI, sep= "_")
all.OT$Spot_trial_shRNA<-paste(all.OT$Animal, all.OT$Spot, all.OT$Trial,all.OT$shRNA, sep= "_")
all.OT$Spot_trial_shRNA_Cond<-paste(all.OT$Spot_trial_shRNA, all.OT$Condition, sep="_")
all.OT$ROIs_Cond_shRNA<-paste(all.OT$ROIs_trial, all.OT$Condition, all.OT$shRNA, sep="_")


# REMOVE duplicate entries from onset time and peak time data frames #needhelp
# only the first entry will be used
all.OT2=all.OT[order(all.OT$OnsetTime),] # sort by ascending onset time

# remove duplicate entries (in theory only the first and therefore fastest onset times will remain) #needhelp
all.OT<-distinct(all.OT2, ROIs_Cond_shRNA, .keep_all = TRUE)

#remove entries with no onset time

all.OT<-subset(all.OT, OnsetTime!="NaN")
all.OT$ROIType= "none"
all.OT.astrocyte<- subset(all.OT, Channel=="GCaMP")
all.OT.neuron<- subset(all.OT, Channel=="RCaMP")



all.OT<-rbind(all.OTA, all.OTB)
all.OT$ROIType<- as.factor(all.OT$ROIType)

rm(all.OT2, all.OTA, all.OTB)