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
baseline1_65626 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 65626/Diameter/2019_11_25 Baseline/2019_11_25_SPOT1_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline2_65626 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 65626/Diameter/2019_11_25 Baseline/2019_11_25_SPOT2_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline3_65626<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 65626/Diameter/2019_11_25 Baseline/2019_11_25_SPOT2_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline4_65626<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 65626/Diameter/2019_11_25 Baseline/2019_11_25_SPOT2_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline5_65626 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 65626/Diameter/2019_11_29 Baseline/2019_11_29_SPOT4_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline6_65626 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 65626/Diameter/2019_11_29 Baseline/2019_11_29_SPOT4_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline7_65626<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 65626/Diameter/2019_11_29 Baseline/2019_11_29_SPOT4_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline8_65626<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 65626/Diameter/2019_11_29 Baseline/2019_11_29_SPOT5_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline9_65626 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 65626/Diameter/2019_11_29 Baseline/2019_11_29_SPOT5_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline10_65626 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 65626/Diameter/2019_11_29 Baseline/2019_11_29_SPOT5_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline11_65626<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 65626/Diameter/2019_11_29 Baseline/2019_11_29_SPOT6_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")

baseline1_79355 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2019_11_05/2019_11_05_SPOT1_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline2_79355 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2019_11_05/2019_11_05_SPOT1_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline3_79355<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2019_11_05/2019_11_05_SPOT1_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline4_79355<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2019_11_05/2019_11_05_SPOT2_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline5_79355 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2019_11_13/2019_11_13_SPOT3_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline6_79355 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2019_11_13/2019_11_13_SPOT3_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline7_79355<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2019_11_13/2019_11_13_SPOT3_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline8_79355<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2019_11_13/2019_11_13_SPOT3_v4_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline9_79355 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2019_11_13/2019_11_13_SPOT4_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline10_79355 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2019_11_13/2019_11_13_SPOT4_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
baseline11_79355 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2019_12_06/2019_12_06_SPOT5_v1_baseline_diameterFWHM.csv", header=TRUE, sep = ",")

baseline1_82093 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 82093/Diameter/2020_01_07/spot1_v1_diameterFWHM.csv", header=TRUE, sep = ",")
baseline2_82093 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 82093/Diameter/2020_01_07/spot2_v1_diameterFWHM.csv", header=TRUE, sep = ",")
baseline3_82093<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 82093/Diameter/2020_01_07/spot2_v2_diameterFWHM.csv", header=TRUE, sep = ",")
baseline4_82093<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 82093/Diameter/2020_01_07/spot3_v1_diameterFWHM.csv", header=TRUE, sep = ",")
baseline5_82093 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 82093/Diameter/2020_01_07/spot3_v2_diameterFWHM.csv", header=TRUE, sep = ",")
baseline6_82093 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 82093/Diameter/2020_01_07/spot3_v3_diameterFWHM.csv", header=TRUE, sep = ",")

#aging
month11_1_65626 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 65626/Diameter/2020_01_09/Spot1_v1_diameterFWHM.csv", header=TRUE, sep = ",")
month11_2_65626 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 65626/Diameter/2020_01_09/Spot2_v1_diameterFWHM.csv", header=TRUE, sep = ",")
month11_3_65626<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 65626/Diameter/2020_01_09/Spot2_v2_diameterFWHM.csv", header=TRUE, sep = ",")
month11_4_65626<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 65626/Diameter/2020_01_09/Spot2_v3_diameterFWHM.csv", header=TRUE, sep = ",")

#nimodipine
nim1_79355 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2020_01_17_nimodipine/spot1_v1_diameterFWHM.csv", header=TRUE, sep = ",")
nim2_79355 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2020_01_17_nimodipine/spot1_v2_diameterFWHM.csv", header=TRUE, sep = ",")
nim3_79355<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2020_01_17_nimodipine/spot1_v3_diameterFWHM.csv", header=TRUE, sep = ",")
nim4_79355<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2020_01_17_nimodipine/spot2_v1_diameterFWHM.csv", header=TRUE, sep = ",")
nim5_79355 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2020_01_17_nimodipine/spot5_v1_diameterFWHM.csv", header=TRUE, sep = ",")

nim1_82093 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 82093/Diameter/2020_01_17_nim_noTR/spot2_v2_short_diameterFWHM.csv", header=TRUE, sep = ",")
nim2_82093 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 82093/Diameter/2020_01_17_nim_noTR/spot3_v3_short_movement_diameterFWHM.csv", header=TRUE, sep = ",")

#Pyr3
Pyr1_79355 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2020_01_14_Pyr3/spot1_v1_diameterFWHM.csv", header=TRUE, sep = ",")
Pyr2_79355 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2020_01_14_Pyr3/spot1_v2_diameterFWHM.csv", header=TRUE, sep = ",")
Pyr3_79355<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2020_01_14_Pyr3/spot1_v3_diameterFWHM.csv", header=TRUE, sep = ",")
Pyr4_79355<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2020_01_14_Pyr3/spot2_v1_diameterFWHM.csv", header=TRUE, sep = ",")
Pyr5_79355 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 79335/Diameter/2020_01_14_Pyr3/spot5_v1-bad_diameterFWHM.csv", header=TRUE, sep = ",")

Pyr1_82093 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 82093/Diameter/2020_01_14_Pyr3/spot1_v1_diameterFWHM.csv", header=TRUE, sep = ",")
Pyr2_82093 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 82093/Diameter/2020_01_14_Pyr3/spot2_v1_diameterFWHM.csv", header=TRUE, sep = ",")
Pyr3_82093<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 82093/Diameter/2020_01_14_Pyr3/spot2_v2_diameterFWHM.csv", header=TRUE, sep = ",")
Pyr4_82093<- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/GCaMP Mice/Results/Bloodflow Results/GCaMP 82093/Diameter/2020_01_14_Pyr3/spot3_v1_movement_diameterFWHM.csv", header=TRUE, sep = ",")

#####
BL1_65<-summarySE(baseline1_65626, measurevar="diameter", na.rm = TRUE)
BL1_65$Animal<-"65626"
BL1_65$Vessel<-"Spot1_v1"
BL1_65$Drug<-"baseline"

BL2_65<-summarySE(baseline2_65626, measurevar="diameter", na.rm = TRUE)
BL2_65$Animal<-"65626"
BL2_65$Vessel<-"Spot2_v1"
BL2_65$Drug<-"baseline"

BL3_65<-summarySE(baseline3_65626, measurevar="diameter", na.rm = TRUE)
BL3_65$Animal<-"65626"
BL3_65$Vessel<-"Spot2_v2"
BL3_65$Drug<-"baseline"

BL4_65<-summarySE(baseline4_65626, measurevar="diameter", na.rm = TRUE)
BL4_65$Animal<-"65626"
BL4_65$Vessel<-"Spot2_v3"
BL4_65$Drug<-"baseline"

BL5_65<-summarySE(baseline5_65626, measurevar="diameter", na.rm = TRUE)
BL5_65$Animal<-"65626"
BL5_65$Vessel<-"Spot4_v1"
BL5_65$Drug<-"baseline"

BL6_65<-summarySE(baseline6_65626, measurevar="diameter", na.rm = TRUE)
BL6_65$Animal<-"65626"
BL6_65$Vessel<-"Spot4_v2"
BL6_65$Drug<-"baseline"

BL7_65<-summarySE(baseline7_65626, measurevar="diameter", na.rm = TRUE)
BL7_65$Animal<-"65626"
BL7_65$Vessel<-"Spot4_v3"
BL7_65$Drug<-"baseline"

BL8_65<-summarySE(baseline8_65626, measurevar="diameter", na.rm = TRUE)
BL8_65$Animal<-"65626"
BL8_65$Vessel<-"Spot5_v1"
BL8_65$Drug<-"baseline"

BL9_65<-summarySE(baseline9_65626, measurevar="diameter", na.rm = TRUE)
BL9_65$Animal<-"65626"
BL9_65$Vessel<-"Spot5_v2"
BL9_65$Drug<-"baseline"

BL10_65<-summarySE(baseline10_65626, measurevar="diameter", na.rm = TRUE)
BL10_65$Animal<-"65626"
BL10_65$Vessel<-"Spot5_v3"
BL10_65$Drug<-"baseline"

BL11_65<-summarySE(baseline11_65626, measurevar="diameter", na.rm = TRUE)
BL11_65$Animal<-"65626"
BL11_65$Vessel<-"Spot6_v1"
BL11_65$Drug<-"baseline"

BL1_79<-summarySE(baseline1_79355, measurevar="diameter", na.rm = TRUE)
BL1_79$Animal<-"79355"
BL1_79$Vessel<-"Spot1_v1"
BL1_79$Drug<-"baseline"

BL2_79<-summarySE(baseline2_79355, measurevar="diameter", na.rm = TRUE)
BL2_79$Animal<-"79355"
BL2_79$Vessel<-"Spot1_v2"
BL2_79$Drug<-"baseline"

BL3_79<-summarySE(baseline3_79355, measurevar="diameter", na.rm = TRUE)
BL3_79$Animal<-"79355"
BL3_79$Vessel<-"Spot1_v3"
BL3_79$Drug<-"baseline"

BL4_79<-summarySE(baseline4_79355, measurevar="diameter", na.rm = TRUE)
BL4_79$Animal<-"79355"
BL4_79$Vessel<-"Spot2_v1"
BL4_79$Drug<-"baseline"

BL5_79<-summarySE(baseline5_79355, measurevar="diameter", na.rm = TRUE)
BL5_79$Animal<-"79355"
BL5_79$Vessel<-"Spot3_v1"
BL5_79$Drug<-"baseline"

BL6_79<-summarySE(baseline6_79355, measurevar="diameter", na.rm = TRUE)
BL6_79$Animal<-"79355"
BL6_79$Vessel<-"Spot3_v2"
BL6_79$Drug<-"baseline"

BL7_79<-summarySE(baseline7_79355, measurevar="diameter", na.rm = TRUE)
BL7_79$Animal<-"79355"
BL7_79$Vessel<-"Spot3_v3"
BL7_79$Drug<-"baseline"

BL8_79<-summarySE(baseline8_79355, measurevar="diameter", na.rm = TRUE)
BL8_79$Animal<-"79355"
BL8_79$Vessel<-"Spot3_v4"
BL8_79$Drug<-"baseline"

BL9_79<-summarySE(baseline9_79355, measurevar="diameter", na.rm = TRUE)
BL9_79$Animal<-"79355"
BL9_79$Vessel<-"Spot4_v1"
BL9_79$Drug<-"baseline"

BL10_79<-summarySE(baseline10_79355, measurevar="diameter", na.rm = TRUE)
BL10_79$Animal<-"79355"
BL10_79$Vessel<-"Spot4_v2"
BL10_79$Drug<-"baseline"

BL11_79<-summarySE(baseline11_79355, measurevar="diameter", na.rm = TRUE)
BL11_79$Animal<-"79355"
BL11_79$Vessel<-"Spot5_v1"
BL11_79$Drug<-"baseline"

BL1_82<-summarySE(baseline1_82093, measurevar="diameter", na.rm = TRUE)
BL1_82$Animal<-"82093"
BL1_82$Vessel<-"Spot1_v1"
BL1_82$Drug<-"baseline"

BL2_82<-summarySE(baseline2_82093, measurevar="diameter", na.rm = TRUE)
BL2_82$Animal<-"82093"
BL2_82$Vessel<-"Spot2_v1"
BL2_82$Drug<-"baseline"

BL3_82<-summarySE(baseline3_82093, measurevar="diameter", na.rm = TRUE)
BL3_82$Animal<-"82093"
BL3_82$Vessel<-"Spot2_v2"
BL3_82$Drug<-"baseline"

BL4_82<-summarySE(baseline4_82093, measurevar="diameter", na.rm = TRUE)
BL4_82$Animal<-"82093"
BL4_82$Vessel<-"Spot3_v1"
BL4_82$Drug<-"baseline"

BL5_82<-summarySE(baseline5_82093, measurevar="diameter", na.rm = TRUE)
BL5_82$Animal<-"82093"
BL5_82$Vessel<-"Spot3_v2"
BL5_82$Drug<-"baseline"

BL6_82<-summarySE(baseline6_82093, measurevar="diameter", na.rm = TRUE)
BL6_82$Animal<-"82093"
BL6_82$Vessel<-"Spot3_v3"
BL6_82$Drug<-"baseline"

#nimodipine

nim1_82<-summarySE(nim2_82093, measurevar="diameter", na.rm = TRUE)
nim1_82$Animal<-"82093"
nim1_82$Vessel<-"Spot3_v3"
nim1_82$Drug<-"nimodipine"

nim2_82<-summarySE(nim1_82093, measurevar="diameter", na.rm = TRUE)
nim2_82$Animal<-"82093"
nim2_82$Vessel<-"Spot2_v2"
nim2_82$Drug<-"nimodipine"

nim1_79<-summarySE(nim1_79355, measurevar="diameter", na.rm = TRUE)
nim1_79$Animal<-"79355"
nim1_79$Vessel<-"Spot1_v1"
nim1_79$Drug<-"nimodipine"

nim2_79<-summarySE(nim2_79355, measurevar="diameter", na.rm = TRUE)
nim2_79$Animal<-"79355"
nim2_79$Vessel<-"Spot1_v2"
nim2_79$Drug<-"nimodipine"

nim3_79<-summarySE(nim3_79355, measurevar="diameter", na.rm = TRUE)
nim3_79$Animal<-"79355"
nim3_79$Vessel<-"Spot1_v3"
nim3_79$Drug<-"nimodipine"

nim4_79<-summarySE(nim4_79355, measurevar="diameter", na.rm = TRUE)
nim4_79$Animal<-"79355"
nim4_79$Vessel<-"Spot2_v1"
nim4_79$Drug<-"nimodipine"

nim5_79<-summarySE(nim5_79355, measurevar="diameter", na.rm = TRUE)
nim5_79$Animal<-"79355"
nim5_79$Vessel<-"Spot5_v1"
nim5_79$Drug<-"nimodipine"


# Pyr3

Pyr1_79<-summarySE(Pyr1_79355, measurevar="diameter", na.rm = TRUE)
Pyr1_79$Animal<-"79355"
Pyr1_79$Vessel<-"Spot1_v1"
Pyr1_79$Drug<-"Pyr3"

Pyr2_79<-summarySE(Pyr2_79355, measurevar="diameter", na.rm = TRUE)
Pyr2_79$Animal<-"79355"
Pyr2_79$Vessel<-"Spot1_v2"
Pyr2_79$Drug<-"Pyr3"

Pyr3_79<-summarySE(Pyr3_79355, measurevar="diameter", na.rm = TRUE)
Pyr3_79$Animal<-"79355"
Pyr3_79$Vessel<-"Spot1_v3"
Pyr3_79$Drug<-"Pyr3"

Pyr4_79<-summarySE(Pyr4_79355, measurevar="diameter", na.rm = TRUE)
Pyr4_79$Animal<-"79355"
Pyr4_79$Vessel<-"Spot2_v1"
Pyr4_79$Drug<-"Pyr3"

Pyr5_79<-summarySE(Pyr5_79355, measurevar="diameter", na.rm = TRUE)
Pyr5_79$Animal<-"79355"
Pyr5_79$Vessel<-"Spot5_v1"
Pyr5_79$Drug<-"Pyr3"


Pyr1_82<-summarySE(Pyr1_82093, measurevar="diameter", na.rm = TRUE)
Pyr1_82$Animal<-"82093"
Pyr1_82$Vessel<-"Spot1_v1"
Pyr1_82$Drug<-"Pyr3"


Pyr2_82<-summarySE(Pyr2_82093, measurevar="diameter", na.rm = TRUE)
Pyr2_82$Animal<-"82093"
Pyr2_82$Vessel<-"Spot2_v1"
Pyr2_82$Drug<-"Pyr3"

Pyr3_82<-summarySE(Pyr3_82093, measurevar="diameter", na.rm = TRUE)
Pyr3_82$Animal<-"82093"
Pyr3_82$Vessel<-"Spot2_v2"
Pyr3_82$Drug<-"Pyr3"

Pyr4_82<-summarySE(Pyr4_82093, measurevar="diameter", na.rm = TRUE)
Pyr4_82$Animal<-"82093"
Pyr4_82$Vessel<-"Spot3_v3"
Pyr4_82$Drug<-"Pyr3"



#
MT1_65<-summarySE(month11_1_65626, measurevar="diameter", na.rm = TRUE)
MT1_65$Animal<-"65626"
MT1_65$Vessel<-"Spot1_v1"
MT1_65$Drug<-"Month11"

MT2_65<-summarySE(month11_2_65626, measurevar="diameter", na.rm = TRUE)
MT2_65$Animal<-"65626"
MT2_65$Vessel<-"Spot2_v1"
MT2_65$Drug<-"Month11"

MT3_65<-summarySE(month11_3_65626, measurevar="diameter", na.rm = TRUE)
MT3_65$Animal<-"65626"
MT3_65$Vessel<-"Spot2_v2"
MT3_65$Drug<-"Month11"

MT4_65<-summarySE(month11_4_65626, measurevar="diameter", na.rm = TRUE)
MT4_65$Animal<-"65626"
MT4_65$Vessel<-"Spot2_v3"
MT4_65$Drug<-"Month11"

#combine data together for each treatment
baseline1<-rbind(BL1_65,BL2_65,BL3_65,BL4_65,BL5_65,BL6_65,BL7_65,BL8_65,BL9_65,BL10_65,BL11_65)
baseline2<-rbind(BL1_79,BL2_79,BL3_79,BL4_79,BL5_79,BL6_79,BL7_79,BL8_79,BL9_79,BL10_79,BL11_79,
                BL1_82,BL2_82,BL3_82,BL4_82,BL5_82,BL6_82) 
nimodipine<-rbind(nim1_79, nim2_79, nim3_79, nim4_79, nim5_79, nim1_82, nim2_82)#L_nimodipine1)
pyr3<-rbind(Pyr1_79, Pyr2_79, Pyr3_79, Pyr4_79, Pyr5_79, Pyr1_82, Pyr2_82, Pyr3_82, Pyr4_82)
month11<-rbind(MT1_65,MT2_65,MT3_65,MT4_65)

#Drug data
DrugData<-rbind(baseline2,nimodipine,pyr3)
DrugData$Drug<-as.factor(DrugData$Drug)

#aging data
AgeData<-rbind(baseline1,month11)
AgeData$Drug<-as.factor(AgeData$Drug)

#pull out the Spot name
pos= regexpr('Spot', DrugData$Vessel)
DrugData$SpotName<-substr(DrugData$Vessel,pos, pos+4)

pos2= regexpr('Spot', AgeData$Vessel)
AgeData$SpotName<-substr(AgeData$Vessel,pos2, pos2+4)

###############################
#histograms
ggplot(DrugData, aes(x=diameter, fill=Drug))+ geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("Distribution of Diameter")

#########

#diameter
df1A <- summarySE(DrugData, measurevar="diameter", groupvars=c("Drug"))
df1B <- summarySE(DrugData, measurevar="diameter", groupvars=c("Drug","Animal"))
df1C <- summarySE(AgeData, measurevar="diameter", groupvars=c("Drug"))

ggplot(data=df1A, aes(x=Drug, y=diameter, fill=Drug)) +
  geom_errorbar(aes(ymin=diameter-se, ymax=diameter+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("diameter") +
  scale_fill_manual(
    values=c("black", "red", "blue", "green")) +
  max.theme

ggplot(data=df1B, aes(x=Animal, y=diameter, fill=Drug)) +
  geom_errorbar(aes(ymin=diameter-se, ymax=diameter+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("diameter") +
  scale_fill_manual(
    values=c("black", "red", "blue", "green")) +
  max.theme

ggplot(data=df1C, aes(x=Drug, y=diameter, fill=Drug)) +
  geom_errorbar(aes(ymin=diameter-se, ymax=diameter+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("diameter") +
  scale_fill_manual(
    values=c("black", "red", "blue", "green")) +
  max.theme

######
#stats

#diameter
diameter.null = lmer(diameter ~ (1|Animal) + (1|SpotName) + (1|Vessel), DrugData,REML=FALSE)
diameter.model1 = lmer(diameter~ Drug + (1|Animal) + (1|SpotName) + (1|Vessel), DrugData,REML=FALSE)
diameter.anova <- anova(diameter.null, diameter.model1)
print(diameter.anova)

# p values
diam.drug <- lsmeans(diameter.model1, pairwise ~ Drug, glhargs=list())
summary(diam.drug)

