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
RR_baseline1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_11/2019_07_11_SPOT3_v1_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_baseline2 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_11/2019_07_11_SPOT3_v4_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_baseline3 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_22/2019_07_22_SPOT4_v1_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_baseline4 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_22/2019_07_22_SPOT4_v2_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_baseline5 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_22/2019_07_22_SPOT4_v3_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_baseline6 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_22/2019_07_22_SPOT4_v4_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")

MJ_baseline1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_06_11 - KX/2019_06_11_SPOT2_v1_vel_Linescan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline2 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_06_11 - KX/2019_06_11_SPOT2_v2_vel_095_Linescan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline3 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_06_11 - KX/2019_06_11_SPOT2_v2_vel_096_Linescan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline4 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_06_11 - KX/2019_06_11_SPOT2_v3_vel_099_Linescan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline5 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_06_11 - KX/2019_06_11_SPOT2_v3_vel_100_Linescan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline6 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_06_11 - KX/2019_06_11_SPOT4_v1_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline7 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_09/2019_07_09_SPOT1_v1_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline8 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_09/2019_07_09_SPOT1_v2_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline9 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_09/2019_07_09_SPOT1_v3_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline10 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_11/2019_07_11_SPOT3_v1_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline11 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_22/2019_07_22_SPOT4_v1_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline12 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_22/2019_07_22_SPOT4_v2_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline13 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_22/2019_07_22_SPOT5_v1_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline14 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_22/2019_07_22_SPOT5_v2_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline15 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_22/2019_07_22_SPOT5_v3_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline16 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_22/2019_07_22_SPOT5_v4_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline17 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_22/2019_07_22_SPOT5_v5_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline18 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_22/2019_07_22_SPOT5_v6_vel_104_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline19 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_22/2019_07_22_SPOT5_v6_vel_106_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline20 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_19/2019_08_19_SPOT1_v1_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline21 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_19/2019_08_19_SPOT1_v2_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline22 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_19/2019_08_19_SPOT1_v3_vel_098_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline23 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_19/2019_08_19_SPOT3_v1_vel_077_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline24 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_19/2019_08_19_SPOT3_v1_vel_078_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_baseline25 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_19/2019_08_19_SPOT6_v1_vel_LineScan_velocityLSPIV.csv", header=TRUE, sep = ",")

#nimodipine
RR_nimodipine1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_16 Nimodipine/2019_07_16_SPOT1_v1_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_nimodipine2 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_16 Nimodipine/2019_07_16_SPOT2_v1_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_nimodipine3 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_16 Nimodipine/2019_07_16_SPOT2_v2_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_nimodipine4 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_16 Nimodipine/2019_07_16_SPOT2_v3_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_nimodipine5 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_16 Nimodipine/2019_07_16_SPOT3_v1_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_nimodipine6 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_16 Nimodipine/2019_07_16_SPOT3_v2_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_nimodipine7 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_16 Nimodipine/2019_07_16_SPOT3_v3_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_nimodipine8 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT4_v1_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_nimodipine9 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT4_v2_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_nimodipine10 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT4_v3_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_nimodipine11 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT4_v4_vel_098_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_nimodipine12 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT5_v1_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_nimodipine13 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT5_v2_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_nimodipine14 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT5_v3_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_nimodipine15 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT5_v4_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_nimodipine16 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_16 Nimodipine/2019_07_16_SPOT3_v4_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")

MJ_nimodipine1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_16 Nimodipine/2019_07_16_SPOT3_v1_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_nimodipine2 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT3_v1_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_nimodipine3 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT4_v1_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_nimodipine4 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT4_v2_vel_115_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_nimodipine5 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT4_v2_vel_116_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_nimodipine6 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT5_v1_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_nimodipine7 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT5_v2_vel_120_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_nimodipine8 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT5_v2_vel_122_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_nimodipine9 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT5_v3_vel_123_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_nimodipine10 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT5_v3_vel_124_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_nimodipine11 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT5_v4_vel_128_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_nimodipine12 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT5_v5_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_nimodipine13 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT5_v6_vel_141_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_nimodipine14 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_02 Nimodipine/2019_08_02_SPOT5_v6_vel_142_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_nimodipine15 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_16 Nimodipine/2019_07_16_SPOT6_v1_vel_LineScan_nimodipine_velocityLSPIV.csv", header=TRUE, sep = ",")

#pyr3
RR_pyr3_1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_19 PYR3/2019_07_19_SPOT1_v1_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_pyr3_2 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_19 PYR3/2019_07_19_SPOT1_v2_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_pyr3_3 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_19 PYR3/2019_07_19_SPOT1_v3_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_pyr3_4 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_19 PYR3/2019_07_19_SPOT1_v4_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_pyr3_5 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_19 PYR3/2019_07_19_SPOT2_v1_vel_096_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_pyr3_6 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_19 PYR3/2019_07_19_SPOT2_v1_vel_097_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_pyr3_7 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_19 PYR3/2019_07_19_SPOT2_v2_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_pyr3_8 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_19 PYR3/2019_07_19_SPOT3_v1_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_pyr3_9 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_19 PYR3/2019_07_19_SPOT3_v2_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
RR_pyr3_10 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Velocity/2019_07_19 PYR3/2019_07_19_SPOT3_v3_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")

MJ_pyr3_1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_19 PYR3/2019_07_19_SPOT1_v1_vel_080_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_pyr3_2 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_19 PYR3/2019_07_19_SPOT1_v1_vel_081_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_pyr3_3 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_19 PYR3/2019_07_19_SPOT1_v2_vel_083_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_pyr3_4 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_19 PYR3/2019_07_19_SPOT1_v2_vel_084_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_pyr3_5 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_07_19 PYR3/2019_07_19_SPOT1_v3_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_pyr3_6 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_14 PYR3/2019_08_14_SPOT3_v1_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_pyr3_7 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_14 PYR3/2019_08_14_SPOT4_v1_vel_074_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_pyr3_8 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_14 PYR3/2019_08_14_SPOT4_v1_vel_075_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_pyr3_9 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_14 PYR3/2019_08_14_SPOT4_v2_vel_079_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_pyr3_10 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_14 PYR3/2019_08_14_SPOT5_v1_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_pyr3_11 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_14 PYR3/2019_08_14_SPOT5_v2_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_pyr3_12 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_14 PYR3/2019_08_14_SPOT5_v3_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_pyr3_13 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_14 PYR3/2019_08_14_SPOT5_v4_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_pyr3_14 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_14 PYR3/2019_08_14_SPOT5_v5_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_pyr3_15 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_14 PYR3/2019_08_14_SPOT5_v6_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")
MJ_pyr3_16 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Velocity/2019_08_14 PYR3/2019_08_14_SPOT6_v1_vel_LineScan_pyr3_velocityLSPIV.csv", header=TRUE, sep = ",")


RR_BL1<-summarySE(RR_baseline1, measurevar="velocity", na.rm = TRUE)
RR_BL1$Animal<-"RR"
RR_BL1$Vessel<-"SPOT3_v1"
RR_BL1$Drug<-"baseline"


RR_BL2<-summarySE(RR_baseline2, measurevar="velocity", na.rm = TRUE)
RR_BL2$Animal<-"RR"
RR_BL2$Vessel<-"SPOT3_v4"
RR_BL2$Drug<-"baseline"

RR_BL3<-summarySE(RR_baseline3, measurevar="velocity", na.rm = TRUE)
RR_BL3$Animal<-"RR"
RR_BL3$Vessel<-"SPOT4_v1"
RR_BL3$Drug<-"baseline"

RR_BL4<-summarySE(RR_baseline4, measurevar="velocity", na.rm = TRUE)
RR_BL4$Animal<-"RR"
RR_BL4$Vessel<-"SPOT4_v2"
RR_BL4$Drug<-"baseline"

RR_BL5<-summarySE(RR_baseline5, measurevar="velocity", na.rm = TRUE)
RR_BL5$Animal<-"RR"
RR_BL5$Vessel<-"SPOT4_v3"
RR_BL5$Drug<-"baseline"

RR_BL6<-summarySE(RR_baseline6, measurevar="velocity", na.rm = TRUE)
RR_BL6$Animal<-"RR"
RR_BL6$Vessel<-"SPOT4_v4"
RR_BL6$Drug<-"baseline"

RR_NP1<-summarySE(RR_nimodipine1, measurevar="velocity", na.rm = TRUE)
RR_NP1$Animal<-"RR"
RR_NP1$Vessel<-"SPOT1_v1"
RR_NP1$Drug<-"nimodipine"

RR_NP2<-summarySE(RR_nimodipine2, measurevar="velocity", na.rm = TRUE)
RR_NP2$Animal<-"RR"
RR_NP2$Vessel<-"SPOT2_v1"
RR_NP2$Drug<-"nimodipine"

RR_NP3<-summarySE(RR_nimodipine3, measurevar="velocity", na.rm = TRUE)
RR_NP3$Animal<-"RR"
RR_NP3$Vessel<-"SPOT2_v2"
RR_NP3$Drug<-"nimodipine"

RR_NP4<-summarySE(RR_nimodipine4, measurevar="velocity", na.rm = TRUE)
RR_NP4$Animal<-"RR"
RR_NP4$Vessel<-"SPOT2_v3"
RR_NP4$Drug<-"nimodipine"

RR_NP5<-summarySE(RR_nimodipine5, measurevar="velocity", na.rm = TRUE)
RR_NP5$Animal<-"RR"
RR_NP5$Vessel<-"SPOT3_v1"
RR_NP5$Drug<-"nimodipine"

RR_NP6<-summarySE(RR_nimodipine6, measurevar="velocity", na.rm = TRUE)
RR_NP6$Animal<-"RR"
RR_NP6$Vessel<-"SPOT3_v2"
RR_NP6$Drug<-"nimodipine"

RR_NP7<-summarySE(RR_nimodipine7, measurevar="velocity", na.rm = TRUE)
RR_NP7$Animal<-"RR"
RR_NP7$Vessel<-"SPOT3_v3"
RR_NP7$Drug<-"nimodipine"

RR_NP8<-summarySE(RR_nimodipine8, measurevar="velocity", na.rm = TRUE)
RR_NP8$Animal<-"RR"
RR_NP8$Vessel<-"SPOT4_v1"
RR_NP8$Drug<-"nimodipine"

RR_NP9<-summarySE(RR_nimodipine9, measurevar="velocity", na.rm = TRUE)
RR_NP9$Animal<-"RR"
RR_NP9$Vessel<-"SPOT4_v2"
RR_NP9$Drug<-"nimodipine"

RR_NP10<-summarySE(RR_nimodipine10, measurevar="velocity", na.rm = TRUE)
RR_NP10$Animal<-"RR"
RR_NP10$Vessel<-"SPOT4_v3"
RR_NP10$Drug<-"nimodipine"

RR_NP11<-summarySE(RR_nimodipine11, measurevar="velocity", na.rm = TRUE)
RR_NP11$Animal<-"RR"
RR_NP11$Vessel<-"SPOT4_v4"
RR_NP11$Drug<-"nimodipine"

RR_NP12<-summarySE(RR_nimodipine12, measurevar="velocity", na.rm = TRUE)
RR_NP12$Animal<-"RR"
RR_NP12$Vessel<-"SPOT5_v1"
RR_NP12$Drug<-"nimodipine"

RR_NP13<-summarySE(RR_nimodipine13, measurevar="velocity", na.rm = TRUE)
RR_NP13$Animal<-"RR"
RR_NP13$Vessel<-"SPOT5_v2"
RR_NP13$Drug<-"nimodipine"

RR_NP14<-summarySE(RR_nimodipine14, measurevar="velocity", na.rm = TRUE)
RR_NP14$Animal<-"RR"
RR_NP14$Vessel<-"SPOT5_v3"
RR_NP14$Drug<-"nimodipine"

RR_NP15<-summarySE(RR_nimodipine15, measurevar="velocity", na.rm = TRUE)
RR_NP15$Animal<-"RR"
RR_NP15$Vessel<-"SPOT5_v4"
RR_NP15$Drug<-"nimodipine"

RR_NP16<-summarySE(RR_nimodipine16, measurevar="velocity", na.rm = TRUE)
RR_NP16$Animal<-"RR"
RR_NP16$Vessel<-"SPOT3_v4"
RR_NP16$Drug<-"nimodipine"

RR_P1<-summarySE(RR_pyr3_1, measurevar="velocity", na.rm = TRUE)
RR_P1$Animal<-"RR"
RR_P1$Vessel<-"SPOT1_v1"
RR_P1$Drug<-"pyr3"

RR_P2<-summarySE(RR_pyr3_2, measurevar="velocity", na.rm = TRUE)
RR_P2$Animal<-"RR"
RR_P2$Vessel<-"SPOT1_v2"
RR_P2$Drug<-"pyr3"

RR_P3<-summarySE(RR_pyr3_3, measurevar="velocity", na.rm = TRUE)
RR_P3$Animal<-"RR"
RR_P3$Vessel<-"SPOT1_v3"
RR_P3$Drug<-"pyr3"

RR_P4<-summarySE(RR_pyr3_4, measurevar="velocity", na.rm = TRUE)
RR_P4$Animal<-"RR"
RR_P4$Vessel<-"SPOT1_v4"
RR_P4$Drug<-"pyr3"

RR_P5<-summarySE(RR_pyr3_5, measurevar="velocity", na.rm = TRUE)
RR_P5$Animal<-"RR"
RR_P5$Vessel<-"SPOT2_v1"
RR_P5$Drug<-"pyr3"

RR_P6<-summarySE(RR_pyr3_6, measurevar="velocity", na.rm = TRUE)
RR_P6$Animal<-"RR"
RR_P6$Vessel<-"SPOT2_v1"
RR_P6$Drug<-"pyr3"

RR_P7<-summarySE(RR_pyr3_7, measurevar="velocity", na.rm = TRUE)
RR_P7$Animal<-"RR"
RR_P7$Vessel<-"SPOT2_v2"
RR_P7$Drug<-"pyr3"

RR_P8<-summarySE(RR_pyr3_8, measurevar="velocity", na.rm = TRUE)
RR_P8$Animal<-"RR"
RR_P8$Vessel<-"SPOT3_v1"
RR_P8$Drug<-"pyr3"

RR_P9<-summarySE(RR_pyr3_9, measurevar="velocity", na.rm = TRUE)
RR_P9$Animal<-"RR"
RR_P9$Vessel<-"SPOT3_v2"
RR_P9$Drug<-"pyr3"

RR_P10<-summarySE(RR_pyr3_10, measurevar="velocity", na.rm = TRUE)
RR_P10$Animal<-"RR"
RR_P10$Vessel<-"SPOT3_v3"
RR_P10$Drug<-"pyr3"

MJ_BL1<-summarySE(MJ_baseline1, measurevar="velocity", na.rm = TRUE)
MJ_BL1$Animal<-"MJ"
MJ_BL1$Vessel<-"SPOT2_v1"
MJ_BL1$Drug<-"baseline"

MJ_BL2<-summarySE(MJ_baseline2, measurevar="velocity", na.rm = TRUE)
MJ_BL2$Animal<-"MJ"
MJ_BL2$Vessel<-"SPOT2_v2"
MJ_BL2$Drug<-"baseline"

MJ_BL3<-summarySE(MJ_baseline3, measurevar="velocity", na.rm = TRUE)
MJ_BL3$Animal<-"MJ"
MJ_BL3$Vessel<-"SPOT2_v2"
MJ_BL3$Drug<-"baseline"

MJ_BL4<-summarySE(MJ_baseline4, measurevar="velocity", na.rm = TRUE)
MJ_BL4$Animal<-"MJ"
MJ_BL4$Vessel<-"SPOT2_v3"
MJ_BL4$Drug<-"baseline"

MJ_BL5<-summarySE(MJ_baseline5, measurevar="velocity", na.rm = TRUE)
MJ_BL5$Animal<-"MJ"
MJ_BL5$Vessel<-"SPOT2_v3"
MJ_BL5$Drug<-"baseline"

MJ_BL6<-summarySE(MJ_baseline6, measurevar="velocity", na.rm = TRUE)
MJ_BL6$Animal<-"MJ"
MJ_BL6$Vessel<-"SPOT4_v1"
MJ_BL6$Drug<-"baseline"

MJ_BL7<-summarySE(MJ_baseline7, measurevar="velocity", na.rm = TRUE)
MJ_BL7$Animal<-"MJ"
MJ_BL7$Vessel<-"SPOT1_v1"
MJ_BL7$Drug<-"baseline"

MJ_BL8<-summarySE(MJ_baseline8, measurevar="velocity", na.rm = TRUE)
MJ_BL8$Animal<-"MJ"
MJ_BL8$Vessel<-"SPOT1_v2"
MJ_BL8$Drug<-"baseline"

MJ_BL9<-summarySE(MJ_baseline9, measurevar="velocity", na.rm = TRUE)
MJ_BL9$Animal<-"MJ"
MJ_BL9$Vessel<-"SPOT1_v3"
MJ_BL9$Drug<-"baseline"

MJ_BL10<-summarySE(MJ_baseline10, measurevar="velocity", na.rm = TRUE)
MJ_BL10$Animal<-"MJ"
MJ_BL10$Vessel<-"SPOT3_v1"
MJ_BL10$Drug<-"baseline"

MJ_BL11<-summarySE(MJ_baseline11, measurevar="velocity", na.rm = TRUE)
MJ_BL11$Animal<-"MJ"
MJ_BL11$Vessel<-"SPOT4_v1"
MJ_BL11$Drug<-"baseline"

MJ_BL12<-summarySE(MJ_baseline12, measurevar="velocity", na.rm = TRUE)
MJ_BL12$Animal<-"MJ"
MJ_BL12$Vessel<-"SPOT4_v2"
MJ_BL12$Drug<-"baseline"

MJ_BL13<-summarySE(MJ_baseline13, measurevar="velocity", na.rm = TRUE)
MJ_BL13$Animal<-"MJ"
MJ_BL13$Vessel<-"SPOT5_v1"
MJ_BL13$Drug<-"baseline"

MJ_BL14<-summarySE(MJ_baseline14, measurevar="velocity", na.rm = TRUE)
MJ_BL14$Animal<-"MJ"
MJ_BL14$Vessel<-"SPOT5_v2"
MJ_BL14$Drug<-"baseline"

MJ_BL15<-summarySE(MJ_baseline15, measurevar="velocity", na.rm = TRUE)
MJ_BL15$Animal<-"MJ"
MJ_BL15$Vessel<-"SPOT5_v3"
MJ_BL15$Drug<-"baseline"

MJ_BL16<-summarySE(MJ_baseline16, measurevar="velocity", na.rm = TRUE)
MJ_BL16$Animal<-"MJ"
MJ_BL16$Vessel<-"SPOT5_v4"
MJ_BL16$Drug<-"baseline"

MJ_BL17<-summarySE(MJ_baseline17, measurevar="velocity", na.rm = TRUE)
MJ_BL17$Animal<-"MJ"
MJ_BL17$Vessel<-"SPOT5_v5"
MJ_BL17$Drug<-"baseline"

MJ_BL18<-summarySE(MJ_baseline18, measurevar="velocity", na.rm = TRUE)
MJ_BL18$Animal<-"MJ"
MJ_BL18$Vessel<-"SPOT5_v6"
MJ_BL18$Drug<-"baseline"

MJ_BL19<-summarySE(MJ_baseline19, measurevar="velocity", na.rm = TRUE)
MJ_BL19$Animal<-"MJ"
MJ_BL19$Vessel<-"SPOT5_v6"
MJ_BL19$Drug<-"baseline"

MJ_BL20<-summarySE(MJ_baseline20, measurevar="velocity", na.rm = TRUE)
MJ_BL20$Animal<-"MJ"
MJ_BL20$Vessel<-"SPOT1_v1"
MJ_BL20$Drug<-"baseline"

MJ_BL21<-summarySE(MJ_baseline21, measurevar="velocity", na.rm = TRUE)
MJ_BL21$Animal<-"MJ"
MJ_BL21$Vessel<-"SPOT1_v2"
MJ_BL21$Drug<-"baseline"

MJ_BL22<-summarySE(MJ_baseline22, measurevar="velocity", na.rm = TRUE)
MJ_BL22$Animal<-"MJ"
MJ_BL22$Vessel<-"SPOT1_v3"
MJ_BL22$Drug<-"baseline"

MJ_BL23<-summarySE(MJ_baseline23, measurevar="velocity", na.rm = TRUE)
MJ_BL23$Animal<-"MJ"
MJ_BL23$Vessel<-"SPOT3_v1"
MJ_BL23$Drug<-"baseline"

MJ_BL24<-summarySE(MJ_baseline24, measurevar="velocity", na.rm = TRUE)
MJ_BL24$Animal<-"MJ"
MJ_BL24$Vessel<-"SPOT3_v1"
MJ_BL24$Drug<-"baseline"

MJ_BL25<-summarySE(MJ_baseline25, measurevar="velocity", na.rm = TRUE)
MJ_BL25$Animal<-"MJ"
MJ_BL25$Vessel<-"SPOT6_v1"
MJ_BL25$Drug<-"baseline"

MJ_NP1<-summarySE(MJ_nimodipine1, measurevar="velocity", na.rm = TRUE)
MJ_NP1$Animal<-"MJ"
MJ_NP1$Vessel<-"SPOT3_v1"
MJ_NP1$Drug<-"nimodipine"

MJ_NP2<-summarySE(MJ_nimodipine2, measurevar="velocity", na.rm = TRUE)
MJ_NP2$Animal<-"MJ"
MJ_NP2$Vessel<-"SPOT3_v1"
MJ_NP2$Drug<-"nimodipine"

MJ_NP3<-summarySE(MJ_nimodipine3, measurevar="velocity", na.rm = TRUE)
MJ_NP3$Animal<-"MJ"
MJ_NP3$Vessel<-"SPOT4_v1"
MJ_NP3$Drug<-"nimodipine"

MJ_NP4<-summarySE(MJ_nimodipine4, measurevar="velocity", na.rm = TRUE)
MJ_NP4$Animal<-"MJ"
MJ_NP4$Vessel<-"SPOT4_v2"
MJ_NP4$Drug<-"nimodipine"

MJ_NP5<-summarySE(MJ_nimodipine5, measurevar="velocity", na.rm = TRUE)
MJ_NP5$Animal<-"MJ"
MJ_NP5$Vessel<-"SPOT4_v2"
MJ_NP5$Drug<-"nimodipine"

MJ_NP6<-summarySE(MJ_nimodipine6, measurevar="velocity", na.rm = TRUE)
MJ_NP6$Animal<-"MJ"
MJ_NP6$Vessel<-"SPOT5_v1"
MJ_NP6$Drug<-"nimodipine"

MJ_NP7<-summarySE(MJ_nimodipine7, measurevar="velocity", na.rm = TRUE)
MJ_NP7$Animal<-"MJ"
MJ_NP7$Vessel<-"SPOT5_v2"
MJ_NP7$Drug<-"nimodipine"

MJ_NP8<-summarySE(MJ_nimodipine8, measurevar="velocity", na.rm = TRUE)
MJ_NP8$Animal<-"MJ"
MJ_NP8$Vessel<-"SPOT5_v2"
MJ_NP8$Drug<-"nimodipine"

MJ_NP9<-summarySE(MJ_nimodipine9, measurevar="velocity", na.rm = TRUE)
MJ_NP9$Animal<-"MJ"
MJ_NP9$Vessel<-"SPOT5_v3"
MJ_NP9$Drug<-"nimodipine"

MJ_NP10<-summarySE(MJ_nimodipine10, measurevar="velocity", na.rm = TRUE)
MJ_NP10$Animal<-"MJ"
MJ_NP10$Vessel<-"SPOT5_v3"
MJ_NP10$Drug<-"nimodipine"

MJ_NP11<-summarySE(MJ_nimodipine11, measurevar="velocity", na.rm = TRUE)
MJ_NP11$Animal<-"MJ"
MJ_NP11$Vessel<-"SPOT5_v4"
MJ_NP11$Drug<-"nimodipine"

MJ_NP12<-summarySE(MJ_nimodipine12, measurevar="velocity", na.rm = TRUE)
MJ_NP12$Animal<-"MJ"
MJ_NP12$Vessel<-"SPOT5_v5"
MJ_NP12$Drug<-"nimodipine"

MJ_NP13<-summarySE(MJ_nimodipine13, measurevar="velocity", na.rm = TRUE)
MJ_NP13$Animal<-"MJ"
MJ_NP13$Vessel<-"SPOT5_v6"
MJ_NP13$Drug<-"nimodipine"

MJ_NP14<-summarySE(MJ_nimodipine14, measurevar="velocity", na.rm = TRUE)
MJ_NP14$Animal<-"MJ"
MJ_NP14$Vessel<-"SPOT5_v6"
MJ_NP14$Drug<-"nimodipine"

MJ_NP15<-summarySE(MJ_nimodipine15, measurevar="velocity", na.rm = TRUE)
MJ_NP15$Animal<-"MJ"
MJ_NP15$Vessel<-"SPOT6_v1"
MJ_NP15$Drug<-"nimodipine"

MJ_P1<-summarySE(MJ_pyr3_1, measurevar="velocity", na.rm = TRUE)
MJ_P1$Animal<-"MJ"
MJ_P1$Vessel<-"SPOT1_v1"
MJ_P1$Drug<-"pyr3"

MJ_P2<-summarySE(MJ_pyr3_2, measurevar="velocity", na.rm = TRUE)
MJ_P2$Animal<-"MJ"
MJ_P2$Vessel<-"SPOT1_v1"
MJ_P2$Drug<-"pyr3"

MJ_P3<-summarySE(MJ_pyr3_3, measurevar="velocity", na.rm = TRUE)
MJ_P3$Animal<-"MJ"
MJ_P3$Vessel<-"SPOT1_v2"
MJ_P3$Drug<-"pyr3"

MJ_P4<-summarySE(MJ_pyr3_4, measurevar="velocity", na.rm = TRUE)
MJ_P4$Animal<-"MJ"
MJ_P4$Vessel<-"SPOT1_v2"
MJ_P4$Drug<-"pyr3"

MJ_P5<-summarySE(MJ_pyr3_5, measurevar="velocity", na.rm = TRUE)
MJ_P5$Animal<-"MJ"
MJ_P5$Vessel<-"SPOT1_v3"
MJ_P5$Drug<-"pyr3"

MJ_P6<-summarySE(MJ_pyr3_6, measurevar="velocity", na.rm = TRUE)
MJ_P6$Animal<-"MJ"
MJ_P6$Vessel<-"SPOT3_v1"
MJ_P6$Drug<-"pyr3"

MJ_P7<-summarySE(MJ_pyr3_7, measurevar="velocity", na.rm = TRUE)
MJ_P7$Animal<-"MJ"
MJ_P7$Vessel<-"SPOT4_v1"
MJ_P7$Drug<-"pyr3"

MJ_P8<-summarySE(MJ_pyr3_8, measurevar="velocity", na.rm = TRUE)
MJ_P8$Animal<-"MJ"
MJ_P8$Vessel<-"SPOT4_v1"
MJ_P8$Drug<-"pyr3"

MJ_P9<-summarySE(MJ_pyr3_9, measurevar="velocity", na.rm = TRUE)
MJ_P9$Animal<-"MJ"
MJ_P9$Vessel<-"SPOT4_v2"
MJ_P9$Drug<-"pyr3"

MJ_P10<-summarySE(MJ_pyr3_10, measurevar="velocity", na.rm = TRUE)
MJ_P10$Animal<-"MJ"
MJ_P10$Vessel<-"SPOT5_v1"
MJ_P10$Drug<-"pyr3"

MJ_P11<-summarySE(MJ_pyr3_11, measurevar="velocity", na.rm = TRUE)
MJ_P11$Animal<-"MJ"
MJ_P11$Vessel<-"SPOT5_v2"
MJ_P11$Drug<-"pyr3"

MJ_P12<-summarySE(MJ_pyr3_12, measurevar="velocity", na.rm = TRUE)
MJ_P12$Animal<-"MJ"
MJ_P12$Vessel<-"SPOT5_v3"
MJ_P12$Drug<-"pyr3"

MJ_P13<-summarySE(MJ_pyr3_13, measurevar="velocity", na.rm = TRUE)
MJ_P13$Animal<-"MJ"
MJ_P13$Vessel<-"SPOT5_v4"
MJ_P13$Drug<-"pyr3"

MJ_P14<-summarySE(MJ_pyr3_14, measurevar="velocity", na.rm = TRUE)
MJ_P14$Animal<-"MJ"
MJ_P14$Vessel<-"SPOT5_v5"
MJ_P14$Drug<-"pyr3"

MJ_P15<-summarySE(MJ_pyr3_15, measurevar="velocity", na.rm = TRUE)
MJ_P15$Animal<-"MJ"
MJ_P15$Vessel<-"SPOT5_v6"
MJ_P15$Drug<-"pyr3"

MJ_P16<-summarySE(MJ_pyr3_16, measurevar="velocity", na.rm = TRUE)
MJ_P16$Animal<-"MJ"
MJ_P16$Vessel<-"SPOT6_v1"
MJ_P16$Drug<-"pyr3"

#combine data together for each treatment
baseline<-rbind(RR_BL1,RR_BL2,RR_BL3,RR_BL4,RR_BL5,RR_BL6,MJ_BL1,MJ_BL2,MJ_BL3,MJ_BL4,MJ_BL5,MJ_BL6,MJ_BL7,MJ_BL8,MJ_BL9,MJ_BL10,MJ_BL11,MJ_BL12,MJ_BL13,MJ_BL14,MJ_BL15,MJ_BL16,MJ_BL17,MJ_BL18,MJ_BL19,MJ_BL20,MJ_BL21,MJ_BL22,MJ_BL23,MJ_BL24,MJ_BL25) #L_baseline1,L_baseline2,L_baseline3,L_baseline4,L_baseline5,L_baseline6,L_baseline7,L_baseline8,L_baseline9,L_baseline10,L_baseline11,JJ_baseline1,JJ_baseline2,JJ_baseline3,JJ_baseline4,JJ_baseline5,JJ_baseline6,OM_baseline1,OM_baseline2,OM_baseline3,OM_baseline4,OM_baseline5,OM_baseline6,OM_baseline7,NL_baseline1,NL_baseline2,NL_baseline3,NL_baseline4)
nimodipine<-rbind(RR_NP1,RR_NP2,RR_NP3,RR_NP4,RR_NP5,RR_NP6,RR_NP7,RR_NP8,RR_NP9,RR_NP10,RR_NP11,RR_NP12,RR_NP13,RR_NP14,RR_NP15,RR_NP16,MJ_NP1,MJ_NP2,MJ_NP3,MJ_NP4,MJ_NP5,MJ_NP6,MJ_NP7,MJ_NP8,MJ_NP9,MJ_NP10,MJ_NP11,MJ_NP12,MJ_NP13,MJ_NP14,MJ_NP15)#L_nimodipine1
pyr3<-rbind(RR_P1,RR_P2,RR_P3,RR_P4,RR_P5,RR_P6,RR_P7,RR_P8,RR_P9,RR_P10,MJ_P1,MJ_P2,MJ_P3,MJ_P4,MJ_P5,MJ_P6,MJ_P7,MJ_P8,MJ_P9,MJ_P10,MJ_P11,MJ_P12,MJ_P13,MJ_P14,MJ_P15,MJ_P16)

AllData1<-rbind(baseline,nimodipine,pyr3)
AllData1$Drug<-as.factor(AllData1$Drug)

#pull out the spot name
pos= regexpr('SPOT', AllData1$Vessel)
AllData1$SpotName<-substr(AllData1$Vessel,pos, pos+4)

AllData1$velocity<-abs(AllData1$velocity)

# load files
#baseline
RR_baseline1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_08/2019_07_08_SPOT1_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
RR_baseline2 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_08/2019_07_08_SPOT1_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
RR_baseline3 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_08/2019_07_08_SPOT1_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
RR_baseline4 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_11/2019_07_11_SPOT3_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
RR_baseline5 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_11/2019_07_11_SPOT3_v4_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
RR_baseline6 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_22/2019_07_22_SPOT4_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
RR_baseline7 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_22/2019_07_22_SPOT4_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
RR_baseline8 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_22/2019_07_22_SPOT4_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
RR_baseline9 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_22/2019_07_22_SPOT4_v4_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")

MJ_baseline1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_06_11_KX/2019_06_11_SPOT2_v1_diam_LineScan_KX_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline2 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_06_11_KX/2019_06_11_SPOT2_v2_diam_LineScan_KX_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline3 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_06_11_KX/2019_06_11_SPOT2_v3_diam_LineScan_KX_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline4 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_06_11_KX/2019_06_11_SPOT4_v1_diam_LineScan_KX_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline5 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_09/2019_07_09_SPOT1_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline6 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_09/2019_07_09_SPOT1_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline7 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_09/2019_07_09_SPOT1_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline8 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_11/2019_07_11_SPOT3_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline9 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_11/2019_07_11_SPOT6_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline10 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_22/2019_07_22_SPOT4_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline11 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_22/2019_07_22_SPOT4_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline12 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_22/2019_07_22_SPOT5_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline13 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_22/2019_07_22_SPOT5_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline14 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_22/2019_07_22_SPOT5_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline15 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_22/2019_07_22_SPOT5_v4_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline16 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_22/2019_07_22_SPOT5_v5_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline17 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_22/2019_07_22_SPOT5_v6_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline18 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_19/2019_08_19_SPOT1_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline19 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_19/2019_08_19_SPOT1_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline20 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_19/2019_08_19_SPOT1_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline21 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_19/2019_08_19_SPOT3_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline22 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_19/2019_08_19_SPOT6_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")

#nimodipine
RR_nimodipine1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_16 Nimodipine/2019_07_16_SPOT3_v1_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
RR_nimodipine2 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_16 Nimodipine/2019_07_16_SPOT3_v2_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
RR_nimodipine3 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_16 Nimodipine/2019_07_16_SPOT3_v3_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
RR_nimodipine4 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_16 Nimodipine/2019_07_16_SPOT1_v1_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
RR_nimodipine5 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_16 Nimodipine/2019_07_16_SPOT1_v2_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
RR_nimodipine6 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_16 Nimodipine/2019_07_16_SPOT1_v3_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
RR_nimodipine7 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_16 Nimodipine/2019_07_16_SPOT1_v4_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
RR_nimodipine8 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT4_v1_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
RR_nimodipine9 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT4_v2_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
RR_nimodipine10 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT4_v3_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
RR_nimodipine11 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT4_v4_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
RR_nimodipine12 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT5_v1_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
RR_nimodipine13 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT5_v2_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
RR_nimodipine14 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT5_v3_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
RR_nimodipine15 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT5_v4_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")

MJ_nimodipine1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_16 Nimodipine/2019_07_16_SPOT1_v1_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_nimodipine2 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_16 Nimodipine/2019_07_16_SPOT3_v1_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_nimodipine3 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_16 Nimodipine/2019_07_16_SPOT3_v2_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_nimodipine4 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_16 Nimodipine/2019_07_16_SPOT3_v3_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_nimodipine5 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT3_v1_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_nimodipine6 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT3_v2_diam_145_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_nimodipine7 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT4_v1_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_nimodipine8 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT4_v2_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_nimodipine9 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT5_v1_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_nimodipine10 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT5_v2_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_nimodipine11 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT5_v3_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_nimodipine12 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT5_v4_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_nimodipine13 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT5_v5_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_nimodipine14 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_02 Nimodipine/2019_08_02_SPOT5_v6_diam_LineScan_nimodipine_diameterFWHM.csv", header=TRUE, sep = ",")

#pyr3
RR_pyr3_1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_19 PYR3/2019_07_19_SPOT1_v1_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
RR_pyr3_2 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_19 PYR3/2019_07_19_SPOT1_v2_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
RR_pyr3_3 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_19 PYR3/2019_07_19_SPOT1_v3_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
RR_pyr3_4 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_19 PYR3/2019_07_19_SPOT1_v4_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
RR_pyr3_5 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_19 PYR3/2019_07_19_SPOT2_v1_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
RR_pyr3_6 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_19 PYR3/2019_07_19_SPOT2_v2_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
RR_pyr3_7 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_19 PYR3/2019_07_19_SPOT3_v1_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
RR_pyr3_8 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_19 PYR3/2019_07_19_SPOT3_v2_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
RR_pyr3_9 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/RR/Linescan - Diameter/2019_07_19 PYR3/2019_07_19_SPOT3_v3_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")

MJ_pyr3_1 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_19 PYR3/2019_07_19_SPOT1_v1_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_pyr3_2 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_19 PYR3/2019_07_19_SPOT1_v2_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_pyr3_3 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_19 PYR3/2019_07_19_SPOT1_v3_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_pyr3_4 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_19 PYR3/2019_07_19_SPOT3_v1_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_pyr3_5 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_07_19 PYR3/2019_07_19_SPOT6_v1_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_pyr3_6 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_14 PYR3/2019_08_14_SPOT3_v1_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_pyr3_7 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_14 PYR3/2019_08_14_SPOT3_v3_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_pyr3_8 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_14 PYR3/2019_08_14_SPOT4_v1_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_pyr3_9 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_14 PYR3/2019_08_14_SPOT4_v2_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_pyr3_10 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_14 PYR3/2019_08_14_SPOT5_v1_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_pyr3_11 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_14 PYR3/2019_08_14_SPOT5_v2_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_pyr3_12 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_14 PYR3/2019_08_14_SPOT5_v3_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_pyr3_13 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_14 PYR3/2019_08_14_SPOT5_v4_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_pyr3_14 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_14 PYR3/2019_08_14_SPOT5_v5_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_pyr3_15 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_14 PYR3/2019_08_14_SPOT5_v6_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_pyr3_16 <- read.table("C:/Users/rodrigmc-INS/Desktop/Analysis In Vivo/MJ/Line Scan - Diameter/2019_08_14 PYR3/2019_08_14_SPOT6_v1_diam_LineScan_pyr3_diameterFWHM.csv", header=TRUE, sep = ",")

#####
RR_BL1<-summarySE(RR_baseline1, measurevar="diameter", na.rm = TRUE)
RR_BL1$Animal<-"RR"
RR_BL1$Vessel<-"SPOT1_v1"
RR_BL1$Drug<-"baseline"


RR_BL2<-summarySE(RR_baseline2, measurevar="diameter", na.rm = TRUE)
RR_BL2$Animal<-"RR"
RR_BL2$Vessel<-"SPOT1_v2"
RR_BL2$Drug<-"baseline"

RR_BL3<-summarySE(RR_baseline3, measurevar="diameter", na.rm = TRUE)
RR_BL3$Animal<-"RR"
RR_BL3$Vessel<-"SPOT1_v3"
RR_BL3$Drug<-"baseline"

RR_BL4<-summarySE(RR_baseline4, measurevar="diameter", na.rm = TRUE)
RR_BL4$Animal<-"RR"
RR_BL4$Vessel<-"SPOT3_v1"
RR_BL4$Drug<-"baseline"

RR_BL5<-summarySE(RR_baseline5, measurevar="diameter", na.rm = TRUE)
RR_BL5$Animal<-"RR"
RR_BL5$Vessel<-"SPOT3_v4"
RR_BL5$Drug<-"baseline"

RR_BL6<-summarySE(RR_baseline6, measurevar="diameter", na.rm = TRUE)
RR_BL6$Animal<-"RR"
RR_BL6$Vessel<-"SPOT4_v1"
RR_BL6$Drug<-"baseline"

RR_BL7<-summarySE(RR_baseline7, measurevar="diameter", na.rm = TRUE)
RR_BL7$Animal<-"RR"
RR_BL7$Vessel<-"SPOT4_v2"
RR_BL7$Drug<-"baseline"

RR_BL8<-summarySE(RR_baseline8, measurevar="diameter", na.rm = TRUE)
RR_BL8$Animal<-"RR"
RR_BL8$Vessel<-"SPOT4_v3"
RR_BL8$Drug<-"baseline"

RR_BL9<-summarySE(RR_baseline9, measurevar="diameter", na.rm = TRUE)
RR_BL9$Animal<-"RR"
RR_BL9$Vessel<-"SPOT4_v3"
RR_BL9$Drug<-"baseline"

RR_NP1<-summarySE(RR_nimodipine1, measurevar="diameter", na.rm = TRUE)
RR_NP1$Animal<-"RR"
RR_NP1$Vessel<-"SPOT3_v1"
RR_NP1$Drug<-"nimodipine"

RR_NP2<-summarySE(RR_nimodipine2, measurevar="diameter", na.rm = TRUE)
RR_NP2$Animal<-"RR"
RR_NP2$Vessel<-"SPOT3_v2"
RR_NP2$Drug<-"nimodipine"

RR_NP3<-summarySE(RR_nimodipine3, measurevar="diameter", na.rm = TRUE)
RR_NP3$Animal<-"RR"
RR_NP3$Vessel<-"SPOT3_v3"
RR_NP3$Drug<-"nimodipine"

RR_NP4<-summarySE(RR_nimodipine4, measurevar="diameter", na.rm = TRUE)
RR_NP4$Animal<-"RR"
RR_NP4$Vessel<-"SPOT1_v1"
RR_NP4$Drug<-"nimodipine"

RR_NP5<-summarySE(RR_nimodipine5, measurevar="diameter", na.rm = TRUE)
RR_NP5$Animal<-"RR"
RR_NP5$Vessel<-"SPOT1_v2"
RR_NP5$Drug<-"nimodipine"

RR_NP6<-summarySE(RR_nimodipine6, measurevar="diameter", na.rm = TRUE)
RR_NP6$Animal<-"RR"
RR_NP6$Vessel<-"SPOT1_v3"
RR_NP6$Drug<-"nimodipine"

RR_NP7<-summarySE(RR_nimodipine7, measurevar="diameter", na.rm = TRUE)
RR_NP7$Animal<-"RR"
RR_NP7$Vessel<-"SPOT1_v4"
RR_NP7$Drug<-"nimodipine"

RR_NP8<-summarySE(RR_nimodipine8, measurevar="diameter", na.rm = TRUE)
RR_NP8$Animal<-"RR"
RR_NP8$Vessel<-"SPOT4_v1"
RR_NP8$Drug<-"nimodipine"

RR_NP9<-summarySE(RR_nimodipine9, measurevar="diameter", na.rm = TRUE)
RR_NP9$Animal<-"RR"
RR_NP9$Vessel<-"SPOT4_v2"
RR_NP9$Drug<-"nimodipine"

RR_NP10<-summarySE(RR_nimodipine10, measurevar="diameter", na.rm = TRUE)
RR_NP10$Animal<-"RR"
RR_NP10$Vessel<-"SPOT4_v3"
RR_NP10$Drug<-"nimodipine"

RR_NP11<-summarySE(RR_nimodipine11, measurevar="diameter", na.rm = TRUE)
RR_NP11$Animal<-"RR"
RR_NP11$Vessel<-"SPOT4_v4"
RR_NP11$Drug<-"nimodipine"

RR_NP12<-summarySE(RR_nimodipine12, measurevar="diameter", na.rm = TRUE)
RR_NP12$Animal<-"RR"
RR_NP12$Vessel<-"SPOT5_v1"
RR_NP12$Drug<-"nimodipine"

RR_NP13<-summarySE(RR_nimodipine13, measurevar="diameter", na.rm = TRUE)
RR_NP13$Animal<-"RR"
RR_NP13$Vessel<-"SPOT5_v2"
RR_NP13$Drug<-"nimodipine"

RR_NP14<-summarySE(RR_nimodipine14, measurevar="diameter", na.rm = TRUE)
RR_NP14$Animal<-"RR"
RR_NP14$Vessel<-"SPOT5_v3"
RR_NP14$Drug<-"nimodipine"

RR_NP15<-summarySE(RR_nimodipine15, measurevar="diameter", na.rm = TRUE)
RR_NP15$Animal<-"RR"
RR_NP15$Vessel<-"SPOT5_v4"
RR_NP15$Drug<-"nimodipine"

RR_P1<-summarySE(RR_pyr3_1, measurevar="diameter", na.rm = TRUE)
RR_P1$Animal<-"RR"
RR_P1$Vessel<-"SPOT1_v1"
RR_P1$Drug<-"pyr3"

RR_P2<-summarySE(RR_pyr3_2, measurevar="diameter", na.rm = TRUE)
RR_P2$Animal<-"RR"
RR_P2$Vessel<-"SPOT1_v2"
RR_P2$Drug<-"pyr3"

RR_P3<-summarySE(RR_pyr3_3, measurevar="diameter", na.rm = TRUE)
RR_P3$Animal<-"RR"
RR_P3$Vessel<-"SPOT1_v3"
RR_P3$Drug<-"pyr3"

RR_P4<-summarySE(RR_pyr3_4, measurevar="diameter", na.rm = TRUE)
RR_P4$Animal<-"RR"
RR_P4$Vessel<-"SPOT1_v4"
RR_P4$Drug<-"pyr3"

RR_P5<-summarySE(RR_pyr3_5, measurevar="diameter", na.rm = TRUE)
RR_P5$Animal<-"RR"
RR_P5$Vessel<-"SPOT2_v1"
RR_P5$Drug<-"pyr3"

RR_P6<-summarySE(RR_pyr3_6, measurevar="diameter", na.rm = TRUE)
RR_P6$Animal<-"RR"
RR_P6$Vessel<-"SPOT2_v2"
RR_P6$Drug<-"pyr3"

RR_P7<-summarySE(RR_pyr3_7, measurevar="diameter", na.rm = TRUE)
RR_P7$Animal<-"RR"
RR_P7$Vessel<-"SPOT3_v1"
RR_P7$Drug<-"pyr3"

RR_P8<-summarySE(RR_pyr3_8, measurevar="diameter", na.rm = TRUE)
RR_P8$Animal<-"RR"
RR_P8$Vessel<-"SPOT3_v2"
RR_P8$Drug<-"pyr3"

RR_P9<-summarySE(RR_pyr3_9, measurevar="diameter", na.rm = TRUE)
RR_P9$Animal<-"RR"
RR_P9$Vessel<-"SPOT3_v3"
RR_P9$Drug<-"pyr3"

MJ_BL1<-summarySE(MJ_baseline1, measurevar="diameter", na.rm = TRUE)
MJ_BL1$Animal<-"MJ"
MJ_BL1$Vessel<-"SPOT2_v1"
MJ_BL1$Drug<-"baseline"

MJ_BL2<-summarySE(MJ_baseline2, measurevar="diameter", na.rm = TRUE)
MJ_BL2$Animal<-"MJ"
MJ_BL2$Vessel<-"SPOT2_v2"
MJ_BL2$Drug<-"baseline"

MJ_BL3<-summarySE(MJ_baseline3, measurevar="diameter", na.rm = TRUE)
MJ_BL3$Animal<-"MJ"
MJ_BL3$Vessel<-"SPOT2_v3"
MJ_BL3$Drug<-"baseline"

MJ_BL4<-summarySE(MJ_baseline4, measurevar="diameter", na.rm = TRUE)
MJ_BL4$Animal<-"MJ"
MJ_BL4$Vessel<-"SPOT4_v1"
MJ_BL4$Drug<-"baseline"

MJ_BL5<-summarySE(MJ_baseline5, measurevar="diameter", na.rm = TRUE)
MJ_BL5$Animal<-"MJ"
MJ_BL5$Vessel<-"SPOT1_v1"
MJ_BL5$Drug<-"baseline"

MJ_BL6<-summarySE(MJ_baseline6, measurevar="diameter", na.rm = TRUE)
MJ_BL6$Animal<-"MJ"
MJ_BL6$Vessel<-"SPOT1_v2"
MJ_BL6$Drug<-"baseline"

MJ_BL7<-summarySE(MJ_baseline7, measurevar="diameter", na.rm = TRUE)
MJ_BL7$Animal<-"MJ"
MJ_BL7$Vessel<-"SPOT1_v3"
MJ_BL7$Drug<-"baseline"

MJ_BL8<-summarySE(MJ_baseline8, measurevar="diameter", na.rm = TRUE)
MJ_BL8$Animal<-"MJ"
MJ_BL8$Vessel<-"SPOT3_v1"
MJ_BL8$Drug<-"baseline"

MJ_BL9<-summarySE(MJ_baseline9, measurevar="diameter", na.rm = TRUE)
MJ_BL9$Animal<-"MJ"
MJ_BL9$Vessel<-"SPOT6_v1"
MJ_BL9$Drug<-"baseline"

MJ_BL10<-summarySE(MJ_baseline10, measurevar="diameter", na.rm = TRUE)
MJ_BL10$Animal<-"MJ"
MJ_BL10$Vessel<-"SPOT4_v1"
MJ_BL10$Drug<-"baseline"

MJ_BL11<-summarySE(MJ_baseline11, measurevar="diameter", na.rm = TRUE)
MJ_BL11$Animal<-"MJ"
MJ_BL11$Vessel<-"SPOT4_v2"
MJ_BL11$Drug<-"baseline"

MJ_BL12<-summarySE(MJ_baseline12, measurevar="diameter", na.rm = TRUE)
MJ_BL12$Animal<-"MJ"
MJ_BL12$Vessel<-"SPOT5_v1"
MJ_BL12$Drug<-"baseline"

MJ_BL13<-summarySE(MJ_baseline13, measurevar="diameter", na.rm = TRUE)
MJ_BL13$Animal<-"MJ"
MJ_BL13$Vessel<-"SPOT5_v2"
MJ_BL13$Drug<-"baseline"

MJ_BL14<-summarySE(MJ_baseline14, measurevar="diameter", na.rm = TRUE)
MJ_BL14$Animal<-"MJ"
MJ_BL14$Vessel<-"SPOT5_v3"
MJ_BL14$Drug<-"baseline"

MJ_BL15<-summarySE(MJ_baseline15, measurevar="diameter", na.rm = TRUE)
MJ_BL15$Animal<-"MJ"
MJ_BL15$Vessel<-"SPOT5_v4"
MJ_BL15$Drug<-"baseline"

MJ_BL16<-summarySE(MJ_baseline16, measurevar="diameter", na.rm = TRUE)
MJ_BL16$Animal<-"MJ"
MJ_BL16$Vessel<-"SPOT5_v5"
MJ_BL16$Drug<-"baseline"

MJ_BL17<-summarySE(MJ_baseline17, measurevar="diameter", na.rm = TRUE)
MJ_BL17$Animal<-"MJ"
MJ_BL17$Vessel<-"SPOT5_v6"
MJ_BL17$Drug<-"baseline"

MJ_BL18<-summarySE(MJ_baseline18, measurevar="diameter", na.rm = TRUE)
MJ_BL18$Animal<-"MJ"
MJ_BL18$Vessel<-"SPOT1_v1"
MJ_BL18$Drug<-"baseline"

MJ_BL19<-summarySE(MJ_baseline19, measurevar="diameter", na.rm = TRUE)
MJ_BL19$Animal<-"MJ"
MJ_BL19$Vessel<-"SPOT1_v2"
MJ_BL19$Drug<-"baseline"

MJ_BL20<-summarySE(MJ_baseline20, measurevar="diameter", na.rm = TRUE)
MJ_BL20$Animal<-"MJ"
MJ_BL20$Vessel<-"SPOT1_v3"
MJ_BL20$Drug<-"baseline"

MJ_BL21<-summarySE(MJ_baseline21, measurevar="diameter", na.rm = TRUE)
MJ_BL21$Animal<-"MJ"
MJ_BL21$Vessel<-"SPOT3_v1"
MJ_BL21$Drug<-"baseline"

MJ_BL22<-summarySE(MJ_baseline22, measurevar="diameter", na.rm = TRUE)
MJ_BL22$Animal<-"MJ"
MJ_BL22$Vessel<-"SPOT6_v1"
MJ_BL22$Drug<-"baseline"

MJ_NP1<-summarySE(MJ_nimodipine1, measurevar="diameter", na.rm = TRUE)
MJ_NP1$Animal<-"MJ"
MJ_NP1$Vessel<-"SPOT1_v1"
MJ_NP1$Drug<-"nimodipine"

MJ_NP2<-summarySE(MJ_nimodipine2, measurevar="diameter", na.rm = TRUE)
MJ_NP2$Animal<-"MJ"
MJ_NP2$Vessel<-"SPOT3_v1"
MJ_NP2$Drug<-"nimodipine"

MJ_NP3<-summarySE(MJ_nimodipine3, measurevar="diameter", na.rm = TRUE)
MJ_NP3$Animal<-"MJ"
MJ_NP3$Vessel<-"SPOT3_v2"
MJ_NP3$Drug<-"nimodipine"

MJ_NP4<-summarySE(MJ_nimodipine4, measurevar="diameter", na.rm = TRUE)
MJ_NP4$Animal<-"MJ"
MJ_NP4$Vessel<-"SPOT3_v3"
MJ_NP4$Drug<-"nimodipine"

MJ_NP5<-summarySE(MJ_nimodipine5, measurevar="diameter", na.rm = TRUE)
MJ_NP5$Animal<-"MJ"
MJ_NP5$Vessel<-"SPOT3_v1"
MJ_NP5$Drug<-"nimodipine"

MJ_NP6<-summarySE(MJ_nimodipine6, measurevar="diameter", na.rm = TRUE)
MJ_NP6$Animal<-"MJ"
MJ_NP6$Vessel<-"SPOT3_v2"
MJ_NP6$Drug<-"nimodipine"

MJ_NP7<-summarySE(MJ_nimodipine7, measurevar="diameter", na.rm = TRUE)
MJ_NP7$Animal<-"MJ"
MJ_NP7$Vessel<-"SPOT4_v1"
MJ_NP7$Drug<-"nimodipine"

MJ_NP8<-summarySE(MJ_nimodipine8, measurevar="diameter", na.rm = TRUE)
MJ_NP8$Animal<-"MJ"
MJ_NP8$Vessel<-"SPOT4_v2"
MJ_NP8$Drug<-"nimodipine"

MJ_NP9<-summarySE(MJ_nimodipine9, measurevar="diameter", na.rm = TRUE)
MJ_NP9$Animal<-"MJ"
MJ_NP9$Vessel<-"SPOT5_v1"
MJ_NP9$Drug<-"nimodipine"

MJ_NP10<-summarySE(MJ_nimodipine10, measurevar="diameter", na.rm = TRUE)
MJ_NP10$Animal<-"MJ"
MJ_NP10$Vessel<-"SPOT5_v2"
MJ_NP10$Drug<-"nimodipine"

MJ_NP11<-summarySE(MJ_nimodipine11, measurevar="diameter", na.rm = TRUE)
MJ_NP11$Animal<-"MJ"
MJ_NP11$Vessel<-"SPOT5_v3"
MJ_NP11$Drug<-"nimodipine"

MJ_NP12<-summarySE(MJ_nimodipine12, measurevar="diameter", na.rm = TRUE)
MJ_NP12$Animal<-"MJ"
MJ_NP12$Vessel<-"SPOT5_v4"
MJ_NP12$Drug<-"nimodipine"

MJ_NP13<-summarySE(MJ_nimodipine13, measurevar="diameter", na.rm = TRUE)
MJ_NP13$Animal<-"MJ"
MJ_NP13$Vessel<-"SPOT5_v5"
MJ_NP13$Drug<-"nimodipine"

MJ_NP14<-summarySE(MJ_nimodipine14, measurevar="diameter", na.rm = TRUE)
MJ_NP14$Animal<-"MJ"
MJ_NP14$Vessel<-"SPOT5_v6"
MJ_NP14$Drug<-"nimodipine"

MJ_P1<-summarySE(MJ_pyr3_1, measurevar="diameter", na.rm = TRUE)
MJ_P1$Animal<-"MJ"
MJ_P1$Vessel<-"SPOT1_v1"
MJ_P1$Drug<-"pyr3"

MJ_P2<-summarySE(MJ_pyr3_2, measurevar="diameter", na.rm = TRUE)
MJ_P2$Animal<-"MJ"
MJ_P2$Vessel<-"SPOT1_v2"
MJ_P2$Drug<-"pyr3"

MJ_P3<-summarySE(MJ_pyr3_3, measurevar="diameter", na.rm = TRUE)
MJ_P3$Animal<-"MJ"
MJ_P3$Vessel<-"SPOT1_v3"
MJ_P3$Drug<-"pyr3"

MJ_P4<-summarySE(MJ_pyr3_4, measurevar="diameter", na.rm = TRUE)
MJ_P4$Animal<-"MJ"
MJ_P4$Vessel<-"SPOT3_v1"
MJ_P4$Drug<-"pyr3"

MJ_P5<-summarySE(MJ_pyr3_5, measurevar="diameter", na.rm = TRUE)
MJ_P5$Animal<-"MJ"
MJ_P5$Vessel<-"SPOT6_v1"
MJ_P5$Drug<-"pyr3"

MJ_P6<-summarySE(MJ_pyr3_6, measurevar="diameter", na.rm = TRUE)
MJ_P6$Animal<-"MJ"
MJ_P6$Vessel<-"SPOT3_v1"
MJ_P6$Drug<-"pyr3"

MJ_P7<-summarySE(MJ_pyr3_7, measurevar="diameter", na.rm = TRUE)
MJ_P7$Animal<-"MJ"
MJ_P7$Vessel<-"SPOT3_v3"
MJ_P7$Drug<-"pyr3"

MJ_P8<-summarySE(MJ_pyr3_8, measurevar="diameter", na.rm = TRUE)
MJ_P8$Animal<-"MJ"
MJ_P8$Vessel<-"SPOT4_v1"
MJ_P8$Drug<-"pyr3"

MJ_P9<-summarySE(MJ_pyr3_9, measurevar="diameter", na.rm = TRUE)
MJ_P9$Animal<-"MJ"
MJ_P9$Vessel<-"SPOT4_v2"
MJ_P9$Drug<-"pyr3"

MJ_P10<-summarySE(MJ_pyr3_10, measurevar="diameter", na.rm = TRUE)
MJ_P10$Animal<-"MJ"
MJ_P10$Vessel<-"SPOT5_v1"
MJ_P10$Drug<-"pyr3"

MJ_P11<-summarySE(MJ_pyr3_11, measurevar="diameter", na.rm = TRUE)
MJ_P11$Animal<-"MJ"
MJ_P11$Vessel<-"SPOT5_v2"
MJ_P11$Drug<-"pyr3"

MJ_P12<-summarySE(MJ_pyr3_12, measurevar="diameter", na.rm = TRUE)
MJ_P12$Animal<-"MJ"
MJ_P12$Vessel<-"SPOT5_v3"
MJ_P12$Drug<-"pyr3"

MJ_P13<-summarySE(MJ_pyr3_13, measurevar="diameter", na.rm = TRUE)
MJ_P13$Animal<-"MJ"
MJ_P13$Vessel<-"SPOT5_v4"
MJ_P13$Drug<-"pyr3"

MJ_P14<-summarySE(MJ_pyr3_14, measurevar="diameter", na.rm = TRUE)
MJ_P14$Animal<-"MJ"
MJ_P14$Vessel<-"SPOT5_v5"
MJ_P14$Drug<-"pyr3"

MJ_P15<-summarySE(MJ_pyr3_15, measurevar="diameter", na.rm = TRUE)
MJ_P15$Animal<-"MJ"
MJ_P15$Vessel<-"SPOT5_v6"
MJ_P15$Drug<-"pyr3"

MJ_P16<-summarySE(MJ_pyr3_16, measurevar="diameter", na.rm = TRUE)
MJ_P16$Animal<-"MJ"
MJ_P16$Vessel<-"SPOT6_v1"
MJ_P16$Drug<-"pyr3"

#combine data together for each treatment
baseline<-rbind(RR_BL1,RR_BL2,RR_BL3,RR_BL4,RR_BL5,RR_BL6,RR_BL7,RR_BL8,RR_BL9,MJ_BL1,MJ_BL2,MJ_BL3,MJ_BL4,MJ_BL5,MJ_BL6,MJ_BL7,MJ_BL8,MJ_BL9,MJ_BL10,MJ_BL11,MJ_BL12,MJ_BL13,MJ_BL14,MJ_BL15,MJ_BL16,MJ_BL17,MJ_BL18,MJ_BL19,MJ_BL20,MJ_BL21,MJ_BL22) #L_baseline1,L_baseline2,L_baseline3,L_baseline4,L_baseline5,L_baseline6,L_baseline7,L_baseline8,L_baseline9,L_baseline10,L_baseline11,JJ_baseline1,JJ_baseline2,JJ_baseline3,JJ_baseline4,JJ_baseline5,JJ_baseline6,OM_baseline1,OM_baseline2,OM_baseline3,OM_baseline4,OM_baseline5,OM_baseline6,OM_baseline7,NL_baseline1,NL_baseline2,NL_baseline3,NL_baseline4)
nimodipine<-rbind(RR_NP1,RR_NP2,RR_NP3,RR_NP4,RR_NP5,RR_NP6,RR_NP7,RR_NP8,RR_NP9,RR_NP10,RR_NP11,RR_NP12,RR_NP13,RR_NP14,RR_NP15,MJ_NP1,MJ_NP2,MJ_NP3,MJ_NP4,MJ_NP5,MJ_NP6,MJ_NP7,MJ_NP8,MJ_NP9,MJ_NP10,MJ_NP11,MJ_NP12,MJ_NP13,MJ_NP14)#L_nimodipine1)
pyr3<-rbind(RR_P1,RR_P2,RR_P3,RR_P4,RR_P5,RR_P6,RR_P7,RR_P8,RR_P9,MJ_P1,MJ_P2,MJ_P3,MJ_P4,MJ_P5,MJ_P6,MJ_P7,MJ_P8,MJ_P9,MJ_P10,MJ_P11,MJ_P12,MJ_P13,MJ_P14,MJ_P15,MJ_P16)

AllData2<-rbind(baseline,nimodipine,pyr3)
AllData2$Drug<-as.factor(AllData2$Drug)

#pull out the spot name
pos= regexpr('SPOT', AllData2$Vessel)
AllData2$SpotName<-substr(AllData2$Vessel,pos, pos+4)

AllData1$Spot_drug<-paste(AllData1$Animal,AllData1$Vessel,AllData1$Drug, sep="_")
AllData2$Spot_drug<-paste(AllData2$Animal, AllData2$Vessel,AllData2$Drug, sep="_")

uniqueVessels<-intersect(AllData2$Spot_drug,AllData1$Spot_drug)

AllData1<-subset(AllData1, Spot_drug %in% uniqueVessels)
AllData2<-subset(AllData2, Spot_drug %in% uniqueVessels)

Vel.means<- ddply(AllData1, c("Animal","Vessel","Drug","Spot_drug"), 
                  summarise, velocity_mean = mean(velocity, na.rm=TRUE))

Diam.means<- ddply(AllData2, c("Animal","Vessel","Drug","Spot_drug"), 
                   summarise, diameter_mean = mean(diameter, na.rm=TRUE))

Vel.means<-merge(Vel.means,Diam.means,id = ~ Spot_drug)

Vel.means$SpotName<-paste(Vel.means$Animal, Vel.means$Vessel, sep="_")

## plot of velocity vs. diameter


ggplot(data=Vel.means, aes(x=diameter_mean, y=velocity_mean, colour=Drug)) +
  geom_point(stat="identity",  size= 3) +
  ylab("Velocity [um/s]") +
  xlab("Diameter [um]") +
  scale_colour_manual(
    values=c("black", "red", "blue")) +
  max.theme


#line plots
ggplot(Vel.means[Vel.means$Drug!="nimodipine",],aes(x=Drug,y=diameter_mean)) +
  geom_point(aes(x=Drug,y=diameter_mean),stat="identity", size=3)+
  geom_line(aes(x=Drug,y=diameter_mean,group=SpotName))+
  #geom_errorbar(aes(ymin=diameter_mean-se, ymax=diameter_mean+se), colour="black", width=.1) +
  ylab("Diameter [um]") +
  max.theme


ggplot(Vel.means[Vel.means$Drug!="pyr3",],aes(x=Drug,y=diameter_mean)) +
  geom_point(aes(x=Drug,y=diameter_mean),stat="identity", size=3)+
  geom_line(aes(x=Drug,y=diameter_mean,group=SpotName))+
  #geom_errorbar(aes(ymin=diameter_mean-se, ymax=diameter_mean+se), colour="black", width=.1) +
  ylab("Diameter [um]") +
  max.theme

ggplot(Vel.means[Vel.means$Drug!="baseline",],aes(x=Drug,y=diameter_mean)) +
  geom_point(aes(x=Drug,y=diameter_mean),stat="identity", size=3)+
  geom_line(aes(x=Drug,y=diameter_mean,group=SpotName))+
  #geom_errorbar(aes(ymin=diameter_mean-se, ymax=diameter_mean+se), colour="black", width=.1) +
  ylab("Diameter [um]") +
  max.theme

ggplot(Vel.means[Vel.means$Drug!="nimodipine",],aes(x=Drug,y=velocity_mean)) +
  geom_point(aes(x=Drug,y=velocity_mean),stat="identity", size=3)+
  geom_line(aes(x=Drug,y=velocity_mean,group=SpotName))+
  #geom_errorbar(aes(ymin=diameter_mean-se, ymax=diameter_mean+se), colour="black", width=.1) +
  ylab("Velocity [um/s]") +
  max.theme

ggplot(Vel.means[Vel.means$Drug!="pyr3",],aes(x=Drug,y=velocity_mean)) +
  geom_point(aes(x=Drug,y=velocity_mean),stat="identity", size=3)+
  geom_line(aes(x=Drug,y=velocity_mean,group=SpotName))+
  #geom_errorbar(aes(ymin=diameter_mean-se, ymax=diameter_mean+se), colour="black", width=.1) +
  ylab("Velocity [um/s]") +
  max.theme

ggplot(Vel.means[Vel.means$Drug!="baseline",],aes(x=Drug,y=velocity_mean)) +
  geom_point(aes(x=Drug,y=velocity_mean),stat="identity", size=3)+
  geom_line(aes(x=Drug,y=velocity_mean,group=SpotName))+
  #geom_errorbar(aes(ymin=diameter_mean-se, ymax=diameter_mean+se), colour="black", width=.1) +
  ylab("Velocity [um/s]") +
  max.theme
