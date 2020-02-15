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
# 5 months
MJ_baseline1 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_06_11_KX/2019_06_11_SPOT2_v1_diam_LineScan_KX_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline2 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_06_11_KX/2019_06_11_SPOT2_v2_diam_LineScan_KX_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline3 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_06_11_KX/2019_06_11_SPOT2_v3_diam_LineScan_KX_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline4 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_06_11_KX/2019_06_11_SPOT4_v1_diam_LineScan_KX_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline5 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_07_09/2019_07_09_SPOT1_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline6 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_07_09/2019_07_09_SPOT1_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline7 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_07_09/2019_07_09_SPOT1_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline8 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_07_11/2019_07_11_SPOT3_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline9 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_07_11/2019_07_11_SPOT6_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline10 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_07_22/2019_07_22_SPOT4_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline11 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_07_22/2019_07_22_SPOT4_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline12 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_07_22/2019_07_22_SPOT5_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline13 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_07_22/2019_07_22_SPOT5_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline14 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_07_22/2019_07_22_SPOT5_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline15 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_07_22/2019_07_22_SPOT5_v4_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline16 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_07_22/2019_07_22_SPOT5_v5_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline17 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_07_22/2019_07_22_SPOT5_v6_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")

# 7 months
MJ_baseline18 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_08_19/2019_08_19_SPOT1_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline19 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_08_19/2019_08_19_SPOT1_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline20 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_08_19/2019_08_19_SPOT1_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline21 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_08_19/2019_08_19_SPOT3_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline22 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_08_19/2019_08_19_SPOT6_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_baseline23 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_08_19/2019_08_19_SPOT3_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")

# 9 months
MJ_9month_1 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_11_07/2019_11_07_SPOT1_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_9month_2 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_11_07/2019_11_07_SPOT1_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_9month_3 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_11_07/2019_11_07_SPOT1_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_9month_4 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_11_07/2019_11_07_SPOT3_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_9month_5 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_11_07/2019_11_07_SPOT6_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_9month_6 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_11_07/2019_11_07_SPOT3_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_9month_7 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_11_07/2019_11_07_SPOT4_v1_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_9month_8 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_11_07/2019_11_07_SPOT4_v2_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_9month_9 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2019_11_07/2019_11_07_SPOT5_v3_diam_LineScan_diameterFWHM.csv", header=TRUE, sep = ",")

# 11 months
MJ_11month_1 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2020_01_06/spot3_v1_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_11month_2 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2020_01_06/spot3_v3_diameterFWHM.csv", header=TRUE, sep = ",")
MJ_11month_3 <- read.table("E:/Jill/Data/Winnipeg/Pericytes/In vivo 2P/RCaMP Mice/Results/Bloodflow Results/RCaMP 59721 MJ/Diameter/2020_01_06/spot6_v1_diameterFWHM.csv", header=TRUE, sep = ",")

#####
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
MJ_BL18$Drug<-"7months"

MJ_BL19<-summarySE(MJ_baseline19, measurevar="diameter", na.rm = TRUE)
MJ_BL19$Animal<-"MJ"
MJ_BL19$Vessel<-"SPOT1_v2"
MJ_BL19$Drug<-"7months"

MJ_BL20<-summarySE(MJ_baseline20, measurevar="diameter", na.rm = TRUE)
MJ_BL20$Animal<-"MJ"
MJ_BL20$Vessel<-"SPOT1_v3"
MJ_BL20$Drug<-"7months"

MJ_BL21<-summarySE(MJ_baseline21, measurevar="diameter", na.rm = TRUE)
MJ_BL21$Animal<-"MJ"
MJ_BL21$Vessel<-"SPOT3_v1"
MJ_BL21$Drug<-"7months"

MJ_BL22<-summarySE(MJ_baseline22, measurevar="diameter", na.rm = TRUE)
MJ_BL22$Animal<-"MJ"
MJ_BL22$Vessel<-"SPOT6_v1"
MJ_BL22$Drug<-"7months"

MJ_BL23<-summarySE(MJ_baseline23, measurevar="diameter", na.rm = TRUE)
MJ_BL23$Animal<-"MJ"
MJ_BL23$Vessel<-"SPOT3_v3"
MJ_BL23$Drug<-"7months"

MJ_9m_1<-summarySE(MJ_9month_1, measurevar="diameter", na.rm = TRUE)
MJ_9m_1$Animal<-"MJ"
MJ_9m_1$Vessel<-"SPOT1_v1"
MJ_9m_1$Drug<-"9months"

MJ_9m_2<-summarySE(MJ_9month_2, measurevar="diameter", na.rm = TRUE)
MJ_9m_2$Animal<-"MJ"
MJ_9m_2$Vessel<-"SPOT1_v2"
MJ_9m_2$Drug<-"9months"

MJ_9m_3<-summarySE(MJ_9month_3, measurevar="diameter", na.rm = TRUE)
MJ_9m_3$Animal<-"MJ"
MJ_9m_3$Vessel<-"SPOT1_v3"
MJ_9m_3$Drug<-"9months"

MJ_9m_4<-summarySE(MJ_9month_4, measurevar="diameter", na.rm = TRUE)
MJ_9m_4$Animal<-"MJ"
MJ_9m_4$Vessel<-"SPOT3_v1"
MJ_9m_4$Drug<-"9months"

MJ_9m_5<-summarySE(MJ_9month_5, measurevar="diameter", na.rm = TRUE)
MJ_9m_5$Animal<-"MJ"
MJ_9m_5$Vessel<-"SPOT6_v1"
MJ_9m_5$Drug<-"9months"

MJ_9m_6<-summarySE(MJ_9month_6, measurevar="diameter", na.rm = TRUE)
MJ_9m_6$Animal<-"MJ"
MJ_9m_6$Vessel<-"SPOT3_v3"
MJ_9m_6$Drug<-"9months"

MJ_9m_7<-summarySE(MJ_9month_7, measurevar="diameter", na.rm = TRUE)
MJ_9m_7$Animal<-"MJ"
MJ_9m_7$Vessel<-"SPOT4_v1"
MJ_9m_7$Drug<-"9months"

MJ_9m_8<-summarySE(MJ_9month_8, measurevar="diameter", na.rm = TRUE)
MJ_9m_8$Animal<-"MJ"
MJ_9m_8$Vessel<-"SPOT4_v2"
MJ_9m_8$Drug<-"9months"

MJ_9m_9<-summarySE(MJ_9month_9, measurevar="diameter", na.rm = TRUE)
MJ_9m_9$Animal<-"MJ"
MJ_9m_9$Vessel<-"SPOT5_v3"
MJ_9m_9$Drug<-"9months"

# month 11
MJ_11m_1<-summarySE(MJ_11month_1, measurevar="diameter", na.rm = TRUE)
MJ_11m_1$Animal<-"MJ"
MJ_11m_1$Vessel<-"SPOT3_v1"
MJ_11m_1$Drug<-"11months"

MJ_11m_2<-summarySE(MJ_11month_2, measurevar="diameter", na.rm = TRUE)
MJ_11m_2$Animal<-"MJ"
MJ_11m_2$Vessel<-"SPOT3_v3"
MJ_11m_2$Drug<-"11months"

MJ_11m_3<-summarySE(MJ_11month_3, measurevar="diameter", na.rm = TRUE)
MJ_11m_3$Animal<-"MJ"
MJ_11m_3$Vessel<-"SPOT6_v1"
MJ_11m_3$Drug<-"11months"

#combine data together for each treatment
baseline<-rbind(MJ_BL1,MJ_BL2,MJ_BL3,MJ_BL4,MJ_BL5,MJ_BL6,MJ_BL7,MJ_BL8,MJ_BL9,MJ_BL10,MJ_BL11,MJ_BL12,MJ_BL13,MJ_BL14,MJ_BL15,MJ_BL16,MJ_BL17)
month7<- rbind(MJ_BL18,MJ_BL19,MJ_BL20,MJ_BL21,MJ_BL22, MJ_BL23)
month9<-rbind(MJ_9m_1,MJ_9m_2,MJ_9m_3,MJ_9m_4,MJ_9m_5,MJ_9m_6,MJ_9m_7,MJ_9m_8,MJ_9m_9)
month11<-rbind(MJ_11m_1,MJ_11m_2,MJ_11m_3)

AllData<-rbind(baseline,month7, month9, month11)
AllData$Drug<-as.factor(AllData$Drug)
AllData$Drug<-factor(AllData$Drug, levels =c("baseline", "7months", "9months", "11months"))

#pull out the spot name
pos= regexpr('SPOT', AllData$Vessel)
AllData$SpotName<-substr(AllData$Vessel,pos, pos+4)


###############################
#histograms
ggplot(AllData, aes(x=diameter, fill=Drug))+ geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("Distribution of Diameter")

#########

#diameter
df1A <- summarySE(AllData, measurevar="diameter", groupvars=c("Drug"))

ggplot(data=df1A, aes(x=Drug, y=diameter, fill=Drug)) +
  geom_errorbar(aes(ymin=diameter-se, ymax=diameter+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("diameter") +
  scale_fill_manual(
    values=c("black", "red", "blue", "green")) +
  max.theme


# only the same vessels over time

AllData2<-subset(AllData, Vessel %in% month11$Vessel)


df1B <- summarySE(AllData2, measurevar="diameter", groupvars=c("Drug"))

ggplot(data=df1B, aes(x=Drug, y=diameter, fill=Drug)) +
  geom_errorbar(aes(ymin=diameter-se, ymax=diameter+se), colour="black", width=.5, size= 1, position=position_dodge(0.9)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width=1, size= 1) +
  ylab("diameter") +
  scale_fill_manual(
    values=c("black", "red", "blue", "green")) +
  max.theme



######
#stats

#diameter
diameter.null = lmer(diameter ~ (1|SpotName) + (1|Vessel), AllData,REML=FALSE)
diameter.model1 = lmer(diameter~ Drug + (1|SpotName) + (1|Vessel), AllData,REML=FALSE)
diameter.anova <- anova(diameter.null, diameter.model1)
print(diameter.anova)

# p values
diam.drug <- lsmeans(diameter.model1, pairwise ~ Drug, glhargs=list())
summary(diam.drug)
