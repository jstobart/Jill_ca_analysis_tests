
library("lme4")
library("lmerTest")
library("lattice")
library("plyr")
library("dplyr")
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
ret_t <- read.csv("E:/Data/Pericyte_project/Two-photon-data/Ret_ret_Mice/O2/Extracted values + Results/Jill-Extracted-values+ret_wt.csv", header=TRUE, sep = ",")
ret_ret <- read.csv("E:/Data/Pericyte_project/Two-photon-data/Ret_ret_Mice/O2/Extracted values + Results/Jill-Extracted-values+ret_ret.csv", header=TRUE, sep = ",")

#str(ret_t)
#str(ret_ret)

ret_t$Genotype<-"Ret+"
ret_ret$Genotype<-"RetRet"

#options(digits=9)
ret_t$pO2<-as.numeric(as.character(ret_t$pO2))
#ret_ret$pO2<-as.numeric(ret_ret$pO2)

#combine data sets
allData<-rbind(ret_t,ret_ret)
#str(allData)
allData$Genotype<- as.factor(allData$Genotype)

# remove points from tissue or with noisy data
allData<-subset(allData, ExcludeData==0)



#########
# unique animal and spot name
allData$PointName <-paste(allData$Animal, allData$Spot, allData$AdjustedDepth, allData$Point, sep= "_")
allData$BranchName <-paste(allData$Animal, allData$Spot, allData$VesselType, sep= "_")

#allData$Genotype<- factor(allData$Genotype,levels = c("Ret+", "RetRet"))
allData$BranchOrder<- as.factor(as.character(allData$BranchOrder))
allData$Animal<- as.character(allData$Animal)
allData$Spot<- as.character(allData$Spot)


# group branch order to simplify analysis
#eGFPpos$BranchGroup[eGFPpos$BranchOrder<=0]<-"arteriole"
#eGFPpos$BranchGroup[eGFPpos$BranchOrder<=4]<-"sm-ensheathing"
#eGFPpos$BranchGroup[eGFPpos$BranchOrder>4]<-"capillary_PC"
#eGFPneg$BranchGroup[eGFPneg$BranchOrder>3]<-"capillary_PC"
#eGFPneg$BranchGroup[eGFPneg$BranchOrder<=3]<-"venule_PC"

#allData$BranchGroup<- as.factor(allData$BranchGroup)


# vessel classification based on O2 levels

# classify as high and low pO2 based on the median value
medianpO2=median(allData$pO2)

allData$O2_type="artery"
allData$O2_type[allData$pO2<medianpO2]="vein"

surfacevessels<-subset(allData, AdjustedDepth<100)

#surfacevessels$O2_type="artery"
#surfacevessels$O2_type[surfacevessels$pO2<medianpO2]="vein"
#surfacevessels$O2_type<-as.factor(surfacevessels$O2_type)

#allData$O2_type=
allData.vesselType<-merge(allData, surfacevessels[, c("BranchName", "O2_type")], by="BranchName", all.x=TRUE)
allData.vesselType<-distinct(allData.vesselType, PointName, .keep_all = TRUE)

###############################

ggplot(allData.vesselType, aes(x = O2_type, y = pO2, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("pO2 [mmHg]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

#line graphs for each vessel
ggplot(allData.vesselType[allData.vesselType$Genotype=="RetRet"& allData.vesselType$O2_type=="artery",], aes(x=AdjustedDepth, y=pO2)) +
  geom_point(aes(colour = BranchName), size=2)+
  geom_line(aes(colour = BranchName), size=1)+
  ggtitle("Retret-pO2 vs depth for all artery") +
  xlab("Depth") + 
  ylab("pO2 [mmHg]") + 
  ylim(0, 70) +
  #scale_colour_manual(values=c("black", "red"), guide=FALSE) + 
  max.theme

ggplot(allData.vesselType[allData.vesselType$Genotype=="RetRet"& allData.vesselType$O2_type=="vein",], aes(x=AdjustedDepth, y=pO2)) +
  geom_point(aes(colour = BranchName), size=2)+
  geom_line(aes(colour = BranchName), size=1)+
  ggtitle("Retret-pO2 vs depth for all veins") +
  xlab("Depth") + 
  ylab("pO2 [mmHg]") + 
  ylim(0, 70) +
  #scale_colour_manual(values=c("black", "red"), guide=FALSE) + 
  max.theme

ggplot(allData.vesselType[allData.vesselType$Genotype=="Ret+"& allData.vesselType$O2_type=="artery",], aes(x=AdjustedDepth, y=pO2)) +
  geom_point(aes(colour = BranchName), size=2)+
  geom_line(aes(colour = BranchName), size=1)+
  ggtitle("Ret+-pO2 vs depth for all artery") +
  xlab("Depth") + 
  ylab("pO2 [mmHg]") + 
  ylim(0, 70) +
  #scale_colour_manual(values=c("black", "red"), guide=FALSE) + 
  max.theme

ggplot(allData.vesselType[allData.vesselType$Genotype=="Ret+"& allData.vesselType$O2_type=="vein",], aes(x=AdjustedDepth, y=pO2)) +
  geom_point(aes(colour = BranchName), size=2)+
  geom_line(aes(colour = BranchName), size=1)+
  ggtitle("Ret+-pO2 vs depth for all veins") +
  xlab("Depth") + 
  ylab("pO2 [mmHg]") + 
  ylim(0, 70) +
  #scale_colour_manual(values=c("black", "red"), guide=FALSE) + 
  max.theme


#######
# mean pO2

df1A<- summarySE(allData.vesselType, measurevar="pO2", groupvars=c("Genotype"))
df1B<- summarySE(allData.vesselType, measurevar="pO2", groupvars=c("Genotype","BranchOrder"))
#df1C<- summarySE(allData, measurevar="pO2", groupvars=c("Genotype","BranchGroup"))
df1D<- summarySE(allData.vesselType, measurevar="pO2", groupvars=c("Genotype","O2_type"))


ggplot(data=df1A, aes(x=Genotype, y=pO2, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=pO2-se, ymax=pO2+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("pO2 [mmHg]") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df1B, aes(x=BranchOrder, y=pO2, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=pO2-se, ymax=pO2+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("pO2[mmHg]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme


ggplot(data=df1C, aes(x=BranchGroup, y=pO2, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=pO2-se, ymax=pO2+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("pO2[mmHg]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme

ggplot(data=df1D, aes(x=O2_type, y=pO2, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=pO2-se, ymax=pO2+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("pO2[mmHg]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme


######
#Stats
## pO2 and genotype or branch order
pO2.null = lmer(pO2 ~ (1|Animal) + (1|Spot), allData,REML=FALSE)
#pO2.model1 = lmer(pO2~ Genotype + (1|Animal) + (1|Spot), allData,REML=FALSE)
pO2.model2A = lmer(pO2~ BranchOrder + (1|Animal) + (1|Spot), allData,REML=FALSE)
pO2.model2B = lmer(pO2~ BranchGroup + (1|Animal) + (1|Spot), allData,REML=FALSE)
#pO2.model3 = lmer(pO2~ Genotype + BranchOrder + (1|Animal) + (1|Spot), allData,REML=FALSE)
#pO2.model4 = lmer(pO2~ Genotype * BranchOrder + (1|Animal) + (1|Spot), allData,REML=FALSE)
pO2.anova <- anova(pO2.null, pO2.model2A,pO2.model2B)
#pO2.anova <- anova(pO2.null, pO2.model1,pO2.model2,pO2.model3,pO2.model4)
print(pO2.anova)
# p values
pO2.BranchGroup <- lsmeans(pO2.model2B, pairwise ~ BranchGroup, glhargs=list())
summary(pO2.BranchGroup)

#pO2.Genotype <- lsmeans(pO2.model1, pairwise ~ Genotype, glhargs=list())
#summary(pO2.Genotype)

#pO2.Genotype_BO <- lsmeans(pO2.model4, pairwise ~ Genotype*BranchOrder, glhargs=list())
#summary(pO2.Genotype_BO)

# Look at the model residuals
plot(pO2.model1)
plot(pO2.model2)
plot(pO2.model2B)

# pO2eter, genotype and depth
pO2.nullB = lmer(pO2 ~ (1|Animal) + (1|Spot), allData,REML=FALSE)
#pO2.model1B = lmer(pO2~ Genotype + (1|AnimalName) + (1|Spot), allData,REML=FALSE)
pO2.model2C = lmer(pO2~ AdjustedDepth + (1|Animal) + (1|Spot), allData,REML=FALSE)
#pO2.model3B = lmer(pO2~ Genotype + Depth + (1|AnimalName) + (1|Spot), allData,REML=FALSE)
#pO2.model4B = lmer(pO2~ Genotype * Depth + (1|AnimalName) + (1|Spot), allData,REML=FALSE)
pO2.anovaB <- anova(pO2.nullB, pO2.model2C)
#pO2.anovaB <- anova(pO2.nullB, pO2.model1B,pO2.model2B,pO2.model3B,pO2.model4B)
print(pO2.anovaB)


##############
