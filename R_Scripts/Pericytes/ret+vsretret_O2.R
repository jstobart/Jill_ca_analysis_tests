
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
allData$Genotype<- factor(allData$Genotype,levels = c("Ret+", "RetRet"))

# remove points from tissue or with noisy data
allData<-subset(allData, ExcludeData==0)


# unique animal and spot name
allData$FOVName <-paste(allData$Animal, allData$Spot, sep= "_")
allData$PointName <-paste(allData$Animal, allData$Spot, allData$AdjustedDepth, allData$Point, sep= "_")
allData$BranchName <-paste(allData$Animal, allData$Spot, allData$VesselType, sep= "_")

#allData$BranchOrder<- as.factor(as.character(allData$BranchOrder))
allData$Animal<- as.character(allData$Animal)
allData$Spot<- as.character(allData$Spot)


#change one branch name that is messed up
allData$BranchName[allData$BranchName=="Prr2_2015_10_27_V3-V4"]="Prr2_2015_10_27_V4"


# outlier test
source("http://goo.gl/UUyEzD")
outlierKD(allData, pO2)
y
allData<-allData[complete.cases(allData$pO2),]


##############

# classify as high and low pO2 based on the mean value
surfacevessels<-subset(allData, AdjustedDepth<100 & BranchOrder<3)

surfacevessels.mean<- ddply(surfacevessels, c("Animal","Spot","FOVName","BranchName","Genotype"), 
                            summarise, meanpO2=mean(pO2))

deepervessels<-subset(allData, AdjustedDepth>150)

deepervessels.mean<- ddply(deepervessels, c("Animal","Spot","FOVName","BranchName","Genotype"), 
                           summarise, meanpO2=mean(pO2))

# vessel classification based on O2 levels
surface.meanpO2=mean(surfacevessels$pO2)
deeper.meanpO2=mean(deepervessels$pO2)

surfacevessels.mean$O2_type="artery"
surfacevessels.mean$O2_type[surfacevessels.mean$meanpO2<surface.meanpO2]="vein"


#############
# combine averages from surface and deep layers
surface.vs.deep<-merge(surfacevessels.mean[, c("BranchName", "O2_type", "meanpO2")], deepervessels.mean, by="BranchName", all.x=TRUE)
surface.vs.deep<-surface.vs.deep[complete.cases(surface.vs.deep$meanpO2.y),]
surface.vs.deep$O2_type<-as.factor(surface.vs.deep$O2_type)

surface.vs.deep2<-melt(surface.vs.deep) #, id.vars=c("BranchName","Genotype"), variable.name=c(""))


df1A.means<- summarySE(surface.vs.deep2, measurevar="value", groupvars=c("Genotype"))
df1B.means<- summarySE(surface.vs.deep2, measurevar="value", groupvars=c("Genotype","variable"))
df1C.means<- summarySE(surface.vs.deep2, measurevar="value", groupvars=c("Genotype","variable","O2_type"))

ggplot(surface.vs.deep2, aes(x=interaction(variable,Genotype), y=value, group = BranchName)) +
  geom_point(aes(colour = O2_type), size=2)+
  geom_line(aes(colour = O2_type), size=0.5)+
  ggtitle("mean pO2 shallow and deep") +
  xlab("shallow or deep") + 
  ylab("pO2 [mmHg]") + 
  ylim(15, 65) +
  #scale_colour_manual(values=c("black", "red"), guide=FALSE) + 
  max.theme


ggplot(data=df1B.means, aes(x=variable, y=value, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("mean pO2 [mmHg]") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df1C.means, aes(x=interaction(variable,O2_type), y=value, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("mean pO2 [mmHg]") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

#############
# difference between deep and shallow layers for matched vessels
surface.vs.deep$difference<-surface.vs.deep$meanpO2.y-surface.vs.deep$meanpO2.x

df.pO2difference<-summarySE(surface.vs.deep, measurevar="difference", groupvars=c("Genotype","O2_type"))

ggplot(data=df.pO2difference, aes(x=O2_type, y=difference, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=difference-se, ymax=difference+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("pO2 difference- deeper to surface [mmHg]") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme


#Stats
## pO2 and genotype or branch group
pO2.diff.null = lmer(difference ~ (1|Animal) + (1|FOVName), surface.vs.deep,REML=FALSE)
pO2.diff.model1 = lmer(difference ~ Genotype + (1|Animal) + (1|FOVName), surface.vs.deep,REML=FALSE)
pO2.diff.model2B = lmer(difference ~ O2_type + (1|Animal) + (1|FOVName), surface.vs.deep,REML=FALSE)
pO2.diff.model3B = lmer(difference ~ Genotype + O2_type + (1|Animal) + (1|FOVName), surface.vs.deep,REML=FALSE)
pO2.diff.model4B = lmer(difference ~ Genotype * O2_type + (1|Animal) + (1|FOVName), surface.vs.deep,REML=FALSE)
pO2.diff.anova <- anova(pO2.diff.null, pO2.diff.model1, pO2.diff.model2B,
                   pO2.diff.model3B,pO2.diff.model4B)
print(pO2.diff.anova)

# p values
pO2.diff.O2type <- lsmeans(pO2.diff.model4B, pairwise ~ Genotype*O2_type, glhargs=list())
summary(pO2.diff.O2type )


######
# apply surface vessel classification to the vessels in the all data table
allData.vesselType<-merge(allData, surfacevessels.mean[, c("BranchName", "O2_type")], by="BranchName", all=TRUE)
#remove vessels that do not appear at surface
allData.vesselType<-allData.vesselType[complete.cases(allData.vesselType$O2_type),]

# group branch order to simplify analysis

# arteriolar side
allData.vesselType$BranchGroup="arteriole"
allData.vesselType$BranchGroup[allData.vesselType$BranchOrder>0 & allData.vesselType$BranchOrder<=4 & allData.vesselType$O2_type=="artery"]<-"ensheathing_PC"
allData.vesselType$BranchGroup[allData.vesselType$BranchOrder>4 & allData.vesselType$O2_type=="artery"]<-"capillary_PC"

#venous side
allData.vesselType$BranchGroup[allData.vesselType$BranchOrder<=3 & allData.vesselType$O2_type=="vein"]<-"venule_PC"
allData.vesselType$BranchGroup[allData.vesselType$BranchOrder==0 & allData.vesselType$O2_type=="vein"]<-"vein"
allData.vesselType$BranchGroup[allData.vesselType$BranchOrder>3 & allData.vesselType$O2_type=="vein"]<-"capillary_PC"

allData.vesselType$BranchGroup<- as.factor(allData.vesselType$BranchGroup)
allData.vesselType$BranchGroup<- factor(allData.vesselType$BranchGroup,levels = c("arteriole", "ensheathing_PC", "capillary_PC", "venule_PC", "vein"))


###############################
# plots

ggplot(allData.vesselType, aes(x=pO2, fill=Genotype)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of pO2")

ggplot(surfacevessels.mean, aes(x=meanpO2, fill=Genotype)) + geom_histogram(binwidth=2, position="dodge") +
  ggtitle("Distribution of mean pO2 per vessel branch")

ggplot(allData.vesselType, aes(x = O2_type, y = pO2, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("pO2 [mmHg]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(allData.vesselType, aes(x = BranchGroup, y = pO2, fill = Genotype)) + 
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
df1C<- summarySE(allData.vesselType, measurevar="pO2", groupvars=c("Genotype","BranchGroup"))

allData.vesselType$AdjustedDepth<-as.factor(allData.vesselType$AdjustedDepth)
df1D<- summarySE(allData.vesselType, measurevar="pO2", groupvars=c("Genotype","BranchGroup","AdjustedDepth"))

ggplot(data=df1A, aes(x=Genotype, y=pO2, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=pO2-se, ymax=pO2+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("pO2 [mmHg]") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df1C, aes(x=BranchGroup, y=pO2, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=pO2-se, ymax=pO2+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("pO2[mmHg]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme


#line graphs for each type of vessel branch
ggplot(df1D, aes(x=AdjustedDepth, y=pO2)) +
  geom_point(aes(colour = interaction(Genotype,BranchGroup)), size=2)+
  geom_line(aes(colour = interaction(Genotype,BranchGroup)), size=1)+
  ggtitle("mean pO2 at different depths") +
  xlab("Depth") + 
  ylab("pO2 [mmHg]") + 
  ylim(0, 70) +
  #scale_colour_manual(values=c("black", "red"), guide=FALSE) + 
  max.theme


######
#Stats
## pO2 and genotype or branch group
pO2.null = lmer(pO2 ~ (1|Animal) + (1|Spot), allData.vesselType,REML=FALSE)
pO2.model1 = lmer(pO2~ Genotype + (1|Animal) + (1|Spot), allData.vesselType,REML=FALSE)
pO2.model2B = lmer(pO2~ BranchGroup + (1|Animal) + (1|Spot), allData.vesselType,REML=FALSE)
pO2.model3B = lmer(pO2~ Genotype + BranchGroup + (1|Animal) + (1|Spot), allData.vesselType,REML=FALSE)
pO2.model4B = lmer(pO2~ Genotype * BranchGroup + (1|Animal) + (1|Spot), allData.vesselType,REML=FALSE)
pO2.anova <- anova(pO2.null, pO2.model1, pO2.model2B,
                   pO2.model3B,pO2.model4B)
print(pO2.anova)

# p values
pO2.BranchGroup <- lsmeans(pO2.model4B, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(pO2.BranchGroup)


pO2.null = lmer(pO2 ~ (1|Animal) + (1|Spot), allData.vesselType,REML=FALSE)
pO2.model1 = lmer(pO2~ Genotype + (1|Animal) + (1|Spot), allData.vesselType,REML=FALSE)
pO2.model2A = lmer(pO2~ O2_type + (1|Animal) + (1|Spot), allData.vesselType,REML=FALSE)
pO2.model2B = lmer(pO2~ BranchGroup + (1|Animal) + (1|Spot), allData.vesselType,REML=FALSE)
pO2.model3A = lmer(pO2~ Genotype + O2_type + (1|Animal) + (1|Spot), allData.vesselType,REML=FALSE)
pO2.model3B = lmer(pO2~ Genotype + BranchGroup + (1|Animal) + (1|Spot), allData.vesselType,REML=FALSE)
pO2.model4A = lmer(pO2~ Genotype * O2_type + (1|Animal) + (1|Spot), allData.vesselType,REML=FALSE)
pO2.model4B = lmer(pO2~ Genotype * BranchGroup + (1|Animal) + (1|Spot), allData.vesselType,REML=FALSE)
pO2.anova <- anova(pO2.null, pO2.model1, pO2.model2A,pO2.model2B,
                   pO2.model3A,pO2.model3B,pO2.model4A,pO2.model4B)
print(pO2.anova)


pO2.O2type <- lsmeans(pO2.model4A, pairwise ~ Genotype*O2_type, glhargs=list())
summary(pO2.BranchGroup)


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
