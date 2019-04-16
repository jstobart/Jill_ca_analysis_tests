
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
#ret_t <- read.csv("E:/Data/Pericyte_project/Two-photon-data/Ret_ret_Mice/O2/Extracted values + Results/Jill-Extracted-values+ret_wt.csv", header=TRUE, sep = ",")
#ret_ret <- read.csv("E:/Data/Pericyte_project/Two-photon-data/Ret_ret_Mice/O2/Extracted values + Results/Jill-Extracted-values+ret_ret.csv", header=TRUE, sep = ",")

ret_t <- read.csv("D:/Data/Pericytes/Results/Jill-Extracted-values+ret_wt.csv", header=TRUE, sep = ",")
ret_ret <- read.csv("D:/Data/Pericytes/Results/Jill-Extracted-values+ret_ret.csv", header=TRUE, sep = ",")


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
#source("http://goo.gl/UUyEzD")
#outlierKD(allData, pO2)
#y
#allData<-allData[complete.cases(allData$pO2),]


##############

# classify as high and low pO2 based on the mean value
surfacevessels<-subset(allData, AdjustedDepth<100)

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
    values=c("blue", "red"),guide=FALSE)+
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

ggplot() +
  geom_jitter(data=surface.vs.deep, aes(y=difference, x=interaction(Genotype,O2_type), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df.pO2difference, aes(x=interaction(Genotype,O2_type), y=difference,ymin=difference, ymax=difference, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df.pO2difference, aes(x=interaction(Genotype,O2_type), ymin=difference-se, ymax=difference+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("O2 gradient between upper and lower cortical layers")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot() +
  geom_bar(data=df.pO2difference, aes(x=interaction(Genotype,O2_type), y=difference, colour=Genotype), fill="white" ,stat="identity", width=0.5, position=position_dodge()) +
  geom_jitter(data=surface.vs.deep, aes(y=difference, x=interaction(Genotype,O2_type), colour=Genotype), 
              position=position_jitter(width=0.12), size=4)+
  geom_errorbar(data=df.pO2difference, aes(x=interaction(Genotype,O2_type), ymin=difference-se, ymax=difference+se), size=2,colour="black", width=0.1,  position=position_dodge()) +
  geom_hline(yintercept = 0)+
  ggtitle("difference all vessels")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
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


ggplot(allData.vesselType, aes(x=pO2, y=..density..,fill=Genotype)) +
  geom_histogram(binwidth=2, position="dodge")+
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme

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
df1B<- summarySE(allData.vesselType[allData.vesselType$BranchOrder>4,], measurevar="pO2", groupvars=c("Genotype", "O2_type","BranchOrder"))
df1C<- summarySE(allData.vesselType, measurevar="pO2", groupvars=c("Genotype","BranchGroup"))
df1D<- summarySE(allData.vesselType[allData.vesselType$O2_type=="artery",], measurevar="pO2", groupvars=c("Genotype","BranchGroup"))

allData.vesselType$AdjustedDepth<-as.factor(allData.vesselType$AdjustedDepth)
df1D<- summarySE(allData.vesselType, measurevar="pO2", groupvars=c("Genotype","BranchGroup","AdjustedDepth"))

ggplot(data=df1A, aes(x=Genotype, y=pO2, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=pO2-se, ymax=pO2+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("pO2 [mmHg]") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot() +
  geom_jitter(data=allData.vesselType, aes(y=pO2, x=Genotype, colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df1A, aes(x=Genotype, y=pO2,ymin=pO2, ymax=pO2, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df1A, aes(x=Genotype, ymin=pO2-se, ymax=pO2+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("pO2 all vessels")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot() +
  geom_bar(data=df1A, aes(x=Genotype, y=pO2, colour=Genotype), fill="white" ,stat="identity", width=0.5, position=position_dodge()) +
  geom_jitter(data=allData.vesselType, aes(y=pO2, x=Genotype, colour=Genotype), 
              position=position_jitter(width=0.12), size=4)+
  geom_errorbar(data=df1A, aes(x=Genotype, ymin=pO2-se, ymax=pO2+se), size=2,colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("pO2 all vessels")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

ggplot(data=df1B, aes(x=interaction(O2_type,BranchOrder), y=pO2, fill=Genotype)) +
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

ggplot() +
  geom_jitter(data=allData.vesselType, aes(y=pO2, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df1C, aes(x=interaction(Genotype,BranchGroup), y=pO2,ymin=pO2, ymax=pO2, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df1C, aes(x=interaction(Genotype,BranchGroup), ymin=pO2-se, ymax=pO2+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("pO2 all vessels")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
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


#####
library(zoo)
# depth vs pO2

ggplot(allData.vesselType, aes(x=AdjustedDepth, y=pO2, colour=interaction(O2_type,Genotype)))+
  geom_jitter(position=position_jitter(width=0.12), size=4)+
  geom_line(aes(y=rollmean(pO2,50, na.pad = TRUE)))+
  ggtitle("pO2 vs depth by vessel type") +
  xlab("depth") + 
  ylab("pO2") + 
  #scale_colour_manual(values=c("blue", "red")) + 
  max.theme

#plots like Chris Schaffer's paper from 2012

#adjust depth
allData.vesselType$AdjustedDepth<-as.numeric(as.character(allData.vesselType$AdjustedDepth))
arteries<-subset(allData.vesselType, O2_type=="artery")
veins<-subset(allData.vesselType, O2_type=="vein")

veins$AdjustedDepth=-(veins$AdjustedDepth)
veins$BranchOrder=-(veins$BranchOrder)
veins$AdjustedpO2=veins$pO2
arteries$AdjustedpO2=-(arteries$pO2)

allVessels.adjusted<-rbind(arteries,veins)

# plots

#pO2 vs branch order
ggplot(data=allVessels.adjusted, aes(x=BranchOrder, y=pO2, colour=Genotype)) +
  geom_jitter(position=position_jitter(width=0.12), size=3, shape=1)+
  #geom_line(aes(y=rollmean(pO2,7, na.pad = TRUE)))+
  facet_grid(. ~O2_type,scales = "free" , space="free")+
  ggtitle("pO2 vs branch order by vessel type") +
  xlab("Branch order") + 
  ylab("pO2") + 
  scale_colour_manual(values=c("blue", "red")) + 
  max.theme

#pO2 vs depth
ggplot(data=allVessels.adjusted, aes(x=AdjustedDepth, y=pO2, colour=Genotype)) +
  geom_jitter(position=position_jitter(width=2), size=3, shape=1)+
  #geom_line(aes(y=rollmean(pO2,7, na.pad = TRUE)))+
  facet_grid(. ~O2_type,scales = "free" , space="free")+
  ggtitle("pO2 vs depth by vessel type") +
  xlab("depth") + 
  ylab("pO2") + 
  scale_colour_manual(values=c("blue", "red")) + 
  max.theme


# calculate the mean at each depth

df2.depth<-summarySE(allVessels.adjusted, measurevar="pO2", groupvars=c("Genotype","AdjustedDepth","O2_type"))

ggplot(df2.depth, aes(x=AdjustedDepth, y=pO2, colour=Genotype))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=pO2-se, ymax=pO2+se), colour="black", width=1,  position=position_dodge()) +
  geom_line()+
  facet_grid(. ~O2_type,scales = "free" , space="free")+
  ggtitle("mean pO2 per depth by vessel type") +
  xlab("depth") + 
  ylab("mean pO2") + 
  scale_colour_manual(values=c("blue", "red")) + 
  max.theme
  

######
#Stats
pO2.null = lmer(pO2 ~ (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType,REML=FALSE)
pO2.model1B = lmer(pO2~ Genotype + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType,REML=FALSE)
pO2.model2B = lmer(pO2~ O2_type + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType,REML=FALSE)
pO2.model2C = lmer(pO2~ BranchGroup + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType,REML=FALSE)
pO2.model3B = lmer(pO2~ Genotype + O2_type + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType,REML=FALSE)
pO2.model3C = lmer(pO2~ Genotype + BranchGroup + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType,REML=FALSE)
pO2.model4B = lmer(pO2~ Genotype * O2_type + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType,REML=FALSE)
pO2.model4C = lmer(pO2~ Genotype * BranchGroup + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType,REML=FALSE)
pO2.anova2 <- anova(pO2.null, pO2.model1B, pO2.model2B,pO2.model2C,
                   pO2.model3B,pO2.model3C,pO2.model4B,pO2.model4C)
print(pO2.anova2)

# p values
pO2.genotype <- lsmeans(pO2.model1B, pairwise ~ Genotype, glhargs=list())
summary(pO2.genotype)

pO2.O2type <- lsmeans(pO2.model4B, pairwise ~ Genotype*O2_type, glhargs=list())
summary(pO2.O2type)

pO2.BranchGroup <- lsmeans(pO2.model4C, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(pO2.BranchGroup)

# Look at the model residuals
plot(pO2.model1B)
plot(pO2.model3B)
plot(pO2.model3C)
plot(pO2.model4B)
plot(pO2.model4C)

# pO2, genotype and depth
allData.vesselType$AdjustedDepth<-as.factor(allData.vesselType$AdjustedDepth)
pO2.nullB = lmer(pO2 ~ (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType,REML=FALSE)
pO2.model1D = lmer(pO2~ Genotype + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType,REML=FALSE)
pO2.model2D = lmer(pO2~ AdjustedDepth + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType,REML=FALSE)
pO2.model3D = lmer(pO2~ Genotype + AdjustedDepth + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType,REML=FALSE)
pO2.model3E = lmer(pO2~ Genotype + AdjustedDepth + O2_type + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType,REML=FALSE)
pO2.model4D = lmer(pO2~ Genotype * AdjustedDepth + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType,REML=FALSE)
pO2.model4E = lmer(pO2~ Genotype * AdjustedDepth*O2_type + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType,REML=FALSE)
pO2.anova.depth <- anova(pO2.nullB, pO2.model1D, pO2.model2D,pO2.model3D,pO2.model3E,
                         pO2.model4D,pO2.model4E)
print(pO2.anova.depth)

pO2.depth <- lsmeans(pO2.model4D, pairwise ~ Genotype*AdjustedDepth, glhargs=list())
summary(pO2.depth)

pO2.depth.type <- lsmeans(pO2.model4E, pairwise ~ Genotype*AdjustedDepth*O2_type, glhargs=list())
summary(pO2.depth.type)


#analyze arteries and veins separately
# arteries
pO2.null.art = lmer(pO2 ~ (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType[allData.vesselType$O2_type=="artery",],REML=FALSE)
pO2.model1.art = lmer(pO2~ Genotype + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType[allData.vesselType$O2_type=="artery",],REML=FALSE)
pO2.model2.art = lmer(pO2~ AdjustedDepth + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType[allData.vesselType$O2_type=="artery",],REML=FALSE)
pO2.model3.art = lmer(pO2~ Genotype + AdjustedDepth + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType[allData.vesselType$O2_type=="artery",],REML=FALSE)
pO2.model4.art = lmer(pO2~ Genotype * AdjustedDepth + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType[allData.vesselType$O2_type=="artery",],REML=FALSE)
pO2.anova.depth.art <- anova(pO2.null.art, pO2.model1.art, pO2.model2.art,pO2.model3.art,pO2.model4.art)
print(pO2.anova.depth.art)

pO2.depth.art <- lsmeans(pO2.model4.art, pairwise ~ Genotype*AdjustedDepth, glhargs=list())
summary(pO2.depth.art)

# veins
pO2.null.vein = lmer(pO2 ~ (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType[allData.vesselType$O2_type=="vein",],REML=FALSE)
pO2.model1.vein = lmer(pO2~ Genotype + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType[allData.vesselType$O2_type=="vein",],REML=FALSE)
pO2.model2.vein = lmer(pO2~ AdjustedDepth + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType[allData.vesselType$O2_type=="vein",],REML=FALSE)
pO2.model3.vein = lmer(pO2~ Genotype + AdjustedDepth + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType[allData.vesselType$O2_type=="vein",],REML=FALSE)
pO2.model4.vein = lmer(pO2~ Genotype * AdjustedDepth + (1|Animal) + (1|Spot) + (1|BranchName), allData.vesselType[allData.vesselType$O2_type=="vein",],REML=FALSE)
pO2.anova.depth.vein <- anova(pO2.null.vein, pO2.model1.vein, pO2.model2.vein,pO2.model3.vein,pO2.model4.vein)
print(pO2.anova.depth.vein)

pO2.depth.vein <- lsmeans(pO2.model4.vein, pairwise ~ Genotype*AdjustedDepth, glhargs=list())
summary(pO2.depth.vein)

##############
