
library("lme4")
library("lmerTest")
library("lattice")
library("plyr")
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
baseline1 <- read.delim("E:/Data/Pericyte_project/Two-photon-data/Ret_ret_Mice/Baseline_BloodFlow/Results/Diam_velocity_03_11_2017.csv", header=TRUE, sep = "\t")
baseline1 <- read.delim("D:/Data/Pericytes/Results/Diam_velocity_27_10_2017.csv", header=TRUE, sep = "\t")


#########
# unique animal and Branchname name
baseline1$Vesselname <-paste(baseline1$AnimalName, baseline1$Spot,sep= "_")
baseline1$Branchname <-paste(baseline1$AnimalName, baseline1$Branch,sep= "_")

baseline1$Genotype<- factor(baseline1$Genotype,levels = c("Ret+", "RetRet"))
#baseline$BranchOrder<- as.factor(baseline$BranchOrder)


# REMOVE DATA !!! from penetrating arterioles, ascending venules or surface vessels
baseline1<-baseline1[!baseline1$BranchOrder<=0,]

baseline1$BranchGroup<-"0"

eGFPpos<-subset(baseline1, eGFP==1)
eGFPneg<-subset(baseline1, eGFP==0)

#eGFPpos$BranchGroup[eGFPpos$BranchOrder<=0]<-"arteriole"
eGFPpos$BranchGroup[eGFPpos$BranchOrder<=4]<-"ensheathing_PC"
eGFPpos$BranchGroup[eGFPpos$BranchOrder>4]<-"capillary_PC_A"
eGFPneg$BranchGroup[eGFPneg$BranchOrder>3]<-"capillary_PC_V"
eGFPneg$BranchGroup[eGFPneg$BranchOrder<=3]<-"venule_PC"

eGFPpos$BranchGroup<- factor(eGFPpos$BranchGroup,levels = c("ensheathing_PC", "capillary_PC_A"))

baseline<-rbind(eGFPpos,eGFPneg)
#baseline<-eGFPpos

baseline$BranchGroup<- as.factor(baseline$BranchGroup)
baseline$BranchGroup<- factor(baseline$BranchGroup,levels = c("ensheathing_PC", "capillary_PC_A",
                                                              "capillary_PC_V","venule_PC"))

baseline.withNoFlow<-baseline
baseline<-baseline[!baseline$Velocity>10,]
baseline<-baseline[complete.cases(baseline$Velocity),]

# only consider data from vessels were we have all the info (flux, etc.)
#baseline<-baseline[complete.cases(baseline$Flux),]



# CONSIDER each Branch

# velocity for each branch at different branch orders
ggplot(eGFPpos[eGFPpos$Genotype=="Ret+",], aes(x=BranchOrder, y=Velocity)) +
  geom_point(aes(colour = Branchname, fill=Branchname), size=2)+
  geom_line(aes(colour = Branchname, fill=Branchname), size=1)+
  ggtitle("Ret+-Velocity vs order for all vessels") +
  xlab("BranchOrder") + 
  ylab("Velocity [mm/s]") + 
  ylim(0,6)+
  #scale_colour_manual(values=c("blue", "red"), guide=FALSE) + 
  max.theme

ggplot(eGFPpos[eGFPpos$Genotype=="RetRet",], aes(x=BranchOrder, y=Velocity)) +
  geom_point(aes(colour = Branchname, fill=Branchname), size=2)+
  geom_line(aes(colour = Branchname, fill=Branchname), size=1)+
  ggtitle("RetRet-Velocity vs order for all vessels") +
  xlab("BranchOrder") + 
  ylab("Velocity [mm/s]") + 
  ylim(0,6)+
  #scale_colour_manual(values=c("blue", "red"), guide=FALSE) + 
  max.theme

###############################
#histograms

ggplot(baseline, aes(x=Diameter, y=..density..,fill=Genotype)) +
  geom_histogram(binwidth=1, position="dodge")+
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme

ggplot(baseline, aes(x=Velocity, y=..density..,fill=Genotype)) +
  geom_histogram(binwidth=0.25, position="dodge")+
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme

ggplot(baseline, aes(x=Flux, y=..density..,fill=Genotype)) +
  geom_histogram(binwidth=15, position="dodge")+
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme



ggplot(baseline, aes(x=PulsatilityIndex, fill=Genotype)) + geom_histogram(binwidth=0.1, position="dodge") +
  ggtitle("Distribution of Pulsatility")



#ggpairs(baseline, columns = 7:17, aes(colour=Genotype), alpha=0.4)


#########
# diameter
df1A<- summarySE(baseline, measurevar="Diameter", groupvars=c("Genotype"))
df1B<- summarySE(eGFPpos, measurevar="Diameter", groupvars=c("Genotype","BranchOrder"))
df1C<- summarySE(baseline, measurevar="Diameter", groupvars=c("Genotype","BranchGroup"))
df1D<- summarySE(eGFPpos, measurevar="Diameter", groupvars=c("Genotype","BranchGroup"))

ggplot(data=df1A, aes(x=Genotype, y=Diameter, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Diameter-se, ymax=Diameter+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Diameter [um]") +
  scale_fill_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

#scatterplot of genotypes- all Data
ggplot() +
  geom_jitter(data=baseline, aes(y=Diameter, x=Genotype, colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df1A, aes(x=Genotype, y=Diameter,ymin=Diameter, ymax=Diameter, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df1A, aes(x=Genotype, ymin=Diameter-se, ymax=Diameter+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Diameter, all Data")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot() +
  geom_bar(data=df1A, aes(x=Genotype, y=Diameter, colour=Genotype), fill="white" ,stat="identity", width=0.5, position=position_dodge()) +
  geom_jitter(data=baseline, aes(y=Diameter, x=Genotype, colour=Genotype), 
              position=position_jitter(width=0.12), size=4)+
  geom_errorbar(data=df1A, aes(x=Genotype, ymin=Diameter-se, ymax=Diameter+se), size=2,colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Diameter all vessels")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

ggplot() +
  geom_bar(data=df1C, aes(x=interaction(Genotype,BranchGroup), y=Diameter, colour=Genotype), fill="white" ,stat="identity", width=0.5, position=position_dodge()) +
  geom_jitter(data=baseline, aes(y=Diameter, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4)+
  geom_errorbar(data=df1C, aes(x=interaction(Genotype,BranchGroup), ymin=Diameter-se, ymax=Diameter+se), size=2,colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Diameter branch groups")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot(data=df1B, aes(x=BranchOrder, y=Diameter, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Diameter-se, ymax=Diameter+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Diameter [um]") +
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme



ggplot(data=df1C, aes(x=BranchGroup, y=Diameter, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Diameter-se, ymax=Diameter+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Diameter [um]") +
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme

#scatterplot of genotypes- all Data
ggplot() +
  geom_jitter(data=baseline, aes(y=Diameter, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df1C, aes(x=interaction(Genotype,BranchGroup), y=Diameter,ymin=Diameter, ymax=Diameter, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df1C, aes(x=interaction(Genotype,BranchGroup), ymin=Diameter-se, ymax=Diameter+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Diameter, all Data")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme




#scatterplot of genotypes- arteriole
ggplot() +
  geom_jitter(data=eGFPpos, aes(y=Diameter, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df1D, aes(x=interaction(Genotype,BranchGroup), y=Diameter,ymin=Diameter, ymax=Diameter, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df1D, aes(x=interaction(Genotype,BranchGroup), ymin=Diameter-se, ymax=Diameter+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Diameter, arteriole")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

## boxplots
ggplot(baseline, aes(x = Genotype, y = Diameter, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Diameter [um]") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchOrder), y = Diameter, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Diameter [um]") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchGroup), y = Diameter, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Diameter [um]") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme




######
#Stats
## Diameter and genotype or branch order
diam.all.null = lmer(Diameter ~ (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
diam.all.model1 = lmer(Diameter~ Genotype + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
diam.all.model2 = lmer(Diameter~ BranchGroup + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
diam.all.model3 = lmer(Diameter~ Genotype + BranchGroup + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
diam.all.model4 = lmer(Diameter~ Genotype * BranchGroup + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
diam.all.anova <- anova(diam.all.null, diam.all.model1,diam.all.model2,diam.all.model3,diam.all.model4)
print(diam.all.anova)
# p values
diam.all.Genotype <- lsmeans(diam.all.model1, pairwise ~ Genotype, glhargs=list())
summary(diam.all.Genotype)

diam.all.Genotype_BO <- lsmeans(diam.all.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(diam.all.Genotype_BO)


# Look at the per-animal difference
dotplot(ranef(diam.all.null, condVar = TRUE))$AnimalName
dotplot(ranef(diam.all.model1, condVar = TRUE))$AnimalName
dotplot(ranef(diam.all.model2, condVar = TRUE))$AnimalName
dotplot(ranef(diam.all.model3, condVar = TRUE))$AnimalName
dotplot(ranef(diam.all.model4, condVar = TRUE))$AnimalName

# Look at the model residuals
plot(diam.all.model1)
plot(diam.all.model2)
plot(diam.all.model3)


# arteriole side only
diam.art.null = lmer(Diameter ~ (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
diam.art.model1 = lmer(Diameter~ Genotype + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
diam.art.model2 = lmer(Diameter~ BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
diam.art.model3 = lmer(Diameter~ Genotype + BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
diam.art.model4 = lmer(Diameter~ Genotype * BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
diam.art.anova <- anova(diam.art.null, diam.art.model1,diam.art.model2,diam.art.model3,diam.art.model4)
print(diam.art.anova)
# p values
diam.art.Genotype <- lsmeans(diam.art.model1, pairwise ~ Genotype, glhargs=list())
summary(diam.art.Genotype)

diam.art.Genotype_BO <- lsmeans(diam.art.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(diam.art.Genotype_BO)


# diameter, genotype and depth
diam.nullB = lmer(Diameter ~ (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
diam.model1B = lmer(Diameter~ Genotype + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
diam.model2B = lmer(Diameter~ Depth + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
diam.model3B = lmer(Diameter~ Genotype + Depth + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
diam.model4B = lmer(Diameter~ Genotype * Depth + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
diam.anovaB <- anova(diam.nullB, diam.model1B,diam.model2B,diam.model3B,diam.model4B)
print(diam.anovaB)
summary(diam.model3B)


# no interaction between genotype and depth or diameter


##############
# depth 

#depth plots

# only arteriole Side
#eGFPpos$BranchOrder<-as.factor(eGFPpos$BranchOrder)
eGFPpos$BranchGroup<- factor(eGFPpos$BranchGroup,levels = c("ensheathing_PC", "capillary_PC_A"))
#eGFPneg$BranchOrder<-as.factor(eGFPneg$BranchOrder)
eGFPneg$BranchGroup<- factor(eGFPneg$BranchGroup,levels = c("capillary_PC_V","venule_PC"))

ggplot(eGFPpos, aes(x = BranchOrder, y = Depth, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("depth [um]") +
  ggtitle("arteriole side")+
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(eGFPneg, aes(x = BranchOrder, y = Depth, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("depth [um]") +
  ggtitle("venous side")+
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

# scatterplot- diameter vs depth
ggplot(eGFPpos, aes(x=Diameter, y=Depth)) +
  geom_point(aes(colour = Genotype), shape = 1, size=2)+
  ggtitle("diameter vs depth for all vessels") +
  xlab("Diameter [um]") + 
  ylab("Depth [um]") + 
  scale_colour_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

# depth per branch order

ggplot(eGFPpos, aes(x=BranchOrder, y=Depth)) +
  geom_jitter(aes(colour = Genotype), shape = 1, size=2)+
  ggtitle("branch order vs depth for arteriole side") +
  xlab("Branch Order") + 
  ylab("Depth [um]") + 
  scale_colour_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

df.depth.art<- summarySE(eGFPpos, measurevar="Depth", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)
df.depth.ven<- summarySE(eGFPneg, measurevar="Depth", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)

#scatterplot of genotypes- all Data
ggplot() +
  geom_jitter(data=eGFPpos, aes(y=Depth, x=interaction(Genotype, BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df.depth.art, aes(x=interaction(Genotype, BranchGroup), y=Depth,ymin=Depth, ymax=Depth, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df.depth.art, aes(x=interaction(Genotype, BranchGroup), ymin=Depth-se, ymax=Depth+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle(" Depth for each branch group arteriole side")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot() +
  geom_bar(data=df.depth.art, aes(x=interaction(Genotype,BranchGroup), y=Depth, colour=Genotype), fill="white" ,stat="identity", width=0.5, position=position_dodge()) +
  geom_jitter(data=eGFPpos, aes(y=Depth, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4)+
  geom_errorbar(data=df.depth.art, aes(x=interaction(Genotype,BranchGroup), ymin=Depth-se, ymax=Depth+se), size=2,colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Diameter depth arteriole")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot(eGFPpos, aes(x = BranchGroup, y = Depth, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("depth [um]") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme



ggplot() +
  geom_jitter(data=eGFPneg, aes(y=Depth, x=interaction(Genotype, BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df.depth.ven, aes(x=interaction(Genotype, BranchGroup), y=Depth,ymin=Depth, ymax=Depth, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df.depth.ven, aes(x=interaction(Genotype, BranchGroup), ymin=Depth-se, ymax=Depth+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle(" Depth for each branch group venule side")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot(eGFPneg, aes(x = BranchGroup, y = Depth, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("depth [um]") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

# depth and branch group
#arteriole side
Adepth.null = lmer(Depth ~ (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
Adepth.model1 = lmer(Depth~ BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
Adepth.model2 = lmer(Depth~ Genotype + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
Adepth.model3 = lmer(Depth~ Genotype + BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
Adepth.model4 = lmer(Depth~ Genotype * BranchGroup+ (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
Adepth.anova <- anova(Adepth.null, Adepth.model1, Adepth.model2, Adepth.model3, Adepth.model4)
print(Adepth.anova)

Adepth.Genotype_BG <- lsmeans(Adepth.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(Adepth.Genotype_BG)

# ALMOST SIGNIFICANT

#venous side
Vdepth.null = lmer(Depth ~ (1|AnimalName) + (1|Branchname), eGFPneg,REML=FALSE)
Vdepth.model1 = lmer(Depth~ BranchGroup + (1|AnimalName) + (1|Branchname), eGFPneg,REML=FALSE)
Vdepth.model2 = lmer(Depth~ Genotype + (1|AnimalName) + (1|Branchname), eGFPneg,REML=FALSE)
Vdepth.model3 = lmer(Depth~ Genotype + BranchGroup + (1|AnimalName) + (1|Branchname), eGFPneg,REML=FALSE)
Vdepth.model4 = lmer(Depth~ Genotype * BranchGroup+ (1|AnimalName) + (1|Branchname), eGFPneg,REML=FALSE)
Vdepth.anova <- anova(Vdepth.null, Vdepth.model1, Vdepth.model2, Vdepth.model3, Vdepth.model4)
print(Vdepth.anova)

Vdepth.Genotype_BG <- lsmeans(Vdepth.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(Vdepth.Genotype_BG)


#########
# velocity

#all data, including those where there is no flux measurements

df2A<- summarySE(baseline, measurevar="Velocity", groupvars=c("Genotype"), na.rm=TRUE)
df2A2<- summarySE(eGFPpos, measurevar="Velocity", groupvars=c("Genotype"), na.rm=TRUE)
df2B<- summarySE(eGFPpos, measurevar="Velocity", groupvars=c("Genotype","BranchOrder"), na.rm=TRUE)
df2C<- summarySE(baseline, measurevar="Velocity", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)
df2D<- summarySE(eGFPpos, measurevar="Velocity", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)

ggplot(data=df2A, aes(x=Genotype, y=Velocity, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Velocity-se, ymax=Velocity+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Velocity [mm/s]") +
  scale_fill_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

#scatterplot of genotypes- all Data
ggplot() +
  geom_jitter(data=baseline, aes(y=Velocity, x=Genotype, colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df2A, aes(x=Genotype, y=Velocity,ymin=Velocity, ymax=Velocity, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df2A, aes(x=Genotype, ymin=Velocity-se, ymax=Velocity+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle(" velocities, all Data")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

#scatterplot of genotypes- arteriole
ggplot() +
  geom_jitter(data=eGFPpos, aes(y=Velocity, x=Genotype, colour=Genotype), position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df2A2, aes(x=Genotype, y=Velocity,ymin=Velocity, ymax=Velocity, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df2A2, aes(x=Genotype, ymin=Velocity-se, ymax=Velocity+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle(" velocities, arteriole")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot() +
  geom_bar(data=df2A, aes(x=Genotype, y=Velocity, colour=Genotype), fill="white" ,stat="identity", width=0.5, position=position_dodge()) +
  geom_jitter(data=baseline, aes(y=Velocity, x=Genotype, colour=Genotype), 
              position=position_jitter(width=0.12), size=4)+
  geom_errorbar(data=df2A, aes(x=Genotype, ymin=Velocity-se, ymax=Velocity+se), size=2,colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Velocity all vessels")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

ggplot() +
  geom_bar(data=df2C, aes(x=interaction(Genotype,BranchGroup), y=Velocity, colour=Genotype), fill="white" ,stat="identity", width=0.5, position=position_dodge()) +
  geom_jitter(data=baseline, aes(y=Velocity, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4)+
  geom_errorbar(data=df2C, aes(x=interaction(Genotype,BranchGroup), ymin=Velocity-se, ymax=Velocity+se), size=2,colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Velocity branch groups")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot(data=df2B, aes(x=BranchOrder, y=Velocity, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Velocity-se, ymax=Velocity+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Velocity [mm/s]") +
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme


ggplot(data=df2C, aes(x=BranchGroup, y=Velocity, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Velocity-se, ymax=Velocity+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Velocity [mm/s]") +
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme

ggplot() +
  geom_jitter(data=baseline, aes(y=Velocity, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df2C, aes(x=interaction(Genotype,BranchGroup), y=Velocity,ymin=Velocity, ymax=Velocity, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df2C, aes(x=interaction(Genotype,BranchGroup), ymin=Velocity-se, ymax=Velocity+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle(" velocities, all Data")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot(data=df2D, aes(x=BranchGroup, y=Velocity, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Velocity-se, ymax=Velocity+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Velocity [mm/s]") +
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme

ggplot() +
  geom_jitter(data=eGFPpos, aes(y=Velocity, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df2D, aes(x=interaction(Genotype,BranchGroup), y=Velocity,ymin=Velocity, ymax=Velocity, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df2D, aes(x=interaction(Genotype,BranchGroup), ymin=Velocity-se, ymax=Velocity+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle(" velocities, arteriole side")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

ggplot() +
  geom_bar(data=df2D, aes(x=interaction(Genotype,BranchGroup), y=Velocity, colour=Genotype), fill="white" ,stat="identity", width=0.5, position=position_dodge()) +
  geom_jitter(data=eGFPpos, aes(y=Velocity, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4)+
  geom_errorbar(data=df2D, aes(x=interaction(Genotype,BranchGroup), ymin=Velocity-se, ymax=Velocity+se), size=2,colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Velocity branch groups")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

###########
## boxplots
ggplot(baseline, aes(x = Genotype, y = Velocity, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Velocity [mm/s]") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(eGFPpos, aes(x = BranchOrder, y = Velocity, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Velocity [mm/s]") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = BranchGroup, y = Velocity, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Velocity [mm/s]") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(eGFPpos, aes(x = BranchGroup, y = Velocity, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Velocity [mm/s]") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme



#depth plots
# scatterplot- Velocity vs depth
ggplot(baseline, aes(x=Velocity, y=Depth)) +
  geom_point(aes(colour = Genotype), shape = 1, size=4)+
  ggtitle("Velocity vs depth for all vessels") +
  xlab("Velocity [mm/s]") + 
  ylab("Depth [mm/s]") + 
  scale_colour_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

# scatterplot- Diameter vs depth
ggplot(baseline, aes(x=Diameter, y=Depth)) +
  geom_point(aes(colour = Genotype), shape = 1, size=4)+
  ggtitle("Diameter vs depth for all vessels") +
  xlab("Diameter]") + 
  ylab("Depth [mm/s]") + 
  scale_colour_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme


######
#Stats
## Velocity and genotype or branch group
vel.all.null = lmer(Velocity ~ (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
vel.all.model1 = lmer(Velocity~ Genotype + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
vel.all.model2 = lmer(Velocity~ BranchGroup + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
vel.all.model3 = lmer(Velocity~ Genotype + BranchGroup + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
vel.all.model4 = lmer(Velocity~ Genotype * BranchGroup + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
vel.all.anova <- anova(vel.all.null, vel.all.model1,vel.all.model2,vel.all.model3,vel.all.model4)
print(vel.all.anova)

# p values
vel.all.Genotype <- lsmeans(vel.all.model1, pairwise ~ Genotype, glhargs=list())
summary(vel.all.Genotype)

vel.all.Genotype_BG <- lsmeans(vel.all.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(vel.all.Genotype_BG)


# Look at the per-animal difference
dotplot(ranef(vel.all.null, condVar = TRUE))$AnimalName
dotplot(ranef(vel.all.model1, condVar = TRUE))$AnimalName
dotplot(ranef(vel.all.model2, condVar = TRUE))$AnimalName
dotplot(ranef(vel.all.model3, condVar = TRUE))$AnimalName
dotplot(ranef(vel.all.model4, condVar = TRUE))$AnimalName

# Look at the model residuals
plot(vel.all.model1)
plot(vel.all.model2)
plot(vel.all.model3)


#arteriole side only

## Velocity and genotype or branch group
vel.art.null = lmer(Velocity ~ (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
vel.art.model1 = lmer(Velocity~ Genotype + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
vel.art.model2 = lmer(Velocity~ BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
vel.art.model3 = lmer(Velocity~ Genotype + BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
vel.art.model4 = lmer(Velocity~ Genotype * BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
vel.art.anova <- anova(vel.art.null, vel.art.model1,vel.art.model2,vel.art.model3,vel.art.model4)
print(vel.art.anova)

# p values
vel.art.Genotype <- lsmeans(vel.art.model1, pairwise ~ Genotype, glhargs=list())
summary(vel.art.Genotype)

vel.art.Genotype_BG <- lsmeans(vel.art.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(vel.art.Genotype_BG)


# Look at the per-animal difference
dotplot(ranef(vel.art.null, condVar = TRUE))$AnimalName
dotplot(ranef(vel.art.model1, condVar = TRUE))$AnimalName
dotplot(ranef(vel.art.model2, condVar = TRUE))$AnimalName
dotplot(ranef(vel.art.model3, condVar = TRUE))$AnimalName
dotplot(ranef(vel.art.model4, condVar = TRUE))$AnimalName

# Look at the model residuals
plot(vel.art.model1)
plot(vel.art.model2)
plot(vel.art.model3)


# Velocity, genotype and depth
vel.nullB = lmer(Velocity ~ (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
vel.model1B = lmer(Velocity~ Genotype + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
vel.model2B = lmer(Velocity~ Depth + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
vel.model3B = lmer(Velocity~ Genotype + Depth + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
vel.model4B = lmer(Velocity~ Genotype * Depth + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
vel.anovaB <- anova(vel.nullB, vel.model1B,vel.model2B,vel.model3B,vel.model4B)
print(vel.anovaB)


#########


# Pulsatility Index

df6A<- summarySE(baseline, measurevar="PulsatilityIndex", groupvars=c("Genotype"), na.rm=TRUE)
df6A2<- summarySE(eGFPpos, measurevar="PulsatilityIndex", groupvars=c("Genotype"), na.rm=TRUE)
df6B<- summarySE(eGFPpos, measurevar="PulsatilityIndex", groupvars=c("Genotype","BranchOrder"), na.rm=TRUE)
df6C<- summarySE(baseline, measurevar="PulsatilityIndex", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)
df6D<- summarySE(eGFPpos, measurevar="PulsatilityIndex", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)



ggplot(data=df6A, aes(x=Genotype, y=PulsatilityIndex, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=PulsatilityIndex-se, ymax=PulsatilityIndex+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("PulsatilityIndex") +
  scale_fill_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


#scatterplot of genotypes- all Data
ggplot() +
  geom_jitter(data=baseline, aes(y=PulsatilityIndex, x=Genotype, colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df6A, aes(x=Genotype, y=PulsatilityIndex,ymin=PulsatilityIndex, ymax=PulsatilityIndex, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df6A, aes(x=Genotype, ymin=PulsatilityIndex-se, ymax=PulsatilityIndex+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("PulsatilityIndex, all Data")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

ggplot() +
  geom_bar(data=df6A, aes(x=Genotype, y=PulsatilityIndex, colour=Genotype), fill="white" ,stat="identity", width=0.5, position=position_dodge()) +
  geom_jitter(data=baseline, aes(y=PulsatilityIndex, x=Genotype, colour=Genotype), 
              position=position_jitter(width=0.12), size=4)+
  geom_errorbar(data=df6A, aes(x=Genotype, ymin=PulsatilityIndex-se, ymax=PulsatilityIndex+se), size=2,colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("PulsatilityIndex all vessels")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme



ggplot(data=df6B, aes(x=BranchOrder, y=PulsatilityIndex, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=PulsatilityIndex-se, ymax=PulsatilityIndex+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("PulsatilityIndex") +
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme


ggplot(data=df6C, aes(x=BranchGroup, y=PulsatilityIndex, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=PulsatilityIndex-se, ymax=PulsatilityIndex+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("PulsatilityIndex") +
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme


#scatterplot of genotypes- all Data
ggplot() +
  geom_jitter(data=baseline, aes(y=PulsatilityIndex, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df6C, aes(x=interaction(Genotype,BranchGroup), y=PulsatilityIndex,ymin=PulsatilityIndex, ymax=PulsatilityIndex, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df6C, aes(x=interaction(Genotype,BranchGroup), ymin=PulsatilityIndex-se, ymax=PulsatilityIndex+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("PulsatilityIndex, all Data")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

#scatterplot of genotypes- all Data
ggplot() +
  geom_jitter(data=eGFPpos, aes(y=PulsatilityIndex, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df6D, aes(x=interaction(Genotype,BranchGroup), y=PulsatilityIndex,ymin=PulsatilityIndex, ymax=PulsatilityIndex, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df6D, aes(x=interaction(Genotype,BranchGroup), ymin=PulsatilityIndex-se, ymax=PulsatilityIndex+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("PulsatilityIndex, arteriole")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot() +
  geom_bar(data=df6D, aes(x=interaction(Genotype,BranchGroup), y=PulsatilityIndex, colour=Genotype), fill="white" ,stat="identity", width=0.5, position=position_dodge()) +
  geom_jitter(data=eGFPpos, aes(y=PulsatilityIndex, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4)+
  geom_errorbar(data=df6D, aes(x=interaction(Genotype,BranchGroup), ymin=PulsatilityIndex-se, ymax=PulsatilityIndex+se), size=2,colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("PulsatilityIndex branch groups")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

######
#boxplots
ggplot(baseline, aes(x = Genotype, y = PulsatilityIndex, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("PulsatilityIndex") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchOrder), y = PulsatilityIndex, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("PulsatilityIndex") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchGroup), y = PulsatilityIndex, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("PulsatilityIndex") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

######
#depth plots
# scatterplot- PulsatilityIndex vs depth
ggplot(baseline, aes(x=PulsatilityIndex, y=Depth)) +
  geom_point(aes(colour = Genotype), shape = 1, size=2)+
  ggtitle("PulsatilityIndex vs depth for all vessels") +
  xlab("PulsatilityIndex") + 
  ylab("Depth [um]") + 
  scale_colour_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme


# pulsatility vs velocity
ggplot(baseline, aes(x=PulsatilityIndex, y=Velocity)) +
  geom_point(aes(colour = Genotype), shape = 1, size=3)+
  ggtitle("PulsatilityIndex vs velocity for all vessels") +
  xlab("PulsatilityIndex") + 
  ylab("velocity") + 
  scale_colour_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme


ggplot(baseline, aes(x=PulsatilityIndex, y=Diameter)) +
  geom_point(aes(colour = Genotype), shape = 1, size=3)+
  ggtitle("PulsatilityIndex vs velocity for all vessels") +
  xlab("PulsatilityIndex") + 
  ylab("diameter") + 
  scale_colour_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

##########
# slow vs fast vessels
medianVelocity=median(baseline$Velocity)

baseline$speed<-"fast"
baseline$speed[baseline$Velocity<medianVelocity]<-"slow"

df6E<- summarySE(baseline, measurevar="PulsatilityIndex", groupvars=c("Genotype","speed"), na.rm=TRUE)


ggplot() +
  geom_jitter(data=baseline, aes(y=PulsatilityIndex, x=interaction(Genotype,speed), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df6E, aes(x=interaction(Genotype,speed), y=PulsatilityIndex,ymin=PulsatilityIndex, ymax=PulsatilityIndex, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df6E, aes(x=interaction(Genotype,speed), ymin=PulsatilityIndex-se, ymax=PulsatilityIndex+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("PulsatilityIndex, speed")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

ggplot() +
  geom_bar(data=df6E, aes(x=interaction(Genotype,speed), y=PulsatilityIndex, colour=Genotype), fill="white" ,stat="identity", width=0.5, position=position_dodge()) +
  geom_jitter(data=baseline, aes(y=PulsatilityIndex, x=interaction(Genotype,speed), colour=Genotype), 
              position=position_jitter(width=0.12), size=4)+
  geom_errorbar(data=df6E, aes(x=interaction(Genotype,speed), ymin=PulsatilityIndex-se, ymax=PulsatilityIndex+se), size=2,colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("PulsatilityIndex branch groups")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

######
#Stats
## PulsatilityIndex and genotype or BranchGroup
PulsatilityIndex.all.null = lmer(PulsatilityIndex ~ (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.all.model1 = lmer(PulsatilityIndex~ Genotype + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.all.model2 = lmer(PulsatilityIndex~ BranchGroup + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.all.model3 = lmer(PulsatilityIndex~ Genotype + BranchGroup + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.all.model4 = lmer(PulsatilityIndex~ Genotype * BranchGroup + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.all.anova <- anova(PulsatilityIndex.all.null, PulsatilityIndex.all.model1,
                                    PulsatilityIndex.all.model2,PulsatilityIndex.all.model3,
                                    PulsatilityIndex.all.model4)
print(PulsatilityIndex.all.anova)
# p values
PulsatilityIndex.all.Genotype <- lsmeans(PulsatilityIndex.all.model1, pairwise ~ Genotype, glhargs=list())
summary(PulsatilityIndex.all.Genotype)

PulsatilityIndex.all.Genotype_BO <- lsmeans(PulsatilityIndex.all.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(PulsatilityIndex.all.Genotype_BO)


# Look at the per-animal difference
dotplot(ranef(PulsatilityIndex.all.null, condVar = TRUE))$AnimalName
dotplot(ranef(PulsatilityIndex.all.model1, condVar = TRUE))$AnimalName
dotplot(ranef(PulsatilityIndex.all.model2, condVar = TRUE))$AnimalName
dotplot(ranef(PulsatilityIndex.all.model3, condVar = TRUE))$AnimalName
dotplot(ranef(PulsatilityIndex.all.model4, condVar = TRUE))$AnimalName

# Look at the model residuals
plot(PulsatilityIndex.all.model1)
plot(PulsatilityIndex.all.model2)
plot(PulsatilityIndex.all.model3)
plot(PulsatilityIndex.all.model4)



## PulsatilityIndex and genotype or branch group
PulsatilityIndex.art.null = lmer(PulsatilityIndex ~ (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
PulsatilityIndex.art.model1 = lmer(PulsatilityIndex~ Genotype + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
PulsatilityIndex.art.model2 = lmer(PulsatilityIndex~ BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
PulsatilityIndex.art.model3 = lmer(PulsatilityIndex~ Genotype + BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
PulsatilityIndex.art.model4 = lmer(PulsatilityIndex~ Genotype * BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
PulsatilityIndex.art.anova <- anova(PulsatilityIndex.art.null, PulsatilityIndex.art.model1,PulsatilityIndex.art.model2,
                                    PulsatilityIndex.art.model3,PulsatilityIndex.art.model4)
print(PulsatilityIndex.art.anova)
# p values
PulsatilityIndex.art.Genotype <- lsmeans(PulsatilityIndex.art.model1, pairwise ~ Genotype, glhargs=list())
summary(PulsatilityIndex.art.Genotype)

PulsatilityIndex.art.Genotype_BO <- lsmeans(PulsatilityIndex.art.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(PulsatilityIndex.art.Genotype_BO)


# PulsatilityIndex, genotype and depth
PulsatilityIndex.nullB = lmer(PulsatilityIndex ~ (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.model1B = lmer(PulsatilityIndex~ Genotype + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.model2B = lmer(PulsatilityIndex~ Depth + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.model3B = lmer(PulsatilityIndex~ Genotype + Depth + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.model4B = lmer(PulsatilityIndex~ Genotype * Depth + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.anovaB <- anova(PulsatilityIndex.nullB, PulsatilityIndex.model1B,PulsatilityIndex.model2B,PulsatilityIndex.model3B,PulsatilityIndex.model4B)
print(PulsatilityIndex.anovaB)

# pulsatility and velocity
PulsatilityIndex.nullC = lmer(PulsatilityIndex ~ (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.model1C = lmer(PulsatilityIndex~ Genotype + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.model2C = lmer(PulsatilityIndex~ Velocity + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.model3C = lmer(PulsatilityIndex~ Genotype + Velocity  + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.model4C = lmer(PulsatilityIndex~ Genotype * Velocity  + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.anovaC <- anova(PulsatilityIndex.nullC, PulsatilityIndex.model1C,PulsatilityIndex.model2C,PulsatilityIndex.model3C,PulsatilityIndex.model4C)
print(PulsatilityIndex.anovaC)


# pulsatility and speed
PulsatilityIndex.speed.null = lmer(PulsatilityIndex ~ (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.speed.model1 = lmer(PulsatilityIndex~ Genotype + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.speed.model2 = lmer(PulsatilityIndex~ speed + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.speed.model3 = lmer(PulsatilityIndex~ Genotype + speed  + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.speed.model4 = lmer(PulsatilityIndex~ Genotype * speed  + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.speed.model5 = lmer(PulsatilityIndex~ Genotype + speed +BranchGroup  + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.speed.model6 = lmer(PulsatilityIndex~ Genotype * speed *BranchGroup  + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
PulsatilityIndex.speed.anova <- anova(PulsatilityIndex.speed.null, PulsatilityIndex.speed.model1,
                                      PulsatilityIndex.speed.model2,PulsatilityIndex.speed.model3,
                                      PulsatilityIndex.speed.model4,PulsatilityIndex.speed.model5,
                                      PulsatilityIndex.speed.model6)
print(PulsatilityIndex.speed.anova)

PulsatilityIndex.speed.Genotype<- lsmeans(PulsatilityIndex.speed.model4, pairwise ~ Genotype*speed, glhargs=list())
summary(PulsatilityIndex.speed.Genotype)

##########################
# vessel data with NO FLOW
baseline.withNoFlow$flow<-"yes"
baseline.withNoFlow$flow[baseline.withNoFlow$NoFlow==1]<-"no"

df.diam.noflow1<- summarySE(baseline.withNoFlow, measurevar="Diameter", groupvars=c("Genotype", "flow"), na.rm=TRUE)
df.diam.noflow3<- summarySE(baseline.withNoFlow, measurevar="Diameter", groupvars=c("Genotype","BranchGroup","flow"), na.rm=TRUE)

ggplot(data=df.diam.noflow1, aes(x=flow, y=Diameter, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Diameter-se, ymax=Diameter+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Diameter") +
  scale_fill_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

ggplot(data=df.diam.noflow3, aes(x=interaction(Genotype,BranchGroup), y=Diameter, fill=flow)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Diameter-se, ymax=Diameter+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Diameter") +
  scale_fill_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


# total non-flowing vessels
noFlow<- ddply(baseline.withNoFlow, c("Genotype","flow"), summarise, nVessels=length(Diameter))

# percentages
noFlow.retplus<-subset(noFlow, Genotype=="Ret+")
noFlow.retret<-subset(noFlow, Genotype=="RetRet")

total_Ret_plus=sum(noFlow$nVessels[noFlow$Genotype=="Ret+"])
total_Ret_Ret=sum(noFlow$nVessels[noFlow$Genotype=="RetRet"])

percent_Ret_plus=(noFlow.retplus$nVessels[noFlow.retplus$flow=="no"]/total_Ret_plus)*100
percent_Ret_Ret=(noFlow.retret$nVessels[noFlow.retret$flow=="no"]/total_Ret_Ret)*100


###################
# Flux

# only data that has 

df3A<- summarySE(baseline, measurevar="Flux", groupvars=c("Genotype"), na.rm=TRUE)
df3A2<- summarySE(eGFPpos, measurevar="Flux", groupvars=c("Genotype"), na.rm=TRUE)
df3B<- summarySE(eGFPpos, measurevar="Flux", groupvars=c("Genotype","BranchOrder"), na.rm=TRUE)
df3C<- summarySE(baseline, measurevar="Flux", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)
df3D<- summarySE(eGFPpos, measurevar="Flux", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)


ggplot(data=df3A, aes(x=Genotype, y=Flux, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Flux-se, ymax=Flux+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Flux [RBCs/s]") +
  scale_fill_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


#scatterplot of genotypes- all Data
ggplot() +
  geom_jitter(data=baseline, aes(y=Flux, x=Genotype, colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df3A, aes(x=Genotype, y=Flux,ymin=Flux, ymax=Flux, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df3A, aes(x=Genotype, ymin=Flux-se, ymax=Flux+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Flux, all Data")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot() +
  geom_bar(data=df3A, aes(x=Genotype, y=Flux, colour=Genotype), fill="white" ,stat="identity", width=0.5, position=position_dodge()) +
  geom_jitter(data=baseline, aes(y=Flux, x=Genotype, colour=Genotype), 
              position=position_jitter(width=0.12), size=4)+
  geom_errorbar(data=df3A, aes(x=Genotype, ymin=Flux-se, ymax=Flux+se), size=2,colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Flux all vessels")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

ggplot() +
  geom_bar(data=df3D, aes(x=interaction(Genotype,BranchGroup), y=Flux, colour=Genotype), fill="white" ,stat="identity", width=0.5, position=position_dodge()) +
  geom_jitter(data=eGFPpos, aes(y=Flux, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4)+
  geom_errorbar(data=df3D, aes(x=interaction(Genotype,BranchGroup), ymin=Flux-se, ymax=Flux+se), size=2,colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Flux branch groups")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot(data=df3B, aes(x=BranchOrder, y=Flux, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Flux-se, ymax=Flux+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Flux [RBCs/s]") +
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme


ggplot(data=df3C, aes(x=BranchGroup, y=Flux, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Flux-se, ymax=Flux+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Flux [RBCs/s]") +
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme

ggplot() +
  geom_jitter(data=baseline, aes(y=Flux, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df3C, aes(x=interaction(Genotype,BranchGroup), y=Flux,ymin=Flux, ymax=Flux, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df3C, aes(x=interaction(Genotype,BranchGroup), ymin=Flux-se, ymax=Flux+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Flux, all Data")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot(data=df3D, aes(x=BranchGroup, y=Flux, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Flux-se, ymax=Flux+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Flux [RBCs/s]") +
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme

ggplot() +
  geom_jitter(data=eGFPpos, aes(y=Flux, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df3D, aes(x=interaction(Genotype,BranchGroup), y=Flux,ymin=Flux, ymax=Flux, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df3D, aes(x=interaction(Genotype,BranchGroup), ymin=Flux-se, ymax=Flux+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Flux, arteriole")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

######
## boxplots
ggplot(baseline, aes(x = Genotype, y = Flux, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Flux [RBCs/s]") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchGroup), y = Flux, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Flux [RBCs/s]") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme


#depth plots
# scatterplot- Flux vs depth
ggplot(baseline, aes(x=Flux, y=Depth)) +
  geom_point(aes(colour = Genotype), shape = 1, size=2)+
  ggtitle("Flux vs depth for all vessels") +
  xlab("Flux [RBCs/s]") + 
  ylab("Depth [um]") + 
  scale_colour_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme


######
#Stats
## Flux and genotype or branch order
flux.all.null = lmer(Flux ~ (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
flux.all.model1 = lmer(Flux~ Genotype + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
flux.all.model2 = lmer(Flux~ BranchGroup + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
flux.all.model3 = lmer(Flux~ Genotype + BranchGroup + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
flux.all.model4 = lmer(Flux~ Genotype * BranchGroup + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
flux.all.anova <- anova(flux.all.null, flux.all.model1,flux.all.model2,flux.all.model3,flux.all.model4)
print(flux.all.anova)
# p values
flux.all.Genotype <- lsmeans(flux.all.model1, pairwise ~ Genotype, glhargs=list())
summary(flux.all.Genotype)

flux.all.Genotype_BO <- lsmeans(flux.all.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(flux.all.Genotype_BO)


# Look at the per-animal difference
dotplot(ranef(flux.all.null, condVar = TRUE))$AnimalName
dotplot(ranef(flux.all.model1, condVar = TRUE))$AnimalName
dotplot(ranef(flux.all.model2, condVar = TRUE))$AnimalName
dotplot(ranef(flux.all.model3, condVar = TRUE))$AnimalName
dotplot(ranef(flux.all.model4, condVar = TRUE))$AnimalName

# Look at the model residuals
plot(flux.all.model1)
plot(flux.all.model2)
plot(flux.all.model3)

## Flux and genotype or branch group
flux.art.null = lmer(Flux ~ (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
flux.art.model1 = lmer(Flux~ Genotype + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
flux.art.model2 = lmer(Flux~ BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
flux.art.model3 = lmer(Flux~ Genotype + BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
flux.art.model4 = lmer(Flux~ Genotype * BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
flux.art.anova <- anova(flux.art.null, flux.art.model1,flux.art.model2,flux.art.model3,flux.art.model4)
print(flux.art.anova)
# p values
flux.art.Genotype <- lsmeans(flux.art.model1, pairwise ~ Genotype, glhargs=list())
summary(flux.art.Genotype)

flux.art.Genotype_BO <- lsmeans(flux.art.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(flux.art.Genotype_BO)


# Flux, genotype and depth
flux.nullB = lmer(Flux ~ (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
flux.model1B = lmer(Flux~ Genotype + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
flux.model2B = lmer(Flux~ Depth + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
flux.model3B = lmer(Flux~ Genotype + Depth + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
flux.model4B = lmer(Flux~ Genotype * Depth + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
flux.anovaB <- anova(flux.nullB, flux.model1B,flux.model2B,flux.model3B,flux.model4B)
print(flux.anovaB)


#########
# linearDensity

# only data that has flux measurements

df4A<- summarySE(baseline, measurevar="linearDensity", groupvars=c("Genotype"), na.rm=TRUE)
df4A2<- summarySE(eGFPpos, measurevar="linearDensity", groupvars=c("Genotype"), na.rm=TRUE)
df4B<- summarySE(eGFPpos, measurevar="linearDensity", groupvars=c("Genotype","BranchOrder"), na.rm=TRUE)
df4C<- summarySE(baseline, measurevar="linearDensity", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)
df4D<- summarySE(eGFPpos, measurevar="linearDensity", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)

ggplot(data=df4A, aes(x=Genotype, y=linearDensity, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=linearDensity-se, ymax=linearDensity+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("linearDensity [RBCs/mm]") +
  scale_fill_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


#scatterplot of genotypes- all Data
ggplot() +
  geom_jitter(data=baseline, aes(y=linearDensity, x=Genotype, colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df4A, aes(x=Genotype, y=linearDensity,ymin=linearDensity, ymax=linearDensity, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df4A, aes(x=Genotype, ymin=linearDensity-se, ymax=linearDensity+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("linearDensity, all Data")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

ggplot(data=df4B, aes(x=BranchOrder, y=linearDensity, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=linearDensity-se, ymax=linearDensity+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("linearDensity [RBCs/mm]") +
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme


ggplot(data=df4C, aes(x=BranchGroup, y=linearDensity, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=linearDensity-se, ymax=linearDensity+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("linearDensity [RBCs/mm]") +
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme

#scatterplot of genotypes- all Data
ggplot() +
  geom_jitter(data=baseline, aes(y=linearDensity, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df4C, aes(x=interaction(Genotype,BranchGroup), y=linearDensity,ymin=linearDensity, ymax=linearDensity, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df4C, aes(x=interaction(Genotype,BranchGroup), ymin=linearDensity-se, ymax=linearDensity+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("linearDensity, all Data")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

ggplot(data=df4D, aes(x=BranchGroup, y=linearDensity, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=linearDensity-se, ymax=linearDensity+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("linearDensity [RBCs/mm]") +
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme

#scatterplot of genotypes- all Data
ggplot() +
  geom_jitter(data=eGFPpos, aes(y=linearDensity, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df4D, aes(x=interaction(Genotype,BranchGroup), y=linearDensity,ymin=linearDensity, ymax=linearDensity, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df4D, aes(x=interaction(Genotype,BranchGroup), ymin=linearDensity-se, ymax=linearDensity+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("linearDensity, all Data")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

#######
## boxplots
ggplot(baseline, aes(x = Genotype, y = linearDensity, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("linearDensity [RBCs/mm]") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchGroup), y = linearDensity, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("linearDensity [RBCs/mm]") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme


# scatterplot- linearDensity vs BO
ggplot(baseline, aes(x=BranchOrder, y=linearDensity)) +
  geom_point(aes(colour = Genotype, fill=Genotype),position="dodge", shape = 1, size=2)+
  ggtitle("linearDensity vs order for all vessels") +
  xlab("BranchOrder") + 
  ylab("linearDensity [RBCs/mm]") + 
  scale_colour_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

#depth plots
# scatterplot- linearDensity vs depth
ggplot(baseline, aes(x=linearDensity, y=Depth)) +
  geom_point(aes(colour = Genotype), shape = 1, size=2)+
  ggtitle("linearDensity vs depth for all vessels") +
  xlab("linearDensity [RBCs/mm]") + 
  ylab("Depth [um]") + 
  scale_colour_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme


######
#Stats
## linearDensity and genotype or branch order
linearDensity.all.null = lmer(linearDensity ~ (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
linearDensity.all.model1 = lmer(linearDensity~ Genotype + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
linearDensity.all.model2 = lmer(linearDensity~ BranchGroup + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
linearDensity.all.model3 = lmer(linearDensity~ Genotype + BranchGroup + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
linearDensity.all.model4 = lmer(linearDensity~ Genotype * BranchGroup + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
linearDensity.all.anova <- anova(linearDensity.all.null, linearDensity.all.model1,linearDensity.all.model2,linearDensity.all.model3,linearDensity.all.model4)
print(linearDensity.all.anova)
# p values
linearDensity.all.Genotype <- lsmeans(linearDensity.all.model1, pairwise ~ Genotype, glhargs=list())
summary(linearDensity.all.Genotype)

linearDensity.all.Genotype_BO <- lsmeans(linearDensity.all.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(linearDensity.all.Genotype_BO)


# Look at the per-animal difference
dotplot(ranef(linearDensity.all.null, condVar = TRUE))$AnimalName
dotplot(ranef(linearDensity.all.model1, condVar = TRUE))$AnimalName
dotplot(ranef(linearDensity.all.model2, condVar = TRUE))$AnimalName
dotplot(ranef(linearDensity.all.model3, condVar = TRUE))$AnimalName
dotplot(ranef(linearDensity.all.model4, condVar = TRUE))$AnimalName

# Look at the model residuals
plot(linearDensity.all.model1)
plot(linearDensity.all.model2)
plot(linearDensity.all.model3)

## linearDensity and genotype or branch group
linearDensity.art.null = lmer(linearDensity ~ (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
linearDensity.art.model1 = lmer(linearDensity~ Genotype + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
linearDensity.art.model2 = lmer(linearDensity~ BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos,REML=FALSE)
linearDensity.art.model3 = lmer(linearDensity~ Genotype + BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos ,REML=FALSE)
linearDensity.art.model4 = lmer(linearDensity~ Genotype * BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos ,REML=FALSE)
linearDensity.art.anova <- anova(linearDensity.art.null, linearDensity.art.model1,linearDensity.art.model2,linearDensity.art.model3,linearDensity.art.model4)
print(linearDensity.art.anova)
# p values
linearDensity.art.Genotype <- lsmeans(linearDensity.art.model1, pairwise ~ Genotype, glhargs=list())
summary(linearDensity.art.Genotype)

linearDensity.art.Genotype_BO <- lsmeans(linearDensity.art.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(linearDensity.art.Genotype_BO)

lsmip(linearDensity.art.model4, Genotype~BranchGroup)


# linearDensity, genotype and depth
linearDensity.nullB = lmer(linearDensity ~ (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
linearDensity.model1B = lmer(linearDensity~ Genotype + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
linearDensity.model2B = lmer(linearDensity~ Depth + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
linearDensity.model3B = lmer(linearDensity~ Genotype + Depth + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
linearDensity.model4B = lmer(linearDensity~ Genotype * Depth + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
linearDensity.anovaB <- anova(linearDensity.nullB, linearDensity.model1B,linearDensity.model2B,linearDensity.model3B,linearDensity.model4B)
print(linearDensity.anovaB)

#########
# Hematocrit

# only data that has flux

# NOTE!!! One weird number in retret data!

# outlier test
source("http://goo.gl/UUyEzD")
outlierKD(baseline, Hematocrit)
y
baseline.HC<-baseline[complete.cases(baseline$Hematocrit),]

source("http://goo.gl/UUyEzD")
outlierKD(eGFPpos, Hematocrit)
y
eGFPpos.HC<-eGFPpos[complete.cases(eGFPpos$Hematocrit),]

df5A<- summarySE(baseline.HC, measurevar="Hematocrit", groupvars=c("Genotype"), na.rm=TRUE)
df5A2<- summarySE(eGFPpos.HC, measurevar="Hematocrit", groupvars=c("Genotype"), na.rm=TRUE)
df5B<- summarySE(eGFPpos.HC, measurevar="Hematocrit", groupvars=c("Genotype","BranchOrder"), na.rm=TRUE)
df5C<- summarySE(baseline.HC, measurevar="Hematocrit", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)
df5D<- summarySE(eGFPpos.HC, measurevar="Hematocrit", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)
df5E<- summarySE(baseline.HC, measurevar="Hematocrit", groupvars=c("Genotype","speed"), na.rm=TRUE)

ggplot(data=df5A, aes(x=Genotype, y=Hematocrit, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Hematocrit-se, ymax=Hematocrit+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Hematocrit [RBCs/mm]") +
  scale_fill_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

#scatterplot of genotypes- all Data
ggplot() +
  geom_jitter(data=baseline.HC, aes(y=Hematocrit, x=Genotype, colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df5A, aes(x=Genotype, y=Hematocrit,ymin=Hematocrit, ymax=Hematocrit, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df5A, aes(x=Genotype, ymin=Hematocrit-se, ymax=Hematocrit+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Hematocrit, all Data")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot() +
  geom_bar(data=df5A, aes(x=Genotype, y=Hematocrit, colour=Genotype), fill="white" ,stat="identity", width=0.5, position=position_dodge()) +
  geom_jitter(data=baseline, aes(y=Hematocrit, x=Genotype, colour=Genotype), 
              position=position_jitter(width=0.12), size=4)+
  geom_errorbar(data=df5A, aes(x=Genotype, ymin=Hematocrit-se, ymax=Hematocrit+se), size=2,colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Hematocrit all vessels")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

ggplot() +
  geom_bar(data=df5D, aes(x=interaction(Genotype,BranchGroup), y=Hematocrit, colour=Genotype), fill="white" ,stat="identity", width=0.5, position=position_dodge()) +
  geom_jitter(data=eGFPpos, aes(y=Hematocrit, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4)+
  geom_errorbar(data=df5D, aes(x=interaction(Genotype,BranchGroup), ymin=Hematocrit-se, ymax=Hematocrit+se), size=2,colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("Hematocrit branch groups")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot(data=df5B, aes(x=BranchOrder, y=Hematocrit, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Hematocrit-se, ymax=Hematocrit+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Hematocrit [RBCs/mm]") +
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme


ggplot(data=df5C, aes(x=BranchGroup, y=Hematocrit, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Hematocrit-se, ymax=Hematocrit+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Hematocrit [RBCs/mm]") +
  scale_fill_manual(
    values=c("blue", "red"))+
  max.theme

ggplot() +
  geom_jitter(data=baseline.HC, aes(y=Hematocrit, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df5C, aes(x=interaction(Genotype,BranchGroup), y=Hematocrit,ymin=Hematocrit, ymax=Hematocrit, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df5C, aes(x=interaction(Genotype,BranchGroup), ymin=Hematocrit-se, ymax=Hematocrit+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("hematocrit all groups")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme


ggplot() +
  geom_jitter(data=eGFPpos.HC, aes(y=Hematocrit, x=interaction(Genotype,BranchGroup), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df5D, aes(x=interaction(Genotype,BranchGroup), y=Hematocrit,ymin=Hematocrit, ymax=Hematocrit, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df5D, aes(x=interaction(Genotype,BranchGroup), ymin=Hematocrit-se, ymax=Hematocrit+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("hematocrit, all Data")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

ggplot() +
  geom_jitter(data=baseline.HC, aes(y=Hematocrit, x=interaction(Genotype,speed), colour=Genotype), 
              position=position_jitter(width=0.12), size=4, shape=21)+
  geom_crossbar(data=df5E, aes(x=interaction(Genotype,speed), y=Hematocrit,ymin=Hematocrit, ymax=Hematocrit, colour=Genotype),  width=0.2,position=position_dodge()) +
  geom_errorbar(data=df5E, aes(x=interaction(Genotype,speed), ymin=Hematocrit-se, ymax=Hematocrit+se), colour="black", width=0.1,  position=position_dodge()) +
  ggtitle("hematocrit, all Data")+
  scale_colour_manual(
    values=c("blue", "red"),guide=FALSE)+
  max.theme

#######
## boxplots
ggplot(baseline, aes(x = Genotype, y = Hematocrit, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Hematocrit [RBCs/mm]") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchOrder), y = Hematocrit, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Hematocrit [RBCs/mm]") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchGroup), y = Hematocrit, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Hematocrit [RBCs/mm]") + 
  scale_fill_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme


# scatterplot- Hematocrit vs BO
ggplot(baseline, aes(x=BranchOrder, y=Hematocrit)) +
  geom_point(aes(colour = Genotype, fill=Genotype),position="dodge", shape = 1, size=2)+
  ggtitle("Hematocrit vs order for all vessels") +
  xlab("BranchOrder") + 
  ylab("Hematocrit [RBCs/mm]") + 
  scale_colour_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme

#depth plots
# scatterplot- Hematocrit vs depth
ggplot(baseline, aes(x=Hematocrit, y=Depth)) +
  geom_point(aes(colour = Genotype), shape = 1, size=2)+
  ggtitle("Hematocrit vs depth for all vessels") +
  xlab("Hematocrit [RBCs/mm]") + 
  ylab("Depth [um]") + 
  scale_colour_manual(
    values=c("blue", "red"), 
    guide=FALSE) + 
  max.theme


######
#Stats
## Hematocrit and genotype or branch order
Hematocrit.null = lmer(Hematocrit ~ (1|AnimalName) + (1|Branchname), baseline.HC,REML=FALSE)
Hematocrit.model1 = lmer(Hematocrit~ Genotype + (1|AnimalName) + (1|Branchname), baseline.HC,REML=FALSE)
Hematocrit.model2 = lmer(Hematocrit~ BranchGroup + (1|AnimalName) + (1|Branchname), baseline.HC,REML=FALSE)
Hematocrit.model3 = lmer(Hematocrit~ Genotype + BranchGroup + (1|AnimalName) + (1|Branchname), baseline.HC,REML=FALSE)
Hematocrit.model4 = lmer(Hematocrit~ Genotype * BranchGroup + (1|AnimalName) + (1|Branchname), baseline.HC,REML=FALSE)
Hematocrit.anova <- anova(Hematocrit.null, Hematocrit.model1,Hematocrit.model2,Hematocrit.model3,Hematocrit.model4)
print(Hematocrit.anova)
# p values
Hematocrit.Genotype <- lsmeans(Hematocrit.model1, pairwise ~ Genotype, glhargs=list())
summary(Hematocrit.Genotype)

Hematocrit.Genotype_BO <- lsmeans(Hematocrit.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(Hematocrit.Genotype_BO)


# Look at the per-animal difference
dotplot(ranef(Hematocrit.null, condVar = TRUE))$AnimalName
dotplot(ranef(Hematocrit.model1, condVar = TRUE))$AnimalName
dotplot(ranef(Hematocrit.model2, condVar = TRUE))$AnimalName
dotplot(ranef(Hematocrit.model3, condVar = TRUE))$AnimalName
dotplot(ranef(Hematocrit.model4, condVar = TRUE))$AnimalName

# Look at the model residuals
plot(Hematocrit.model1)
plot(Hematocrit.model2)
plot(Hematocrit.model3)

## Hematocrit and genotype or branch group
Hematocrit.art.null = lmer(Hematocrit ~ (1|AnimalName) + (1|Branchname), eGFPpos.HC,REML=FALSE)
Hematocrit.art.model1 = lmer(Hematocrit~ Genotype + (1|AnimalName) + (1|Branchname), eGFPpos.HC,REML=FALSE)
Hematocrit.art.model2 = lmer(Hematocrit~ BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos.HC,REML=FALSE)
Hematocrit.art.model3 = lmer(Hematocrit~ Genotype + BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos.HC,REML=FALSE)
Hematocrit.art.model4 = lmer(Hematocrit~ Genotype * BranchGroup + (1|AnimalName) + (1|Branchname), eGFPpos.HC,REML=FALSE)
Hematocrit.art.anova <- anova(Hematocrit.art.null, Hematocrit.art.model1,Hematocrit.art.model2,
                              Hematocrit.art.model3,Hematocrit.art.model4)
print(Hematocrit.art.anova)
# p values
Hematocrit.art.Genotype <- lsmeans(Hematocrit.art.model1, pairwise ~ Genotype, glhargs=list())
summary(Hematocrit.art.Genotype)

Hematocrit.art.Genotype_BO <- lsmeans(Hematocrit.art.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(Hematocrit.art.Genotype_BO)

#Hematocrit.Genotype_BO2 <- lsmeans(Hematocrit.model3, pairwise ~ Genotype+BranchGroup, glhargs=list())
#summary(Hematocrit.Genotype_BO2)

# Hematocrit, genotype and depth
Hematocrit.nullB = lmer(Hematocrit ~ (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
Hematocrit.model1B = lmer(Hematocrit~ Genotype + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
Hematocrit.model2B = lmer(Hematocrit~ Depth + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
Hematocrit.model3B = lmer(Hematocrit~ Genotype + Depth + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
Hematocrit.model4B = lmer(Hematocrit~ Genotype * Depth + (1|AnimalName) + (1|Branchname), baseline,REML=FALSE)
Hematocrit.anovaB <- anova(Hematocrit.nullB, Hematocrit.model1B,Hematocrit.model2B,Hematocrit.model3B,Hematocrit.model4B)
print(Hematocrit.anovaB)


######################



#plots (like Chris Schaffer's paper from 2012)

eGFPpos$AdjDiam=-(eGFPpos$Diameter)
eGFPpos$AdjBO=eGFPpos$BranchOrder
eGFPpos$AdjVel=-(eGFPpos$Velocity)
eGFPpos$vessel="artery"
eGFPneg$AdjDiam=eGFPneg$Diameter
eGFPneg$AdjBO=-(eGFPneg$BranchOrder)
eGFPneg$AdjVel=eGFPneg$Velocity
eGFPneg$vessel="vein"

baseline.adjustedDiameter<-rbind(eGFPpos, eGFPneg)
baseline.adjustedDiameter<-baseline.adjustedDiameter[!baseline.adjustedDiameter$Velocity>10,]
baseline.adjustedDiameter<-baseline.adjustedDiameter[complete.cases(baseline.adjustedDiameter$Velocity),]


# plots
library(zoo)

#velocity vs diameter
ggplot(data=baseline.adjustedDiameter, aes(x=AdjDiam, y=Velocity, colour=Genotype)) +
  geom_point(size=3)+
  #geom_line(aes(y=rollmean(Velocity,31, na.pad = TRUE)))+
  facet_grid(. ~vessel,scales = "free" , space="free")+
  ggtitle("diam vs velocity by vessel type") +
  xlab("Diameter") + 
  ylab("Velocity [mm/s]") + 
  scale_colour_manual(values=c("blue", "red")) + 
  max.theme


ggplot(data=baseline.adjustedDiameter[baseline.adjustedDiameter$AdjDiam>-6 & baseline.adjustedDiameter$AdjDiam<6,], aes(x=AdjDiam, y=Velocity, colour=Genotype)) +
  geom_point(size=3)+
  #geom_line(aes(y=rollmean(Velocity,31, na.pad = TRUE)))+
  facet_grid(. ~vessel , scales="free", space="free")+
  ggtitle("diam vs velocity capillaries") +
  xlab("Diameter") + 
  ylab("Velocity [mm/s]") + 
  scale_colour_manual(values=c("blue", "red"), guide=FALSE) + 
  max.theme


# velocity vs branch order

AdjustedBO<-baseline.adjustedDiameter[baseline.adjustedDiameter$BranchGroup=="capillary_PC_A" | baseline.adjustedDiameter$BranchGroup=="capillary_PC_V",]
ggplot(data=AdjustedBO, aes(x=AdjBO, y=Velocity, colour=Genotype)) +
  geom_point(size=3)+
  geom_line(aes(y=rollmean(Velocity,7, na.pad = TRUE)))+
  facet_grid(. ~vessel , scales="free", space="free")+
  ggtitle("branchvs velocity- capillaries") +
  xlab("Branch Order") + 
  ylab("Velocity [mm/s]") + 
  scale_colour_manual(values=c("blue", "red")) + 
  max.theme


# diameter vs pulsatility
ggplot(data=baseline.adjustedDiameter, aes(x=AdjDiam, y=PulsatilityIndex, colour=Genotype)) +
  geom_point(size=3)+
  #geom_line(aes(y=rollmean(Velocity,31, na.pad = TRUE)))+
  facet_grid(. ~vessel,scales = "free" , space="free")+
  ggtitle("diam vs pulsatility by vessel type") +
  xlab("Diameter") + 
  ylab("PI") + 
  scale_colour_manual(values=c("blue", "red")) + 
  max.theme


ggplot(data=baseline.adjustedDiameter, aes(x=AdjVel, y=PulsatilityIndex, colour=Genotype)) +
  geom_point(size=3)+
  #geom_line(aes(y=rollmean(Velocity,31, na.pad = TRUE)))+
  facet_grid(. ~vessel,scales = "free" , space="free")+
  ggtitle("diam vs pulsatility by vessel type") +
  xlab("velocity") + 
  ylab("PI") + 
  scale_colour_manual(values=c("blue", "red")) + 
  max.theme



# linear density vs diameter
ggplot(data=baseline.adjustedDiameter, aes(x=AdjDiam, y=linearDensity, colour=Genotype)) +
  geom_point(size=3)+
  #geom_line(aes(y=rollmean(Velocity,31, na.pad = TRUE)))+
  facet_grid(. ~vessel,scales = "free" , space="free")+
  ggtitle("diam vs LD by vessel type") +
  xlab("diameter") + 
  ylab("LD") + 
  scale_colour_manual(values=c("blue", "red")) + 
  max.theme

# linear density vs branch order
ggplot(data=AdjustedBO, aes(x=AdjBO, y=Hematocrit, colour=Genotype)) +
  geom_point(size=3)+
  #geom_line(aes(y=rollmean(Velocity,31, na.pad = TRUE)))+
  facet_grid(. ~vessel,scales = "free" , space="free")+
  ggtitle("BO vs LD by vessel type") +
  xlab("BO") + 
  ylab("LD") + 
  scale_colour_manual(values=c("blue", "red")) + 
  max.theme

