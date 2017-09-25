
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
baseline <- read.delim("E:/Data/Pericyte_project/Two-photon-data/Ret_ret_Mice/Baseline_BloodFlow/Results/Ret+_diam_velocity_25_09_2017.csv", header=TRUE, sep = "\t")


#########
# unique animal and spot name
baseline$Vesselname <-paste(baseline$AnimalName, baseline$Spot,sep= "_")

baseline$Genotype<- factor(baseline$Genotype,levels = c("Ret+", "RetRet"))
baseline$BranchOrder<- as.factor(baseline$BranchOrder)

###############################
#histograms
ggplot(baseline, aes(x=Diameter, fill=Genotype)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of diameters")

ggplot(baseline, aes(x=Velocity, fill=Genotype))+ geom_histogram(binwidth=0.5, position="dodge") +
  ggtitle("Distribution of velocites")

ggplot(baseline, aes(x=Flux, fill=Genotype)) + geom_histogram(binwidth=10, position="dodge") +
  ggtitle("Distribution of Flux")

ggplot(baseline, aes(x=PulsatilityIndex, fill=Genotype)) + geom_histogram(binwidth=0.1, position="dodge") +
  ggtitle("Distribution of Pulsatility")

ggpairs(baseline, columns = 7:17, aes(colour=Genotype), alpha=0.4)


#########
# diameter
df1A<- summarySE(baseline, measurevar="Diameter", groupvars=c("Genotype"))
df1B<- summarySE(baseline, measurevar="Diameter", groupvars=c("Genotype","BranchOrder"))

ggplot(data=df1A, aes(x=Genotype, y=Diameter, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Diameter-se, ymax=Diameter+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Diameter [um]") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df1B, aes(x=BranchOrder, y=Diameter, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Diameter-se, ymax=Diameter+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Diameter [um]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme

ggplot(data=df1B, aes(x=BranchOrder, y=Diameter, colour=Genotype)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=Diameter-se, ymax=Diameter+se),width=.2) +
  ylab("Diameter [um]") +
  scale_colour_manual(
    values=c("black", "red"))+
  max.theme

## boxplots
ggplot(baseline, aes(x = Genotype, y = Diameter, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Diameter [um]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchOrder), y = Diameter, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Diameter [um]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

#depth plots
# scatterplot- diameter vs BO
ggplot(baseline, aes(x=BranchOrder, y=Diameter)) +
  geom_point(aes(colour = Genotype, fill=Genotype),position="dodge", shape = 1, size=2)+
  ggtitle("diameter vs order for all vessels") +
  xlab("BranchOrder") + 
  ylab("Diameter [um]") + 
  scale_colour_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

# scatterplot- diameter vs depth
ggplot(baseline, aes(x=Diameter, y=Depth)) +
  geom_point(aes(colour = Genotype), shape = 1, size=2)+
  ggtitle("diameter vs depth for all vessels") +
  xlab("Diameter [um]") + 
  ylab("Depth [um]") + 
  scale_colour_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


######
#Stats
## Diameter and genotype or branch order
diam.null = lmer(Diameter ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
diam.model1 = lmer(Diameter~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
diam.model2 = lmer(Diameter~ BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
diam.model3 = lmer(Diameter~ Genotype + BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
diam.model4 = lmer(Diameter~ Genotype * BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
diam.anova <- anova(diam.null, diam.model1,diam.model2,diam.model3,diam.model4)
print(diam.anova)
# p values
diam.Genotype <- lsmeans(diam.model1, pairwise ~ Genotype, glhargs=list())
summary(diam.Genotype)

diam.Genotype_BO <- lsmeans(diam.model4, pairwise ~ Genotype*BranchOrder, glhargs=list())
summary(diam.Genotype_BO)


# Look at the per-animal difference
dotplot(ranef(diam.null, condVar = TRUE))$AnimalName
dotplot(ranef(diam.model1, condVar = TRUE))$AnimalName
dotplot(ranef(diam.model2, condVar = TRUE))$AnimalName
dotplot(ranef(diam.model3, condVar = TRUE))$AnimalName
dotplot(ranef(diam.model4, condVar = TRUE))$AnimalName

# Look at the model residuals
plot(diam.model1)
plot(diam.model2)
plot(diam.model3)

# diameter, genotype and depth
diam.nullB = lmer(Diameter ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
diam.model1B = lmer(Diameter~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
diam.model2B = lmer(Diameter~ Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
diam.model3B = lmer(Diameter~ Genotype + Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
diam.model4B = lmer(Diameter~ Genotype * Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
diam.anovaB <- anova(diam.nullB, diam.model1B,diam.model2B,diam.model3B,diam.model4B)
print(diam.anovaB)


# NOTE: an effect of depth and genotype on diameter!





#########
# velocity
df2A<- summarySE(baseline, measurevar="Velocity", groupvars=c("Genotype"), na.rm=TRUE)
df2B<- summarySE(baseline, measurevar="Velocity", groupvars=c("Genotype","BranchOrder"), na.rm=TRUE)

ggplot(data=df2A, aes(x=Genotype, y=Velocity, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Velocity-se, ymax=Velocity+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Velocity [mm/s]") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df2B, aes(x=BranchOrder, y=Velocity, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Velocity-se, ymax=Velocity+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Velocity [mm/s]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme

ggplot(data=df2B, aes(x=BranchOrder, y=Velocity, colour=Genotype)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=Velocity-se, ymax=Velocity+se),width=.2) +
  ylab("Velocity [mm/s]") +
  scale_colour_manual(
    values=c("black", "red"))+
  max.theme

## boxplots
ggplot(baseline, aes(x = Genotype, y = Velocity, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Velocity [mm/s]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchOrder), y = Velocity, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Velocity [mm/s]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

#depth plots
# scatterplot- Velocity vs BO
ggplot(baseline, aes(x=BranchOrder, y=Velocity)) +
  geom_point(aes(colour = Genotype, fill=Genotype),position="dodge", shape = 1, size=2)+
  ggtitle("Velocity vs order for all vessels") +
  xlab("BranchOrder") + 
  ylab("Velocity [mm/s]") + 
  scale_colour_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

# scatterplot- Velocity vs depth
ggplot(baseline, aes(x=Velocity, y=Depth)) +
  geom_point(aes(colour = Genotype), shape = 1, size=2)+
  ggtitle("Velocity vs depth for all vessels") +
  xlab("Velocity [mm/s]") + 
  ylab("Depth [mm/s]") + 
  scale_colour_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


######
#Stats
## Velocity and genotype or branch order
vel.null = lmer(Velocity ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.model1 = lmer(Velocity~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.model2 = lmer(Velocity~ BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.model3 = lmer(Velocity~ Genotype + BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.model4 = lmer(Velocity~ Genotype * BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.anova <- anova(vel.null, vel.model1,vel.model2,vel.model3,vel.model4)
print(vel.anova)
# p values
vel.Genotype <- lsmeans(vel.model1, pairwise ~ Genotype, glhargs=list())
summary(vel.Genotype)

vel.Genotype_BO <- lsmeans(vel.model4, pairwise ~ Genotype*BranchOrder, glhargs=list())
summary(vel.Genotype_BO)


# Look at the per-animal difference
dotplot(ranef(vel.null, condVar = TRUE))$AnimalName
dotplot(ranef(vel.model1, condVar = TRUE))$AnimalName
dotplot(ranef(vel.model2, condVar = TRUE))$AnimalName
dotplot(ranef(vel.model3, condVar = TRUE))$AnimalName
dotplot(ranef(vel.model4, condVar = TRUE))$AnimalName

# Look at the model residuals
plot(vel.model1)
plot(vel.model2)
plot(vel.model3)

# Velocity, genotype and depth
vel.nullB = lmer(Velocity ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.model1B = lmer(Velocity~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.model2B = lmer(Velocity~ Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.model3B = lmer(Velocity~ Genotype + Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.model4B = lmer(Velocity~ Genotype * Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.anovaB <- anova(vel.nullB, vel.model1B,vel.model2B,vel.model3B,vel.model4B)
print(vel.anovaB)


# NOTE: an effect of depth and genotype on Velocity!




# Group data with NO FLOW
