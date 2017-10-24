
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
baseline <- read.delim("E:/Data/Pericyte_project/Two-photon-data/Ret_ret_Mice/Baseline_BloodFlow/Results/Diam_velocity_24_10_2017.csv", header=TRUE, sep = "\t")


#########
# unique animal and spot name
baseline$Vesselname <-paste(baseline$AnimalName, baseline$Spot,sep= "_")
baseline$Branchname <-paste(baseline$AnimalName, baseline$Branch,sep= "_")

baseline$Genotype<- factor(baseline$Genotype,levels = c("Ret+", "RetRet"))
baseline$BranchOrder<- as.factor(baseline$BranchOrder)


# group branch order to simplify analysis
# groups 1-3, 4-6, 7-9
baseline$BranchGroup<-"0"
baseline$BranchGroup[baseline$BranchOrder==1]<-"1-3"
baseline$BranchGroup[baseline$BranchOrder==2]<-"1-3"
baseline$BranchGroup[baseline$BranchOrder==3]<-"1-3"
baseline$BranchGroup[baseline$BranchOrder==4]<-"4-6"
baseline$BranchGroup[baseline$BranchOrder==5]<-"4-6"
baseline$BranchGroup[baseline$BranchOrder==6]<-"4-6"
baseline$BranchGroup[baseline$BranchOrder==7]<-"7-9"
baseline$BranchGroup[baseline$BranchOrder==8]<-"7-9"
baseline$BranchGroup[baseline$BranchOrder==9]<-"7-9"

baseline$BranchGroup<- as.factor(baseline$BranchGroup)


# only consider data from vessels were we have all the info (flux, etc.)
baseline<-baseline[complete.cases(baseline$Flux),]


# only consider branch order greater than 4

baseline<- baseline[!baseline$BranchGroup=="1-3",]

###############################
#histograms
ggplot(baseline, aes(x=Diameter, fill=Genotype)) + geom_histogram(binwidth=1, position="dodge") +
  ggtitle("Distribution of diameters")

ggplot(baseline, aes(x=Velocity, fill=Genotype))+ geom_histogram(binwidth=0.2, position="dodge") +
  ggtitle("Distribution of velocites")

ggplot(baseline, aes(x=Velocity, fill=Genotype)) +
  geom_histogram(binwidth=0.2, position="dodge", aes(y=..count../sum(..count..)))

ggplot(baseline, aes(x=Velocity, y=..density..,fill=Genotype)) +
  geom_histogram(binwidth=0.2, position="dodge")

ggplot(baseline, aes(x=Flux, fill=Genotype)) + geom_histogram(binwidth=10, position="dodge") +
  ggtitle("Distribution of Flux")

ggplot(baseline, aes(x=Flux, y=..density..,fill=Genotype)) +
  geom_histogram(bidwidth=10, position="dodge")

ggplot(baseline, aes(x=PulsatilityIndex, fill=Genotype)) + geom_histogram(binwidth=0.1, position="dodge") +
  ggtitle("Distribution of Pulsatility")

#ggpairs(baseline, columns = 7:17, aes(colour=Genotype), alpha=0.4)


#########
# diameter
df1A<- summarySE(baseline, measurevar="Diameter", groupvars=c("Genotype"))
df1B<- summarySE(baseline, measurevar="Diameter", groupvars=c("Genotype","BranchOrder"))
df1C<- summarySE(baseline, measurevar="Diameter", groupvars=c("Genotype","BranchGroup"))

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

ggplot(data=df1C, aes(x=BranchGroup, y=Diameter, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Diameter-se, ymax=Diameter+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Diameter [um]") +
  scale_fill_manual(
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

ggplot(baseline, aes(x = interaction(Genotype,BranchGroup), y = Diameter, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Diameter [um]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

# categorical scatterplots
ggplot(baseline, aes(x=BranchOrder, y=Diameter, colour=Genotype))+
  geom_jitter() +
  geom_crossbar(data=df1B,aes(x=BranchOrder,ymin=Diameter, ymax=Diameter,y=Diameter,group=BranchOrder), width = 0.5) +
  scale_colour_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x=BranchGroup, y=Diameter, colour=Genotype))+
  geom_jitter() +
  geom_crossbar(data=df1C,aes(x=BranchGroup,ymin=Diameter, ymax=Diameter,y=Diameter,group=BranchGroup), width = 0.5) +
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


##############
# depth 

#depth plots

ggplot(baseline, aes(x = interaction(Genotype,BranchOrder), y = Depth, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("depth [um]") + 
  scale_fill_manual(
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

# depth per branch order
# branch order vs depth
ggplot(baseline, aes(x=BranchOrder, y=Depth)) +
  geom_jitter(aes(colour = Genotype), shape = 1, size=2)+
  ggtitle("branch order vs depth for all vessels") +
  xlab("Branch Order") + 
  ylab("Depth [um]") + 
  scale_colour_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


ggplot(baseline, aes(x = interaction(Genotype,BranchGroup), y = Depth, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("depth [um]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


# depth and branch order
depth.null = lmer(Depth ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
depth.model1 = lmer(Depth~ BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
depth.model2 = lmer(Depth~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
depth.model3 = lmer(Depth~ Genotype + BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
depth.model4 = lmer(Depth~ Genotype * BranchOrder+ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
depth.anova <- anova(depth.null, depth.model1, depth.model2, depth.model3, depth.model4)
print(depth.anova)

depth.Genotype_BO <- lsmeans(depth.model4, pairwise ~ Genotype*BranchOrder, glhargs=list())
summary(depth.Genotype_BO)

# depth and branch group
depth.null = lmer(Depth ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
depth.model1 = lmer(Depth~ BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
depth.model2 = lmer(Depth~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
depth.model3 = lmer(Depth~ Genotype + BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
depth.model4 = lmer(Depth~ Genotype * BranchGroup+ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
depth.anova <- anova(depth.null, depth.model1, depth.model2, depth.model3, depth.model4)
print(depth.anova)

depth.Genotype_BG <- lsmeans(depth.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(depth.Genotype_BG)


# no significant different between genotypes and branch order at certain depths,
# ie. the ret ret do not seem to branch more at earlier depths

#########
# velocity

# only data that has 

df2A<- summarySE(baseline, measurevar="Velocity", groupvars=c("Genotype"), na.rm=TRUE)
df2B<- summarySE(baseline, measurevar="Velocity", groupvars=c("Genotype","BranchOrder"), na.rm=TRUE)
df2C<- summarySE(baseline, measurevar="Velocity", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)

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

ggplot(data=df2C, aes(x=BranchGroup, y=Velocity, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Velocity-se, ymax=Velocity+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Velocity [mm/s]") +
  scale_fill_manual(
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

ggplot(baseline, aes(x = interaction(Genotype,BranchGroup), y = Velocity, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Velocity [mm/s]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


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

#depth plots
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

## Velocity and genotype or branch group
vel.null = lmer(Velocity ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.model1 = lmer(Velocity~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.model2 = lmer(Velocity~ BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.model3 = lmer(Velocity~ Genotype + BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.model4 = lmer(Velocity~ Genotype * BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.anova <- anova(vel.null, vel.model1,vel.model2,vel.model3,vel.model4)
print(vel.anova)
# p values
vel.Genotype <- lsmeans(vel.model1, pairwise ~ Genotype, glhargs=list())
summary(vel.Genotype)

vel.Genotype_BO <- lsmeans(vel.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(vel.Genotype_BO)


# Velocity, genotype and depth
vel.nullB = lmer(Velocity ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.model1B = lmer(Velocity~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.model2B = lmer(Velocity~ Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.model3B = lmer(Velocity~ Genotype + Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.model4B = lmer(Velocity~ Genotype * Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
vel.anovaB <- anova(vel.nullB, vel.model1B,vel.model2B,vel.model3B,vel.model4B)
print(vel.anovaB)


#########
# Flux

# only data that has 

df3A<- summarySE(baseline, measurevar="Flux", groupvars=c("Genotype"), na.rm=TRUE)
df3B<- summarySE(baseline, measurevar="Flux", groupvars=c("Genotype","BranchOrder"), na.rm=TRUE)
df3C<- summarySE(baseline, measurevar="Flux", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)

ggplot(data=df3A, aes(x=Genotype, y=Flux, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Flux-se, ymax=Flux+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Flux [RBCs/s]") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df3B, aes(x=BranchOrder, y=Flux, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Flux-se, ymax=Flux+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Flux [RBCs/s]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme

ggplot(data=df3B, aes(x=BranchOrder, y=Flux, colour=Genotype)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=Flux-se, ymax=Flux+se),width=.2) +
  ylab("Flux [RBCs/s]") +
  scale_colour_manual(
    values=c("black", "red"))+
  max.theme

ggplot(data=df3C, aes(x=BranchGroup, y=Flux, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Flux-se, ymax=Flux+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Flux [RBCs/s]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme

## boxplots
ggplot(baseline, aes(x = Genotype, y = Flux, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Flux [RBCs/s]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchOrder), y = Flux, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Flux [RBCs/s]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchGroup), y = Flux, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Flux [RBCs/s]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


# scatterplot- Flux vs BO
ggplot(baseline, aes(x=BranchOrder, y=Flux)) +
  geom_point(aes(colour = Genotype, fill=Genotype),position="dodge", shape = 1, size=2)+
  ggtitle("Flux vs order for all vessels") +
  xlab("BranchOrder") + 
  ylab("Flux [RBCs/s]") + 
  scale_colour_manual(
    values=c("black", "red"), 
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
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


######
#Stats
## Flux and genotype or branch order
flux.null = lmer(Flux ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
flux.model1 = lmer(Flux~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
flux.model2 = lmer(Flux~ BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
flux.model3 = lmer(Flux~ Genotype + BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
flux.model4 = lmer(Flux~ Genotype * BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
flux.anova <- anova(flux.null, flux.model1,flux.model2,flux.model3,flux.model4)
print(flux.anova)
# p values
flux.Genotype <- lsmeans(flux.model1, pairwise ~ Genotype, glhargs=list())
summary(flux.Genotype)

flux.Genotype_BO <- lsmeans(flux.model4, pairwise ~ Genotype*BranchOrder, glhargs=list())
summary(flux.Genotype_BO)


# Look at the per-animal difference
dotplot(ranef(flux.null, condVar = TRUE))$AnimalName
dotplot(ranef(flux.model1, condVar = TRUE))$AnimalName
dotplot(ranef(flux.model2, condVar = TRUE))$AnimalName
dotplot(ranef(flux.model3, condVar = TRUE))$AnimalName
dotplot(ranef(flux.model4, condVar = TRUE))$AnimalName

# Look at the model residuals
plot(flux.model1)
plot(flux.model2)
plot(flux.model3)

## Flux and genotype or branch group
flux.null = lmer(Flux ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
flux.model1 = lmer(Flux~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
flux.model2 = lmer(Flux~ BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
flux.model3 = lmer(Flux~ Genotype + BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
flux.model4 = lmer(Flux~ Genotype * BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
flux.anova <- anova(flux.null, flux.model1,flux.model2,flux.model3,flux.model4)
print(flux.anova)
# p values
flux.Genotype <- lsmeans(flux.model1, pairwise ~ Genotype, glhargs=list())
summary(flux.Genotype)

flux.Genotype_BO <- lsmeans(flux.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(flux.Genotype_BO)


# Flux, genotype and depth
flux.nullB = lmer(Flux ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
flux.model1B = lmer(Flux~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
flux.model2B = lmer(Flux~ Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
flux.model3B = lmer(Flux~ Genotype + Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
flux.model4B = lmer(Flux~ Genotype * Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
flux.anovaB <- anova(flux.nullB, flux.model1B,flux.model2B,flux.model3B,flux.model4B)
print(flux.anovaB)


#########
# linearDensity

# only data that has 

df4A<- summarySE(baseline, measurevar="linearDensity", groupvars=c("Genotype"), na.rm=TRUE)
df4B<- summarySE(baseline, measurevar="linearDensity", groupvars=c("Genotype","BranchOrder"), na.rm=TRUE)
df4C<- summarySE(baseline, measurevar="linearDensity", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)

ggplot(data=df4A, aes(x=Genotype, y=linearDensity, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=linearDensity-se, ymax=linearDensity+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("linearDensity [RBCs/mm]") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df4B, aes(x=BranchOrder, y=linearDensity, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=linearDensity-se, ymax=linearDensity+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("linearDensity [RBCs/mm]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme

ggplot(data=df4B, aes(x=BranchOrder, y=linearDensity, colour=Genotype)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=linearDensity-se, ymax=linearDensity+se),width=.2) +
  ylab("linearDensity [RBCs/mm]") +
  scale_colour_manual(
    values=c("black", "red"))+
  max.theme

ggplot(data=df4C, aes(x=BranchGroup, y=linearDensity, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=linearDensity-se, ymax=linearDensity+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("linearDensity [RBCs/mm]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme

## boxplots
ggplot(baseline, aes(x = Genotype, y = linearDensity, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("linearDensity [RBCs/mm]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchOrder), y = linearDensity, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("linearDensity [RBCs/mm]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchGroup), y = linearDensity, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("linearDensity [RBCs/mm]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


# scatterplot- linearDensity vs BO
ggplot(baseline, aes(x=BranchOrder, y=linearDensity)) +
  geom_point(aes(colour = Genotype, fill=Genotype),position="dodge", shape = 1, size=2)+
  ggtitle("linearDensity vs order for all vessels") +
  xlab("BranchOrder") + 
  ylab("linearDensity [RBCs/mm]") + 
  scale_colour_manual(
    values=c("black", "red"), 
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
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


######
#Stats
## linearDensity and genotype or branch order
linearDensity.null = lmer(linearDensity ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
linearDensity.model1 = lmer(linearDensity~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
linearDensity.model2 = lmer(linearDensity~ BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
linearDensity.model3 = lmer(linearDensity~ Genotype + BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
linearDensity.model4 = lmer(linearDensity~ Genotype * BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
linearDensity.anova <- anova(linearDensity.null, linearDensity.model1,linearDensity.model2,linearDensity.model3,linearDensity.model4)
print(linearDensity.anova)
# p values
linearDensity.Genotype <- lsmeans(linearDensity.model1, pairwise ~ Genotype, glhargs=list())
summary(linearDensity.Genotype)

linearDensity.Genotype_BO <- lsmeans(linearDensity.model4, pairwise ~ Genotype*BranchOrder, glhargs=list())
summary(linearDensity.Genotype_BO)


# Look at the per-animal difference
dotplot(ranef(linearDensity.null, condVar = TRUE))$AnimalName
dotplot(ranef(linearDensity.model1, condVar = TRUE))$AnimalName
dotplot(ranef(linearDensity.model2, condVar = TRUE))$AnimalName
dotplot(ranef(linearDensity.model3, condVar = TRUE))$AnimalName
dotplot(ranef(linearDensity.model4, condVar = TRUE))$AnimalName

# Look at the model residuals
plot(linearDensity.model1)
plot(linearDensity.model2)
plot(linearDensity.model3)

## linearDensity and genotype or branch group
linearDensity.null = lmer(linearDensity ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
linearDensity.model1 = lmer(linearDensity~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
linearDensity.model2 = lmer(linearDensity~ BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
linearDensity.model3 = lmer(linearDensity~ Genotype + BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
linearDensity.model4 = lmer(linearDensity~ Genotype * BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
linearDensity.anova <- anova(linearDensity.null, linearDensity.model1,linearDensity.model2,linearDensity.model3,linearDensity.model4)
print(linearDensity.anova)
# p values
linearDensity.Genotype <- lsmeans(linearDensity.model1, pairwise ~ Genotype, glhargs=list())
summary(linearDensity.Genotype)

linearDensity.Genotype_BO <- lsmeans(linearDensity.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(linearDensity.Genotype_BO)

#linearDensity.Genotype_BO2 <- lsmeans(linearDensity.model3, pairwise ~ Genotype+BranchGroup, glhargs=list())
#summary(linearDensity.Genotype_BO2)

# linearDensity, genotype and depth
linearDensity.nullB = lmer(linearDensity ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
linearDensity.model1B = lmer(linearDensity~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
linearDensity.model2B = lmer(linearDensity~ Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
linearDensity.model3B = lmer(linearDensity~ Genotype + Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
linearDensity.model4B = lmer(linearDensity~ Genotype * Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
linearDensity.anovaB <- anova(linearDensity.nullB, linearDensity.model1B,linearDensity.model2B,linearDensity.model3B,linearDensity.model4B)
print(linearDensity.anovaB)

#########
# Hematocrit

# only data that has flux

df5A<- summarySE(baseline, measurevar="Hematocrit", groupvars=c("Genotype"), na.rm=TRUE)
df5B<- summarySE(baseline, measurevar="Hematocrit", groupvars=c("Genotype","BranchOrder"), na.rm=TRUE)
df5C<- summarySE(baseline, measurevar="Hematocrit", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)

ggplot(data=df5A, aes(x=Genotype, y=Hematocrit, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Hematocrit-se, ymax=Hematocrit+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Hematocrit [RBCs/mm]") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df5B, aes(x=BranchOrder, y=Hematocrit, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Hematocrit-se, ymax=Hematocrit+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Hematocrit [RBCs/mm]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme

ggplot(data=df5B, aes(x=BranchOrder, y=Hematocrit, colour=Genotype)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=Hematocrit-se, ymax=Hematocrit+se),width=.2) +
  ylab("Hematocrit [RBCs/mm]") +
  scale_colour_manual(
    values=c("black", "red"))+
  max.theme

ggplot(data=df5C, aes(x=BranchGroup, y=Hematocrit, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Hematocrit-se, ymax=Hematocrit+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Hematocrit [RBCs/mm]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme

## boxplots
ggplot(baseline, aes(x = Genotype, y = Hematocrit, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Hematocrit [RBCs/mm]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchOrder), y = Hematocrit, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Hematocrit [RBCs/mm]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchGroup), y = Hematocrit, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Hematocrit [RBCs/mm]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


# scatterplot- Hematocrit vs BO
ggplot(baseline, aes(x=BranchOrder, y=Hematocrit)) +
  geom_point(aes(colour = Genotype, fill=Genotype),position="dodge", shape = 1, size=2)+
  ggtitle("Hematocrit vs order for all vessels") +
  xlab("BranchOrder") + 
  ylab("Hematocrit [RBCs/mm]") + 
  scale_colour_manual(
    values=c("black", "red"), 
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
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


######
#Stats
## Hematocrit and genotype or branch order
Hematocrit.null = lmer(Hematocrit ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
Hematocrit.model1 = lmer(Hematocrit~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
Hematocrit.model2 = lmer(Hematocrit~ BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
Hematocrit.model3 = lmer(Hematocrit~ Genotype + BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
Hematocrit.model4 = lmer(Hematocrit~ Genotype * BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
Hematocrit.anova <- anova(Hematocrit.null, Hematocrit.model1,Hematocrit.model2,Hematocrit.model3,Hematocrit.model4)
print(Hematocrit.anova)
# p values
Hematocrit.Genotype <- lsmeans(Hematocrit.model1, pairwise ~ Genotype, glhargs=list())
summary(Hematocrit.Genotype)

Hematocrit.Genotype_BO <- lsmeans(Hematocrit.model4, pairwise ~ Genotype*BranchOrder, glhargs=list())
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
Hematocrit.null = lmer(Hematocrit ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
Hematocrit.model1 = lmer(Hematocrit~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
Hematocrit.model2 = lmer(Hematocrit~ BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
Hematocrit.model3 = lmer(Hematocrit~ Genotype + BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
Hematocrit.model4 = lmer(Hematocrit~ Genotype * BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
Hematocrit.anova <- anova(Hematocrit.null, Hematocrit.model1,Hematocrit.model2,Hematocrit.model3,Hematocrit.model4)
print(Hematocrit.anova)
# p values
Hematocrit.Genotype <- lsmeans(Hematocrit.model1, pairwise ~ Genotype, glhargs=list())
summary(Hematocrit.Genotype)

Hematocrit.Genotype_BO <- lsmeans(Hematocrit.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(Hematocrit.Genotype_BO)

#Hematocrit.Genotype_BO2 <- lsmeans(Hematocrit.model3, pairwise ~ Genotype+BranchGroup, glhargs=list())
#summary(Hematocrit.Genotype_BO2)

# Hematocrit, genotype and depth
Hematocrit.nullB = lmer(Hematocrit ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
Hematocrit.model1B = lmer(Hematocrit~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
Hematocrit.model2B = lmer(Hematocrit~ Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
Hematocrit.model3B = lmer(Hematocrit~ Genotype + Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
Hematocrit.model4B = lmer(Hematocrit~ Genotype * Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
Hematocrit.anovaB <- anova(Hematocrit.nullB, Hematocrit.model1B,Hematocrit.model2B,Hematocrit.model3B,Hematocrit.model4B)
print(Hematocrit.anovaB)


######################
# Pulsatility Index

df6A<- summarySE(baseline, measurevar="PulsatilityIndex", groupvars=c("Genotype"), na.rm=TRUE)
df6B<- summarySE(baseline, measurevar="PulsatilityIndex", groupvars=c("Genotype","BranchOrder"), na.rm=TRUE)
df6C<- summarySE(baseline, measurevar="PulsatilityIndex", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)

ggplot(data=df6A, aes(x=Genotype, y=PulsatilityIndex, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=PulsatilityIndex-se, ymax=PulsatilityIndex+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("PulsatilityIndex") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df6B, aes(x=BranchOrder, y=PulsatilityIndex, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=PulsatilityIndex-se, ymax=PulsatilityIndex+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("PulsatilityIndex") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme

ggplot(data=df6B, aes(x=BranchOrder, y=PulsatilityIndex, colour=Genotype)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=PulsatilityIndex-se, ymax=PulsatilityIndex+se),width=.2) +
  ylab("PulsatilityIndex") +
  scale_colour_manual(
    values=c("black", "red"))+
  max.theme

ggplot(data=df6C, aes(x=BranchGroup, y=PulsatilityIndex, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=PulsatilityIndex-se, ymax=PulsatilityIndex+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("PulsatilityIndex") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme

## boxplots
ggplot(baseline, aes(x = Genotype, y = PulsatilityIndex, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("PulsatilityIndex") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchOrder), y = PulsatilityIndex, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("PulsatilityIndex") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(baseline, aes(x = interaction(Genotype,BranchGroup), y = PulsatilityIndex, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("PulsatilityIndex") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


# scatterplot- PulsatilityIndex vs BO
ggplot(baseline, aes(x=BranchOrder, y=PulsatilityIndex)) +
  geom_point(aes(colour = Genotype, fill=Genotype),position="dodge", shape = 1, size=2)+
  ggtitle("PulsatilityIndex vs order for all vessels") +
  xlab("BranchOrder") + 
  ylab("PulsatilityIndex") + 
  scale_colour_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

#depth plots
# scatterplot- PulsatilityIndex vs depth
ggplot(baseline, aes(x=PulsatilityIndex, y=Depth)) +
  geom_point(aes(colour = Genotype), shape = 1, size=2)+
  ggtitle("PulsatilityIndex vs depth for all vessels") +
  xlab("PulsatilityIndex [RBCs/mm]") + 
  ylab("Depth [um]") + 
  scale_colour_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


######
#Stats
## PulsatilityIndex and genotype or branch order
PulsatilityIndex.null = lmer(PulsatilityIndex ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
PulsatilityIndex.model1 = lmer(PulsatilityIndex~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
PulsatilityIndex.model2 = lmer(PulsatilityIndex~ BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
PulsatilityIndex.model3 = lmer(PulsatilityIndex~ Genotype + BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
PulsatilityIndex.model4 = lmer(PulsatilityIndex~ Genotype * BranchOrder + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
PulsatilityIndex.anova <- anova(PulsatilityIndex.null, PulsatilityIndex.model1,PulsatilityIndex.model2,PulsatilityIndex.model3,PulsatilityIndex.model4)
print(PulsatilityIndex.anova)
# p values
PulsatilityIndex.Genotype <- lsmeans(PulsatilityIndex.model1, pairwise ~ Genotype, glhargs=list())
summary(PulsatilityIndex.Genotype)

PulsatilityIndex.Genotype_BO <- lsmeans(PulsatilityIndex.model4, pairwise ~ Genotype*BranchOrder, glhargs=list())
summary(PulsatilityIndex.Genotype_BO)


# Look at the per-animal difference
dotplot(ranef(PulsatilityIndex.null, condVar = TRUE))$AnimalName
dotplot(ranef(PulsatilityIndex.model1, condVar = TRUE))$AnimalName
dotplot(ranef(PulsatilityIndex.model2, condVar = TRUE))$AnimalName
dotplot(ranef(PulsatilityIndex.model3, condVar = TRUE))$AnimalName
dotplot(ranef(PulsatilityIndex.model4, condVar = TRUE))$AnimalName

# Look at the model residuals
plot(PulsatilityIndex.model1)
plot(PulsatilityIndex.model2)
plot(PulsatilityIndex.model3)
plot(PulsatilityIndex.model4)

## PulsatilityIndex and genotype or branch group
PulsatilityIndex.null = lmer(PulsatilityIndex ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
PulsatilityIndex.model1 = lmer(PulsatilityIndex~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
PulsatilityIndex.model2 = lmer(PulsatilityIndex~ BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
PulsatilityIndex.model3 = lmer(PulsatilityIndex~ Genotype + BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
PulsatilityIndex.model4 = lmer(PulsatilityIndex~ Genotype * BranchGroup + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
PulsatilityIndex.anova <- anova(PulsatilityIndex.null, PulsatilityIndex.model1,PulsatilityIndex.model2,PulsatilityIndex.model3,PulsatilityIndex.model4)
print(PulsatilityIndex.anova)
# p values
PulsatilityIndex.Genotype <- lsmeans(PulsatilityIndex.model1, pairwise ~ Genotype, glhargs=list())
summary(PulsatilityIndex.Genotype)

PulsatilityIndex.Genotype_BO <- lsmeans(PulsatilityIndex.model4, pairwise ~ Genotype*BranchGroup, glhargs=list())
summary(PulsatilityIndex.Genotype_BO)

#PulsatilityIndex.Genotype_BO2 <- lsmeans(PulsatilityIndex.model3, pairwise ~ Genotype+BranchGroup, glhargs=list())
#summary(PulsatilityIndex.Genotype_BO2)

# PulsatilityIndex, genotype and depth
PulsatilityIndex.nullB = lmer(PulsatilityIndex ~ (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
PulsatilityIndex.model1B = lmer(PulsatilityIndex~ Genotype + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
PulsatilityIndex.model2B = lmer(PulsatilityIndex~ Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
PulsatilityIndex.model3B = lmer(PulsatilityIndex~ Genotype + Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
PulsatilityIndex.model4B = lmer(PulsatilityIndex~ Genotype * Depth + (1|AnimalName) + (1|Spot), baseline,REML=FALSE)
PulsatilityIndex.anovaB <- anova(PulsatilityIndex.nullB, PulsatilityIndex.model1B,PulsatilityIndex.model2B,PulsatilityIndex.model3B,PulsatilityIndex.model4B)
print(PulsatilityIndex.anovaB)


######
df7A<- summarySE(baseline, measurevar="SD_LD", groupvars=c("Genotype"), na.rm=TRUE)
df7B<- summarySE(baseline, measurevar="linearDensity", groupvars=c("Genotype","BranchOrder"), na.rm=TRUE)
df7C<- summarySE(baseline, measurevar="linearDensity", groupvars=c("Genotype","BranchGroup"), na.rm=TRUE)

ggplot(data=df7A, aes(x=Genotype, y=SD_LD, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=SD_LD-se, ymax=SD_LD+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("SD_LD") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df4B, aes(x=BranchOrder, y=linearDensity, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=linearDensity-se, ymax=linearDensity+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("linearDensity [RBCs/mm]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme

ggplot(data=df4B, aes(x=BranchOrder, y=linearDensity, colour=Genotype)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin=linearDensity-se, ymax=linearDensity+se),width=.2) +
  ylab("linearDensity [RBCs/mm]") +
  scale_colour_manual(
    values=c("black", "red"))+
  max.theme

ggplot(data=df4C, aes(x=BranchGroup, y=linearDensity, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=linearDensity-se, ymax=linearDensity+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("linearDensity [RBCs/mm]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme

##########################
# vessel data with NO FLOW
baseline$flow<-"yes"
baseline$flow[baseline$NoFlow==1]<-"no"

df.diam.noflow1<- summarySE(baseline, measurevar="Diameter", groupvars=c("Genotype", "flow"), na.rm=TRUE)
df.diam.noflow2<- summarySE(baseline, measurevar="Diameter", groupvars=c("Genotype","BranchOrder","flow"), na.rm=TRUE)
df.diam.noflow3<- summarySE(baseline, measurevar="Diameter", groupvars=c("Genotype","BranchGroup","flow"), na.rm=TRUE)

ggplot(data=df.diam.noflow1, aes(x=Genotype, y=Diameter, fill=flow)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Diameter-se, ymax=Diameter+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Diameter") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df.diam.noflow3, aes(x=interaction(Genotype,BranchGroup), y=Diameter, fill=flow)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Diameter-se, ymax=Diameter+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Diameter") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme


# total non-flowing vessels
noFlow<- ddply(baseline, c("Genotype","flow"), summarise, nVessels=length(Diameter))

# percentages
noFlow.retplus<-subset(noFlow, Genotype=="Ret+")
noFlow.retret<-subset(noFlow, Genotype=="RetRet")

total_Ret_plus=sum(noFlow$nVessels[noFlow$Genotype=="Ret+"])
total_Ret_Ret=sum(noFlow$nVessels[noFlow$Genotype=="RetRet"])

percent_Ret_plus=(noFlow.retplus$nVessels[noFlow.retplus$flow=="no"]/total_Ret_plus)*100
percent_Ret_Ret=(noFlow.retret$nVessels[noFlow.retret$flow=="no"]/total_Ret_Ret)*100
