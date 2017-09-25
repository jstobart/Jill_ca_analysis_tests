# Load the required packages.  You may need to install these.
library("lme4")
library("lmerTest")
library("lattice")
library("plyr")
library("ggplot2")
library("gplots")
library("lsmeans")
library("bear")
library("MASS")
library("multcomp")
library("reshape2")
library("data.table")
library("Hmisc")

max.theme <- theme_classic() + 
  theme(
    text=element_text(size=12),
    axis.text=element_text(size=12),
    axis.title=element_text(size=14, face="bold"),
    axis.title.y=element_text(vjust=1),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14, face="bold"))
labs.gt <- c("Ret/+", "Ret/Ret")

#########
# diameter, velocity data
fnCSV <- "E:/Data/Pericyte project/Two-photon data/Max/Baseline_BloodFlow/Results/Diam_vel_ret_ret.csv"
prr.raw <- read.csv(fnCSV)
prr.raw$Genotype="Ret/Ret"

#fix name of one entry
pr3=prr.raw$animal=="p"
prr.raw$animal[pr3]="prr3"

fnCSV <- "E:/Data/Pericyte project/Two-photon data/Max/Baseline_BloodFlow/Results/Diam_vel_ret_+.csv"
prt.raw <- read.csv(fnCSV)
prt.raw$Genotype="Ret/+"

# pulsatility index data
PICSV <- "E:/Data/Pericyte project/Two-photon data/Max/Baseline_BloodFlow/Results/PulsatilityIndex_ret_ret.csv"
PIret.raw <- read.csv(PICSV)
PIret.raw$Genotype="Ret/Ret"

#fix name of one entry
pr3=PIret.raw$animal=="p"
PIret.raw$animal[pr3]="prr3"

PICSV <- "E:/Data/Pericyte project/Two-photon data/Max/Baseline_BloodFlow/Results/PulsatilityIndex_ret_+.csv"
PIwt.raw <- read.csv(PICSV)
PIwt.raw$Genotype="Ret/+"

data<- rbind(prr.raw,prt.raw)
PIdata<- rbind(PIret.raw,PIwt.raw)


# make genotype an ordered factor
data$Genotype<- factor(data$Genotype,levels = c("Ret/+", "Ret/Ret"))
PIdata$Genotype<- factor(PIdata$Genotype,levels = c("Ret/+", "Ret/Ret"))
data$branch_order<- as.factor(data$branch_order)
PIdata$branch_order<- as.factor(PIdata$branch_order)

#unqiue vessel name
data$VesselName<- paste (data$animal, data$vessel, sep = "_", collapse = NULL)
PIdata$VesselName<- paste (PIdata$animal, PIdata$vessel, sep = "_", collapse = NULL)

# group branch orders into 3 groups (1st-3rd, 4th-6th, 7-9th)
data$branch_group=0
group1=subset(data, branch_order==1)
group1$branch_group="1"
group2=subset(data, branch_order==2)
group2$branch_group="2-3"
group3=subset(data, branch_order==3)
group3$branch_group="2-3"

group4=subset(data, branch_order==4)
group4$branch_group="4-5"
group5=subset(data, branch_order==5)
group5$branch_group="4-5"
group6=subset(data, branch_order==6)
group6$branch_group="6-7"
group7=subset(data, branch_order==7)
group7$branch_group="6-7"
group8=subset(data, branch_order==8)
group8$branch_group="8-9"
group9=subset(data, branch_order==9)
group9$branch_group="8-9"

all.data<-rbind(group1,group2,group3,group4,group5,group6,group7,group8,group9)
all.data$branch_group<- as.factor(all.data$branch_group)


PIdata$branch_group=0
group1=subset(PIdata, branch_order==1)
group1$branch_group="1"
group2=subset(PIdata, branch_order==2)
group2$branch_group="2-3"
group3=subset(PIdata, branch_order==3)
group3$branch_group="2-3"

group4=subset(PIdata, branch_order==4)
group4$branch_group="4-5"
group5=subset(PIdata, branch_order==5)
group5$branch_group="4-5"
group6=subset(PIdata, branch_order==6)
group6$branch_group="6-7"

group7=subset(PIdata, branch_order==7)
group7$branch_group="6-7"
group8=subset(PIdata, branch_order==8)
group8$branch_group="8-9"
group9=subset(PIdata, branch_order==9)
group9$branch_group="8-9"

all.PIdata<-rbind(group1,group2,group3,group4,group5,group6,group7,group8,group9)
all.PIdata$branch_group<- as.factor(all.PIdata$branch_group)

# remove weird data point with Infinte pulsatility index
all.PIdata<- subset(all.PIdata, PI!=Inf)

# outlier test
#source("http://goo.gl/UUyEzD")
#outlierKD(all.data, vel_radon)
#y


#outlierKD(all.data, diameter)
#outlierKD(all.data, flux)


# remove same outliers from PI index data
#outliers<-subset(all.data, is.na(vel_radon))#="NA")
#outliername<-outliers$VesselName
#all.PIdata2<-subset(all.PIdata, VesselName %!in% outliername)

# remove all data of velocity outliers
#all.data<-subset(all.data, vel_radon!="NA")


# Flux subset (only accurate fluxes included)
Flux<-subset(all.data, flux_note==1)



#####
#-----------------------------------------------------------------------------------------------------------------#
# GRAPHS
# diameter, velocity, flux, linear density

#histograms
ggplot(all.data, aes(x=diameter, fill=Genotype)) + geom_histogram(binwidth=1, position="dodge", xmin=1) +
  ggtitle("Distribution of diameters")
ggplot(all.data, aes(x=vel_radon, fill=Genotype)) + geom_histogram(binwidth=0.5, position="dodge", xmin=0.5) +
  ggtitle("Distribution of velocities")

ggplot(Flux, aes(x=flux, fill=Genotype)) + geom_histogram(binwidth=10, position="dodge", xmin=10) +
  ggtitle("Distribution of fluxes")

ggplot(Flux, aes(x=linear_density, fill=Genotype)) + geom_histogram(binwidth=10, position="dodge", xmin=10) +
  ggtitle("Distribution of linear density")

ggplot(Flux, aes(x=Hc, fill=Genotype)) + geom_histogram(binwidth=0.01, position="dodge", xmin=0.01) +
  ggtitle("Distribution of hematocrit")

ggplot(all.PIdata, aes(x=PI, fill=Genotype)) + geom_histogram(binwidth=0.2, position="dodge", xmin=0.2) +
  ggtitle("Distribution of pulsatility index")



## boxplots
#diameter
ggplot(all.data, aes(x = Genotype, y = diameter, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Diameter [um]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

#velocity
ggplot(all.data, aes(x = Genotype, y = vel_radon, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Velocity [mm/s]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

# flux
ggplot(Flux, aes(x = Genotype, y = flux, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Flux [RBCs/s]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


# linear density
ggplot(Flux, aes(x = Genotype, y = linear_density, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Linear Density [RBCs/mm]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

# hematocrit
ggplot(Flux, aes(x = Genotype, y = Hc, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Local hematocrit") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


# pulsatility index
ggplot(all.PIdata, aes(x = Genotype, y = PI, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Pulsatility Index") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


## scatterplot
# diameter vs velocity
ggplot(all.data, aes(x=diameter, y=vel_radon)) +
  geom_point(aes(fill = Genotype, colour=Genotype), shape = 1, size=3)+
  ggtitle("diameter vs velocity for all vessels") +
  xlab("Diameter [um]") + 
  ylab("Velocity [mm/s]") + 
  scale_colour_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

#diameter vs flux
ggplot(Flux, aes(x=diameter, y=flux)) +
  geom_point(aes(fill = Genotype, colour=Genotype), shape = 1, size=3)+
  ggtitle("diameter vs velocity for all vessels") +
  xlab("Diameter [um]") + 
  ylab("Flux [RBCs/s]") + 
  scale_colour_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

#depth plots
# scatterplot- diameter vs depth
ggplot(all.data, aes(x=diameter, y=depth)) +
  geom_point(aes(fill = Genotype, colour=Genotype), shape = 1, size=3)+
  ggtitle("diameter vs depth for all vessels") +
  xlab("Diameter [um]") + 
  ylab("Depth [um]") + 
  scale_colour_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme

ggplot(all.data, aes(x=vel_radon, y=depth)) +
  geom_point(aes(fill = Genotype, colour=Genotype), shape = 1, size=3)+
  ggtitle("velocity vs depth for all vessels") +
  xlab("Velocity [mm/s]") + 
  ylab("Depth [um]") + 
  scale_colour_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  max.theme


#####
# grouping by branch order

# drop group 1
all.data.groups<-subset(all.data, branch_group!="1")
Flux.groups<-subset(Flux, branch_group!="1")
all.PIdata.groups<-subset(all.PIdata, branch_group!="1")


# diameter
df1A<- summarySE(all.data.groups, measurevar="diameter", groupvars=c("Genotype"))
df1B<- summarySE(all.data.groups, measurevar="diameter", groupvars=c("Genotype","branch_group"))

ggplot(data=df1A, aes(x=Genotype, y=diameter, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=diameter-se, ymax=diameter+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Diameter [um]") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df1B, aes(x=branch_group, y=diameter, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=diameter-se, ymax=diameter+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Diameter [um]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme

# velocity
df2A<- summarySE(all.data.groups, measurevar="vel_radon", groupvars=c("Genotype"))
df2B<- summarySE(all.data.groups, measurevar="vel_radon", groupvars=c("Genotype","branch_group"))

ggplot(data=df2A, aes(x=Genotype, y=vel_radon, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=vel_radon-se, ymax=vel_radon+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Velocity [mm/s]") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df2B, aes(x=branch_group, y=vel_radon, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=vel_radon-se, ymax=vel_radon+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Velocity [mm/s]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme


# flux
df3A<- summarySE(Flux.groups, measurevar="flux", groupvars=c("Genotype"))
df3B<- summarySE(Flux.groups, measurevar="flux", groupvars=c("Genotype","branch_group"))

ggplot(data=df3A, aes(x=Genotype, y=flux, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=flux-se, ymax=flux+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Flux [RBCs/s]") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df3B, aes(x=branch_group, y=flux, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=flux-se, ymax=flux+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Flux [RBCs/s]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme


# linear density
df4A<- summarySE(Flux.groups, measurevar="linear_density", groupvars=c("Genotype"))
df4B<- summarySE(Flux.groups, measurevar="linear_density", groupvars=c("Genotype","branch_group"))

ggplot(data=df4A, aes(x=Genotype, y=linear_density, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=linear_density-se, ymax=linear_density+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Linear Density [RBCs/mm]") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df4B, aes(x=branch_group, y=linear_density, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=linear_density-se, ymax=linear_density+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Linear Density [RBCs/mm]") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme


# hematocrit
df5A<- summarySE(Flux.groups, measurevar="Hc", groupvars=c("Genotype"))
df5B<- summarySE(Flux.groups, measurevar="Hc", groupvars=c("Genotype","branch_group"))

ggplot(data=df5A, aes(x=Genotype, y=Hc, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Hc-se, ymax=Hc+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Local Hematocrit") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df5B, aes(x=branch_group, y=Hc, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Hc-se, ymax=Hc+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Local Hematocrit") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme


# pulsatiliy index (Vmax-Vmin)/Vmean
df6A<- summarySE(all.PIdata.groups, measurevar="PI", groupvars=c("Genotype"))
df6B<- summarySE(all.PIdata.groups, measurevar="PI", groupvars=c("Genotype","branch_group"))

ggplot(data=df6A, aes(x=Genotype, y=PI, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=PI-se, ymax=PI+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Pulsatility Index") +
  scale_fill_manual(
    values=c("black", "red"),guide=FALSE)+
  max.theme

ggplot(data=df6B, aes(x=branch_group, y=PI, fill=Genotype)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=PI-se, ymax=PI+se), colour="black", width=.1,  position=position_dodge(.9)) +
  ylab("Pulsatility Index") +
  scale_fill_manual(
    values=c("black", "red"))+
  max.theme



#####
# stats
# --------------------------------------------------------------------------- #
## Diameter alone

# Fit an empty linear model, and look at the output
diam.null<- lmer(diameter ~ (1|animal), all.data.groups, REML=FALSE)
summary(diam.null)
diam.model1 <- lmer(diameter ~ Genotype + (1|animal), all.data.groups, REML=FALSE)
summary(diam.model1)
diam.model2 <- lmer(diameter ~ Genotype + branch_group + (1|animal), all.data.groups, REML=FALSE)
summary(diam.model2)
diam.model3 <- lmer(diameter ~ Genotype*branch_group + (1|animal), all.data.groups, REML=FALSE)
summary(diam.model3)

# Look at the per-animal difference
dotplot(ranef(diam.null, condVar = TRUE))$animal
dotplot(ranef(diam.model1, condVar = TRUE))$animal
dotplot(ranef(diam.model2, condVar = TRUE))$animal
dotplot(ranef(diam.model3, condVar = TRUE))$animal

# Look at the model residuals
plot(diam.model2)
plot(diam.model3)

# Compare the difference between the null model with the alternative
# to see if there's an effect of Genotype
an01 <- anova(diam.null, diam.model1,diam.model2,diam.model3)
print(an01)

#p values
lsmeans1 <- lsmeans(diam.model1, pairwise ~ Genotype, glhargs=list())
summary(lsmeans1)
lsmeans2 <- lsmeans(diam.model3, pairwise ~ Genotype*branch_group, glhargs=list())
summary(lsmeans2)


# --------------------------------------------------------------------------- #
## Velocity alone

# Fit an empty linear model, and look at the output
vel.null<- lmer(vel_radon ~ (1|animal), all.data.groups, REML=FALSE)
summary(vel.null)
vel.model1 <- lmer(vel_radon ~ Genotype + (1|animal), all.data.groups, REML=FALSE)
summary(vel.model1)
vel.model2 <- lmer(vel_radon ~ Genotype + branch_group + (1|animal), all.data.groups, REML=FALSE)
summary(vel.model2)
vel.model3 <- lmer(vel_radon ~ Genotype*branch_group + (1|animal), all.data.groups, REML=FALSE)
summary(vel.model3)

# Look at the per-animal difference
dotplot(ranef(vel.null, condVar = TRUE))$animal
dotplot(ranef(vel.model1, condVar = TRUE))$animal
dotplot(ranef(vel.model2, condVar = TRUE))$animal
dotplot(ranef(vel.model3, condVar = TRUE))$animal

# Look at the model residuals
plot(vel.model2)
plot(vel.model3)

# Compare the difference between the null model with the alternative
# to see if there's an effect of Genotype
an02 <- anova(vel.null, vel.model1, vel.model2, vel.model3)
print(an02)

#p values
lsmeans3 <- lsmeans(vel.model1, pairwise ~ Genotype, glhargs=list())
summary(lsmeans3)
lsmeans4 <- lsmeans(vel.model3, pairwise ~ Genotype*branch_group, glhargs=list())
summary(lsmeans4)


# --------------------------------------------------------------------------- #
## flux alone

# Fit an empty linear model, and look at the output
flux.null<- lmer(flux ~ (1|animal), Flux.groups, REML=FALSE)
summary(flux.null)
flux.model1 <- lmer(flux ~ Genotype + (1|animal), Flux.groups, REML=FALSE)
summary(flux.model1)
flux.model2 <- lmer(flux ~ Genotype + branch_group + (1|animal), Flux.groups, REML=FALSE)
summary(flux.model2)
flux.model3 <- lmer(flux ~ Genotype*branch_group + (1|animal), Flux.groups, REML=FALSE)
summary(flux.model3)

# Look at the per-animal difference
dotplot(ranef(flux.null, condVar = TRUE))$animal
dotplot(ranef(flux.model1, condVar = TRUE))$animal
dotplot(ranef(flux.model2, condVar = TRUE))$animal
dotplot(ranef(flux.model3, condVar = TRUE))$animal

# Look at the model residuals
plot(flux.model2)
plot(flux.model3)

# Compare the difference between the null model with the alternative
# to see if there's an effect of Genotype
an03 <- anova(flux.null, flux.model1, flux.model2, flux.model3)
print(an03)

#p values
lsmeans5 <- lsmeans(flux.model1, pairwise ~ Genotype, glhargs=list())
summary(lsmeans5)
lsmeans6 <- lsmeans(flux.model3, pairwise ~ Genotype*branch_group, glhargs=list())
summary(lsmeans6)

# --------------------------------------------------------------------------- #
## linear density alone

# Fit an empty linear model, and look at the output
LD.null<- lmer(linear_density ~ (1|animal), Flux.groups, REML=FALSE)
summary(LD.null)
LD.model1 <- lmer(linear_density ~ Genotype + (1|animal), Flux.groups, REML=FALSE)
summary(LD.model1)
LD.model2 <- lmer(linear_density  ~ Genotype + branch_group + (1|animal), Flux.groups, REML=FALSE)
summary(LD.model2)
LD.model3 <- lmer(linear_density  ~ Genotype*branch_group + (1|animal), Flux.groups, REML=FALSE)
summary(LD.model3)

# Look at the per-animal difference
dotplot(ranef(LD.null, condVar = TRUE))$animal
dotplot(ranef(LD.model1, condVar = TRUE))$animal
dotplot(ranef(LD.model2, condVar = TRUE))$animal
dotplot(ranef(LD.model3, condVar = TRUE))$animal

# Look at the model residuals
plot(LD.model2)
plot(LD.model3)

# Compare the difference between the null model with the alternative
# to see if there's an effect of Genotype
an04 <- anova(LD.null, LD.model1, LD.model2, LD.model3)
print(an04)

#p values
lsmeans7 <- lsmeans(LD.model1, pairwise ~ Genotype, glhargs=list())
summary(lsmeans7)
lsmeans8 <- lsmeans(LD.model3, pairwise ~ Genotype*branch_group, glhargs=list())
summary(lsmeans8)

# --------------------------------------------------------------------------- #
## hematocrit alone

# Fit an empty linear model, and look at the output
Hc.null<- lmer(Hc ~ (1|animal), Flux.groups, REML=FALSE)
summary(Hc.null)
Hc.model1 <- lmer(Hc ~ Genotype + (1|animal), Flux.groups, REML=FALSE)
summary(Hc.model1)
Hc.model2 <- lmer(Hc~ Genotype + branch_group + (1|animal), Flux.groups, REML=FALSE)
summary(Hc.model2)
Hc.model3 <- lmer(Hc ~ Genotype*branch_group + (1|animal), Flux.groups, REML=FALSE)
summary(Hc.model3)

# Look at the per-animal difference
dotplot(ranef(Hc.null, condVar = TRUE))$animal
dotplot(ranef(Hc.model1, condVar = TRUE))$animal
dotplot(ranef(Hc.model2, condVar = TRUE))$animal
dotplot(ranef(Hc.model3, condVar = TRUE))$animal

# Look at the model residuals
plot(Hc.model2)
plot(Hc.model3)

# Compare the difference between the null model with the alternative
# to see if there's an effect of Genotype
an05 <- anova(Hc.null, Hc.model1, Hc.model2, Hc.model3)
print(an05)

#p values
lsmeans9 <- lsmeans(Hc.model1, pairwise ~ Genotype, glhargs=list())
summary(lsmeans9)
lsmeans10 <- lsmeans(Hc.model3, pairwise ~ Genotype*branch_group, glhargs=list())
summary(lsmeans10)


# --------------------------------------------------------------------------- #
## pulsatility index
# (Vmax-Vmin)/Vmean
# all data re-processed with 1000ms window time

# Fit an empty linear model, and look at the output
PI.null<- lmer(PI ~ (1|animal), all.PIdata.groups, REML=FALSE)
summary(PI.null)
PI.model1 <- lmer(PI ~ Genotype + (1|animal), all.PIdata.groups, REML=FALSE)
summary(PI.model1)
PI.model2 <- lmer(PI~ Genotype + branch_group + (1|animal), all.PIdata.groups, REML=FALSE)
summary(PI.model2)
PI.model3 <- lmer(PI ~ Genotype*branch_group + (1|animal), all.PIdata.groups, REML=FALSE)
summary(PI.model3)

# Look at the per-animal difference
dotplot(ranef(PI.null, condVar = TRUE))$animal
dotplot(ranef(PI.model1, condVar = TRUE))$animal
dotplot(ranef(PI.model2, condVar = TRUE))$animal
dotplot(ranef(PI.model3, condVar = TRUE))$animal

# Look at the model residuals
plot(PI.model2)
plot(PI.model3)

# Compare the difference between the null model with the alternative
# to see if there's an effect of Genotype
an05 <- anova(PI.null, PI.model1, PI.model2, PI.model3)
print(an05)

#p values
lsmeans11 <- lsmeans(PI.model1, pairwise ~ Genotype, glhargs=list())
summary(lsmeans11)
lsmeans12 <- lsmeans(PI.model3, pairwise ~ Genotype*branch_group, glhargs=list())
summary(lsmeans12)


#####
# branch order per depth???
# does it branch more in superficial layers?

# order vs depth
ggplot(all.data, aes(x=branch_order, y=depth)) +
  geom_point(aes(fill = Genotype, colour=Genotype), shape = 1, size=3)+
  xlab("branch order") + 
  ylab("depth [um]") + 
  scale_colour_manual(
    values=c("black", "red"))+ 
  max.theme

#boxplot
ggplot(all.data, aes(x = branch_order, y = depth, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("depth [um]") + 
  scale_fill_manual(
    values=c("black", "red"))
  max.theme

## branch order vs depth

# Fit an empty linear model, and look at the output
depth.b_o.Null <- lmer(depth ~ (1|animal), all.data, REML=FALSE)
summary(depth.b_o.Null)

# Fit an linear model with an effect , and look at the output
depth.b_o.model1 <- lmer(depth ~ branch_order + (1|animal), all.data, REML=FALSE)
summary(depth.b_o.model1)

# Fit an linear model with an effect , and look at the output
depth.b_o.model2 <- lmer(depth ~ branch_order + Genotype + (1|animal), all.data, REML=FALSE)
summary(depth.b_o.model2)

# Fit an linear model with an effect , and look at the output
depth.b_o.model3 <- lmer(depth ~ branch_order*Genotype + (1|animal), all.data, REML=FALSE)
summary(depth.b_o.model3)

# Look at the per-animal difference
dotplot(ranef(depth.b_o.model3, condVar = TRUE))$animal
# qqmath(ranef(fm.vel.diam.Alt06, condVar = TRUE))$animal # an alternative view

# Look at the model residuals
plot(depth.b_o.model3)

# Compare the difference between the null model with the alternative
# to see if there's an effect of cell type
depth.b_o.anova <- anova(depth.b_o.Null,depth.b_o.model1,depth.b_o.model2,depth.b_o.model3)
print(depth.b_o.anova)

lsmeans13 <- lsmeans(depth.b_o.model3, pairwise ~ Genotype*branch_order, glhargs=list())
summary(lsmeans13)


##########

# Matt's old code......
ggplot(dat.raw, aes(x = Genotype, y = Diameter, fill = Genotype)) + 
  geom_boxplot() + 
  ylab("Diameter [um]") + 
  scale_fill_manual(
    values=c("black", "red"), 
    guide=FALSE) + 
  scale_x_discrete(name="", labels=labs.gt) + 
  max.theme
ggsave(filename=file.path(dirFig, "diameter_box.png"), 
       dpi = 300, width=figWidthBox, height=figHeightBox)


ggplot(dat.raw, aes(x = Diameter, y = Velocity, colour = Genotype, shape = Genotype)) + 
  geom_point() + 
  ylim(0, yLim.vel) + 
  ylab("Velocity [mm/s]") + 
  xlab("Diameter [um]") + 
  scale_colour_manual(values=c("black", "red"), labels=labs.gt) +
  scale_shape_discrete(labels=labs.gt) + 
  max.theme + 
  theme(legend.justification=c(1,1), legend.position=c(1,1))
ggsave(filename=file.path(dirFig, "velocity_point.png"), 
       dpi = 300, width=figWidthPoint, height=figHeightBox)


# --------------------------------------------------------------------------- #
## Velocity vs Diameter

# Fit an empty linear model, and look at the output
fm.vel.diam.Null <- lmer(log(Velocity) ~ (1|Animal.ID), dat.raw, REML=FALSE)
summary(fm.vel.diam.Null)

# Fit an linear model with an effect , and look at the output
fm.vel.diam.Al01 <- lmer(log(Velocity) ~ Diameter + (1|Animal.ID), dat.raw, REML=FALSE)
summary(fm.vel.diam.Al01)

# Fit an linear model with an effect , and look at the output
fm.vel.diam.Alt02 <- lmer(log(Velocity) ~ Diameter + (Diameter|Animal.ID), dat.raw, REML=FALSE)
summary(fm.vel.diam.Alt02)

# Fit an linear model with an effect , and look at the output
fm.vel.diam.Alt03 <- lmer(log(Velocity) ~ Genotype + Diameter + (1|Animal.ID), dat.raw, REML=FALSE)
summary(fm.vel.diam.Alt03)

# Fit an linear model with an effect , and look at the output
fm.vel.diam.Alt04 <- lmer(log(Velocity) ~ Genotype*Diameter + (1|Animal.ID), dat.raw, REML=FALSE)
summary(fm.vel.diam.Alt04)

# Fit an linear model with an effect , and look at the output
fm.vel.diam.Alt05 <- lmer(log(Velocity) ~ Genotype + Diameter + (Diameter|Animal.ID), dat.raw, REML=FALSE)
summary(fm.vel.diam.Alt05)

# Fit an linear model with an effect , and look at the output
fm.vel.diam.Alt06 <- lmer(log(Velocity) ~ Genotype*Diameter + (Diameter|Animal.ID), dat.raw, REML=FALSE)
summary(fm.vel.diam.Alt06)

# Look at the per-animal difference
dotplot(ranef(fm.vel.diam.Alt06, condVar = TRUE))$Animal.ID
# qqmath(ranef(fm.vel.diam.Alt06, condVar = TRUE))$animal # an alternative view

# Look at the model residuals
plot(fm.vel.diam.Alt05)

# Compare the difference between the null model with the alternative
# to see if there's an effect of cell type
an.vel.diam.01 <- anova(fm.vel.diam.Null,  fm.vel.diam.Al01, fm.vel.diam.Alt02, 
                        fm.vel.diam.Alt03, fm.vel.diam.Alt04, fm.vel.diam.Alt05,
                        fm.vel.diam.Alt06)
print(an.vel.diam.01)

an.vel.diam.02 <- anova(fm.vel.diam.Null, fm.vel.diam.Al01, fm.vel.diam.Alt02, 
                        fm.vel.diam.Alt05)
print(an.vel.diam.02)


# --------------------------------------------------------------------------- #
## Summary

direct.comparisons <- c(
  sprintf("Diameter: %6.5f", an.diam.01$`Pr(>Chisq)`[2]),
  sprintf("Velocity: %6.5f", an.vel.01$`Pr(>Chisq)`[2]),
  sprintf("Flux: %6.5f", an.flux.01$`Pr(>Chisq)`[2]),
  sprintf("Linear Density: %6.5f", an.ld.01$`Pr(>Chisq)`[2]))

print(direct.comparisons)


print(an.vel.diam.02)

print(an.flux.diam.02)

print(an.ld.diam.01)
