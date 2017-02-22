
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
    legend.title=element_text(size=14, face="bold"))

########################
#Correlation of RCaMP and GCaMP Signals during long stim
CorrData <- read.table("D:/Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/LongStim_Correlations.csv", header=TRUE, sep = ",")

CorrData$CompChannel<-paste(CorrData$ChannelX, CorrData$ChannelY, sep= "_")
CorrData$CompType<-paste(CorrData$ROI_X_type, CorrData$ROI_Y_type, sep= "_")

ggplot(CorrData, aes(x=Long_Corr, fill=CompType)) + geom_histogram(binwidth=0.05, position="dodge") +
  ggtitle("Distribution of Correlation- All comparisons")

ggplot(CorrData, aes(x=Long_Corr, fill=CompChannel)) + geom_histogram(binwidth=0.05, position="dodge") +
  ggtitle("Distribution of Correlation- All comparisons")

ggplot(CorrData, aes(x=Short_Corr, fill=CompChannel)) + geom_histogram(binwidth=0.05, position="dodge") +
  ggtitle("Distribution of stim window Correlation- All comparisons")

ggplot(CorrData, aes(x=Lag, fill=CompChannel)) + geom_histogram(binwidth=0.5, position="dodge") +
  ggtitle("Distribution of lags- All comparisons")

ggplot(CorrData, aes(x=xCorr, fill=CompChannel)) + geom_histogram(binwidth=0.1, position="dodge") +
  ggtitle("Distribution of CrossCorrelations- All comparisons")


df1A<- summarySE(CorrData, measurevar="Short_Corr", groupvars=c("CompChannel"))
df1B<- summarySE(CorrData, measurevar="xCorr", groupvars=c("CompChannel"))
df1C<- summarySE(CorrData, measurevar="Lag", groupvars=c("CompChannel"))


GCaMP_RCaMP<-subset(CorrData, CompChannel=="GCaMP_RCaMP")

with(GCaMP_RCaMP, xyplot(GCaMP_RCaMP$xCorr ~ GCaMP_RCaMP$minDistance, 
                                   main="post30 vs post30 ",xlab="ROI Distance ", ylab="xCorrelation",))

##################################
#higher correlated data during stimulation window

highcorr<- subset(CorrData, Short_Corr>0.5)

df1A<- summarySE(highcorr, measurevar="Short_Corr", groupvars=c("CompType"))

ggplot(highcorr, aes(x=xCorr, fill=CompChannel)) + geom_histogram(binwidth=0.1, position="dodge") +
  ggtitle("Distribution of CrossCorrelations- All comparisons")

##################
df1A<- summarySE(WF.ROIs.Corr, measurevar="Corr_mean", groupvars=c("Condition"))

#####
ggplot(data=df1A, aes(x=Condition, y=Corr_mean, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Corr_mean-se, ymax=Corr_mean+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Condition") +
  ylab("Correlation") +
  ggtitle("Correlation for all ROIs in all trials") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

# box plot
ggplot(data=WF.ROIs.Corr, aes(x=Corr_mean, fill=Condition)) +
  geom_density(alpha=.2) 

ggplot(data=WF.ROIs.Corr, aes(x=Condition, y=Corr_mean, fill=Condition)) +
  geom_boxplot() 

# cumulative fraction plot
ggplot(WF.ROIs.Corr, aes(x=Corr_mean, colour = Condition)) + stat_ecdf()

######
#time lag
df2A<- summarySE(WF.ROIs.Corr, measurevar="Lag_mean", groupvars=c("Condition"))

#####
ggplot(data=df2A, aes(x=Condition, y=Lag_mean, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag_mean-se, ymax=Lag_mean+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Condition") +
  ylab("Time Lag (s)") +
  ggtitle("Time Lag for all ROIs in all trials") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

ggplot(data=WF.ROIs.Corr, aes(x=Condition, y=Lag_mean, fill=Condition)) +
  geom_boxplot() +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

#####
# Likelihood-ratio test 
# WF ROIs
Corr.WF.null = lmer(Corr_mean ~ (1|Animal) + (1|Spot), WF.ROIs.Corr,REML=FALSE)
Corr.WF.model1 = lmer(Corr_mean ~ Condition + (1|Animal) + (1|Spot), WF.ROIs.Corr,REML=FALSE)
Corr.WF.anova1 <- anova(Corr.WF.null, Corr.WF.model1)
print(Corr.WF.anova1)

# check residuals for linearity
plot(fitted(Corr.WF.model1), residuals(Corr.WF.model1),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(Corr.WF.model1), residuals(Corr.WF.model1)), col=46, lwd=2.5)

# p values
Corr.WF.stimpv <- glht(Corr.WF.model1, mcp(Condition= "Tukey"))
summary(Corr.WF.stimpv)


# Likelihood-ratio test 
# WF ROIs
Lag.WF.null = lmer(Lag_mean ~ (1|Animal) + (1|Spot), WF.ROIs.Corr,REML=FALSE)
Lag.WF.model1 = lmer(Lag_mean ~ Condition + (1|Animal) + (1|Spot), WF.ROIs.Corr,REML=FALSE)
Lag.WF.anova1 <- anova(Lag.WF.null, Lag.WF.model1)
print(Lag.WF.anova1)

# check residuals for linearity
plot(fitted(Lag.WF.model1), residuals(Lag.WF.model1),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(Lag.WF.model1), residuals(Lag.WF.model1)), col=46, lwd=2.5)

# p values
Lag.WF.stimpv <- glht(Lag.WF.model1, mcp(Condition= "Tukey"))
summary(Lag.WF.stimpv)

##############
#Correlation of RCaMP and GCaMP Signals
AvsA<- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/AC_Corr.csv", header=TRUE, sep = ",")
NvsN <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/Neuron_Corr.csv", header=TRUE, sep = ",")
NvsA <- read.table("E:/Data/Two_Photon_Data/GCaMP_RCaMP/cyto_GCaMP6s/Results/NeuronAC_Corr.csv", header=TRUE, sep = ",")

####################
# population plots
AvsA$type<-"AvsA"

NvsA$type<-"NvsA"

NvsN$type<-"NvsN"

All = rbind(NvsA,AvsA,NvsN)
ggplot(All, aes(x=Corr, fill=Condition)) + geom_histogram(binwidth=0.05, position="dodge", xmin=1) +
  ggtitle("Distribution of Correlation- All comparisons")

ggplot(All, aes(x=TimeLag, fill=Condition)) + geom_histogram(binwidth=5, position="dodge", xmin=1) +
  ggtitle("Distribution of TimeLags- All comparisons")

##################################
# aggregate ROIs
all.ROIs.Corr<- ddply(All, c("Animal", "Spot", "Layer", "Condition","ROI_X", "ROI_Y","type"), summarise, 
                      Corr_mean = mean(Corr), Corr_SD = sd(Corr),
                      Dis_mean = mean(Distance), Dis_SD = sd(Distance),
                      Lag_mean = mean(TimeLag), Lag_SD = sd(TimeLag),
                      N= length(Corr))
all.ROIs.Corr$combROI<-paste(all.ROIs.Corr$ROI_X, all.ROIs.Corr$ROI_Y, sep= "_")

# compare the 90Hz and Nostim peak AUC for significantly greater 90Hz
#stim= subset(all.ROIs.Corr, Condition=="Stim")
#nostim= subset(all.ROIs.Corr, Condition=="Nostim")

#all.ROIs.diff<- stim
#all.ROIs.diff$Corr_diff <-stim$Corr_mean-nostim$Corr_mean

##################
df1A<- summarySE(all.ROIs.Corr, measurevar="Corr_mean", groupvars=c("Condition"))
df1B<- summarySE(all.ROIs.Corr, measurevar="Corr_mean", groupvars=c("Condition","type"))


#####
ggplot(data=df1A, aes(x=Condition, y=Corr_mean, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Corr_mean-se, ymax=Corr_mean+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Condition") +
  ylab("Correlation") +
  ggtitle("Correlation for all ROIs in all trials") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

ggplot(data=df1B, aes(x=type, y=Corr_mean, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Corr_mean-se, ymax=Corr_mean+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Condition") +
  ylab("Correlation") +
  ggtitle("Correlation for all ROIs in all trials") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

# box plot
ggplot(data=all.ROIs.Corr, aes(x=Corr_mean, fill=Condition)) +
  geom_density(alpha=.2) 

ggplot(data=all.ROIs.Corr, aes(x=Condition, y=Corr_mean, fill=Condition)) +
  geom_boxplot() 

# cumulative fraction plot
ggplot(all.ROIs.Corr, aes(x=Corr_mean, colour = Condition)) + stat_ecdf()
 
#####
# Likelihood-ratio test 
# All ROIs
all.ROIs.Corr$Cond_dis<- interaction(all.ROIs.Corr$Condition,all.ROIs.Corr$Dis_mean)
all.ROIs.Corr$Cond_type<- interaction(all.ROIs.Corr$Condition,all.ROIs.Corr$type)

Corr.all.null = lmer(Corr_mean ~ (1|Animal) + (1|Spot), all.ROIs.Corr,REML=FALSE)
Corr.all.model1 = lmer(Corr_mean ~ Condition + (1|Animal) + (1|Spot), all.ROIs.Corr,REML=FALSE)
Corr.all.model1B = lmer(Corr_mean ~ Dis_mean + (1|Animal) + (1|Spot), all.ROIs.Corr,REML=FALSE)
Corr.all.model2 = lmer(Corr_mean ~ Condition+Dis_mean + (1|Animal) + (1|Spot), all.ROIs.Corr,REML=FALSE)
Corr.all.model2B = lmer(Corr_mean ~ Cond_dis + (1|Animal) + (1|Spot), all.ROIs.Corr,REML=FALSE)
Corr.all.anova1 <- anova(Corr.all.null, Corr.all.model1B, Corr.all.model1, Corr.all.model2,Corr.all.model2B)
print(Corr.all.anova1)

Corr.all.model4= lmer(Corr_mean ~ Condition*type + (1|Animal) + (1|Spot), all.ROIs.Corr,REML=FALSE)
Corr.all.model4B= lmer(Corr_mean ~ type + (1|Animal) + (1|Spot), all.ROIs.Corr,REML=FALSE)
Corr.all.anova3 <- anova(Corr.all.null, Corr.all.model4B, Corr.all.model1, Corr.all.model4)
print(Corr.all.anova3)

#####
# p values
Corr.all.stimpv <- glht(Corr.all.model1, mcp(Condition= "Tukey"))
summary(Corr.all.stimpv)

Corr.all.stim_dis <- glht(Corr.all.model1, mcp(Condition+Dis_mean= "Tukey"))
summary(Corr.all.stim_dis)

Corr.all.stimpv  <- lsmeans(Corr.all.model1, pairwise ~ Condition, glhargs=list())
summary(Corr.all.stimpv)

######
#time lag
df2A<- summarySE(all.ROIs.Corr, measurevar="Lag_mean", groupvars=c("Condition"))
df2B<- summarySE(all.ROIs.Corr, measurevar="Lag_mean", groupvars=c("Condition","type"))

#####
ggplot(data=df2A, aes(x=Condition, y=Lag_mean, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag_mean-se, ymax=Lag_mean+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Condition") +
  ylab("Time Lag (s)") +
  ggtitle("Time Lag for all ROIs in all trials") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

ggplot(data=df2B, aes(x=type, y=Lag_mean, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Lag_mean-se, ymax=Lag_mean+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Condition") +
  ylab("Time Lag (s)") +
  ggtitle("Time Lag for all ROIs in all trials") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

ggplot(data=all.ROIs.Corr, aes(x=type, y=Lag_mean, fill=Condition)) +
  geom_boxplot() +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

#######
# Likelihood-ratio test 
# WF ROIs
Lag.null = lmer(Lag_mean ~ (1|Animal) + (1|Spot), all.ROIs.Corr,REML=FALSE)
Lag.model1 = lmer(Lag_mean ~ Condition + (1|Animal) + (1|Spot), all.ROIs.Corr,REML=FALSE)
Lag.anova1 <- anova(Lag.null, Lag.model1)
print(Lag.WF.anova1)

# check residuals for linearity
plot(fitted(Lag.model1), residuals(Lag.model1),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(Lag.model1), residuals(Lag.model1)), col=46, lwd=2.5)

# p values
Lag.stimpv <- glht(Lag.model1, mcp(Condition= "Tukey"))
summary(Lag.stimpv)


# scatter plots of 
# scatter plots of distance vs lag
# scatter plots of 


#########################
#scatter plot of stimulation data

stim <- subset(all.ROIs.Corr, Condition=="Stim")
shortstim <- subset(all.ROIs.Corr, Condition=="ShortStim")

#distance vs correlation
with(stim, xyplot(stim$Corr_mean ~ stim$Dis_mean, group=type,
                  main="LongStim",xlab="ROI Distance ", ylab="Correlation mean ",  
                  auto.key=T))

with(shortstim, xyplot(shortstim$Corr_mean ~ shortstim$Dis_mean, group=type,
                  main="ShortStim",xlab="ROI Distance ", ylab="Correlation mean ",  
                  auto.key=T))


#lag vs correlation
with(stim, xyplot(stim$Corr_mean ~ stim$Lag_mean, group=type,
                  main="LongStim",xlab="Lag mean ", ylab="Correlation mean ",
                  auto.key=T))

with(shortstim, xyplot(shortstim$Corr_mean ~ shortstim$Lag_mean, group=type,
                  main="ShortStim",xlab="Lag mean ", ylab="Correlation mean ",
                  auto.key=T))

#lag vs distance
with(stim, xyplot(stim$Dis_mean ~ stim$Lag_mean, group=type,
                  main="LongStim",xlab="Lag mean ", ylab="Distance mean ",
                  auto.key=T))

with(shortstim, xyplot(shortstim$Dis_mean ~ shortstim$Lag_mean, group=type,
                  main="ShortStim",xlab="Lag mean ", ylab="Distance mean ",
                  auto.key=T))


###############
# sort data by highest correlation and closest distance

close<-subset(stim, Dis_mean<=20)
mid<- subset(stim, Dis_mean>20)
far<- subset(mid, Dis_mean>=40)
mid2<- subset(mid, Dis_mean<40)

close$Distance="Close"
mid2$Distance="Mid"
far$Distance="Far"

stim.dis<- rbind(close,mid2,far)

df3A1<- summarySE(stim.dis, measurevar="Corr_mean", groupvars=c("Distance"))
df3A2<- summarySE(stim.dis, measurevar="Corr_mean", groupvars=c("Distance","type"))

df3B1<- summarySE(stim.dis, measurevar="Lag_mean", groupvars=c("Distance"))
df3B2<- summarySE(stim.dis, measurevar="Lag_mean", groupvars=c("Distance","type"))




#
NvsA.close<- subset(stim, type=="NvsA" & Dis_mean<20 & Corr_mean>0.5)

with(NvsA.close, xyplot(NvsA.close$Dis_mean ~ NvsA.close$Lag_mean, group=type,
                       main="close NvsA",xlab="Lag mean ", ylab="Distance mean ",
                       auto.key=T))

with(NvsA.close, xyplot(NvsA.close$Corr_mean ~ NvsA.close$Lag_mean, group=type,
                        main="close NvsA",xlab="Lag mean ", ylab="Corr mean ",
                        auto.key=T))

with(NvsA.close, xyplot(NvsA.close$Corr_mean ~ NvsA.close$Dis_mean, group=type,
                        main="close NvsA",xlab="Dis mean ", ylab="Corr mean ",
                        auto.key=T))

df4A<- summarySE(NvsA.close, measurevar="Corr_mean")


# all post30 ROIs (with condition) that increase correlation with stimulation
high.diffpost30 <- unique(high.ROIs.diff.post30$combROI)
post30.ROIs.high <- subset(post30.ROIs, combROI %in% high.diffpost30)

#scatter plot
post30.ROIs.high.stim<-subset(post30.ROIs.high, Condition =="Stim")
with(post30.ROIs.high.stim, xyplot(post30.ROIs.high.stim$Corr_mean ~ post30.ROIs.high.stim$Dis_mean, group= type,
                                   main="post30 vs post30 ",xlab="ROI Distance ", ylab="Mean Correlation  ",))
####
df3A<- summarySE(post30.ROIs.high, measurevar="Corr_mean", groupvars=c("Condition"))
df3B<- summarySE(post30.ROIs.high, measurevar="Corr_mean", groupvars=c("Condition","type"))
df3C<- summarySE(post30.ROIs.high, measurevar="Corr_mean", groupvars=c("Condition","Layer"))
df3D<- summarySE(post30.ROIs.high, measurevar="Corr_mean", groupvars=c("Condition","Layer","type"))


#####
ggplot(data=df3A, aes(x=Condition, y=Corr_mean, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Corr_mean-se, ymax=Corr_mean+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Condition") +
  ylab("Correlation") +
  ggtitle("Correlation for post30 ROIs in positive") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

ggplot(data=df3B, aes(x=type, y=Corr_mean, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Corr_mean-se, ymax=Corr_mean+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Condition") +
  ylab("Correlation") +
  ggtitle("Correlation for post30 ROIs in positive") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme

ggplot(data=df3C, aes(x=Layer, y=Corr_mean, fill=Condition)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  geom_errorbar(aes(ymin=Corr_mean-se, ymax=Corr_mean+se), colour="black", width=.1,  position=position_dodge(.9)) +
  xlab("Condition") +
  ylab("Correlation") +
  ggtitle("Correlation for post30 ROIs in positive") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  science_theme





SomavsProc <-subset(post30.ROIs.high.stim, type=="SvsP")
with(SomavsProc, xyplot(SomavsProc$Corr_mean ~ SomavsProc$Dis_mean, group= type,
                                   main="SomavsProc ",xlab="ROI Distance ", ylab="Mean Correlation  ",))

# unique
test1 = setdiff(unique(high.ROIs.diff.post30x$ROI_X), unique(high.ROIs.diff.post30y$ROI_Y))
test2 = intersect(unique(high.ROIs.diff.post30x$ROI_X), unique(high.ROIs.diff.post30y$ROI_Y))
test3 = setdiff(unique(high.ROIs.diff.post30y$ROI_Y),unique(high.ROIs.diff.post30x$ROI_X))



post30vspost30 =length(high.ROIs.diff.post30xy$ROI_X)
numpost30 = length(unique(high.ROIs.diff.post30y$ROI_Y))

high.ROIs.diff.post30x_no30y <-subset(high.ROIs.diff.post30x, !(ROI_Y %in% post30names))  #post 30 X vs other ROIs y
post30vsother =length(high.ROIs.diff.post30x_no30y$ROI_X)

high.ROIs.diff.no30x <-subset(high.ROIs.diff, !(ROI_X %in% post30names))
high.ROIs.diff.no30x_30y <-subset(high.ROIs.diff.no30x,ROI_Y %in% post30names) #other ROIs X vs post 30 y
othervspost30 =length(high.ROIs.diff.no30x_30y$ROI_X)
high.ROIs.diff.no30x_no30y <-subset(high.ROIs.diff.no30x,!(ROI_Y %in% post30names)) #other ROIs X vs other ROIs y
othervsother =length(high.ROIs.diff.no30x_no30y $ROI_X)

# percentage of each ROI type in low NS, high S population
slices <- c(post30vspost30,post30vsother, othervspost30,othervsother) 
pct <- round(slices/sum(slices)*100)







# all ROIs
with(all.ROIs, xyplot(all.ROIs$Corr_mean ~ all.ROIs$Dis_mean, group=type,
                      main="all.ROIs",xlab="ROI Distance ", ylab="Correlation ",))

#######
# distance between somas and all other signals
# distance between all other signals

# ROIs subsetted by comparison type and considering correlation difference

all.mean.diff$Disgroup<-0

all.mean.diff$Disgroup[all.mean.diff$Dis_mean<100]<-"near"
all.mean.diff$Disgroup[all.mean.diff$Dis_mean>=100]<-"far"

df.dis1 <- summarySE(all.mean.diff, measurevar="Corr_diff", groupvars=c("Disgroup"))
df.dis2 <- summarySE(all.mean.diff, measurevar="Corr_diff", groupvars=c("Disgroup","type"))
df.dis3 <- summarySE(all.mean.diff, measurevar="Corr_diff", groupvars=c("Disgroup","type","Layer"))

# Likelihood-ratio test 
# All ROIs
all.mean.diff$type<- as.factor(all.mean.diff$type)
all.mean.diff$Disgroup<- as.factor(all.mean.diff$Disgroup)
all.mean.diff$Disgroup_type <- interaction(all.mean.diff$Disgroup, all.mean.diff$type)
all.mean.diff$Disgroup_type_Lay <- interaction(all.mean.diff$Disgroup, all.mean.diff$type,all.mean.diff$Layer)
all.mean.diff$Disgroup_type_type2 <- interaction(all.mean.diff$Disgroup, all.mean.diff$type,all.mean.diff$type2)

Corr.diff.dis.null = lmer(Corr_diff ~ (1|Animal) + (1|Spot), all.mean.diff,REML=FALSE)
Corr.diff.dis.model1 = lmer(Corr_diff ~ Disgroup + (1|Animal) + (1|Spot), all.mean.diff,REML=FALSE)
Corr.diff.dis.model2A = lmer(Corr_diff ~ Disgroup + type + (1|Animal) + (1|Spot), all.mean.diff,REML=FALSE)
Corr.diff.dis.model2B = lmer(Corr_diff ~ Disgroup_type  + (1|Animal) + (1|Spot), all.mean.diff,REML=FALSE)
Corr.diff.dis.model3A = lmer(Corr_diff ~ Disgroup+type+Layer + (1|Animal) + (1|Spot), all.mean.diff,REML=FALSE)
Corr.diff.dis.model3B = lmer(Corr_diff ~ Disgroup_type_Lay + (1|Animal) + (1|Spot), all.mean.diff,REML=FALSE)
Corr.diff.dis.anova1 <- anova(Corr.diff.dis.null, Corr.diff.dis.model1,Corr.diff.dis.model2A,Corr.diff.dis.model2B,
                              Corr.diff.dis.model3A,Corr.diff.dis.model3B)
print(Corr.diff.dis.anova1)

# p values
#Corr.diff.dis.pvalue <- lsmeans(Corr.diff.dis.model2B, pairwise ~ Disgroup*type, glhargs=list())
#summary(Corr.diff.dis.pvalue)

posthoc3 <- glht(Corr.diff.dis.model2B, mcp(Disgroup_type= "Tukey"))
summary(posthoc3)

posthoc4 <- glht(Corr.diff.dis.model3B, mcp(Disgroup_type_Lay= "Tukey"))
summary(posthoc4)

ggplot(data=all.mean.diff, aes(x=Disgroup, y=Corr_diff, fill=type)) +
  geom_boxplot() 

ggplot(data=all.mean.diff, aes(x=interaction(Disgroup,Layer), y=Corr_diff, fill=type)) +
  geom_boxplot() 


#####
# check residuals for linearity
plot(fitted(Corr.diff.all.model1), residuals(Corr.diff.all.model1),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(Corr.diff.all.model1), residuals(Corr.diff.all.model1)), col=46, lwd=2.5)


#####
# only responding ROIs
post30.mean.diff<- subset(all.mean.diff, type2=="P30vsP30")

df.dis4 <- summarySE(post30.mean.diff, measurevar="Corr_diff", groupvars=c("Disgroup"))
df.dis5 <- summarySE(post30.mean.diff, measurevar="Corr_diff", groupvars=c("Disgroup","type"))
df.dis6 <- summarySE(post30.mean.diff, measurevar="Corr_diff", groupvars=c("Disgroup","type","Layer"))

# Likelihood-ratio test 
# All ROIs
post30.mean.diff$type<- as.factor(post30.mean.diff$type)
post30.mean.diff$Disgroup<- as.factor(post30.mean.diff$Disgroup)
post30.mean.diff$Disgroup_type <- interaction(post30.mean.diff$Disgroup, post30.mean.diff$type)
post30.mean.diff$Disgroup_type_Lay <- interaction(post30.mean.diff$Disgroup, post30.mean.diff$type,post30.mean.diff$Layer)

post30.diff.dis.null = lmer(Corr_diff ~ (1|Animal) + (1|Spot), post30.mean.diff,REML=FALSE)
post30.diff.dis.model1 = lmer(Corr_diff ~ Disgroup + (1|Animal) + (1|Spot), post30.mean.diff,REML=FALSE)
post30.diff.dis.model2A = lmer(Corr_diff ~ Disgroup + type + (1|Animal) + (1|Spot), post30.mean.diff,REML=FALSE)
post30.diff.dis.model2B = lmer(Corr_diff ~ Disgroup_type  + (1|Animal) + (1|Spot), post30.mean.diff,REML=FALSE)
post30.diff.dis.model3A = lmer(Corr_diff ~ Disgroup+type+Layer + (1|Animal) + (1|Spot), post30.mean.diff,REML=FALSE)
post30.diff.dis.model3B = lmer(Corr_diff ~ Disgroup_type_Lay + (1|Animal) + (1|Spot), post30.mean.diff,REML=FALSE)
post30.diff.dis.anova1 <- anova(post30.diff.dis.null, post30.diff.dis.model1,post30.diff.dis.model2A,post30.diff.dis.model2B,
                                post30.diff.dis.model3A,post30.diff.dis.model3B)
print(post30.diff.dis.anova1)

posthoc5B <- glht(post30.diff.dis.model1, mcp(Disgroup= "Tukey"))
summary(posthoc5B)

posthoc5 <- glht(post30.diff.dis.model2B, mcp(Disgroup_type= "Tukey"))
summary(posthoc5)

posthoc6 <- glht(post30.diff.dis.model3B, mcp(Disgroup_type_Lay= "Tukey"))
summary(posthoc6)

ggplot(data=post30.mean.diff, aes(x=Disgroup, y=Corr_diff, fill=type)) +
  geom_boxplot() 

ggplot(data=post30.mean.diff, aes(x=interaction(Disgroup,Layer), y=Corr_diff, fill=type)) +
  geom_boxplot() 

############################
#load data with info about active ROIs and trials


sigROIs<- read.table("E:/Data/Two_Photon_Data/GFAP_GCaMP6/Whisker_stim/Longstim/Results/sigROIs.csv", header=TRUE, sep = ",")



# the ROI_Y that are NOT post30 ROIs
nonpost30 <- setdiff(post30.all$ROI_Y,post30.all2$ROI_Y)
nonpost30.all<- data.frame()
for (xx in 1:length(nonpost30))
{
  name2B = nonpost30[xx]
  ROIsubset2B = subset(post30.all, ROI_Y == name2B)
  nonpost30.all<- rbind(nonpost30.all,ROIsubset2B)
}

# aggregate trials
post30.ROIs<- ddply(post30.all2, c("Animal", "Spot", "Layer", "Condition","ROI_X", "ROI_Y","type"), summarise, 
                    Corr_mean = mean(Corr), Corr_SD = sd(Corr),
                    Dis_mean = mean(Distance), Dis_SD = sd(Distance),
                    N= length(Corr))

nonpost30.ROIs<- ddply(nonpost30.all, c("Animal", "Spot", "Layer", "Condition","ROI_X", "ROI_Y","type"), summarise, 
                       Corr_mean = mean(Corr), Corr_SD = sd(Corr),
                       Dis_mean = mean(Distance), Dis_SD = sd(Distance),
                       N= length(Corr))

#histogram of post30 ROIs
post30.all2<-subset(post30.all2, Corr<1)  #remove comparisons to themselves

ggplot(post30.all2, aes(x=Corr, fill=Condition)) + geom_histogram(binwidth=0.05, position="dodge", xmin=1) +
  ggtitle("Distribution of Correlation- All post30 ROIs and all trials")

post30.ROIs<-subset(post30.ROIs, Corr_mean<1)  #remove comparisons to themselves
ggplot(post30.ROIs, aes(x=Corr_mean, fill=Condition)) + geom_histogram(binwidth=0.05, position="dodge", xmin=1) +
  ggtitle("Distribution of Mean Correlation- All post30 ROIs and all trials")

nonpost30.ROIs<-subset(nonpost30.ROIs, Corr_mean<1)  #remove comparisons to themselves
ggplot(nonpost30.ROIs, aes(x=Corr_mean, fill=Condition)) + geom_histogram(binwidth=0.05, position="dodge", xmin=1) +
  ggtitle("Distribution of Mean Correlation-non post30 ROIs")

######
# sort out the significant ROIs
sigROI<-as.character(unique(sigROIs$ROI))

sigROIs.all<- data.frame()
# find all entries compared to other post30 ROIs
for (ii in 1:length(sigROI))
{
  name4 = sigROI[ii]
  ROIsubset4 = subset(All, ROI_X == name4)
  sigROIs.all<- rbind(sigROIs.all,ROIsubset4)
}

sigROIs.all2 <- data.frame()
# find all entries compared to other post30 ROIs
for (xx in 1:length(post30names))
{
  name5 = sigROI[xx]
  ROIsubset6 = subset(sigROIs.all, ROI_Y == name5)
  sigROIs.all2<- rbind(sigROIs.all2,ROIsubset6)
}

#histogram of sig ROIs
sigROIs.all2<-subset(sigROIs.all2, Corr<1)  #remove comparisons to themselves

ggplot(sigROIs.all2, aes(x=Corr, fill=Condition)) + geom_histogram(binwidth=0.05, position="dodge", xmin=1) +
  ggtitle("Distribution of Correlation- All sigROIs and all trials")

# aggregate trials
sig.ROIs<- ddply(sigROIs.all2, c("Animal", "Spot", "Layer", "Condition","ROI_X", "ROI_Y","type"), summarise, 
                 Corr_mean = mean(Corr), Corr_SD = sd(Corr),
                 Dis_mean = mean(Distance), Dis_SD = sd(Distance),
                 N= length(Corr))

#histogram of sig ROIs
sig.ROIs<-subset(sig.ROIs, Corr_Mean<1)  #remove comparisons to themselves

ggplot(sig.ROIs, aes(x=Corr_mean, fill=Condition)) + geom_histogram(binwidth=0.05, position="dodge", xmin=1) +
  ggtitle("Distribution of Mean Correlation- All sigROIs and all trials")

###################
#sort out the active trials
post30_ActTrials<- read.table("E:/Data/Two_Photon_Data/GFAP_GCaMP6/Whisker_stim/Longstim/Results/Post30ROIs_activeTrials.csv", header=TRUE, sep = ",")
active.trials<-subset(post30_ActTrials, N2>0)

post30.activetrials<- data.frame()
# find all entries compared to other post30 ROIs
for (x in 1:length(active.trials$Trial))
{
  name3 = as.character(active.trials$ROI[x])
  trials = as.character(active.trials$Trial[x])
  ROIsubset3 = subset(post30.all2, ROI_X == name3 & Trial==trials)
  post30.activetrials<- rbind(post30.activetrials,ROIsubset3)
}

# aggregate trials
post30.ROIs.ActTrials<- ddply(post30.activetrials, c("Animal", "Spot", "Layer", "Condition","ROI_X", "ROI_Y","type"), summarise, 
                              Corr_mean = mean(Corr), Corr_SD = sd(Corr),
                              Dis_mean = mean(Distance), Dis_SD = sd(Distance),
                              N= length(Corr))

#histogram of post30 ROIs(as seeds for correlation) & active trials
ggplot(post30.ROIs.ActTrials, aes(x=Corr_mean, fill=Condition)) + geom_histogram(binwidth=0.05, position="dodge", xmin=1) +
  ggtitle("Distribution of Mean Correlation- post30 ROIs(as seeds for correlation) & active trials")

