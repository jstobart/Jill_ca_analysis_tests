
# Load some extra 'packages' that we need for the analysis
library("ggplot2") # plotting
library("lme4") # linear mixed effects models (multilevel analysis)
library("lattice") # plotting
library("multcomp") # post hoc testing
library("lsmeans") # post hoc testing
library("plyr")
library("Rmisc")
library("reshape2")
library("data.table")

# Load the raw data
data.raw <- read.table("D:/Data/YC_data/WTKI_raw.csv", header=TRUE, sep=",") #
#data.raw <- read.table("E:/Data/Two_Photon_Data/YC_data/WTDN_allareas_matched_ROIs_all.csv", header=TRUE, sep=",") #
#data.raw <- read.table("D:/Data/YC_data/WTDN_allareas_matched_ROIs_all.csv", header=TRUE, sep=",") #

# make NaNs zero
#data.raw[is.na(data.raw)] <- 0

# Factors
data.raw$trial <- factor(data.raw$trial)
data.raw$area <- factor(data.raw$area)
data.raw$session <- factor(data.raw$session)


data.raw$ROI_area<-paste(data.raw$roiName, data.raw$area, sep="_")

data.raw.short<-subset(data.raw, peakTime<5)
#data.raw.short<-data.raw
##########################################

# aggregate the data (first 5 s)

# How many peaks do rois give for each trial
sum.data.short.raw<- ddply(data.raw.short, c("animal","ROI_area","condition","genotype","session"), summarise, 
                   PA_mean = mean(peakAUC, na.rm = TRUE), Amp_mean =mean(amplitude,na.rm = TRUE), nEvents = sum(numPeaks),
                   Dur_mean = mean(halfWidth,na.rm = TRUE), PT_mean = mean(peakTime,na.rm = TRUE), ntrials=40)

sum.data.short.raw$Freq<-sum.data.short.raw$nEvents/sum.data.short.raw$ntrials

# sum across trials
sum.data.raw.roi<- ddply(sum.data.short.raw, c("animal","ROI_area","condition","genotype"), summarise, 
                         PA_mean2 = mean(PA_mean), Amp_mean2 =mean(Amp_mean), nEvents2 = sum(nEvents),
                         Dur_mean2 = mean(Dur_mean), PT_mean2 = mean(PT_mean), Freq_mean=mean(Freq))

##########
# --------------------------------------------------------------------------- #
## Do some exploratory plotting  ALL DATA

# Data summarised (boxplot)
ggplot(data.raw.short, aes(x = animal, y = halfWidth, fill =session )) + 
  geom_boxplot() + 
  facet_grid(area~condition) + 
  theme_bw()+
  ggtitle("boxplot halfWidth duration short")

ggplot(data.raw.short, aes(x = animal, y = peakAUC, fill =session )) + 
  geom_boxplot() + 
  facet_grid(area~condition) + 
  theme_bw()+
  ggtitle("boxplot peakAUC short")

ggplot(data.raw.short, aes(x = area, y = amplitude, fill = condition)) + 
  geom_boxplot() + 
  facet_grid(genotype ~ condition) + 
  theme_bw()+
  ggtitle("boxplot amp short")

ggplot(data.raw.short, aes(x = area, y = halfWidth, fill = condition)) + 
  geom_boxplot() + 
  facet_grid(genotype ~ condition) + 
  theme_bw()+
  ggtitle("boxplot duration halfwidth short")

ggplot(data.raw.short, aes(x = area, y = fullWidth, fill = condition)) + 
  geom_boxplot() + 
  facet_grid(genotype~ condition) + 
  theme_bw()+
  ggtitle("boxplot duration fullWidth short")



ggplot(data.raw.short, aes(x = condition, y = amplitude, colour = condition)) + 
  geom_boxplot() + 
  facet_grid(area ~ genotype) + 
  theme_bw()


######
# plots of means of all ROIs
sum.data.raw.roi$genotype<- factor(sum.data.raw.roi$genotype, levels= c("WT","KI"))
ggplot(sum.data.raw.roi, aes(x = genotype, y = Amp_mean2, colour = condition)) + 
  geom_boxplot() + 
  ylab("Mean Amplitude") + 
  theme_bw()

ggplot(sum.data.raw.roi, aes(x = genotype, y = PA_mean2, colour = condition)) + 
  geom_boxplot() + 
  ylab("Mean Peak AUC") + 
  theme_bw()

ggplot(sum.data.raw.roi, aes(x = genotype, y = Dur_mean2, colour = condition)) + 
  geom_boxplot() + 
  ylab("Mean Duration (s)") + 
  theme_bw()


ggplot(sum.data.raw.roi, aes(x = genotype, y = Freq_mean, colour = condition)) + 
  geom_boxplot() + 
  ylab("Mean Number of Events/session") + 
  theme_bw()

###########
# Stats

#amplitude
amp.null = lmer(Amp_mean2 ~ (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)
amp.model1 = lmer(Amp_mean2 ~ condition + (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)
amp.model2 = lmer(Amp_mean2 ~ genotype + (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)
amp.model3 = lmer(Amp_mean2 ~ condition + genotype + (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)
amp.model4 = lmer(Amp_mean2 ~ condition * genotype + (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)

amp.anova <- anova(amp.null, amp.model1, amp.model2, amp.model3, amp.model4)
print(amp.anova)

# p values
amp.pvalues <- lsmeans(amp.model4, pairwise ~ condition*genotype, glhargs=list())
summary(amp.pvalues)


#peak AUC
peakauc.null = lmer(PA_mean2 ~ (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)
peakauc.model1 = lmer(PA_mean2 ~ condition + (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)
peakauc.model2 = lmer(PA_mean2 ~ genotype + (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)
peakauc.model3 = lmer(PA_mean2 ~ condition + genotype + (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)
peakauc.model4 = lmer(PA_mean2 ~ condition * genotype + (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)

peakauc.anova <- anova(peakauc.null, peakauc.model1, peakauc.model2, peakauc.model3, peakauc.model4)
print(peakauc.anova)

# p values
peakauc.pvalues <- lsmeans(peakauc.model4, pairwise ~ condition*genotype, glhargs=list())
summary(peakauc.pvalues)


#duration
dur.null = lmer(Dur_mean2 ~ (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)
dur.model1 = lmer(Dur_mean2 ~ condition + (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)
dur.model2 = lmer(Dur_mean2 ~ genotype + (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)
dur.model3 = lmer(Dur_mean2 ~ condition + genotype + (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)
dur.model4 = lmer(Dur_mean2 ~ condition * genotype + (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)

dur.anova <- anova(dur.null, dur.model1, dur.model2, dur.model3, dur.model4)
print(dur.anova)

# p values
dur.pvalues <- lsmeans(dur.model4, pairwise ~ condition*genotype, glhargs=list())
summary(dur.pvalues)


#mean number of events per session (divided by the number of trials per session)
neps.null = lmer(Freq_mean ~ (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)
neps.model1 = lmer(Freq_mean~ condition + (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)
neps.model2 = lmer(Freq_mean ~ genotype + (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)
neps.model3 = lmer(Freq_mean ~ condition + genotype + (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)
neps.model4 = lmer(Freq_mean ~ condition * genotype + (1|animal) + (1|ROI_area), sum.data.raw.roi,REML=FALSE)

neps.anova <- anova(neps.null, neps.model1, neps.model2, neps.model3, neps.model4)
print(neps.anova)

# p values
neps.pvalues <- lsmeans(neps.model4, pairwise ~ condition*genotype, glhargs=list())
summary(neps.pvalues)




