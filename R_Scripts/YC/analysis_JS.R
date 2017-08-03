
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
data.raw <- read.table("E:/Data/Two_Photon_Data/YC_data/WTDN_allareas_matched_ROIs_all.csv", header=TRUE, sep=",") #
data.raw <- read.table("D:/Data/YC_data/WTDN_allareas_matched_ROIs_all.csv", header=TRUE, sep=",") #

# make NaNs zero
data.raw[is.na(data.raw)] <- 0

# Factors
data.raw$trial <- factor(data.raw$trial)
data.raw$area <- factor(data.raw$area)
data.raw$session <- factor(data.raw$session)

### data subsets
# by animals
data.raw.YC20 <- subset(data.raw, animal=="YC20")
data.raw.YC21 <- subset(data.raw, animal=="YC21")
data.raw.YC22 <- subset(data.raw, animal=="YC22")
data.raw.YC33 <- subset(data.raw, animal=="YC33")

# by area
data.raw.short.area1 <- subset(data.raw.short, area=="1")
data.raw.short.area2 <- subset(data.raw.short, area=="2")
data.raw.short.area3 <- subset(data.raw.short, area=="3")

data.raw.short.area2_3 <- rbind(data.raw.short.area2,data.raw.short.area3)

data.raw.short.area1.YC20 <- subset(data.raw.short.area1, animal=="YC20")
data.raw.short.area1.YC21 <- subset(data.raw.short.area1, animal=="YC21")
data.raw.short.area1.YC22 <- subset(data.raw.short.area1, animal=="YC22")

data.raw.short.area2.YC20 <- subset(data.raw.short.area2, animal=="YC20")
data.raw.short.area2.YC21 <- subset(data.raw.short.area2, animal=="YC21")
data.raw.short.area2.YC22 <- subset(data.raw.short.area2, animal=="YC22")

data.raw.short.area3.YC20 <- subset(data.raw.short.area3, animal=="YC20")
data.raw.short.area3.YC21 <- subset(data.raw.short.area3, animal=="YC21")
data.raw.short.area3.YC22 <- subset(data.raw.short.area3, animal=="YC22")


####################################
# correcting peakTime for YC20, YC21 and YC22 (take away all the baseline period)
data.raw.YC20.normPT <- subset(data.raw.YC20, session=="1" |session=="2" |session=="3"|session=="4")
data.raw.YC20.longPT <- subset(data.raw.YC20, session=="5" |session=="6")
data.raw.YC20.normPT$peakTimecor <- data.raw.YC20.normPT$peakTime-2.5
data.raw.YC20.longPT$peakTimecor <- data.raw.YC20.longPT$peakTime-3
data.raw.YC20 <- rbind(data.raw.YC20.normPT,data.raw.YC20.longPT)

data.raw.YC21.normPT <- subset(data.raw.YC21, session=="1" |session=="2" |session=="3"|session=="4")
data.raw.YC21.longPT <- subset(data.raw.YC21, session=="5" |session=="6")
data.raw.YC21.normPT$peakTimecor <- data.raw.YC21.normPT$peakTime-2.5
data.raw.YC21.longPT$peakTimecor <- data.raw.YC21.longPT$peakTime-3
data.raw.YC21 <- rbind(data.raw.YC21.normPT,data.raw.YC21.longPT)

data.raw.YC22$peakTimecor <- data.raw.YC22$peakTime-3
data.raw.YC33$peakTimecor <- data.raw.YC33$peakTime-2.5
#data.raw$peakTimecor <- data.raw$peakTime-2.5


data.raw2 <- rbind(data.raw.YC20,data.raw.YC21,data.raw.YC22,data.raw.YC33)

##########################################

# stimulating period plus 1s for each animal
data.raw.short.YC20 <- subset(data.raw.YC20, peakTimecor>=0 & peakTimecor<=2)
data.raw.short.YC21 <- subset(data.raw.YC21, peakTimecor>=0 & peakTimecor<=2)
data.raw.short.YC22 <- subset(data.raw.YC22, peakTimecor>=0 & peakTimecor<=2)
data.raw.short.YC33 <- subset(data.raw.YC33, peakTimecor>=0 & peakTimecor<=2)

#data.raw.short <- subset(data.raw, peakTimecor>=0 & peakTimecor<=2)
# combine subsets
data.raw.short <- rbind(data.raw.short.YC20,data.raw.short.YC21,data.raw.short.YC22,data.raw.short.YC33)


# Unique session name and roi name
#data.raw.short$sessionU <- interaction(data.raw.short$animal,data.raw.short$area, data.raw.short$session)
#data.raw.short$ROInameU <- interaction(data.raw.short$animal,data.raw.short$area, data.raw.short$ROIname)


# rough grouping based on amplitude
#data.raw.short.smallamp <- subset(data.raw.short, amplitude<=1)
#data.raw.short.midlamp <- subset(data.raw.short, amplitude>1 & amplitude<=3)
#data.raw.short.largelamp <- subset(data.raw.short, amplitude>3)



# aggregate the data

# How many peaks do rois give for each trial
sum.data.short.raw<- ddply(data.raw.short, c("animal","ROInameU","condition","treatment","area","session","trial"), summarise, 
                   PA_mean = mean(peakAUC), Amp_mean =mean(amplitude), nEvents = sum(numPeaks),
                   Dur_mean = mean(halfWidth), PT_mean = mean(peakTimecor))
# number of events per trial (use this for whole data)
sum.data.raw.Ev.whole<- ddply(data.raw, c("animal","ROInameU","condition","treatment","area","sessionU"), summarise, 
                           nEvents = sum(numPeaks))
# number of events per session (use this for short data)
sum.data.raw.Ev.short<- ddply(data.raw.short, c("animal","ROInameU","condition","treatment","area","sessionU"), summarise, 
                        nEvents = sum(numPeaks))

# sum across trials
sum.data.raw.roi<- ddply(sum.data.short.raw, c("animal","ROInameU","condition","treatment","area"), summarise, 
                     PA_mean2 = mean(PA_mean), Amp_mean2 =mean(Amp_mean), nEvents2 = sum(nEvents),
                     Dur_mean2 = mean(Dur_mean), PT_mean2 = mean(PT_mean))



##########
# --------------------------------------------------------------------------- #
## Do some exploratory plotting

# Data summarised (boxplot)
ggplot(data.raw.short, aes(x = animal, y = halfWidth, fill =session )) + 
  geom_boxplot() + 
  facet_grid(area~condition) + 
  theme_bw()+
  ggtitle("boxplot halfWidth duration short")

ggplot(data.raw.short, aes(x = area, y = amplitude, fill = condition)) + 
  geom_boxplot() + 
  facet_grid(treatment ~ condition) + 
  theme_bw()+
  ggtitle("boxplot amp short")

ggplot(data.raw.short, aes(x = area, y = halfWidth, fill = condition)) + 
  geom_boxplot() + 
  facet_grid(treatment ~ condition) + 
  theme_bw()+
  ggtitle("boxplot duration halfwidth short")

ggplot(data.raw.short, aes(x = area, y = fullWidth, fill = condition)) + 
  geom_boxplot() + 
  facet_grid(treatment~ condition) + 
  theme_bw()+
  ggtitle("boxplot duration fullWidth short")

ggplot(sum.data.raw.Ev.short, aes(x = sessionU, y = nEvents, fill = animal)) + 
  geom_boxplot() + 
  facet_grid(area~ condition) + 
  theme_bw()+
  ggtitle("boxplot nEvents short")


ggplot(data.raw.short, aes(x = condition, y = amplitude, colour = condition)) + 
  geom_boxplot() + 
  facet_grid(area ~ treatment) + 
  theme_bw()


######
# plots of means of all ROIs

ggplot(sum.data.raw.roi, aes(x = treatment, y = Amp_mean2, colour = condition)) + 
  geom_boxplot() + 
  facet_grid(~ area) + 
  theme_bw()

ggplot(sum.data.raw.roi, aes(x = Amp_mean2, fill = treatment)) + 
  geom_histogram(binwidth=.2, position="dodge") + 
  facet_grid(condition~ area) + 
  xlim(0,5)+
  theme_bw()


# plots of amplitudes

#before.stim<-subset(data.raw.short, treatment=="None"& condition=="Stim")

# one number for each ROI for each trial
#before.stim$trialU<-paste(before.stim$sessionU, before.stim$trial, sep=".")

#before.stim.ROI<- ddply(before.stim, c("animal","ROInameU","trialU","condition","treatment","area","session"), summarise, 
# Amp_mean =max(amplitude), nEvents = sum(numPeaks))


#write.csv(before.stim.ROI, file = "E:/Data/Two_Photon_Data/YC_data/before_stim_ROI.csv")

library(RColorBrewer)
jBuPuFun <- colorRampPalette(brewer.pal(n = 9, "BuPu"))
paletteSize <- 256
jBuPuPalette <- jBuPuFun(paletteSize)

## 'jet.colors' is "as in Matlab"
## (and hurting the eyes by over-saturation)
jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
jJet.colors <- jet.colors(paletteSize)

sum.data.raw.roi<- sum.data.raw.roi[order(sum.data.raw.roi$Amp_mean2),]

ggplot(sum.data.raw.roi[sum.data.raw.roi$area==1,], aes(x = condition, y = ROInameU, fill = Amp_mean2)) +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_tile() +
  facet_grid(~treatment) + 
  scale_fill_gradient2(low = jJet.colors[1],
                       mid = jJet.colors[paletteSize/2],
                       high = jJet.colors[paletteSize],
                       midpoint = max(sum.data.raw.roi$Amp_mean2) + min(sum.data.raw.roi$Amp_mean2)) / 2,
                       name = "Amplitude")



###########
# neurons that are high, mid, low responding in first 3 sessions (before tamoxifen)

# consider the groups for each Area independently

# area 1
area1.beforeTam.stim<-subset(sum.data.raw.roi, area==1 & treatment=="None"& condition=="Stim")

area1_percentiles<-quantile(area1.beforeTam.stim$Amp_mean2, prob = seq(0, 1, length = 21), type = 5)

area1.high<-area1.beforeTam.stim[area1.beforeTam.stim$Amp_mean2>area1_percentiles[20],]
area1.low<-area1.beforeTam.stim[area1.beforeTam.stim$Amp_mean2<area1_percentiles[11],]
area1.mid<-area1.beforeTam.stim[area1.beforeTam.stim$Amp_mean2<=area1_percentiles[20] & area1.beforeTam.stim$Amp_mean2>=area1_percentiles[11],]

# area 2
area2.beforeTam.stim<-subset(sum.data.raw.roi, area==2 & treatment=="None"& condition=="Stim")

area2_percentiles<-quantile(area2.beforeTam.stim$Amp_mean2, prob = seq(0, 1, length = 21), type = 5)

area2.high<-area2.beforeTam.stim[area2.beforeTam.stim$Amp_mean2>area2_percentiles[20],]
area2.low<-area2.beforeTam.stim[area2.beforeTam.stim$Amp_mean2<area2_percentiles[11],]
area2.mid<-area2.beforeTam.stim[area2.beforeTam.stim$Amp_mean2<=area2_percentiles[20] & area2.beforeTam.stim$Amp_mean2>=area2_percentiles[11],]


# area 3
area3.beforeTam.stim<-subset(sum.data.raw.roi, area==3 & treatment=="None"& condition=="Stim")

area3_percentiles<-quantile(area3.beforeTam.stim$Amp_mean2, prob = seq(0, 1, length = 21), type = 5)

area3.high<-area3.beforeTam.stim[area3.beforeTam.stim$Amp_mean2>area3_percentiles[20],]
area3.low<-area3.beforeTam.stim[area3.beforeTam.stim$Amp_mean2<area3_percentiles[11],]
area3.mid<-area3.beforeTam.stim[area3.beforeTam.stim$Amp_mean2<=area3_percentiles[20] & area3.beforeTam.stim$Amp_mean2>=area3_percentiles[11],]


high<-rbind(area1.high,area2.high,area3.high)
mid<-rbind(area1.mid,area2.mid,area3.mid)
low<-rbind(area1.low,area2.low,area3.low)


sum.data.raw.roi$responders=0

ROIdata<-data.frame()
for (ii in 1:nrow(high))
{
  xROI=high$ROInameU[ii]
  subset1=subset(sum.data.raw.roi, ROInameU==xROI)
  subset1$responders="high"
  ROIdata<-rbind(ROIdata, subset1) 
}

for (ii in 1:nrow(mid))
{
  xROI=mid$ROInameU[ii]
  subset1=subset(sum.data.raw.roi, ROInameU==xROI)
  subset1$responders="mid"
  ROIdata<-rbind(ROIdata, subset1) 
}

for (ii in 1:nrow(low))
{
  xROI=low$ROInameU[ii]
  subset1=subset(sum.data.raw.roi, ROInameU==xROI)
  subset1$responders="low"
  ROIdata<-rbind(ROIdata, subset1) 
}


ROIdata$responders<-factor(ROIdata$responders, levels=c("low","mid","high"))

ggplot(ROIdata, aes(x = responders, y = Amp_mean2, colour = condition, shape = condition)) + 
  geom_point(position = position_jitter(width = 0.2), size=3) + 
  facet_grid(area ~ treatment) + 
  theme_bw()

ggplot(ROIdata, aes(x =interaction(condition,responders), y = Amp_mean2, colour = condition, shape = condition)) + 
  geom_boxplot() + 
  facet_grid(area ~ treatment) + 
  theme_bw()

ggplot(ROIdata, aes(x =condition, y = Amp_mean2, colour = condition)) +   
  geom_point(size = 2) +
  geom_line(aes(group=ROInameU)) +
  facet_grid(area ~ interaction(treatment,responders)) + 
  theme_bw()

ggplot(ROIdata[ROIdata$responders=="high",], aes(x =treatment, y = Amp_mean2, colour = treatment)) +   
  geom_point(size = 2) +
  geom_line(aes(group=ROInameU)) +
  facet_grid(area ~ condition) + 
  ggtitle("high responders")+
  theme_bw()

ggplot(ROIdata[ROIdata$responders=="mid",], aes(x =treatment, y = Amp_mean2, colour = treatment)) +   
    geom_point(size = 2) +
    geom_line(aes(group=ROInameU)) +
    facet_grid(area ~ condition) + 
    ggtitle("mid responders")+
  theme_bw()
  
ggplot(ROIdata[ROIdata$responders=="low",], aes(x =treatment, y = Amp_mean2, colour = treatment)) +   
    geom_point(size = 2) +
    geom_line(aes(group=ROInameU)) +
    facet_grid(area ~ condition) + 
    ggtitle("low responders")+
  theme_bw()




######
# normalize based on no stim, no tamoxifen
diff.data.short.amp <- dcast(ROIdata, animal + ROInameU + area + responders   ~ treatment+condition, value.var="Amp_mean2")

diff.data.short.amp  <- diff.data.short.amp[complete.cases(diff.data.short.amp),]

diff.data.short.amp$norm_base_noTAM <- (diff.data.short.amp$None_Nostim/diff.data.short.amp$None_Nostim)
diff.data.short.amp$norm_Stim_noTAM <- (diff.data.short.amp$None_Stim/diff.data.short.amp$None_Nostim)
diff.data.short.amp$norm_nostim_TAM <- (diff.data.short.amp$TAM_Nostim/diff.data.short.amp$None_Nostim)
diff.data.short.amp$norm_Stim_TAM <- (diff.data.short.amp$TAM_Stim/diff.data.short.amp$None_Nostim)


diff.data.short.amp_1 <- diff.data.short.amp[c("animal","ROInameU","area","responders","norm_base_noTAM","norm_Stim_noTAM","norm_nostim_TAM","norm_Stim_TAM")]
diff.data.short.amp_2 <- melt(diff.data.short.amp_1)
diff.data.short.amp_2$logvalue <- log10(diff.data.short.amp_2$value)
  

# Normalized TAMStim to NoneStim
diff.data.short.amp$norm_nostim<- (diff.data.short.amp$TAM_Nostim/diff.data.short.amp$None_Nostim)
diff.data.short.amp$norm_stim <-(diff.data.short.amp$TAM_Stim/diff.data.short.amp$None_Stim)
  
diff.data.short.amp_3 <- diff.data.short.amp[c("animal","ROInameU","area","responders","norm_nostim","norm_stim")]
diff.data.short.amp_3 <- melt(diff.data.short.amp_3)
diff.data.short.amp_3$logvalue <- log10(diff.data.short.amp_3$value)


# plots
ggplot(diff.data.short.amp_2, aes(x = responders, y = value, fill = variable)) + 
  geom_boxplot() + 
  facet_grid(~area) + 
  theme_bw()+
  ggtitle("boxplot Amp")


ggplot(diff.data.short.amp_3, aes(x = responders, y = value, fill = variable)) + 
  geom_boxplot() + 
  facet_grid(~area) + 
  theme_bw()+
  ggtitle("boxplot Amp")


df.short.amp1<-summarySE(data=diff.data.short.amp_2, measurevar="value", groupvars=c("responders","area","variable"), na.rm=TRUE)
df.short.amp2<-summarySE(data=diff.data.short.amp_3, measurevar="value", groupvars=c("responders","area","variable"), na.rm=TRUE)

ggplot(df.short.amp1, aes(x = responders, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour="black", width=.6,  position=position_dodge(.9)) + 
  facet_grid(~ area) + 
  theme_bw()+
  ggtitle("bargraph amp- normalized to no stim, no TAM")

ggplot(df.short.amp2, aes(x = responders, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), colour="black", width=.6,  position=position_dodge(.9)) + 
  facet_grid(~ area) + 
  theme_bw()+
  ggtitle("bargraph amp- normalized to no TAM")



#########################
# Analysis and Stats

diff.data.short.amp_2$CondTreatResp <- interaction(diff.data.short.amp_2$variable, diff.data.short.amp_2$area, diff.data.short.amp_2$responders)

# Fit an empty linear model, and look at the output
diff.amp.Null <- lmer(value ~ (1|animal) + (1|ROInameU),diff.data.short.amp_2, REML=FALSE)
summary(diff.amp.Null)

# Fit an linear model with an effect of cond only, and look at the output
diff.amp.model1 <- lmer(value ~ variable + (1|animal) + (1|ROInameU), diff.data.short.amp_2, REML=FALSE)
summary(diff.amp.model1)

# Fit an linear model with an effect of cond only, and look at the output
diff.amp.model2 <- lmer(value ~ variable + area + (1|animal) + (1|ROInameU), diff.data.short.amp_2, REML=FALSE)
summary(diff.amp.model2)

# Fit an linear model with an effect of cond and area, and look at the output
diff.amp.model3<- lmer(value ~ variable*area + (1|animal) + (1|ROInameU), diff.data.short.amp_2, REML=FALSE)
summary(diff.amp.model3)

# Fit an linear model with an effect of cond and area and responders, and look at the output
diff.amp.model4<- lmer(value ~ variable+area+responders + (1|animal) + (1|ROInameU), diff.data.short.amp_2, REML=FALSE)
summary(diff.amp.model4)

# Fit an linear model with an effect of cond and area and responders, and look at the output
diff.amp.model5<- lmer(value ~ CondTreatResp + (1|animal) + (1|ROInameU), diff.data.short.amp_2, REML=FALSE)
summary(diff.amp.model5)

# ANOVA comparisons of models
an.diff.all <- anova(diff.amp.Null, diff.amp.model1, diff.amp.model2, diff.amp.model3,diff.amp.model4, diff.amp.model5)
print(an.diff.all)




# consider each area separately because there are too many comparisons for post hoc tests
area1.diff.amp<-subset(diff.data.short.amp_2, area==1)
area2.diff.amp<-subset(diff.data.short.amp_2, area==2)
area3.diff.amp<-subset(diff.data.short.amp_2, area==3)

area1.diff.amp$CondResp<-interaction(area1.diff.amp$variable, area1.diff.amp$responders)
area2.diff.amp$CondResp<-interaction(area2.diff.amp$variable, area2.diff.amp$responders)
area3.diff.amp$CondResp<-interaction(area3.diff.amp$variable, area3.diff.amp$responders)

# Fit an linear model with an effect of cond and responders, and look at the output
diff.amp.model6a<- lmer(value ~ CondResp + (1|animal) + (1|ROInameU), area1.diff.amp, REML=FALSE)
summary(diff.amp.model6a)

diff.amp.model6b<- lmer(value ~ CondResp + (1|animal) + (1|ROInameU), area2.diff.amp, REML=FALSE)
summary(diff.amp.model6b)

diff.amp.model6c<- lmer(value ~ CondResp + (1|animal) + (1|ROInameU), area3.diff.amp, REML=FALSE)
summary(diff.amp.model6c)


# Post-hoc test
diff.amp.pv.CondResp.A1 <- glht(diff.amp.model6a, mcp(CondResp="Tukey"))
summary(diff.amp.pv.CondResp.A1)


diff.amp.pv.CondResp.A2 <- glht(diff.amp.model6b, mcp(CondResp="Tukey"))
summary(diff.amp.pv.CondResp.A2)

diff.amp.pv.CondResp.A3 <- glht(diff.amp.model6c, mcp(CondResp="Tukey"))
summary(diff.amp.pv.CondResp.A3)



