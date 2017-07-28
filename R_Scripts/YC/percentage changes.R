##########################
# TAM-None for each ROI

data.raw.short$area <- factor(data.raw.short$area)

# Summary for other variables other than nEvents
diff.data.short<- ddply(data.raw.short, c("animal","ROInameU","condition","treatment","area"), summarise, 
                           PA_mean = mean(peakAUC), Amp_mean =mean(amplitude),
                           Dur_mean = mean(halfWidth), PT_mean = mean(peakTimecor))
## Summary nEvents only
sum.data.raw.Ev<- ddply(data.raw.short, c("animal","ROInameU","condition","treatment","area","sessionU"), summarise, 
                        nEvents = sum(numPeaks))
diff.data.short1 <- sum.data.raw.Ev
diff.data.short1 <- ddply(diff.data.short1, c("animal","ROInameU","condition","treatment","area"), summarise, 
                      Events_mean = mean(nEvents))

# Reshape
# Normalized to Nostim+None
diff.data.short.PT <- dcast(diff.data.short, animal + ROInameU + area   ~ treatment+condition, value.var="PT_mean")
diff.data.short.PT <- diff.data.short.PT[complete.cases(diff.data.short.PT),]
diff.data.short.PT$norm_base <- (diff.data.short.PT$None_Nostim/diff.data.short.PT$None_Nostim)
diff.data.short.PT$norm_Stim <- (diff.data.short.PT$None_Stim/diff.data.short.PT$None_Nostim)
diff.data.short.PT$norm_TAM <- (diff.data.short.PT$TAM_Nostim/diff.data.short.PT$None_Nostim)
diff.data.short.PT$norm_StimTAM <- (diff.data.short.PT$TAM_Stim/diff.data.short.PT$None_Nostim)
diff.data.short.PT_1 <- diff.data.short.PT[c("animal","ROInameU","area","norm_base","norm_Stim","norm_TAM","norm_StimTAM")]
diff.data.short.PT_2 <- melt(diff.data.short.PT_1)
diff.data.short.PT_2$logvalue <- log10(diff.data.short.PT_2$value)

# Normalized TAMStim to NoneStim
diff.data.short.PT_3 <- dcast(diff.data.short.PT_2, animal + ROInameU + area ~ variable, value.var="value")
diff.data.short.PT_3$diff_response <- (diff.data.short.PT_3$norm_StimTAM/diff.data.short.PT_3$norm_Stim)
diff.data.short.PT_3 <- diff.data.short.PT_3[c("animal","ROInameU","area","diff_response")]
diff.data.short.PT_3$log_diff_response <- log10(diff.data.short.PT_3$diff_response)

# Normalize to Nostim separately for before and after TAM
diff.data.short.PT_d1 <- diff.data.short.PT[c("animal","ROInameU","area","None_Nostim","None_Stim","TAM_Nostim","TAM_Stim")]
diff.data.short.PT_d1$norm_baseNone <- (diff.data.short.PT_d1$None_Nostim/diff.data.short.PT_d1$None_Nostim)
diff.data.short.PT_d1$norm_Stim <- (diff.data.short.PT_d1$None_Stim/diff.data.short.PT_d1$None_Nostim)
diff.data.short.PT_d1$norm_baseTAM <- (diff.data.short.PT_d1$TAM_Nostim/diff.data.short.PT_d1$TAM_Nostim)
diff.data.short.PT_d1$norm_StimTAM <- (diff.data.short.PT_d1$TAM_Stim/diff.data.short.PT_d1$TAM_Nostim)
diff.data.short.PT_d1 <- diff.data.short.PT_d1[c("animal","ROInameU","area","norm_baseNone","norm_Stim","norm_baseTAM","norm_StimTAM")]
diff.data.short.PT_d2 <- melt(diff.data.short.PT_d1)
diff.data.short.PT_d2$logvalue <- log10(diff.data.short.PT_d2$value)

# Find the ROIs with biggest difference

sort(diff.data.short.Amp_3$diff_response, decreasing = TRUE)

# Remove ROIs with na
diff.data.short.Ev_4 <- diff.data.short.Ev_3[complete.cases(diff.data.short.Ev_3),]



############################
## Check if there is any extreme outliers with boxplot

ggplot(diff.data.short.Amp_2, aes(x = area, y = value, fill = variable)) + 
  geom_boxplot() + 
  #facet_grid(treatment~condition) + 
  theme_bw()+
  ggtitle("boxplot Amp")

### Remove outliers

library(outliers)
outlier_tf = outlier(diff.data.short.Ev_5$diff2,logical=TRUE)
#This gives an array with all values False, except for the outlier (as defined in the package documentation "Finds value with largest difference between it and sample mean, which can be an outlier")
#That value isreturned as True. 
find_outlier = which(outlier_tf==TRUE,arr.ind=TRUE)
#This finds the location of the outlier by finding that "True" value within the "outlier_tf" array.
diff.data.short.Ev_5 = diff.data.short.Ev_5[-find_outlier,]
#This creates a new dataset based on the old data, removing the one row that contains the outlier

### Remove non-finite values
diff.data.short.Ev_5 <- diff.data.short.Ev_4
diff.data.short.Ev_5 <- diff.data.short.Ev_5[Reduce(`&`, lapply(diff.data.short.Ev_5, is.finite)),]

##########################
# Graph

df.short.PT <- summarySE(diff.data.short.PT_3, measurevar="log_diff_response", groupvars=c("area"), na.rm=TRUE)

ggplot(df.short.PT, aes(x = area, y = log_diff_response, fill = area)) + 
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_errorbar(aes(ymin=log_diff_response-se, ymax=log_diff_response+se), colour="black", width=.6,  position=position_dodge(.9)) + 
  #facet_grid(~ area) + 
  theme_bw()+
  ggtitle("bargraph log_PT difference in response after TAM")

#########################
# Analysis and Stats

diff.data.short.PT_d2$CondTreatAll <- interaction(diff.data.short.PT_d2$variable, diff.data.short.PT_d2$area)

# Fit an empty linear model, and look at the output
fm.diff.Null <- lmer(value ~ (1|animal) + (1|ROInameU),diff.data.short.PT_d2, REML=FALSE)
summary(fm.diff.Null)

# Fit an linear model with an effect of cond only, and look at the output
fm.diff.Alt01 <- lmer(value ~ variable + (1|animal) + (1|ROInameU), diff.data.short.PT_d2, REML=FALSE)
summary(fm.diff.Alt01)

# Fit an linear model with an effect of cond only, and look at the output
fm.diff.Alt02 <- lmer(value ~ variable + area + (1|animal) + (1|ROInameU), diff.data.short.PT_d2, REML=FALSE)
summary(fm.diff.Alt02)

# Fit an linear model with an effect of cond and area, and look at the output
fm.diff.Alt03 <- lmer(value ~ variable*area + (1|animal) + (1|ROInameU), diff.data.short.PT_d2, REML=FALSE)
summary(fm.diff.Alt03)

fm.diff.Alt03a <- lmer(value ~ CondTreatAll + (1|animal) + (1|ROInameU), diff.data.short.PT_d2, REML=FALSE)
summary(fm.diff.Alt03a)

# ANOVA comparisons of models
an.diff.all <- anova(fm.diff.Null, fm.diff.Alt01, fm.diff.Alt02, fm.diff.Alt03)
print(an.diff.all)

an.diff.03a <- anova(fm.diff.Alt03a, fm.diff.Alt03)
print(an.diff.03a)

# Post-hoc test
ph.diff.03a <- glht(fm.diff.Alt03a, mcp(CondTreatAll="Tukey"))
summary(ph.diff.03a)

### For double normalization #################
# Fit an empty linear model, and look at the output  + (1|ROInameU)
fm.diff.Null <- lmer(diff_response ~ (1|animal),diff.data.short.PT_3, REML=FALSE)
summary(fm.diff.Null)

fm.diff.Alt01 <- lmer(diff_response ~ area + (1|animal),diff.data.short.PT_3, REML=FALSE)
summary(fm.diff.Alt01)

an.diff.01 <- anova(fm.diff.Null, fm.diff.Alt01)
print(an.diff.01)

ph.diff.01 <- glht(fm.diff.Alt01, mcp(area="Tukey"))
summary(ph.diff.01)
