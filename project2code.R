#Clear environment
rm(list=ls())

################################################################################
################################################################################

#Load libraries 

################################################################################
################################################################################

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(onewaytests) #For brown-forsythe test
library(forcats)
library(ggpubr)
library(gridExtra)
library(lmPerm)
library(jmuOutlier)
library(xtable)

################################################################################
################################################################################

#Load data

################################################################################
################################################################################

##Load interferon data
interferon=read_excel("Project 2-1.xlsx",
                     sheet=1)

##Load ohd data
ohd=read_excel("Project 2-1.xlsx",
               sheet=2)

##Load lactoferrin data
lactoferrin=read_excel("Project 2-1.xlsx",
                       sheet=3)

##Load immunoglobulyn data
immunoglobulyn=read_excel("Project 2-1.xlsx",
                          sheet=4)

##Load level data
level=read_excel("Project 2-2.xlsx",
                 sheet=4)

################################################################################
################################################################################

#Exploratory Data Analysis

################################################################################
################################################################################

#Interferon

#Convert data from wide to long format
interferon=interferon %>%
  pivot_longer(cols=TIME1:TIME4,
               names_to = "TimeType",
               values_to = "TimeValue")

#Box plot
ggsummarystats(
  interferon,
  x = "TRT",
  y = "TimeValue",
  ggfunc = ggboxplot,
  add = "jitter",
  color = "TRT",
  palette = "npg",
  ylab="Interferon Level",
  main="Box plots of Interferon Levels by Group"
)


#Histograms; bin width of 5
x1 = ggplot(data=interferon[interferon$TRT=="CONT",])+
  geom_histogram(mapping=aes(x=TimeValue,
                             y=..density..),
                 binwidth=5,
                 color="black",
                 fill="red")+
  geom_density(mapping=aes(x=TimeValue))+
  labs(x="Interferon Levels",
       y="Frequency")+
  ggtitle("Control Group")+
  theme_minimal()

x2 = ggplot(data=interferon[interferon$TRT=="EXPT",])+
  geom_histogram(mapping=aes(x=TimeValue,
                             y=..density..),
                 binwidth=5,
                 color="black",
                 fill="blue")+
  geom_density(mapping=aes(x=TimeValue))+
  labs(x="Interferon Levels",
       y="Frequency")+
  ggtitle("Experimental Group")+
  theme_minimal()

grid.arrange(x1,
             x2,
             nrow=1,
             top="Histograms of Interferon Levels by Group")


#Calculate summary statistics
#Calculate means and standard deviation
intcontmean = mean(interferon[interferon$TRT=="CONT","TimeValue"]$TimeValue)
intexptmean =  mean(interferon[interferon$TRT=="EXPT","TimeValue"]$TimeValue)
intcontmean
intexptmean

intcontsd = sd(interferon[interferon$TRT=="CONT","TimeValue"]$TimeValue)
intexptsd =  sd(interferon[interferon$TRT=="EXPT","TimeValue"]$TimeValue)
intcontsd
intexptsd

intcontsum = summary(interferon[interferon$TRT=="CONT","TimeValue"]$TimeValue)
intcontsum
intexptsum = summary(interferon[interferon$TRT=="EXPT","TimeValue"]$TimeValue)
intexptsum

#Brown-forsythe test for equal variances
bf.test(TimeValue~TRT,
        data=interferon)

#Reject hypothesis that variances are equal

#Shapiro Wilk test for normality
shapiro.test(interferon[interferon$TRT=="CONT","TimeValue"]$TimeValue)

shapiro.test(interferon[interferon$TRT=="EXPT","TimeValue"]$TimeValue)

#Cannot reject null hypotheses that distributions are normal

#Use t-test with unequal variances
t.test(interferon[interferon$TRT=="CONT","TimeValue"]$TimeValue,
       interferon[interferon$TRT=="EXPT","TimeValue"]$TimeValue,
       var.equal=FALSE)

#Reject null hypothesis that means are equal

################################################################################
################################################################################

##25OHD

################################################################################
################################################################################

#Convert data from wide to long format
ohd=ohd %>%
  pivot_longer(cols=TIME1:TIME4,
               names_to = "TimeType",
               values_to = "TimeValue")

#Box plot
ggsummarystats(
  ohd,
  x = "TRT",
  y = "TimeValue",
  ggfunc = ggboxplot,
  add = "jitter",
  color = "TRT",
  palette = "npg",
  ylab="25-OHD Level",
  main="Box plots of 25-OHD Levels by Group"
)


#Histograms; bin width of 5
x1=ggplot(data=ohd[ohd$TRT=="CONT",])+
  geom_histogram(mapping=aes(x=TimeValue,
                             y=..density..),
                 binwidth=5,
                 color="black",
                 fill="red")+
  geom_density(mapping=aes(x=TimeValue))+
  labs(x="Times",
       y="Frequency")+
  ggtitle("Control Group")+
  theme_minimal()

x2=ggplot(data=ohd[ohd$TRT=="EXPT",])+
  geom_histogram(mapping=aes(x=TimeValue,
                             y=..density..),
                 binwidth=5,
                 color="black",
                 fill="blue")+
  geom_density(mapping=aes(x=TimeValue))+
  labs(x="Times",
       y="Frequency")+
  ggtitle("Experimental Group")+
  theme_minimal()

grid.arrange(x1,
             x2,
             nrow=1,
             top="Histograms of 25-OHD Levels by Group")

#Calculate means and standard deviation
ohdcontmean = mean(ohd[ohd$TRT=="CONT","TimeValue"]$TimeValue)
ohdexptmean =  mean(ohd[ohd$TRT=="EXPT","TimeValue"]$TimeValue)
ohdcontmean
ohdexptmean

ohdcontsd = sd(ohd[ohd$TRT=="CONT","TimeValue"]$TimeValue)
ohdexptsd =  sd(ohd[ohd$TRT=="EXPT","TimeValue"]$TimeValue)
ohdcontsd
ohdexptsd

ohdcontsum = summary(ohd[ohd$TRT=="CONT","TimeValue"]$TimeValue)
ohdcontsum
intexptsum = summary(ohd[ohd$TRT=="EXPT","TimeValue"]$TimeValue)
intexptsum

#Brown-forsythe test for equal variances
bf.test(TimeValue~TRT,
        data=ohd)

#Reject hypothesis that variances are equal

#Shapiro Wilk test for normality
shapiro.test(ohd[ohd$TRT=="CONT","TimeValue"]$TimeValue)

shapiro.test(ohd[ohd$TRT=="EXPT","TimeValue"]$TimeValue)

#Cannot reject null hypotheses that distributions are normal

#Use t-test with unequal variances
t.test(ohd[ohd$TRT=="CONT","TimeValue"]$TimeValue,
       ohd[ohd$TRT=="EXPT","TimeValue"]$TimeValue,
       var.equal=FALSE)

#Reject null hypothesis that means are equal

################################################################################
################################################################################

##Lactoferrin

################################################################################
################################################################################

#Convert data from wide to long format
lactoferrin=lactoferrin %>%
  pivot_longer(cols=TIME1:TIME4,
               names_to = "TimeType",
               values_to = "TimeValue")

#Box plot
ggsummarystats(
  lactoferrin,
  x = "TRT",
  y = "TimeValue",
  ggfunc = ggboxplot,
  add = "jitter",
  color = "TRT",
  palette = "npg",
  ylab="Lactoferrin Level",
  main="Box plots of Lactoferrin Levels by Group"
)

#Histograms; bin width of 5
x1=ggplot(data=lactoferrin[lactoferrin$TRT=="CONT",])+
  geom_histogram(mapping=aes(x=TimeValue,
                             y=..density..),
                 binwidth=5,
                 color="black",
                 fill="red")+
  geom_density(mapping=aes(x=TimeValue))+
  labs(x="Times",
       y="Frequency")+
  ggtitle("Control Group")+
  theme_minimal()

x2=ggplot(data=lactoferrin[lactoferrin$TRT=="EXPT",])+
  geom_histogram(mapping=aes(x=TimeValue,
                             y=..density..),
                 binwidth=5,
                 color="black",
                 fill="blue")+
  geom_density(mapping=aes(x=TimeValue))+
  labs(x="Times",
       y="Frequency")+
  ggtitle("Experimental Group")+
  theme_minimal()

grid.arrange(x1,
             x2,
             nrow=1,
             top="Histograms of Lactoferrin Levels by Group")

#Calculate means and standard deviation
laccontmean = mean(lactoferrin[lactoferrin$TRT=="CONT","TimeValue"]$TimeValue)
lacexptmean =  mean(lactoferrin[lactoferrin$TRT=="EXPT","TimeValue"]$TimeValue)
laccontmean
lacexptmean

laccontsd = sd(lactoferrin[lactoferrin$TRT=="CONT","TimeValue"]$TimeValue)
lacexptsd =  sd(lactoferrin[lactoferrin$TRT=="EXPT","TimeValue"]$TimeValue)
laccontsd
lacexptsd

laccontsum = summary(lactoferrin[lactoferrin$TRT=="CONT","TimeValue"]$TimeValue)
laccontsum
lacexptsum = summary(lactoferrin[lactoferrin$TRT=="EXPT","TimeValue"]$TimeValue)
lacexptsum

#Brown-forsythe test for equal variances
bf.test(TimeValue~TRT,
        data=lactoferrin)

#Reject hypothesis that variances are equal

#Shapiro Wilk test for normality
shapiro.test(lactoferrin[lactoferrin$TRT=="CONT","TimeValue"]$TimeValue)

shapiro.test(lactoferrin[lactoferrin$TRT=="EXPT","TimeValue"]$TimeValue)

#Treatment group does not appear to be normal

#Perform permutation test
set.seed(123)
x = perm.test(x=lactoferrin[lactoferrin$TRT=="CONT","TimeValue"]$TimeValue,lactoferrin[lactoferrin$TRT=="EXPT","TimeValue"]$TimeValue) 
perm.test(x=lactoferrin[lactoferrin$TRT=="CONT","TimeValue"]$TimeValue,lactoferrin[lactoferrin$TRT=="EXPT","TimeValue"]$TimeValue,plot=TRUE)

mean(lactoferrin[lactoferrin$TRT=="CONT","TimeValue"]$TimeValue)-mean(lactoferrin[lactoferrin$TRT=="EXPT","TimeValue"]$TimeValue)
#Reject null hypothesis that means/medians are equal

################################################################################
################################################################################

##Immunoglobulyn

################################################################################
################################################################################

#Convert data from wide to long format
immunoglobulyn=immunoglobulyn %>%
  pivot_longer(cols=TIME1:TIME4,
               names_to = "TimeType",
               values_to = "TimeValue")

#Box plot
ggsummarystats(
  immunoglobulyn,
  x = "TRT",
  y = "TimeValue",
  ggfunc = ggboxplot,
  add = "jitter",
  color = "TRT",
  palette = "npg",
  ylab="Immunoglobulyn Level",
  main="Box plots of Immunoglobulyn Levels by Group"
)


#Histograms; bin width of 5
x1 = ggplot(data=immunoglobulyn[immunoglobulyn$TRT=="CONT",])+
  geom_histogram(mapping=aes(x=TimeValue,
                             y=..density..),
                 binwidth=5,
                 color="black",
                 fill="red")+
  geom_density(mapping=aes(x=TimeValue))+
  labs(x="Times",
       y="Frequency")+
  ggtitle("Control Group")+
  theme_minimal()

ggplot(data=immunoglobulyn[immunoglobulyn$TRT=="EXPT",])+
  geom_histogram(mapping=aes(x=TimeValue,
                             y=..density..),
                 binwidth=5,
                 color="black",
                 fill="blue")+
  geom_density(mapping=aes(x=TimeValue))+
  labs(x="Times",
       y="Frequency")+
  ggtitle("Treatment Group")+
  theme_minimal()

grid.arrange(x1,
             x2,
             nrow=1,
             top="Histograms of Immunoglobulyn Levels by Group")

#Calculate means and standard deviation
immcontmean = mean(immunoglobulyn[immunoglobulyn$TRT=="CONT","TimeValue"]$TimeValue)
immexptmean =  mean(immunoglobulyn[immunoglobulyn$TRT=="EXPT","TimeValue"]$TimeValue)
immcontmean
immexptmean

immcontsd = sd(immunoglobulyn[immunoglobulyn$TRT=="CONT","TimeValue"]$TimeValue)
immexptsd =  sd(immunoglobulyn[immunoglobulyn$TRT=="EXPT","TimeValue"]$TimeValue)
immcontsd
immexptsd

immcontsum = summary(immunoglobulyn[immunoglobulyn$TRT=="CONT","TimeValue"]$TimeValue)
immcontsum
immexptsum = summary(immunoglobulyn[immunoglobulyn$TRT=="EXPT","TimeValue"]$TimeValue)
immexptsum

#Brown-forsythe test for equal variances
bf.test(TimeValue~TRT,
        data=immunoglobulyn)

#Reject hypothesis that variances are equal

#Shapiro Wilk test for normality
shapiro.test(immunoglobulyn[immunoglobulyn$TRT=="CONT","TimeValue"]$TimeValue)

shapiro.test(immunoglobulyn[immunoglobulyn$TRT=="EXPT","TimeValue"]$TimeValue)


#Permutation test
set.seed(123)
x = perm.test(x=immunoglobulyn[immunoglobulyn$TRT=="CONT","TimeValue"]$TimeValue,immunoglobulyn[immunoglobulyn$TRT=="EXPT","TimeValue"]$TimeValue) 
x
perm.test(x=immunoglobulyn[immunoglobulyn$TRT=="CONT","TimeValue"]$TimeValue,immunoglobulyn[immunoglobulyn$TRT=="EXPT","TimeValue"]$TimeValue,plot=TRUE)


mean(immunoglobulyn[immunoglobulyn$TRT=="CONT","TimeValue"]$TimeValue)-mean(immunoglobulyn[immunoglobulyn$TRT=="EXPT","TimeValue"]$TimeValue)
#Reject null hypothesis that means/medians are equal
#Reject null hypothesis that means/medians are equal



################################################################################
################################################################################

##Level/Viability data

################################################################################
################################################################################

level = level %>%
  mutate(LEVEL=fct_relevel(LEVEL,
                           "L0",
                           "L5",
                           "L10",
                           "L15",
                           "L20",
                           "L25"))
#Summary stats
##L0
summary(level[level$LEVEL=="L5" & level$VIABILITY=="EXPT","VALUE"]$VALUE)
sd(level[level$LEVEL=="L5" & level$VIABILITY=="EXPT","VALUE"]$VALUE)

summary(level[level$LEVEL=="L0" & level$VIABILITY=="CONT","VALUE"]$VALUE)
sd(level[level$LEVEL=="L0" & level$VIABILITY=="CONT","VALUE"]$VALUE)

##L5
summary(level[level$LEVEL=="L5" & level$VIABILITY=="EXPT","VALUE"]$VALUE)
sd(level[level$LEVEL=="L5" & level$VIABILITY=="EXPT","VALUE"]$VALUE)

summary(level[level$LEVEL=="L5" & level$VIABILITY=="CONT","VALUE"]$VALUE)
sd(level[level$LEVEL=="L5" & level$VIABILITY=="CONT","VALUE"]$VALUE)

##L10
summary(level[level$LEVEL=="L10" & level$VIABILITY=="EXPT","VALUE"]$VALUE)
sd(level[level$LEVEL=="L10" & level$VIABILITY=="EXPT","VALUE"]$VALUE)

summary(level[level$LEVEL=="L10" & level$VIABILITY=="CONT","VALUE"]$VALUE)
sd(level[level$LEVEL=="L10" & level$VIABILITY=="CONT","VALUE"]$VALUE)



##L15
summary(level[level$LEVEL=="L15" & level$VIABILITY=="EXPT","VALUE"]$VALUE)
sd(level[level$LEVEL=="L15" & level$VIABILITY=="EXPT","VALUE"]$VALUE)

summary(level[level$LEVEL=="L15" & level$VIABILITY=="CONT","VALUE"]$VALUE)
sd(level[level$LEVEL=="L15" & level$VIABILITY=="CONT","VALUE"]$VALUE)

##L20
summary(level[level$LEVEL=="L20" & level$VIABILITY=="EXPT","VALUE"]$VALUE)
sd(level[level$LEVEL=="L20" & level$VIABILITY=="EXPT","VALUE"]$VALUE)

summary(level[level$LEVEL=="L20" & level$VIABILITY=="CONT","VALUE"]$VALUE)
sd(level[level$LEVEL=="L20" & level$VIABILITY=="CONT","VALUE"]$VALUE)

##L25
summary(level[level$LEVEL=="L25" & level$VIABILITY=="EXPT","VALUE"]$VALUE)
sd(level[level$LEVEL=="L25" & level$VIABILITY=="EXPT","VALUE"]$VALUE)

summary(level[level$LEVEL=="L25" & level$VIABILITY=="CONT","VALUE"]$VALUE)
sd(level[level$LEVEL=="L25" & level$VIABILITY=="CONT","VALUE"]$VALUE)


#Box plots
ggsummarystats(
  level,
  x = "LEVEL",
  y = "VALUE",
  ggfunc = ggboxplot,
  add = "jitter",
  color = "VIABILITY",
  palette = "npg",
  ylab="Values",
  main="Box plots of Values by Level and Group"
)

#Fit nested anova model
mod = aov(VALUE~LEVEL*VIABILITY,
          data=level)
xtable(summary(mod))
par(mfrow=c(2,2))
plot(mod)
par(mfrow=c(1,1))

#Test normality of residuals
qqnorm(mod$residuals)
qqline(mod$residuals)

shapiro.test(mod$residuals)

#Plot and extract Tukey's HSD simultaneous comparisons
tukey = TukeyHSD(mod,
         conf.level=0.95)
xtable(tukey$`LEVEL:VIABILITY`)

