## This code performs the statistical analyses described in the paper : Synchronization between keyboard typing
## and neural oscillations
## It performs linear mixed models on RT, peak frequency extracted from Kernel density estimation and 
## clustering value at each investigated frequency

## Version 2 19/11/2020
## duprez.joan@gmail.com
## https://duprezjoan.wixsite.com/joanduprez

# Needs the nlme, MuMIn, and multcomp and forecast packages


lapply(c("nlme", "MuMIn", "multcomp", "forecast"), require, character.only = TRUE)
inpath = 'your data folder'
outpath = 'your result folder'
## ---
## RT ANALYSES
## ---

dataRT = read.csv(paste(inpath, "/dataRT4stat.csv", sep = ""))
dataRT$Subject = as.factor(dataRT$Subject)
dataRT$Precision = as.factor(dataRT$Precision)
dataRT$Condition = as.factor(dataRT$Condition)


modRT = lme(BoxCox(RT, lambda = "auto")~Condition+Precision, random=~1|Subject, data=dataRT)
anova(modRT)


residus<-residuals(modRT)
qqnorm(residus)
qqline(residus)

plot(fitted(modRT), residuals(modRT),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modRT), residuals(modRT)))

# conditional R squared instead of effect size
r.squaredGLMM(modRT)
# Posthoc tests
# for experimental conditions
require(multcomp)
postRT = glht(modRT,linfct=mcp(Condition="Tukey"))
summary(postRT)
# Get condition specific mean and sd to interpret
with(dataRT, aggregate(RT, list(dataRT2$Condition), mean))
with(dataRT, aggregate(RT, list(dataRT2$Condition), sd))

## ---
## ACCURACY ANALYSES
##

dataacc = read.table(paste(inpath, "/err_det.txt"), sep = "", header=TRUE)
dataacc$n = as.factor(dataacc$n)
dataacc$condition = as.factor(dataacc$condition)

modacc = lme(log(acc)~condition, random=~1|n, data=dataacc)
anova(modacc)

residus<-residuals(modacc)
qqnorm(residus)
qqline(residus)

plot(fitted(modacc), residuals(modacc),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modacc), residuals(modacc)))

r.squaredGLMM(modacc)

postacc = glht(modacc,linfct=mcp(condition="Tukey"))
summary(postacc)

with(dataacc, aggregate(acc, list(dataacc$condition), mean))
with(dataacc, aggregate(acc, list(dataacc$condition), sd))

## ---
## IKI ANALYSES
## ---

dataIKI = read.csv(paste(inpath, "/dataIKI4stat.csv", sep = ""))

dataIKI$sub = as.factor(dataIKI2$sub) # Otherwise sub will be considered as a continuous quantitative variable
dataIKI$Condition = as.factor(dataIKI2$Condition) # Otherwise sub will be considered as a continuous quantitative variable
dataIKI$outcome = as.factor(dataIKI2$outcome) # Otherwise sub will be considered as a continuous quantitative variable


modIKI = lme(log(IKI)~Condition+outcome, random=~1|sub, data=dataIKI)
anova(modIKI)

# graphical check of assumptions
residus<-residuals(modIKI)
qqnorm(residus)
qqline(residus)

plot(fitted(modIKI), residuals(modIKI),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modIKI), residuals(modIKI)))

# conditional R squared instead of effect size
r.squaredGLMM(modIKI)
# Posthoc tests
# for experimental conditions
require(multcomp)
postIKI = glht(modIKI,linfct=mcp(Condition="Tukey"))
summary(postIKI)
# Get condition specific mean and sd to interpret
with(dataIKI, aggregate(IKI, list(dataIKI$Condition), mean))
with(dataIKI, aggregate(IKI, list(dataIKI$Condition), sd))

# for outcome
postIKI2 = glht(modIKI,linfct=mcp(outcome="Tukey"))
summary(postIKI2)
# Get condition specific mean and sd to interpret
with(dataIKI, aggregate(IKI, list(dataIKI$outcome), mean))
with(dataIKI, aggregate(IKI, list(dataIKI$outcome), sd))

## Same analyses but separate linguistic vs non-linguistic stimuli
dataIKI = read.csv(paste(inpath, "/dataIKI4stat.csv", sep = ""))

dataIKI$ling = dataIKI$Condition

dataIKI$ling [which(dataIKI$ling =="Words")]="L"
dataIKI$ling [which(dataIKI$ling =="Sentences")]="L"
dataIKI$ling [which(dataIKI$ling =="Pseudo-words")]="NL"
dataIKI$ling [which(dataIKI$ling =="Pseudo-sentences")]="NL"

dataIKI$sub = as.factor(dataIKI$sub)
dataIKI$ling = as.factor(dataIKI$ling) 
dataIKI$outcome = as.factor(dataIKI$outcome) 

modIKI = lme(log(IKI)~ling+outcome, random=~1|sub, data=dataIKI)
anova(modIKI)

# graphical check of assumptions
residus<-residuals(modIKI)
qqnorm(residus)
qqline(residus)

plot(fitted(modIKI), residuals(modIKI),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modIKI), residuals(modIKI)))

# conditional R squared instead of effect size
r.squaredGLMM(modIKI)
# Posthoc tests
# for experimental conditions
postIKI = glht(modIKI,linfct=mcp(outcome="Tukey"))
summary(postIKI)
# Get condition specific mean and sd to interpret
with(dataIKI, aggregate(IKI, list(dataIKI$ling), mean))
with(dataIKI, aggregate(IKI, list(dataIKI$ling), sd))

# Comp IKI of corrected errors with and withou BS IKI

dataIKI = read.table(paste(inpath, "/comp_correctedERRIKI_BS_wBS.txt"), sep = "")
colnames(dataIKI)=c("n","IKI", "cond","backspace")
## stats
dataIKI$n = as.factor(dataIKI$n) # Otherwise sub will be considered as a continuous quantitative variable
dataIKI$cond = as.factor(dataIKI$cond) # Otherwise sub will be considered as a continuous quantitative variable
dataIKI$backspace = as.factor(dataIKI$backspace) # Otherwise sub will be considered as a continuous quantitative variable


modIKI = lme(log(IKI)~cond+backspace, random=~1|n, data=dataIKI)
anova(modIKI)

# graphical check of assumptions
residus<-residuals(modIKI)
qqnorm(residus)
qqline(residus)

plot(fitted(modIKI), residuals(modIKI),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modIKI), residuals(modIKI)))

# conditional R squared instead of effect size
r.squaredGLMM(modIKI)
# Posthoc tests
# for experimental conditions
require(multcomp)
postIKI = glht(modIKI,linfct=mcp(Condition="Tukey"))
summary(postIKI)
# Get condition specific mean and sd to interpret
with(dataIKI, aggregate(IKI, list(dataIKI$backspace), mean))
with(dataIKI, aggregate(IKI, list(dataIKI$backspace), sd))

# for outcome
postIKI2 = glht(modIKI,linfct=mcp(outcome="Tukey"))
summary(postIKI2)
# Get outcome specific mean to interpret
with(dataIKI, aggregate(IKI, list(dataIKI$outcome), mean))
with(dataIKI, aggregate(IKI, list(dataIKI$outcome), sd))


## ---
## KDE ANALYSES
## ---

dataKDE = read.csv(paste(inpath, "/dataKDE4stat.csv", sep = ""))
dataKDE$sub = as.factor(dataKDE$sub)
dataKDE$Condition = as.factor(dataKDE$Condition)
dataKDE$outcome = as.factor(dataKDE$outcome)

modKDE = lme(BoxCox(peak.frequency, lambda = "auto")~Condition+outcome, random=~1|sub, data=dataKDE)
anova(modKDE)
# graphical check of assumptions
residus<-residuals(modKDE)
qqnorm(residus)
qqline(residus)

plot(fitted(modKDE), residuals(modKDE),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modKDE), residuals(modKDE)))

# conditional R squared instead of effect size
r.squaredGLMM(modKDE)

# Posthoc tests
# for experimental conditions
require(multcomp)
postKDE = glht(modKDE,linfct=mcp(Condition="Tukey"))
summary(postKDE)
# Get condition specific mean and to interpret
with(dataKDE, aggregate(peak.frequency, list(dataKDE$Condition), mean))
with(dataKDE, aggregate(peak.frequency, list(dataKDE$Condition), sd))

# for outcome
postKDE2 = glht(modKDE,linfct=mcp(outcome="Tukey"))
summary(postKDE2)
# Get outcome specific mean and sd to interpret
with(dataKDE, aggregate(peak.frequency, list(dataKDE$outcome), mean))
with(dataKDE, aggregate(peak.frequency, list(dataKDE$outcome), sd))

# Same analyses but separate linguistic versus non-linguistic stimuli
dataKDE = read.csv(paste(inpath, "/dataKDE4stat.csv", sep = ""))

dataKDE$ling=dataKDE$Condition

dataKDE$ling [which(dataKDE$ling =="W")]="L"
dataKDE$ling [which(dataKDE$ling =="S")]="L"
dataKDE$ling [which(dataKDE$ling =="pW")]="NL"
dataKDE$ling [which(dataKDE$ling =="pS")]="NL"

dataKDE$sub = as.factor(dataKDE$sub)
dataKDE$ling = as.factor(dataKDE$ling)
dataKDE$outcome = as.factor(dataKDE$outcome)

modKDE = lme(BoxCox(peak.frequency, lambda = "auto")~ling+outcome, random=~1|sub, data=dataKDE)
anova(modKDE)
# graphical check of assumptions
residus<-residuals(modKDE)
qqnorm(residus)
qqline(residus)

plot(fitted(modKDE), residuals(modKDE),
     xlab = "Fitted Values", ylab = "Residuals")
abline(h=0, lty=2)
lines(smooth.spline(fitted(modKDE), residuals(modKDE)))

# conditional R squared instead of effect size
r.squaredGLMM(modKDE)

# Posthoc tests
# for experimental conditions
require(multcomp)
postKDE = glht(modKDE,linfct=mcp(outcome="Tukey"))
summary(postKDE)
# Get condition specific mean and sd to interpret
with(dataKDE, aggregate(peak.frequency, list(dataKDE$ling), mean))
with(dataKDE, aggregate(peak.frequency, list(dataKDE$ling), sd))


## ---
## PHASE CLUSTERING
## ---

dataCLUST = read.csv(paste(inpath, "/dataCLUST4stat_perm.csv", sep = ""))

dataCLUST$Subject = as.factor(dataCLUST$Subject)
dataCLUST$Frequency = as.factor(dataCLUST$Frequency)

# linear mixed model at each freq
# output vector
resCLUST <- vector()
resCLUST
require(nlme)
require(MuMIn)
resCLUST <- matrix(, nrow = 28, ncol = 10)
for(freqi in 2:28){
  
  tempdat = droplevels(dataCLUST[which(dataCLUST$Frequency == freqi),])
  modCLUST = anova(lme(ITPCz~Precision+Condition , random = ~1|Subject, data=tempdat))
  resCLUST[freqi, 1] = modCLUST[2,1]
  resCLUST[freqi, 2] = modCLUST[2,2]
  resCLUST[freqi, 3] = modCLUST[2,3]
  resCLUST[freqi, 4] = modCLUST[2,4]
  resCLUST[freqi, 5] = modCLUST[3,1]
  resCLUST[freqi, 6] = modCLUST[3,2]
  resCLUST[freqi, 7] = modCLUST[3,3]
  resCLUST[freqi, 8] = modCLUST[3,4]
  resCLUST[freqi, 9] = r.squaredGLMM(lme(ITPCz~Precision+Condition , random = ~1|Subject, data=tempdat))[1]
  resCLUST[freqi, 10] = r.squaredGLMM(lme(ITPCz~Precision+Condition , random = ~1|Subject, data=tempdat))[2]
  
}

write.csv(resCLUST, paste(outpath,"/results_clust_anova.csv", sep = ""))
write.csv(resCLUST, paste(outpath,"/CLUSTsignif_inter.csv", sep = ""))

# random graphical check of assumptions

# take a random freqi from 1 to 25
randfreq = sample(1:25, 1)
tempdat = droplevels(dataCLUST[which(dataCLUST$Frequency == randfreq-1),])
modCLUST = lme(ITPCz~Precision*Condition, random=~1|Subject, data=tempdat) # cube root is the best transform compared to log, 1/n and squared root
# graphical check of assumptions
residus<-residuals(modCLUST)
qqnorm(residus, main=randfreq)
qqline(residus)

plot(fitted(modCLUST), residuals(modCLUST),
     xlab = "Fitted Values", ylab = "Residuals", main = randfreq)
abline(h=0, lty=2)
lines(smooth.spline(fitted(modCLUST), residuals(modCLUST)))

# Same analyses but separating linguistic vs. non-linguistic stimuli and word vs. sentences
dataCLUST = read.csv(paste(inpath, "/dataCLUST4stat_perm.csv", sep = ""))

dataCLUST$ling = dataCLUST$Condition
dataCLUST$WS = dataCLUST$Condition

dataCLUST$ling [which(dataCLUST$ling =="Words")]="L"
dataCLUST$ling [which(dataCLUST$ling =="Sentences")]="L"
dataCLUST$ling [which(dataCLUST$ling =="Pseudowords")]="NL"
dataCLUST$ling [which(dataCLUST$ling =="Pseudosentences")]="NL"

dataCLUST$WS [which(dataCLUST$WS =="Words")]="W"
dataCLUST$WS [which(dataCLUST$WS =="Sentences")]="S"
dataCLUST$WS [which(dataCLUST$WS =="Pseudowords")]="W"
dataCLUST$WS [which(dataCLUST$WS =="Pseudosentences")]="S"

dataCLUST$Subject = as.factor(dataCLUST$Subject)
dataCLUST$Frequency = as.factor(dataCLUST$Frequency)
dataCLUST$ling = as.factor(dataCLUST$ling)
dataCLUST$WS = as.factor(dataCLUST$WS)
dataCLUST$Precision= as.factor(dataCLUST$Precision)

resCLUST <- matrix(, nrow = 28, ncol = 14)
for(freqi in 2:28){
  
  tempdat = droplevels(dataCLUST[which(dataCLUST$Frequency == freqi),])
  modCLUST = anova(lme(ITPCz~Precision+ling+WS , random = ~1|Subject, data=tempdat))
  resCLUST[freqi, 1] = modCLUST[2,1]
  resCLUST[freqi, 2] = modCLUST[2,2]
  resCLUST[freqi, 3] = modCLUST[2,3]
  resCLUST[freqi, 4] = modCLUST[2,4]
  resCLUST[freqi, 5] = modCLUST[3,1]
  resCLUST[freqi, 6] = modCLUST[3,2]
  resCLUST[freqi, 7] = modCLUST[3,3]
  resCLUST[freqi, 8] = modCLUST[3,4]
  resCLUST[freqi, 9] = modCLUST[4,1]
  resCLUST[freqi, 10] = modCLUST[4,2]
  resCLUST[freqi, 11] = modCLUST[4,3]
  resCLUST[freqi, 12] = modCLUST[4,4]
  resCLUST[freqi, 13] = r.squaredGLMM(lme(ITPCz~Precision+ling+WS , random = ~1|Subject, data=tempdat))[1]
  resCLUST[freqi, 14] = r.squaredGLMM(lme(ITPCz~Precision+ling+WS , random = ~1|Subject, data=tempdat))[2]
  
}

colnames(resCLUST)=c("ndf_Precision","ddf_Precision","F_Precision","p_precision","ndf_ling","ddf_ling","F_ling","p_ling", "ndf_WS","ddf_WS","F_WS","p_WS", "R2m", "R2c")
write.csv(resCLUST, paste(outpath,"/results_clust_anova_ling_words.csv", sep = ""))
write.csv(resCLUST, paste(outpath,"/CLUSTsignif_inter.csv", sep = ""))


lingITPCz = with(dataCLUST, aggregate(ITPCz, list(dataCLUST$ling, dataCLUST$Frequency), mean))
write.csv(lingITPCz, paste(bpath,"Dropbox/SINS/Keyboard/Keysync/lingITPCz.csv", sep = ""))

WSITPCz = with(dataCLUST, aggregate(ITPCz, list(dataCLUST$WS, dataCLUST$Frequency), mean))
write.csv(WSITPCz, paste(bpath,"Dropbox/SINS/Keyboard/Keysync/WSITPCz.csv", sep = ""))



# END