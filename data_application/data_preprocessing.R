
# set working directory
rm(list=ls())
library(nnet)
library("survival")
library("dplyr")
library("PSweight")
library("splines")

#library("mstate")
#library("tidyr")
#library("reshape2")
#library("gbm")

dat<- read.csv("dat_lymph.csv",header=T)

#some patients are treated at diagnosis
# rx0<-filter(dat,days2trx==0) #548 patients
# dat[dat$days2trx==0,]$days2trx<-.5


#check columns with missing values
#sapply(dat, function(x) sum(is.na(x)))
colSums(is.na(dat))
dat<-na.omit(dat)
# dat<-dat[dat$RACE!="99",]
# dat<-dat[dat$SPANISH_HISPANIC_ORIGIN!="9",]


# sum(is.na(dat$NO_HSD_QUAR_00)&is.na(dat$MED_INC_QUAR_00)&is.na(dat$UR_CD_03)
#     &is.na(dat$CROWFLY)&is.na(dat$CS_SITESPECIFIC_FACTOR_8))

dat[dat$CS_SITESPECIFIC_FACTOR_8<=6,]$CS_SITESPECIFIC_FACTOR_8<-6
dat[!dat$RACE %in% c(1,2),]$RACE<-3 #combine non-white & non-black race group
dat[dat$SPANISH_HISPANIC_ORIGIN %in% c(1:8),]$SPANISH_HISPANIC_ORIGIN<- 1
dat<-dat %>% mutate (tstage =ifelse(TNM_CLIN_T %in% c("3","3A","3B","4"),2,1))
dat$TNM_CLIN_T=NULL
#trx=1:"RP", trx=2: "EBRT+AD", trx=3: "EBRT+brachy±AD"
names(dat)<-c("age","race","spanish","insurance","income","education",
              "deyo","dyear","psa","gs", "surgsite", "regdose", "boostdose","surgseq","hormone",
              "fupmonth","id","trx","days2trx","death","totdose","tstage")

dat<- dat %>% mutate(dyearcat=ifelse(dyear %in% c(2004:2007), 1, ifelse(dyear %in% c(2008:2010),2,3)))
dat$dyearcat2<- as.factor(dat$dyear)

which(sapply(dat,is.factor))
numindex<-which(sapply(dat,is.numeric)) #age and psa are continous
numindex

ind2<-numindex[c("race","spanish","insurance","income","education",
                 "deyo","gs","trx","tstage","dyearcat")]

facVars = lapply(dat[ , ind2],
                 function(x) {
                   x = as.factor(x)
                   x
                 })

dat1<-cbind(facVars, dat[ , - ind2])
which(sapply(dat1,is.factor))
dat<-dat1

#-----------------------#
#     11/13/2017        #
#  psa rescale by 1/10  #
#-----------------------#
# dat$psa <- dat$psa/10

#==========================================#
# 10/11/2020 for overlap weights project #
#==========================================#
# use average days per month (365/12 = 30.42) to obtain time-to-treatment 
dat<- subset(dat, gs %in% c(6,7,8,9,10), drop=TRUE)
dat$gs <- factor(dat$gs, levels = c(6,7,8,9,10))
dat$month2trx <- dat$days2trx/30.42
dat$time <- dat$fupmonth - dat$month2trx

# In the dataset, the covariates are specified in the PS model as follows
# the outcome is the vector (time, death), where death is the censoring indicator (1 = death/event, 0 = censoring)

# the follwing code is commented out, but is used to generate the analysis dataset

# write.table(dat, file = "dat_survivalOW.csv", row.names=FALSE, na="",col.names=FALSE, sep=",")

# use the largest category as reference
dat$race <- relevel(factor(dat$race),ref='1')
dat$spanish <- relevel(factor(dat$spanish),ref='0')
dat$insurance <- relevel(factor(dat$insurance),ref='3')
dat$income <- relevel(factor(dat$income),ref='4')
dat$education <- relevel(factor(dat$education),ref='4')
dat$deyo <-relevel(factor(dat$deyo),ref='0')
dat$gs <- relevel(factor(dat$gs),ref='8')
dat$tstage <- relevel(factor(dat$tstage),ref='1')
dat$dyearcat <- relevel(factor(dat$dyearcat),ref='3')


### Covariates
gps.model <- multinom(trx ~ ns(dat$age,df=4) + ns(psa,df=4) + as.factor(race) + as.factor(spanish) +
                        as.factor(insurance) + as.factor(income) + as.factor(education) + as.factor(deyo) +
                        as.factor(gs) + as.factor(tstage) + as.factor(dyearcat), data = dat,
                      maxit = 500,
                      Hess = TRUE,
                      trace = FALSE)
# Extract Design Matrix
X = model.matrix(gps.model)
# Survival Time
Y = dat$time
DELTA = dat$death
# adjust for treatment indicator to be 0:(J-1)
Z = as.numeric(dat$trx)-1

# ## Some naive analysis
# # Estimate generalized propensity scores with natural cubic splines on cont vars
# gps.model <- trx ~ ns(age,df=4) + ns(psa,df=4) + as.factor(race) + as.factor(spanish) +
#   as.factor(insurance) + as.factor(income) + as.factor(education) + as.factor(deyo) +
#   as.factor(gs) + as.factor(tstage) + as.factor(dyearcat)


gps.model <- trx ~ age + psa + as.factor(race) + as.factor(spanish) +
  as.factor(insurance) + as.factor(income) + as.factor(education) + as.factor(deyo) +
  as.factor(gs) + as.factor(tstage) + as.factor(dyearcat)
# 
bal.mult<-SumStat(ps.formula=gps.model, weight=c("IPW","overlap"), data=dat)

pdf("GPS_plot.pdf",width=12,height=4)
par(mfrow=c(1,3))
plot(density(bal.mult$propensity[Z==0,1]),col="red",xlab = "Estimated GPS",bty="n",
     main="Propensity for RP", lwd=2,yaxt='n')
lines(density(bal.mult$propensity[Z==1,1]),col="blue", lwd=2)
lines(density(bal.mult$propensity[Z==2,1]),col="green",lwd=2)

plot(density(bal.mult$propensity[Z==0,2]),col="red",xlab = "Estimated GPS",bty="n",
     main="Propensity for EBRT+AD", lwd=2,yaxt='n')
lines(density(bal.mult$propensity[Z==1,2]),col="blue", lwd=2)
lines(density(bal.mult$propensity[Z==2,2]),col="green",lwd=2)

plot(density(bal.mult$propensity[Z==0,3]),col="red",xlab = "Estimated GPS",bty="n",
     main="Propensity for EBRT+brachy±AD", lwd=2, yaxt='n',ylim=c(0,20))
lines(density(bal.mult$propensity[Z==1,3]),col="blue", lwd=2)
lines(density(bal.mult$propensity[Z==2,3]),col="green",lwd=2)
legend("topright",legend=c("RP","EBRT+AD","EBRT+brachy±AD"),col=c("red","blue","green"),lwd=2)
dev.off()

# pdf("binary_ps.pdf",width=6,height=5)
# plot(density(bal.mult$propensity[Z==0,2]),col="red",xlab = "Estimated PS",bty="n",
#      main="Propensity for receving EBRT+AD", lwd=2,yaxt='n')
# lines(density(bal.mult$propensity[Z==1,2]),col="blue", lwd=2)
# legend("topright",legend=c("RP","EBRT+AD"),col=c("red","blue"),lwd=2)
# dev.off()

# output_table = NULL
output_table = cbind(bal.mult$unweighted.sumstat$mres%*%(table(Z)/sum(table(Z))),
  bal.mult$unweighted.sumstat$mres,
apply(bal.mult$unweighted.sumstat$ASD.unweighted.var,1,FUN = function(x){max(abs(x))}),
apply(bal.mult$IPW.sumstat$ASD.unweighted.var,1,FUN = function(x){max(abs(x))}),
apply(bal.mult$overlap.sumstat$ASD.unweighted.var,1,FUN = function(x){max(abs(x))}))

output_table = cbind(bal.mult$unweighted.sumstat$mres%*%(table(Z)),
                     bal.mult$unweighted.sumstat$mres[,1]*table(Z)[1],
                     bal.mult$unweighted.sumstat$mres[,2]*table(Z)[2],
                     bal.mult$unweighted.sumstat$mres[,3]*table(Z)[3],
                     apply(bal.mult$unweighted.sumstat$ASD.unweighted.var,1,FUN = function(x){max(abs(x))}),
                     apply(bal.mult$IPW.sumstat$ASD.unweighted.var,1,FUN = function(x){max(abs(x))}),
                     apply(bal.mult$overlap.sumstat$ASD.unweighted.var,1,FUN = function(x){max(abs(x))}))

library(xtable)
xtable(output_table,digits = c(0,0,0,0,0,3,3,3))
# 
# # check GPS distribution
pdf("GPS_distributions.pdf", width = 7, height = 7, paper = "special")
par(mar=c(5.1, 6.1, 4.1, 2.1))
par(mfrow=c(3,1))
plot(bal.mult, type = "density")

dev.off()
# 
# # check balance: it looks like IPW does a good job in balancing covariates based on the 0.1 threshold
# #                but OW should provide more efficient results
# pdf("Cov_balance.pdf", width = 13, height = 7, paper = "special")
# 
# plot(bal.mult, metric = "ASD")
# plot(bal.mult, metric = "PSD")
# 
# dev.off()
# 
# #possible subsetting
# dats<- filter(dat, (trx==3 & hormone==1) | (trx==2 & totdose>=7920 & regdose != 99999 & boostdose !=99999) )
# #
# # somewhat naive analysis; i.e. fitting weighted Cox models
# library(survey)
# 
# # design objects
# des.IPW <- svydesign(ids=~1, weights=~bal.mult$ps.weights$IPW, data=dat)
# des.Overlap <- svydesign(ids=~1, weights=~bal.mult$ps.weights$overlap, data=dat)
# 
# # unweighted analysis - subject to large bias
# cox <- coxph(Surv(time, death) ~ factor(trx), data = dat)
# summary(cox)
# 
# # weighted analysis: OW will lead to significant results for both comparisons,
# #                    while IPW only leads to significant results for one comparison
# 
# cox_Overlap <- svycoxph(Surv(time, death)~factor(trx),design=des.Overlap)
# cox_IPW <- svycoxph(Surv(time, death)~factor(trx),design=des.IPW)
# 
# summary(cox_Overlap)
# summary(cox_IPW)
# 
# # All pairwise CIs exclude zero with OW on the HR scale!
# # Estimand on Hazard Ratio
# x <- summary(cox)$coef
# out <- matrix(NA,3,3)
# rownames(out) <- c("2 vs 1", "3 vs 1", "2 vs 3")
# colnames(out) <- c("HR", "LCL", "UCL")
# out[1,] <- exp(c(x[1,1], x[1,1]-qnorm(0.975)*x[1,3], x[1,1]+qnorm(0.975)*x[1,3]))
# out[2,] <- exp(c(x[2,1], x[2,1]-qnorm(0.975)*x[2,3], x[2,1]+qnorm(0.975)*x[2,3]))
# out[3,] <- exp(c(x[1,1]-x[2,1],
#                  (x[1,1]-x[2,1])-qnorm(0.975)*sqrt(t(c(1,-1))%*%vcov(cox_Overlap)%*%c(1,-1)),
#                  (x[1,1]-x[2,1])+qnorm(0.975)*sqrt(t(c(1,-1))%*%vcov(cox_Overlap)%*%c(1,-1))))
# round(out,2)
# 
# x <- summary(cox_Overlap)$coef
# out <- matrix(NA,3,3)
# rownames(out) <- c("2 vs 1", "3 vs 1", "2 vs 3")
# colnames(out) <- c("HR", "LCL", "UCL")
# out[1,] <- exp(c(x[1,1], x[1,1]-qnorm(0.975)*x[1,3], x[1,1]+qnorm(0.975)*x[1,3]))
# out[2,] <- exp(c(x[2,1], x[2,1]-qnorm(0.975)*x[2,3], x[2,1]+qnorm(0.975)*x[2,3]))
# out[3,] <- exp(c(x[1,1]-x[2,1],
#                  (x[1,1]-x[2,1])-qnorm(0.975)*sqrt(t(c(1,-1))%*%vcov(cox_Overlap)%*%c(1,-1)),
#                  (x[1,1]-x[2,1])+qnorm(0.975)*sqrt(t(c(1,-1))%*%vcov(cox_Overlap)%*%c(1,-1))))
# round(out,2)
# 
# x <- summary(cox_IPW)$coef
# out <- matrix(NA,3,3)
# rownames(out) <- c("2 vs 1", "3 vs 1", "2 vs 3")
# colnames(out) <- c("HR", "LCL", "UCL")
# out[1,] <- exp(c(x[1,1], x[1,1]-qnorm(0.975)*x[1,3], x[1,1]+qnorm(0.975)*x[1,3]))
# out[2,] <- exp(c(x[2,1], x[2,1]-qnorm(0.975)*x[2,3], x[2,1]+qnorm(0.975)*x[2,3]))
# out[3,] <- exp(c(x[1,1]-x[2,1],
#                  (x[1,1]-x[2,1])-qnorm(0.975)*sqrt(t(c(1,-1))%*%vcov(cox_IPW)%*%c(1,-1)),
#                  (x[1,1]-x[2,1])+qnorm(0.975)*sqrt(t(c(1,-1))%*%vcov(cox_IPW)%*%c(1,-1))))
# round(out,2)

