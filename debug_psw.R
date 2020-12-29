setwd("C:/Users/Shuxi ZENG/Dropbox/Fourth Year/OW_Survival/codebase")

rm(list=ls())

good_overlap=1; sample_size=200;truncate=40

# args=commandArgs(trailingOnly = TRUE)
# if(length(args)==0){
#   print("No arguments supplied.")
# }else{
#   for(i in 1:length(args)){
#     eval(parse(text=args[[i]]))
#   }
# }
# good_overlap=TRUE
# estimand_id=1
# sample_size=150
# install.packages("pseudo_1.4.3.tar.gz", repos = NULL, type="source")
# estimand_grid=c("mean","survprob","restricted_mean")
# estimand=estimand_grid[estimand_id]

library(MASS)
library(mvtnorm)
library(pseudo)
library(flexsurv)
library(survival)
source("IPWC.R")
source("OW.R")

set.seed(2020)


## Lower dimension of covariates, same from Lunceford setting
X4=sample(c(0,1),sample_size,replace=T,prob=c(0.5,0.5))
X3=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.6-0.2*X4[x],0.2*X4[x]+0.4))}))
X=cbind(1,mvrnorm(n=sample_size,mu=c(0,0),Sigma=matrix(c(2,0.5,0.5,2),2,2)),X3,X4)
colnames(X)=c("X0","X1","X2","X3","X4")
# Propensity Score Model
if (good_overlap==1){
  beta=c(0,0,0,0,0)
}else if (good_overlap==2){
  beta=c(-1.2,0.85,0.9,0.45,-0.25)
}else{
  beta=c(0.4,3,2,-3.5,-2)
}

true_e=1/(1+exp(-X%*%beta))
Z=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(1-true_e[x],true_e[x]))}))


# Survival time
gamma=1;alpha=c(2,1.5,-1,1);lambda=0.0001;v=3
hazard=as.vector(gamma*Z+X[,-1]%*%alpha)
random_simu=runif(sample_size)
Survival_time=(-log(random_simu)/(lambda*exp(hazard)))^(1/v)


# Censor time
Censor_time=runif(sample_size,min=0,max=115)
# Censor_time=runif(sample_size,min=2000,max=2115)



# Observed survival
Censored=Survival_time>Censor_time
Y=Survival_time
Y[Censored]=Censor_time[Censored]
DELTA=1-Censored

weight.type="OW"
estimand.type="SPCE"
ps.threshold = NA
evaluate.time = 50

var.method=2;weight.type="IPW";
estimand.type="ASCE";alpha=c(-1,1)
evaluate.time = 50; ps.threshold = 0.05


var.method=2;weight.type="IPW";
estimand.type="ASCE";alpha=c(-1,0,1)
evaluate.time = 50; ps.threshold = 0.05

source("PSW_pseudo.R")
A=PSW.pseudo(Y, Z, DELTA, X,var_method=1,weight.type="OW",
           estimand.type="SPCE",
           evaluate.time = 40)
pseudo_obs=pseudomean(Y,event=DELTA)
OW(y=pseudo_obs, z=Z, W=X)
