setwd("C:/Users/Shuxi ZENG/Dropbox/Fourth Year/OW_Survival/codebase")

rm(list=ls())

good_overlap=2; sample_size=200;truncate=50;
n_simu=200;n_mc=10000

multi.arm = F

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
# library(pseudo)
library(flexsurv)
library(survival)
source("IPWC.R")
source("OW.R")
source("PSW_pseudo.R")
# source("Mao_Method_func.R")
set.seed(2020)
coverage_rate_calc<-function(point_est,se_est,true)
{
  non.na.id = which(!is.na(point_est))
  point_est = point_est[non.na.id]
  se_est = se_est[non.na.id]
  if(length(true)>2){
    true = true[non.na.id]
  }

  return (mean((true<point_est+se_est*qnorm(0.975))&(true>=point_est-se_est*qnorm(0.975))))
}





ow_mse=NULL
ipw_mse=NULL
ow_bias=NULL
ipw_bias=NULL
ow_coverage=NULL
ipw_coverage=NULL
cox_q_bias=NULL
cox_q_mse=NULL

ow_est=matrix(NA,n_simu,3)
ipw_est=matrix(NA,n_simu,3)
unadj_est=matrix(NA,n_simu,3)
cox_q_est=matrix(NA,n_simu,3)

ow_est_mao=matrix(NA,n_simu,3)
mw_est_mao=matrix(NA,n_simu,3)
ipw_est_mao=matrix(NA,n_simu,3)
uw_est_mao=matrix(NA,n_simu,3)

ow_se=matrix(NA,n_simu,3)
ipw_se=matrix(NA,n_simu,3)

ow_se_mao=matrix(NA,n_simu,3)
mw_se_mao=matrix(NA,n_simu,3)
ipw_se_mao=matrix(NA,n_simu,3)
uw_se_mao=matrix(NA,n_simu,3)

true_est_finite=matrix(NA,n_simu,3)
true_est_ow_finite=matrix(NA,n_simu,3)
true_est=matrix(NA,n_mc,3)
true_est_ow=matrix(NA,n_mc,3)

ow_est_by_group=matrix(NA,n_simu,3)
ipw_est_by_group=matrix(NA,n_simu,3)
ow_se_by_group=matrix(NA,n_simu,3)
ipw_se_by_group=matrix(NA,n_simu,3)

#########
######### Use MC to approximate true
for (i in 1:n_mc){
  ## Lower dimension of covariates, same from Lunceford setting
  X4=sample(c(0,1),sample_size,replace=T,prob=c(0.5,0.5))
  X3=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.6-0.2*X4[x],0.2*X4[x]+0.4))}))
  X=cbind(1,mvrnorm(n=sample_size,mu=c(0,0),Sigma=matrix(c(2,0.5,0.5,2),2,2)),X3,X4)
  colnames(X)=c("X0","X1","X2","X3","X4")
  # Propensity Score Model
  beta.0 = c(0,0,0,0,0)
  if (good_overlap==1){
    beta.1=c(0,0,0,0,0)
    beta.2=c(0,0,0,0,0)
  }else if (good_overlap==2){
    beta.1 = c(-0.4,0.85,0.9,0.45,-0.25)
    beta.2 = beta.1*0.2
  }else{
    beta.1 = c(1.2,1.5,1,-1.5,-1)
    beta.2 = beta.1*0.2
  }
  
  if(multi.arm){
    beta = cbind(beta.0,beta.1,beta.2)
  }else{
    beta = cbind(beta.0,beta.1)
  }
  exp.part = exp(X%*%beta)
  soft.max = exp.part/rowSums(exp.part)
  
  # true_e=1/(1+exp(-X%*%beta))
  Z=unlist(lapply(1:sample_size,FUN=function(x){sample(0:(ncol(soft.max)-1),1,prob=soft.max[x,])}))
  
  
  # Survival time
  gamma.1=0.2;gamma.2=0.5;alpha=c(2,1.5,-1,1);lambda=0.0001;v=3
  hazard=as.vector(gamma.1*as.numeric(Z==1)+gamma.2*as.numeric(Z==2)+X[,-1]%*%alpha)
  random_simu=runif(sample_size)
  Survival_time=(-log(random_simu)/(lambda*exp(hazard)))^(1/v)
  
  # get true est
  hazard_0=as.vector(X[,-1]%*%alpha)
  hazard_1=as.vector(gamma.1+X[,-1]%*%alpha)
  hazard_2=as.vector(gamma.2+X[,-1]%*%alpha)
  Survival_time_0=(-log(random_simu)/(lambda*exp(hazard_0)))^(1/v)
  Survival_time_1=(-log(random_simu)/(lambda*exp(hazard_1)))^(1/v)
  Survival_time_2=(-log(random_simu)/(lambda*exp(hazard_2)))^(1/v)
  
  # truncate=min(quantile(Survival_time_1[Z==1],0.75),quantile(Survival_time_0[Z==0],0.75))
  tilt.h = 1/rowSums(1/soft.max)
  if(!multi.arm){
    true_est[i,1]=mean(Survival_time_1)-mean(Survival_time_0)
    true_est_ow[i,1]=sum((Survival_time_1-Survival_time_0)*tilt.h)/sum(tilt.h)
    
    true_est[i,3]=mean(Survival_time_1>truncate)-mean(Survival_time_0>truncate)
    true_est_ow[i,3]=sum((as.numeric(Survival_time_1>truncate)-as.numeric(Survival_time_0>truncate))*tilt.h)/sum(tilt.h)

    Survival_time_1[Survival_time_1>truncate]=truncate
    Survival_time_0[Survival_time_0>truncate]=truncate
    true_est[i,2]=mean(Survival_time_1)-mean(Survival_time_0)
    true_est_ow[i,2]=sum((Survival_time_1-Survival_time_0)*tilt.h)/sum(tilt.h)
  
  }else{
    true_est[i,1] = mean(Survival_time_2)-mean(Survival_time_0)
    true_est_ow[i,1] = sum((Survival_time_2-Survival_time_0)*tilt.h)/sum(tilt.h)
    
    true_est[i,3]= mean(Survival_time_2>truncate)-mean(Survival_time_0>truncate)
    true_est_ow[i,3] = sum((as.numeric(Survival_time_2>truncate)-as.numeric(Survival_time_0>truncate))*tilt.h)/sum(tilt.h)

    Survival_time_1[Survival_time_2>truncate] = truncate
    Survival_time_0[Survival_time_0>truncate] = truncate
    true_est[i,2] = mean(Survival_time_2)-mean(Survival_time_0)
    true_est_ow[i,2] = sum((Survival_time_2-Survival_time_0)*tilt.h)/sum(tilt.h)
  }



  if (i%%5000==0)
  {
    print(paste("== MC",i,"=="))
  }
}
#########

for (i in (1:n_simu)){
  ## Lower dimension of covariates, same from Lunceford setting
  X4=sample(c(0,1),sample_size,replace=T,prob=c(0.5,0.5))
  X3=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.6-0.2*X4[x],0.2*X4[x]+0.4))}))
  X=cbind(1,mvrnorm(n=sample_size,mu=c(0,0),Sigma=matrix(c(2,0.5,0.5,2),2,2)),X3,X4)
  colnames(X)=c("X0","X1","X2","X3","X4")
  # Propensity Score Model
  beta.0 = c(0,0,0,0,0)
  if (good_overlap==1){
    beta.1=c(0,0,0,0,0)
    beta.2=c(0,0,0,0,0)
  }else if (good_overlap==2){
    beta.1 = c(-0.4,0.85,0.9,0.45,-0.25)
    beta.2 = beta.1*0.2
  }else{
    beta.1 = c(1.2,1.5,1,-1.5,-1)
    beta.2 = beta.1*0.2
  }

  if(multi.arm){
    beta = cbind(beta.0,beta.1,beta.2)
  }else{
    beta = cbind(beta.0,beta.1)
  }
  exp.part = exp(X%*%beta)
  soft.max = exp.part/rowSums(exp.part)

  # true_e=1/(1+exp(-X%*%beta))
  Z=unlist(lapply(1:sample_size,FUN=function(x){sample(0:(ncol(soft.max)-1),1,prob=soft.max[x,])}))

  
  # Survival time
  gamma.1=0.2;gamma.2=0.5;alpha=c(2,1.5,-1,1);lambda=0.0001;v=3
  hazard=as.vector(gamma.1*as.numeric(Z==1)+gamma.2*as.numeric(Z==2)+X[,-1]%*%alpha)
  random_simu=runif(sample_size)
  Survival_time=(-log(random_simu)/(lambda*exp(hazard)))^(1/v)
  
  # get true est
  hazard_0=as.vector(X[,-1]%*%alpha)
  hazard_1=as.vector(gamma.1+X[,-1]%*%alpha)
  hazard_2=as.vector(gamma.2+X[,-1]%*%alpha)
  Survival_time_0=(-log(random_simu)/(lambda*exp(hazard_0)))^(1/v)
  Survival_time_1=(-log(random_simu)/(lambda*exp(hazard_1)))^(1/v)
  Survival_time_2=(-log(random_simu)/(lambda*exp(hazard_2)))^(1/v)
  
  # truncate=min(quantile(Survival_time_1[Z==1],0.75),quantile(Survival_time_0[Z==0],0.75))
  tilt.h = 1/rowSums(1/soft.max)
  if(!multi.arm){
    true_est_finite[i,1]=mean(Survival_time_1)-mean(Survival_time_0)
    true_est_ow_finite[i,1]=sum((Survival_time_1-Survival_time_0)*tilt.h)/sum(tilt.h)
    
    # true_est_finite[i,3]=mean(Survival_time_1>truncate)-mean(Survival_time_0>truncate)
    # true_est_ow_finite[i,3]=sum((as.numeric(Survival_time_1>truncate)-as.numeric(Survival_time_0>truncate))*tilt.h)/sum(tilt.h)
    # 
    # Survival_time_1[Survival_time_1>truncate]=truncate
    # Survival_time_0[Survival_time_0>truncate]=truncate
    # true_est_finite[i,2]=mean(Survival_time_1)-mean(Survival_time_0)
    # true_est_ow_finite[i,2]=sum((Survival_time_1-Survival_time_0)*tilt.h)/sum(tilt.h)
  }else{
    # alpha = c(-1,0,1)
    true_est_finite[i,1] = mean(Survival_time_2)-mean(Survival_time_0)
    true_est_ow_finite[i,1] = sum((Survival_time_2-Survival_time_0)*tilt.h)/sum(tilt.h)
    
    # true_est_finite[i,3]= mean(Survival_time_2>truncate)-mean(Survival_time_0>truncate)
    # true_est_ow_finite[i,3] = sum((as.numeric(Survival_time_2>truncate)-as.numeric(Survival_time_0>truncate))*tilt.h)/sum(tilt.h)
    # 
    # Survival_time_1[Survival_time_2>truncate] = truncate
    # Survival_time_0[Survival_time_0>truncate] = truncate
    # true_est_finite[i,2] = mean(Survival_time_2)-mean(Survival_time_0)
    # true_est_ow_finite[i,2] = sum((Survival_time_2-Survival_time_0)*tilt.h)/sum(tilt.h)
  }
  
  
  # Censor time
  Censor_time=runif(sample_size,min=0,max=115)
  # Censor_time=runif(sample_size,min=2000,max=2115)
  
  
  
  # Observed survival
  Censored=Survival_time>Censor_time
  Y=Survival_time
  Y[Censored]=Censor_time[Censored]
  DELTA=1-Censored
  if(multi.arm){
    alpha.arg = c(-1,0,1)
  }else{
    alpha.arg = c(-1,1) 
  }
  res.IPWC=PSW.pseudo(Y, Z, DELTA, X,var.method=2,weight.type="IPW",
             estimand.type="ASCE", ps.threshold = NA, alpha = alpha.arg)
  ipw_est[i,1] <- res.IPWC$tau
  ipw_se[i,1]<- res.IPWC$se

  res.OW=PSW.pseudo(Y, Z, DELTA, X,var.method=2,weight.type="OW",
                      estimand.type="ASCE", alpha = alpha.arg)

  ow_est[i,1] <- res.OW$tau
  ow_se[i,1]<- res.OW$se

  res.IPWC=PSW.pseudo(Y, Z, DELTA, X,var.method=2,weight.type="IPW",
                      estimand.type="RACE", ps.threshold = 0.03, alpha = alpha.arg,
                      evaluate.time = truncate)
  ipw_est[i,2] <- res.IPWC$tau
  ipw_se[i,2]<- res.IPWC$se

  res.OW=PSW.pseudo(Y, Z, DELTA, X,var.method=2,weight.type="OW",
                    estimand.type="RACE", alpha = alpha.arg,
                    evaluate.time = truncate)
  ow_est[i,2] <- res.OW$tau
  ow_se[i,2]<- res.OW$se

  res.IPWC=PSW.pseudo(Y, Z, DELTA, X,var.method=2,weight.type="IPW",
                      estimand.type="SPCE", ps.threshold = 0.03, alpha = alpha.arg,
                      evaluate.time = truncate)
  ipw_est[i,3] <- res.IPWC$tau
  ipw_se[i,3]<- res.IPWC$se

  res.OW=PSW.pseudo(Y, Z, DELTA, X,var.method=2,weight.type="OW",
                    estimand.type="SPCE", alpha = alpha.arg,
                    evaluate.time = truncate)
  ow_est[i,3] <- res.OW$tau
  ow_se[i,3]<- res.OW$se
  
  # ASCE
  # pseudo_obs=pseudomean(Y,event=DELTA)
  # # # IPW
  # res.IPWC <- IPWC(y.all=pseudo_obs, z.all=Z, W.all=X, q.all=0)
  # ipw_est[i,1] <- res.IPWC$TAU
  # ipw_se[i,1]<- res.IPWC$SE
  # # # OW
  # res.OW <- OW(y=pseudo_obs, z=Z, W=X)
  # ow_est[i,1] <- res.OW$tau
  # ow_se[i,1] <- res.OW$se
  # # 
  # unadj_est[i,1]=mean(pseudo_obs[Z==1])-mean(pseudo_obs[Z==0])
  # 
  # # RACE
  # pseudo_obs=pseudomean(Y,event=delta,tmax=truncate)
  # # IPW
  # res.IPWC <- IPWC(y.all=pseudo_obs, z.all=Z, W.all=X, q.all=0)
  # ipw_est[i,2] <- res.IPWC$TAU
  # ipw_se[i,2]<- res.IPWC$SE
  # # OW
  # res.OW <- OW(y=pseudo_obs, z=Z, W=X)
  # ow_est[i,2] <- res.OW$tau
  # ow_se[i,2] <- res.OW$se
  # 
  # unadj_est[i,2]=mean(pseudo_obs[Z==1])-mean(pseudo_obs[Z==0])
  # 
  # # SPCE
  # pseudo_obs=pseudosurv(Y,event=delta,tmax=truncate)
  # pseudo_obs=pseudo_obs$pseudo
  # # IPW
  # res.IPWC <- IPWC(y.all=pseudo_obs, z.all=Z, W.all=X, q.all=0)
  # ipw_est[i,3] <- res.IPWC$TAU
  # ipw_se[i,3]<- res.IPWC$SE
  # # OW
  # res.OW <- OW(y=pseudo_obs, z=Z, W=X)
  # ow_est[i,3] <- res.OW$tau
  # ow_se[i,3] <- res.OW$se
  
  # unadj_est[i,3]=mean(pseudo_obs[Z==1])-mean(pseudo_obs[Z==0])
  
  
  # {
  #   # ASCE
  #   pseudo_obs=numeric(length(Y))
  #   pseudo_obs[Z==1]=pseudomean(Y[Z==1],event=delta[Z==1])
  #   pseudo_obs[Z==0]=pseudomean(Y[Z==0],event=delta[Z==0])
  #   # IPW
  #   res.IPWC <- IPWC(y.all=pseudo_obs, z.all=Z, W.all=X, q.all=0)
  #   ipw_est_by_group[i,1] <- res.IPWC$TAU
  #   ipw_se_by_group[i,1]<- res.IPWC$SE
  #   # OW
  #   res.OW <- OW(y=pseudo_obs, z=Z, W=X)
  #   ow_est_by_group[i,1] <- res.OW$tau
  #   ow_se_by_group[i,1] <- res.OW$se
  #   
  #   # RACE
  #   pseudo_obs=numeric(length(Y))
  #   pseudo_obs[Z==1]=pseudomean(Y[Z==1],event=delta[Z==1],tmax=truncate)
  #   pseudo_obs[Z==0]=pseudomean(Y[Z==0],event=delta[Z==0],tmax=truncate)
  #   # IPW
  #   res.IPWC <- IPWC(y.all=pseudo_obs, z.all=Z, W.all=X, q.all=0)
  #   ipw_est_by_group[i,2] <- res.IPWC$TAU
  #   ipw_se_by_group[i,2]<- res.IPWC$SE
  #   # OW
  #   res.OW <- OW(y=pseudo_obs, z=Z, W=X)
  #   ow_est_by_group[i,2] <- res.OW$tau
  #   ow_se_by_group[i,2] <- res.OW$se
  #   
  #   # SPCE
  #   pseudo_obs=numeric(length(Y))
  #   pseudo_obs_obj=pseudosurv(Y[Z==1],event=delta[Z==1],tmax=truncate)
  #   pseudo_obs[Z==1]=pseudo_obs_obj$pseudo
  #   pseudo_obs_obj=pseudosurv(Y[Z==0],event=delta[Z==0],tmax=truncate)
  #   pseudo_obs[Z==0]=pseudo_obs_obj$pseudo
  #   
  #   # IPW
  #   res.IPWC <- IPWC(y.all=pseudo_obs, z.all=Z, W.all=X, q.all=0)
  #   ipw_est_by_group[i,3] <- res.IPWC$TAU
  #   ipw_se_by_group[i,3]<- res.IPWC$SE
  #   # OW
  #   res.OW <- OW(y=pseudo_obs, z=Z, W=X)
  #   ow_est_by_group[i,3] <- res.OW$tau
  #   ow_se_by_group[i,3] <- res.OW$se
  #   
  # 
  # }
  # 

  ### Mao's method
  # res.IPW = estimand_analysis(X = X[,-1],Z = Z,Y = Y,delta = delta, weight.type = "IPW",
  #                             t.trunc=truncate, tmax = 120)
  # res.OW = estimand_analysis(X = X[,-1],Z = Z,Y = Y,delta = delta, weight.type = "OVERLAP",
  #                            t.trunc=truncate, tmax = 120)
  # res.MW = estimand_analysis(X = X[,-1],Z = Z,Y = Y,delta = delta, weight.type = "MW",
  #                            t.trunc=truncate, tmax = 120)
  # res.UW = estimand_analysis(X = X[,-1],Z = Z,Y = Y,delta = delta, weight.type = "UNWEIGHT",
  #                            t.trunc=truncate, tmax = 120)
  # 
  # ipw_est_mao[i,]=res.IPW$ans[,2]
  # ow_est_mao[i,]=res.OW$ans[,2]
  # mw_est_mao[i,]=res.MW$ans[,2]
  # uw_est_mao[i,]=res.UW$ans[,2]
  # 
  # ipw_se_mao[i,]=res.IPW$ans[,3]
  # ow_se_mao[i,]=res.OW$ans[,3]
  # mw_se_mao[i,]=res.MW$ans[,3]
  # uw_se_mao[i,]=res.UW$ans[,3]
  
  # ## Cox Q model
  # simu_data=data.frame(cbind(X,Z))
  # cox_q_model=flexsurvreg(Surv(Observed_Survival_time,Death)~X1+X2+X3+X4+Z,data=simu_data,
  #                       dist="weibull")
  # treated_data=simu_data
  # treated_data$Z=1
  # control_data=simu_data
  # control_data$Z=0
  # ## G-Formula
  # if(estimand=="mean"){
  # cox_q_treated=summary(cox_q_model,type="mean",newdata=treated_data)
  # cox_q_control=summary(cox_q_model,type="mean",newdata=control_data)
  # }else if(estimand=="survprob"){
  # cox_q_treated=summary(cox_q_model,type="survival",newdata=treated_data,t=c(60))
  # cox_q_control=summary(cox_q_model,type="survival",newdata=control_data,t=c(60))
  # }else{
  # cox_q_treated=summary(cox_q_model,type="rmst",newdata=treated_data,t=c(60))
  # cox_q_control=summary(cox_q_model,type="rmst",newdata=control_data,t=c(60))
  # }
  # treated_surv=unlist(lapply(cox_q_treated,FUN=function(x){x[2]}))
  # control_surv=unlist(lapply(cox_q_control,FUN=function(x){x[2]}))
  # cox_q_est[i]=mean(treated_surv-control_surv)
  ## Bootstrap for SE
  
  print(paste("==",i,"=="))
  
  # if(i==1){
  # setwd("results_1217/")
  # }
  
  
  # save(i,ow_est,ipw_est,unadj_est,
  #      cox_q_est,ow_se,ipw_se,
  #      true_est,true_est_ow,
  #      true_est_finite,true_est_ow_finite,
  #      ipw_est_mao,ipw_se_mao,
  #      uw_est_mao,uw_se_mao,
  #      ow_est_mao,ow_se_mao,
  #      mw_est_mao,mw_se_mao,
  #      ow_est_by_group,
  #      ipw_est_by_group,
  #      ow_se_by_group,
  #      ipw_se_by_group,
  #      file=paste(good_overlap,sample_size,multi.arm,"new_result.RData",sep="_"))
}

print(paste("==Sample size",sample_size,"=="))
print(paste("==overlap",good_overlap,"=="))
#  }
# dir.create(file.path("results"), showWarnings = FALSE)


#  }
#}
# ow_mse=sqrt(ow_mse)
# ipw_mse=sqrt(ipw_mse)
# ow_mse_table=matrix(ow_mse,ncol=2,byrow = T)
# ow_mse_table=cbind(matrix(ow_mse_table[,1],ncol=3),matrix(ow_mse_table[,2],ncol=3))
# 
# ipw_mse_table=matrix(ipw_mse,ncol=2,byrow = T)
# ipw_mse_table=cbind(matrix(ipw_mse_table[,1],ncol=3),matrix(ipw_mse_table[,2],ncol=3))
# 
# ow_bias_table=matrix(ow_bias,ncol=2,byrow = T)
# ow_bias_table=cbind(matrix(ow_bias_table[,1],ncol=3),matrix(ow_bias_table[,2],ncol=3))
# 
# 
# ipw_bias_table=matrix(ipw_bias,ncol=2,byrow = T)
# ipw_bias_table=cbind(matrix(ipw_bias_table[,1],ncol=3),matrix(ipw_bias_table[,2],ncol=3))
# 
# ow_coverage_table=matrix(ow_coverage,ncol=2,byrow = T)
# ow_coverage_table=cbind(matrix(ow_coverage_table[,1],ncol=3),matrix(ow_coverage_table[,2],ncol=3))
# 
# ipw_coverage_table=matrix(ipw_coverage,ncol=2,byrow = T)
# ipw_coverage_table=cbind(matrix(ipw_coverage_table[,1],ncol=3),matrix(ipw_coverage_table[,2],ncol=3))
# 
# 
# plot_size=2.5
# pdf("Summary.pdf",height=3.3*plot_size,width=3*plot_size)
# m <- matrix(c(1,2,3,4,5,6,7,8,9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)
# layout(mat = m,heights = c(0.5,0.5,0.5,0.15))
# plot(sample_size_grid,ow_mse_table[,1],type='o',
#      xlab="Sample Size",
#      ylab='RMSE',
#      main="RMSE comparison, ASCE",
#      ylim=range(ow_mse_table,ipw_mse_table),
#      col="black",pch=15)
# lines(sample_size_grid,ow_mse_table[,4],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_mse_table[,1],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_mse_table[,4],type='o',col="blue",pch=19)
# 
# plot(sample_size_grid,ow_bias_table[,1],type='o',
#      xlab="Sample Size",
#      ylab='Absolute Bias',
#      main="Absolute bias comparison, ASCE",
#      ylim=range(ow_bias_table,ipw_bias_table),
#      col="black",pch=15)
# lines(sample_size_grid,ow_bias_table[,4],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_bias_table[,1],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_bias_table[,4],type='o',col="blue",pch=19)
# 
# 
# plot(sample_size_grid,ow_coverage_table[,1],type='o',
#      xlab="Sample Size",
#      ylab='Coverage Rate',
#      main="Coverage Rate for 95% CI, ASCE",
#      ylim=range(ow_coverage_table,ipw_coverage_table),
#      col="black",pch=15)
# abline(h=0.95,lty=2)
# lines(sample_size_grid,ow_coverage_table[,4],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_coverage_table[,1],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_coverage_table[,4],type='o',col="blue",pch=19)
# 
# plot(sample_size_grid,ow_mse_table[,3],type='o',
#      xlab="Sample Size",
#      ylab='RMSE',
#      main="RMSE comparison, RACE",
#      ylim=range(ow_mse_table,ipw_mse_table),
#      col="black",pch=15)
# lines(sample_size_grid,ow_mse_table[,6],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_mse_table[,3],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_mse_table[,6],type='o',col="blue",pch=19)
# 
# plot(sample_size_grid,ow_bias_table[,3],type='o',
#      xlab="Sample Size",
#      ylab='Absolute Bias',
#      main="Absolute bias comparison, RACE",
#      ylim=range(ow_bias_table,ipw_bias_table),
#      col="black",pch=15)
# lines(sample_size_grid,ow_bias_table[,6],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_bias_table[,3],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_bias_table[,6],type='o',col="blue",pch=19)
# 
# 
# plot(sample_size_grid,ow_coverage_table[,3],type='o',
#      xlab="Sample Size",
#      ylab='Coverage Rate',
#      main="Coverage Rate for 95% CI, RACE",
#      ylim=range(ow_coverage_table,ipw_coverage_table),
#      col="black",pch=15)
# abline(h=0.95,lty=2)
# lines(sample_size_grid,ow_coverage_table[,6],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_coverage_table[,3],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_coverage_table[,6],type='o',col="blue",pch=19)
# 
# 
# plot(sample_size_grid,ow_mse_table[,2],type='o',
#      xlab="Sample Size",
#      ylab='RMSE',
#      main="RMSE comparison, SPCE",
#      ylim=range(ow_mse_table[,2],ipw_mse_table[,5]),
#      col="black",pch=15)
# lines(sample_size_grid,ow_mse_table[,5],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_mse_table[,2],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_mse_table[,5],type='o',col="blue",pch=19)
# 
# plot(sample_size_grid,ow_bias_table[,2],type='o',
#      xlab="Sample Size",
#      ylab='Absolute Bias',
#      main="Absolute bias comparison, SPCE",
#      ylim=range(ow_bias_table[,2],ipw_bias_table[,5]),
#      col="black",pch=15)
# lines(sample_size_grid,ow_bias_table[,5],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_bias_table[,2],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_bias_table[,5],type='o',col="blue",pch=19)
# 
# 
# plot(sample_size_grid,ow_coverage_table[,2],type='o',
#      xlab="Sample Size",
#      ylab='Coverage Rate',
#      main="Coverage Rate for 95% CI, SPCE",
#      ylim=range(ow_coverage_table,ipw_coverage_table),
#      col="black",pch=15)
# abline(h=0.95,lty=2)
# lines(sample_size_grid,ow_coverage_table[,5],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_coverage_table[,2],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_coverage_table[,5],type='o',col="blue",pch=19)
# 
# par(mar = c(0.1,0.1,1,0.1))
# plot(1, type = "n", axes=FALSE, xlab="", ylab="")
# legend("top",inset=0,title="Method, Overlaps Condition", col=c("black","red","purple","blue"),lwd=1.5,
#        lty=1,pch=c(15,17,18,19),legend=c("OW, Good Overlaps","OW, Poor Overlaps","IPW, Good Overlaps","IPW, Poor Overlaps"),horiz=TRUE)
# dev.off()
coverage_rate_calc(ow_est[,1],ow_se[,1],mean(true_est_ow[,1]))
coverage_rate_calc(ow_est[,2],ow_se[,2],mean(true_est_ow[,2]))
coverage_rate_calc(ow_est[,3],ow_se[,3],mean(true_est_ow[,3]))

coverage_rate_calc(ow_est[,1],ow_se[,1],mean(true_est_ow[,1]))
coverage_rate_calc(ow_est[,2],ow_se[,2],mean(true_est_ow[,2]))
coverage_rate_calc(ow_est[,3],ow_se[,3],mean(true_est_ow[,3]))

sqrt(mean((unadj_est[,1]-true_est_finite[,1])^2))/mean(abs(true_est_finite[,1]))

sqrt(mean((ipw_est[,1]-true_est_finite[,1])^2))/mean(abs(true_est_finite[,1]))
sqrt(mean((ipw_est_by_group[,1]-true_est_finite[,1])^2))/mean(abs(true_est_finite[,1]))
sqrt(mean((ipw_est_mao[,1]-true_est_finite[,1])^2))/mean(abs(true_est_finite[,1]))

sqrt(mean((ow_est[,1]-true_est_ow_finite[,1])^2))/mean(abs(true_est_ow_finite[,1]))
sqrt(mean((ow_est_by_group[,1]-true_est_ow_finite[,1])^2))/mean(abs(true_est_ow_finite[,1]))
sqrt(mean((ow_est_mao[,1]-true_est_ow_finite[,1])^2))/mean(abs(true_est_ow_finite[,1]))

sqrt(mean((unadj_est[,2]-true_est_finite[,2])^2))/mean(abs(true_est_finite[,2]))

sqrt(mean((ipw_est[,2]-true_est_finite[,2])^2))/mean(abs(true_est_finite[,2]))
sqrt(mean((ipw_est_by_group[,2]-true_est_finite[,2])^2))/mean(abs(true_est_finite[,2]))
sqrt(mean((ipw_est_mao[,2]-true_est_finite[,2])^2))/mean(abs(true_est_finite[,2]))


sqrt(mean((ow_est[,2]-true_est_ow_finite[,2])^2))/mean(abs(true_est_ow_finite[,2]))
sqrt(mean((ow_est_by_group[,2]-true_est_ow_finite[,2])^2))/mean(abs(true_est_ow_finite[,2]))
sqrt(mean((ow_est_mao[,2]-true_est_ow_finite[,2])^2))/mean(abs(true_est_finite[,2]))



sqrt(mean((ipw_est[id,3]-true_est_finite[id,3])^2))/mean(abs(true_est_finite[id,3]))
sqrt(mean((ipw_est_by_group[,3]-true_est_finite[,3])^2))/mean(abs(true_est_finite[,3]))
sqrt(mean((ipw_est_mao[,3]-true_est_finite[,3])^2))/mean(abs(true_est_finite[,3]))


sqrt(mean((ow_est[,3]-true_est_ow_finite[,3])^2))/mean(abs(true_est_ow_finite[,3]))
sqrt(mean((ow_est_by_group[,3]-true_est_ow_finite[,3])^2))/mean(abs(true_est_ow_finite[,3]))
sqrt(mean((ow_est_mao[,3]-true_est_ow_finite[,3])^2))/mean(abs(true_est_finite[,3]))


# coverage_rate_calc(ow_est[,1],ow_se[,1],true_est_ow_finite[,1])
# coverage_rate_calc(ow_est[,2],ow_se[,2],true_est_ow_finite[,2])
# coverage_rate_calc(ow_est[,3],ow_se[,3],true_est_ow_finite[,3])
# 
# coverage_rate_calc(ipw_est[,1],ipw_se[,1],true_est_finite[,1])
# coverage_rate_calc(ipw_est[,2],ipw_se[,2],true_est_finite[,2])
# coverage_rate_calc(ipw_est[,3],ipw_se[,3],true_est_finite[,3])

