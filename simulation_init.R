setwd("C:/Users/Shuxi ZENG/Dropbox/Fourth Year/OW_Survival/codebase")
rm(list=ls())
library(MASS)
library(mvtnorm)
library(pseudo)
library(flexsurv)
library(survival)
source("IPWC.R")
source("OW.R")
source("Mao_Method_func.R")
set.seed(2020)

n_simu=50;n_mc=1000
sample_size_grid=c(150,250,350,500)
good_overlap_grid=c(TRUE,FALSE)
estimand_grid=c("mean","survprob","restricted_mean")


ow_mse=NULL
ipw_mse=NULL
ow_bias=NULL
ipw_bias=NULL
ow_coverage=NULL
ipw_coverage=NULL
cox_q_bias=NULL
cox_q_mse=NULL

for (estimand in estimand_grid){
for (sample_size in sample_size_grid){
for (good_overlap in good_overlap_grid){
  
ow_est=numeric(n_simu)
ipw_est=numeric(n_simu)
unadj_est=numeric(n_simu)
cox_q_est=numeric(n_simu)

ow_se=numeric(n_simu)
ipw_se=numeric(n_simu)

true_est=numeric(n_mc)
true_est_ow=numeric(n_mc)

## Use MC to approximate true
for (i in 1:n_mc){
  X4=sample(c(0,1),sample_size,replace=T,prob=c(0.5,0.5))
  X3=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.6-0.2*X4[x],0.2*X4[x]+0.4))}))
  X=cbind(1,mvrnorm(n=sample_size,mu=c(0,0),Sigma=matrix(c(2,0.5,0.5,2),2,2)),X3,X4)
  colnames(X)=c("X0","X1","X2","X3","X4")
  # Propensity Score Model
  if (good_overlap){
    # beta=c(-1.2,0.85,0.9,0.45,-0.25)
    beta=c(0,0,0,0,0)
  }else{
    beta=c(-1.2,0.85,0.9,0.45,-0.25)
    # beta=c(0.4,3,2,-3.5,-2)
  }
  
  true_e=1/(1+exp(-X%*%beta))
  Z=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(1-true_e[x],true_e[x]))}))
  
  # Plot PS histogram
  # a<-hist(true_e[Z==1],breaks=50,freq=FALSE)
  # a$counts<-a$counts/sum(a$counts)
  # 
  # b<-hist(true_e[Z==0],breaks=50,add=T)
  # b$counts<-b$counts/sum(b$counts)
  # 
  # pdf("good.pdf",width=8,height=6)
  # plot(a,col=rgb(1,0,0,0.5),main="Good Overlaps",cex.main=0.9,xlab="Propensity Score",ylab="Density",xlim=c(0,1),ylim=c(0,0.1))
  # plot(b,add=TRUE,col=rgb(0,0,1,0.5),main="Good Overlaps",xlab="Propensity Score",ylab="Density")
  # legend("top",legend=c("Treated Group","Control Group"),col=c("red","blue"),pch=15,cex=0.7)
  # dev.off()
  
  # Survival time
  gamma=1;alpha=c(2,1.5,-1,1);lambda=0.0001;v=3
  hazard=gamma*Z+X[,-1]%*%alpha
  random_simu=runif(sample_size)
  Survival_time=(-log(random_simu)/(lambda*exp(hazard)))^(1/v)
  
  # get true est
  hazard_0=X[,-1]%*%alpha
  hazard_1=gamma+X[,-1]%*%alpha
  Survival_time_0=(-log(random_simu)/(lambda*exp(hazard_0)))^(1/v)
  Survival_time_1=(-log(random_simu)/(lambda*exp(hazard_1)))^(1/v)
  if (estimand=="mean"){
    true_est[i]=mean(Survival_time_1)-mean(Survival_time_0)
    true_est_ow[i]=sum((Survival_time_1-Survival_time_0)*true_e*(1-true_e))/sum(true_e*(1-true_e))
  }else if (estimand=="survprob"){
    true_est[i]=mean(Survival_time_1>60)-mean(Survival_time_0>60)
    true_est_ow[i]=sum((as.numeric(Survival_time_1>60)-as.numeric(Survival_time_0>60))*true_e*(1-true_e))/sum(true_e*(1-true_e))
  }else{
    Survival_time_1[Survival_time_1>100]=100
    Survival_time_0[Survival_time_0>100]=100
    true_est[i]=mean(Survival_time_1)-mean(Survival_time_0)
    true_est_ow[i]=sum((Survival_time_1-Survival_time_0)*true_e*(1-true_e))/sum(true_e*(1-true_e))
    
  }
  if (i%%5000==0)
  {
    print(paste("== MC",i,"=="))
  }
}


for (i in 1:n_simu){
## Lower dimension of covariates, same from Lunceford setting
X4=sample(c(0,1),sample_size,replace=T,prob=c(0.5,0.5))
X3=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.6-0.2*X4[x],0.2*X4[x]+0.4))}))
X=cbind(1,mvrnorm(n=sample_size,mu=c(0,0),Sigma=matrix(c(2,0.5,0.5,2),2,2)),X3,X4)
colnames(X)=c("X0","X1","X2","X3","X4")
# Propensity Score Model
if (good_overlap){
# beta=c(-1.2,0.85,0.9,0.45,-0.25)
beta=c(0,0,0,0,0)
}else{
beta=c(-1.2,0.85,0.9,0.45,-0.25)
# beta=c(0.4,3,2,-3.5,-2)
}

true_e=1/(1+exp(-X%*%beta))
Z=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(1-true_e[x],true_e[x]))}))

# Plot PS histogram
# ps_model=glm(Z~X,family = binomial(link="logit"))
# ps_e=fitted(ps_model)
# a<-hist(ps_e[Z==1],breaks=50,freq=FALSE)
# a$counts<-a$counts/sum(a$counts)
# 
# b<-hist(ps_e[Z==0],breaks=50,add=T)
# b$counts<-b$counts/sum(b$counts)
# if(good_overlap){
# pdf("good_new.pdf",width=8,height=6)
# plot(a,col=rgb(1,0,0,0.5),main="Good Overlaps",cex.main=0.9,xlab="Propensity Score",ylab="Density",xlim=c(0.3,0.7),ylim=c(0,0.1))
# plot(b,add=TRUE,col=rgb(0,0,1,0.5),main="Good Overlaps",xlab="Propensity Score",ylab="Density")
# legend("top",legend=c("Treated Group","Control Group"),col=c("red","blue"),pch=15,cex=0.7)
# dev.off()}else{
# pdf("poor_new.pdf",width=8,height=6)
# plot(a,col=rgb(1,0,0,0.5),main="Poor Overlaps",cex.main=0.9,xlab="Propensity Score",ylab="Density",xlim=c(0,1),ylim=c(0,0.1))
# plot(b,add=TRUE,col=rgb(0,0,1,0.5),main="Poor Overlaps",xlab="Propensity Score",ylab="Density")
# legend("top",legend=c("Treated Group","Control Group"),col=c("red","blue"),pch=15,cex=0.7)
# dev.off()}

# Survival time
gamma=1;alpha=c(2,1.5,-1,1);lambda=0.0001;v=3
hazard=as.vector(gamma*Z+X[,-1]%*%alpha)
random_simu=runif(sample_size)
Survival_time=(-log(random_simu)/(lambda*exp(hazard)))^(1/v)

# get true est
hazard_0=X[,-1]%*%alpha
hazard_1=gamma+X[,-1]%*%alpha
Survival_time_0=(-log(random_simu)/(lambda*exp(hazard_0)))^(1/v)
Survival_time_1=(-log(random_simu)/(lambda*exp(hazard_1)))^(1/v)
# if (estimand=="mean"){
# true_est[i]=mean(Survival_time_1)-mean(Survival_time_0)
# true_est_ow[i]=sum((Survival_time_1-Survival_time_0)*true_e*(1-true_e))/sum(true_e*(1-true_e))
# }else if (estimand=="survprob"){
# true_est[i]=mean(Survival_time_1>60)-mean(Survival_time_0>60)
# true_est_ow[i]=sum((as.numeric(Survival_time_1>60)-as.numeric(Survival_time_0>60))*true_e*(1-true_e))/sum(true_e*(1-true_e))
# }else{
#   Survival_time_1[Survival_time_1>100]=100
#   Survival_time_0[Survival_time_0>100]=100
#   true_est[i]=mean(Survival_time_1)-mean(Survival_time_0)
#   true_est_ow[i]=sum((Survival_time_1-Survival_time_0)*true_e*(1-true_e))/sum(true_e*(1-true_e))
#   
# }

# Censor time
Censor_time=runif(sample_size,min=0,max=115)
# Censor_time=runif(sample_size,min=2000,max=2115)



# Observed survival
Censored=Survival_time>Censor_time
Y=Survival_time
Y[Censored]=Censor_time[Censored]
delta=1-Censored

### Mao's method
res.IPW = estimand_analysis(X = X[,-1],Z = Z,Y = Y,delta = delta, weight.type = "UNWEIGHT",
                        t.trunc=60, tmax = 120)

if (estimand=="mean"){
  pseudo_obs=pseudomean(Y,event=delta)
}else if(estimand=="survprob"){
  pseudo_obs=pseudosurv(Y,event=delta,tmax=60)
  pseudo_obs=pseudo_obs$pseudo
}else{
  pseudo_obs=pseudomean(Y,event=delta,tmax=60)
}
# ps_model=glm(Z~X,family = binomial(link="logit"))
# ps=fitted(ps_model)
# ow_weight=ps*(1-Z)+(1-ps)*Z
# ow_weight[Z==1]=ow_weight[Z==1]/mean(ow_weight[Z==1])
# ow_weight[Z==0]=ow_weight[Z==0]/mean(ow_weight[Z==0])
# 
# ipw_weight=(1-Z)/(1-ps)+Z/ps
# ipw_weight[Z==1]=ipw_weight[Z==1]/mean(ipw_weight[Z==1])
# ipw_weight[Z==0]=ipw_weight[Z==0]/mean(ipw_weight[Z==0])


# IPW
res.IPWC <- IPWC(y.all=pseudo_obs, z.all=Z, W.all=X, q.all=0)
ipw_est[i] <- res.IPWC$TAU
ipw_se[i]<- res.IPWC$SE



res.OW <- OW(y=pseudo_obs, z=Z, W=X)
ow_est[i] <- res.OW$tau
ow_se[i] <- res.OW$se

unadj_est[i]=mean(pseudo_obs[Z==1])-mean(pseudo_obs[Z==0])


# ## Cox Q model
# simu_data=data.frame(cbind(X,Z))
# cox_q_model=flexsurvreg(Surv(Observed_Survival_time,Death)~X1+X2+X3+X4+Z,data=simu_data,
#                         dist="weibull")
# treated_data=simu_data
# treated_data$Z=1
# control_data=simu_data
# control_data$Z=0
# ## G-Formula
# if(estimand=="mean"){
#   cox_q_treated=summary(cox_q_model,type="mean",newdata=treated_data)
#   cox_q_control=summary(cox_q_model,type="mean",newdata=control_data)
# }else if(estimand=="survprob"){
#   cox_q_treated=summary(cox_q_model,type="survival",newdata=treated_data,t=c(60))
#   cox_q_control=summary(cox_q_model,type="survival",newdata=control_data,t=c(60))
# }else{
#   cox_q_treated=summary(cox_q_model,type="rmst",newdata=treated_data,t=c(60))
#   cox_q_control=summary(cox_q_model,type="rmst",newdata=control_data,t=c(60))
# }
# treated_surv=unlist(lapply(cox_q_treated,FUN=function(x){x[2]}))
# control_surv=unlist(lapply(cox_q_control,FUN=function(x){x[2]}))
# cox_q_est[i]=mean(treated_surv-control_surv)
## Bootstrap for SE



if (i%%10==0)
{
  print(paste("==",i,"=="))
}
}

# sample_size=1000
# ## Lower dimension of covariates, same from Lunceford setting
# X4=sample(c(0,1),sample_size,replace=T,prob=c(0.5,0.5))
# X3=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.6-0.2*X4[x],0.2*X4[x]+0.4))}))
# X=cbind(1,mvrnorm(n=sample_size,mu=c(0,0),Sigma=matrix(c(2,0.5,0.5,2),2,2)),X3,X4)
# 
# # Propensity Score Model
# 
# # Survival time
# gamma=1;alpha=c(2,1.5,-1,1);lambda=0.0001;v=3
# hazard=gamma*1+X[,-1]%*%alpha
# Survival_time_1=(-log(runif(sample_size))/(lambda*exp(hazard)))^(1/v)
# mean(Survival_time_1>=10)
# 
# hazard=gamma*0+X[,-1]%*%alpha
# Survival_time_0=(-log(runif(sample_size))/(lambda*exp(hazard)))^(1/v)
# 
# 
# Survival_time_1[Survival_time_1>100]=100
# Survival_time_0[Survival_time_0>100]=100
# mean(Survival_time_1)-mean(Survival_time_0)

# Summerize SE
# ow_mse=c(ow_mse,mean((ow_est-true_est_ow)^2))
# ipw_mse=c(ipw_mse,mean((ipw_est-true_est)^2))

# Summarize RMSE/true_tau

ow_est=ow_est[ow_est!=0]
ipw_est=ipw_est[ipw_est!=0]
cox_q_est=cox_q_est[cox_q_est!=0]
ow_se=ow_se[ow_se!=0]
ipw_se=ipw_se[ipw_se!=0]

ow_mse=c(ow_mse,sqrt(mean((ow_est-mean(true_est_ow))^2))/abs(mean(true_est_ow)))
ipw_mse=c(ipw_mse,sqrt(mean((ipw_est-mean(true_est))^2))/abs(mean(true_est)))
cox_q_mse=c(cox_q_mse,sqrt(mean((cox_q_est-mean(true_est))^2))/abs(mean(true_est)))

# unadj_mse=c(unadj_mse,mean(abs(unadj_est-mean(true_est)))/abs(mean(true_est)))

# Bias/true_tau
ow_bias=c(ow_bias,abs(mean(ow_est-mean(true_est_ow)))/abs(mean(true_est_ow)))
ipw_bias=c(ipw_bias,abs(mean(ipw_est-mean(true_est)))/abs(mean(true_est)))
cox_q_bias=c(cox_q_bias,abs(mean(cox_q_est-mean(true_est)))/abs(mean(true_est)))
# unadj_bias=c(unadj_bias,abs(mean(unadj_est-mean(true_est))))


# Coverage
ow_coverage=c(ow_coverage,mean(mean(true_est_ow) <= ow_est + qnorm(0.975)*ow_se & 
           mean(true_est_ow) >= ow_est - qnorm(0.975)*ow_se))

ipw_coverage=c(ipw_coverage,mean(mean(true_est) <= ipw_est + qnorm(0.975)*ipw_se & 
       mean(true_est) >= ipw_est - qnorm(0.975)*ipw_se))

print(paste("==Sample size",sample_size,"=="))
print(paste("==overlap",good_overlap,"=="))
print(paste("==estimand",estimand,"=="))
}
}
}
# ow_mse=sqrt(ow_mse)
# ipw_mse=sqrt(ipw_mse)
ow_mse_table=matrix(ow_mse,ncol=2,byrow = T)
ow_mse_table=cbind(matrix(ow_mse_table[,1],ncol=3),matrix(ow_mse_table[,2],ncol=3))

ipw_mse_table=matrix(ipw_mse,ncol=2,byrow = T)
ipw_mse_table=cbind(matrix(ipw_mse_table[,1],ncol=3),matrix(ipw_mse_table[,2],ncol=3))

ow_bias_table=matrix(ow_bias,ncol=2,byrow = T)
ow_bias_table=cbind(matrix(ow_bias_table[,1],ncol=3),matrix(ow_bias_table[,2],ncol=3))


ipw_bias_table=matrix(ipw_bias,ncol=2,byrow = T)
ipw_bias_table=cbind(matrix(ipw_bias_table[,1],ncol=3),matrix(ipw_bias_table[,2],ncol=3))

ow_coverage_table=matrix(ow_coverage,ncol=2,byrow = T)
ow_coverage_table=cbind(matrix(ow_coverage_table[,1],ncol=3),matrix(ow_coverage_table[,2],ncol=3))

ipw_coverage_table=matrix(ipw_coverage,ncol=2,byrow = T)
ipw_coverage_table=cbind(matrix(ipw_coverage_table[,1],ncol=3),matrix(ipw_coverage_table[,2],ncol=3))


plot_size=2.5
pdf("Summary.pdf",height=3.3*plot_size,width=3*plot_size)
m <- matrix(c(1,2,3,4,5,6,7,8,9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(0.5,0.5,0.5,0.15))
plot(sample_size_grid,ow_mse_table[,1],type='o',
     xlab="Sample Size",
     ylab='RMSE',
     main="RMSE comparison, ASCE",
     ylim=range(ow_mse_table,ipw_mse_table),
     col="black",pch=15)
lines(sample_size_grid,ow_mse_table[,4],type='o',col="red",pch=17)
lines(sample_size_grid,ipw_mse_table[,1],type='o',col="purple",pch=18)
lines(sample_size_grid,ipw_mse_table[,4],type='o',col="blue",pch=19)

plot(sample_size_grid,ow_bias_table[,1],type='o',
     xlab="Sample Size",
     ylab='Absolute Bias',
     main="Absolute bias comparison, ASCE",
     ylim=range(ow_bias_table,ipw_bias_table),
     col="black",pch=15)
lines(sample_size_grid,ow_bias_table[,4],type='o',col="red",pch=17)
lines(sample_size_grid,ipw_bias_table[,1],type='o',col="purple",pch=18)
lines(sample_size_grid,ipw_bias_table[,4],type='o',col="blue",pch=19)


plot(sample_size_grid,ow_coverage_table[,1],type='o',
     xlab="Sample Size",
     ylab='Coverage Rate',
     main="Coverage Rate for 95% CI, ASCE",
     ylim=range(ow_coverage_table,ipw_coverage_table),
     col="black",pch=15)
abline(h=0.95,lty=2)
lines(sample_size_grid,ow_coverage_table[,4],type='o',col="red",pch=17)
lines(sample_size_grid,ipw_coverage_table[,1],type='o',col="purple",pch=18)
lines(sample_size_grid,ipw_coverage_table[,4],type='o',col="blue",pch=19)

plot(sample_size_grid,ow_mse_table[,3],type='o',
     xlab="Sample Size",
     ylab='RMSE',
     main="RMSE comparison, RACE",
     ylim=range(ow_mse_table,ipw_mse_table),
     col="black",pch=15)
lines(sample_size_grid,ow_mse_table[,6],type='o',col="red",pch=17)
lines(sample_size_grid,ipw_mse_table[,3],type='o',col="purple",pch=18)
lines(sample_size_grid,ipw_mse_table[,6],type='o',col="blue",pch=19)

plot(sample_size_grid,ow_bias_table[,3],type='o',
     xlab="Sample Size",
     ylab='Absolute Bias',
     main="Absolute bias comparison, RACE",
     ylim=range(ow_bias_table,ipw_bias_table),
     col="black",pch=15)
lines(sample_size_grid,ow_bias_table[,6],type='o',col="red",pch=17)
lines(sample_size_grid,ipw_bias_table[,3],type='o',col="purple",pch=18)
lines(sample_size_grid,ipw_bias_table[,6],type='o',col="blue",pch=19)


plot(sample_size_grid,ow_coverage_table[,3],type='o',
     xlab="Sample Size",
     ylab='Coverage Rate',
     main="Coverage Rate for 95% CI, RACE",
     ylim=range(ow_coverage_table,ipw_coverage_table),
     col="black",pch=15)
abline(h=0.95,lty=2)
lines(sample_size_grid,ow_coverage_table[,6],type='o',col="red",pch=17)
lines(sample_size_grid,ipw_coverage_table[,3],type='o',col="purple",pch=18)
lines(sample_size_grid,ipw_coverage_table[,6],type='o',col="blue",pch=19)


plot(sample_size_grid,ow_mse_table[,2],type='o',
     xlab="Sample Size",
     ylab='RMSE',
     main="RMSE comparison, SPCE",
     ylim=range(ow_mse_table[,2],ipw_mse_table[,5]),
     col="black",pch=15)
lines(sample_size_grid,ow_mse_table[,5],type='o',col="red",pch=17)
lines(sample_size_grid,ipw_mse_table[,2],type='o',col="purple",pch=18)
lines(sample_size_grid,ipw_mse_table[,5],type='o',col="blue",pch=19)

plot(sample_size_grid,ow_bias_table[,2],type='o',
     xlab="Sample Size",
     ylab='Absolute Bias',
     main="Absolute bias comparison, SPCE",
     ylim=range(ow_bias_table[,2],ipw_bias_table[,5]),
     col="black",pch=15)
lines(sample_size_grid,ow_bias_table[,5],type='o',col="red",pch=17)
lines(sample_size_grid,ipw_bias_table[,2],type='o',col="purple",pch=18)
lines(sample_size_grid,ipw_bias_table[,5],type='o',col="blue",pch=19)


plot(sample_size_grid,ow_coverage_table[,2],type='o',
     xlab="Sample Size",
     ylab='Coverage Rate',
     main="Coverage Rate for 95% CI, SPCE",
     ylim=range(ow_coverage_table,ipw_coverage_table),
     col="black",pch=15)
abline(h=0.95,lty=2)
lines(sample_size_grid,ow_coverage_table[,5],type='o',col="red",pch=17)
lines(sample_size_grid,ipw_coverage_table[,2],type='o',col="purple",pch=18)
lines(sample_size_grid,ipw_coverage_table[,5],type='o',col="blue",pch=19)

par(mar = c(0.1,0.1,1,0.1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top",inset=0,title="Method, Overlaps Condition", col=c("black","red","purple","blue"),lwd=1.5,
       lty=1,pch=c(15,17,18,19),legend=c("OW, Good Overlaps","OW, Poor Overlaps","IPW, Good Overlaps","IPW, Poor Overlaps"),horiz=TRUE)
dev.off()
