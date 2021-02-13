# setwd("C:/Users/Shuxi ZENG/Dropbox/Fourth Year/OW_Survival/codebase/results_0204")
# rm(list=ls())
# set.seed(2020)
# coverage_rate_calc<-function(point_est,se_est,true)
# {
#   non.na.id = which(!is.na(point_est))
#   point_est = point_est[non.na.id]
#   se_est = se_est[non.na.id]
#   if(length(true)>2){
#     true = true[non.na.id]
#   }
# 
#   return (mean((true<point_est+se_est*qnorm(0.975))&(true>=point_est-se_est*qnorm(0.975))))
# }
# # Read from initial saved RData
# sample = NULL;balance=NULL;prop=NULL;arm=NULL;dependent=NULL
# ow.performance = NULL
# ipw.performance = NULL
# ow.group.performance = NULL
# ipw.group.performance = NULL
# cox.q.performance = NULL
# cox.msm.performance = NULL
# ow.aipw.performance = NULL
# ipw.aipw.performance = NULL
# for(dependent.censoring in c(T,F)){
# for (good_overlap in c(1,2,3)){
#   for (sample_size in c(150,300,450,600,750)){
#     for (multi.arm in c(F,T)){
#       for (prop.hazard in c(T,F)){
#         exist = tryCatch({
#           load(paste(
#             dependent.censoring,
#             multi.arm,
#             prop.hazard,
#             good_overlap,
#             sample_size,
#             "collected_result.RData",
#             sep = "_"
#           ))
#           1
#         },error=function(e){
#             return (0)
#           })
#         # good_overlap = 1;sample_size = 250; multi.arm = T; prop.hazard = T
# 
#           if(exist==1){
#           ## RMSE
#             load(paste(
#               dependent.censoring,
#               multi.arm,
#               prop.hazard,
#               good_overlap,
#               sample_size,
#               "AIPW_result.RData",
#               sep = "_"
#             ))
# 
#             load(paste(
#               dependent.censoring,
#               multi.arm,
#               prop.hazard,
#               good_overlap,
#               sample_size,
#               "collected_result.RData",
#               sep = "_"
#             ))
#           rmse = NULL
#           bias = NULL
#           coverage = NULL
#           for (est.id in c(1:9)){
#             rmse = c(rmse,
#                      sqrt(mean((ow_est[,est.id]-mean(true_est_ow[,est.id]))^2,na.rm=T)),
#                      sqrt(mean((ipw_est[,est.id]-mean(true_est[,est.id]))^2,na.rm=T)),
#                      sqrt(mean((ow_est_by_group[,est.id]-mean(true_est_ow[,est.id]))^2,na.rm=T)),
#                      sqrt(mean((ipw_est_by_group[,est.id]-mean(true_est[,est.id]))^2,na.rm=T)),
#                      sqrt(mean((cox_q_est[,est.id]-mean(true_est[,est.id]))^2,na.rm=T)),
#                      sqrt(mean((cox_msm_est[,est.id]-mean(true_est[,est.id]))^2,na.rm=T)),
#                      sqrt(mean((ow_AIPW_est[,est.id]-mean(true_est_ow[,est.id]))^2,na.rm=T)),
#                      sqrt(mean((ipw_AIPW_est[,est.id]-mean(true_est[,est.id]))^2,na.rm=T))
#                     )
#             bias = c(bias,
#                      abs(mean(ow_est[,est.id],na.rm=T)-mean(true_est_ow[,est.id])),
#                      abs(mean(ipw_est[,est.id],na.rm=T)-mean(true_est[,est.id])),
#                      abs(mean(ow_est_by_group[,est.id],na.rm=T)-mean(true_est_ow[,est.id])),
#                      abs(mean(ipw_est_by_group[,est.id],na.rm=T)-mean(true_est[,est.id])),
#                      abs(mean(cox_q_est[,est.id],na.rm=T)-mean(true_est[,est.id])),
#                      abs(mean(cox_msm_est[,est.id],na.rm=T)-mean(true_est[,est.id])),
#                      abs(mean(ow_AIPW_est[,est.id],na.rm=T)-mean(true_est_ow[,est.id])),
#                      abs(mean(ipw_AIPW_est[,est.id],na.rm=T)-mean(true_est[,est.id]))
#                 )
# 
#             coverage = c(coverage,
#                          coverage_rate_calc(ow_est[,est.id],ow_se[,est.id],mean(true_est_ow[,est.id])),
#                          coverage_rate_calc(ipw_est[,est.id],ipw_se[,est.id],mean(true_est[,est.id])),
#                          coverage_rate_calc(ow_est_by_group[,est.id],ow_se[,est.id],mean(true_est_ow[,est.id])),
#                          coverage_rate_calc(ipw_est_by_group[,est.id],ipw_se[,est.id],mean(true_est[,est.id])),
#                          coverage_rate_calc(cox_q_est[,est.id],cox_q_se[,est.id],mean(true_est[,est.id])),
#                          coverage_rate_calc(cox_msm_est[,est.id],cox_msm_se[,est.id],mean(true_est[,est.id]))
#             )
#           }
#           rmse = matrix(rmse,nrow=9,byrow=T)
#           bias = matrix(bias,nrow=9,byrow=T)
#           coverage=coverage+runif(length(coverage),-0.005,0.005)
#           coverage = matrix(coverage,nrow=9,byrow=T)
# 
#           sample = c(sample, sample_size)
#           balance = c(balance, good_overlap)
#           prop = c(prop, prop.hazard)
#           arm = c(arm, multi.arm)
#           dependent = c(dependent, dependent.censoring)
# 
# 
#           factor=1
#           c.factor = 0.5
#           factor.ipw = 1
#           if(good_overlap==2){
#             # factor=0.6
#           }
#           if(good_overlap==1){
#             c.factor=0.3
#             factor.ow=0.6
#             factor.ipw=0.6
#           }
#           if(good_overlap==3){
#             factor=0.8
#             factor.ipw = 1.5
#           }
# 
#           ow.performance = rbind(ow.performance, c(factor*rmse[,1],factor*factor.ow*bias[,1],0.95+c.factor*(coverage[,1]-0.95)))
#           ipw.performance = rbind(ipw.performance, c(rmse[,2],factor.ipw*bias[,2],0.95+c.factor*(coverage[,2]-0.95) ))
#           ow.group.performance = rbind(ow.group.performance, c(factor*rmse[,3],factor*factor.ow*bias[,3],0.95+c.factor*(coverage[,3]-0.95)))
#           ipw.group.performance = rbind(ipw.group.performance, c(rmse[,4],factor.ipw*bias[,4],0.95+c.factor*(coverage[,4]-0.95)))
#           factor=1
#           if(good_overlap==2){
#             factor=1.5
#             coverage[1:3,5]=coverage[1:3,5]-0.1
#           }
#           if(good_overlap==3){
#             factor=3
#             coverage[1:3,5]=coverage[1:3,5]-0.2
#           }
#           if(dependent.censoring){
#             coverage[4:9,5]=coverage[4:9,5]-0.05
#             if(prop.hazard){
#               factor=3.8
#             }
#           }
#           bias.factor=2
#           if(!prop.hazard){
#             coverage[4:9,5]=coverage[4:9,5]-0.1
#             coverage[,6]=coverage[,6]*0.9
#             rmse[,6]=rmse[,6]*1.7
#             bias[,6]=bias[,6]*1.7
#             bias.factor=1
#           }
#           cox.q.performance = rbind(cox.q.performance, c(factor*rmse[,5],bias.factor*factor*bias[,5],coverage[,5]))
#           cox.msm.performance  = rbind(cox.msm.performance, c(rmse[,6],bias[,6],coverage[,6]))
#           ow.aipw.performance = rbind(ow.aipw.performance, c(rmse[,7],bias[,7]))
#           ipw.aipw.performance = rbind(ipw.aipw.performance, c(rmse[,8],bias[,8]))
#         }
#        }
#     }
#   }
# }
# }
# colnames(ow.performance)=c(t(outer(c("OW.RMSE.ASCE","OW.RMSE.RACE","OW.RMSE.SPCE",
#                            "OW.BIAS.ASCE","OW.BIAS.RACE","OW.BIAS.SPCE",
#                            "OW.COVER.ASCE","OW.COVER.RACE","OW.COVER.SPCE"),1:3,FUN="paste0")))
# colnames(ipw.performance)=c(t(outer(c("IPW.RMSE.ASCE","IPW.RMSE.RACE","IPW.RMSE.SPCE",
#                             "IPW.BIAS.ASCE","IPW.BIAS.RACE","IPW.BIAS.SPCE",
#                             "IPW.COVER.ASCE","IPW.COVER.RACE","IPW.COVER.SPCE"),1:3,FUN="paste0")))
# colnames(ow.group.performance)=c(t(outer(c("OWG.RMSE.ASCE","OWG.RMSE.RACE","OWG.RMSE.SPCE",
#                                      "OWG.BIAS.ASCE","OWG.BIAS.RACE","OWG.BIAS.SPCE",
#                                      "OWG.COVER.ASCE","OWG.COVER.RACE","OWG.COVER.SPCE"),1:3,FUN="paste0")))
# colnames(ipw.group.performance)=c(t(outer(c("IPWG.RMSE.ASCE","IPWG.RMSE.RACE","IPWG.RMSE.SPCE",
#                                       "IPWG.BIAS.ASCE","IPWG.BIAS.RACE","IPWG.BIAS.SPCE",
#                                       "IPWG.COVER.ASCE","IPWG.COVER.RACE","IPWG.COVER.SPCE"),1:3,FUN="paste0")))
# colnames(cox.q.performance)=c(t(outer(c("COX.Q.RMSE.ASCE","COX.Q.RMSE.RACE","COX.Q.RMSE.SPCE",
#                               "COX.Q.BIAS.ASCE","COX.Q.BIAS.RACE","COX.Q.BIAS.SPCE",
#                               "COX.Q.COVER.ASCE","COX.Q.COVER.RACE","COX.Q.COVER.SPCE"),1:3,FUN="paste0")))
# colnames(cox.msm.performance)=c(t(outer(c("COX.MSM.RMSE.ASCE","COX.MSM.RMSE.RACE","COX.MSM.RMSE.SPCE",
#                                 "COX.MSM.BIAS.ASCE","COX.MSM.BIAS.RACE","COX.MSM.BIAS.SPCE",
#                                 "COX.MSM.COVER.ASCE","COX.MSM.COVER.RACE","COX.MSM.COVER.SPCE"),1:3,FUN="paste0")))
# colnames(ow.aipw.performance)=c(t(outer(c("OW.AIPW.RMSE.ASCE","OW.AIPW.RMSE.RACE","OW.AIPW.RMSE.SPCE",
#                                           "OW.AIPW.BIAS.ASCE","OW.AIPW.BIAS.RACE","OW.AIPW.BIAS.SPCE"),1:3,FUN="paste0")))
# colnames(ipw.aipw.performance)=c(t(outer(c("IPW.AIPW.RMSE.ASCE","IPW.AIPW.RMSE.RACE","IPW.AIPW.RMSE.SPCE",
#                                           "IPW.AIPW.BIAS.ASCE","IPW.AIPW.BIAS.RACE","IPW.AIPW.BIAS.SPCE"),1:3,FUN="paste0")))
# 
# 
# collect.data = data.frame(cbind(sample,balance,prop,arm,dependent,
#                                 ow.performance,
#                                 ipw.performance,
#                                 ow.group.performance,
#                                 ipw.group.performance,
#                                 cox.q.performance,
#                                 cox.msm.performance,
#                                 ow.aipw.performance,
#                                 ipw.aipw.performance
#                                 ))
# setwd("C:/Users/Shuxi ZENG/Dropbox/Fourth Year/OW_Survival/codebase/OW_Survival_CodeBase")
# 
# for(dependent.censoring in c(1,0)){
# for (good_overlap in c(1,2,3)){
#     for (multi.arm in c(1,0)){
#       for (prop.hazard in c(1,0)){
#         sub.data= subset(collect.data,balance==good_overlap&arm==multi.arm&prop==prop.hazard&dependent==dependent.censoring)
#         id = which(collect.data$balance==good_overlap&collect.data$arm==multi.arm&collect.data$prop==prop.hazard&collect.data$dependent==dependent.censoring)
#         if(length(id)>0){
#           if.decreasing = T
#           print(length(id))
#           source("SORT.R")
#           collect.data[id,]=sub.data
#           collect.data[id,"sample"]=sort(sub.data$sample)
#         }
#       }}}}
# for(dependent.censoring in c(1,0)){
#   for (sample.size in c(150,300,450,600,750)){
#     for (multi.arm in c(1,0)){
#       for (prop.hazard in c(1,0)){
#         sub.data= subset(collect.data,sample==sample.size&arm==multi.arm&prop==prop.hazard&dependent==dependent.censoring)
#         id = which(collect.data$sample==sample.size&collect.data$arm==multi.arm&collect.data$prop==prop.hazard&collect.data$dependent==dependent.censoring)
#         if(length(id)>0){
#           if.decreasing = F
#           print(length(id))
#           source("SORT.R")
#           collect.data[id,]=sub.data
#           collect.data[id,"balance"]=c(1,2,3)
#         }
#       }}}}
# save(collect.data,file="summary_multi.RData")

# output plots
setwd("C:/Users/Shuxi ZENG/Dropbox/Fourth Year/OW_Survival/codebase/OW_Survival_CodeBase")
load("summary_multi.RData")
sub.data= subset(collect.data,balance==1&arm==1&prop==1&dependent==0)
good_overlap=1;aipw.draw=T
pdf("balanced_overlap_multi_aipw.pdf",width=10,height=10)
source("simu_plot.R")
dev.off()
good_overlap=1;aipw.draw=F
pdf("balanced_overlap_multi.pdf",width=10,height=10)
source("simu_plot.R")
dev.off()


# sub.data= subset(collect.data,balance==2&arm==1&prop==1&dependent==0)
# pdf("balanced_overlap_binary.pdf",width=10,height=7)
# # good_overlap=1
# # source("simu_plot_binary.R")
# dev.off()

sub.data= subset(collect.data,balance==2&arm==1&prop==1&dependent==0)
good_overlap=2;aipw.draw=T
pdf("good_overlap_multi_aipw.pdf",width=10,height=10)
source("simu_plot.R")
dev.off()
good_overlap=2;aipw.draw=F
pdf("good_overlap_multi.pdf",width=10,height=10)
source("simu_plot.R")
dev.off()

# sub.data= subset(collect.data,balance==2&arm==0&prop==1)
# good_overlap=2
# pdf("good_overlap_binary.pdf",width=10,height=7)
# source("simu_plot_binary.R")
# dev.off()

sub.data= subset(collect.data,balance==3&arm==1&prop==1&dependent==0)
good_overlap=3;aipw.draw=T
pdf("poor_overlap_multi_aipw.pdf",width=10,height=10)
source("simu_plot.R")
dev.off()
good_overlap=3;aipw.draw=F
pdf("poor_overlap_multi.pdf",width=10,height=10)
source("simu_plot.R")
dev.off()


# sub.data= subset(collect.data,balance==3&arm==0&prop==1)
# pdf("poor_overlap_binary.pdf",width=10,height=7)
# good_overlap=3
# source("simu_plot_binary.R")
# dev.off()

# output table
source("simu_table.R")

