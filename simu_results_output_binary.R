# # setwd("results_0103/")
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
# Read from initial saved RData
# sample = NULL;balance=NULL;prop=NULL;arm=NULL
# ow.performance = NULL
# ipw.performance = NULL
# cox.q.performance = NULL
# ipw.mao.performance = NULL
# ow.mao.performance = NULL
# mw.mao.performance = NULL
# cox.msm.performance = NULL
# 
# for (good_overlap in c(1,2,3)){
#   for (sample_size in c(150,250,350,400,500)){
#     for (multi.arm in c(F,T)){
#       for (prop.hazard in c(T,F)){
#         exist = tryCatch({
#           load(paste(
#             good_overlap,
#             sample_size,
#             multi.arm,
#             prop.hazard,
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
#               good_overlap,
#               sample_size,
#               multi.arm,
#               prop.hazard,
#               "collected_result.RData",
#               sep = "_"
#             ))
#           rmse = NULL
#           bias = NULL
#           coverage = NULL
#           for (est.id in c(1,2,3)){
#             rmse = c(rmse,sqrt(mean((ow_est[,est.id]-mean(true_est_ow[,est.id]))^2,na.rm=T))/abs(mean(true_est_ow[,est.id])),
#                      sqrt(mean((ipw_est[,est.id]-mean(true_est[,est.id]))^2,na.rm=T))/abs(mean(true_est[,est.id])),
#                      sqrt(mean((cox_q_est[,est.id]-mean(true_est[,est.id]))^2,na.rm=T))/abs(mean(true_est[,est.id])),
#                      sqrt(mean((ipw_est_mao[,est.id]-mean(true_est[,est.id]))^2,na.rm=T))/abs(mean(true_est[,est.id])),
#                      sqrt(mean((ow_est_mao[,est.id]-mean(true_est[,est.id]))^2,na.rm=T))/abs(mean(true_est[,est.id])),
#                      sqrt(mean((mw_est_mao[,est.id]-mean(true_est[,est.id]))^2,na.rm=T))/abs(mean(true_est[,est.id])),
#                      sqrt(mean((cox_msm_est[,est.id]-mean(true_est[,est.id]))^2,na.rm=T))/abs(mean(true_est[,est.id]))
#             )
#             bias = c(bias, abs(mean(ow_est[,est.id],na.rm=T)-mean(true_est_ow[,est.id]))/abs(mean(true_est_ow[,est.id])),
#                      abs(mean(ipw_est[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id])),
#                      abs(mean(cox_q_est[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id])),
#                      abs(mean(ipw_est_mao[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id])),
#                      abs(mean(ow_est_mao[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id])),
#                      abs(mean(mw_est_mao[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id])),
#                      abs(mean(cox_msm_est[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id]))
#             )
#             
#             coverage = c(coverage,coverage_rate_calc(ow_est[,est.id],ow_se[,est.id],mean(true_est_ow[,est.id])),
#                          coverage_rate_calc(ipw_est[,est.id],ipw_se[,est.id],mean(true_est[,est.id])),
#                          coverage_rate_calc(cox_q_est[,est.id],cox_q_se[,est.id],mean(true_est[,est.id])),
#                          coverage_rate_calc(ipw_est_mao[,est.id],ipw_se_mao[,est.id],mean(true_est[,est.id])),
#                          coverage_rate_calc(ow_est_mao[,est.id],ow_se_mao[,est.id],mean(true_est[,est.id])),
#                          coverage_rate_calc(mw_est_mao[,est.id],mw_se_mao[,est.id],mean(true_est[,est.id])),
#                          coverage_rate_calc(cox_msm_est[,est.id],cox_msm_se[,est.id],mean(true_est[,est.id]))
#             )
#           }
#           rmse = matrix(rmse,nrow=3,byrow=T)
#           bias = matrix(bias,nrow=3,byrow=T)
#           coverage = matrix(coverage,nrow=3,byrow=T)
#           
#           sample = c(sample, sample_size)
#           balance = c(balance, good_overlap)
#           prop = c(prop, prop.hazard)
#           arm = c(arm, multi.arm)
#           
#           ow.performance = rbind(ow.performance, c(rmse[,1],bias[,1],coverage[,1]))
#           ipw.performance = rbind(ipw.performance, c(rmse[,2],bias[,2],coverage[,2]))
#           cox.q.performance = rbind(cox.q.performance, c(rmse[,3],bias[,3],coverage[,3]))
#           ipw.mao.performance = rbind(ipw.mao.performance, c(rmse[,4],bias[,4],coverage[,4]))
#           ow.mao.performance = rbind(ow.mao.performance, c(rmse[,5],bias[,5],coverage[,5]))
#           mw.mao.performance  = rbind(mw.mao.performance, c(rmse[,6],bias[,6],coverage[,6]))
#           cox.msm.performance  = rbind(cox.msm.performance, c(rmse[,7],bias[,7],coverage[,7]))
#         }
#        }
#     }
#   }
# }
# colnames(ow.performance)=c("OW.RMSE.ASCE","OW.RMSE.RACE","OW.RMSE.SPCE",
#                            "OW.BIAS.ASCE","OW.BIAS.RACE","OW.BIAS.SPCE",
#                            "OW.COVER.ASCE","OW.COVER.RACE","OW.COVER.SPCE")
# colnames(ipw.performance)=c("IPW.RMSE.ASCE","IPW.RMSE.RACE","IPW.RMSE.SPCE",
#                             "IPW.BIAS.ASCE","IPW.BIAS.RACE","IPW.BIAS.SPCE",
#                             "IPW.COVER.ASCE","IPW.COVER.RACE","IPW.COVER.SPCE")
# colnames(cox.q.performance)=c("COX.Q.RMSE.ASCE","COX.Q.RMSE.RACE","COX.Q.RMSE.SPCE",
#                               "COX.Q.BIAS.ASCE","COX.Q.BIAS.RACE","COX.Q.BIAS.SPCE",
#                               "COX.Q.COVER.ASCE","COX.Q.COVER.RACE","COX.Q.COVER.SPCE")
# colnames(ipw.mao.performance)=c("IPW.MAO.RMSE.ASCE","IPW.MAO.RMSE.RACE","IPW.MAO.RMSE.SPCE",
#                                 "IPW.MAO.BIAS.ASCE","IPW.MAO.BIAS.RACE","IPW.MAO.BIAS.SPCE",
#                                 "IPW.MAO.COVER.ASCE","IPW.MAO.COVER.RACE","IPW.MAO.COVER.SPCE")
# colnames(ow.mao.performance)=c("OW.MAO.RMSE.ASCE","OW.MAO.RMSE.RACE","OW.MAO.RMSE.SPCE",
#                                "OW.MAO.BIAS.ASCE","OW.MAO.BIAS.RACE","OW.MAO.BIAS.SPCE",
#                                "OW.MAO.COVER.ASCE","OW.MAO.COVER.RACE","OW.MAO.COVER.SPCE")
# colnames(mw.mao.performance)=c("MW.MAO.RMSE.ASCE","MW.MAO.RMSE.RACE","MW.MAO.RMSE.SPCE",
#                                "MW.MAO.BIAS.ASCE","MW.MAO.BIAS.RACE","MW.MAO.BIAS.SPCE",
#                                "MW.MAO.COVER.ASCE","MW.MAO.COVER.RACE","MW.MAO.COVER.SPCE")
# colnames(cox.msm.performance)=c("COX.MSM.RMSE.ASCE","COX.MSM.RMSE.RACE","COX.MSM.RMSE.SPCE",
#                                 "COX.MSM.BIAS.ASCE","COX.MSM.BIAS.RACE","COX.MSM.BIAS.SPCE",
#                                 "COX.MSM.COVER.ASCE","COX.MSM.COVER.RACE","COX.MSM.COVER.SPCE")
# 
# 
# collect.data = data.frame(cbind(sample,balance,prop,arm,
#                                 ow.performance,
#                                 ipw.performance,
#                                 cox.q.performance,
#                                 cox.msm.performance,
#                                 ipw.mao.performance,
#                                 ow.mao.performance,
#                                 mw.mao.performance))
# 
# save(collect.data,file="summary.RData")

# output plots
load("summary.RData")
sub.data= subset(collect.data,balance==1&arm==0&prop==1)
pdf("balanced_overlap.pdf",width=10,height=10)
source("simu_plot.R")
dev.off()

sub.data= subset(collect.data,balance==1&arm==0&prop==1)
pdf("balanced_overlap_binary.pdf",width=10,height=7)
good_overlap=1
source("simu_plot_binary.R")
dev.off()

sub.data= subset(collect.data,balance==2&arm==0&prop==1)
pdf("good_overlap.pdf",width=10,height=10)
source("simu_plot.R")
dev.off()

# sub.data= subset(collect.data,balance==2&arm==0&prop==1)
# good_overlap=2
# pdf("good_overlap_binary.pdf",width=10,height=7)
# source("simu_plot_binary.R")
# dev.off()

sub.data= subset(collect.data,balance==3&arm==0&prop==1)
pdf("poor_overlap.pdf",width=10,height=10)
source("simu_plot.R")
dev.off()

sub.data= subset(collect.data,balance==3&arm==0&prop==1)
pdf("poor_overlap_binary.pdf",width=10,height=7)
good_overlap=3
source("simu_plot_binary.R")
dev.off()

# output table
source("simu_table.R")

