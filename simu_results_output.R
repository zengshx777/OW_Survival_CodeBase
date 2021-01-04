# setwd("results_0103/")
good_overlap = 3;sample_size = 250; multi.arm = F; prop.hazard = T
load(paste(
  good_overlap,
  sample_size,
  multi.arm,
  prop.hazard,
  "collected_result.RData",
  sep = "_"
))

## RMSE
est.id = 1
sqrt(mean((ow_est[,est.id]-mean(true_est_ow[,est.id])),na.rm=T)^2)/abs(mean(true_est_ow[,est.id]))
sqrt(mean((ipw_est[,est.id]-mean(true_est[,est.id])),na.rm=T)^2)/abs(mean(true_est[,est.id]))
sqrt(mean((cox_q_est[,est.id]-mean(true_est[,est.id]))^2))/abs(mean(true_est[,est.id]))
sqrt(mean((ipw_est_mao[,est.id]-mean(true_est[,est.id]))^2))/abs(mean(true_est[,est.id]))
sqrt(mean((ow_est_mao[,est.id]-mean(true_est[,est.id]))^2))/abs(mean(true_est[,est.id]))
sqrt(mean((mw_est_mao[,est.id]-mean(true_est[,est.id]))^2))/abs(mean(true_est[,est.id]))
est.id = 2
sqrt(mean((ow_est[,est.id]-mean(true_est_ow[,est.id])),na.rm=T)^2)/abs(mean(true_est_ow[,est.id]))
sqrt(mean((ipw_est[,est.id]-mean(true_est[,est.id]))^2))/abs(mean(true_est[,est.id]))
sqrt(mean((cox_q_est[,est.id]-mean(true_est[,est.id]))^2))/abs(mean(true_est[,est.id]))
sqrt(mean((ipw_est_mao[,est.id]-mean(true_est[,est.id]))^2))/abs(mean(true_est[,est.id]))
sqrt(mean((ow_est_mao[,est.id]-mean(true_est[,est.id]))^2))/abs(mean(true_est[,est.id]))
sqrt(mean((mw_est_mao[,est.id]-mean(true_est[,est.id]))^2))/abs(mean(true_est[,est.id]))
est.id = 3
sqrt(mean((ow_est[,est.id]-mean(true_est_ow[,est.id])),na.rm=T)^2)/abs(mean(true_est_ow[,est.id]))
sqrt(mean((ipw_est[,est.id]-mean(true_est[,est.id])),na.rm=T)^2)/abs(mean(true_est[,est.id]))
sqrt(mean((cox_q_est[,est.id]-mean(true_est[,est.id]))^2))/abs(mean(true_est[,est.id]))
sqrt(mean((ipw_est_mao[,est.id]-mean(true_est[,est.id]))^2))/abs(mean(true_est[,est.id]))
sqrt(mean((ow_est_mao[,est.id]-mean(true_est[,est.id]))^2))/abs(mean(true_est[,est.id]))
sqrt(mean((mw_est_mao[,est.id]-mean(true_est[,est.id]))^2))/abs(mean(true_est[,est.id]))

## BIAS
est.id = 1
abs(mean(ow_est[,est.id],na.rm=T)-mean(true_est_ow[,est.id]))/abs(mean(true_est_ow[,est.id]))
abs(mean(ipw_est[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id]))
abs(mean(cox_q_est[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id]))
abs(mean(ipw_est_mao[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id]))
abs(mean(ow_est_mao[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id]))
abs(mean(mw_est_mao[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id]))
est.id = 2
abs(mean(ow_est[,est.id],na.rm=T)-mean(true_est_ow[,est.id]))/abs(mean(true_est_ow[,est.id]))
abs(mean(ipw_est[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id]))
abs(mean(cox_q_est[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id]))
abs(mean(ipw_est_mao[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id]))
abs(mean(ow_est_mao[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id]))
abs(mean(mw_est_mao[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id]))
est.id = 3
abs(mean(ow_est[,est.id],na.rm=T)-mean(true_est_ow[,est.id]))/abs(mean(true_est_ow[,est.id]))
abs(mean(ipw_est[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id]))
abs(mean(cox_q_est[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id]))
abs(mean(ipw_est_mao[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id]))
abs(mean(ow_est_mao[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id]))
abs(mean(mw_est_mao[,est.id],na.rm=T)-mean(true_est[,est.id]))/abs(mean(true_est[,est.id]))
## Covarage rate
est.id = 1
coverage_rate_calc(ow_est[,est.id],ow_se[,est.id],mean(true_est_ow[,est.id]))
coverage_rate_calc(ipw_est[,est.id],ipw_se[,est.id],mean(true_est[,est.id]))
coverage_rate_calc(cox_q_est[,est.id],cox_q_se[,est.id],mean(true_est[,est.id]))
coverage_rate_calc(ipw_est_mao[,est.id],ipw_se_mao[,est.id],mean(true_est[,est.id]))
coverage_rate_calc(ow_est_mao[,est.id],ipw_se_mao[,est.id],mean(true_est[,est.id]))
coverage_rate_calc(mw_est_mao[,est.id],ipw_se_mao[,est.id],mean(true_est[,est.id]))
est.id = 2
coverage_rate_calc(ow_est[,est.id],ow_se[,est.id],mean(true_est_ow[,est.id]))
coverage_rate_calc(ipw_est[,est.id],ipw_se[,est.id],mean(true_est[,est.id]))
coverage_rate_calc(cox_q_est[,est.id],cox_q_se[,est.id],mean(true_est[,est.id]))
coverage_rate_calc(ipw_est_mao[,est.id],ipw_se_mao[,est.id],mean(true_est[,est.id]))
coverage_rate_calc(ow_est_mao[,est.id],ipw_se_mao[,est.id],mean(true_est[,est.id]))
coverage_rate_calc(mw_est_mao[,est.id],ipw_se_mao[,est.id],mean(true_est[,est.id]))
est.id = 3
coverage_rate_calc(ow_est[,est.id],ow_se[,est.id],mean(true_est_ow[,est.id]))
coverage_rate_calc(ipw_est[,est.id],ipw_se[,est.id],mean(true_est[,est.id]))
coverage_rate_calc(cox_q_est[,est.id],cox_q_se[,est.id],mean(true_est[,est.id]))
coverage_rate_calc(ipw_est_mao[,est.id],ipw_se_mao[,est.id],mean(true_est[,est.id]))
coverage_rate_calc(ow_est_mao[,est.id],ipw_se_mao[,est.id],mean(true_est[,est.id]))
coverage_rate_calc(mw_est_mao[,est.id],ipw_se_mao[,est.id],mean(true_est[,est.id]))
