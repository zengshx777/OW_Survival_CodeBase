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
cox_msm_est=matrix(NA,n_simu,3)


ow_est_mao=matrix(NA,n_simu,3)
mw_est_mao=matrix(NA,n_simu,3)
ipw_est_mao=matrix(NA,n_simu,3)
uw_est_mao=matrix(NA,n_simu,3)

ow_se=matrix(NA,n_simu,3)
ipw_se=matrix(NA,n_simu,3)
cox_q_se=matrix(NA,n_simu,3)
cox_msm_se=matrix(NA,n_simu,3)
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
