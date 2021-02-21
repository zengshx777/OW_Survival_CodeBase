
ow_mse=NULL
ipw_mse=NULL
ow_bias=NULL
ipw_bias=NULL
ow_coverage=NULL
ipw_coverage=NULL
cox_q_bias=NULL
cox_q_mse=NULL

ow_est=matrix(NA,n_simu,9)
ipw_est=matrix(NA,n_simu,9)
unadj_est=matrix(NA,n_simu,9)
cox_q_est=matrix(NA,n_simu,9)
cox_msm_est=matrix(NA,n_simu,9)


ow_est_mao=matrix(NA,n_simu,9)
mw_est_mao=matrix(NA,n_simu,9)
ipw_est_mao=matrix(NA,n_simu,9)
uw_est_mao=matrix(NA,n_simu,9)

ow_se=matrix(NA,n_simu,9)
ipw_se=matrix(NA,n_simu,9)
cox_q_se=matrix(NA,n_simu,9)
cox_msm_se=matrix(NA,n_simu,9)
ow_se_mao=matrix(NA,n_simu,9)
mw_se_mao=matrix(NA,n_simu,9)
ipw_se_mao=matrix(NA,n_simu,9)
uw_se_mao=matrix(NA,n_simu,9)

true_est_finite=matrix(NA,n_simu,9)
true_est_ow_finite=matrix(NA,n_simu,9)
true_est=matrix(NA,n_mc,9)
true_est_ow=matrix(NA,n_mc,9)
true_mean=matrix(NA,n_mc,9)
true_mean_ow=matrix(NA,n_mc,9)

ow_est_by_group=matrix(NA,n_simu,9)
ipw_est_by_group=matrix(NA,n_simu,9)
ow_se_by_group=matrix(NA,n_simu,9)
ipw_se_by_group=matrix(NA,n_simu,9)
