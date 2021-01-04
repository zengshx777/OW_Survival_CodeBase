# setwd("C:/Users/Shuxi ZENG/Dropbox/Fourth Year/OW_Survival/codebase/OW_Survival_CodeBase")
# setwd("~/OW_Survival/Codebase/OW_Survival_CodeBase")

rm(list = ls())
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  print("No arguments supplied.")
} else{
  for (i in 1:length(args)) {
    eval(parse(text = args[[i]]))
  }
}
# good_overlap = 1; sample_size = 200; multi.arm = F; prop.hazard = T
truncate.ph = 50
truncate.aft = 60

n_simu = 200
n_mc = 100000

mao.method = T
cox.q.method = T
cox.msm.method = T

library(MASS)
library(mvtnorm)
source("cox_model.R")
source("PSW_pseudo.R")
source("Mao_Method_func.R")
source("simu_utils.R")
set.seed(2020)

if (multi.arm) {
  alpha.arg = c(-1, 0, 1)
} else{
  alpha.arg = c(-1, 1)
}

#########
######### Use MC to approximate true
for (i in 1:n_mc) {
  ## Lower dimension of covariates, same from Lunceford setting
  X4 = sample(c(0, 1),
              sample_size,
              replace = T,
              prob = c(0.5, 0.5))
  X3 = unlist(lapply(
    1:sample_size,
    FUN = function(x) {
      sample(c(0, 1), 1, prob = c(0.6 - 0.2 * X4[x], 0.2 * X4[x] + 0.4))
    }
  ))
  X = cbind(1, mvrnorm(
    n = sample_size,
    mu = c(0, 0),
    Sigma = matrix(c(2, 0.5, 0.5, 2), 2, 2)
  ), X3, X4)
  colnames(X) = c("X0", "X1", "X2", "X3", "X4")
  # Propensity Score Model
  beta.0 = c(0, 0, 0, 0, 0)
  if (good_overlap == 1) {
    beta.1 = c(0, 0, 0, 0, 0)
    beta.2 = c(0, 0, 0, 0, 0)
  } else if (good_overlap == 2) {
    beta.1 = c(-0.4, 0.85, 0.9, 0.45, -0.25)
    beta.2 = beta.1 * 0.2
  } else{
    beta.1 = c(1.2, 1.5, 1, -1.5, -1)
    beta.2 = beta.1 * 0.2
  }
  
  if (multi.arm) {
    beta = cbind(beta.0, beta.1, beta.2)
  } else{
    beta = cbind(beta.0, beta.1)
  }
  exp.part = exp(X %*% beta)
  soft.max = exp.part / rowSums(exp.part)
  
  # true_e=1/(1+exp(-X%*%beta))
  Z = unlist(lapply(
    1:sample_size,
    FUN = function(x) {
      sample(0:(ncol(soft.max) - 1), 1, prob = soft.max[x, ])
    }
  ))
  
  
  # Survival time
  
  if (prop.hazard) {
    gamma.1 = 0.5
    gamma.2 = 0.5
    alpha = c(2, 1.5, -1, 1)
    lambda = 0.0001
    v = 3
    truncate = truncate.ph
    # The covariate and the treatment is defined on proportional hazard
    random_simu = runif(sample_size)
    # get true est
    hazard_0 = as.vector(X[, -1] %*% alpha)
    hazard_1 = as.vector(gamma.1 + X[, -1] %*% alpha)
    hazard_2 = as.vector(gamma.2 + X[, -1] %*% alpha)
    Survival_time_0 = (-log(random_simu) / (lambda * exp(hazard_0))) ^ (1 /
                                                                          v)
    Survival_time_1 = (-log(random_simu) / (lambda * exp(hazard_1))) ^ (1 /
                                                                          v)
    Survival_time_2 = (-log(random_simu) / (lambda * exp(hazard_2))) ^ (1 /
                                                                          v)
  } else{
    gamma.1 = -1
    gamma.2 = -1
    alpha = -c(2, 1.5, -1, 1)
    lambda = 0.0001
    v = 3
    truncate = truncate.aft
    # AFT but not PH; log-normal
    random_simu = rnorm(sample_size, sd = 0.81)
    Survival_time_0 = exp(3.5 + 0.2 * (X[, -1] %*% alpha) + random_simu)
    Survival_time_1 = exp(3.5 + 0.2 * (gamma.1 + X[, -1] %*% alpha) + random_simu)
    Survival_time_2 = exp(3.5 + 0.2 * (gamma.2 + X[, -1] %*% alpha) + random_simu)
  }
  
  
  # truncate=min(quantile(Survival_time_1[Z==1],0.75),quantile(Survival_time_0[Z==0],0.75))
  tilt.h = 1 / rowSums(1 / soft.max)
  if (!multi.arm) {
    true_est[i, 1] = mean(Survival_time_1) - mean(Survival_time_0)
    true_est_ow[i, 1] = sum((Survival_time_1 - Survival_time_0) * tilt.h) /
      sum(tilt.h)
    true_est[i, 3] = mean(Survival_time_1 > truncate) - mean(Survival_time_0 >
                                                               truncate)
    true_est_ow[i, 3] = sum((
      as.numeric(Survival_time_1 > truncate) - as.numeric(Survival_time_0 > truncate)
    ) * tilt.h) / sum(tilt.h)
    
    Survival_time_1[Survival_time_1 > truncate] = truncate
    Survival_time_0[Survival_time_0 > truncate] = truncate
    true_est[i, 2] = mean(Survival_time_1) - mean(Survival_time_0)
    true_est_ow[i, 2] = sum((Survival_time_1 - Survival_time_0) * tilt.h) /
      sum(tilt.h)
  } else{
    true_est[i, 1] = mean(Survival_time_2) - mean(Survival_time_0)
    true_est_ow[i, 1] = sum((Survival_time_2 - Survival_time_0) * tilt.h) /
      sum(tilt.h)
    
    true_est[i, 3] = mean(Survival_time_2 > truncate) - mean(Survival_time_0 >
                                                               truncate)
    true_est_ow[i, 3] = sum((
      as.numeric(Survival_time_2 > truncate) - as.numeric(Survival_time_0 > truncate)
    ) * tilt.h) / sum(tilt.h)
    
    Survival_time_2[Survival_time_2 > truncate] = truncate
    Survival_time_0[Survival_time_0 > truncate] = truncate
    true_est[i, 2] = mean(Survival_time_2) - mean(Survival_time_0)
    true_est_ow[i, 2] = sum((Survival_time_2 - Survival_time_0) * tilt.h) /
      sum(tilt.h)
  }
  
  
  
  if (i %% 5000 == 0)
  {
    print(paste("== MC", i, "=="))
  }
}
#########

for (i in (1:n_simu)) {
  source("simu_data_gen.R")
  
  ### PS Pseudo approach
  res.IPWC = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    var.method = 2,
    weight.type = "IPW",
    estimand.type = "ASCE",
    ps.threshold = NA,
    alpha = alpha.arg
  )
  ipw_est[i, 1] <- res.IPWC$tau
  ipw_se[i, 1] <- res.IPWC$se
  
  res.OW = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    var.method = 2,
    weight.type = "OW",
    estimand.type = "ASCE",
    alpha = alpha.arg
  )
  
  ow_est[i, 1] <- res.OW$tau
  ow_se[i, 1] <- res.OW$se
  
  res.IPWC = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    var.method = 2,
    weight.type = "IPW",
    estimand.type = "RACE",
    ps.threshold = 0.03,
    alpha = alpha.arg,
    evaluate.time = truncate
  )
  ipw_est[i, 2] <- res.IPWC$tau
  ipw_se[i, 2] <- res.IPWC$se
  
  res.OW = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    var.method = 2,
    weight.type = "OW",
    estimand.type = "RACE",
    alpha = alpha.arg,
    evaluate.time = truncate
  )
  ow_est[i, 2] <- res.OW$tau
  ow_se[i, 2] <- res.OW$se
  
  res.IPWC = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    var.method = 2,
    weight.type = "IPW",
    estimand.type = "SPCE",
    ps.threshold = 0.03,
    alpha = alpha.arg,
    evaluate.time = truncate
  )
  ipw_est[i, 3] <- res.IPWC$tau
  ipw_se[i, 3] <- res.IPWC$se
  
  res.OW = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    var.method = 2,
    weight.type = "OW",
    estimand.type = "SPCE",
    alpha = alpha.arg,
    evaluate.time = truncate
  )
  ow_est[i, 3] <- res.OW$tau
  ow_se[i, 3] <- res.OW$se
  
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
  
  ### Mao's method
  if (mao.method) {
    if (!multi.arm) {
      tmax = min(max(Y[Z == 1]), max(Y[Z == 0]))
      ## IPW based
      res.mao.IPW = estimand_analysis(
        X = X[, -1],
        Z = Z,
        Y = Y,
        delta = DELTA,
        weight.type = "IPW",
        t.trunc = truncate,
        tmax = tmax
      )
      ## OW based
      res.mao.OW = estimand_analysis(
        X = X[, -1],
        Z = Z,
        Y = Y,
        delta = DELTA,
        weight.type = "OVERLAP",
        t.trunc = truncate,
        tmax = tmax
      )
      ## MW based
      res.mao.MW = estimand_analysis(
        X = X[, -1],
        Z = Z,
        Y = Y,
        delta = DELTA,
        weight.type = "MW",
        t.trunc = truncate,
        tmax = tmax
      )
      ## Unweighted Cox
      # res.mao.UW = estimand_analysis(
      #   X = X[, -1],
      #   Z = Z,
      #   Y = Y,
      #   delta = DELTA,
      #   weight.type = "UNWEIGHT",
      #   t.trunc = truncate,
      #   tmax = tmax
      # )
    } else{
      tmax = min(max(Y[Z == 0]), max(Y[Z == 2]))
      X.aug = cbind(X, as.numeric(Z == 1))
      colnames(X.aug) = c("X0", "X1", "X2", "X3", "X4", "Z01")
      Z.aug = as.numeric(Z == 2)
      ## IPW based
      res.mao.IPW = estimand_analysis(
        X = X.aug[, -1],
        Z = Z.aug,
        Y = Y,
        delta = DELTA,
        weight.type = "IPW",
        t.trunc = truncate,
        tmax = tmax
      )
      ## OW based
      res.mao.OW = estimand_analysis(
        X = X.aug[, -1],
        Z = Z.aug,
        Y = Y,
        delta = DELTA,
        weight.type = "OVERLAP",
        t.trunc = truncate,
        tmax = tmax
      )
      ## MW based
      res.mao.MW = estimand_analysis(
        X = X.aug[, -1],
        Z = Z.aug,
        Y = Y,
        delta = DELTA,
        weight.type = "MW",
        t.trunc = truncate,
        tmax = tmax
      )
      ## Unweighted Cox
      # res.mao.UW = estimand_analysis(
      #   X = X.aug[, -1],
      #   Z = Z.aug,
      #   Y = Y,
      #   delta = DELTA,
      #   weight.type = "UNWEIGHT",
      #   t.trunc = truncate,
      #   tmax = tmax
      # )
    }
    ipw_est_mao[i, ] = res.mao.IPW$ans[, 2]
    ow_est_mao[i, ] = res.mao.OW$ans[, 2]
    mw_est_mao[i, ] = res.mao.MW$ans[, 2]
    # uw_est_mao[i, ] = res.mao.UW$ans[, 2]
    
    ipw_se_mao[i, ] = res.mao.IPW$ans[, 3]
    ow_se_mao[i, ] = res.mao.OW$ans[, 3]
    mw_se_mao[i, ] = res.mao.MW$ans[, 3]
    # uw_se_mao[i, ] = res.mao.UW$ans[, 3]
  }
  
  ### Cox Q model
  if (cox.q.method) {
    res.cox.q = cox.q.model.fit(
      Y,
      DELTA,
      X,
      Z,
      alpha = alpha.arg,
      truncate =  truncate,
      boot.time = 250
    )
    cox_q_est[i, ] = res.cox.q$est
    cox_q_se[i, ] = res.cox.q$se
  }
  
  ## MSM model
  if (cox.msm.method) {
    res.cox.msm = cox.msm.model.fit(
      Y,
      DELTA,
      X,
      Z,
      alpha = alpha.arg,
      truncate =  truncate,
      boot.time = 250
    )
    cox_msm_est[i, ] = res.cox.msm$est
    cox_msm_se[i, ] = res.cox.msm$se
  }
  
  print(paste("==", i, "=="))
  
  save(
    i,
    ow_est,
    ipw_est,
    unadj_est,
    ow_se,
    ipw_se,
    true_est,
    true_est_ow,
    true_est_finite,
    true_est_ow_finite,
    ipw_est_mao,
    ipw_se_mao,
    uw_est_mao,
    uw_se_mao,
    ow_est_mao,
    ow_se_mao,
    mw_est_mao,
    mw_se_mao,
    cox_q_est,
    cox_q_se,
    cox_msm_est,
    cox_msm_se,
    ow_est_by_group,
    ipw_est_by_group,
    ow_se_by_group,
    ipw_se_by_group,
    file = paste(
      good_overlap,
      sample_size,
      multi.arm,
      prop.hazard,
      "collected_result.RData",
      sep = "_"
    )
  )
}

print(paste("==Sample size", sample_size, "=="))
print(paste("==overlap", good_overlap, "=="))
#  }
# dir.create(file.path("results"), showWarnings = FALSE)


#source("simu_results_analysis.R")
