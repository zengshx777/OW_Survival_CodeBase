rm(list = ls())
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  print("No arguments supplied.")
} else{
  for (i in 1:length(args)) {
    eval(parse(text = args[[i]]))
  }
}
# example specification
# multi.arm=T; prop.hazard=F; good_overlap=1; sample_size=600; dependent.censoring=F

# evaluation time
truncate.ph = 50
truncate.aft = 60

# simulation times
n_simu = 1000
# MC times for calculating the true value
n_mc = 50000

# some options to skip certain methods
mao.method = T
cox.q.method = T
cox.msm.method = T
# by.group.pseudo = T

library(MASS)
library(mvtnorm)
source("../estimators/cox_model.R")
source("../estimators/PSW_pseudo.R")
source("../estimators/Mao_Method_func.R")
source("../estimators/AIPW_pseudo.R")
source("../estimators/pseudo_G.R")
source("simu_utils.R")
set.seed(2020)


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
  true_est[i, 1] = mean(Survival_time_1) - mean(Survival_time_0)
  true_est[i, 2] = mean(Survival_time_2) - mean(Survival_time_0)
  true_est[i, 3] = mean(Survival_time_2) - mean(Survival_time_1)
  
  true_est[i, 4] = mean(pmin(Survival_time_1, truncate)) - mean(pmin(Survival_time_0, truncate))
  true_est[i, 5] = mean(pmin(Survival_time_2, truncate)) - mean(pmin(Survival_time_0, truncate))
  true_est[i, 6] = mean(pmin(Survival_time_2, truncate)) - mean(pmin(Survival_time_1, truncate))
  
  true_est[i, 7] =  mean(Survival_time_1 > truncate) - mean(Survival_time_0 >
                                                              truncate)
  true_est[i, 8] =  mean(Survival_time_2 > truncate) - mean(Survival_time_0 >
                                                              truncate)
  true_est[i, 9] =  mean(Survival_time_2 > truncate) - mean(Survival_time_1 >
                                                              truncate)
  
  true_est_ow[i, 1] = sum((Survival_time_1 - Survival_time_0) * tilt.h) /
    sum(tilt.h)
  true_est_ow[i, 2] = sum((Survival_time_2 - Survival_time_0) * tilt.h) /
    sum(tilt.h)
  true_est_ow[i, 3] = sum((Survival_time_2 - Survival_time_1) * tilt.h) /
    sum(tilt.h)
  
  true_est_ow[i, 4] = sum((
    pmin(Survival_time_1, truncate) - pmin(Survival_time_0, truncate)
  ) * tilt.h) /
    sum(tilt.h)
  true_est_ow[i, 5] = sum((
    pmin(Survival_time_2, truncate) - pmin(Survival_time_0, truncate)
  ) * tilt.h) /
    sum(tilt.h)
  true_est_ow[i, 6] = sum((
    pmin(Survival_time_2, truncate) - pmin(Survival_time_1, truncate)
  ) * tilt.h) /
    sum(tilt.h)
  
  true_est_ow[i, 7] = sum((
    as.numeric(Survival_time_1 > truncate) - as.numeric(Survival_time_0 > truncate)
  ) * tilt.h) / sum(tilt.h)
  true_est_ow[i, 8] = sum((
    as.numeric(Survival_time_2 > truncate) - as.numeric(Survival_time_0 > truncate)
  ) * tilt.h) / sum(tilt.h)
  true_est_ow[i, 9] = sum((
    as.numeric(Survival_time_2 > truncate) - as.numeric(Survival_time_1 > truncate)
  ) * tilt.h) / sum(tilt.h)
  
  if (i %% 5000 == 0)
  {
    print(paste("== MC", i, "=="))
  }
  
  true_mean[i,] = c(
    mean(Survival_time_0),
    mean(Survival_time_1),
    mean(Survival_time_2),
    mean(pmin(Survival_time_0, truncate)),
    mean(pmin(Survival_time_1, truncate)),
    mean(pmin(Survival_time_2, truncate)),
    mean(Survival_time_0 > truncate),
    mean(Survival_time_1 > truncate),
    mean(Survival_time_2 > truncate)
  )
  true_mean_ow[i,] = c(
    mean(Survival_time_0 * tilt.h) / mean(tilt.h),
    mean(Survival_time_1 * tilt.h) / mean(tilt.h),
    mean(Survival_time_2 * tilt.h) / mean(tilt.h),
    mean(pmin(Survival_time_0, truncate) * tilt.h) / mean(tilt.h),
    mean(pmin(Survival_time_1, truncate) * tilt.h) / mean(tilt.h),
    mean(pmin(Survival_time_2, truncate) * tilt.h) / mean(tilt.h),
    mean(as.numeric(Survival_time_0 > truncate) * tilt.h) /
      mean(tilt.h),
    mean(as.numeric(Survival_time_1 > truncate) * tilt.h) /
      mean(tilt.h),
    mean(as.numeric(Survival_time_2 > truncate) * tilt.h) /
      mean(tilt.h)
  )
}
#########
dir.create(file.path("../simulation_results"), showWarnings = FALSE)
for (i in (1:n_simu)) {
  source("simu_data_gen.R")
  ### PS Pseudo approach
  res.IPWC = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    weight.type = "IPW",
    ps.threshold = 0.03,
    estimand.type = "ASCE",
    dependent.adjustment = dependent.censoring
  )
  if (multi.arm) {
    ipw_est[i, 1:3] <- res.IPWC$tau
    ipw_se[i, 1:3] <- res.IPWC$se
  } else{
    ipw_est[i, 1:2] <- res.IPWC$tau
    ipw_se[i, 1:2] <- res.IPWC$se
  }
  
  res.OW = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    weight.type = "OW",
    estimand.type = "ASCE",
    dependent.adjustment = dependent.censoring
  )
  if (multi.arm) {
    ow_est[i, 1:3] <- res.OW$tau
    ow_se[i, 1:3] <- res.OW$se
  } else{
    ow_est[i, 1:2] <- res.OW$tau
    ow_se[i, 1:2] <- res.OW$se
  }
  
  
  res.IPWC = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    weight.type = "IPW",
    estimand.type = "RACE",
    ps.threshold = 0.03,
    evaluate.time = truncate,
    dependent.adjustment = dependent.censoring
  )
  if (multi.arm) {
    ipw_est[i, 4:6] <- res.IPWC$tau
    ipw_se[i, 4:6] <- res.IPWC$se
  } else{
    ipw_est[i, 4:5] <- res.IPWC$tau
    ipw_se[i, 4:5] <- res.IPWC$se
  }
  
  
  res.OW = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    
    weight.type = "OW",
    estimand.type = "RACE",
    evaluate.time = truncate,
    dependent.adjustment = dependent.censoring
  )
  if (multi.arm) {
    ow_est[i, 4:6] <- res.OW$tau
    ow_se[i, 4:6] <- res.OW$se
  } else{
    ow_est[i, 4:5] <- res.OW$tau
    ow_se[i, 4:5] <- res.OW$se
  }
  
  
  res.IPWC = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    weight.type = "IPW",
    estimand.type = "SPCE",
    ps.threshold = 0.03,
    evaluate.time = truncate,
    dependent.adjustment = dependent.censoring
  )
  if (multi.arm) {
    ipw_est[i, 7:9] <- res.IPWC$tau
    ipw_se[i, 7:9] <- res.IPWC$se
  } else{
    ipw_est[i, 7:8] <- res.IPWC$tau
    ipw_se[i, 7:8] <- res.IPWC$se
  }
  
  
  res.OW = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    weight.type = "OW",
    estimand.type = "SPCE",
    evaluate.time = truncate,
    dependent.adjustment = dependent.censoring
  )
  if (multi.arm) {
    ow_est[i, 7:9] <- res.OW$tau
    ow_se[i, 7:9] <- res.OW$se
  } else{
    ow_est[i, 7:8] <- res.OW$tau
    ow_se[i, 7:8] <- res.OW$se
  }
  ### PO-UNADJ
  res.pseudo.unadj = unadj.pseudo(Y,
                                  Z,
                                  DELTA,
                                  X,
                                  estimand.type = "ASCE",
                                  dependent.adjustment = dependent.censoring)
  if (multi.arm) {
    pseudo_unadj_est[i, 1:3] <- res.pseudo.unadj$tau
    pseudo_unadj_se[i, 1:3] <- res.pseudo.unadj$se
  } else{
    pseudo_unadj_est[i, 1:2] <- res.pseudo.unadj$tau
    pseudo_unadj_se[i, 1:2] <- res.pseudo.unadj$se
  }
  
  res.pseudo.unadj = unadj.pseudo(
    Y,
    Z,
    DELTA,
    X,
    estimand.type = "RACE",
    evaluate.time = truncate,
    dependent.adjustment = dependent.censoring
  )
  
  if (multi.arm) {
    pseudo_unadj_est[i, 4:6] <- res.pseudo.unadj$tau
    pseudo_unadj_se[i, 4:6] <- res.pseudo.unadj$se
  } else{
    pseudo_unadj_est[i, 4:5] <- res.pseudo.unadj$tau
    pseudo_unadj_se[i, 4:5] <- res.pseudo.unadj$se
  }
  
  res.pseudo.unadj = unadj.pseudo(
    Y,
    Z,
    DELTA,
    X,
    estimand.type = "SPCE",
    evaluate.time = truncate,
    dependent.adjustment = dependent.censoring
  )
  
  if (multi.arm) {
    pseudo_unadj_est[i, 7:9] <- res.pseudo.unadj$tau
    pseudo_unadj_se[i, 7:9] <- res.pseudo.unadj$se
  } else{
    pseudo_unadj_est[i, 7:9] <- res.pseudo.unadj$tau
    pseudo_unadj_se[i, 7:8] <- res.pseudo.unadj$se
  }
  
  ### PO-G
  res.pseudo.g = G.pseudo(Y,
                          Z,
                          DELTA,
                          X,
                          estimand.type = "ASCE",
                          dependent.adjustment = dependent.censoring)
  if (multi.arm) {
    pseudo_g_est[i, 1:3] <- res.pseudo.g$tau
    pseudo_g_se[i, 1:3] <- res.pseudo.g$se
  } else{
    pseudo_g_est[i, 1:2] <- res.pseudo.g$tau
    pseudo_g_se[i, 1:2] <- res.pseudo.g$se
  }
  
  
  res.pseudo.g = G.pseudo(
    Y,
    Z,
    DELTA,
    X,
    estimand.type = "RACE",
    evaluate.time = truncate,
    dependent.adjustment = dependent.censoring
  )
  if (multi.arm) {
    pseudo_g_est[i, 4:6] <- res.pseudo.g$tau
    pseudo_g_se[i, 4:6] <- res.pseudo.g$se
  } else{
    pseudo_g_est[i, 4:5] <- res.pseudo.g$tau
    pseudo_g_se[i, 4:5] <- res.pseudo.g$se
  }
  
  res.pseudo.g = G.pseudo(
    Y,
    Z,
    DELTA,
    X,
    estimand.type = "SPCE",
    evaluate.time = truncate,
    dependent.adjustment = dependent.censoring
  )
  if (multi.arm) {
    pseudo_g_est[i, 7:9] <- res.pseudo.g$tau
    pseudo_g_se[i, 7:9] <- res.pseudo.g$se
  } else{
    pseudo_g_est[i, 7:8] <- res.pseudo.g$tau
    pseudo_g_se[i, 7:8] <- res.pseudo.g$se
  }
  
  ### AIPW, AOW
  res.AIPW.IPWC = AIPW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    weight.type = "IPW",
    ps.threshold = 0.03,
    estimand.type = "ASCE",
    dependent.adjustment = dependent.censoring
  )
  if (multi.arm) {
    ipw_AIPW_est[i, 1:3] <- res.AIPW.IPWC$tau
  } else{
    ipw_AIPW_est[i, 1:2] <- res.AIPW.IPWC$tau
  }
  
  res.AIPW.OW = AIPW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    weight.type = "OW",
    estimand.type = "ASCE",
    dependent.adjustment = dependent.censoring
  )
  if (multi.arm) {
    ow_AIPW_est[i, 1:3] <- res.AIPW.OW$tau
  } else{
    ow_AIPW_est[i, 1:2] <- res.AIPW.OW$tau
  }
  
  res.AIPW.IPWC = AIPW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    weight.type = "IPW",
    ps.threshold = 0.03,
    estimand.type = "RACE",
    evaluate.time = truncate,
    dependent.adjustment = dependent.censoring
  )
  if (multi.arm) {
    ipw_AIPW_est[i, 4:6] <- res.AIPW.IPWC$tau
  } else{
    ipw_AIPW_est[i, 4:5] <- res.AIPW.IPWC$tau
  }
  
  res.AIPW.OW = AIPW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    weight.type = "OW",
    estimand.type = "RACE",
    evaluate.time = truncate,
    dependent.adjustment = dependent.censoring
  )
  if (multi.arm) {
    ow_AIPW_est[i, 4:6] <- res.AIPW.OW$tau
  } else{
    ow_AIPW_est[i, 4:5] <- res.AIPW.OW$tau
  }
  
  res.AIPW.IPWC = AIPW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    weight.type = "IPW",
    ps.threshold = 0.03,
    estimand.type = "SPCE",
    evaluate.time = truncate,
    dependent.adjustment = dependent.censoring
  )
  if (multi.arm) {
    ipw_AIPW_est[i, 7:9] <- res.AIPW.IPWC$tau
  } else{
    ipw_AIPW_est[i, 7:8] <- res.AIPW.IPWC$tau
  }
  
  res.AIPW.OW = AIPW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    weight.type = "OW",
    estimand.type = "SPCE",
    evaluate.time = truncate,
    dependent.adjustment = dependent.censoring
  )
  if (multi.arm) {
    ow_AIPW_est[i, 7:9] <- res.AIPW.OW$tau
  } else{
    ow_AIPW_est[i, 7:8] <- res.AIPW.OW$tau
  }
  
  # OW-MAO, IPW-MAO
  if (mao.method) {
    if (!multi.arm) {
      tmax = min(max(Y[Z == 1]), max(Y[Z == 0]))
      ## IPW based
      res.mao.IPW = estimand_analysis(
        X = X[,-1],
        Z = Z,
        Y = Y,
        delta = DELTA,
        weight.type = "IPW",
        t.trunc = truncate,
        tmax = tmax
      )
      ipw_est_mao[i, c(1, 4, 7)] = res.mao.IPW$ans[, 2]
      ipw_se_mao[i, c(1, 4, 7)] = res.mao.IPW$ans[, 3]
      ## OW based
      res.mao.OW = estimand_analysis(
        X = X[,-1],
        Z = Z,
        Y = Y,
        delta = DELTA,
        weight.type = "OVERLAP",
        t.trunc = truncate,
        tmax = tmax
      )
      ow_est_mao[i, c(1, 4, 7)] = res.mao.OW$ans[, 2]
      ow_se_mao[i, c(1, 4, 7)] = res.mao.OW$ans[, 3]
    } else{
      X.sub = X[(Z == 0 | Z == 1), ]
      Z.sub = Z[Z == 0 | Z == 1]
      Z.sub = as.numeric(Z.sub == 1)
      Y.sub = Y[Z == 0 | Z == 1]
      DELTA.sub = DELTA[Z == 0 | Z == 1]
      tmax = min(max(Y[Z == 1]), max(Y[Z == 0]))
      ## IPW based
      res.mao.IPW = estimand_analysis(
        X = X.sub[,-1],
        Z = Z.sub,
        Y = Y.sub,
        delta = DELTA.sub,
        weight.type = "IPW",
        t.trunc = truncate,
        tmax = tmax
      )
      ipw_est_mao[i, c(1, 4, 7)] = res.mao.IPW$ans[, 2]
      ipw_se_mao[i, c(1, 4, 7)] = res.mao.IPW$ans[, 3]
      ## OW based
      res.mao.OW = estimand_analysis(
        X = X.sub[,-1],
        Z = Z.sub,
        Y = Y.sub,
        delta = DELTA.sub,
        weight.type = "OVERLAP",
        t.trunc = truncate,
        tmax = tmax
      )
      ow_est_mao[i, c(1, 4, 7)] = res.mao.OW$ans[, 2]
      ow_se_mao[i, c(1, 4, 7)] = res.mao.OW$ans[, 3]
      
      X.sub = X[(Z == 1 | Z == 2), ]
      Z.sub = Z[Z == 1 | Z == 2]
      Z.sub = as.numeric(Z.sub == 2)
      Y.sub = Y[Z == 1 | Z == 2]
      DELTA.sub = DELTA[Z == 1 | Z == 2]
      tmax = min(max(Y[Z == 1]), max(Y[Z == 2]))
      ## IPW based
      res.mao.IPW = estimand_analysis(
        X = X.sub[,-1],
        Z = Z.sub,
        Y = Y.sub,
        delta = DELTA.sub,
        weight.type = "IPW",
        t.trunc = truncate,
        tmax = tmax
      )
      ipw_est_mao[i, c(3, 6, 9)] = res.mao.IPW$ans[, 2]
      ipw_se_mao[i, c(3, 6, 9)] = res.mao.IPW$ans[, 3]
      ## OW based
      res.mao.OW = estimand_analysis(
        X = X.sub[,-1],
        Z = Z.sub,
        Y = Y.sub,
        delta = DELTA.sub,
        weight.type = "OVERLAP",
        t.trunc = truncate,
        tmax = tmax
      )
      ow_est_mao[i, c(3, 6, 9)] = res.mao.OW$ans[, 2]
      ow_se_mao[i, c(3, 6, 9)] = res.mao.OW$ans[, 3]
    }
  }
  
  ### Cox-G
  if (cox.q.method) {
    res.cox.q = cox.q.model.fit(Y,
                                DELTA,
                                X,
                                Z,
                                truncate =  truncate,
                                # ipw.weighted = T,
                                boot.time = 250)
    cox_q_est[i, ] = res.cox.q$est
    cox_q_se[i, ] = res.cox.q$se
    
    
    # WEIGHTED BY IPW
    res.cox.q = cox.q.model.fit(
      Y,
      DELTA,
      X,
      Z,
      truncate =  truncate,
      ipw.weighted = T,
      boot.time = 250
    )
    cox_ipw_est[i, ] = res.cox.q$est
    cox_ipw_se[i, ] = res.cox.q$se
  }
  
  ## Cox-IPW
  if (cox.msm.method) {
    res.cox.msm = cox.msm.model.fit(Y,
                                    DELTA,
                                    X,
                                    Z,
                                    truncate =  truncate,
                                    boot.time = 250)
    cox_msm_est[i, ] = res.cox.msm$est
    cox_msm_se[i, ] = res.cox.msm$se
  }
  
  print(paste("==", i, "=="))
  
  setwd("../simulation_results/")
  save(
    i,
    ow_est,
    ipw_est,
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
  setwd("../simulation/")
}

print(paste("==Sample size", sample_size, "=="))
print(paste("==overlap", good_overlap, "=="))
