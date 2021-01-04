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
  hazard = as.vector(gamma.1 * as.numeric(Z == 1) + gamma.2 * as.numeric(Z ==
                                                                           2) + X[, -1] %*% alpha)
  random_simu = runif(sample_size)
  Survival_time = (-log(random_simu) / (lambda * exp(hazard))) ^ (1 / v)
  
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
  # AFT but not PH; log-normal
  gamma.1 = -1
  gamma.2 = -1
  alpha = -c(2, 1.5, -1, 1)
  lambda = 0.0001
  v = 3
  truncate = truncate.aft
  random_simu = rnorm(sample_size, sd = 0.81)
  Survival_time = exp(3.5 + 0.2 * (
    as.vector(
      gamma.1 * as.numeric(Z == 1) + gamma.2 * as.numeric(Z == 2) + X[, -1] %*%
        alpha
    )
  ) + random_simu)
  Survival_time_0 = exp(3.5 + 0.2 * (X[, -1] %*% alpha) + random_simu)
  Survival_time_1 = exp(3.5 + 0.2 * (gamma.1 + X[, -1] %*% alpha) + random_simu)
  Survival_time_2 = exp(3.5 + 0.2 * (gamma.2 + X[, -1] %*% alpha) + random_simu)
}
# Censor time
Censor_time = runif(sample_size, min = 0, max = 115)
# Censor_time=runif(sample_size,min=2000,max=2115)

# Observed survival
Censored = Survival_time > Censor_time
Y = Survival_time
Y[Censored] = Censor_time[Censored]
DELTA = 1 - Censored

# Store true value on finite sample
# truncate=min(quantile(Survival_time_1[Z==1],0.75),quantile(Survival_time_0[Z==0],0.75))
tilt.h = 1 / rowSums(1 / soft.max)
if (!multi.arm) {
  true_est_finite[i, 1] = mean(Survival_time_1) - mean(Survival_time_0)
  true_est_ow_finite[i, 1] = sum((Survival_time_1 - Survival_time_0) * tilt.h) /
    sum(tilt.h)
  
  true_est_finite[i, 3] = mean(Survival_time_1 > truncate) - mean(Survival_time_0 >
                                                                    truncate)
  true_est_ow_finite[i, 3] = sum((
    as.numeric(Survival_time_1 > truncate) - as.numeric(Survival_time_0 > truncate)
  ) * tilt.h) / sum(tilt.h)
  
  Survival_time_1[Survival_time_1 > truncate] = truncate
  Survival_time_0[Survival_time_0 > truncate] = truncate
  true_est_finite[i, 2] = mean(Survival_time_1) - mean(Survival_time_0)
  true_est_ow_finite[i, 2] = sum((Survival_time_1 - Survival_time_0) * tilt.h) /
    sum(tilt.h)
} else{
  # alpha = c(-1,0,1)
  true_est_finite[i, 1] = mean(Survival_time_2) - mean(Survival_time_0)
  true_est_ow_finite[i, 1] = sum((Survival_time_2 - Survival_time_0) * tilt.h) /
    sum(tilt.h)
  
  true_est_finite[i, 3] = mean(Survival_time_2 > truncate) - mean(Survival_time_0 >
                                                                    truncate)
  true_est_ow_finite[i, 3] = sum((
    as.numeric(Survival_time_2 > truncate) - as.numeric(Survival_time_0 > truncate)
  ) * tilt.h) / sum(tilt.h)
  
  Survival_time_2[Survival_time_2 > truncate] = truncate
  Survival_time_0[Survival_time_0 > truncate] = truncate
  true_est_finite[i, 2] = mean(Survival_time_2) - mean(Survival_time_0)
  true_est_ow_finite[i, 2] = sum((Survival_time_2 - Survival_time_0) * tilt.h) /
    sum(tilt.h)
}
