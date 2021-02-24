library(mice)
library(survival)
source("../estimators/fast_pseudo_calculation.R")

G.pseudo <- function(Y,
                     Z,
                     DELTA,
                     X,
                     estimand.type = "ASCE",
                     evaluate.time = NA,
                     by.group.version = F,
                     dependent.adjustment = F,
                     if.cloglog.link = F) {
  ############
  # Y: (vec) Vector of end time min(T,C)
  # Z: (vec) Vector of treatment arm, encoding from 0 to (J-1)
  # DELTA: (vec) Whether dead within the period
  # X: (matrix N*p) design matrix
  # estimand.type: (str) type of survival estimand SPCE, ASCE, RACE
  # evaluate.time: (float) time point for evaluation (SPCE, RACE)
  # by.group.version: (bool) if we calculate PO in different arms.
  # dependent.censoring: (bool): whether use the model counted for dependent censoring.
  # if.cloglog.link (bool): whether we use cloglog link
  ############
  
  ### Add intercept if not for the design matrix
  if (!all(X[, 1] == 1)) {
    X = cbind(1, X)
  }
  
  ### Summary statistics
  N = length(Y)
  
  if (estimand.type == "ASCE") {
    evaluate.time = (max(Y) + max(Y[Y != max(Y)])) / 2
  }
  
  ### Handle error for truncation
  if (max(Y) < evaluate.time) {
    return (list(tau = NA, se = NA))
  }
  if (dependent.adjustment) {
    Censoring.DELTA = 1 - DELTA
    censoring.model = coxph(Surv(Y, Censoring.DELTA) ~ X)
    km.object = survfit(censoring.model, newdata = data.frame(X = X))
    km.time = km.object$time
    ind = unlist(lapply(
      pmin(Y, evaluate.time),
      FUN = function(x) {
        max(which(x >= km.time))
      }
    ))
    weight.prob = unlist(lapply(
      1:N,
      FUN = function(x) {
        km.object$surv[ind[x], x]
      }
    ))
    
    if (estimand.type == "ASCE") {
      V = Y
    } else if (estimand.type == "RACE") {
      V = pmin(Y, evaluate.time)
    } else if (estimand.type == "SPCE") {
      V = as.numeric(Y >= evaluate.time)
    } else{
      stop("Undefined estimand")
    }
    # truncation
    keep.ind = which(weight.prob >= 0.03)
    V = V[keep.ind]
    weight.prob = weight.prob[keep.ind]
    Y = Y[keep.ind]
    Z = Z[keep.ind]
    X = X[keep.ind,]
    DELTA = DELTA[keep.ind]
    N = length(Y)
    indicator = as.numeric((DELTA == 1) | (Y >= evaluate.time))
    pseudo.obs = V * indicator / weight.prob
    var.method = 3
  } else{
    ### Pseudo Obs
    if (estimand.type == "ASCE") {
      if (by.group.version) {
        pseudo.obs = rep(NA, length(Y))
        for (k in 0:(max(Z))) {
          pseudo.obs[Z == k] = fast.pseudo(Y[Z == k], event = DELTA[Z == k], type = "mean")
        }
      } else{
        pseudo.obs = fast.pseudo(Y, event = DELTA, type = "mean")
      }
    } else if (estimand.type == "RACE") {
      if (by.group.version) {
        for (k in 0:max(Z)) {
          if (max(Y[Z == k][DELTA[Z == k] == 1]) < evaluate.time) {
            print("not valid for by group pseudo obs calculation")
            return (list(
              tau = NA,
              mu = NA,
              se = NA,
              alpha.matrix = NA,
              Sigma.matrix = NA
            ))
          }
        }
        pseudo.obs = rep(NA, length(Y))
        for (k in 0:(max(Z))) {
          pseudo.obs[Z == k] = fast.pseudo(Y[Z == k],
                                           event = DELTA[Z == k],
                                           tmax = evaluate.time,
                                           type = "mean")
        }
      } else{
        pseudo.obs = fast.pseudo(Y,
                                 event = DELTA,
                                 tmax = evaluate.time,
                                 type = "mean")
      }
      
    } else if (estimand.type == "SPCE") {
      if (by.group.version) {
        pseudo.obs = rep(NA, length(Y))
        for (k in 0:(max(Z))) {
          pseudo.obs[Z == k] = fast.pseudo(Y[Z == k],
                                           event = DELTA[Z == k],
                                           tmax = evaluate.time,
                                           type = "surv")
        }
      } else{
        pseudo.obs = fast.pseudo(Y,
                                 event = DELTA,
                                 tmax = evaluate.time,
                                 type = "surv")
      }
    } else{
      stop("Undefined estimand")
    }
  }
  
  J = max(Z) + 1   # Number of arms
  A = matrix(unlist(lapply(
    Z,
    FUN = function(x) {
      as.numeric(x == 0:(J - 1))
    }
  )), byrow = T, ncol = J)  # A_ij matrix
  # muhat.h = as.numeric(colSums(A * pseudo.obs * w) / colSums(A * w))  # \mu^{h}_{j}
  
  ind.alpha = combn(J, 2)
  alpha.matrix = matrix(0, nrow = ncol(ind.alpha), ncol = J)
  alpha.matrix[cbind(1:nrow(alpha.matrix), ind.alpha[1, ])] = -1
  alpha.matrix[cbind(1:nrow(alpha.matrix), ind.alpha[2, ])] = 1
  
  
  LINK.TEXT = "log"
  if (estimand.type == "SPCE") {
    # maybe subject to numeric error if use cloglog
    if (if.cloglog.link) {
      LINK.TEXT = "cloglog"
    } else{
      LINK.TEXT = "identity"
    }
  }
  data.temp = cbind(X[,-1], Z)
  colnames(data.temp) = c(paste("X", 1:(ncol(X) - 1), sep = "."), "Z")
  data.temp = data.frame(data.temp)
  data.temp$Z = as.factor(data.temp$Z)
  formula.text = paste("pseudo.obs~", paste(
    paste(paste("X", 1:(ncol(
      X
    ) - 1), sep = "."), collapse = "+"),
    "Z",
    paste(paste(paste(
      "X", 1:(ncol(X) - 1), sep = "."
    ), "*Z"), collapse = "+"),
    sep = "+"
  ))
  outcome.model  = glm(as.formula(formula.text),
                       data = data.temp,
                       family = gaussian(link = LINK.TEXT))
  mu.fit = matrix(NA, nrow = length(Z), ncol = J)
  X.m = NULL
  Mean.Matrix = NULL
  for (j in 0:(J - 1)) {
    data.use = data.temp
    data.use$Z = rep(factor(j, levels = 0:(J - 1)))
    mu.fit[, (j + 1)] = predict(outcome.model, newdata = data.use, type =
                                  "response")
    # create some design matrix
    X.m = rbind(X.m, model.matrix.default(as.formula(formula.text), data.use))
    zero = matrix(0, nrow = N, ncol = J)
    if (LINK.TEXT == "log") {
      zero[, (j + 1)] = mu.fit[, (j + 1)]
    } else if (LINK.TEXT == "identity") {
      zero[, (j + 1)] = 1
    } else{
      zero[, (j + 1)] = exp(mu.fit[, (j + 1)] - exp(mu.fit[, (j + 1)]))
    }
    Mean.Matrix = rbind(Mean.Matrix, as.numeric(zero) / N)
  }
  mu.model = as.numeric(colMeans(mu.fit))
  muhat.h = mu.model
  ind.alpha = combn(J, 2)
  alpha.matrix = matrix(0, nrow = ncol(ind.alpha), ncol = J)
  alpha.matrix[cbind(1:nrow(alpha.matrix), ind.alpha[1, ])] = -1
  alpha.matrix[cbind(1:nrow(alpha.matrix), ind.alpha[2, ])] = 1
  tau.h = as.matrix(muhat.h %*% t(alpha.matrix))  # Point estimation
  
  # Empirical variance
  COV = vcov(outcome.model, method = "white1", type = "HC1")
  # Covariance for mu1,mu0..
  COV.m = Mean.Matrix %*% X.m %*% COV %*% t(X.m) %*% t(Mean.Matrix)
  # variance for tau12,tau23
  Sigmah = diag(alpha.matrix %*% COV.m %*% t(alpha.matrix))
  return (
    list(
      tau = tau.h,
      mu = muhat.h,
      se = sqrt(Sigmah),
      alpha.matrix = alpha.matrix,
      Sigma.matrix = COV.m
    )
  )
}


unadj.pseudo <- function(Y,
                         Z,
                         DELTA,
                         X,
                         estimand.type = "ASCE",
                         evaluate.time = NA,
                         by.group.version = F,
                         dependent.adjustment = F) {
  ############
  # Y: (vec) Vector of end time min(T,C)
  # Z: (vec) Vector of treatment arm, encoding from 0 to (J-1)
  # DELTA: (vec) Whether dead within the period
  # X: (matrix N*p) design matrix
  # estimand.type: (str) type of survival estimand SPCE, ASCE, RACE
  # evaluate.time: (float) time point for evaluation (SPCE, RACE)
  # by.group.version: (bool) if we calculate PO in different arms.
  # dependent.censoring: (bool): whether use the model counted for dependent censoring.
  ############
  
  ### Add intercept if not for the design matrix
  if (!all(X[, 1] == 1)) {
    X = cbind(1, X)
  }
  
  ### Summary statistics
  N = length(Y)
  
  if (estimand.type == "ASCE") {
    evaluate.time = (max(Y) + max(Y[Y != max(Y)])) / 2
  }
  
  ### Handle error for truncation
  if (max(Y) < evaluate.time) {
    return (list(tau = NA, se = NA))
  }
  if (dependent.adjustment) {
    Censoring.DELTA = 1 - DELTA
    censoring.model = coxph(Surv(Y, Censoring.DELTA) ~ X)
    km.object = survfit(censoring.model, newdata = data.frame(X = X))
    km.time = km.object$time
    ind = unlist(lapply(
      pmin(Y, evaluate.time),
      FUN = function(x) {
        max(which(x >= km.time))
      }
    ))
    weight.prob = unlist(lapply(
      1:N,
      FUN = function(x) {
        km.object$surv[ind[x], x]
      }
    ))
    
    if (estimand.type == "ASCE") {
      V = Y
    } else if (estimand.type == "RACE") {
      V = pmin(Y, evaluate.time)
    } else if (estimand.type == "SPCE") {
      V = as.numeric(Y >= evaluate.time)
    } else{
      stop("Undefined estimand")
    }
    # truncation
    keep.ind = which(weight.prob >= 0.03)
    V = V[keep.ind]
    weight.prob = weight.prob[keep.ind]
    Y = Y[keep.ind]
    Z = Z[keep.ind]
    X = X[keep.ind,]
    DELTA = DELTA[keep.ind]
    N = length(Y)
    indicator = as.numeric((DELTA == 1) | (Y >= evaluate.time))
    pseudo.obs = V * indicator / weight.prob
    var.method = 3
  } else{
    ### Pseudo Obs
    if (estimand.type == "ASCE") {
      if (by.group.version) {
        pseudo.obs = rep(NA, length(Y))
        for (k in 0:(max(Z))) {
          pseudo.obs[Z == k] = fast.pseudo(Y[Z == k], event = DELTA[Z == k], type = "mean")
        }
      } else{
        pseudo.obs = fast.pseudo(Y, event = DELTA, type = "mean")
      }
    } else if (estimand.type == "RACE") {
      if (by.group.version) {
        for (k in 0:max(Z)) {
          if (max(Y[Z == k][DELTA[Z == k] == 1]) < evaluate.time) {
            print("not valid for by group pseudo obs calculation")
            return (list(
              tau = NA,
              mu = NA,
              se = NA,
              alpha.matrix = NA,
              Sigma.matrix = NA
            ))
          }
        }
        pseudo.obs = rep(NA, length(Y))
        for (k in 0:(max(Z))) {
          pseudo.obs[Z == k] = fast.pseudo(Y[Z == k],
                                           event = DELTA[Z == k],
                                           tmax = evaluate.time,
                                           type = "mean")
        }
      } else{
        pseudo.obs = fast.pseudo(Y,
                                 event = DELTA,
                                 tmax = evaluate.time,
                                 type = "mean")
      }
      
    } else if (estimand.type == "SPCE") {
      if (by.group.version) {
        pseudo.obs = rep(NA, length(Y))
        for (k in 0:(max(Z))) {
          pseudo.obs[Z == k] = fast.pseudo(Y[Z == k],
                                           event = DELTA[Z == k],
                                           tmax = evaluate.time,
                                           type = "surv")
        }
      } else{
        pseudo.obs = fast.pseudo(Y,
                                 event = DELTA,
                                 tmax = evaluate.time,
                                 type = "surv")
      }
    } else{
      stop("Undefined estimand")
    }
  }
  
  J = max(Z) + 1   # Number of arms
  A = matrix(unlist(lapply(
    Z,
    FUN = function(x) {
      as.numeric(x == 0:(J - 1))
    }
  )), byrow = T, ncol = J)  # A_ij matrix
  muhat.h = as.numeric(colSums(A * pseudo.obs) / colSums(A))
  
  ind.alpha = combn(J, 2)
  alpha.matrix = matrix(0, nrow = ncol(ind.alpha), ncol = J)
  alpha.matrix[cbind(1:nrow(alpha.matrix), ind.alpha[1, ])] = -1
  alpha.matrix[cbind(1:nrow(alpha.matrix), ind.alpha[2, ])] = 1
  
  var.list = unlist(lapply(
    0:(J - 1),
    FUN = function(x) {
      var(pseudo.obs[Z == x]) / sum(Z == x)
    }
  ))
  
  tau.h = as.matrix(muhat.h %*% t(alpha.matrix))  # Point estimation
  
  COV.m = diag(var.list)
  # variance for tau12,tau23
  Sigmah = diag(alpha.matrix %*% COV.m %*% t(alpha.matrix))
  return (
    list(
      tau = tau.h,
      mu = muhat.h,
      se = sqrt(Sigmah),
      alpha.matrix = alpha.matrix,
      Sigma.matrix = COV.m
    )
  )
}
