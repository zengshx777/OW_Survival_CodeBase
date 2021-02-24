#
library(nnet)
library(mice)
library(numDeriv)
library(survival)
source("../estimators/fast_pseudo_calculation.R")

AIPW.pseudo <- function(Y,
                        Z,
                        DELTA,
                        X,
                        weight.type = "IPW",
                        estimand.type = "ASCE",
                        evaluate.time = NA,
                        var.method = 2,
                        target.j = NA,
                        ps.threshold = NA,
                        by.group.version = F,
                        dependent.adjustment = F) {
  ############
  # Y: (vec) Vector of end time min(T,C)
  # Z: (vec) Vector of treatment arm, encoding from 0 to (J-1)
  # DELTA: (vec) Whether dead within the period
  # X: (matrix N*p) design matrix
  # weight.type: (str) weighting schemes
  # estimand.type: (str) type of survival estimand SPCE, ASCE, RACE
  # evaluate.time: (float) time point for evaluation (SPCE, RACE)
  # var.method: (int) type of variance calculation
  # target.j: (int) index for arm, ranging from 0 to (J-1), only used for ATT
  # ps.threshold: threshold for ps truncation 0.1 e.g., usually for IPW
  # by.group.version: (bool) if we calculate PO in different arms.
  # dependent.censoring: (bool): whether use the model counted for dependent censoring.
  ############
  
  ### Add intercept if not for the design matrix
  if (!all(X[, 1] == 1)) {
    X = cbind(1, X)
  }
  
  ### Summary statistics
  N = length(Y)
  
  ### Estimate ps
  # fit <- glm(Z ~ X, family = binomial(link = "logit"))
  # e.h <- as.numeric(fit$fitted.values)
  gps.fit = multinom(Z ~ -1 + X,
                     maxit = 500,
                     Hess = TRUE,
                     trace = FALSE)
  Thetah = t(coef(gps.fit))   # matrix coefficient
  thetah = c(Thetah)          # vec operator - long vector
  IthetaInv = N * vcov(gps.fit) # Information Matrix
  e = gps.fit$fitted.values   # Predicted PS
  ### Adjust for the case with two arms
  if (length(unique(Z)) == 2) {
    e = cbind(1 - e, e)
  }
  
  
  ### Truncation for poor overlapping units
  if (!is.na(ps.threshold))
  {
    if (weight.type == "OW") {
      ### Generally not used for OW
      warning("truncate propensity score for overlap weights")
    }
    keep.ind = which(apply(
      e,
      1,
      FUN = function(x) {
        (min(x) > ps.threshold) && (max(x) < (1 - ps.threshold))
      }
    ))
    keep.ratio = length(keep.ind) / N
    if (keep.ratio < 0.8)
    {
      warning(paste((1 - keep.ratio) * 100,
                    "%(>=20%) observations are disgarded due to PS truncation",
                    sep = ""
      ))
    }
    Y = Y[keep.ind]
    e = e[keep.ind,]
    Z = Z[keep.ind]
    X = X[keep.ind,]
    DELTA = DELTA[keep.ind]
    ### Summary statistics recalculate
    N = length(Y)
  }
  
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
    e = e[keep.ind,]
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
  eInv = 1 / e   # Inverse of PS
  # print("calculate pseudo finished")
  ### Weighting Type
  if (weight.type == "IPW") {
    tilt.h = rep(1, length(Z))
  } else if (weight.type == "OW") {
    tilt.h = as.numeric(1 / rowSums(eInv))
  } else if (weight.type == "ATT") {
    tilt.h = e[, (target.j + 1)]
  } else{
    stop("Undefined weight type")
  }
  
  
  w = eInv * tilt.h # weight matrix w_ij
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
  residual.matrix = matrix(NA, nrow = length(Z), ncol = J)
  for (j in 0:(J - 1)) {
    data.use = data.temp
    data.use$Z = rep(factor(j, levels = 0:(J - 1)))
    mu.fit[, (j + 1)] = predict(outcome.model, newdata = data.use, type =
                                  "response")
    residual.matrix[, (j + 1)] = pseudo.obs - mu.fit[, (j + 1)]
  }
  mu.weight = as.numeric(colSums(A * residual.matrix * w) / colSums(A * w))
  mu.model = as.numeric(colSums(mu.fit * tilt.h) / sum(tilt.h))
  muhat.h = mu.model + mu.weight
  ind.alpha = combn(J, 2)
  alpha.matrix = matrix(0, nrow = ncol(ind.alpha), ncol = J)
  alpha.matrix[cbind(1:nrow(alpha.matrix), ind.alpha[1, ])] = -1
  alpha.matrix[cbind(1:nrow(alpha.matrix), ind.alpha[2, ])] = 1
  tau.h = as.matrix(muhat.h %*% t(alpha.matrix))  # Point estimation
  return(list(tau = tau.h, mu = muhat.h))
  
  # ### Calculate influence function for pseudo observations
  # if (var.method == 3) {
  #   obs.approxi = pseudo.obs
  # } else{
  #   if (by.group.version) {
  #     obs.approxi = rep(NA, length(Y))
  #     for (k in 0:max(Z)) {
  #       source("integral_utils.R")
  #       obs.approxi[Z == k] = calculate.obs.approximate (Y[Z == k], DELTA[Z ==
  #                                                                           k], evaluate.time, estimand.type, var.method)
  #     }
  #   } else{
  #     source("integral_utils.R")
  #     obs.approxi = calculate.obs.approximate (Y, DELTA, evaluate.time, estimand.type, var.method)
  #   }
  # }
  #
  #
  # ### Calculate the uncertainty from PS model
  # ### Calculate gradient of weights
  # wj = function(theta) {
  #   Theta = matrix(theta, ncol(X), ncol(A) - 1)
  #   Eta = X %*% Theta
  #   if (weight.type == "OW") {
  #     return(1 / (1 + as.numeric(rowSums(exp(
  #       -Eta
  #     )))))
  #   } else if (weight.type == "IPW") {
  #     return(1 + as.numeric(rowSums(exp(Eta))))
  #   } else if (weight.type == "ATT") {
  #     if (target.j >= 1) {
  #       return(as.numeric(exp(Eta[, target.j])))
  #     } else{
  #       return (1)
  #     }
  #   } else{
  #     stop("Undefined weight type")
  #   }
  # }
  # wdotj = jacobian(wj, thetah)
  # Hmat = c(colMeans(A[, 1] * (pseudo.obs - muhat.h[1]) * wdotj))
  #
  # for (j in 2:J) {
  #   wj = function(theta) {
  #     Theta = matrix(theta, ncol(X), ncol(A) - 1)
  #     Eta = X %*% Theta
  #     if (weight.type == "OW") {
  #       return(as.numeric(exp(-Eta[, j - 1])) / (1 + as.numeric(rowSums(exp(
  #         -Eta
  #       )))))
  #     } else if (weight.type == "IPW") {
  #       return((1 + as.numeric(rowSums(exp(
  #         Eta
  #       )))) / exp(Eta[, j - 1]))
  #     } else if (weight.type == "ATT") {
  #       if (target.j >= 1) {
  #         return (exp(Eta[, target.j] - Eta[, j - 1]))
  #       } else{
  #         return (exp(-Eta[, j - 1]))
  #       }
  #     } else{
  #       stop("Undefined weight type")
  #     }
  #   }
  #   wdotj = jacobian(wj, thetah)
  #   Hmat = rbind(Hmat, c(colMeans(A[, j] * (
  #     pseudo.obs - muhat.h[j]
  #   ) * wdotj)))
  # }
  #
  # ### Score function
  # loglik = function(theta) {
  #   Theta = matrix(theta, ncol(X), ncol(A) - 1)
  #   Eta = X %*% Theta
  #   ltheta = as.numeric(rowSums(A[,-1] * Eta) - log(1 + rowSums(exp(Eta))))
  #   return(ltheta)
  # }
  # Sthetah = jacobian(loglik, thetah)
  #
  #
  # ### Calculate influence function
  # yMat = matrix(rep(obs.approxi, J), N, J)
  # mhatMat = matrix(rep(muhat.h, each = N), N, J)
  # Psi = t(A * (yMat - mhatMat) * w) + Hmat %*% IthetaInv %*% t(Sthetah)
  # omega = mean(tilt.h)
  #
  # Sigmah = diag(alpha.matrix %*% tcrossprod(Psi) %*% t(alpha.matrix) / (N * omega) ^ 2)
  # Sigma.matrix = tcrossprod(Psi) / (N * omega) ^ 2
  #
  # return (
  #   list(
  #     tau = tau.h,
  #     mu = muhat.h,
  #     se = sqrt(Sigmah),
  #     alpha.matrix = alpha.matrix,
  #     Sigma.matrix = Sigma.matrix
  #   )
  # )
}