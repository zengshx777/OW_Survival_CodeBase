library(pseudo)
library(nnet)
library(mice)
# library(rootSolve)
library(numDeriv)
library(survival)
source("fast_pseudo_calculation.R")

PSW.pseudo <- function(Y,
                       Z,
                       DELTA,
                       X,
                       alpha = c(-1, 1),
                       weight.type = "IPW",
                       estimand.type = "ASCE",
                       evaluate.time = NA,
                       var.method = 1,
                       target.j = NA,
                       ps.threshold = NA) {
  ############
  # Y: (vec) Vector of end time min(T,C)
  # Z: (vec) Vector of treatment arm, encoding from 0 to (J-1)
  # DELTA: (vec) Whether dead within the period
  # X: (matrix N*p) design matrix
  # alpha: (vec) setup for estimand
  # weight.type: (str) weighting schemes
  # estimand.type: (str) type of survival estimand SPCE, ASCE, RACE
  # evaluate.time: (float) time point for evaluation (SPCE, RACE)
  # var.method: (int) type of variance calculation
  # target.j: (int) index for arm, ranging from 0 to (J-1), only used for ATT
  # ps.threshold: threshold for ps truncation 0.1 e.g., usually for IPW
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
  
  eInv = 1 / e   # Inverse of PS
  
  if (estimand.type == "ASCE") {
    evaluate.time = (max(Y) + max(Y[Y != max(Y)])) / 2
  }
  
  ### Handle error for truncation
  if (max(Y) < evaluate.time) {
    return (list(tau = NA, se = NA))
  }
  
  ### Pseudo Obs
  if (estimand.type == "ASCE") {
    pseudo.obs = fast.pseudo(Y, event = DELTA, type = "mean")
  } else if (estimand.type == "RACE") {
    pseudo.obs = fast.pseudo(Y, event = DELTA, tmax = evaluate.time, type = "mean")
  } else if (estimand.type == "SPCE") {
    pseudo.obs = fast.pseudo(Y, event = DELTA, tmax = evaluate.time, type = "surv")
  } else{
    stop("Undefined estimand")
  }
  print("calculate pseudo finished")
  ### Weighting Type
  if (weight.type == "IPW") {
    tilt.h = 1
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
  muhat.h = as.numeric(colSums(A * pseudo.obs * w) / colSums(A * w))  # \mu^{h}_{j}
  tau.h = as.matrix(alpha %*% muhat.h)  # Point estimation
  
  ### Calculate Cumulative Hazard
  temp.data = data.frame(Y = Y, DELTA = DELTA)
  Cumu.hazard = nelsonaalen(temp.data, Y, DELTA)
  CH.list = list(Time = sort(Y), CH = sort(Cumu.hazard))
  
  # Jump for cumulative hazard function
  first.diff.CH = c(0, CH.list$CH[1], diff(CH.list$CH))
  eval.time = c(0, CH.list$Time)
  
  H.left = ecdf(Y) # Calculate the Probability of (not) being at risk
  
  KM.object = survfit(Surv(Y, DELTA) ~ 1, data = temp.data)
  S = summary(KM.object)   # KM Estimate of S
  
  time.grid = c(0, S$time) # Time grids for evaluation
  diff.time.grid = c(diff(time.grid), 0) # Difference of time
  final.index = max(which(time.grid <= evaluate.time))
  diff.time.grid.temp = diff.time.grid
  diff.time.grid.temp[final.index] = evaluate.time - time.grid[final.index]
  
  ### source("FCLT_util.R")
  ### Wrap surv prob into a function
  Surv.fun <- function(time.input)
  {
    if (time.input < min(S$time)) {
      return (1)
    }
    return(S$surv[max(which(S$time <= time.input))])
  }
  
  # Utility function to calculate the integral from 0 to tmax
  # first order \int_{0}^{tmax}1/H(s)d CH(s)
  # second order \int_{0}^{tmax}1/H^2(s)d CH(s)
  integral.H.CH <- function(tmax, order = 1) {
    end.index = max(which(eval.time <= tmax))
    if (order == 1) {
      jumps = unlist(lapply(
        1:end.index,
        FUN = function(i) {
          first.diff.CH[i] / (1 - H.left(eval.time[i]))
        }
      ))
    } else{
      jumps = unlist(lapply(
        1:end.index,
        FUN = function(i) {
          first.diff.CH[i] / ((1 - H.left(eval.time[i])) ^ 2)
        }
      ))
    }
    return (sum(jumps))
  }
  
  ### Tau'(Xi,F) for cumulative hazard functional
  Tau.1 <- function(i, end.time) {
    Tau.1.order = -integral.H.CH(tmax = min(end.time, Y[i]), order = 1)
    if (DELTA[i] == 1 && Y[i] <= end.time)
    {
      Tau.1.order = Tau.1.order + 1 / (1 - H.left(Y[i]))
    }
    return (Tau.1.order)
  }
  
  ### Tau''(Xj,Xi,F) for cumulative hazard functional
  Tau.2 <- function(i, j, end.time) {
    Tau.2.order =   Tau.1(i, end.time) + Tau.1(j, end.time) +
      2 * integral.H.CH(tmax = min(end.time, Y[i], Y[j]), order = 2)
    if (DELTA[i] == 1 && Y[i] <= min(Y[j], end.time))
    {
      Tau.2.order = Tau.2.order - 1 / ((1 - H.left(Y[i])) ^ 2)
    }
    if (DELTA[j] == 1 && Y[j] <= min(Y[i], end.time))
    {
      Tau.2.order = Tau.2.order - 1 / ((1 - H.left(Y[j])) ^ 2)
    }
    return (Tau.2.order)
  }
  
  ### Tau''(Xj,Xi,F)
  # Tau.2.matrix<-function(end.time){
  #
  #   second.order.matrix=matrix(NA,N,N)
  #   for (i in 1:N){
  #     for (j in 1:i){
  #       second.order.matrix[i,j] =  Tau.1(i, end.time)+Tau.1(j, end.time)+
  #         2*integral.H.CH(tmax = min(end.time,Y[i],Y[j]), order = 2)
  #     }
  #     if(DELTA[i]==1&&Y[i]<=min(Y[j],end.time))
  #     {
  #       Tau.2.order = Tau.2.order - 1/((1-H.left(Y[i]))^2)
  #     }
  #     if(DELTA[j]==1&&Y[j]<=min(Y[i],end.time))
  #     {
  #       Tau.2.order = Tau.2.order - 1/((1-H.left(Y[j]))^2)
  #     }
  #   }
  #   return (second.order.matrix)
  # }
  #
  
  ### T.1 for Surv Prob
  T.1.order <- function(i, time.point)
  {
    return (-Surv.fun(time.point) * Tau.1(i, end.time = time.point))
  }
  
  ### T.2 for Surv Prob
  T.2.order <- function(i, j, time.point)
  {
    T.2.order.value = -Surv.fun(time.point) * (
      Tau.2(i, j, end.time = time.point) - Tau.1(i, end.time = time.point) * Tau.1(j, end.time = time.point)
    )
    if (i == j && DELTA[i] == 1 && Y[i] <= time.point)
    {
      T.2.order.value = T.2.order.value - Surv.fun(time.point) / ((1 - H.left(Y[i])) ^
                                                                    2)
    }
    return (T.2.order.value)
  }
  
  ### Functional value, zero order derivative
  zero.derivative <- function() {
    if (estimand.type == "SPCE") {
      zero.derivative.value = Surv.fun(evaluate.time)
    } else {
      zero.derivative.value = sum(unlist(lapply(
        1:final.index,
        FUN = function(x) {
          Surv.fun(time.grid[x]) * diff.time.grid.temp[x]
        }
      )))
    }
    return (zero.derivative.value)
  }
  
  ### Functional Derivative, first order
  first.derivative <- function(i)
  {
    if (estimand.type == "SPCE")
    {
      first.derivative.value = T.1.order(i, time.point = evaluate.time)
    } else {
      first.derivative.value = sum(unlist(lapply(
        1:final.index,
        FUN = function(x) {
          T.1.order(i, time.point = time.grid[x]) * diff.time.grid.temp[x]
        }
      )))
    }
    return (first.derivative.value)
  }
  
  ### Functional Derivative, second order
  second.derivative <- function(i, j)
  {
    if (estimand.type == "SPCE")
    {
      second.derivative.value = T.2.order(i, j, time.point = evaluate.time)
    } else {
      ### Calculate second order derivative through integral
      second.derivative.value = sum(unlist(lapply(
        1:final.index,
        FUN = function(x) {
          T.2.order(i, j, time.point = time.grid[x]) * diff.time.grid.temp[x]
        }
      )))
    }
    return (second.derivative.value)
  }
  
  ### Calculate influence function
  ### Calculate the uncertainty from PS model
  ### Calculate gradient of weights
  wj = function(theta) {
    Theta = matrix(theta, ncol(X), ncol(A) - 1)
    Eta = X %*% Theta
    if (weight.type == "OW") {
      return(1 / (1 + as.numeric(rowSums(exp(
        -Eta
      )))))
    } else if (weight.type == "IPW") {
      return(1 + as.numeric(rowSums(exp(Eta))))
    } else if (weight.type == "ATT") {
      if (target.j >= 1) {
        return(as.numeric(exp(Eta[, target.j])))
      } else{
        return (1)
      }
    } else{
      stop("Undefined weight type")
    }
  }
  wdotj = jacobian(wj, thetah)
  Hmat = c(colMeans(A[, 1] * (pseudo.obs - muhat.h[1]) * wdotj))
  
  for (j in 2:J) {
    wj = function(theta) {
      Theta = matrix(theta, ncol(X), ncol(A) - 1)
      Eta = X %*% Theta
      if (weight.type == "OW") {
        return(as.numeric(exp(-Eta[, j - 1])) / (1 + as.numeric(rowSums(exp(
          -Eta
        )))))
      } else if (weight.type == "IPW") {
        return((1 + as.numeric(rowSums(exp(
          Eta
        )))) / exp(Eta[, j - 1]))
      } else if (weight.type == "ATT") {
        if (target.j >= 1) {
          return (exp(Eta[, target.j] - Eta[, j - 1]))
        } else{
          return (exp(-Eta[, j - 1]))
        }
      } else{
        stop("Undefined weight type")
      }
    }
    wdotj = jacobian(wj, thetah)
    Hmat = rbind(Hmat, c(colMeans(A[, j] * (
      pseudo.obs - muhat.h[j]
    ) * wdotj)))
  }
  
  ### Score function
  loglik = function(theta) {
    Theta = matrix(theta, ncol(X), ncol(A) - 1)
    Eta = X %*% Theta
    ltheta = as.numeric(rowSums(A[,-1] * Eta) - log(1 + rowSums(exp(Eta))))
    return(ltheta)
  }
  Sthetah = jacobian(loglik, thetah)
  
  ### Construct observation for influence function
  if (var.method == 1) {
    ### Most accurate, taking account of second derivative
    adjust.matrix = as.vector((w * A) %*% (alpha))
    obs.approxi = zero.derivative() + unlist(lapply(1:N, first.derivative))
    ### Only applied for SPCE, computationnally expensive for ASCE, RACE
    ### TODO: Problem to be solved ###
    if (estimand.type == "SPCE")
    {
      obs.approxi = obs.approxi + unlist(lapply(
        1:N,
        FUN = function(j) {
          mean(unlist(lapply(
            1:N,
            FUN = function(i) {
              second.derivative(i, j) * adjust.matrix[i]
            }
          )))
        }
      ))
    }
  } else if (var.method == 2) {
    ### Less accurate, ignoring the second order derivative, overestimate
    obs.approxi = zero.derivative() + unlist(lapply(1:N, first.derivative))
  } else{
    ### Naive Calculation,ignoring the correlation between pseudos
    obs.approxi = pseudo.obs
  }
  
  ### Calculate influence function
  yMat = matrix(rep(obs.approxi, J), N, J)
  mhatMat = matrix(rep(muhat.h, each = N), N, J)
  Psi = t(A * (yMat - mhatMat) * w) + Hmat %*% IthetaInv %*% t(Sthetah)
  omega = mean(tilt.h)
  Sigmah = diag(t(alpha) %*% tcrossprod(Psi) %*% alpha / (N * omega) ^ 2)
  return (list(tau = tau.h, se = sqrt(Sigmah)))
}