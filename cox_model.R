######## Cox-q model and MSM model
library(survival)
library(nnet)

g.calculation <- function(km.object, end.time, type = 1) {
  time.grid = km.object$time
  surv.prob = km.object$surv
  n.size = ncol(surv.prob)
  time.grid = c(0, time.grid) # Time grids for evaluation
  diff.time.grid = c(diff(time.grid)) # Difference of time
  if (type == 1) {
    final.index = max(which(time.grid <= end.time))
    return (mean(surv.prob[final.index, ]))
  } else{
    if (!is.na(end.time)) {
      final.index = max(which(time.grid <= end.time))
      diff.time.grid.temp = diff.time.grid
      diff.time.grid.temp[final.index] = end.time - time.grid[final.index]
    } else{
      final.index = length(diff.time.grid)
      diff.time.grid.temp = diff.time.grid
    }
    output = mean(unlist(lapply(
      1:n.size,
      FUN = function(y) {
        sum(unlist(lapply(
          1:final.index,
          FUN = function(x) {
            surv.prob[x, y] * diff.time.grid.temp[x]
          }
        )))
      }
    )))
    return (output)
  }
}

cox.q.model.est <-
  function(Y, DELTA, X, Z, alpha = c(-1, 1), truncate) {
    J = max(Z) + 1   # Number of arms
    A = matrix(unlist(lapply(
      Z,
      FUN = function(x) {
        as.numeric(x == 0:(J - 1))
      }
    )), byrow = T, ncol = J)
    # Delete intercept term
    if (all(X[, 1] == 1)) {
      X = X[, -1]
    }
    for (k in 1:(J - 1)) {
      nam <- paste("A", k, sep = ".")
      assign(nam, A[, (k + 1)])
    }
    cox.formula = as.formula(paste("Surv(Y,DELTA)~X", paste(paste("A", 1:(
      J - 1
    ), sep = "."), collapse = "+"), sep = "+"))
    cox_q_model = coxph(cox.formula, method = "breslow")
    A.treat = matrix(0, ncol = J, nrow = length(Y))
    A.treat[, which(alpha == 1)] = 1
    treated_data = data.frame(cbind(X, A.treat[, -1]))
    colnames(treated_data) = names(cox_q_model$coefficients)
    A.control = matrix(0, ncol = J, nrow = length(Y))
    A.control[, which(alpha == -1)] = 1
    control_data = data.frame(cbind(X, A.control[, -1]))
    colnames(control_data) = names(cox_q_model$coefficients)
    # G computation
    km.object.treated = survfit(cox_q_model, newdata = treated_data)
    km.object.control = survfit(cox_q_model, newdata = control_data)
    
    point.est = rep(NA, 3)
    point.est[1] =
      g.calculation(km.object = km.object.treated,
                    end.time = NA,
                    type = 2) -
      g.calculation(km.object = km.object.control,
                    end.time = NA,
                    type = 2)
    
    point.est[2] =
      g.calculation(km.object = km.object.treated,
                    end.time = truncate,
                    type = 2) -
      g.calculation(km.object = km.object.control,
                    end.time = truncate,
                    type = 2)
    
    point.est[3] =
      g.calculation(km.object = km.object.treated,
                    end.time = truncate,
                    type = 1) -
      g.calculation(km.object = km.object.control,
                    end.time = truncate,
                    type = 1)
    
    return (point.est)
  }

cox.q.model.fit <-
  function(Y,
           DELTA,
           X,
           Z,
           alpha = c(-1, 1),
           truncate,
           boot.time = 250) {
    point.est =  cox.q.model.est(Y, DELTA, X, Z, alpha, truncate)
    boot.est = matrix(NA, ncol = 3, nrow = boot.time)
    n.size = length(Y)
    for (b in 1:boot.time) {
      b.id = sample(1:n.size, n.size, replace = T)
      b.Y = Y[b.id]
      b.DELTA = DELTA[b.id]
      b.X = X[b.id,]
      b.Z = Z[b.id]
      if (var(b.Z) == 0) {
        next
      }
      boot.est[b,] = cox.q.model.est(b.Y, b.DELTA, b.X, b.Z, alpha, truncate)
    }
    se = apply(boot.est, 2, sd, na.rm = T)
    names(point.est) = c("ASCE", "RACE", "SPCE")
    names(se) = c("ASCE", "RACE", "SPCE")
    return (list(est = point.est, se = se))
  }

msm.calculation <- function(km.object, end.time, type = 1) {
  time.grid = km.object$time
  surv.prob = km.object$surv
  levels = ncol(surv.prob)
  time.grid = c(0, time.grid) # Time grids for evaluation
  diff.time.grid = c(diff(time.grid)) # Difference of time
  if (type == 1) {
    final.index = max(which(time.grid <= end.time))
    return (surv.prob[final.index, ])
  } else{
    if (!is.na(end.time)) {
      final.index = max(which(time.grid <= end.time))
      diff.time.grid.temp = diff.time.grid
      diff.time.grid.temp[final.index] = end.time - time.grid[final.index]
    } else{
      final.index = length(diff.time.grid)
      diff.time.grid.temp = diff.time.grid
    }
    output = unlist(lapply(
      1:levels,
      FUN = function(y) {
        sum(unlist(lapply(
          1:final.index,
          FUN = function(x) {
            surv.prob[x, y] * diff.time.grid.temp[x]
          }
        )))
      }
    ))
    return (output)
  }
}

cox.msm.model.est <- function(Y, DELTA, X, Z, alpha = alpha, truncate) {
  # Delete intercept term
  if (all(X[, 1] == 1)) {
    X = X[, -1]
  }
  gps.model = multinom(Z ~ X,
                       maxit = 500,
                       Hess = TRUE,
                       trace = FALSE)
  e = gps.model$fitted.values
  if (length(unique(Z)) == 2) {
    e = cbind(1 - e, e)
  }
  ## Calculate Stablized Weights
  J = max(Z) + 1   # Number of arms
  w = matrix(rep(colMeans(e), length(Y)), ncol = J, byrow = T) / e
  A = matrix(unlist(lapply(
    Z,
    FUN = function(x) {
      as.numeric(x == 0:(J - 1))
    }
  )), byrow = T, ncol = J)
  s.weights = rowSums(A * w)
  for (k in 1:(J - 1)) {
    nam <- paste("A", k, sep = ".")
    assign(nam, A[, (k + 1)])
  }
  cox.formula = as.formula(paste("Surv(Y,DELTA)~", paste(paste("A", 1:(
    J - 1
  ), sep = "."), collapse = "+"), sep = ""))
  cox.msm.model = coxph(cox.formula, weights = s.weights, method = "breslow")
  A.unique = diag(rep(1, J))
  A.unique = data.frame(A.unique[, -1])
  colnames(A.unique) = names(cox.msm.model$coefficients)
  km.object = survfit(cox.msm.model, newdata = A.unique)
  
  point.est = rep(NA, 3)
  group.mean = msm.calculation(km.object, end.time = NA, type = 2)
  point.est[1] = sum(group.mean * alpha)
  group.mean = msm.calculation(km.object, end.time = truncate, type = 2)
  point.est[2] = sum(group.mean * alpha)
  group.mean = msm.calculation(km.object, end.time = truncate, type = 1)
  point.est[3] = sum(group.mean * alpha)
  return (point.est)
}

cox.msm.model.fit <-
  function(Y,
           DELTA,
           X,
           Z,
           alpha = c(-1, 1),
           truncate,
           boot.time = 250) {
    point.est =  cox.msm.model.est(Y, DELTA, X, Z, alpha, truncate)
    boot.est = matrix(NA, ncol = 3, nrow = boot.time)
    n.size = length(Y)
    for (b in 1:boot.time) {
      b.id = sample(1:n.size, n.size, replace = T)
      b.Y = Y[b.id]
      b.DELTA = DELTA[b.id]
      b.X = X[b.id, ]
      b.Z = Z[b.id]
      if (var(b.Z) == 0) {
        next
      }
      boot.est[b, ] = cox.msm.model.est(b.Y, b.DELTA, b.X, b.Z, alpha, truncate)
    }
    se = apply(boot.est, 2, sd, na.rm = T)
    names(point.est) = c("ASCE", "RACE", "SPCE")
    names(se) = c("ASCE", "RACE", "SPCE")
    return (list(est = point.est, se = se))
  }