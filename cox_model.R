######## Cox-q model and MSM model
library(survival)
library(nnet)
### Apply G-computation based on the results from Cox model
g.calculation <- function(km.object, end.time, type = 1) {
  time.grid = km.object$time
  surv.prob = km.object$surv
  n.size = ncol(surv.prob)
  time.grid = c(0, time.grid) # Time grids for evaluation
  diff.time.grid = c(diff(time.grid)) # Difference of time
  if (type == 1) {
    final.index = max(which(time.grid <= end.time))
    final.index = min(final.index, nrow(surv.prob))
    return (mean(surv.prob[final.index,]))
  } else{
    if (!is.na(end.time)) {
      final.index = max(which(time.grid <= end.time))
      diff.time.grid.temp = diff.time.grid
      diff.time.grid.temp[final.index] = end.time - time.grid[final.index]
    } else{
      final.index = length(diff.time.grid)
      diff.time.grid.temp = diff.time.grid
    }
    final.index = min(nrow(surv.prob), final.index)
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
### Obtain point estimation from fitting Cox model
cox.q.model.est <-
  function(Y, DELTA, X, Z, truncate, ipw.weighted) {
    J = max(Z) + 1   # Number of arms
    mu.est.asce <- mu.est.race <- mu.est.spce <- rep(NA, J)
    
    
    A = matrix(unlist(lapply(
      Z,
      FUN = function(x) {
        as.numeric(x == 0:(J - 1))
      }
    )), byrow = T, ncol = J)
    # Delete intercept term
    if (all(X[, 1] == 1)) {
      X = X[,-1]
    }
    
    if (ipw.weighted) {
      gps.model = multinom(Z ~ X,
                           maxit = 500,
                           Hess = TRUE,
                           trace = FALSE)
      e = gps.model$fitted.values
      if (length(unique(Z)) == 2) {
        e = cbind(1 - e, e)
      }
      w = matrix(rep(colMeans(e), length(Y)),
                 ncol = J,
                 byrow = T) / e
      # Stablized weights
      s.weights = rowSums(A * w)
    }
    for (k in 1:(J - 1)) {
      nam <- paste("A", k, sep = ".")
      assign(nam, A[, (k + 1)])
    }
    cox.formula = as.formula(paste("Surv(Y,DELTA)~X", paste(paste("A", 1:(
      J - 1
    ), sep = "."), collapse = "+"), sep = "+"))
    if (ipw.weighted) {
      cox_q_model = coxph(cox.formula, weights = s.weights, method = "breslow")
    } else{
      cox_q_model = coxph(cox.formula, method = "breslow")
    }
    
    alpha.matrix = matrix(0, ncol = (J - 1), nrow = length(Y))
    
    temp.data = data.frame(cbind(X, alpha.matrix))
    colnames(temp.data) = names(cox_q_model$coefficients)
    km.object = survfit(cox_q_model, newdata = temp.data)
    mu.est.asce[1] = g.calculation(km.object = km.object,
                                   end.time = NA,
                                   type = 2)
    
    mu.est.race[1] = g.calculation(km.object = km.object,
                                   end.time = truncate,
                                   type = 2)
    
    mu.est.spce[1] = g.calculation(km.object = km.object,
                                   end.time = truncate,
                                   type = 1)
    
    for (k in 1:(J - 1)) {
      arm.matrix = alpha.matrix
      arm.matrix[, k] = 1
      temp.data = data.frame(cbind(X, arm.matrix))
      colnames(temp.data) = names(cox_q_model$coefficients)
      km.object = survfit(cox_q_model, newdata = temp.data)
      mu.est.asce[k + 1] = g.calculation(km.object = km.object,
                                         end.time = NA,
                                         type = 2)
      
      mu.est.race[k + 1] = g.calculation(km.object = km.object,
                                         end.time = truncate,
                                         type = 2)
      
      mu.est.spce[k + 1] = g.calculation(km.object = km.object,
                                         end.time = truncate,
                                         type = 1)
    }
    ind.alpha = combn(J, 2)
    point.est = matrix(NA, ncol = ncol(ind.alpha), nrow = 3)
    for (k in 1:ncol(point.est)) {
      point.est[1, k] = mu.est.asce[ind.alpha[2, k]] - mu.est.asce[ind.alpha[1, k]]
      point.est[2, k] = mu.est.race[ind.alpha[2, k]] - mu.est.race[ind.alpha[1, k]]
      point.est[3, k] = mu.est.spce[ind.alpha[2, k]] - mu.est.spce[ind.alpha[1, k]]
    }
    res = c(point.est[1, ], point.est[2, ], point.est[3, ])
    return (res)
  }
### Obtain point estimation and CI (bootstrap) from fitting Cox model
cox.q.model.fit <-
  function(Y,
           DELTA,
           X,
           Z,
           truncate,
           ipw.weighted = FALSE,
           boot.time = 250) {
    point.est =  cox.q.model.est(Y, DELTA, X, Z, truncate, ipw.weighted)
    alpha.ind = combn(max(Z) + 1, 2)
    comb.num = ncol(alpha.ind)
    boot.est = matrix(NA, ncol = 3 * comb.num, nrow = boot.time)
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
      boot.est[b, ] = cox.q.model.est(b.Y, b.DELTA, b.X, b.Z, truncate, ipw.weighted)
      # print(paste("== Bootstrap time", b, "=="))
    }
    se = apply(boot.est, 2, sd, na.rm = T)
    names(point.est) = rep(c("ASCE", "RACE", "SPCE"), each = comb.num)
    names(se) = rep(c("ASCE", "RACE", "SPCE"), each = comb.num)
    return (list(
      est = point.est,
      se = se,
      alpha.ind = alpha.ind
    ))
  }
### Apply MSM computation from the Cox model
msm.calculation <- function(km.object, end.time, type = 1) {
  time.grid = km.object$time
  surv.prob = km.object$surv
  levels = ncol(surv.prob)
  time.grid = c(0, time.grid) # Time grids for evaluation
  diff.time.grid = c(diff(time.grid)) # Difference of time
  if (type == 1) {
    final.index = max(which(time.grid <= end.time))
    final.index = min(final.index, nrow(surv.prob))
    return (surv.prob[final.index,])
  } else{
    if (!is.na(end.time)) {
      final.index = max(which(time.grid <= end.time))
      diff.time.grid.temp = diff.time.grid
      diff.time.grid.temp[final.index] = end.time - time.grid[final.index]
    } else{
      final.index = length(diff.time.grid)
      diff.time.grid.temp = diff.time.grid
    }
    final.index = min(nrow(surv.prob), final.index)
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
### Obtain MSM point estimation from the Cox model
cox.msm.model.est <-
  function(Y, DELTA, X, Z, truncate) {
    # Delete intercept term
    if (all(X[, 1] == 1)) {
      X = X[,-1]
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
    # Stablized weights
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
    A.unique = data.frame(A.unique[,-1])
    colnames(A.unique) = names(cox.msm.model$coefficients)
    km.object = survfit(cox.msm.model, newdata = A.unique)
    
    mu.est.asce = msm.calculation(km.object, end.time = NA, type = 2)
    mu.est.race = msm.calculation(km.object, end.time = truncate, type = 2)
    mu.est.spce = msm.calculation(km.object, end.time = truncate, type = 1)
    
    ind.alpha = combn(J, 2)
    point.est = matrix(NA, ncol = ncol(ind.alpha), nrow = 3)
    for (k in 1:ncol(point.est)) {
      point.est[1, k] = mu.est.asce[ind.alpha[2, k]] - mu.est.asce[ind.alpha[1, k]]
      point.est[2, k] = mu.est.race[ind.alpha[2, k]] - mu.est.race[ind.alpha[1, k]]
      point.est[3, k] = mu.est.spce[ind.alpha[2, k]] - mu.est.spce[ind.alpha[1, k]]
    }
    res = c(point.est[1, ], point.est[2, ], point.est[3, ])
    return (res)
  }
### Obtain MSM point estimation and CI from the Cox model
cox.msm.model.fit <-
  function(Y,
           DELTA,
           X,
           Z,
           truncate,
           boot.time = 250) {
    point.est =  cox.msm.model.est(Y, DELTA, X, Z, truncate)
    alpha.ind = combn((max(Z) + 1), 2)
    comb.num = ncol(alpha.ind)
    boot.est = matrix(NA, ncol = 3 * comb.num, nrow = boot.time)
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
      boot.est[b,] = cox.msm.model.est(b.Y, b.DELTA, b.X, b.Z, truncate)
      # print(paste("== Bootstrap time", b, "=="))
    }
    se = apply(boot.est, 2, sd, na.rm = T)
    names(point.est) = rep(c("ASCE", "RACE", "SPCE"), each = comb.num)
    names(se) = rep(c("ASCE", "RACE", "SPCE"), each = comb.num)
    return (list(
      est = point.est,
      se = se,
      alpha.ind = alpha.ind
    ))
  }