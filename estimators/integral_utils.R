### Calculate Cumulative Hazard
calculate.obs.approximate <-
  function(Y,
           DELTA,
           evaluate.time,
           estimand.type = "ASCE",
           var.method = 2) {
    if(estimand.type == "ASCE"){
      evaluate.time = (max(Y) + max(Y[Y != max(Y)])) / 2
    }
    N = length(Y)
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
    
    
    ### Construct observation for influence function
    if (var.method == 1) {
      ### Most accurate, taking account of second derivative
      # adjust.matrix = as.vector((w.est * A) %*% t(alpha.matrix))
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
    } else {
      ### Less accurate, ignoring the second order derivative, overestimate
      obs.approxi = zero.derivative() + unlist(lapply(1:N, first.derivative))
    } 
    return (obs.approxi)
  }