#setwd("C:/Users/Shuxi ZENG/Dropbox/Fourth Year/OW_Survival/codebase/OW_Survival_CodeBase")
source("data_preprocessing.R")
source("cox_model.R")
source("PSW_pseudo.R")
source("Mao_Method_func.R")

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  print("No arguments supplied.")
} else{
  for (i in 1:length(args)) {
    eval(parse(text = args[[i]]))
  }
}
truncate.time=60
# Y = round(Y)
# alpha.arg = c(-1,1,0);truncate.time=60;
### PS Pseudo approach
if (first.diff) {
  alpha.arg = c(-1, 1, 0)
  if(approx.var==3){
    setwd("results_application_2/")
  }else{
    setwd("results_application_2_var/")
  }
} else{
  alpha.arg = c(-1, 0, 1)
  if(approx.var==3){
    setwd("results_application_3/")
  }else{
    setwd("results_application_3_var/")
  }
}

if (id == 1) {
  res.IPWC.ASCE = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    var.method = approx.var,
    weight.type = "IPW",
    estimand.type = "ASCE",
    ps.threshold = NA,
    alpha = alpha.arg
  )
  save(res.IPWC.ASCE, file = "IPWC_ASCE_Application.RData")
}

if (id == 2) {
  res.OW.ASCE = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    var.method = approx.var,
    weight.type = "OW",
    estimand.type = "ASCE",
    ps.threshold = NA,
    alpha = alpha.arg
  )
  save(res.OW.ASCE, file = "OW_ASCE_Application.RData")
}

if (id == 3) {
  res.IPWC.RACE = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    var.method = approx.var,
    weight.type = "IPW",
    estimand.type = "RACE",
    ps.threshold = 0.03,
    alpha = alpha.arg,
    evaluate.time = truncate.time
  )
  save(res.IPWC.RACE, file = "IPW_RACE_Application.RData")
}

if (id == 4) {
  res.OW.RACE = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    var.method = approx.var,
    weight.type = "OW",
    estimand.type = "RACE",
    alpha = alpha.arg,
    evaluate.time = truncate.time
  )
  save(res.OW.RACE, file = "OW_RACE_Application.RData")
}

if (id == 5) {
  res.IPWC.SPCE = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    var.method = 2,
    weight.type = "IPW",
    estimand.type = "SPCE",
    ps.threshold = 0.03,
    alpha = alpha.arg,
    evaluate.time = truncate.time
  )
  save(res.IPWC.SPCE, file = "IPW_SPCE_Application.RData")
}

if (id == 6) {
  res.OW.SPCE = PSW.pseudo(
    Y,
    Z,
    DELTA,
    X,
    var.method = 2,
    weight.type = "OW",
    estimand.type = "SPCE",
    alpha = alpha.arg,
    evaluate.time = truncate.time
  )
  save(res.OW.SPCE, file = "OW_SPCE_Application.RData")
}

if (id == 7){
  tmax = min(max(Y[Z == 0]), max(Y[Z == 2]))
  X.aug = cbind(X, as.numeric(Z == 1))
  colnames(X.aug) = c(paste("X", 0:(ncol(X) - 1), sep = "."), "Z.01")
  Z.aug = as.numeric(Z == 2)
  res.mao.IPW = estimand_analysis(
    X = X.aug[,-1],
    Z = Z.aug,
    Y = Y,
    delta = DELTA,
    weight.type = "IPW",
    t.trunc = truncate.time,
    tmax = tmax
  )
  
  # res.cox.q = cox.q.model.fit(
  #   Y,
  #   DELTA,
  #   X,
  #   Z,
  #   alpha = alpha.arg,
  #   truncate =  truncate.time,
  #   boot.time = 250
  # )
  
  res.cox.msm = cox.msm.model.fit(
    Y,
    DELTA,
    X,
    Z,
    alpha = alpha.arg,
    truncate =  truncate.time,
    boot.time = 250
  )
  
  save(res.cox.msm, res.mao.IPW, file = "Other_Application.RData")
}