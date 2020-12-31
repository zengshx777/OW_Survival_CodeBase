setwd("C:/Users/Shuxi ZENG/Dropbox/Fourth Year/OW_Survival/codebase/OW_Survival_CodeBase")
source("data_preprocessing.R")
source("cox_model.R")
source("PSW_pseudo.R")
source("Mao_Method_func.R")

### PS Pseudo approach
alpha.arg = c(-1,1,0)
res.IPWC.ASCE = PSW.pseudo(
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

res.OW.ASCE = PSW.pseudo(
  Y,
  Z,
  DELTA,
  X,
  var.method = 2,
  weight.type = "OW",
  estimand.type = "ASCE",
  alpha = alpha.arg
)

res.IPWC.RACE = PSW.pseudo(
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

res.OW.RACE = PSW.pseudo(
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
  evaluate.time = truncate
)

res.OW.SPCE = PSW.pseudo(
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



tmax = min(max(Y[Z == 0]), max(Y[Z == 2]))
X.aug = cbind(X, as.numeric(Z == 1))
colnames(X.aug) = c(paste("X",0:(ncol(X)-1),sep="."),"Z.01")
Z.aug = as.numeric(Z == 2)
res.mao.IPW = estimand_analysis(
  X = X.aug[, -1],
  Z = Z.aug,
  Y = Y,
  delta = DELTA,
  weight.type = "IPW",
  t.trunc = truncate,
  tmax = tmax
)

res.cox.q = cox.q.model.fit(
  Y,
  DELTA,
  X,
  Z,
  alpha = alpha.arg,
  truncate =  truncate,
  boot.time = 250
)

res.cox.msm = cox.msm.model.fit(
  Y,
  DELTA,
  X,
  Z,
  alpha = alpha.arg,
  truncate =  truncate,
  boot.time = 250
)
