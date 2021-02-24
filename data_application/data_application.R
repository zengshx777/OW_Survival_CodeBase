source("data_preprocessing.R")
source("../estimators/cox_model.R")
source("../estimators/PSW_pseudo.R")
# source("Mao_Method_func.R")
# change time points for IPW-PO and OW-PO
grids = seq(1,by=0.5,length=220)
for(truncate.time in grids){
res.weighting = PSW.pseudo.OW.IPW(
  Y,
  Z,
  DELTA,
  X,
  var.method = 2,
  estimand.type = "SPCE",
  evaluate.time = grids[i]
)
save(res.weighting, file = paste(grids[i],"PSW_SPCE_Application.RData",sep="_"))
}
truncate.time=60

# COX-G
res.cox.q = cox.q.model.fit(
  Y,
  DELTA,
  X,
  Z,
  truncate =  truncate.time,
  boot.time = 250
)
save(res.cox.msm, file = "Cox_G_application_results.RData")

# COX-IPW
res.cox.msm = cox.msm.model.fit(
  Y,
  DELTA,
  X,
  Z,
  truncate =  truncate.time,
  boot.time = 250
)

save(res.cox.msm, file = "Cox_IPW_application_results.RData")