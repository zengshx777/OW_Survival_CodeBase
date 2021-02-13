setwd("C:/Users/Shuxi ZENG/Dropbox/Fourth Year/OW_Survival/codebase/OW_Survival_CodeBase")

##Helper function to output all results
output_helper<-function(est,se,level=0.05)
{
  res=c(est,se,est-qnorm(1-level)*se,est+qnorm(1-level)*se,
        2*(1-pnorm(abs(est/se))))
  names(res)=c("Estimate","SE","95% CI","95% CI","p-value")
  return (res)
}
load("Application_12_comparison.RData")
results_12 = NULL
# results_12=rbind(results_12,
#   output_helper(res.OW.ASCE$tau,res.OW.ASCE$se),
#   output_helper(-res.IPWC.ASCE$tau,res.IPWC.ASCE$se),
#   output_helper(res.cox.q$est[1],res.cox.q$se[1]),
#   output_helper(res.cox.msm$est[1],res.cox.msm$se[1])
# #  output_helper(res.mao.IPW$ans$est[1],res.mao.IPW$ans$std[1])
# )

results_12=rbind(results_12,
  output_helper(res.OW.RACE$tau,res.OW.RACE$se),
  output_helper(res.IPWC.RACE$tau,res.IPWC.RACE$se),
  output_helper(res.cox.q$est[2],res.cox.q$se[2]),
  output_helper(res.cox.msm$est[2],res.cox.msm$se[2])
#  output_helper(res.mao.IPW$ans$est[2],res.mao.IPW$ans$std[2])
)

results_12=rbind(results_12,
                 output_helper(res.OW.SPCE$tau,res.OW.SPCE$se),
                 output_helper(res.IPWC.SPCE$tau,res.IPWC.SPCE$se),
                 output_helper(res.cox.q$est[3],res.cox.q$se[3]),
                 output_helper(res.cox.msm$est[3],res.cox.msm$se[3])
#                 output_helper(res.mao.IPW$ans$est[3],res.mao.IPW$ans$std[3])
)
load("Application_13_comparison.RData")
results_13 = NULL
# results_13=rbind(results_13,
#   output_helper(res.OW.ASCE$tau,res.OW.ASCE$se),
#   output_helper(-res.IPWC.ASCE$tau,res.IPWC.ASCE$se),
#   output_helper(res.cox.q$est[1],res.cox.q$se[1]),
#   output_helper(res.cox.msm$est[1],res.cox.msm$se[1])
#   #  output_helper(res.mao.IPW$ans$est[1],res.mao.IPW$ans$std[1])
# )

results_13=rbind(results_13,
                 output_helper(res.OW.RACE$tau,res.OW.RACE$se),
                 output_helper(res.IPWC.RACE$tau,res.IPWC.RACE$se),
                 output_helper(res.cox.q$est[2],res.cox.q$se[2]),
                 output_helper(res.cox.msm$est[2],res.cox.msm$se[2])
                 #  output_helper(res.mao.IPW$ans$est[2],res.mao.IPW$ans$std[2])
)

results_13=rbind(results_13,
                 output_helper(res.OW.SPCE$tau,res.OW.SPCE$se),
                 output_helper(res.IPWC.SPCE$tau,res.IPWC.SPCE$se),
                 output_helper(res.cox.q$est[3],res.cox.q$se[3]),
                 output_helper(res.cox.msm$est[3],res.cox.msm$se[3])
                 #                 output_helper(res.mao.IPW$ans$est[3],res.mao.IPW$ans$std[3])
)

load("Application_23_comparison.RData")
results_23 = NULL
# results_23=rbind(results_23,
#   output_helper(res.OW.ASCE$tau,res.OW.ASCE$se),
#   output_helper(-res.IPWC.ASCE$tau,res.IPWC.ASCE$se),
#   output_helper(res.cox.q$est[1],res.cox.q$se[1]),
#   output_helper(res.cox.msm$est[1],res.cox.msm$se[1])
#   #  output_helper(res.mao.IPW$ans$est[1],res.mao.IPW$ans$std[1])
# )

results_23=rbind(results_23,
                 output_helper(res.OW.RACE$tau,res.OW.RACE$se),
                 output_helper(res.IPWC.RACE$tau,res.IPWC.RACE$se),
                 output_helper(res.cox.q$est[2],res.cox.q$se[2]),
                 output_helper(res.cox.msm$est[2],res.cox.msm$se[2])
                 #  output_helper(res.mao.IPW$ans$est[2],res.mao.IPW$ans$std[2])
)

results_23=rbind(results_23,
                 output_helper(res.OW.SPCE$tau,res.OW.SPCE$se),
                 output_helper(res.IPWC.SPCE$tau,res.IPWC.SPCE$se),
                 output_helper(res.cox.q$est[3],res.cox.q$se[3]),
                 output_helper(res.cox.msm$est[3],res.cox.msm$se[3])
                 #                 output_helper(res.mao.IPW$ans$est[3],res.mao.IPW$ans$std[3])
)
library(xtable)
results = rbind(results_12,results_13, results_23)
row.names(results)=rep(c("OW","IPW","COX-Q","MSM"),nrow(results)/4)
print(xtable(results,digits=3))

