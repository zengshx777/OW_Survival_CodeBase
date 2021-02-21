##
setwd("C:/Users/Shuxi ZENG/Dropbox/Fourth Year/OW_Survival/codebase/application_PSW_curve_drawing/results_application_surv_prob_PSW")
# setwd("C:/Users/Shuxi ZENG/Dropbox/Fourth Year/OW_Survival/codebase/results_application_surv_prob/results_application_surv_prob")
TIME.LENGTH = 220
time.grid = seq(1,by=0.5,length=TIME.LENGTH)
ow.race.diff = matrix(NA,ncol=3,nrow=TIME.LENGTH)
ow.race.se = matrix(NA,ncol=3,nrow=TIME.LENGTH)
ipw.race.diff = matrix(NA,ncol=3,nrow=TIME.LENGTH)
ipw.race.se = matrix(NA,ncol=3,nrow=TIME.LENGTH)
ow.surv.diff = matrix(NA,ncol=3,nrow=TIME.LENGTH)
ow.surv.se = matrix(NA,ncol=3,nrow=TIME.LENGTH)
ipw.surv.diff = matrix(NA,ncol=3,nrow=TIME.LENGTH)
ipw.surv.se = matrix(NA,ncol=3,nrow=TIME.LENGTH)
ipw.surv.mu = matrix(NA,ncol=3,nrow=TIME.LENGTH)
ipw.surv.mu.se = matrix(NA,ncol=3,nrow=TIME.LENGTH)
ow.surv.mu = matrix(NA,ncol=3,nrow=TIME.LENGTH)
ow.surv.mu.se = matrix(NA,ncol=3,nrow=TIME.LENGTH)
for (i in 1:TIME.LENGTH){
  load(paste(time.grid[i],"PSW_RACE_Application.RData",sep="_"))

  ipw.race.diff[i,] = res.weighting$tau.ipw
  ipw.race.se[i,] = res.weighting$se.ipw
  ow.race.diff[i,] = res.weighting$tau.ow
  ow.race.se[i,] = res.weighting$se.ow

  exist = tryCatch({
    load(paste(time.grid[i],"PSW_SPCE_Application.RData",sep="_"))
              1
            },error=function(e){
                return (0)
              })
  if(exist==1){
    # ipw.surv.diff[i,] = res.weighting$tau.ipw
    # ipw.surv.se[i,] = res.weighting$se.ipw
    # ow.surv.diff[i,] = res.weighting$tau.ow
    # ow.surv.se[i,] = res.weighting$se.ow

    ipw.surv.mu[i,] = res.weighting$mu.ipw
    ipw.surv.mu.se[i,] = sqrt(diag(res.weighting$Sigma.matrix.ipw))
    ow.surv.mu[i,] = res.weighting$mu.ow
    ow.surv.mu.se[i,] = sqrt(diag(res.weighting$Sigma.matrix.ow))
  }
  load(paste(time.grid[i],"OW_SPCE_Application.RData",sep="_"))
  ow.surv.diff[i,] = res.OW.SPCE$tau
  ow.surv.se[i,] = res.OW.SPCE$se
  load(paste(time.grid[i],"IPW_SPCE_Application.RData",sep="_"))
  ipw.surv.diff[i,] = res.IPWC.SPCE$tau
  ipw.surv.se[i,] = res.IPWC.SPCE$se
}
for (j in 1:3){
  order = order(ipw.surv.mu[,j],decreasing = T)
  ipw.surv.mu[,j] = ipw.surv.mu[order,j]
  ipw.surv.mu.se[order,j] = ipw.surv.mu.se[order,j]
  # ipw.surv.diff[,j] = ipw.surv.diff[order,j]
  # ipw.surv.se[,j] = ipw.surv.se[order,j]
  order = order(ow.surv.mu[,j],decreasing = T)
  ow.surv.mu[,j] = ow.surv.mu[order,j]
  ow.surv.mu.se[,j] = ow.surv.mu.se[order,j]
  # ow.surv.diff[,j] = ow.surv.diff[order,j]
  # ow.surv.se[,j] = ow.surv.se[order,j]
}
temp = ipw.surv.diff
ipw.surv.diff = ow.surv.diff
ow.surv.diff = temp
temp = ipw.surv.se
ipw.surv.se= ow.surv.se
ow.surv.se = temp
{
setwd("C:/Users/Shuxi ZENG/Dropbox/Fourth Year/OW_Survival/codebase/OW_Survival_CodeBase")
output_helper<-function(est,se,level=0.05)
{
  res=c(est,se,est-qnorm(1-level)*se,est+qnorm(1-level)*se,
        2*(1-pnorm(abs(est/se))))
  names(res)=c("Estimate","SE","95% CI","95% CI","p-value")
  return (res)
}

load("Application_12_comparison.RData")
results_12=rbind(
                 output_helper(ow.race.diff[which(time.grid==60),1],ow.race.se[which(time.grid==60),1]),
                 output_helper(ipw.race.diff[which(time.grid==60),1],ipw.race.se[which(time.grid==60),1]),
                 output_helper(res.cox.q$est[2],res.cox.q$se[2]),
                 output_helper(res.cox.msm$est[2],res.cox.msm$se[2])
                 #  output_helper(res.mao.IPW$ans$est[2],res.mao.IPW$ans$std[2])
)
results_12=rbind(results_12,
  output_helper(ow.surv.diff[which(time.grid==60),1],ow.surv.se[which(time.grid==60),1]),
  output_helper(ipw.surv.diff[which(time.grid==60),1],ipw.surv.se[which(time.grid==60),1]),
  output_helper(res.cox.q$est[3],res.cox.q$se[3]),
  output_helper(res.cox.msm$est[3],res.cox.msm$se[3])
  #  output_helper(res.mao.IPW$ans$est[2],res.mao.IPW$ans$std[2])
)

load("Application_13_comparison.RData")
results_13=rbind(
  output_helper(ow.race.diff[which(time.grid==60),2],ow.race.se[which(time.grid==60),2]),
  output_helper(ipw.race.diff[which(time.grid==60),2],ipw.race.se[which(time.grid==60),2]),
  output_helper(res.cox.q$est[2],res.cox.q$se[2]),
  output_helper(res.cox.msm$est[2],res.cox.msm$se[2])
  #  output_helper(res.mao.IPW$ans$est[2],res.mao.IPW$ans$std[2])
)
results_13=rbind(results_13,
                 output_helper(ow.surv.diff[which(time.grid==60),2],ow.surv.se[which(time.grid==60),2]),
                 output_helper(ipw.surv.diff[which(time.grid==60),2],ipw.surv.se[which(time.grid==60),2]),
                 output_helper(res.cox.q$est[3],res.cox.q$se[3]),
                 output_helper(res.cox.msm$est[3],res.cox.msm$se[3])
                 #  output_helper(res.mao.IPW$ans$est[2],res.mao.IPW$ans$std[2])
)

load("Application_23_comparison.RData")
results_23=rbind(
  output_helper(ow.race.diff[which(time.grid==60),3],ow.race.se[which(time.grid==60),3]),
  output_helper(ipw.race.diff[which(time.grid==60),3],ipw.race.se[which(time.grid==60),3]),
  output_helper(res.cox.q$est[2],res.cox.q$se[2]),
  output_helper(res.cox.msm$est[2],res.cox.msm$se[2])
  #  output_helper(res.mao.IPW$ans$est[2],res.mao.IPW$ans$std[2])
)
results_23=rbind(results_23,
                 output_helper(ow.surv.diff[which(time.grid==60),3],ow.surv.se[which(time.grid==60),3]),
                 output_helper(ipw.surv.diff[which(time.grid==60),3],ipw.surv.se[which(time.grid==60),3]),
                 output_helper(res.cox.q$est[3],res.cox.q$se[3]),
                 output_helper(res.cox.msm$est[3],res.cox.msm$se[3])
                 #  output_helper(res.mao.IPW$ans$est[2],res.mao.IPW$ans$std[2])
)

library(xtable)
results = rbind(results_12,results_13, results_23)
row.names(results)=rep(c("OW","IPW","COX-Q","MSM"),nrow(results)/4)
print(xtable(results,digits=3))
}


pdf("race_diff.pdf",width=12,height=4)
par(mfrow=c(1,3))
plot(time.grid,ow.race.diff[,1],type='s',xlab="Months after treatment",
     ylab="EBRT+AD vs RP Comparison, RACE",ylim=c(-5.5,0.5),lwd=2)
lines(time.grid,ow.race.diff[,1]+qnorm(0.975)*ow.race.se[,1],lty=2,lwd=2,type='s')
lines(time.grid,ow.race.diff[,1]-qnorm(0.975)*ow.race.se[,1],lty=2,lwd=2,type='s')
lines(time.grid,ipw.race.diff[,1],type='s')
lines(time.grid,ipw.race.diff[,1]+qnorm(0.975)*ipw.race.se[,1],lty=2,type='s')
lines(time.grid,ipw.race.diff[,1]-qnorm(0.975)*ipw.race.se[,1],lty=2,type='s')
abline(h=0,lty=2)

plot(time.grid,ow.race.diff[,2],type='s',xlab="Months after treatment",
     ylab="EBRT+brachy±AD vs RP Comparison, RACE",ylim=c(-3,2.5),lwd=2)
lines(time.grid,ow.race.diff[,2]+qnorm(0.975)*ow.race.se[,2],lty=2,lwd=2)
lines(time.grid,ow.race.diff[,2]-qnorm(0.975)*ow.race.se[,2],lty=2,lwd=2)
lines(time.grid,ipw.race.diff[,2])
lines(time.grid,ipw.race.diff[,2]+qnorm(0.975)*ipw.race.se[,2],lty=2)
lines(time.grid,ipw.race.diff[,2]-qnorm(0.975)*ipw.race.se[,2],lty=2)
abline(h=0,lty=2)
legend("bottom",legend = c("OW","IPW"),lty=1,lwd=c(2,1),horiz = T)

plot(time.grid,ow.race.diff[,3],type='s',xlab="Months after treatment",
     ylab="EBRT+brachy±AD vs EBRT+AD Comparison, RACE", ylim=c(-0.5,4),lwd=2)
lines(time.grid,ow.race.diff[,3]+qnorm(0.975)*ow.race.se[,3],lty=2,lwd=2)
lines(time.grid,ow.race.diff[,3]-qnorm(0.975)*ow.race.se[,3],lty=2,lwd=2)
lines(time.grid,ipw.race.diff[,3])
lines(time.grid,ipw.race.diff[,3]+qnorm(0.975)*ipw.race.se[,3],lty=2)
lines(time.grid,ipw.race.diff[,3]-qnorm(0.975)*ipw.race.se[,3],lty=2)
abline(h=0,lty=2)
dev.off()

pdf("Survival_prob_application.pdf",width=10,height = 5)
par(mfrow=c(1,2))
non.na.id = which(!is.na(ipw.surv.diff[,1]))
plot(time.grid[non.na.id],ipw.surv.mu[non.na.id,1],type='s',xlab="Months after treatment",
     ylab="Survival function, IPW adjusted",col="red", ylim=c(0,1),lwd=2)
lines(time.grid[non.na.id],ipw.surv.mu[non.na.id,1]+qnorm(0.975)*ipw.surv.mu.se[non.na.id,1],lty=2,col="red",lwd=1.5)
lines(time.grid[non.na.id],ipw.surv.mu[non.na.id,1]-qnorm(0.975)*ipw.surv.mu.se[non.na.id,1],lty=2,col="red",lwd=1.5)
lines(time.grid[non.na.id],ipw.surv.mu[non.na.id,2],col="blue",lwd=2)
lines(time.grid[non.na.id],ipw.surv.mu[non.na.id,2]+qnorm(0.975)*ipw.surv.mu.se[non.na.id,2],lty=2,col="blue",lwd=1.5)
lines(time.grid[non.na.id],ipw.surv.mu[non.na.id,2]-qnorm(0.975)*ipw.surv.mu.se[non.na.id,2],lty=2,col="blue",lwd=1.5)
lines(time.grid[non.na.id],ipw.surv.mu[non.na.id,3],col="green",lwd=2)
lines(time.grid[non.na.id],ipw.surv.mu[non.na.id,3]+qnorm(0.975)*ipw.surv.mu.se[non.na.id,3],lty=2,col="green",lwd=1.5)
lines(time.grid[non.na.id],ipw.surv.mu[non.na.id,3]-qnorm(0.975)*ipw.surv.mu.se[non.na.id,3],lty=2,col="green",lwd=1.5)
legend("bottomleft",legend = c("RP","","EBRT+AD","","EBRT+brachy±AD",""),
       lty=rep(c(1,2),3),col=rep(c("red","blue","green"),each=2),lwd=2)

plot(time.grid[non.na.id],ow.surv.mu[non.na.id,1],type='s',xlab="Months after treatment",
     ylab="Survival function, OW adjusted",col="red", ylim=c(0,1),lwd=2)
lines(time.grid[non.na.id],ow.surv.mu[non.na.id,1]+qnorm(0.975)*ow.surv.mu.se[non.na.id,1],lty=2,col="red",lwd=1.5)
lines(time.grid[non.na.id],ow.surv.mu[non.na.id,1]-qnorm(0.975)*ow.surv.mu.se[non.na.id,1],lty=2,col="red",lwd=1.5)
lines(time.grid[non.na.id],ow.surv.mu[non.na.id,2],col="blue",lwd=2)
lines(time.grid[non.na.id],ow.surv.mu[non.na.id,2]+qnorm(0.975)*ow.surv.mu.se[non.na.id,2],lty=2,col="blue",lwd=1.5)
lines(time.grid[non.na.id],ow.surv.mu[non.na.id,2]-qnorm(0.975)*ow.surv.mu.se[non.na.id,2],lty=2,col="blue",lwd=1.5)
lines(time.grid[non.na.id],ow.surv.mu[non.na.id,3],col="green",lwd=2)
lines(time.grid[non.na.id],ow.surv.mu[non.na.id,3]+qnorm(0.975)*ow.surv.mu.se[non.na.id,3],lty=2,col="green",lwd=1.5)
lines(time.grid[non.na.id],ow.surv.mu[non.na.id,3]-qnorm(0.975)*ow.surv.mu.se[non.na.id,3],lty=2,col="green",lwd=1.5)
# legend("bottomleft",legend = c("RP","","EBRT+AD","","EBRT+brachy±AD",""),
#        lty=rep(c(1,2),3),col=rep(c("red","blue","green"),each=2),lwd=2)
dev.off()


pdf("Surv_prob_diff.pdf",width=12,height=4)
par(mfrow=c(1,3))
non.na.id = which(!is.na(ow.surv.diff[,1]))
plot(time.grid[non.na.id],ow.surv.diff[non.na.id,1],type='s',ylim=c(-0.4,0.35),xlab="Months after treatment",
     ylab="EBRT+AD vs RP Comparison, SPCE",lwd=1)
lines(time.grid[non.na.id],ow.surv.diff[non.na.id,1]+qnorm(0.975)*ow.surv.se[non.na.id,1],lty=2,type='s')
lines(time.grid[non.na.id],ow.surv.diff[non.na.id,1]-qnorm(0.975)*ow.surv.se[non.na.id,1],lty=2,type='s')
lines(time.grid[non.na.id],ipw.surv.diff[non.na.id,1],type='s',lwd=2)
lines(time.grid[non.na.id],ipw.surv.diff[non.na.id,1]+qnorm(0.975)*ipw.surv.se[non.na.id,1],lty=2,type='s',lwd=2)
lines(time.grid[non.na.id],ipw.surv.diff[non.na.id,1]-qnorm(0.975)*ipw.surv.se[non.na.id,1],lty=2,type='s',lwd=2)
abline(h=0,lty=2)


non.na.id = which(!is.na(ow.surv.diff[,2]))
plot(time.grid[non.na.id],ow.surv.diff[non.na.id,2],type='s',ylim=c(-0.4,0.35),xlab="Months after treatment",
     ylab="EBRT+brachy±AD vs RP Comparison, SPCE")
lines(time.grid[non.na.id],ow.surv.diff[non.na.id,2]+qnorm(0.975)*ow.surv.se[non.na.id,2],lty=2,type='s')
lines(time.grid[non.na.id],ow.surv.diff[non.na.id,2]-qnorm(0.975)*ow.surv.se[non.na.id,2],lty=2,type='s')
lines(time.grid[non.na.id],ipw.surv.diff[non.na.id,2],lwd=2,type='s')
lines(time.grid[non.na.id],ipw.surv.diff[non.na.id,2]+qnorm(0.975)*ipw.surv.se[non.na.id,2],lty=2,lwd=2,type='s')
lines(time.grid[non.na.id],ipw.surv.diff[non.na.id,2]-qnorm(0.975)*ipw.surv.se[non.na.id,2],lty=2,lwd=2,type='s')
abline(h=0,lty=2)
legend("bottom",legend = c("OW","IPW"),lty=1,lwd= c(1,2), horiz = T)

non.na.id = which(!is.na(ow.surv.diff[,3]))
plot(time.grid[non.na.id],ow.surv.diff[non.na.id,3],type='s',ylim=c(-0.4,0.35),xlab="Months after treatment",
     ylab="EBRT+brachy±AD vs EBRT+AD Comparison, SPCE")
lines(time.grid[non.na.id],ow.surv.diff[non.na.id,3]+qnorm(0.975)*ow.surv.se[non.na.id,3],lty=2,type='s')
lines(time.grid[non.na.id],ow.surv.diff[non.na.id,3]-qnorm(0.975)*ow.surv.se[non.na.id,3],lty=2,type='s')
lines(time.grid[non.na.id],ipw.surv.diff[non.na.id,3],lwd=2,type='s')
lines(time.grid[non.na.id],ipw.surv.diff[non.na.id,3]+qnorm(0.975)*ipw.surv.se[non.na.id,3],lty=2,lwd=2,type='s')
lines(time.grid[non.na.id],ipw.surv.diff[non.na.id,3]-qnorm(0.975)*ipw.surv.se[non.na.id,3],lty=2,lwd=2,type='s')
abline(h=0,lty=2)
dev.off()

library(ggplot2)
library(gridExtra)
pdf("race_diff_ggplot.pdf",width=12,height=4)
data = cbind(time.grid,ow.race.diff,ow.race.se,ipw.race.diff,ipw.race.se)
colnames(data) = c("time","OW.RACE.DIFF.1","OW.RACE.DIFF.2","OW.RACE.DIFF.3",
                   "OW.RACE.SE.1","OW.RACE.SE.2","OW.RACE.SE.3",
                   "IPW.RACE.DIFF.1","IPW.RACE.DIFF.2","IPW.RACE.DIFF.3",
                   "IPW.RACE.SE.1","IPW.RACE.SE.2","IPW.RACE.SE.3")
data =  as.data.frame(data)

p1 = ggplot(aes(time),data = data) + 
  geom_ribbon(aes(ymin = OW.RACE.DIFF.1 - qnorm(0.975)*OW.RACE.SE.1,
                  ymax = OW.RACE.DIFF.1 + qnorm(0.975)*OW.RACE.SE.1),    # shadowing cnf intervals
              fill = "firebrick",alpha=0.5) + 
  geom_line(aes(x=time, y=OW.RACE.DIFF.1,color="OW"),size = 1)   +
geom_ribbon(aes(ymin = IPW.RACE.DIFF.1 - qnorm(0.975)*IPW.RACE.SE.1,
                ymax = IPW.RACE.DIFF.1 + qnorm(0.975)*IPW.RACE.SE.1),    # shadowing cnf intervals
            fill = "steelblue2",alpha=0.5) + 
  geom_line(aes(time, IPW.RACE.DIFF.1,color = "IPW"),
            size = 1) +xlab("Months after treatment") +
  scale_color_manual(values = c(
    'OW' = 'red',
    'IPW' = 'darkblue'))+labs(color = 'Method')+
  ylab("EBRT+AD vs RP Comparison, restricted mean difference")+geom_hline(yintercept = 0, linetype="dashed")+theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position =  c(0.2, 0.15))

p2 =ggplot(aes(time),data = data) + 
  geom_ribbon(aes(ymin = OW.RACE.DIFF.2 - qnorm(0.975)*OW.RACE.SE.2,
                  ymax = OW.RACE.DIFF.2 + qnorm(0.975)*OW.RACE.SE.2),    # shadowing cnf intervals
              fill = "firebrick",alpha=0.5) + 
  geom_line(aes(x=time, y=OW.RACE.DIFF.2,color="OW"),size = 1)   +
  geom_ribbon(aes(ymin = IPW.RACE.DIFF.2 - qnorm(0.975)*IPW.RACE.SE.2,
                  ymax = IPW.RACE.DIFF.2 + qnorm(0.975)*IPW.RACE.SE.2),    # shadowing cnf intervals
              fill = "steelblue2",alpha=0.5) + 
  geom_line(aes(time, IPW.RACE.DIFF.2,color = "IPW"),
            size = 1) +xlab("Months after treatment") +
  scale_color_manual(values = c(
    'OW' = 'red',
    'IPW' = 'darkblue'))+labs(color = 'Method')+
  ylab("EBRT+brachy±AD vs RP Comparison, RACE")+geom_hline(yintercept = 0, linetype="dashed")+theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")


p3 = ggplot(aes(time),data = data) + 
  geom_ribbon(aes(ymin = OW.RACE.DIFF.3 - qnorm(0.975)*OW.RACE.SE.3,
                  ymax = OW.RACE.DIFF.3 + qnorm(0.975)*OW.RACE.SE.3),    # shadowing cnf intervals
              fill = "firebrick",alpha=0.5) + 
  geom_line(aes(x=time, y=OW.RACE.DIFF.3,color="OW"),size = 1)   +
  geom_ribbon(aes(ymin = IPW.RACE.DIFF.3 - qnorm(0.975)*IPW.RACE.SE.3,
                  ymax = IPW.RACE.DIFF.3 + qnorm(0.975)*IPW.RACE.SE.3),    # shadowing cnf intervals
              fill = "steelblue2",alpha=0.5) + 
  geom_line(aes(time, IPW.RACE.DIFF.3,color = "IPW"),
            size = 1) +xlab("Months after treatment") +
  scale_color_manual(values = c(
    'OW' = 'red',
    'IPW' = 'darkblue'))+labs(color = 'Method')+
  ylab("EBRT+brachy±AD vs EBRT+AD Comparison, RACE")+geom_hline(yintercept = 0, linetype="dashed")+  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")
grid.arrange(p1,p2,p3,ncol=3)
dev.off()



pdf("spce_diff_ggplot.pdf",width=12,height=4)
data = cbind(time.grid,ow.surv.diff,ow.surv.se,ipw.surv.diff,ipw.surv.se)
colnames(data) = c("time","OW.spce.DIFF.1","OW.spce.DIFF.2","OW.spce.DIFF.3",
                   "OW.spce.SE.1","OW.spce.SE.2","OW.spce.SE.3",
                   "IPW.spce.DIFF.1","IPW.spce.DIFF.2","IPW.spce.DIFF.3",
                   "IPW.spce.SE.1","IPW.spce.SE.2","IPW.spce.SE.3")
data =  as.data.frame(data)

q1 = ggplot(aes(time),data = data) + 
  geom_ribbon(aes(ymin = OW.spce.DIFF.1 - qnorm(0.975)*OW.spce.SE.1,
                  ymax = OW.spce.DIFF.1 + qnorm(0.975)*OW.spce.SE.1),    # shadowing cnf intervals
              fill = "firebrick",alpha=0.5) + 
  geom_line(aes(x=time, y=OW.spce.DIFF.1,color="OW"),size = 1)   +
  geom_ribbon(aes(ymin = IPW.spce.DIFF.1 - qnorm(0.975)*IPW.spce.SE.1,
                  ymax = IPW.spce.DIFF.1 + qnorm(0.975)*IPW.spce.SE.1),    # shadowing cnf intervals
              fill = "steelblue2",alpha=0.5) + 
  geom_line(aes(time, IPW.spce.DIFF.1,color = "IPW"),
            size = 1) +xlab("Months after treatment") +
  scale_color_manual(values = c(
    'OW' = 'red',
    'IPW' = 'darkblue'))+labs(color = 'Method')+
  ylab("EBRT+AD vs RP Comparison, survival prob difference")+geom_hline(yintercept = 0, linetype="dashed")+theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = c(0.2, 0.15))

q2 =ggplot(aes(time),data = data) + 
  geom_ribbon(aes(ymin = OW.spce.DIFF.2 - qnorm(0.975)*OW.spce.SE.2,
                  ymax = OW.spce.DIFF.2 + qnorm(0.975)*OW.spce.SE.2),    # shadowing cnf intervals
              fill = "firebrick",alpha=0.5) + 
  geom_line(aes(x=time, y=OW.spce.DIFF.2,color="OW"),size = 1)   +
  geom_ribbon(aes(ymin = IPW.spce.DIFF.2 - qnorm(0.975)*IPW.spce.SE.2,
                  ymax = IPW.spce.DIFF.2 + qnorm(0.975)*IPW.spce.SE.2),    # shadowing cnf intervals
              fill = "steelblue2",alpha=0.5) + 
  geom_line(aes(time, IPW.spce.DIFF.2,color = "IPW"),
            size = 1) +xlab("Months after treatment") +
  scale_color_manual(values = c(
    'OW' = 'red',
    'IPW' = 'darkblue'))+labs(color = 'Method')+
  ylab("EBRT+brachy±AD vs RP Comparison, SPCE")+geom_hline(yintercept = 0, linetype="dashed")+theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")


q3 = ggplot(aes(time),data = data) + 
  geom_ribbon(aes(ymin = OW.spce.DIFF.3 - qnorm(0.975)*OW.spce.SE.3,
                  ymax = OW.spce.DIFF.3 + qnorm(0.975)*OW.spce.SE.3),    # shadowing cnf intervals
              fill = "firebrick",alpha=0.5) + 
  geom_line(aes(x=time, y=OW.spce.DIFF.3,color="OW"),size = 1)   +
  geom_ribbon(aes(ymin = IPW.spce.DIFF.3 - qnorm(0.975)*IPW.spce.SE.3,
                  ymax = IPW.spce.DIFF.3 + qnorm(0.975)*IPW.spce.SE.3),    # shadowing cnf intervals
              fill = "steelblue2",alpha=0.5) + 
  geom_line(aes(time, IPW.spce.DIFF.3,color = "IPW"),
            size = 1) +xlab("Months after treatment") +
  scale_color_manual(values = c(
    'OW' = 'red',
    'IPW' = 'darkblue'))+labs(color = 'Method')+
  ylab("EBRT+brachy±AD vs EBRT+AD Comparison, SPCE")+geom_hline(yintercept = 0, linetype="dashed")+  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")
grid.arrange(q1,q2,q3,ncol=3)
dev.off()

pdf("surv_prob_ggplot.pdf",width=10,height = 5)
# NA.ID = which(is.na(ow.surv.mu[,1]))
data = cbind(time.grid,ow.surv.mu,
             ow.surv.mu.se,
             ipw.surv.mu,
             ipw.surv.mu.se)
colnames(data) = c("time","OW.SURV.1","OW.SURV.2","OW.SURV.3",
                   "OW.SURV.SE.1","OW.SURV.SE.2","OW.SURV.SE.3",
                   "IPW.SURV.1","IPW.SURV.2","IPW.SURV.3",
                   "IPW.SURV.SE.1","IPW.SURV.SE.2","IPW.SURV.SE.3")
data =  as.data.frame(data)

s1 = ggplot(aes(time),data = data) + 
  geom_ribbon(aes(ymin = IPW.SURV.1 - qnorm(0.975)*IPW.SURV.SE.1,
                  ymax = IPW.SURV.1 + qnorm(0.975)*IPW.SURV.SE.1),    # shadIPWing cnf intervals
              fill = "magenta",alpha=0.2) + 
  geom_line(aes(x=time, y=IPW.SURV.1,color="RP"),size = 1)   +
  geom_ribbon(aes(ymin = IPW.SURV.2 - qnorm(0.975)*IPW.SURV.SE.2,
                  ymax = IPW.SURV.2 + qnorm(0.975)*IPW.SURV.SE.2),    # shadIPWing cnf intervals
              fill = "dodgerblue",alpha=0.2) + 
  geom_line(aes(time, IPW.SURV.2,color = "EBRT+AD"),size = 1) + 
  geom_ribbon(aes(ymin = IPW.SURV.3 - qnorm(0.975)*IPW.SURV.SE.3,
                  ymax = IPW.SURV.3 + qnorm(0.975)*IPW.SURV.SE.3),    # shadIPWing cnf intervals
              fill = "orange",alpha=0.2) + 
  geom_line(aes(time, IPW.SURV.3,color = "EBRT+brachy±AD"),size = 1) + 
  xlab("Months after treatment") +
  scale_color_manual(values = c(
    'RP' = 'purple','EBRT+AD' = 'darkblue','EBRT+brachy±AD' = 'darkorange4'))+labs(color = '')+
  ylab("Survival Prob, adjusted by IPW")+theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = c(0.2, 0.15))


s2 = ggplot(aes(time),data = data) + 
  geom_ribbon(aes(ymin = OW.SURV.1 - qnorm(0.975)*OW.SURV.SE.1,
                  ymax = OW.SURV.1 + qnorm(0.975)*OW.SURV.SE.1),    # shadowing cnf intervals
              fill = "magenta",alpha=0.2) + 
  geom_line(aes(x=time, y=OW.SURV.1,color="RP"),size = 1)   +
  geom_ribbon(aes(ymin = OW.SURV.2 - qnorm(0.975)*OW.SURV.SE.2,
                  ymax = OW.SURV.2 + qnorm(0.975)*OW.SURV.SE.2),    # shadowing cnf intervals
              fill = "dodgerblue",alpha=0.2) + 
  geom_line(aes(time, OW.SURV.2,color = "EBRT+AD"),size = 1) + 
  geom_ribbon(aes(ymin = OW.SURV.3 - qnorm(0.975)*OW.SURV.SE.3,
                  ymax = OW.SURV.3 + qnorm(0.975)*OW.SURV.SE.3),    # shadowing cnf intervals
              fill = "orange",alpha=0.2) + 
  geom_line(aes(time, OW.SURV.3,color = "EBRT+brachy±AD"),size = 1) + 
  xlab("Months after treatment") +
  scale_color_manual(values = c(
    'RP' = 'purple','EBRT+AD' = 'darkblue','EBRT+brachy±AD' = 'darkorange4'))+labs(color = '')+
  ylab("Survival Prob, adjusted by OW")+theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none")
grid.arrange(s1,s2,ncol=2)
dev.off()

# pdf("Binary_diff.pdf",width=10,height=5)
# grid.arrange(q1,p1,ncol=2)
# dev.off()
# library(survival)
# km_fit_0 <- survfit(Surv(Y[Z==0], DELTA[Z==0]) ~ 1, data=dat)
# km_fit_1 <- survfit(Surv(Y[Z==1], DELTA[Z==1]) ~ 1, data=dat)
# data_0 = data.frame(cbind(km_fit_0$time,km_fit_0$surv,km_fit_0$std.err))
# data_1 = data.frame(cbind(km_fit_1$time,km_fit_1$surv,km_fit_1$std.err))
# colnames(data_0)=c("time","UW.SURV.1","UW.SURV.SE.1")
# colnames(data_1)=c("time","UW.SURV.2","UW.SURV.SE.2")
# 
# s0 = ggplot(aes(time),data = data_0) +
#   geom_ribbon(aes(ymin = pmax(0,UW.SURV.1 - qnorm(0.975)*UW.SURV.SE.1),
#                   ymax = pmin(1,UW.SURV.1 + qnorm(0.975)*UW.SURV.SE.1)),    # shadIPWing cnf intervals
#               fill = "magenta",alpha=0.2) +
#   geom_line(aes(x=time, y=UW.SURV.1,color="RP"),size = 1)+
#   geom_ribbon(data =data_1, aes(x=time,ymin = pmax(0,UW.SURV.2 - qnorm(0.975)*UW.SURV.SE.2),
#                   ymax = pmin(1,UW.SURV.2 + qnorm(0.975)*UW.SURV.SE.2)),    # shadIPWing cnf intervals
#               fill = "dodgerblue",alpha=0.2) +
#   geom_line(aes(time, UW.SURV.2,color = "EBRT+AD"),data=data_1,size = 1) +
#   xlab("Months after treatment") +
#   scale_color_manual(values = c(
#     'RP' = 'purple','EBRT+AD' = 'darkblue'))+labs(color = '')+
#   ylab("Survival Prob, Unweighted")+theme_bw()+
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   theme(legend.position = c(0.2, 0.15))
# 
# data = cbind(time.grid,ow.surv.mu,
#              ow.surv.mu.se,
#              ipw.surv.mu,
#              ipw.surv.mu.se)
# colnames(data) = c("time","OW.SURV.1","OW.SURV.2","OW.SURV.3",
#                    "OW.SURV.SE.1","OW.SURV.SE.2","OW.SURV.SE.3",
#                    "IPW.SURV.1","IPW.SURV.2","IPW.SURV.3",
#                    "IPW.SURV.SE.1","IPW.SURV.SE.2","IPW.SURV.SE.3")
# data =  as.data.frame(data)
# s1 = ggplot(aes(time),data = data) +
#   geom_ribbon(aes(ymin = IPW.SURV.1 - qnorm(0.975)*IPW.SURV.SE.1,
#                   ymax = IPW.SURV.1 + qnorm(0.975)*IPW.SURV.SE.1),    # shadIPWing cnf intervals
#               fill = "magenta",alpha=0.2) +
#   geom_line(aes(x=time, y=IPW.SURV.1,color="RP"),size = 1)   +
#   geom_ribbon(aes(ymin = IPW.SURV.2 - qnorm(0.975)*IPW.SURV.SE.2,
#                   ymax = IPW.SURV.2 + qnorm(0.975)*IPW.SURV.SE.2),    # shadIPWing cnf intervals
#               fill = "dodgerblue",alpha=0.2) +
#   geom_line(aes(time, IPW.SURV.2,color = "EBRT+AD"),size = 1) +
#   xlab("Months after treatment") +
#   scale_color_manual(values = c(
#     'RP' = 'purple','EBRT+AD' = 'darkblue'))+labs(color = '')+
#   ylab("Survival Prob, adjusted by IPW")+theme_bw()+ylim(c(0,1))+
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   theme(legend.position = "none")
# 
# s2 = ggplot(aes(time),data = data) +
#   geom_ribbon(aes(ymin = OW.SURV.1 - qnorm(0.975)*OW.SURV.SE.1,
#                   ymax = OW.SURV.1 + qnorm(0.975)*OW.SURV.SE.1),    # shadowing cnf intervals
#               fill = "magenta",alpha=0.2) +
#   geom_line(aes(x=time, y=OW.SURV.1,color="RP"),size = 1)   +
#   geom_ribbon(aes(ymin = OW.SURV.2 - qnorm(0.975)*OW.SURV.SE.2,
#                   ymax = OW.SURV.2 + qnorm(0.975)*OW.SURV.SE.2),    # shadowing cnf intervals
#               fill = "dodgerblue",alpha=0.2) +
#   geom_line(aes(time, OW.SURV.2,color = "EBRT+AD"),size = 1) +
#   xlab("Months after treatment") +
#   scale_color_manual(values = c(
#     'RP' = 'purple','EBRT+AD' = 'darkblue'))+labs(color = '')+
#   ylab("Survival Prob, adjusted by OW")+theme_bw()+ylim(c(0,1))+
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank()) +
#   theme(legend.position = "none")
# pdf("Binary_surv.pdf",width=12,height=4)
# grid.arrange(s0,s1,s2,ncol=3)
# dev.off()
