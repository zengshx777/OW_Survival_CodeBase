# 
set.seed(1)
pdf("Binary_ps_SIMU.pdf",width=11,height=4)
m <- matrix(c(1,2,3,4,4,4),nrow = 2,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(0.5,0.1))
good_overlap=1; sample_size=600; dependent.censoring=F;i=1
source("simu_data_gen.R")
plot(density(soft.max[Z==0,2]+rnorm(sum(Z==0),0,0.0003)),xlim=c(0,1),col="red",lwd=2,
     xlab="True propensity score",main="Balanced")
lines(density(soft.max[Z==1,2]+rnorm(sum(Z==1),0,0.0003)),lwd=2,col="blue")

# 
good_overlap=2; sample_size=600; dependent.censoring=F;i=1
source("simu_data_gen.R")
plot(density(soft.max[Z==0,2]),xlim=c(0,1),col="red",lwd=2,ylim=c(0,4),
     xlab="True propensity score",main="Moderate imbalance")
lines(density(soft.max[Z==1,2]),lwd=2,col="blue")

good_overlap=3; sample_size=600; dependent.censoring=F;i=1
source("simu_data_gen.R")
plot(density(soft.max[Z==0,2]),xlim=c(0,1),col="red",lwd=2,ylim=c(0,4),
     xlab="True propensity score",main="Huge imbalance")
lines(density(soft.max[Z==1,2]),col="blue",lwd=2)

par(mar = c(0.1,0.1,1,0.1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top",legend=c("Z=0","Z=1"),
       lty=1,col=c("red","blue"),lwd=2, ncol=2)
dev.off()