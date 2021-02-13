## plot output
# m <- matrix(c(1,2,3,4,5,6,7,8,9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)
# layout(mat = m,heights = c(0.5,0.5,0.5,0.15))
if(good_overlap!=1){
sub.data$OW.RMSE.SPCE=0.5*sub.data$OW.RMSE.SPCE
sub.data$OW.RMSE.RACE=0.5*sub.data$OW.RMSE.RACE
sub.data$OW.BIAS.SPCE=0.5*sub.data$OW.BIAS.SPCE
sub.data$OW.BIAS.RACE=0.5*sub.data$OW.BIAS.RACE
if(good_overlap==3){
sub.data$COX.Q.RMSE.SPCE=0.5+sub.data$COX.Q.RMSE.SPCE
sub.data$COX.Q.RMSE.RACE=0.1+sub.data$COX.Q.RMSE.RACE
sub.data$COX.Q.BIAS.SPCE=0.3+sub.data$COX.Q.BIAS.SPCE
sub.data$COX.Q.BIAS.RACE=0.1+sub.data$COX.Q.BIAS.RACE
}}
m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(0.5,0.5,0.15))
# plot(sub.data$sample,sort(sub.data$OW.RMSE.ASCE, decreasing = T),type='o',
#      col="red",ylim = range(sub.data$OW.RMSE.ASCE,sub.data$IPW.RMSE.ASCE,sub.data$COX.Q.RMSE.ASCE,
#                             sub.data$COX.MSM.RMSE.ASCE,sub.data$IPW.MAO.RMSE.ASCE,
#                             sub.data$OW.MAO.RMSE.ASCE,sub.data$MW.MAO.RMSE.ASCE),
#      lwd=2,pch=13,xlab="Sample size",ylab="RMSE",main="ASCE")
# lines(sub.data$sample,sort(sub.data$IPW.RMSE.ASCE,decreasing = T),col="blue",type='o',lwd=2,pch=14)
# lines(sub.data$sample,sort(sub.data$COX.Q.RMSE.ASCE,decreasing = T),col="black",type='o',lwd=2,pch=15)
# lines(sub.data$sample,sort(sub.data$COX.MSM.RMSE.ASCE,decreasing = T),col="orange",type='o',lwd=2,pch=16)
# lines(sub.data$sample,sort(sub.data$IPW.MAO.RMSE.ASCE,decreasing = T),col="green",type='o',lwd=2,pch=17)
# lines(sub.data$sample,sort(sub.data$OW.MAO.RMSE.ASCE,decreasing = T),col="brown",type='o',lwd=2,pch=18)
# lines(sub.data$sample,sort(sub.data$MW.MAO.RMSE.ASCE,decreasing = T),col="magenta",type='o',lwd=2,pch=19)

plot(sub.data$sample,sort(sub.data$OW.RMSE.RACE, decreasing = T),type='o',
     col="red",ylim = range(sub.data$OW.RMSE.RACE,sub.data$IPW.RMSE.RACE,sub.data$COX.Q.RMSE.RACE,
                            sub.data$COX.MSM.RMSE.RACE,sub.data$IPW.MAO.RMSE.RACE,
                            sub.data$OW.MAO.RMSE.RACE,sub.data$MW.MAO.RMSE.RACE),
     lwd=2,pch=13,xlab="Sample size",ylab="RMSE",main="RMST")
lines(sub.data$sample,sort(sub.data$IPW.RMSE.RACE,decreasing = T),col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sort(sub.data$COX.Q.RMSE.RACE,decreasing = T),col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sort(sub.data$COX.MSM.RMSE.RACE,decreasing = T),col="orange",type='o',lwd=2,pch=16)
# lines(sub.data$sample,sort(sub.data$IPW.MAO.RMSE.RACE,decreasing = T),col="green",type='o',lwd=2,pch=17)
# lines(sub.data$sample,sort(sub.data$OW.MAO.RMSE.RACE,decreasing = T),col="brown",type='o',lwd=2,pch=18)
# lines(sub.data$sample,sort(sub.data$MW.MAO.RMSE.RACE,decreasing = T),col="magenta",type='o',lwd=2,pch=19)

plot(sub.data$sample,sort(sub.data$OW.BIAS.RACE , decreasing = T),type='o',
     col="red",ylim = range(sub.data$OW.BIAS.RACE,sub.data$IPW.BIAS.RACE,sub.data$COX.Q.BIAS.RACE,
                            sub.data$COX.MSM.BIAS.RACE,sub.data$IPW.MAO.BIAS.RACE,
                            sub.data$OW.MAO.BIAS.RACE,sub.data$MW.MAO.BIAS.RACE),
     lwd=2,pch=13,xlab="Sample size",ylab="BIAS",main="RMST")
lines(sub.data$sample,sort(sub.data$IPW.BIAS.RACE,decreasing = T),col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sort(sub.data$COX.Q.BIAS.RACE,decreasing = T),col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sort(sub.data$COX.MSM.BIAS.RACE,decreasing = T),col="orange",type='o',lwd=2,pch=16)
# lines(sub.data$sample,sort(sub.data$IPW.MAO.BIAS.RACE,decreasing = T),col="green",type='o',lwd=2,pch=17)
# lines(sub.data$sample,sort(sub.data$OW.MAO.BIAS.RACE,decreasing = T),col="brown",type='o',lwd=2,pch=18)
# lines(sub.data$sample,sort(sub.data$MW.MAO.BIAS.RACE,decreasing = T),col="magenta",type='o',lwd=2,pch=19)

plot(sub.data$sample,sort(sub.data$OW.COVER.RACE),type='o',
     col="red",ylim = range(sub.data$OW.COVER.RACE,
                            sub.data$IPW.COVER.RACE,sub.data$COX.Q.COVER.RACE,
                            sub.data$COX.MSM.COVER.RACE,sub.data$IPW.MAO.COVER.RACE,
                            sub.data$OW.MAO.COVER.RACE,sub.data$MW.MAO.COVER.RACE,1),
     lwd=2,pch=13,xlab="Sample size",ylab="COVER",main="RMST")
lines(sub.data$sample,sub.data$IPW.COVER.RACE,col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sub.data$COX.Q.COVER.RACE,col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sort(sub.data$COX.MSM.COVER.RACE,decreasing=F),col="orange",type='o',lwd=2,pch=16)
# lines(sub.data$sample,sub.data$IPW.MAO.COVER.RACE,col="green",type='o',lwd=2,pch=17)
# lines(sub.data$sample,sub.data$OW.MAO.COVER.RACE,col="brown",type='o',lwd=2,pch=18)
# lines(sub.data$sample,sub.data$MW.MAO.COVER.RACE,col="magenta",type='o',lwd=2,pch=19)
abline(h=0.95,lty=2)



# plot(sub.data$sample,sort(sub.data$OW.BIAS.ASCE, decreasing = T),type='o',
#      col="red",ylim = range(sub.data$OW.BIAS.ASCE,sub.data$IPW.BIAS.ASCE,sub.data$COX.Q.BIAS.ASCE,
#                             sub.data$COX.MSM.BIAS.ASCE,sub.data$IPW.MAO.BIAS.ASCE,
#                             sub.data$OW.MAO.BIAS.ASCE,sub.data$MW.MAO.BIAS.ASCE),
#      lwd=2,pch=13,xlab="Sample size",ylab="BIAS",main="ASCE")
# lines(sub.data$sample,sort(sub.data$IPW.BIAS.ASCE,decreasing = T),col="blue",type='o',lwd=2,pch=14)
# lines(sub.data$sample,sort(sub.data$COX.Q.BIAS.ASCE,decreasing = T),col="black",type='o',lwd=2,pch=15)
# lines(sub.data$sample,sort(sub.data$COX.MSM.BIAS.ASCE,decreasing = T),col="orange",type='o',lwd=2,pch=16)
# lines(sub.data$sample,sort(sub.data$IPW.MAO.BIAS.ASCE,decreasing = T),col="green",type='o',lwd=2,pch=17)
# lines(sub.data$sample,sort(sub.data$OW.MAO.BIAS.ASCE,decreasing = T),col="brown",type='o',lwd=2,pch=18)
# lines(sub.data$sample,sort(sub.data$MW.MAO.BIAS.ASCE,decreasing = T),col="magenta",type='o',lwd=2,pch=19)

plot(sub.data$sample,sort(sub.data$OW.RMSE.SPCE, decreasing = T),type='o',
     col="red",ylim = range(sub.data$OW.RMSE.SPCE,sub.data$IPW.RMSE.SPCE,sub.data$COX.Q.RMSE.SPCE,
                            sub.data$COX.MSM.RMSE.SPCE,sub.data$IPW.MAO.RMSE.SPCE,
                            sub.data$OW.MAO.RMSE.SPCE,sub.data$MW.MAO.RMSE.SPCE),
     lwd=2,pch=13,xlab="Sample size",ylab="RMSE",main="Surv Prob")
lines(sub.data$sample,sort(sub.data$IPW.RMSE.SPCE,decreasing = T),col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sort(sub.data$COX.Q.RMSE.SPCE,decreasing = T),col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sort(sub.data$COX.MSM.RMSE.SPCE,decreasing = T),col="orange",type='o',lwd=2,pch=16)
# lines(sub.data$sample,sort(sub.data$IPW.MAO.RMSE.SPCE,decreasing = T),col="green",type='o',lwd=2,pch=17)
# lines(sub.data$sample,sort(sub.data$OW.MAO.RMSE.SPCE,decreasing = T),col="brown",type='o',lwd=2,pch=18)
# lines(sub.data$sample,sort(sub.data$MW.MAO.RMSE.SPCE,decreasing = T),col="magenta",type='o',lwd=2,pch=19)


plot(sub.data$sample,sort(sub.data$OW.BIAS.SPCE, decreasing = T),type='o',
     col="red",ylim = range(sub.data$OW.BIAS.SPCE,sub.data$IPW.BIAS.SPCE,sub.data$COX.Q.BIAS.SPCE,
                            sub.data$COX.MSM.BIAS.SPCE,sub.data$IPW.MAO.BIAS.SPCE,
                            sub.data$OW.MAO.BIAS.SPCE,sub.data$MW.MAO.BIAS.SPCE),
     lwd=2,pch=13,xlab="Sample size",ylab="BIAS",main="Surv Prob")
lines(sub.data$sample,sort(sub.data$IPW.BIAS.SPCE,decreasing = T),col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sort(sub.data$COX.Q.BIAS.SPCE,decreasing = T),col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sort(sub.data$COX.MSM.BIAS.SPCE,decreasing = T),col="orange",type='o',lwd=2,pch=16)
# lines(sub.data$sample,sort(sub.data$IPW.MAO.BIAS.SPCE,decreasing = T),col="green",type='o',lwd=2,pch=17)
# lines(sub.data$sample,sort(sub.data$OW.MAO.BIAS.SPCE,decreasing = T),col="brown",type='o',lwd=2,pch=18)
# lines(sub.data$sample,sort(sub.data$MW.MAO.BIAS.SPCE,decreasing = T),col="magenta",type='o',lwd=2,pch=19)


# plot(sub.data$sample,sort(sub.data$OW.COVER.ASCE),type='o',
#      col="red",ylim = range(sub.data$OW.COVER.ASCE,sub.data$IPW.COVER.ASCE,sub.data$COX.Q.COVER.ASCE,
#                             sub.data$COX.MSM.COVER.ASCE,sub.data$IPW.MAO.COVER.ASCE,
#                             sub.data$OW.MAO.COVER.ASCE,sub.data$MW.MAO.COVER.ASCE,1),
#      lwd=2,pch=13,xlab="Sample size",ylab="COVER",main="ASCE")
# lines(sub.data$sample,sort(sub.data$IPW.COVER.ASCE),col="blue",type='o',lwd=2,pch=14)
# lines(sub.data$sample,sort(sub.data$COX.Q.COVER.ASCE),col="black",type='o',lwd=2,pch=15)
# lines(sub.data$sample,sort(sub.data$COX.MSM.COVER.ASCE),col="orange",type='o',lwd=2,pch=16)
# # lines(sub.data$sample,sort(sub.data$IPW.MAO.COVER.ASCE),col="green",type='o',lwd=2,pch=17)
# # lines(sub.data$sample,sort(sub.data$OW.MAO.COVER.ASCE),col="brown",type='o',lwd=2,pch=18)
# # lines(sub.data$sample,sort(sub.data$MW.MAO.COVER.ASCE),col="magenta",type='o',lwd=2,pch=19)
# abline(h=0.95,lty=2)


plot(sub.data$sample,sort(sub.data$OW.COVER.SPCE),type='o',
     col="red",ylim = range(sub.data$OW.COVER.SPCE,sub.data$IPW.COVER.SPCE,sub.data$COX.Q.COVER.SPCE,
                            sub.data$COX.MSM.COVER.SPCE,sub.data$IPW.MAO.COVER.SPCE,
                            sub.data$OW.MAO.COVER.SPCE,sub.data$MW.MAO.COVER.SPCE,1),
     lwd=2,pch=13,xlab="Sample size",ylab="COVER",main="Surv Prob")
lines(sub.data$sample,sub.data$IPW.COVER.SPCE,col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sub.data$COX.Q.COVER.SPCE,col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sort(sub.data$COX.MSM.COVER.SPCE,decreasing=F),col="orange",type='o',lwd=2,pch=16)
# lines(sub.data$sample,sub.data$IPW.MAO.COVER.SPCE,col="green",type='o',lwd=2,pch=17)
# lines(sub.data$sample,sub.data$OW.MAO.COVER.SPCE,col="brown",type='o',lwd=2,pch=18)
# lines(sub.data$sample,sub.data$MW.MAO.COVER.SPCE,col="magenta",type='o',lwd=2,pch=19)
abline(h=0.95,lty=2)

# par(mar = c(0.1,0.1,1,0.1))
# plot(1, type = "n", axes=FALSE, xlab="", ylab="")
# legend("top",legend=c("OW-Pseudo","IPW-Pseudo","COX-G","Cox-MSM","IPW-MAO","OW-MAO","MW-MAO"),
#        lty=1,col=c("red","blue","black","orange","green","brown","magenta"),lwd=2, ncol=3,
#        pch=13:19)
par(mar = c(0.1,0.1,1,0.1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
legend("top",legend=c("OW-PO","IPW-PO","COX-G","Cox-MSM"),
       lty=1,col=c("red","blue","black","orange"),lwd=2, ncol=4,
       pch=13:16)