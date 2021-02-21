## plot output
if(aipw.draw){
  m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)
  layout(mat = m,heights = c(0.5,0.5,0.15))
}else{
  m <- matrix(c(1,2,3,4,5,6,7,8,9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)
  layout(mat = m,heights = c(0.5,0.5,0.5,0.15))
}



if(good_overlap<=2){
  ylim.range = range(sub.data$OW.BIAS.SPCE1,sub.data$IPW.UNTRIM.BIAS.SPCE1,sub.data$COX.Q.BIAS.SPCE1,
                     sub.data$COX.MSM.BIAS.SPCE1,sub.data$IPW.AIPW.BIAS.SPCE1,
                     sub.data$OW.AIPW.BIAS.SPCE1)
}else{
  ylim.range = range(sub.data$OW.BIAS.SPCE1,sub.data$IPW.UNTRIM.BIAS.SPCE1,sub.data$COX.Q.BIAS.SPCE1,
                     sub.data$IPW.AIPW.BIAS.SPCE1,
                     sub.data$OW.AIPW.BIAS.SPCE1)
}
plot(sub.data$sample,sub.data$OW.BIAS.SPCE1,type='o',
     col="red",ylim = ylim.range,
     lwd=2,pch=13,xlab="Sample size",ylab="BIAS",main="SPCE")
lines(sub.data$sample,sub.data$IPW.UNTRIM.BIAS.SPCE1,col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sub.data$COX.Q.BIAS.SPCE1,col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sub.data$COX.MSM.BIAS.SPCE1,col="orange",type='o',lwd=2,pch=16)

if(aipw.draw){
  lines(sub.data$sample,sub.data$OW.AIPW.BIAS.SPCE1,col="green",type='o',lwd=2,pch=17)
  lines(sub.data$sample,sub.data$IPW.AIPW.BIAS.SPCE1,col="magenta",type='o',lwd=2,pch=19)
}
if(pseudo.g.draw){
  lines(sub.data$sample,sub.data$PSEUDO.G.BIAS.SPCE1,col="green",type='o',lwd=2,pch=20)
  lines(sub.data$sample,sub.data$PSEUDO.UNADJ.BIAS.SPCE1,col="magenta",type='o',lwd=2,pch=21)
}
if(mao.draw){
  lines(sub.data$sample,sub.data$IPW.MAO.BIAS.SPCE1,col="green",type='o',lwd=2,pch=18)
  lines(sub.data$sample,sub.data$OW.MAO.BIAS.SPCE1,col="magenta",type='o',lwd=2,pch=22)
}


if(good_overlap<=2){
  ylim.range = range(sub.data$OW.BIAS.RACE1,sub.data$IPW.UNTRIM.BIAS.RACE1,sub.data$COX.Q.BIAS.RACE1,
                     sub.data$COX.MSM.BIAS.RACE1,
                     sub.data$IPW.AIPW.BIAS.RACE1,
                     sub.data$OW.AIPW.BIAS.RACE1)
}else{
  ylim.range = range(sub.data$OW.BIAS.RACE1,sub.data$IPW.UNTRIM.BIAS.RACE1,sub.data$COX.Q.BIAS.RACE1,
                     sub.data$IPW.AIPW.BIAS.RACE1,
                     sub.data$OW.AIPW.BIAS.RACE1)
}
plot(sub.data$sample,sub.data$OW.BIAS.RACE1,type='o',
     col="red",ylim = ylim.range,
     lwd=2,pch=13,xlab="Sample size",ylab="BIAS",main="RACE")
lines(sub.data$sample,sub.data$IPW.UNTRIM.BIAS.RACE1,col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sub.data$COX.Q.BIAS.RACE1,col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sub.data$COX.MSM.BIAS.RACE1,col="orange",type='o',lwd=2,pch=16)
if(aipw.draw){
  lines(sub.data$sample,sub.data$OW.AIPW.BIAS.RACE1,col="green",type='o',lwd=2,pch=17)
  lines(sub.data$sample,sub.data$IPW.AIPW.BIAS.RACE1,col="magenta",type='o',lwd=2,pch=19)
}
if(pseudo.g.draw){
  lines(sub.data$sample,sub.data$PSEUDO.G.BIAS.RACE1,col="green",type='o',lwd=2,pch=20)
  lines(sub.data$sample,sub.data$PSEUDO.UNADJ.BIAS.RACE1,col="magenta",type='o',lwd=2,pch=21)
}
if(mao.draw){
  lines(sub.data$sample,sub.data$IPW.MAO.BIAS.RACE1,col="green",type='o',lwd=2,pch=18)
  lines(sub.data$sample,sub.data$OW.MAO.BIAS.RACE1,col="magenta",type='o',lwd=2,pch=22)
}

if(good_overlap<=2){
  ylim.range = range(sub.data$OW.BIAS.ASCE1,sub.data$IPW.UNTRIM.BIAS.ASCE1,sub.data$COX.Q.BIAS.ASCE1,
                     sub.data$COX.MSM.BIAS.ASCE1,
                     sub.data$IPW.AIPW.BIAS.ASCE1,
                     sub.data$OW.AIPW.BIAS.ASCE1)
}else{
  ylim.range = range(sub.data$OW.BIAS.ASCE1,sub.data$IPW.UNTRIM.BIAS.ASCE1,sub.data$COX.Q.BIAS.ASCE1,
                     sub.data$IPW.AIPW.BIAS.ASCE1,
                     sub.data$OW.AIPW.BIAS.ASCE1)
}
plot(sub.data$sample,sub.data$OW.BIAS.ASCE1,type='o',
     col="red",ylim = ylim.range,
     lwd=2,pch=13,xlab="Sample size",ylab="BIAS",main="ASCE")
lines(sub.data$sample,sub.data$IPW.UNTRIM.BIAS.ASCE1,col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sub.data$COX.Q.BIAS.ASCE1,col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sub.data$COX.MSM.BIAS.ASCE1,col="orange",type='o',lwd=2,pch=16)
if(aipw.draw){
  lines(sub.data$sample,sub.data$OW.AIPW.BIAS.ASCE1,col="green",type='o',lwd=2,pch=17)
  lines(sub.data$sample,sub.data$IPW.AIPW.BIAS.ASCE1,col="magenta",type='o',lwd=2,pch=19)
}
if(pseudo.g.draw){
  lines(sub.data$sample,sub.data$PSEUDO.G.BIAS.ASCE1,col="green",type='o',lwd=2,pch=20)
  lines(sub.data$sample,sub.data$PSEUDO.UNADJ.BIAS.ASCE1,col="magenta",type='o',lwd=2,pch=21)
}
if(mao.draw){
  lines(sub.data$sample,sub.data$IPW.MAO.BIAS.ASCE1,col="green",type='o',lwd=2,pch=18)
  lines(sub.data$sample,sub.data$OW.MAO.BIAS.ASCE1,col="magenta",type='o',lwd=2,pch=22)
}



plot(sub.data$sample,sub.data$OW.RMSE.SPCE1,type='o',
     col="red",ylim = range(sub.data$OW.RMSE.SPCE1,sub.data$IPW.UNTRIM.RMSE.SPCE1,sub.data$COX.Q.RMSE.SPCE1,
                            sub.data$COX.MSM.RMSE.SPCE1,
                            sub.data$IPW.AIPW.RMSE.SPCE1,
                            sub.data$OW.AIPW.RMSE.SPCE1),
     lwd=2,pch=13,xlab="Sample size",ylab="RMSE",main="SPCE")
lines(sub.data$sample,sub.data$IPW.RMSE.SPCE1,col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sub.data$COX.Q.RMSE.SPCE1,col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sub.data$COX.MSM.RMSE.SPCE1,col="orange",type='o',lwd=2,pch=16)
if(aipw.draw){
  lines(sub.data$sample,sub.data$OW.AIPW.RMSE.SPCE1,col="green",type='o',lwd=2,pch=17)
  lines(sub.data$sample,sub.data$IPW.AIPW.RMSE.SPCE1,col="magenta",type='o',lwd=2,pch=19)
}
if(pseudo.g.draw){
  lines(sub.data$sample,sub.data$PSEUDO.G.RMSE.SPCE1,col="green",type='o',lwd=2,pch=20)
  lines(sub.data$sample,sub.data$PSEUDO.UNADJ.RMSE.SPCE1,col="magenta",type='o',lwd=2,pch=21)
}
if(mao.draw){
  lines(sub.data$sample,sub.data$IPW.MAO.RMSE.SPCE1,col="green",type='o',lwd=2,pch=18)
  lines(sub.data$sample,sub.data$OW.MAO.RMSE.SPCE1,col="magenta",type='o',lwd=2,pch=22)
}


plot(sub.data$sample,sub.data$OW.RMSE.RACE1,type='o',
     col="red",ylim = range(sub.data$OW.RMSE.RACE1,sub.data$IPW.UNTRIM.RMSE.RACE1,sub.data$COX.Q.RMSE.RACE1,
                            sub.data$COX.MSM.RMSE.RACE1,
                            sub.data$IPW.AIPW.RMSE.RACE1,
                            sub.data$OW.AIPW.RMSE.RACE1),
     lwd=2,pch=13,xlab="Sample size",ylab="RMSE",main="RACE")
lines(sub.data$sample,sub.data$IPW.UNTRIM.RMSE.RACE1,col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sub.data$COX.Q.RMSE.RACE1,col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sub.data$COX.MSM.RMSE.RACE1,col="orange",type='o',lwd=2,pch=16)
if(aipw.draw){
  lines(sub.data$sample,sub.data$OW.AIPW.RMSE.RACE1,col="green",type='o',lwd=2,pch=17)
  lines(sub.data$sample,sub.data$IPW.AIPW.RMSE.RACE1,col="magenta",type='o',lwd=2,pch=19)
}
if(pseudo.g.draw){
  lines(sub.data$sample,sub.data$PSEUDO.G.RMSE.RACE1,col="green",type='o',lwd=2,pch=20)
  lines(sub.data$sample,sub.data$PSEUDO.UNADJ.RMSE.RACE1,col="magenta",type='o',lwd=2,pch=21)
}
if(mao.draw){
  lines(sub.data$sample,sub.data$IPW.MAO.RMSE.RACE1,col="green",type='o',lwd=2,pch=18)
  lines(sub.data$sample,sub.data$OW.MAO.RMSE.RACE1,col="magenta",type='o',lwd=2,pch=22)
}



plot(sub.data$sample,sub.data$OW.RMSE.ASCE1,type='o',
     col="red",ylim = range(sub.data$OW.RMSE.ASCE1,sub.data$IPW.RMSE.ASCE1,sub.data$COX.Q.RMSE.ASCE1,
                            sub.data$COX.MSM.RMSE.ASCE1,
                            sub.data$IPW.AIPW.RMSE.ASCE1,
                            sub.data$OW.AIPW.RMSE.ASCE1),
     lwd=2,pch=13,xlab="Sample size",ylab="RMSE",main="ASCE")
lines(sub.data$sample,sub.data$IPW.UNTRIM.RMSE.ASCE1,col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sub.data$COX.Q.RMSE.ASCE1,col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sub.data$COX.MSM.RMSE.ASCE1,col="orange",type='o',lwd=2,pch=15)
if(aipw.draw){
  lines(sub.data$sample,sub.data$OW.AIPW.RMSE.ASCE1,col="green",type='o',lwd=2,pch=17)
  lines(sub.data$sample,sub.data$IPW.AIPW.RMSE.ASCE1,col="magenta",type='o',lwd=2,pch=19)
}
if(pseudo.g.draw){
  lines(sub.data$sample,sub.data$PSEUDO.G.RMSE.ASCE1,col="green",type='o',lwd=2,pch=20)
  lines(sub.data$sample,sub.data$PSEUDO.UNADJ.RMSE.ASCE1,col="magenta",type='o',lwd=2,pch=21)
}
if(mao.draw){
  lines(sub.data$sample,sub.data$IPW.MAO.RMSE.ASCE1,col="green",type='o',lwd=2,pch=18)
  lines(sub.data$sample,sub.data$OW.MAO.RMSE.ASCE1,col="magenta",type='o',lwd=2,pch=22)
}

if(!aipw.draw){
  plot(sub.data$sample,sub.data$OW.COVER.SPCE1,type='o',
       col="red",ylim = range(sub.data$OW.COVER.SPCE1,sub.data$IPW.UNTRIM.COVER.SPCE1,sub.data$COX.Q.COVER.SPCE1,
                              sub.data$COX.MSM.COVER.SPCE1,1),
       lwd=2,pch=13,xlab="Sample size",ylab="COVER",main="SPCE")
  lines(sub.data$sample,sub.data$IPW.UNTRIM.COVER.SPCE1,col="blue",type='o',lwd=2,pch=14)
  lines(sub.data$sample,sub.data$COX.Q.COVER.SPCE1,col="black",type='o',lwd=2,pch=15)
  lines(sub.data$sample,sub.data$COX.MSM.COVER.SPCE1,col="orange",type='o',lwd=2,pch=16)
  if(pseudo.g.draw){
    lines(sub.data$sample,sub.data$PSEUDO.G.COVER.SPCE1,col="green",type='o',lwd=2,pch=20)
    lines(sub.data$sample,sub.data$PSEUDO.UNADJ.COVER.SPCE1,col="magenta",type='o',lwd=2,pch=21)
  }
  if(mao.draw){
    lines(sub.data$sample,sub.data$IPW.MAO.COVER.SPCE1,col="green",type='o',lwd=2,pch=18)
    lines(sub.data$sample,sub.data$OW.MAO.COVER.SPCE1,col="magenta",type='o',lwd=2,pch=22)
  }
  abline(h=0.95,lty=2)
  
  
  plot(sub.data$sample,sub.data$OW.COVER.RACE1,type='o',
       col="red",ylim = range(sub.data$OW.COVER.RACE1,
                              sub.data$IPW.UNTRIM.COVER.RACE1,sub.data$COX.Q.COVER.RACE1,
                              sub.data$COX.MSM.COVER.RACE1,1),
       lwd=2,pch=13,xlab="Sample size",ylab="COVER",main="RACE")
  lines(sub.data$sample,sub.data$IPW.UNTRIM.COVER.RACE1,col="blue",type='o',lwd=2,pch=14)
  lines(sub.data$sample,sub.data$COX.Q.COVER.RACE1,col="black",type='o',lwd=2,pch=15)
  lines(sub.data$sample,sub.data$COX.MSM.COVER.RACE1,col="orange",type='o',lwd=2,pch=16)
  if(pseudo.g.draw){
    lines(sub.data$sample,sub.data$PSEUDO.G.COVER.RACE1,col="green",type='o',lwd=2,pch=20)
    lines(sub.data$sample,sub.data$PSEUDO.UNADJ.COVER.RACE1,col="magenta",type='o',lwd=2,pch=21)
  }
  if(mao.draw){
    lines(sub.data$sample,sub.data$IPW.MAO.COVER.RACE1,col="green",type='o',lwd=2,pch=18)
    lines(sub.data$sample,sub.data$OW.MAO.COVER.RACE1,col="magenta",type='o',lwd=2,pch=22)
  }
  abline(h=0.95,lty=2)
  
  if(!mao.draw){
    ylim.range = range(sub.data$OW.COVER.ASCE1,sub.data$IPW.UNTRIM.COVER.ASCE1,sub.data$COX.Q.COVER.ASCE1,
                       sub.data$COX.MSM.COVER.ASCE1,1)
  }else{
    ylim.range = range(sub.data$OW.COVER.ASCE1,sub.data$IPW.UNTRIM.COVER.ASCE1,sub.data$COX.Q.COVER.ASCE1,
                       sub.data$COX.MSM.COVER.ASCE1,sub.data$IPW.MAO.COVER.ASCE1,sub.data$OW.MAO.COVER.ASCE11)
  }
  
  plot(sub.data$sample,sub.data$OW.COVER.ASCE1,type='o',
       col="red",ylim = ylim.range,
       lwd=2,pch=13,xlab="Sample size",ylab="COVER",main="ASCE")
  lines(sub.data$sample,sub.data$IPW.UNTRIM.COVER.ASCE1,col="blue",type='o',lwd=2,pch=14)
  lines(sub.data$sample,sub.data$COX.Q.COVER.ASCE1,col="black",type='o',lwd=2,pch=15)
  lines(sub.data$sample,sub.data$COX.MSM.COVER.ASCE1,col="orange",type='o',lwd=2,pch=16)
  if(pseudo.g.draw){
    lines(sub.data$sample,sub.data$PSEUDO.G.COVER.ASCE1,col="green",type='o',lwd=2,pch=20)
    lines(sub.data$sample,sub.data$PSEUDO.UNADJ.COVER.ASCE1,col="magenta",type='o',lwd=2,pch=21)
  }
  if(mao.draw){
    lines(sub.data$sample,sub.data$IPW.MAO.COVER.ASCE1,col="green",type='o',lwd=2,pch=18)
    lines(sub.data$sample,sub.data$OW.MAO.COVER.ASCE1,col="magenta",type='o',lwd=2,pch=22)
  }
  abline(h=0.95,lty=2)
}

par(mar = c(0.1,0.1,1,0.1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
if(aipw.draw){
  legend("top",legend=c("OW","IPW","Cox-G","Cox-IPW","AOW","AIPW"),
         lty=1,col=c("red","blue","black","orange","green","magenta"),lwd=2, horiz = T,
         pch=c(13,14,15,16,17,19))
}else if(pseudo.g.draw){
  legend("top",legend=c("OW","IPW","Cox-G","Cox-IPW","PO-G","PO-UNADJ"),
         lty=1,col=c("red","blue","black","orange","green","magenta"),lwd=2, horiz = T,
         pch=c(13,14,15,16,20,21))
}else if(mao.draw){
  legend("top",legend=c("OW","IPW","Cox-G","Cox-IPW","IPW-MAO","OW-MAO"),
         lty=1,col=c("red","blue","black","orange","green","magenta"),lwd=2, horiz = T,
         pch=c(13,14,15,16,18,22))
}else{
  legend("top",legend=c("OW","IPW","COX-G","Cox-IPW"),
         lty=1,col=c("red","blue","black","orange"),lwd=2, ncol=4,
         pch=c(13,14,15,16))
}

