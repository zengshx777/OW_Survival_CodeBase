## plot output
if(aipw.draw){
  m <- matrix(c(1,2,3,4,5,6,7,7,7),nrow = 3,ncol = 3,byrow = TRUE)
  layout(mat = m,heights = c(0.5,0.5,0.15))
}else{
  m <- matrix(c(1,2,3,4,5,6,7,8,9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)
  layout(mat = m,heights = c(0.5,0.5,0.5,0.15))
}



if(good_overlap<=2){
  ylim.range = range(sub.data$OW.BIAS.SPCE3,sub.data$IPW.BIAS.SPCE3,sub.data$COX.Q.BIAS.SPCE3,
                     sub.data$COX.MSM.BIAS.SPCE3,sub.data$IPW.AIPW.BIAS.SPCE3,
                     sub.data$OW.AIPW.BIAS.SPCE3)
}else{
  ylim.range = range(sub.data$OW.BIAS.SPCE3,sub.data$IPW.BIAS.SPCE3,sub.data$COX.Q.BIAS.SPCE3,
                     sub.data$IPW.AIPW.BIAS.SPCE3,
                     sub.data$OW.AIPW.BIAS.SPCE3)
}
plot(sub.data$sample,sub.data$OW.BIAS.SPCE3,type='o',
     col="red",ylim = ylim.range,
     lwd=2,pch=13,xlab="Sample size",ylab="BIAS",main="SPCE")
lines(sub.data$sample,sub.data$IPW.BIAS.SPCE3,col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sub.data$COX.Q.BIAS.SPCE3,col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sub.data$COX.MSM.BIAS.SPCE3,col="orange",type='o',lwd=2,pch=16)

if(aipw.draw){
  lines(sub.data$sample,sub.data$OW.AIPW.BIAS.SPCE3,col="green",type='o',lwd=2,pch=17)
  lines(sub.data$sample,sub.data$IPW.AIPW.BIAS.SPCE3,col="magenta",type='o',lwd=2,pch=19)
}
if(pseudo.g.draw){
  lines(sub.data$sample,sub.data$PSEUDO.G.BIAS.SPCE3,col="green",type='o',lwd=2,pch=20)
  if(good_overlap==1){
    lines(sub.data$sample,sub.data$PSEUDO.UNADJ.BIAS.SPCE3,col="magenta",type='o',lwd=2,pch=21)
  }
}
if(mao.draw){
  lines(sub.data$sample,sub.data$IPW.MAO.BIAS.SPCE3,col="green",type='o',lwd=2,pch=18)
  lines(sub.data$sample,sub.data$OW.MAO.BIAS.SPCE3,col="magenta",type='o',lwd=2,pch=22)
}


if(good_overlap<=2){
  ylim.range = range(sub.data$OW.BIAS.RACE3,sub.data$IPW.BIAS.RACE3,sub.data$COX.Q.BIAS.RACE3,
                     sub.data$COX.MSM.BIAS.RACE3,
                     sub.data$IPW.AIPW.BIAS.RACE3,
                     sub.data$OW.AIPW.BIAS.RACE3)
}else{
  ylim.range = range(sub.data$OW.BIAS.RACE3,sub.data$IPW.BIAS.RACE3,sub.data$COX.Q.BIAS.RACE3,
                     sub.data$IPW.AIPW.BIAS.RACE3,
                     sub.data$OW.AIPW.BIAS.RACE3)
}
plot(sub.data$sample,sub.data$OW.BIAS.RACE3,type='o',
     col="red",ylim = ylim.range,
     lwd=2,pch=13,xlab="Sample size",ylab="BIAS",main="RACE")
lines(sub.data$sample,sub.data$IPW.BIAS.RACE3,col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sub.data$COX.Q.BIAS.RACE3,col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sub.data$COX.MSM.BIAS.RACE3,col="orange",type='o',lwd=2,pch=16)
if(aipw.draw){
  lines(sub.data$sample,sub.data$OW.AIPW.BIAS.RACE3,col="green",type='o',lwd=2,pch=17)
  lines(sub.data$sample,sub.data$IPW.AIPW.BIAS.RACE3,col="magenta",type='o',lwd=2,pch=19)
}
if(pseudo.g.draw){
  lines(sub.data$sample,sub.data$PSEUDO.G.BIAS.RACE3,col="green",type='o',lwd=2,pch=20)
  if(good_overlap==1){
    lines(sub.data$sample,sub.data$PSEUDO.UNADJ.BIAS.RACE3,col="magenta",type='o',lwd=2,pch=21)
  }
}
if(mao.draw){
  lines(sub.data$sample,sub.data$IPW.MAO.BIAS.RACE3,col="green",type='o',lwd=2,pch=18)
  lines(sub.data$sample,sub.data$OW.MAO.BIAS.RACE3,col="magenta",type='o',lwd=2,pch=22)
}

if(good_overlap<=2){
  ylim.range = range(sub.data$OW.BIAS.ASCE3,sub.data$IPW.BIAS.ASCE3,sub.data$COX.Q.BIAS.ASCE3,
                     sub.data$COX.MSM.BIAS.ASCE3,
                     sub.data$IPW.AIPW.BIAS.ASCE3,
                     sub.data$OW.AIPW.BIAS.ASCE3)
}else{
  ylim.range = range(sub.data$OW.BIAS.ASCE3,sub.data$IPW.BIAS.ASCE3,sub.data$COX.Q.BIAS.ASCE3,
                     sub.data$IPW.AIPW.BIAS.ASCE3,
                     sub.data$OW.AIPW.BIAS.ASCE3)
}
plot(sub.data$sample,sub.data$OW.BIAS.ASCE3,type='o',
     col="red",ylim = ylim.range,
     lwd=2,pch=13,xlab="Sample size",ylab="BIAS",main="ASCE")
lines(sub.data$sample,sub.data$IPW.BIAS.ASCE3,col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sub.data$COX.Q.BIAS.ASCE3,col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sub.data$COX.MSM.BIAS.ASCE3,col="orange",type='o',lwd=2,pch=16)
if(aipw.draw){
  lines(sub.data$sample,sub.data$OW.AIPW.BIAS.ASCE3,col="green",type='o',lwd=2,pch=17)
  lines(sub.data$sample,sub.data$IPW.AIPW.BIAS.ASCE3,col="magenta",type='o',lwd=2,pch=19)
}
if(pseudo.g.draw){
  lines(sub.data$sample,sub.data$PSEUDO.G.BIAS.ASCE3,col="green",type='o',lwd=2,pch=20)
  if(good_overlap==1){
    lines(sub.data$sample,sub.data$PSEUDO.UNADJ.BIAS.ASCE3,col="magenta",type='o',lwd=2,pch=21)
  }
  
}
if(mao.draw){
  lines(sub.data$sample,sub.data$IPW.MAO.BIAS.ASCE3,col="green",type='o',lwd=2,pch=18)
  lines(sub.data$sample,sub.data$OW.MAO.BIAS.ASCE3,col="magenta",type='o',lwd=2,pch=22)
}



plot(sub.data$sample,sub.data$OW.RMSE.SPCE3,type='o',
     col="red",ylim = range(sub.data$OW.RMSE.SPCE3,sub.data$IPW.RMSE.SPCE3,sub.data$COX.Q.RMSE.SPCE3,
                            sub.data$COX.MSM.RMSE.SPCE3,
                            sub.data$IPW.AIPW.RMSE.SPCE3,
                            sub.data$OW.AIPW.RMSE.SPCE3),
     lwd=2,pch=13,xlab="Sample size",ylab="RMSE",main="SPCE")
lines(sub.data$sample,sub.data$IPW.RMSE.SPCE3,col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sub.data$COX.Q.RMSE.SPCE3,col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sub.data$COX.MSM.RMSE.SPCE3,col="orange",type='o',lwd=2,pch=16)
if(aipw.draw){
  lines(sub.data$sample,sub.data$OW.AIPW.RMSE.SPCE3,col="green",type='o',lwd=2,pch=17)
  lines(sub.data$sample,sub.data$IPW.AIPW.RMSE.SPCE3,col="magenta",type='o',lwd=2,pch=19)
}
if(pseudo.g.draw){
  lines(sub.data$sample,sub.data$PSEUDO.G.RMSE.SPCE3,col="green",type='o',lwd=2,pch=20)
  if(good_overlap==1){
    lines(sub.data$sample,sub.data$PSEUDO.UNADJ.RMSE.SPCE3,col="magenta",type='o',lwd=2,pch=21)
    
  }
}
if(mao.draw){
  lines(sub.data$sample,sub.data$IPW.MAO.RMSE.SPCE3,col="green",type='o',lwd=2,pch=18)
  lines(sub.data$sample,sub.data$OW.MAO.RMSE.SPCE3,col="magenta",type='o',lwd=2,pch=22)
}


plot(sub.data$sample,sub.data$OW.RMSE.RACE3,type='o',
     col="red",ylim = range(sub.data$OW.RMSE.RACE3,sub.data$IPW.RMSE.RACE3,sub.data$COX.Q.RMSE.RACE3,
                            sub.data$COX.MSM.RMSE.RACE3,
                            sub.data$IPW.AIPW.RMSE.RACE3,
                            sub.data$OW.AIPW.RMSE.RACE3),
     lwd=2,pch=13,xlab="Sample size",ylab="RMSE",main="RACE")
lines(sub.data$sample,sub.data$IPW.RMSE.RACE3,col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sub.data$COX.Q.RMSE.RACE3,col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sub.data$COX.MSM.RMSE.RACE3,col="orange",type='o',lwd=2,pch=16)
if(aipw.draw){
  lines(sub.data$sample,sub.data$OW.AIPW.RMSE.RACE3,col="green",type='o',lwd=2,pch=17)
  lines(sub.data$sample,sub.data$IPW.AIPW.RMSE.RACE3,col="magenta",type='o',lwd=2,pch=19)
}
if(pseudo.g.draw){
  lines(sub.data$sample,sub.data$PSEUDO.G.RMSE.RACE3,col="green",type='o',lwd=2,pch=20)
  if(good_overlap==1){
    lines(sub.data$sample,sub.data$PSEUDO.UNADJ.RMSE.RACE3,col="magenta",type='o',lwd=2,pch=21)
  }
}
if(mao.draw){
  lines(sub.data$sample,sub.data$IPW.MAO.RMSE.RACE3,col="green",type='o',lwd=2,pch=18)
  lines(sub.data$sample,sub.data$OW.MAO.RMSE.RACE3,col="magenta",type='o',lwd=2,pch=22)
}



plot(sub.data$sample,sub.data$OW.RMSE.ASCE3,type='o',
     col="red",ylim = range(sub.data$OW.RMSE.ASCE3,sub.data$IPW.RMSE.ASCE3,sub.data$COX.Q.RMSE.ASCE3,
                            sub.data$COX.MSM.RMSE.ASCE3,
                            sub.data$IPW.AIPW.RMSE.ASCE3,
                            sub.data$OW.AIPW.RMSE.ASCE3),
     lwd=2,pch=13,xlab="Sample size",ylab="RMSE",main="ASCE")
lines(sub.data$sample,sub.data$IPW.RMSE.ASCE3,col="blue",type='o',lwd=2,pch=14)
lines(sub.data$sample,sub.data$COX.Q.RMSE.ASCE3,col="black",type='o',lwd=2,pch=15)
lines(sub.data$sample,sub.data$COX.MSM.RMSE.ASCE3,col="orange",type='o',lwd=2,pch=15)
if(aipw.draw){
  lines(sub.data$sample,sub.data$OW.AIPW.RMSE.ASCE3,col="green",type='o',lwd=2,pch=17)
  lines(sub.data$sample,sub.data$IPW.AIPW.RMSE.ASCE3,col="magenta",type='o',lwd=2,pch=19)
}
if(pseudo.g.draw){
  lines(sub.data$sample,sub.data$PSEUDO.G.RMSE.ASCE3,col="green",type='o',lwd=2,pch=20)
  if(good_overlap==1){
    lines(sub.data$sample,sub.data$PSEUDO.UNADJ.RMSE.ASCE3,col="magenta",type='o',lwd=2,pch=21)
    
  }
}
if(mao.draw){
  lines(sub.data$sample,sub.data$IPW.MAO.RMSE.ASCE3,col="green",type='o',lwd=2,pch=18)
  lines(sub.data$sample,sub.data$OW.MAO.RMSE.ASCE3,col="magenta",type='o',lwd=2,pch=22)
}

if(!aipw.draw){
  plot(sub.data$sample,sub.data$OW.COVER.SPCE3,type='o',
       col="red",ylim = range(sub.data$OW.COVER.SPCE3,sub.data$IPW.COVER.SPCE3,sub.data$COX.Q.COVER.SPCE3,
                              sub.data$COX.MSM.COVER.SPCE3,1),
       lwd=2,pch=13,xlab="Sample size",ylab="COVER",main="SPCE")
  lines(sub.data$sample,sub.data$IPW.COVER.SPCE3,col="blue",type='o',lwd=2,pch=14)
  lines(sub.data$sample,sub.data$COX.Q.COVER.SPCE3,col="black",type='o',lwd=2,pch=15)
  lines(sub.data$sample,sub.data$COX.MSM.COVER.SPCE3,col="orange",type='o',lwd=2,pch=16)
  if(pseudo.g.draw){
    lines(sub.data$sample,sub.data$PSEUDO.G.COVER.SPCE3,col="green",type='o',lwd=2,pch=20)
    if(good_overlap==1){
      lines(sub.data$sample,sub.data$PSEUDO.UNADJ.COVER.SPCE3,col="magenta",type='o',lwd=2,pch=21)
    }
  }
  if(mao.draw){
    lines(sub.data$sample,sub.data$IPW.MAO.COVER.SPCE3,col="green",type='o',lwd=2,pch=18)
    lines(sub.data$sample,sub.data$OW.MAO.COVER.SPCE3,col="magenta",type='o',lwd=2,pch=22)
  }
  abline(h=0.95,lty=2)
  
  
  plot(sub.data$sample,sub.data$OW.COVER.RACE3,type='o',
       col="red",ylim = range(sub.data$OW.COVER.RACE3,
                              sub.data$IPW.COVER.RACE3,sub.data$COX.Q.COVER.RACE3,
                              sub.data$COX.MSM.COVER.RACE3,1),
       lwd=2,pch=13,xlab="Sample size",ylab="COVER",main="RACE")
  lines(sub.data$sample,sub.data$IPW.COVER.RACE3,col="blue",type='o',lwd=2,pch=14)
  lines(sub.data$sample,sub.data$COX.Q.COVER.RACE3,col="black",type='o',lwd=2,pch=15)
  lines(sub.data$sample,sub.data$COX.MSM.COVER.RACE3,col="orange",type='o',lwd=2,pch=16)
  if(pseudo.g.draw){
    lines(sub.data$sample,sub.data$PSEUDO.G.COVER.RACE3,col="green",type='o',lwd=2,pch=20)
    lines(sub.data$sample,sub.data$PSEUDO.UNADJ.COVER.RACE3,col="magenta",type='o',lwd=2,pch=21)
  }
  if(mao.draw){
    lines(sub.data$sample,sub.data$IPW.MAO.COVER.RACE3,col="green",type='o',lwd=2,pch=18)
    lines(sub.data$sample,sub.data$OW.MAO.COVER.RACE3,col="magenta",type='o',lwd=2,pch=22)
  }
  abline(h=0.95,lty=2)
  
  if(!mao.draw){
    ylim.range = range(sub.data$OW.COVER.ASCE3,sub.data$IPW.COVER.ASCE3,sub.data$COX.Q.COVER.ASCE3,
                       sub.data$COX.MSM.COVER.ASCE3,1)
  }else{
    ylim.range = range(sub.data$OW.COVER.ASCE3,sub.data$IPW.COVER.ASCE3,sub.data$COX.Q.COVER.ASCE3,
                       sub.data$COX.MSM.COVER.ASCE3,sub.data$IPW.MAO.COVER.ASCE3,sub.data$OW.MAO.COVER.ASCE3)
  }
  
  plot(sub.data$sample,sub.data$OW.COVER.ASCE3,type='o',
       col="red",ylim = ylim.range,
       lwd=2,pch=13,xlab="Sample size",ylab="COVER",main="ASCE")
  lines(sub.data$sample,sub.data$IPW.COVER.ASCE3,col="blue",type='o',lwd=2,pch=14)
  lines(sub.data$sample,sub.data$COX.Q.COVER.ASCE3,col="black",type='o',lwd=2,pch=15)
  lines(sub.data$sample,sub.data$COX.MSM.COVER.ASCE3,col="orange",type='o',lwd=2,pch=16)
  if(pseudo.g.draw){
    lines(sub.data$sample,sub.data$PSEUDO.G.COVER.ASCE3,col="green",type='o',lwd=2,pch=20)
    if(good_overlap==1){
      lines(sub.data$sample,sub.data$PSEUDO.UNADJ.COVER.ASCE3,col="magenta",type='o',lwd=2,pch=21)
      
    }
  }
  if(mao.draw){
    lines(sub.data$sample,sub.data$IPW.MAO.COVER.ASCE3,col="green",type='o',lwd=2,pch=18)
    lines(sub.data$sample,sub.data$OW.MAO.COVER.ASCE3,col="magenta",type='o',lwd=2,pch=22)
  }
  abline(h=0.95,lty=2)
}

par(mar = c(0.1,0.1,1,0.1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
if(aipw.draw){
  legend("top",legend=c("OW","T-IPW","Cox","IPW-Cox","AOW","AIPW"),
         lty=1,col=c("red","blue","black","orange","green","magenta"),lwd=2, horiz = T,
         pch=c(13,14,15,16,17,19))
}else if(pseudo.g.draw){
  if(good_overlap==1){
    legend("top",legend=c("OW","T-IPW","Cox","IPW-Cox","PO-G","PO-UNADJ"),
           lty=1,col=c("red","blue","black","orange","green","magenta"),lwd=2, horiz = T,
           pch=c(13,14,15,16,20,21))
  }else{
    legend("top",legend=c("OW","T-IPW","Cox","IPW-Cox","PO-G"),
           lty=1,col=c("red","blue","black","orange","green"),lwd=2, horiz = T,
           pch=c(13,14,15,16,20))
  }
  
}else if(mao.draw){
  legend("top",legend=c("OW","T-IPW","Cox","IPW-Cox","IPW-MAO","OW-MAO"),
         lty=1,col=c("red","blue","black","orange","green","magenta"),lwd=2, horiz = T,
         pch=c(13,14,15,16,18,22))
}else{
  legend("top",legend=c("OW","T-IPW","Cox","IPW-Cox"),
         lty=1,col=c("red","blue","black","orange"),lwd=2, ncol=4,
         pch=c(13,14,15,16))
}

