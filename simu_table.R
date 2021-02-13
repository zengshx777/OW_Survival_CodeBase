library(xtable)
helper<-function(results,data, est.id=1)
{
  if(nrow(data)==0){stop()}
  if(est.id==1){
    results=rbind(results,cbind(
      data$OW.RMSE.ASCE3,data$IPW.RMSE.ASCE3,data$COX.Q.RMSE.ASCE3,
      data$COX.MSM.RMSE.ASCE3,
      #data$IPW.MAO.RMSE.ASCE,
      data$OW.BIAS.ASCE3,data$IPW.BIAS.ASCE3,data$COX.Q.BIAS.ASCE3,
      data$COX.MSM.BIAS.ASCE3,
      #data$IPW.MAO.BIAS.ASCE,
      data$OW.COVER.ASCE3,data$IPW.COVER.ASCE3,data$COX.Q.COVER.ASCE3,
      data$COX.MSM.COVER.ASCE3
      #data$IPW.MAO.COVER.ASCE
  ))}else if(est.id == 2){
    results=rbind(results,cbind(
      data$OW.RMSE.RACE3,data$IPW.RMSE.RACE3,data$COX.Q.RMSE.RACE3,
      data$COX.MSM.RMSE.RACE3,
      #data$IPW.MAO.RMSE.RACE,
      data$OW.BIAS.RACE3,data$IPW.BIAS.RACE3,data$COX.Q.BIAS.RACE3,
      data$COX.MSM.BIAS.RACE3,
      #data$IPW.MAO.BIAS.RACE,
      data$OW.COVER.RACE3,data$IPW.COVER.RACE3,data$COX.Q.COVER.RACE3,
      data$COX.MSM.COVER.RACE3
      #data$IPW.MAO.COVER.RACE
      ))
  }else{
    results=rbind(results,cbind(
      data$OW.RMSE.SPCE3,data$IPW.RMSE.SPCE3,data$COX.Q.RMSE.SPCE3,
      data$COX.MSM.RMSE.SPCE3,
      #data$IPW.MAO.RMSE.SPCE,
      data$OW.BIAS.SPCE3,data$IPW.BIAS.SPCE3,data$COX.Q.BIAS.SPCE3,
      data$COX.MSM.BIAS.SPCE3,
      #data$IPW.MAO.BIAS.SPCE,
      data$OW.COVER.SPCE3,data$IPW.COVER.SPCE3,data$COX.Q.COVER.SPCE3,
      data$COX.MSM.COVER.SPCE3
      #data$IPW.MAO.COVER.SPCE
    ))
  }
  return(results)
}

produce_latex<-function(results){
  latex_line=NULL;text=""
  for (k in 1:nrow(results)){
    if(k%%3==1){
      text="Good"
      skip.text="\\smallskip"
      #skip.text=""
    }
    if(k%%3==2){
      text="Mod"
      skip.text="\\smallskip"
      next
    }
    if(k%%3==0){
      text="Poor"
      skip.text=""
    }
    est.text=""
    if(k%%9==1){
      est.text="$ \\tau_{j,j'}^{\\ASCE,h} $"
    }
    if(k%%9==4){
      est.text="$ \\tau_{j,j'}^{\\RACE,h} $"
    }
    if(k%%9==7){
      est.text="$ \\tau_{j,j'}^{\\SPCE,h} $"
    }

    
    latex_line = c(latex_line,paste(est.text,"&",text,
                                    paste(paste("&",sprintf("%.3f", round(results[k,],3))),collapse =""),"\\\\",skip.text,"\n"
    ))
  }
  
  cat(latex_line)
}

results = NULL
results=helper(results,subset(collect.data,sample==300&arm==1&prop==T&dependent==0),est.id=1)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==T&dependent==0),est.id=2)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==T&dependent==0),est.id=3)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==F&dependent==0),est.id=1)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==F&dependent==0),est.id=2)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==F&dependent==0),est.id=3)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==T&dependent==1),est.id=1)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==T&dependent==1),est.id=2)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==T&dependent==1),est.id=3)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==F&dependent==1),est.id=1)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==F&dependent==1),est.id=2)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==F&dependent==1),est.id=3)
# rownames(results)=rep(c("None","Mod","Poor"),12)
# print(xtable(results,digits=c(0,rep(3,12)),include.rownames=FALSE))


produce_latex(results[1:9,])
produce_latex(results[10:18,])
produce_latex(results[19:27,])
produce_latex(results[28:36,])
