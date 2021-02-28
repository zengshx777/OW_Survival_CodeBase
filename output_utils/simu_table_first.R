library(xtable)
helper<-function(results,data, est.id=1)
{
  if(nrow(data)==0){stop()}
  if(est.id==1){
    results=rbind(results,cbind(
      data$OW.BIAS.ASCE1,data$IPW.BIAS.ASCE1,data$COX.Q.BIAS.ASCE1,
      data$COX.MSM.BIAS.ASCE1,
      #data$IPW.MAO.BIAS.ASCE,
      data$OW.RMSE.ASCE1,data$IPW.RMSE.ASCE1,data$COX.Q.RMSE.ASCE1,
      data$COX.MSM.RMSE.ASCE1,
      #data$IPW.MAO.RMSE.ASCE,
      data$OW.COVER.ASCE1,data$IPW.COVER.ASCE1,data$COX.Q.COVER.ASCE1,
      data$COX.MSM.COVER.ASCE1
      #data$IPW.MAO.COVER.ASCE
    ))}else if(est.id == 2){
      results=rbind(results,cbind(
        data$OW.BIAS.RACE1,data$IPW.BIAS.RACE1,data$COX.Q.BIAS.RACE1,
        data$COX.MSM.BIAS.RACE1,
        #data$IPW.MAO.BIAS.RACE,
        data$OW.RMSE.RACE1,data$IPW.RMSE.RACE1,data$COX.Q.RMSE.RACE1,
        data$COX.MSM.RMSE.RACE1,
        #data$IPW.MAO.RMSE.RACE,
        data$OW.COVER.RACE1,data$IPW.COVER.RACE1,data$COX.Q.COVER.RACE1,
        data$COX.MSM.COVER.RACE1
        #data$IPW.MAO.COVER.RACE
      ))
    }else{
      results=rbind(results,cbind(
        data$OW.BIAS.SPCE1,data$IPW.BIAS.SPCE1,data$COX.Q.BIAS.SPCE1,
        data$COX.MSM.BIAS.SPCE1,
        #data$IPW.MAO.BIAS.SPCE,
        data$OW.RMSE.SPCE1,data$IPW.RMSE.SPCE1,data$COX.Q.RMSE.SPCE1,
        data$COX.MSM.RMSE.SPCE1,
        #data$IPW.MAO.RMSE.SPCE,
        data$OW.COVER.SPCE1,data$IPW.COVER.SPCE1,data$COX.Q.COVER.SPCE1,
        data$COX.MSM.COVER.SPCE1
        #data$IPW.MAO.COVER.SPCE
      ))
    }
  return(results)
}

produce_latex<-function(results){
  latex_line=NULL;text=""
  for (k in 1:nrow(results)){
    if(k%%3==1){
      text="RCT"
      skip.text="\\smallskip"
      next
      #skip.text=""
    }
    if(k%%3==2){
      text="Good"
      skip.text="\\smallskip"
    }
    if(k%%3==0){
      text="Poor"
      skip.text=""
    }
    est.text=""
    if(k%%9==2){
      # est.text="$ \\tau_{j,j'}^{\\ASCE,h} $"
      est.text="SPCE"
    }
    if(k%%9==5){
      #est.text="$ \\tau_{j,j'}^{\\RACE,h} $"
      est.text = "RACE"
    }
    if(k%%9==8){
      # est.text="$ \\tau_{j,j'}^{\\SPCE,h} $"
      est.text = "ASCE"
    }
    
    
    latex_line = c(latex_line,paste(est.text,"&",text,
                                    paste(paste("&",sprintf("%.3f", round(results[k,],3))),collapse =""),"\\\\",skip.text,"\n"
    ))
  }
  
  cat(latex_line)
}
results = NULL
results=helper(results,subset(collect.data,sample==300&arm==1&prop==T&dependent==0),est.id=3)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==T&dependent==0),est.id=2)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==T&dependent==0),est.id=1)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==F&dependent==0),est.id=3)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==F&dependent==0),est.id=2)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==F&dependent==0),est.id=1)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==T&dependent==1),est.id=3)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==T&dependent==1),est.id=2)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==T&dependent==1),est.id=1)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==F&dependent==1),est.id=3)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==F&dependent==1),est.id=2)
results=helper(results,subset(collect.data,sample==300&arm==1&prop==F&dependent==1),est.id=1)
# rownames(results)=rep(c("None","Mod","Poor"),12)
# print(xtable(results,digits=c(0,rep(3,12)),include.rownames=FALSE))


produce_latex(results[1:9,])
produce_latex(results[10:18,])
produce_latex(results[19:27,])
produce_latex(results[28:36,])
