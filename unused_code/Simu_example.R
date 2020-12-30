sample_size=200
# Z=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.5,0.5))}))
# M1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.55-Z[x]*0.1,0.45+Z[x]*0.1))}))
# M2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.75-Z[x]*0.5,0.25+Z[x]*0.5))}))
# P1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.4+Z[x]*0.1+M1[x]*0.1,0.6-Z[x]*0.1-M1[x]*0.1))}))
# P2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.4+Z[x]*0.1+(M1[x]+M2[x])*0.1,0.6-Z[x]*0.1-(M1[x]+M2[x])*0.1))}))
# Survival=P1+P1*P2
#pm1=0.4;pm2=0.1

effect_on_M_results=NULL
M_on_Y_results=NULL
Indirect_effect_results=NULL

pm0=0.2
# pm1=0.1;pm2=0.1
pm1_grid=seq(0,0.5,length=11)
#pm2_grid=seq(0.25,0.35,length=3)


set.seed(100)
# effect on M
for (early_effect in c(T,F)){
for (pm1 in pm1_grid){
pm2=0.5-pm1
# for (pm2 in pm2_grid){
n_simu=1000
effect_on_M=numeric(n_simu)
M_on_Y=numeric(n_simu)
Indirect_effect=numeric(n_simu)
for (n in 1:n_simu){
Z=rep(0,sample_size)
if(early_effect){
  M1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.55-Z[x]*0.3,0.45+Z[x]*0.3))}))
  M2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.75-Z[x]*0.1,0.25+Z[x]*0.1))}))
}else{
  M1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.55-Z[x]*0.1,0.45+Z[x]*0.1))}))
  M2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.75-Z[x]*0.5,0.25+Z[x]*0.5))}))
}
P1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.4+Z[x]*0.1+M1[x]*pm0,0.6-Z[x]*0.1-M1[x]*pm0))}))
P2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.4+Z[x]*0.1+M1[x]*pm1+M2[x]*pm2,0.6-Z[x]*0.1-(M1[x]*pm1+M2[x]*pm2)))}))
M_0=(M1+M2*P1)/(1+P1)


Z=rep(1,sample_size)
if(early_effect){
  M1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.55-Z[x]*0.3,0.45+Z[x]*0.3))}))
  M2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.75-Z[x]*0.1,0.25+Z[x]*0.1))}))
}else{
  M1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.55-Z[x]*0.1,0.45+Z[x]*0.1))}))
  M2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.75-Z[x]*0.5,0.25+Z[x]*0.5))}))
}
M_1=(M1+M2*P1)/(1+P1)
effect_on_M[n]=mean(M_1)-mean(M_0)

# effect of M on Survival
# effect on M
Z=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.5,0.5))}))
# M1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.55-Z[x]*0.1,0.45+Z[x]*0.1))}))
# M2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.75-Z[x]*0.5,0.25+Z[x]*0.5))}))
M1=rep(0,sample_size)
M2=rep(0,sample_size)
P1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.4+Z[x]*0.1+M1[x]*pm0,0.6-Z[x]*0.1-M1[x]*pm0))}))
P2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.4+Z[x]*0.1+M1[x]*pm1+M2[x]*pm2,0.6-Z[x]*0.1-(M1[x]*pm1+M2[x]*pm2)))}))
Survival_0=P1+P1*P2


# M1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.55-Z[x]*0.1,0.45+Z[x]*0.1))}))
# M2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.75-Z[x]*0.5,0.25+Z[x]*0.5))}))
M1=M1+1;M2=M2+1
P1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.4+Z[x]*0.1+M1[x]*pm0,0.6-Z[x]*0.1-M1[x]*pm0))}))
P2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.4+Z[x]*0.1+M1[x]*pm1+M2[x]*pm2,0.6-Z[x]*0.1-(M1[x]*pm1+M2[x]*pm2)))}))
Survival_1=P1+P1*P2
M_on_Y[n]=mean(Survival_1)-mean(Survival_0)


# Indirect effect
# effect on M
Z=rep(0,sample_size)
if(early_effect){
  M1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.55-Z[x]*0.3,0.45+Z[x]*0.3))}))
  M2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.75-Z[x]*0.1,0.25+Z[x]*0.1))}))
}else{
  M1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.55-Z[x]*0.1,0.45+Z[x]*0.1))}))
  M2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.75-Z[x]*0.5,0.25+Z[x]*0.5))}))
}

P1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.4+M1[x]*pm0,0.6-M1[x]*pm0))}))
P2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.4+M1[x]*pm1+M2[x]*pm2,0.6-M1[x]*pm1-M2[x]*pm2))}))
Survival_0=P1+P1*P2

Z=rep(1,sample_size)
if(early_effect){
  M1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.55-Z[x]*0.3,0.45+Z[x]*0.3))}))
  M2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.75-Z[x]*0.1,0.25+Z[x]*0.1))}))
}else{
  M1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.55-Z[x]*0.1,0.45+Z[x]*0.1))}))
  M2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.75-Z[x]*0.5,0.25+Z[x]*0.5))}))
}

P1=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.4+M1[x]*pm0,0.6-M1[x]*pm0))}))
P2=unlist(lapply(1:sample_size,FUN=function(x){sample(c(0,1),1,prob=c(0.4+M1[x]*pm1+M2[x]*pm2,0.6-M1[x]*pm1-M2[x]*pm2))}))
Survival_1=P1+P1*P2
Indirect_effect[n]=mean(Survival_1)-mean(Survival_0)
if (n%%500==0)
{
  print(paste(pm1,"==",pm2,"== MC",n,"=="))
}
}

effect_on_M_results=rbind(effect_on_M_results,c(pm1,pm2,quantile(effect_on_M,c(0.025,0.5,0.975))))
M_on_Y_results=rbind(M_on_Y_results,c(pm1,pm2,quantile(M_on_Y,c(0.025,0.5,0.975))))
Indirect_effect_results=rbind(Indirect_effect_results,c(pm1,pm2,quantile(Indirect_effect,c(0.025,0.5,0.975))))

#}
}
}
# effect_on_M_results_LATER=effect_on_M_results
# M_on_Y_results_LATER=effect_on_M_results
# Indirect_effect_results_LATER=Indirect_effect_results
# save(effect_on_M_results_LATER,M_on_Y_results_LATER,Indirect_effect_results_LATER,file="later_simu.RDATA")
# pdf("summary.pdf",width=10,height=3)
par(mfrow=c(1,3))
plot(effect_on_M_results[1:11,2],effect_on_M_results[1:11,4],ylim=c(0.05,0.35),type='o',main="Effect on Mediator",
     xlab="Stength of effect on later stage from M",ylab="Effect on M")
lines(effect_on_M_results[1:11,2],effect_on_M_results[1:11,3],lty=2)
lines(effect_on_M_results[1:11,2],effect_on_M_results[1:11,5],lty=2)
points(effect_on_M_results[12:22,2],effect_on_M_results[12:22,4],col="blue",lwd=2)
lines(effect_on_M_results[12:22,2],effect_on_M_results[12:22,4],lty=1,col="blue")
lines(effect_on_M_results[12:22,2],effect_on_M_results[12:22,3],lty=2,col="blue")
lines(effect_on_M_results[12:22,2],effect_on_M_results[12:22,5],lty=2,col="blue")

plot(effect_on_M_results[1:11,2],M_on_Y_results[1:11,4],ylim=c(-0.7,-0.2),type='o',main="Effect of Mediator on Survival",
     xlab="Strength of effect on later stage from M",ylab="Effect on survival")
lines(effect_on_M_results[1:11,2],M_on_Y_results[1:11,3],lty=2)
lines(effect_on_M_results[1:11,2],M_on_Y_results[1:11,5],lty=2)
points(effect_on_M_results[12:22,2],M_on_Y_results[12:22,4],col="blue",lwd=2)
lines(effect_on_M_results[12:22,2],M_on_Y_results[12:22,4],lty=1,col="blue")
lines(effect_on_M_results[12:22,2],M_on_Y_results[12:22,3],lty=2,col="blue")
lines(effect_on_M_results[12:22,2],M_on_Y_results[12:22,5],lty=2,col="blue")
legend("top",legend=c("M change early","M change later"),col=c("black","blue"),lty=1,lwd=2)

plot(effect_on_M_results[1:11,2],Indirect_effect_results[1:11,4],ylim=c(-0.4,0.2),type='o',main="Indirect Effect",
     xlab="Strength of effect on later stage from M",ylab="Indirect effect on survival",lwd=2)
lines(effect_on_M_results[1:11,2],Indirect_effect_results[1:11,3],lty=2)
lines(effect_on_M_results[1:11,2],Indirect_effect_results[1:11,5],lty=2)
points(effect_on_M_results[12:22,2],Indirect_effect_results[12:22,4],col="blue")
lines(effect_on_M_results[12:22,2],Indirect_effect_results[12:22,4],lty=1,col="blue",lwd=2)
lines(effect_on_M_results[12:22,2],Indirect_effect_results[12:22,3],lty=2,col="blue")
lines(effect_on_M_results[12:22,2],Indirect_effect_results[12:22,5],lty=2,col="blue")
abline(h=0,lty=2)
# dev.off()