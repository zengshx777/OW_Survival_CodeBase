# set.seed(201)
# pdf("balanced_condition.pdf",width = 10,height = 4)
# m <- matrix(c(1,2,3,4,4,4),nrow = 2,ncol = 3,byrow = TRUE)
# layout(mat = m,heights = c(0.5,0.06))
# pdf("GPS_plot_simu.pdf",width=12,height=12)
# m <- matrix(c(1,2,3,4,5,6,7,8,9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)
# layout(mat = m,heights = c(0.5,0.5,0.5,0.06))
# for (good_overlap in c(1,2,3)){
# Lower dimension of covariates, same from Lunceford setting
X4 = sample(c(0, 1),
            sample_size,
            replace = T,
            prob = c(0.5, 0.5))
X3 = unlist(lapply(
  1:sample_size,
  FUN = function(x) {
    sample(c(0, 1), 1, prob = c(0.6 - 0.2 * X4[x], 0.2 * X4[x] + 0.4))
  }
))
X = cbind(1, mvrnorm(
  n = sample_size,
  mu = c(0, 0),
  Sigma = matrix(c(2, 0.5, 0.5, 2), 2, 2)
), X3, X4)
colnames(X) = c("X0", "X1", "X2", "X3", "X4")
# Propensity Score Model
beta.0 = c(0, 0, 0, 0, 0)
if (good_overlap == 1) {
  beta.1 = c(0, 0, 0, 0, 0)
  beta.2 = c(0, 0, 0, 0, 0)
} else if (good_overlap == 2) {
  beta.1 = c(-0.4, 0.85, 0.9, 0.45,-0.25)
  beta.2 = -beta.1 * 0.2
} else{
  beta.1 = c(1.2, 1.5, 1,-1.5,-1)
  beta.2 = - beta.1 * 0.2
}
if (multi.arm) {
  beta = cbind(beta.0, beta.1, beta.2)
} else{
  beta = cbind(beta.0, beta.1)
}
exp.part = exp(X %*% beta)
soft.max = exp.part / rowSums(exp.part)


# true_e=1/(1+exp(-X%*%beta))
Z = unlist(lapply(
  1:sample_size,
  FUN = function(x) {
    sample(0:(ncol(soft.max) - 1), 1, prob = soft.max[x,])
  }
))
# if(good_overlap==1){
#   main.text = "\n good overlap"
#   xlab.text = ""
# }else if(good_overlap == 2){
#   main.text = "\n moderate overlap"
#   xlab.text = ""
# }else{
#   main.text = "\n poor overlap"
#   xlab.text = "True GPS"
# }
# ## Snippet for draw PC plot
# # pc.object = prcomp(X[,-1])
# # score = pc.object$x
# # plot(score[Z==0,c(1,2)],col = rgb(red = 0, green = 0, blue = 1, alpha = 0.7),
# #      main = main.text, pch = 17 ,xlim = c(-4,4), ylim = c(-3,3))
# # points(score[Z==1,c(1,2)],col = rgb(red = 1, green = 0, blue = 0, alpha = 0.7), pch = 17)
# # points(score[Z==2,c(1,2)],col = rgb(red = 0, green = 0, blue = 0, alpha = 0.5), pch = 17)
# if(good_overlap==1){ymax=100}else{
#   ymax=max(density(soft.max[Z==0,1])$y,density(soft.max[Z==1,1])$y,density(soft.max[Z==2,1])$y)
# }
# plot(density(soft.max[Z==0,1]+rnorm(sum(Z==0),0,0.005)),col="red",xlab = xlab.text,bty="n",ylab="Density",
#      main=paste("Propensity for being assigned to ARM 1",main.text), lwd=2,yaxt='n',
#      ylim=c(0,ymax),xlim=c(0,1))
# lines(density(soft.max[Z==1,1]+rnorm(sum(Z==1),0,0.005)),col="blue", lwd=2)
# lines(density(soft.max[Z==2,1]+rnorm(sum(Z==2),0,0.005)),col="green",lwd=2)
# if(good_overlap==1){ymax=100}else{
#   ymax=max(density(soft.max[Z==0,2])$y,density(soft.max[Z==1,2])$y,density(soft.max[Z==2,2])$y)
# }
# plot(density(soft.max[Z==0,2]+rnorm(sum(Z==0),0,0.005)),col="red",xlab = xlab.text,bty="n",ylab="",
#      main=paste("Propensity for being assigned to ARM 2",main.text), lwd=2,yaxt='n',
#      ylim=c(0,ymax),xlim=c(0,1))
# lines(density(soft.max[Z==1,2]+rnorm(sum(Z==1),0,0.005)),col="blue", lwd=2)
# lines(density(soft.max[Z==2,2]+rnorm(sum(Z==2),0,0.005)),col="green",lwd=2)
# if(good_overlap==1){ymax=100}else{
#   ymax=max(density(soft.max[Z==0,3])$y,density(soft.max[Z==1,3])$y,density(soft.max[Z==2,3])$y)
# }
# plot(density(soft.max[Z==0,3]+rnorm(sum(Z==0),0,0.005)),col="red",xlab = xlab.text,bty="n",ylab="",
#      main=paste("Propensity for being assigned to ARM 3",main.text), lwd=2,yaxt='n',
#      ylim=c(0,ymax),xlim=c(0,1))
# lines(density(soft.max[Z==1,3]+rnorm(sum(Z==1),0,0.005)),col="blue", lwd=2)
# lines(density(soft.max[Z==2,3]+rnorm(sum(Z==2),0,0.005)),col="green",lwd=2)
# }
# par(mar = c(0.1,0.1,1,0.1))
# plot(1, type = "n", axes=FALSE, xlab="", ylab="")
# # legend("top",legend=c("Z=1","Z=2","Z=3"),
# #        col=c("red","blue","black"), ncol=3,
# #        pch=17)
# legend("top",legend=c("Z=1","Z=2","Z=3"),
#        col=c("red","blue","green"), ncol=3,
#        pch=17)
# dev.off()
# Survival time
if (prop.hazard) {
  gamma.1 = 0.5
  gamma.2 = 0.5
  alpha = c(2, 1.5,-1, 1)
  lambda = 0.0001
  v = 3
  truncate = truncate.ph
  # The covariate and the treatment is defined on proportional hazard
  hazard = as.vector(gamma.1 * as.numeric(Z == 1) + gamma.2 * as.numeric(Z ==
                                                                           2) + X[,-1] %*% alpha)
  random_simu = runif(sample_size)
  Survival_time = (-log(random_simu) / (lambda * exp(hazard))) ^ (1 / v)
  
  # get true est
  hazard_0 = as.vector(X[,-1] %*% alpha)
  hazard_1 = as.vector(gamma.1 + X[,-1] %*% alpha)
  hazard_2 = as.vector(gamma.2 + X[,-1] %*% alpha)
  Survival_time_0 = (-log(random_simu) / (lambda * exp(hazard_0))) ^ (1 /
                                                                        v)
  Survival_time_1 = (-log(random_simu) / (lambda * exp(hazard_1))) ^ (1 /
                                                                        v)
  Survival_time_2 = (-log(random_simu) / (lambda * exp(hazard_2))) ^ (1 /
                                                                        v)
} else{
  # AFT but not PH; log-normal
  gamma.1 = -1
  gamma.2 = -1
  alpha = -c(2, 1.5,-1, 1)
  lambda = 0.0001
  v = 3
  truncate = truncate.aft
  random_simu = rnorm(sample_size, sd = 0.81)
  Survival_time = exp(3.5 + 0.2 * (
    as.vector(
      gamma.1 * as.numeric(Z == 1) + gamma.2 * as.numeric(Z == 2) + X[,-1] %*%
        alpha
    )
  ) + random_simu)
  Survival_time_0 = exp(3.5 + 0.2 * (X[,-1] %*% alpha) + random_simu)
  Survival_time_1 = exp(3.5 + 0.2 * (gamma.1 + X[,-1] %*% alpha) + random_simu)
  Survival_time_2 = exp(3.5 + 0.2 * (gamma.2 + X[,-1] %*% alpha) + random_simu)
}
# Censor time
if (dependent.censoring) {
  lambda.c = 0.0001
  v.c = 2.7
  alpha.c = c(1, 0.5,-0.5, 0.5)  # The covariate and the treatment is defined on proportional hazard
  c.hazard = as.vector(X[,-1] %*% alpha)
  random_simu = runif(sample_size)
  Censor_time = (-log(random_simu) / (lambda.c * exp(c.hazard))) ^ (1 / v.c)
  # mean(Survival_time > Censor_time)
} else{
  Censor_time = runif(sample_size, min = 0, max = 115)
  # mean(Survival_time > Censor_time)
}

# Censor_time=runif(sample_size,min=2000,max=2115)

# Observed survival
Censored = Survival_time > Censor_time
Y = pmin(Survival_time, Censor_time)
DELTA = 1 - Censored

# Store true value on finite sample
# truncate=min(quantile(Survival_time_1[Z==1],0.75),quantile(Survival_time_0[Z==0],0.75))
tilt.h = 1 / rowSums(1 / soft.max)

tilt.h = 1 / rowSums(1 / soft.max)
true_est_finite[i, 1] = mean(Survival_time_1) - mean(Survival_time_0)
true_est_finite[i, 2] = mean(Survival_time_2) - mean(Survival_time_0)
true_est_finite[i, 3] = mean(Survival_time_2) - mean(Survival_time_1)

true_est_finite[i, 4] = mean(pmin(Survival_time_1, truncate)) - mean(pmin(Survival_time_0, truncate))
true_est_finite[i, 5] = mean(pmin(Survival_time_2, truncate)) - mean(pmin(Survival_time_0, truncate))
true_est_finite[i, 6] = mean(pmin(Survival_time_2, truncate)) - mean(pmin(Survival_time_1, truncate))

true_est_finite[i, 7] =  mean(Survival_time_1 > truncate) - mean(Survival_time_0 >
                                                                   truncate)
true_est_finite[i, 8] =  mean(Survival_time_2 > truncate) - mean(Survival_time_0 >
                                                                   truncate)
true_est_finite[i, 9] =  mean(Survival_time_2 > truncate) - mean(Survival_time_1 >
                                                                   truncate)

true_est_ow_finite[i, 1] = sum((Survival_time_1 - Survival_time_0) * tilt.h) /
  sum(tilt.h)
true_est_ow_finite[i, 2] = sum((Survival_time_2 - Survival_time_0) * tilt.h) /
  sum(tilt.h)
true_est_ow_finite[i, 3] = sum((Survival_time_2 - Survival_time_1) * tilt.h) /
  sum(tilt.h)

true_est_ow_finite[i, 4] = sum((
  pmin(Survival_time_1, truncate) - pmin(Survival_time_0, truncate)
) * tilt.h) /
  sum(tilt.h)
true_est_ow_finite[i, 5] = sum((
  pmin(Survival_time_2, truncate) - pmin(Survival_time_0, truncate)
) * tilt.h) /
  sum(tilt.h)
true_est_ow_finite[i, 6] = sum((
  pmin(Survival_time_2, truncate) - pmin(Survival_time_1, truncate)
) * tilt.h) /
  sum(tilt.h)

true_est_ow_finite[i, 7] = sum((
  as.numeric(Survival_time_1 > truncate) - as.numeric(Survival_time_0 > truncate)
) * tilt.h) / sum(tilt.h)
true_est_ow_finite[i, 8] = sum((
  as.numeric(Survival_time_2 > truncate) - as.numeric(Survival_time_0 > truncate)
) * tilt.h) / sum(tilt.h)
true_est_ow_finite[i, 9] = sum((
  as.numeric(Survival_time_2 > truncate) - as.numeric(Survival_time_1 > truncate)
) * tilt.h) / sum(tilt.h)
