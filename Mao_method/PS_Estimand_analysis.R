###########################################################################################################
##
## Example code for propensity score weighting adjusted analysis for estimation of survival estimands 
## 05/25/2018 @ Huzhang Mao
## 
##  WIN 7, R version 3.3.3 (2017-03-06)
###########################################################################################################

rm( list=ls( all = T ) );
date.tag <- "2018_05_25_";

set.seed(123);
setwd("C:/Users/Shuxi ZENG/Dropbox/Fourth Year/OW_Survival/codebase/Mao_method")
library( survival );
source( "functions//func_PS_Estimand_analysis.R" );
source( "functions//GIPW_function_omega.R" );


##################################################################################
## Data generation:
## Note
##  1. Treatment variable musgt be 0 or 1;
##  2. Categorical variable must be dummy coded to 0 or 1 for each level
##################################################################################
n.sample <- 500;
PID <- c( 1 : n.sample );
X0 <- rep( 1, n.sample );
X1 <- rnorm( n.sample, 0, 1 );
X2 <- rnorm( n.sample, 1, 2 );
X3 <- rbinom( n.sample, 1, 0.5 );
X4 <- rbinom( n.sample, 1, 0.3 );
X.sum <- X1 + X2 + X3 + X4;
Z <- rbinom( n.sample, 1, exp( X.sum )/( 1 + exp(X.sum) ) );
Y <- rexp( n.sample, 0.2 );
censor.time <- runif( n.sample, 0, 0.5*max(Y) );
delta <- as.numeric( Y <= censor.time ); # event indicator

X4PS <- c( "X1", "X2", "X3", "X4" );
data <- as.data.frame( cbind( PID, X0, X1, X2, X3, X4, Z, Y, delta ) );
write.csv( data, "Example_data.csv" );

#####################################################################
## Setup for PS_Estimand analysis
#####################################################################
setup <- list( nsample = 0,
               sp.degree = 1,
               nknots = 10,
               trt.type = ' ',
               weight.type = ' ',
               overlap = ' ',
               data.type = ' ',
               race.time = 4,
               spce.time = 4,
               sqe.quantile = 0.5 );

#############################################################################
## Note
## !!! The following variable names are fixed !!!
#############################################################################
IDname <- "PID";
prefix.knots <- "lambda";
prefix.fix <- "time";
knots.name <- paste(set.basis.name(nbasis=setup$nknots, prefix=prefix.knots), sep=',');
var.name <- c( "X0", t( sapply( c(1:setup$sp.degree), function(x) paste( prefix.fix, x, sep='') ) ), knots.name);
out.name <- "Yig";
x.name <- c( "X0", X4PS );
offset.name <- "Cig";
wt.name <- "wt";
censor.name <- "delta";
tname <- "Z";
time.name <- "Y";
alpha1.name <- sapply( 1 : (setup$nknots+setup$sp.degree+1), function(x) paste( "alpha1.", x, sep='') );
alpha0.name <- sapply( 1 : (setup$nknots+setup$sp.degree+1), function(x) paste( "alpha0.", x, sep='') );
beta.name <- sapply( c( 0 : (length(x.name)-1) ) , function(x) paste( "beta", x, sep='') );
ps.form <- as.formula( paste( "Z ~ ", paste( X4PS, collapse = " + " ), sep="" ) );
surv.form <- as.formula( 'Surv(Y, delta) ~ Z' );
D0 <- D1 <- cbind(matrix(0, ncol=(setup$sp.degree+1), nrow=setup$nknots), diag(setup$nknots));
par.start <- rep(0, length = (setup$sp.degree+setup$nknots+1) );


##############################################################################################
############# Estimand analysis ##################
##############################################################################################

write.poisson.data <- function( data, weight.type ) {
  # poissonization for to calculate pentalty parameter (sigma)
  #   data: data to be used
  #   weight.type: PS weighting method
  
  dir.create(file.path("C:/Users/Shuxi ZENG/Dropbox/Fourth Year/OW_Survival/codebase/Mao_method", "poisson_data"), showWarnings = FALSE)
  
  out.ps <- ps.model(dat=data, form=ps.form );
  data$ps <- out.ps$ps.hat;  # estimated ps
  
  if(weight.type!="UNWEIGHT"){
    omega <- sapply( data$ps, calc.omega, weight = weight.type );
  }else{
    omega = data$ps * data$Z + (1-data$ps) * (1-data$Z)
  }
  data$wt <- omega / ( data$ps * data$Z + (1-data$ps) * (1-data$Z) );
  
  # treated group
  Z <- 1;
  tmp <- cox.to.poisson( dat=data[ data$Z==Z, ], nknots=setup$nknots, dgr=setup$sp.degree, basis.name=var.name, out.name=out.name, 
                         offset.name=offset.name, tname=tname, censor.name=censor.name, wt.name=wt.name, ID.name=IDname, time.name=time.name );
  write.csv( tmp$dat.poisson, paste( "poisson_data//data_poisson_estimand",weight.type, "_", ifelse(Z==1, "treated", "control"), ".csv", sep="" ), row.names = F );
  
  # control group
  Z <- 0;
  tmp <- cox.to.poisson( dat=data[ data$Z==Z, ], nknots=setup$nknots, dgr=setup$sp.degree, basis.name=var.name, out.name=out.name, 
                         offset.name=offset.name, tname=tname, censor.name=censor.name, wt.name=wt.name, ID.name=IDname, time.name=time.name );
  write.csv( tmp$dat.poisson, paste( "poisson_data//data_poisson_estimand",weight.type, "_", ifelse(Z==1, "treated", "control"), ".csv", sep="" ), row.names = F );
}

# IPW
write.poisson.data( data = data, weight.type="IPW" );
# MW
write.poisson.data( data = data, weight.type="MW" );
# OW
write.poisson.data( data = data, weight.type="OVERLAP" );
# UNWEIGHT
write.poisson.data( data = data, weight.type="UNWEIGHT" );

estimand.analysis.fun <- function( data, data.setup ) {
  # function for PS_Estimand analysis
  # 
  # data: data to used
  # data.setup: analysis setup up
  
  out.ps <- ps.model(dat=data, form=ps.form );
  beta.hat <- coef( out.ps$fm );
  names( beta.hat ) <- beta.name;
  data$ps <- out.ps$ps.hat;  # estimated ps
  if(setup$weight.type!="UNWEIGHT"){
    omega <- sapply( data$ps, calc.omega, weight=setup$weight.type ) / 2;
  }else{
    omega = data$ps * data$Z + (1-data$ps) * (1-data$Z)
  }
  data$wt <- omega / ( data$ps * data$Z + (1-data$ps) * (1-data$Z) );

  ## Treated group
  Z <- 1;
  data.treat <- data[ data$Z == Z, ];
  res1 <- calc.par( dat=data.treat, par.start=par.start, sigma=sigma1, nknots=data.setup$nknots, Dmat=D1, censor.name=censor.name, dgr=data.setup$sp.degree,
                    weight=data.setup$weight.type, out.name=out.name, offset.name=offset.name, wt.name=wt.name, basis.name=var.name, tname=tname, 
                    ID.name=IDname, time.name=time.name );
  pois.treat <- res1$pois.dat;
  alpha1 <- res1$alpha;
  names( alpha1 ) <- alpha1.name;
  alpha1.var <- res1$alpha.var;
  knots1 <- res1$knots;
  ## Control group
  Z <- 0;
  data.control <- data[ data$Z == Z, ];
  res0 <- calc.par( dat=data.control, par.start=par.start, sigma=sigma0, nknots=data.setup$nknots, Dmat=D0, censor.name=censor.name, dgr=data.setup$sp.degree,
                    weight=data.setup$weight.type, out.name=out.name, offset.name=offset.name, wt.name=wt.name, basis.name=var.name, tname=tname, 
                    ID.name=IDname, time.name=time.name );
  pois.control <- res0$pois.dat;
  alpha0 <- res0$alpha;
  names( alpha0 ) <- alpha0.name;
  alpha0.var <- res0$alpha.var;
  knots0 <- res0$knots;
  
  ## combine data after poissonization
  data.poisson <- rbind( pois.treat, pois.control );
  colnames( data.poisson );
  ## Q matrix 
  Q <- matrix( 0, nrow=2*(data.setup$sp.degree + data.setup$nknots + 1) + length(x.name), 
               ncol=2*(data.setup$sp.degree + data.setup$nknots + 1) + length(x.name) );
  diag(Q) <- c( rep(0, data.setup$sp.degree+1), rep(1/sigma1^2, data.setup$nknots), rep(0, data.setup$sp.degree+1), 
                rep(1/sigma0^2, data.setup$nknots), rep(0, length(x.name) ) );
  param.hat <- c( alpha1, alpha0, beta.hat );
  
  est.var <- calc.var(param=param.hat, dat=data, dat.pois=data.poisson, coef1.name=alpha1.name, coef0.name=alpha0.name, 
                      beta.name=beta.name, out.name=out.name, xname=x.name, basis.name=var.name, tname=tname,
                      offset.name=offset.name, Qmat=Q, weight=data.setup$weight.type, ID.name=IDname );
  est.var.alpha <- est.var[c(1 : (2*(data.setup$nknots+data.setup$sp.degree+1))), c(1 : (2*(data.setup$nknots+data.setup$sp.degree+1)))];  #covariance matrix for spline curves
  zero.matrix <- matrix( 0, ncol=ncol(alpha1.var), nrow=nrow(alpha1.var) );
  est.var.alpha.noPS <- rbind( cbind(alpha1.var, zero.matrix),  cbind(zero.matrix, alpha0.var) );  
  
  ##################### weighting method ########################
  #### RACE
  tmp.race <- calc.race( tstar1 = data.setup$race.time, tstar0 = data.setup$race.time, var.alpha=est.var.alpha, var.alpha.noPS=est.var.alpha.noPS, 
                         alpha1=alpha1, alpha0=alpha0, knots1=knots1, knots0=knots0, dgr=data.setup$sp.degree );
  #print(tmp.race);
  #### SPCE
  tmp.spce <- calc.spce( select.time = data.setup$spce.time, var.alpha=est.var.alpha, var.alpha.noPS=est.var.alpha.noPS, 
                         alpha1=alpha1, alpha0=alpha0, knots1=knots1, knots0=knots0, dgr=data.setup$sp.degree );
  #print(tmp.spce);
  #### SQE
  # tmp.sqe <- calc.sqe( q = data.setup$sqe.quantile, var.alpha=est.var.alpha, var.alpha.noPS=est.var.alpha.noPS, 
  #                      alpha1=alpha1, alpha0=alpha0, knots1=knots1, knots0=knots0, dgr=data.setup$sp.degree );
  # #print( tmp.sqe );
  # ans <- rbind( tmp.race, tmp.spce, tmp.sqe );
  # ans <- cbind( data.frame( estimand = c( "RACE", "SPCE", "SQE" ) ), ans );

  ans <- rbind( tmp.race, tmp.spce );
  ans <- cbind( data.frame( estimand = c( "RACE", "SPCE" ) ), ans );
  return( list( ans=ans, alpha1=alpha1, alpha0=alpha0, knots1=knots1, knots0=knots0 ) );
}

###################################################################################################################
## 1. sigma is the smoothing paramber, which is weighting method dependent
## 2. sigma can be obtained using SAS PROC GLIMMIX on poissonized data generated by function write.poisson.data(),
##    and file "PS_HR_sigma_glimmix.sas" is prepared for GLIMMIX fitting
## 3. The following sigma values are only for illustration purpose and are not obtained from GLIMMIX fitting. 
####################################################################################################################
setup$weight.type <- "IPW";
sigma0 <- 0.2;
sigma1 <- 0.2;
res.estimand.IPW <- estimand.analysis.fun( data = data, data.setup = setup );

setup$weight.type <- "MW";
sigma0 <- 0.2
sigma1 <- 0.2
res.estimand.MW <- estimand.analysis.fun( data = data, data.setup = setup );

setup$weight.type <- "OVERLAP";
sigma0 <- 0.2;
sigma1 <- 0.2
res.estimand.OW <- estimand.analysis.fun( data = data, data.setup = setup );


setup$weight.type <- "UNWEIGHT";
sigma0 <- 0.2;
sigma1 <- 0.2
res.estimand.UNW <- estimand.analysis.fun( data = data, data.setup = setup );


##### Fitted survival estimation after PS weighting adjustment #####
# time0 <- sort( data$Y[ (data$Z == 0) & (data$delta == 1) ] );
# time1 <- sort( data$Y[ (data$Z == 1) & (data$delta == 1) ] );
# 
# surv1.IPW <- sapply( time1, calc.est.surv, knots=res.estimand.IPW$knots1, alpha=res.estimand.IPW$alpha1, dgr=setup$sp.degree );
# surv0.IPW <- sapply( time0, calc.est.surv, knots=res.estimand.IPW$knots0, alpha=res.estimand.IPW$alpha0, dgr=setup$sp.degree );
# 
# surv1.MW <- sapply( time1, calc.est.surv, knots=res.estimand.MW$knots1, alpha=res.estimand.MW$alpha1, dgr=setup$sp.degree );
# surv0.MW <- sapply( time0, calc.est.surv, knots=res.estimand.MW$knots0, alpha=res.estimand.MW$alpha0, dgr=setup$sp.degree );
# 
# surv1.OW <- sapply( time1, calc.est.surv, knots=res.estimand.OW$knots1, alpha=res.estimand.OW$alpha1, dgr=setup$sp.degree );
# surv0.OW <- sapply( time0, calc.est.surv, knots=res.estimand.OW$knots0, alpha=res.estimand.OW$alpha0, dgr=setup$sp.degree )
