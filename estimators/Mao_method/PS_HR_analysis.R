####################################################################################################
##
## Example code for propensity score weighting adjusted analysis for hazard ratio estimation 
## 05/25/2018 @ Huzhang Mao
## 
##  WIN 7, R version 3.3.3 (2017-03-06)
####################################################################################################
rm( list=ls( all = T ) );
date.tag <- "2018_05_25_";

set.seed(123);
setwd( "E:\\PS_surv\\2018_05_4Li" );
source( "functions//func_PS_HR_analysis.R" );
source( "functions//GIPW_function_omega.R" );

##################################################################################
## Data generation:
## Note
##  1. Treatment variable must be 0 or 1;
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

#####################################################################
## setup for spline analysis
#####################################################################
data.setup <- list( nsample = 0,
                    sp.degree = 1,
                    nknots = 10,
                    weight.type = ' '
                  );

IDname <- "PID";
prefix.knots <- "lambda";
prefix.fix <- "time";
knots.name <- paste(set.basis.name(nbasis = data.setup$nknots, prefix=prefix.knots), sep=',');
var.name <- c( "X0", t( sapply( c(1 : data.setup$sp.degree), function(x) paste( prefix.fix, x, sep='') ) ), knots.name);
out.name <- "Yig";
x.name <- c( "X0", X4PS );
offset.name <- "Cig";
wt.name <- "wt";
censor.name <- "delta";
tname <- "Z";
time.name <- "Y";
alpha.name <- c('alpha.Z', sapply( 1 : ( data.setup$nknots + data.setup$sp.degree+1), function(x) paste( 'alpha.', x, sep='') ) );
beta.name <- sapply( c( 0 : (length(x.name)-1) ) , function(x) paste( "beta", x, sep='') );
ps.form <- as.formula( paste( "Z ~ ", paste( X4PS, collapse = " + " ), sep="" ) );
surv.form <- as.formula( 'Surv(Y, delta) ~ Z' );
D <- cbind(matrix( 0, ncol = ( data.setup$sp.degree + length(tname) + 1 ), nrow = data.setup$nknots ), diag( data.setup$nknots ) );
par.start <- rep(0, length = ( data.setup$sp.degree + data.setup$nknots + length(tname)+1) );


########### Poissonization data to calculate penalty parameter (sigma) ################
write.poisson.data <- function( data, weight.type ) {
  # poissonization for to calculate pentalty parameter (sigma)
  #   data: data to be used
  #   weight.type: PS weighting method
  
  out.ps <- ps.model(dat=data, form=ps.form );
  data$ps <- out.ps$ps.hat;  # estimated ps
  
  if(weight.type!="UNWEIGHT"){
    omega <- sapply( data$ps, calc.omega, weight = weight.type );
  }else{
    omega = data$ps * data$Z + (1-data$ps) * (1-data$Z)
  }
  data$wt <- omega / ( data$ps * data$Z + (1-data$ps) * (1-data$Z) );
  
  tmp <- cox.to.poisson( dat=data, nknots=data.setup$nknots, dgr=data.setup$sp.degree, basis.name=var.name, out.name=out.name, 
                         offset.name=offset.name, tname=tname, censor.name=censor.name, wt.name=wt.name, ID.name=IDname, time.name=time.name );
  write.csv( tmp$dat.poisson, paste( "poisson_data//data_poisson_HR_", weight.type, ".csv", sep="" ), row.names = F );
}

# IPW
write.poisson.data( data = data, weight.type="IPW" );
# MW
write.poisson.data( data = data, weight.type="MW" );
# OW
write.poisson.data( data = data, weight.type="OVERLAP" );
# UNWEIGHT
write.poisson.data( data = data, weight.type="UNWEIGHT" );


HR.analysis.fun <- function( data, setup ) {
  # Analysis of hazard ratio using the simulated data
  #   data: data set
  #   setup: analysis setup

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

  #### Weighted Spline Cox regression ####
  res.sp <- calc.HR.sp( data=data, par.start=par.start, beta=beta.hat, sigma=sigma, nknots=setup$nknots, dgr=setup$sp.degree, 
                        weight=setup$weight.type, censor.name=censor.name, out.name=out.name, offset.name=offset.name, wt.name=wt.name, 
                        basis.name=var.name, beta.name=beta.name, tname=tname, alpha.name=alpha.name, ID.name=IDname, time.name=time.name );
  
  ans <- data.frame( est = res.sp[, "est.sp"], std = res.sp[, "std.sp"], lower = res.sp[, "est.sp"] - 1.96*res.sp[, "std.sp"], 
                      upper = res.sp[, "est.sp"] + 1.96*res.sp[, "std.sp"], pvalue = res.sp[, "pvalue.sp"] );

  return( ans );
}


###################################################################################################################
## 1. sigma is the smoothing paramber, which is weighting method dependent
## 2. sigma can be obtained using SAS PROC GLIMMIX on poissonized data generated by function write.poisson.data(),
##    and file "PS_HR_sigma_glimmix.sas" is prepared for GLIMMIX fitting
## 3. The following sigma values are only for illustration purpose and are not obtained from GLIMMIX fitting. 
####################################################################################################################
data.setup$weight.type <- "IPW";
sigma <- 0.1; # mock value
res.HR.IPW <- HR.analysis.fun( data = data, setup = data.setup );


data.setup$weight.type <- "MW";
sigma <- 0.2;  # mock value
res.HR.MW <- HR.analysis.fun( data = data, setup = data.setup );

data.setup$weight.type <- "OVERLAP";
sigma <- 0.25; # mock value
res.HR.OW <- HR.analysis.fun( data = data, setup = data.setup );

data.setup$weight.type <- "UNWEIGHT";
sigma <- 0.25; # mock value
res.HR.UW <- HR.analysis.fun( data = data, setup = data.setup )
