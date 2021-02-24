###########################################################################
##
## !!! Supporting functions for "main_ps_surv_test.R" without censor!!!
##
## WIN 7, R 3.3.1 @ Huzhang Mao
## 12/22/2016
###########################################################################
library(MASS);  # rmvnorm();
library(boot);  # inv.logit();
library(numDeriv);  # grad();
library(rootSolve);  # multiroot();

cox.to.poisson <- function (dat, nknots, dgr, basis.name, out.name, offset.name, tname, censor.name, wt.name, ID.name, time.name ) {
  # poissoinzation of cox data to poisson data
  #
  # Args
  #   dat: data to be poissonized 
  #   sim: simulation number
  #   nknots: number of spline knots for log basline hazard function
  #   dgr: spline degree
  #   basis.name: basis name
  #   out.name: outcome name;
  #   offset.name: offset name
  #   tname: treatment indicator name
  #   censor.name: censor indicator name
  #   wt.name: weight name
  #   ID.name: name for ID
  #   time.name: survival time name
  # Returns
  #   Poissonized data set
  
  Y <- dat[, time.name];
  Delta <- dat[ , censor.name ];
  tau <- seq( 0, max(Y), length = 50 );
  tau <- c(0, tau);

  knots <- quantile( Y[Delta==1], seq(0, 1, length = nknots+2) );
  knots <- knots[-c(1, length(knots))];  # equally spaced internal knots
  
  dat.poisson <- Reduce( "rbind", lapply( c(1:nrow(dat)), cox.to.poisson.core, dat=dat, out.name=out.name, wt.name=wt.name, 
                                        tname=tname, censor.name=censor.name, tau=tau, dgr=dgr, knots=knots, ID.name=ID.name, time.name=time.name ) );
  dat.poisson <- cbind(c(1:nrow(dat.poisson)), dat.poisson);
  
  colnames(dat.poisson) <- c("ID", ID.name, out.name, offset.name, tname, basis.name, censor.name, wt.name );

  return(list(dat.poisson=dat.poisson, knots=knots));
}

cox.to.poisson.core <- function(lnum, dat, out.name, tname, censor.name, wt.name, tau, dgr, knots, ID.name, time.name) {
  # poissonization core function
  # 
  # Args
  #   lnum: lnum-th row of survial data
  #   dat: survival data
  #   out.name: outcome name;
  #   tname: treatment indicator name  
  #   censor.name: censor indicator name
  #   wt.name: weight name
  #   tau: trapzoidal integration cut point
  #   dgr: spline degree
  #   knots: spline knots
  #   ID.name: name for ID
  #   time.name: survival time name
  # Returns
  #   Poissonized data for lnum-th subject
  
  this.dat <- dat[lnum, ];
  Ti <- this.dat[, time.name];
  PIDi <- this.dat[, ID.name];

  Deltai <- this.dat[ , censor.name];
  Zi <- this.dat[ , tname];
  wti <- this.dat[ , wt.name];
  Mi <- findInterval(Ti, tau, rightmost.closed = T)+1;
  
  g <- c(2 : Mi);
  Yig <- Deltai * ifelse(g==Mi, 1, 0);
  cig <- sapply(g, function(x) (min(tau[ifelse(x==length(tau), x, x+1)], Ti) - tau[x-1])/2);
  basis <- t( sapply(pmin(tau[g], Ti), calc.trunc.basis, dgr=dgr, Knots=knots) );

  Uig <- cbind( data.frame( PID=rep(PIDi, length(Yig) ) ), 
                as.data.frame( cbind( Yig, cig, rep(Zi, length(g)), basis, rep(Deltai, length(g)), rep(wti, length(g))) ) );
  Uig;
}


calc.spce <- function(select.time, var.alpha, var.alpha.noPS, alpha1, alpha0, knots1, knots0, dgr){
  # calcuate survival probability causal effect at selected time (SPCE)
  # Args
  #   select.time: time at which survival probability to be calculated
  #   var.alpha: variance-covariance matrix for alpha
  #   var.alpha.noPS: variance-covariance matrix for alpha without PS adjustment
  #   alpha1: fitted alpha for treated
  #   alpha0: fitted alpha for control
  #   knots1: knots for treated
  #   knots0: knots for control
  # Returns
  #   A data frame of true difference, estimated difference, bias, 
  #   estimated standard deviation, 95% confidence interval, and coverage
  
  est.surv1 <- calc.est.surv( time=select.time, knots=knots1, alpha=alpha1, dgr=dgr );  # survival probality
  surv1.deriv <- grad( calc.est.surv, x=alpha1, time=select.time, knots=knots1, dgr=dgr );  # 1st derivative of survival probability

  est.surv0 <- calc.est.surv( time=select.time, knots=knots0, alpha=alpha0, dgr=dgr );
  surv0.deriv <- grad( calc.est.surv, x=alpha0, time=select.time, knots=knots0, dgr=dgr );

  est.diff <- est.surv1 - est.surv0;
  
  surv.var <- t( c(surv1.deriv, -surv0.deriv) ) %*% var.alpha %*% c(surv1.deriv, -surv0.deriv);
  lower <- est.diff - 1.96*sqrt(surv.var);
  upper <- est.diff + 1.96*sqrt(surv.var);
  pvalue <- 2 * ( 1- pnorm( abs( est.diff/sqrt(surv.var) ) ) );
  
  surv.var.noPS <- t( c(surv1.deriv, -surv0.deriv) ) %*% var.alpha.noPS %*% c( surv1.deriv, -surv0.deriv );
  lower.noPS <- est.diff - 1.96*sqrt(surv.var.noPS);
  upper.noPS <- est.diff + 1.96*sqrt(surv.var.noPS);
  pvalue.noPS <- 2 * ( 1- pnorm( abs( est.diff/sqrt(surv.var.noPS) ) ) );
  
  ans <- data.frame( est = est.diff, std=sqrt(surv.var), lower=lower, upper=upper, pvalue=pvalue,
                     std.noPS = sqrt(surv.var.noPS), lower.noPS=lower.noPS, upper.noPS=upper.noPS, pvalue.noPS=pvalue.noPS );
  return( ans );
}



calc.sqe <- function( q, var.alpha, var.alpha.noPS, alpha1, alpha0, knots1, knots0, dgr ){
  # calcuate survival time distribution quantile causal effect at selected quantile (SQE)
  #
  # q: selected survial time dist quantile
  # var.alpha: variance-covariance matrix for alpha
  # var.alpha.noPS: variance-covariance matrix for alpha without PS adjustment
  # alpha1: fitted alpha for treated
  # alpha0: fitted alpha for control
  # knots1: knots for treated
  # knots0: knots for control
  # dgr: spline degree
  #
  # Return: A data from of true difference, estimated difference, bias, 
  # estimated standard deviation, 95% confidence interval, and coverage
  
  est.t1 <- calc.quantile.time( alpha=alpha1, q=q, knots=knots1, dgr=dgr );
  q1.deriv <- grad( calc.quantile.time, x=alpha1, q=q, knots=knots1, dgr=dgr );
  
  est.t0 <- calc.quantile.time( alpha=alpha0, q=q, knots=knots0, dgr=dgr );
  q0.deriv <- grad( calc.quantile.time, x=alpha0, q=q, knots=knots0, dgr=dgr );
  
  est.diff <- est.t1 - est.t0;
  
  q.var <- t( c(q1.deriv, -q0.deriv) ) %*% var.alpha %*% c(q1.deriv, -q0.deriv);
  lower <- est.diff - 1.96*sqrt(q.var);
  upper <- est.diff + 1.96*sqrt(q.var);
  pvalue <- 2 * ( 1- pnorm( abs( est.diff/sqrt(q.var) ) ) );
  
  q.var.noPS <- t( c(q1.deriv, -q0.deriv) ) %*% var.alpha.noPS %*% c(q1.deriv, -q0.deriv);
  lower.noPS <- est.diff - 1.96*sqrt(q.var.noPS);
  upper.noPS <- est.diff + 1.96*sqrt(q.var.noPS);
  pvalue.noPS <- 2 * ( 1- pnorm( abs( est.diff/sqrt(q.var.noPS) ) ) );
  
  ans <- data.frame( est = est.diff, std = sqrt(q.var), lower=lower, upper=upper, pvalue = pvalue,
                     std.noPS = sqrt(q.var.noPS), lower.noPS=lower.noPS, upper.noPS=upper.noPS, pvalue.noPS = pvalue.noPS );
  return( ans );
}



calc.race <- function(tstar1, tstar0, var.alpha, var.alpha.noPS, alpha1, alpha0, knots1, knots0, dgr ) {
  # calculate restricted average survival time causal effect (RACE)
  #
  # Args
  #   tstar1: retricted time for treated
  #   tstar0: retricted time for control
  #   var.alpha: variance-covariance matrix for alpha
  #   var.alpha.noPS: variance-covariance matrix for alpha without PS adjustment
  #   alpha1: fitted alpha for treated
  #   alpha0: fitted alpha for control
  #   knots1: knots for treated
  #   knots0: knots for control
  #   dgr: spline degree
  #
  # Return
  #   A data from of true difference, estimated difference, bias, 
  #   estimated standard deviation, 95% confidence interval 
  
  # exact method
  est.t1 <- calc.avg.time(time=tstar1, knots=knots1, alpha=alpha1, dgr=dgr);
  t1.deriv <- grad( calc.avg.time, x=alpha1, time=tstar1, knots=knots1, dgr=dgr );
  est.t0 <- calc.avg.time(time=tstar0, knots=knots0, alpha=alpha0, dgr=dgr);
  t0.deriv <- grad( calc.avg.time, x=alpha0, time=tstar0, knots=knots0, dgr=dgr);

  est.diff <- est.t1 - est.t0;
  
  t.var <- t( c(t1.deriv, -t0.deriv) ) %*% var.alpha %*% c(t1.deriv, -t0.deriv);
  lower <- est.diff - 1.96*sqrt(t.var);
  upper <- est.diff + 1.96*sqrt(t.var);
  pvalue <- 2 * ( 1- pnorm( abs( est.diff/sqrt(t.var) ) ) );
  
  t.var.noPS <- t( c(t1.deriv, -t0.deriv) ) %*% var.alpha.noPS %*% c(t1.deriv, -t0.deriv);
  lower.noPS <- est.diff - 1.96*sqrt(t.var.noPS);
  upper.noPS <- est.diff + 1.96*sqrt(t.var.noPS);
  pvalue.noPS <- 2 * ( 1- pnorm( abs( est.diff/sqrt(t.var.noPS) ) ) );
  
  ans <- data.frame(est = est.diff, std = sqrt(t.var), lower=lower, upper=upper, pvalue=pvalue, 
                    std.noPS = sqrt(t.var.noPS), lower.noPS=lower.noPS, upper.noPS=upper.noPS, pvalue.noPS = pvalue.noPS );
  return( ans );
}


calc.ace <- function(t1, t0, var.alpha, var.alpha.noPS, alpha1, alpha0, knots1, knots0, dgr) {
  # calculate average survival time causal effect (ACE) using numeric derivative
  #
  # t1: max time for treated 
  # t0: max time for control
  # var.alpha: variance-covariance matrix for alpha
  # var.alpha.noPS: variance-covariance matrix for alpha without PS adjustment
  # alpha1: fitted alpha for treated
  # alpha0: fitted alpha for control
  # knots1: knots for treated
  # knots0: knots for control
  # true.diff: true difference
  # dgr: spline degree
  #
  # Return: A data from of true difference, estimated difference, bias, 
  # estimated standard deviation, 95% confidence interval 
  
  ## exact method
  est.t1 <- calc.avg.time(time=t1, knots=knots1, alpha=alpha1, dgr=dgr );
  t1.deriv <- grad(calc.avg.time, x=alpha1, time=t1, knots=knots1, dgr=dgr );
  est.t0 <- calc.avg.time( time=t0, knots=knots0, alpha=alpha0, dgr=dgr );
  t0.deriv <- grad( calc.avg.time, x=alpha0, time=t0, knots=knots0, dgr=dgr );

  est.diff <- est.t1 - est.t0;
  t.var <- t( c(t1.deriv, -t0.deriv) ) %*% var.alpha %*% c(t1.deriv, -t0.deriv);
  lower <- est.diff - 1.96*sqrt(t.var);
  upper <- est.diff + 1.96*sqrt(t.var);
  pvalue <- 2 * ( 1- pnorm( abs( est.diff/sqrt(t.var) ) ) );
  
  t.var.noPS <- t( c(t1.deriv, -t0.deriv) ) %*% var.alpha.noPS %*% c(t1.deriv, -t0.deriv);
  lower.noPS <- est.diff - 1.96*sqrt(t.var.noPS);
  upper.noPS <- est.diff + 1.96*sqrt(t.var.noPS);
  pvalue.noPS <- 2 * ( 1- pnorm( abs( est.diff/sqrt(t.var.noPS) ) ) );
  
  ans <- data.frame( est = est.diff, std = sqrt(t.var), lower=lower, upper=upper, pvalue = pvalue,
                     std.noPS = sqrt(t.var.noPS), lower.noPS=lower.noPS, upper.noPS=upper.noPS, pvalue.noPS = pvalue.noPS );
  return( ans );
}


calc.par <- function(dat, dgr, par.start, sigma, nknots, basis.name, tname, censor.name, Dmat, weight, 
                     out.name, offset.name, wt.name, ID.name, time.name ) {
  # function to calculate fitted parameter for HR estimation
  #
  #   dat: data to be used
  #   dgr: spline degree
  #   par.start: initial value for N-R parameter calculation
  #   sigma: sigma parameter
  #   basis.name: spline basis name
  #   Dmat: D matrix
  #   weight: weight name
  #   out.name: outcome variable name
  #   offset.name: offset variable name
  #   censor.name: name of censor indicator
  #   tname: treatment indicator name
  #   wt.name: weight variable name
  #   ID.name: name for ID
  #   time.name: survival time name
  #   plot.curve: plot fitted curve? Only valid for causal estimand estimation.
  # Return
  #   Parameter estimation for Poisson model
  
  pois <- cox.to.poisson( dat=dat, dgr=dgr, nknots=nknots, basis.name=basis.name, tname=tname, censor.name=censor.name, 
                          out.name=out.name, offset.name=offset.name, wt.name=wt.name, ID.name=ID.name, time.name=time.name );
  
  pois.dat <- pois$dat.poisson;
  knots <- pois$knots;
  tmp <- newton.par(dat=pois.dat, par.init=par.start, sigma.b=sigma, nknots=nknots, 
                      out.name=out.name, var.name=basis.name, offset.name=offset.name, wt.name=wt.name, Dmat=Dmat);
  
  alpha <- tmp$par;  # coefficient estimation
  alpha.var <- solve( -tmp$hess );  # coefficient variance

  return( list( alpha = alpha, alpha.var = alpha.var, knots = knots, pois.dat = pois.dat ) );
}

calc.var <- function(param, dat, dat.pois, coef1.name, coef0.name, beta.name, out.name, 
                         xname, basis.name, tname, offset.name, Qmat, weight, ID.name ) {
  # calculate var-covar of all parameters
  #
  #   param: A named vector of parameters
  #   dat: data set
  #   dat.pois: poissonized data
  #   coef1.name, coef0.name: coefficient name for the treated  and control groups respectively
  #   beta.name: coefficient name for the PS model
  #   out.come: outcome name after poissonization
  #   xname: covariate name ins the PS model
  #   basis.name: basis name in spline regression
  #   tname: treatment indicator name
  #   offset.name: offset name
  #   Qmat: Q matrix in estimating equation to exace proper parameters
  #   weight: PS weight name
  #   ID.name: name for ID
  
  # Return: A variance-covariance matrix
  
  n <- nrow(dat);
  Amat <- Bmat <- 0;
  for ( i in 1 : n ) {
    
    this.data <- dat[ i, ];  #data for each subject
    this.ID.name <- ( dat[ , ID.name] )[i];
    this.pois.data <- dat.pois[ dat.pois[ , ID.name ] == this.ID.name, ];  # poissonized data for each PID
    
    Xi <- t( this.data[ , xname ] ); 
    Zi <- as.numeric( this.data[, tname ] );
    Yig <- this.pois.data[ , out.name ];
    Uig <- as.matrix( this.pois.data[ , basis.name ] );
    Cig <- this.pois.data[ , offset.name ];
    
    phi <- calc.phi(param=param, coef1.name=coef1.name, coef0.name=coef0.name, beta.name=beta.name, 
                        Xi=Xi, Zi=Zi, Yig=Yig, Uig=Uig, Cig=Cig, Qmat=Qmat, weight=weight, ndata=n );
    Bmat <- Bmat + outer(phi, phi);
    
    phi.deriv <- jacobian( calc.phi, param, coef1.name=coef1.name, coef0.name=coef0.name, beta.name=beta.name, 
                           Xi=Xi, Zi=Zi, Yig=Yig, Uig=Uig, Cig=Cig, Qmat=Qmat, weight=weight, ndata=n );
    Amat <- Amat + phi.deriv;
  }

  Amat <- Amat/n;
  Bmat <- Bmat/n;
  Amat.inv <- solve(Amat);
  var.mat <- ( Amat.inv %*% Bmat %*% t(Amat.inv) ) / n;
  
  colnames(var.mat) <- rownames(var.mat) <- names(param);
  
  return(var.mat);
}

calc.phi <- function( param, coef1.name, coef0.name, beta.name, Xi, Zi, Yig, Uig, Cig, Qmat, weight, ndata ) {
  # calculate the phi in estimating equation
  
  #   param: A named vector of parameters
  #   coef1.name, coef0.name: coefficient name for the treated  and control groups respectively
  #   beta.name: coefficient name for the PS model
  #   out.come: outcome name after poissonization
  #   xname: covariate name ins the PS model
  #   basis.name: basis name in spline regression
  #   tname: treatment indicator name
  #   offset.name: offset name
  #   Qmat: Q matrix in estimating equation to exace proper parameters
  #   weight: PS weight name
  #   ndata: number of observations in data
  
  # Return: a vector

  alpha1 <- param[ coef1.name ];
  alpha0 <- param[ coef0.name ];
  beta <- param[ beta.name ];
  
  #print( Xi );
  
  ei <- calc.ps.Xbeta( Xmat=Xi, beta=beta );
  
  #print( ei );
  
  Wi <- calc.ps.weight( ps=ei, Zi=Zi, weight=weight );
  
  ans1.tmp <-  as.numeric( Yig - exp( Uig %*% alpha1 + log(Cig) ) );
  ans1 <- Zi * Wi * apply( ans1.tmp * Uig, 2, sum );
  
  ans2.tmp <- as.numeric( Yig - exp( Uig %*% alpha0 + log(Cig) ) );
  ans2 <- (1-Zi) * Wi * apply( ans2.tmp * Uig, 2, sum );
  
  ans3 <- (Zi-ei)*Xi;
  
  ans4 <- as.numeric( Qmat %*% param / ndata );
  
  ans <- c( ans1, ans2, ans3 ) - ans4;
  ans;
}

#set.varcovar.name <- function(prefix, ncoef) {
  # set dim name for var-covar matrix
#  tmp <- sapply(c(1:ncoef), function(x) paste(prefix, x, sep=".")); 

#  return(tmp);
#}

calc.trunc.basis <- function(ti, dgr, Knots) {
  # calculate power truncated basis
  #
  #   ti: time points
  #   dgr: degree of power truncated basis
  #   Knots: spline knots
  
  #print( dgr );
  
  tmp0 <- sapply( c( 0 : dgr), function(x) ti^x );
  tmp1 <- (ti-Knots);
  tmp2 <- tmp1^dgr * (tmp1 >= 0);
  ans <- c( tmp0, tmp2);
  ans;  
}

set.basis.name <- function (nbasis, prefix) {
  # Set basis names 
  #
  # Args:
  # n.basis: number of basis used
  #
  # Returns:
  # A character vector of basis names with defined prefix
  
  ans <- vector(mode = "character", length = 0)
  for (i in 1 : nbasis) {
    ans <- c(ans, paste(prefix, i, sep = ""))
  }
  ans
}


calc.est.cumhaz <- function (time, knots, alpha, dgr) {
# fitted cumulative hazard spline curve
#
# time: time at which cumulative hazard to be estimated
# Knots: knots for basline hazard spline
# alpha: fitted coefficients
# dgr: spline degree
  
  #print( dgr );

  trap.length <- 100;
  sub.time <- seq(0, time, length=trap.length);
  tmp1 <- sapply(sub.time, function(x) exp( calc.trunc.basis(ti=x, Knots=knots, dgr=dgr) %*% alpha) );
  cumhaz <- trapz(x=sub.time, y=tmp1);
  
  cumhaz;
}



calc.est.surv <- function( time, knots, alpha, dgr ) {
  # calculate estimated survival probability at select time
  #
  # Args
  #  time: selected time
  #  knots: spline knots
  #  alpha: fitted spline coefficients
  #
  # Return
  #  estimated survival probability

  ans <- exp( -calc.est.cumhaz( time=time, knots=knots, alpha=alpha, dgr=dgr ) );
  ans;
}


calc.quantile.time <- function(alpha, q, knots, dgr) {
  # caulculate survival probability at specified time
  #
  # alpha: fitted spline coefficients
  # q: survival quantile
  # knots: spline knots
  # dgr: spline degree
  
  uniroot( calc.quantile.time.core, c(0, 1e10), q=q, knots=knots, alpha=alpha, dgr=dgr, tol=1e-8 )$root;
}


calc.quantile.time.core <- function(time, q, knots, alpha, dgr) {
  # calculate survival time for selected survival quantile
  #
  # x: survival time for which survival probability is q to be calculated
  # q: survival time quanitle
  # knots: spline knots
  # alpha: spline coefficients
  # dgr: spline degree
  
  # Return: time at survival quantile q
  
  calc.est.surv( time=time, knots=knots, alpha=alpha, dgr=dgr ) - (1-q);
}



newton.par <- function(dat, par.init, sigma.b, nknots, out.name, var.name, offset.name, wt.name, Dmat, tol=1e-6, maxit=10000) {
  # Newton-Raphson maximization for penalized weighted likelihood
  
  # Args:
  # dat.pois: poissonized data.
  # par.init: initial parameter value.
  # sigma.b: value for sigma.b
  # nknots: number of internal knots
  # out.name: outcome name
  # var.name: covariate names
  # offset.name: offset name
  # Dmat: D matrix
  # nknots: number of internal knots

  par.old <- par.init;

  itnum <- 1;
  while (itnum <= maxit) {
    grad <- Reduce( "+", lapply(c(1:nrow(dat)), calc.newton.grad, dat.pois=dat, alpha=par.old, outname=out.name, 
                                varname=var.name, offsetname=offset.name, wtname=wt.name) ) - t(Dmat)%*%Dmat%*%par.old/sigma.b^2;
    
    hess <- Reduce( "+", lapply(c(1:nrow(dat)), calc.newton.hess, dat.pois=dat, alpha=par.old, varname=var.name,
                                offsetname=offset.name, wtname=wt.name) ) - t(Dmat)%*%Dmat/sigma.b^2;

    par.new <- par.old - solve(hess) %*% grad;  # update
    
    converge <- sum(abs(par.new - par.old))/(sum(abs(par.old)));  # relative converge criteria   
    if (converge <= tol) {
      cat("N-R converge successfully at\t", itnum, "-th iteration\n");
      return( list(par = par.new, hess = hess) );
    } else {
      par.old <- par.new;
      itnum <- itnum + 1;
    }
  }
  stop("N-R gamma converge failed");
}


calc.newton.grad <- function(lnum, dat.pois, alpha, outname, varname, offsetname, wtname) {
  # calculate gradient vector for Newton-Raphson method
  
  # lum: line number for each data line 
  # dat.pois: poisson data.
  # alpha: paramete value.
  # outname: outcome name.
  # varname: fixed and random effects name.
  # offsetname: offset name
  # wtname: weight name
  
  Yig <- as.numeric( dat.pois[lnum, ][outname] );
  Uig <- as.numeric( dat.pois[lnum, ][varname] );
  cig <- as.numeric( dat.pois[lnum, ][offsetname] );
  tmp1 <- as.numeric( exp( Uig %*% alpha + log(cig) ) );

  ans <- as.numeric(dat.pois[lnum, ][wtname]) * ( (Yig - tmp1)*colVec(Uig) ); 

  ans;
}

calc.newton.hess <- function(lnum, dat.pois, alpha, varname, offsetname, wtname) {
  # calculate hessian matrix for Newton-Raphson method
  
  # lum: line number for each data line 
  # dat.pois: poisson data.
  # alpha: paramete value.
  # varname: fixed and random effects name.
  # offsetname: offset name
  # wtname: weight name
  
  Uig <- as.numeric( as.vector( dat.pois[lnum, ][varname] ) );
  cig <- as.numeric( dat.pois[lnum, ][offsetname] );
  tmp1 <- as.numeric( exp( Uig %*% alpha + log(cig) ) );
  ans <- as.numeric(dat.pois[lnum, ][wtname]) * ( -tmp1 * colVec(Uig) %*% Uig );
  
  ans;
}


trapz <- function(y, x){
  idx <- 2:length(x);
  return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)  
}

calc.avg.time <- function(time, knots, alpha, dgr) {
  # calculate averaged survial time
  #
  # time: pre-selected time
  # knots: spline knots
  # alpha: spline coefficients
  
  time.seq <- seq(0, time, length=100);
  tmp <- sapply(time.seq, calc.integ1, knots=knots, alpha=alpha, dgr=dgr);
  trapz( y=exp(-tmp), x=time.seq );
}


calc.integ1 <- function(time, knots, alpha, dgr) {
  # calcualate integration exp( U(t)*alpha ) via trapezoidal rule
  #
  # time: pre-specified time
  # knots: time knots
  # alpha: fitted alpha
  # dgr: spline degree
  #
  # Return: integrated value of length 1
  
  time.range <- seq(0, time, length=100);
  #Ut <- t( sapply(time.range, calc.linear.basis, Knots=knots) );
  Ut <- t( sapply(time.range, calc.trunc.basis, Knots=knots, dgr=dgr) );
  ans <- trapz(y=exp(Ut %*% alpha), x=time.range);
  
  ans;
}
