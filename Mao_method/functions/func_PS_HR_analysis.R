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
library(survival);

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


calc.par.HR <- function(dat, dgr, par.start, sigma, nknots, basis.name, tname, censor.name, Dmat, 
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
  #   plot.curve: plot fitted curve? Only valid for causal estimand estimation.
  # Return
  #   Parameter estimation for Poisson model
  
  pois <- cox.to.poisson( dat=dat, dgr=dgr, nknots=nknots, basis.name=basis.name, tname=tname, censor.name=censor.name, 
                          out.name=out.name, offset.name=offset.name, wt.name=wt.name, ID.name=ID.name, time.name=time.name );
  
  pois.dat <- pois$dat.poisson;
  knots <- pois$knots;
  tmp <- newton.par.HR( dat=pois.dat, par.init=par.start, sigma.b=sigma, nknots=nknots, out.name=out.name, 
                        var.name=c(tname, basis.name), offset.name=offset.name, wt.name=wt.name, Dmat=Dmat );
  
  alpha <- tmp$par;  # coefficient estimation
  alpha.var <- solve( -tmp$hess );  # coefficient variance
  
  return( list( alpha = alpha, alpha.var = alpha.var, knots = knots, pois.dat = pois.dat ) );
}

calc.var.HR <- function( param, dat, dat.pois, coef.name, beta.name, out.name, 
                         xname, basis.name, tname, offset.name, Qmat, weight, ID.name ) {
  # calculate var-covar of all parameters given ESTIMATED propensity score 
  #
  #   param: A named vector of parameters whose variances are to be estimated
  #   beta: true coefficient in the PS model
  #   dat: data set
  #   dat.pois: poissonized data
  #   coef.name: spline coefficient names
  #   beta.name: PS model coefficient names
  #   out.come: outcome name after poissonization
  #   xname: covariate name ins the PS model
  #   basis.name: basis name in spline regression
  #   tname: treatment indicator name
  #   offset.name: offset name
  #   Qmat: Q matrix in estimating equation to exace penalized spline coefficients
  #   weight: PS weight name
  
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
    Uig <- as.matrix( this.pois.data[ , c(tname, basis.name) ] );
    Cig <- this.pois.data[ , offset.name ];
    
    phi <- calc.phi.HR(param=param, coef.name=coef.name, beta.name=beta.name, 
                    Xi=Xi, Zi=Zi, Yig=Yig, Uig=Uig, Cig=Cig, Qmat=Qmat, weight=weight, ndata=n );
    Bmat <- Bmat + outer(phi, phi);

    phi.deriv <- jacobian( calc.phi.HR, param, coef.name=coef.name, beta.name=beta.name, 
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

calc.phi.HR <- function( param, coef.name, beta.name, Xi, Zi, Yig, Uig, Cig, Qmat, weight, ndata ) {
  # calculate the score equation in estimating equation given ESTIMATED propensity score
  
  #   param: A named vector of parameters
  #   coef.name: spline coefficient names
  #   beta.name:  PS model coefficient names
  #   Xi: covariates in the PS model
  #   Zi: treatment indicator
  #   Yig: poissonized outcome
  #   Uig: spline basis
  #   Cig: offset
  #   Qmat: Q matrix in estimating equation to exact penalized coefficients
  #   weight: PS weight name
  #   ndata: number of observations in data
  
  # Return: a vector of score equation
  
  alpha <- param[ coef.name ];
  beta <- param[ beta.name ];
  
  ei <- calc.ps.Xbeta( Xmat=Xi, beta=beta );
  Wi <- calc.ps.weight( ps=ei, Zi=Zi, weight=weight );
  
  ans1.tmp <-  as.numeric( Yig - exp( Uig %*% alpha + log(Cig) ) );
  ans1 <- Wi * apply( ans1.tmp * Uig, 2, sum );
  
  ans2 <- (Zi-ei)*Xi;
  
  ans3 <- as.numeric( Qmat %*% param / ndata );
  
  ans <- c( ans1, ans2 ) - ans3;
  ans;
}

calc.HR.sp <- function( data, par.start, beta, sigma, nknots, dgr, weight, censor.name, out.name, offset.name, 
                        wt.name, basis.name, beta.name, tname, alpha.name, ID.name, time.name ) {
  # Calculate HR using spline method
  #
  # Args:
  #   data: data setup
  #   par.start: initial paramter value for Newton-Raphson fitting
  #   beta: PS model coefficients
  #   sigma: standard devation of random effect distribution
  #   nknots: number of spline knots
  #   dgr: degree of penalized spline
  #   weight: PS weight type
  #   censor.name: censor variable name
  #   out.name: poissonized outcome name
  #   offset.name: offset name
  #   wt.name: wighting variable name
  #   basis.name: spline basis name
  #   beta.name: PS model coefficient name
  #   tname: treatment indicator name
  #   alpha.name: fitted spline coefficient name
  #   ps.type: PS type, estimate or true? By default is "est"
  #   true.HR: true marginal hazard ratio
  
  # Returns
  #   A data frame with point estimation, variance, and coverage 
  
  Dmat <- cbind( matrix( 0, ncol=(dgr+length(tname)+1), nrow=nknots), diag(nknots) );  # matrix to extract penalized paratmeters in spline fitting
  
  # estimated PS for Q matrix to extract penalized parameter in estimating equation 
  Q <- matrix( 0, nrow=(dgr + nknots + length(tname) + 1) + length(beta.name), 
               ncol=(dgr + nknots + length(tname) + 1) + length(beta.name) );
  diag(Q) <- c( rep(0, dgr + length(tname) + 1), rep(1/sigma^2, nknots), rep(0, length(beta.name) ) );

  fit.sp <- calc.par.HR( dat=data, par.start=par.start, sigma=sigma, nknots=nknots, Dmat=Dmat, censor.name=censor.name, dgr=dgr,
                         out.name=out.name, offset.name=offset.name, wt.name=wt.name, basis.name=basis.name, tname=tname, 
                         ID.name=ID.name, time.name=time.name );
  data.poisson <- fit.sp$pois.dat;
  alpha <- fit.sp$alpha;
  names( alpha ) <- alpha.name;
  alpha.var <- fit.sp$alpha.var;
  knots <- fit.sp$knots;
  est.sp <- alpha[1];
  
  #### Estimated PS
  param.hat <- c( alpha, beta );
  est.var <- calc.var.HR( param=param.hat, dat=data, dat.pois=data.poisson, coef.name=alpha.name,
                          beta.name=beta.name, out.name=out.name, xname=x.name, basis.name=basis.name, tname=tname,
                          offset.name=offset.name, Qmat=Q, weight=weight, ID.name=ID.name );
  
  std.sp <- sqrt( est.var[1, 1] );
  pvalue.sp <- 2 * ( 1 - pnorm( abs( est.sp/std.sp ) ) );
  
  std.noPS <- sqrt( alpha.var[1, 1] );
  pvalue.noPS <- 2 * ( 1 - pnorm( abs( est.sp/std.noPS ) ) );
  
  ans <- data.frame( est.sp=est.sp, std.sp=std.sp, pvalue.sp=pvalue.sp, std.noPS=std.noPS, pvalue.noPS=pvalue.noPS );
  
  return( ans );
}




calc.trunc.basis <- function(ti, dgr, Knots) {
  # calculate power truncated basis
  #
  #   ti: time points
  #   dgr: degree of power truncated basis
  #   Knots: spline knots
  
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

#calc.est.cumhaz <- function (time, knots, alpha, dgr) {
# fitted cumulative hazard spline curve
#
# time: time at which cumulative hazard to be estimated
# Knots: knots for basline hazard spline
# alpha: fitted coefficients
# dgr: spline degree

#  trap.length <- 100;
#  sub.time <- seq(0, time, length=trap.length);
#  tmp1 <- sapply(sub.time, function(x) exp( calc.trunc.basis(ti=x, Knots=knots, dgr=dgr) %*% alpha) );
#  cumhaz <- trapz(x=sub.time, y=tmp1);
  
#  cumhaz;
#}


#calc.est.surv <- function( time, knots, alpha, dgr ) {
  # calculate estimated survival probability at select time
  #
  # Args
  #  time: selected time
  #  knots: spline knots
  #  alpha: fitted spline coefficients
  #
  # Return
  #  estimated survival probability
  
#  ans <- exp( -calc.est.cumhaz( time=time, knots=knots, alpha=alpha, dgr=dgr ) );
#  ans;
#}

#calc.quantile.time <- function(alpha, q, knots, dgr) {
  # caulculate survival probability at specified time
  #
  # alpha: fitted spline coefficients
  # q: survival quantile
  # knots: spline knots
  # dgr: spline degree
  
#  uniroot( calc.quantile.time.core, c(0, 1e2), q=q, knots=knots, alpha=alpha, dgr=dgr, tol=1e-8 )$root;
#}

#calc.quantile.time.core <- function(time, q, knots, alpha, dgr) {
  # calculate survival time for selected survival quantile
  #
  # x: survival time for which survival probability is q to be calculated
  # q: survival time quanitle
  # knots: spline knots
  # alpha: spline coefficients
  # dgr: spline degree
  
  # Return: time at survival quantile q
  
#  calc.est.surv( time=time, knots=knots, alpha=alpha, dgr=dgr ) - (1-q);
#}


newton.par.HR <- function(dat, par.init, sigma.b, nknots, out.name, var.name, offset.name, wt.name, Dmat, tol=1e-6, maxit=10000) {
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

