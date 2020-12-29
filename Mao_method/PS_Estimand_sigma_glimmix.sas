/*****************************************************************************************************************
 * <May 25, 2018> @ Huzhang Mao
 * ---------------------------------------------------------------------------------------------------------------
 * (1) Estimate penalty parameters (sigma0, sigma1) for propensity score weighting adjusted estimand analysis 
*      using mixed poisson regression with GLIMMIX
******************************************************************************************************************/
%let wkdir = E:\PS_surv\2018_05_4Li;
%let datetag = 2018_05_25_ ;
%let wt_type = MW;
%let spdgr = 1;
%let nknots = 10;
libname BCS "&wkdir";
data _NULL_ ;
	call system('cd "&wkdir"') ; 
run;

data sigma_est;
	length z 8 sigma_est 8;
	stop;
run;

/*****************************************************
Treated group
*****************************************************/
PROC IMPORT OUT=treated
	DATAFILE="&wkdir\poisson_data\data_poisson_estimand_&wt_type._treated.csv" DBMS = CSV REPLACE;
	GETNAMES = YES;
RUN;


/*****************************
Compute offset
******************************/
data treated;
	set treated;
	lCig = log(Cig);
run;

ods output CovParms = cov_est;
proc glimmix data=treated;
	model Yig = X0 time1-time&spdgr / solution noint dist=poisson offset=lCig link=log;
	random lambda1-lambda&nknots / s G type=TOEP(1);
	weight wt;
run;

data cov_est;
	set cov_est;
	z = 1;
	sigma_est = sqrt( estimate );
	keep z sigma_est;
run;
proc append base=sigma_est data=cov_est;
run;

/*****************************************************
Control group
*****************************************************/
PROC IMPORT OUT=control
	DATAFILE="&wkdir\poisson_data\data_poisson_estimand_&wt_type._control.csv" DBMS = CSV REPLACE;
	GETNAMES = YES;
RUN;
/*****************************
Compute offset
******************************/
data control;
	set control;
	lCig = log(Cig);
run;

ods output CovParms = cov_est;
proc glimmix data=control;
	model Yig = X0 time1-time&spdgr / solution noint dist=poisson offset=lCig link=log;
	random lambda1-lambda&nknots / s G type=TOEP(1);
	weight wt;
run;

data cov_est;
	set cov_est;
	z = 0;
	sigma_est = sqrt( estimate );
	keep z sigma_est;
run;
proc append base=sigma_est data=cov_est;
run;

proc print data=sigma_est;
run;

