/**************************************************************************************************************************
 * <May 25, 2018> @ Huzhang Mao
 * ---------------------------------------------------------------------------------------------------------------
 * (1) Estimate penalty parameter (sigma) for propensity score weighting adjusted hazard ratio analysis
 *     using mixed poisson regression with GLIMMIX
****************************************************************************************************************************/
%let wkdir = E:\PS_surv\2018_05_4Li;
%let datetag = 2018_05_25_ ;
%let wt_type = MW;
%let spdgr = 1;
%let nknots = 10;
libname BCS "&wkdir";
data _NULL_ ;
	call system('cd "&wkdir"') ; 
run;

PROC IMPORT OUT=HR
	DATAFILE="&wkdir\poisson_data\data_poisson_HR_&wt_type..csv" DBMS = CSV REPLACE;
	GETNAMES = YES;
RUN;

/*****************************
Compute offset
******************************/
data HR;
	set HR;
	lCig = log(Cig);
run;

ods output CovParms = cov_est;
proc glimmix data=HR;
	model Yig = X0 Z time1-time&spdgr / solution noint dist=poisson offset=lCig link=log;
	random lambda1-lambda&nknots / s G type=TOEP(1);
	weight wt;
run;

data cov_est;
	set cov_est;
	sigma_est = sqrt( estimate );
	keep sigma_est;
run;
proc append base=sigma_est data=cov_est;
run;

proc print data=sigma_est;
run;
