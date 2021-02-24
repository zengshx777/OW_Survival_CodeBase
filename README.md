# OW_Survival_CodeBase
### R Scripts in  folder ```estimators```

* ```fast_pseudo_calculation.R``` : Function for calculating pseudo-observations, which is faster than Rpackage **pseudo**.
* ```PSW_pseudo.R``` : Function for IPW-PO and OW-PO.
* ```cox_model.R```  : Functions for estimator Cox and Cox-IPW.
* ```pseudo_G.R```  : Function for PO-UNADJ and PO-G.
* ```Mao_Method_func.R```  : Function for estimators in Mao's paper.
* ```AIPW_pseudo.R```  : Function for AIPW and AOW.

### R Scripts in folder ```simulation```

* ```simu_main.R``` : Main script for running simulations.
* ```simu_utils.R``` :  Utility functions for simulations.
* ```simu_data_gen.R```  : Utility functions for generating simulated data.
* ```simu_exe.sh```  : Bash script to run simulations in all settings.

### R Scripts in  folder```data_application```

* ```data_preprocessing.R``` : Cleaning for data application.
* ```data_application.R``` :  Analyze function for data application.

### Run simulations

To run the simulation in the paper, you can run following the command and set the ```simulation``` as the working directory.
```
git clone https://github.com/zengshx777/OW\Survival_CodeBase

R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=F good_overlap=1 sample_size=150' simu_main.R R1.out
```


where ```dependent.censoring``` controls whether the censoring is independent of the covariates; ```multi.arm``` controls the number of arms in the data (```T``` for J=3, ```F```} for J=2); ```prop.hazard``` controls whether the proportional hazard assumption is correct; ```good_overlap``` control the overlap conditions (```1``` for RCT, ```2``` for good overlap, ```3``` for poor overlap); ```sample_size``` control the sample size. One simple way to run many simulations in different settings in parallel is to run the ```simu_exe.sh``` directly (you can customize the scenario in this file). The current ```simu_main.R``` will run all estimators mentioned in the paper by default, which might be time-consuming. You can comment out certain estimators to speed up.

The results will be saved in the folder ```simulation_results``` To output the similar Figures and Tables in the paper, please refer to the scripts in folder ```output_utils```

The NCDB data used in the case study is publicly available upon approval of the NCDB Participant User File application.
