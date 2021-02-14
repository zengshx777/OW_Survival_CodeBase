R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=F good_overlap=1 sample_size=150' simu_main.R Lf1.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=F good_overlap=2 sample_size=150' simu_main.R Lf2.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=F good_overlap=3 sample_size=150' simu_main.R Lf3.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=F good_overlap=1 sample_size=300' simu_main.R Mf1.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=F good_overlap=2 sample_size=300' simu_main.R Mf2.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=F good_overlap=3 sample_size=300' simu_main.R Mf3.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=F good_overlap=1 sample_size=450' simu_main.R Lf4.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=F good_overlap=2 sample_size=450' simu_main.R Lf5.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=F good_overlap=3 sample_size=450' simu_main.R Lf6.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=F good_overlap=1 sample_size=600' simu_main.R Lf7.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=F good_overlap=2 sample_size=600' simu_main.R Lf8.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=F good_overlap=3 sample_size=600' simu_main.R Lf9.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=F good_overlap=1 sample_size=750' simu_main.R Lf10.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=F good_overlap=2 sample_size=750' simu_main.R Lf11.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=F good_overlap=3 sample_size=750' simu_main.R Lf12.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=T good_overlap=1 sample_size=150' simu_main.R Lt1.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=T good_overlap=2 sample_size=150' simu_main.R Lt2.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=T good_overlap=3 sample_size=150' simu_main.R Lt3.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=T good_overlap=1 sample_size=300' simu_main.R Mf4.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=T good_overlap=2 sample_size=300' simu_main.R Mf5.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=T good_overlap=3 sample_size=300' simu_main.R Mf6ut &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=T good_overlap=1 sample_size=450' simu_main.R Lt4.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=T good_overlap=2 sample_size=450' simu_main.R Lt5.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=T good_overlap=3 sample_size=450' simu_main.R Lt6.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=T good_overlap=1 sample_size=600' simu_main.R Lt7.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=T good_overlap=2 sample_size=600' simu_main.R Lt8.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=T good_overlap=3 sample_size=600' simu_main.R Lt9.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=T good_overlap=1 sample_size=750' simu_main.R Lt10.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=T good_overlap=2 sample_size=750' simu_main.R Lt11.out &
R CMD BATCH --vanilla '--args dependent.censoring=F multi.arm=T prop.hazard=T good_overlap=3 sample_size=750' simu_main.R Lt12.out &

R CMD BATCH --vanilla '--args dependent.censoring=T multi.arm=T prop.hazard=F good_overlap=1 sample_size=300' simu_main.R df7.out &
R CMD BATCH --vanilla '--args dependent.censoring=T multi.arm=T prop.hazard=F good_overlap=2 sample_size=300' simu_main.R df8.out &
R CMD BATCH --vanilla '--args dependent.censoring=T multi.arm=T prop.hazard=F good_overlap=3 sample_size=300' simu_main.R df9.out &
R CMD BATCH --vanilla '--args dependent.censoring=T multi.arm=T prop.hazard=T good_overlap=1 sample_size=300' simu_main.R df10.out &
R CMD BATCH --vanilla '--args dependent.censoring=T multi.arm=T prop.hazard=T good_overlap=2 sample_size=300' simu_main.R df11.out &
R CMD BATCH --vanilla '--args dependent.censoring=T multi.arm=T prop.hazard=T good_overlap=3 sample_size=300' simu_main.R df12.out &
R CMD BATCH --vanilla '--args dependent.censoring=T multi.arm=T prop.hazard=F good_overlap=1 sample_size=150' simu_main.R df13.out &
R CMD BATCH --vanilla '--args dependent.censoring=T multi.arm=T prop.hazard=F good_overlap=2 sample_size=150' simu_main.R df14.out &
R CMD BATCH --vanilla '--args dependent.censoring=T multi.arm=T prop.hazard=F good_overlap=3 sample_size=150' simu_main.R df15.out &
R CMD BATCH --vanilla '--args dependent.censoring=T multi.arm=T prop.hazard=T good_overlap=1 sample_size=150' simu_main.R df16.out &
R CMD BATCH --vanilla '--args dependent.censoring=T multi.arm=T prop.hazard=T good_overlap=2 sample_size=150' simu_main.R df17.out &
R CMD BATCH --vanilla '--args dependent.censoring=T multi.arm=T prop.hazard=T good_overlap=3 sample_size=150' simu_main.R df18.out &
