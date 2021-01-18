# R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=T good_overlap=1 sample_size=150' simu_main.R M1.out &
# R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=T good_overlap=2 sample_size=150' simu_main.R M2.out &
# R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=T good_overlap=3 sample_size=150' simu_main.R M3.out &
# R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=T good_overlap=1 sample_size=250' simu_main.R M4.out &
# R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=T good_overlap=2 sample_size=250' simu_main.R M5.out &
# R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=T good_overlap=3 sample_size=250' simu_main.R M6.out &
# # R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=T good_overlap=1 sample_size=400' simu_main.R R7.out &
# # R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=T good_overlap=2 sample_size=400' simu_main.R R8.out &
# # R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=T good_overlap=3 sample_size=400' simu_main.R R9.out &
# # R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=T good_overlap=1 sample_size=500' simu_main.R R10.out &
# # R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=T good_overlap=2 sample_size=500' simu_main.R R11.out &
# # R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=T good_overlap=3 sample_size=500' simu_main.R R12.out &
# # R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=F good_overlap=1 sample_size=250' simu_main.R R13.out &
# # R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=F good_overlap=2 sample_size=250' simu_main.R R14.out &
# # R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=F good_overlap=3 sample_size=250' simu_main.R R15.out &
# 
# R CMD BATCH --vanilla '--args multi.arm=T prop.hazard=T good_overlap=1 sample_size=250' simu_main.R R16.out &
# R CMD BATCH --vanilla '--args multi.arm=T prop.hazard=T good_overlap=2 sample_size=250' simu_main.R R17.out &
# R CMD BATCH --vanilla '--args multi.arm=T prop.hazard=T good_overlap=3 sample_size=250' simu_main.R R18.out &
# R CMD BATCH --vanilla '--args multi.arm=T prop.hazard=F good_overlap=1 sample_size=250' simu_main.R R19.out &
# R CMD BATCH --vanilla '--args multi.arm=T prop.hazard=F good_overlap=2 sample_size=250' simu_main.R R20.out &
# R CMD BATCH --vanilla '--args multi.arm=T prop.hazard=F good_overlap=3 sample_size=250' simu_main.R R21.out &

R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=T good_overlap=1 sample_size=250 dependent.censoring=T' simu_main.R M1.out &
R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=T good_overlap=2 sample_size=250 dependent.censoring=T' simu_main.R M2.out &
R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=T good_overlap=3 sample_size=250 dependent.censoring=T' simu_main.R M3.out &
R CMD BATCH --vanilla '--args multi.arm=T prop.hazard=T good_overlap=1 sample_size=250 dependent.censoring=T' simu_main.R M4.out &
R CMD BATCH --vanilla '--args multi.arm=T prop.hazard=T good_overlap=2 sample_size=250 dependent.censoring=T' simu_main.R M5.out &
R CMD BATCH --vanilla '--args multi.arm=T prop.hazard=T good_overlap=3 sample_size=250 dependent.censoring=T' simu_main.R M6.out &
R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=F good_overlap=1 sample_size=250 dependent.censoring=T' simu_main.R M7.out &
R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=F good_overlap=2 sample_size=250 dependent.censoring=T' simu_main.R M8.out &
R CMD BATCH --vanilla '--args multi.arm=F prop.hazard=F good_overlap=3 sample_size=250 dependent.censoring=T' simu_main.R M9.out &
