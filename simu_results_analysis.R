#  }
#}
# ow_mse=sqrt(ow_mse)
# ipw_mse=sqrt(ipw_mse)
# ow_mse_table=matrix(ow_mse,ncol=2,byrow = T)
# ow_mse_table=cbind(matrix(ow_mse_table[,1],ncol=3),matrix(ow_mse_table[,2],ncol=3))
# 
# ipw_mse_table=matrix(ipw_mse,ncol=2,byrow = T)
# ipw_mse_table=cbind(matrix(ipw_mse_table[,1],ncol=3),matrix(ipw_mse_table[,2],ncol=3))
# 
# ow_bias_table=matrix(ow_bias,ncol=2,byrow = T)
# ow_bias_table=cbind(matrix(ow_bias_table[,1],ncol=3),matrix(ow_bias_table[,2],ncol=3))
# 
# 
# ipw_bias_table=matrix(ipw_bias,ncol=2,byrow = T)
# ipw_bias_table=cbind(matrix(ipw_bias_table[,1],ncol=3),matrix(ipw_bias_table[,2],ncol=3))
# 
# ow_coverage_table=matrix(ow_coverage,ncol=2,byrow = T)
# ow_coverage_table=cbind(matrix(ow_coverage_table[,1],ncol=3),matrix(ow_coverage_table[,2],ncol=3))
# 
# ipw_coverage_table=matrix(ipw_coverage,ncol=2,byrow = T)
# ipw_coverage_table=cbind(matrix(ipw_coverage_table[,1],ncol=3),matrix(ipw_coverage_table[,2],ncol=3))
# 
# 
# plot_size=2.5
# pdf("Summary.pdf",height=3.3*plot_size,width=3*plot_size)
# m <- matrix(c(1,2,3,4,5,6,7,8,9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)
# layout(mat = m,heights = c(0.5,0.5,0.5,0.15))
# plot(sample_size_grid,ow_mse_table[,1],type='o',
#      xlab="Sample Size",
#      ylab='RMSE',
#      main="RMSE comparison, ASCE",
#      ylim=range(ow_mse_table,ipw_mse_table),
#      col="black",pch=15)
# lines(sample_size_grid,ow_mse_table[,4],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_mse_table[,1],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_mse_table[,4],type='o',col="blue",pch=19)
# 
# plot(sample_size_grid,ow_bias_table[,1],type='o',
#      xlab="Sample Size",
#      ylab='Absolute Bias',
#      main="Absolute bias comparison, ASCE",
#      ylim=range(ow_bias_table,ipw_bias_table),
#      col="black",pch=15)
# lines(sample_size_grid,ow_bias_table[,4],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_bias_table[,1],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_bias_table[,4],type='o',col="blue",pch=19)
# 
# 
# plot(sample_size_grid,ow_coverage_table[,1],type='o',
#      xlab="Sample Size",
#      ylab='Coverage Rate',
#      main="Coverage Rate for 95% CI, ASCE",
#      ylim=range(ow_coverage_table,ipw_coverage_table),
#      col="black",pch=15)
# abline(h=0.95,lty=2)
# lines(sample_size_grid,ow_coverage_table[,4],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_coverage_table[,1],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_coverage_table[,4],type='o',col="blue",pch=19)
# 
# plot(sample_size_grid,ow_mse_table[,3],type='o',
#      xlab="Sample Size",
#      ylab='RMSE',
#      main="RMSE comparison, RACE",
#      ylim=range(ow_mse_table,ipw_mse_table),
#      col="black",pch=15)
# lines(sample_size_grid,ow_mse_table[,6],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_mse_table[,3],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_mse_table[,6],type='o',col="blue",pch=19)
# 
# plot(sample_size_grid,ow_bias_table[,3],type='o',
#      xlab="Sample Size",
#      ylab='Absolute Bias',
#      main="Absolute bias comparison, RACE",
#      ylim=range(ow_bias_table,ipw_bias_table),
#      col="black",pch=15)
# lines(sample_size_grid,ow_bias_table[,6],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_bias_table[,3],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_bias_table[,6],type='o',col="blue",pch=19)
# 
# 
# plot(sample_size_grid,ow_coverage_table[,3],type='o',
#      xlab="Sample Size",
#      ylab='Coverage Rate',
#      main="Coverage Rate for 95% CI, RACE",
#      ylim=range(ow_coverage_table,ipw_coverage_table),
#      col="black",pch=15)
# abline(h=0.95,lty=2)
# lines(sample_size_grid,ow_coverage_table[,6],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_coverage_table[,3],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_coverage_table[,6],type='o',col="blue",pch=19)
# 
# 
# plot(sample_size_grid,ow_mse_table[,2],type='o',
#      xlab="Sample Size",
#      ylab='RMSE',
#      main="RMSE comparison, SPCE",
#      ylim=range(ow_mse_table[,2],ipw_mse_table[,5]),
#      col="black",pch=15)
# lines(sample_size_grid,ow_mse_table[,5],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_mse_table[,2],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_mse_table[,5],type='o',col="blue",pch=19)
# 
# plot(sample_size_grid,ow_bias_table[,2],type='o',
#      xlab="Sample Size",
#      ylab='Absolute Bias',
#      main="Absolute bias comparison, SPCE",
#      ylim=range(ow_bias_table[,2],ipw_bias_table[,5]),
#      col="black",pch=15)
# lines(sample_size_grid,ow_bias_table[,5],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_bias_table[,2],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_bias_table[,5],type='o',col="blue",pch=19)
# 
# 
# plot(sample_size_grid,ow_coverage_table[,2],type='o',
#      xlab="Sample Size",
#      ylab='Coverage Rate',
#      main="Coverage Rate for 95% CI, SPCE",
#      ylim=range(ow_coverage_table,ipw_coverage_table),
#      col="black",pch=15)
# abline(h=0.95,lty=2)
# lines(sample_size_grid,ow_coverage_table[,5],type='o',col="red",pch=17)
# lines(sample_size_grid,ipw_coverage_table[,2],type='o',col="purple",pch=18)
# lines(sample_size_grid,ipw_coverage_table[,5],type='o',col="blue",pch=19)
# 
# par(mar = c(0.1,0.1,1,0.1))
# plot(1, type = "n", axes=FALSE, xlab="", ylab="")
# legend("top",inset=0,title="Method, Overlaps Condition", col=c("black","red","purple","blue"),lwd=1.5,
#        lty=1,pch=c(15,17,18,19),legend=c("OW, Good Overlaps","OW, Poor Overlaps","IPW, Good Overlaps","IPW, Poor Overlaps"),horiz=TRUE)
# dev.off()
# coverage_rate_calc(ow_est[,1],ow_se[,1],mean(true_est_ow[,1]))
# coverage_rate_calc(ow_est[,2],ow_se[,2],mean(true_est_ow[,2]))
# coverage_rate_calc(ow_est[,3],ow_se[,3],mean(true_est_ow[,3]))
# 
# coverage_rate_calc(ow_est[,1],ow_se[,1],mean(true_est_ow[,1]))
# coverage_rate_calc(ow_est[,2],ow_se[,2],mean(true_est_ow[,2]))
# coverage_rate_calc(ow_est[,3],ow_se[,3],mean(true_est_ow[,3]))
# 
# sqrt(mean((unadj_est[,1]-true_est_finite[,1])^2))/mean(abs(true_est_finite[,1]))
# 
# sqrt(mean((ipw_est[,1]-true_est_finite[,1])^2))/mean(abs(true_est_finite[,1]))
# sqrt(mean((ipw_est_by_group[,1]-true_est_finite[,1])^2))/mean(abs(true_est_finite[,1]))
# sqrt(mean((ipw_est_mao[,1]-true_est_finite[,1])^2))/mean(abs(true_est_finite[,1]))
# 
# sqrt(mean((ow_est[,1]-true_est_ow_finite[,1])^2))/mean(abs(true_est_ow_finite[,1]))
# sqrt(mean((ow_est_by_group[,1]-true_est_ow_finite[,1])^2))/mean(abs(true_est_ow_finite[,1]))
# sqrt(mean((ow_est_mao[,1]-true_est_ow_finite[,1])^2))/mean(abs(true_est_ow_finite[,1]))
# 
# sqrt(mean((unadj_est[,2]-true_est_finite[,2])^2))/mean(abs(true_est_finite[,2]))
# 
# sqrt(mean((ipw_est[,2]-true_est_finite[,2])^2))/mean(abs(true_est_finite[,2]))
# sqrt(mean((ipw_est_by_group[,2]-true_est_finite[,2])^2))/mean(abs(true_est_finite[,2]))
# sqrt(mean((ipw_est_mao[,2]-true_est_finite[,2])^2))/mean(abs(true_est_finite[,2]))
# 
# 
# sqrt(mean((ow_est[,2]-true_est_ow_finite[,2])^2))/mean(abs(true_est_ow_finite[,2]))
# sqrt(mean((ow_est_by_group[,2]-true_est_ow_finite[,2])^2))/mean(abs(true_est_ow_finite[,2]))
# sqrt(mean((ow_est_mao[,2]-true_est_ow_finite[,2])^2))/mean(abs(true_est_finite[,2]))
# 
# 
# 
# sqrt(mean((ipw_est[id,3]-true_est_finite[id,3])^2))/mean(abs(true_est_finite[id,3]))
# sqrt(mean((ipw_est_by_group[,3]-true_est_finite[,3])^2))/mean(abs(true_est_finite[,3]))
# sqrt(mean((ipw_est_mao[,3]-true_est_finite[,3])^2))/mean(abs(true_est_finite[,3]))
# 
# 
# sqrt(mean((ow_est[,3]-true_est_ow_finite[,3])^2))/mean(abs(true_est_ow_finite[,3]))
# sqrt(mean((ow_est_by_group[,3]-true_est_ow_finite[,3])^2))/mean(abs(true_est_ow_finite[,3]))
# sqrt(mean((ow_est_mao[,3]-true_est_ow_finite[,3])^2))/mean(abs(true_est_finite[,3]))


# coverage_rate_calc(ow_est[,1],ow_se[,1],true_est_ow_finite[,1])
# coverage_rate_calc(ow_est[,2],ow_se[,2],true_est_ow_finite[,2])
# coverage_rate_calc(ow_est[,3],ow_se[,3],true_est_ow_finite[,3])
# 
# coverage_rate_calc(ipw_est[,1],ipw_se[,1],true_est_finite[,1])
# coverage_rate_calc(ipw_est[,2],ipw_se[,2],true_est_finite[,2])
# coverage_rate_calc(ipw_est[,3],ipw_se[,3],true_est_finite[,3])
