rm(list=ls())
graphics.off()
cat("\014")

setwd("C:/Users/parfait/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/codes/training_models/")
setwd("C:/Users/eparfmu/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/codes/training_models/")

data_file_path <- 'C:/Users/parfait/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/codes/training_models/'
data_file_path <- 'C:/Users/eparfmu/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/codes/training_models/'

# load the data 
source("/data_cleaner.R")
source("interval_predictions.R")

# load posterior output
load("selected_model_posterior.R")
load("posterior_for_model_efficiency_new.R")
# setting plot input
exp.col<- c("ChangedModule","FileComplexity")
mix.col<- c("Commits")
n_comp <- 2
pred_intervals <- c(1, 9, 20)
y.max <- 10
L<-length(time.intervals) - 1 # length of the time intervals

# prediction based on the reference posterior
posterior_ref = list(ref = model_posterior$poisson_mix)
predictions_ref <- interval_predictions(Data, Test.data, time.intervals, posterior_ref,
                                    pred_intervals, exp.col, mix.col, n_comp, y.max)

observed_pred <- predictions$observed
pred_ref <- predictions_ref$predicted

# prediction based on small particle samples

predictions_small_pf <- interval_predictions(Data, Test.data, time.intervals, output_efficiency,
                                        pred_intervals, exp.col, mix.col, n_comp, y.max)

pred_small_pf <- predictions_small_pf$predicted
## plot the predictive distributions for different intervals
windows(12,4)
par(mfrow=c(1,length(pred_intervals)), font=2,family="sans",
    font.axis=1, cex.lab=.8, cex.axis=.8, cex.main=1)
i <- 1

for (j in pred_intervals){
  pred_quantiles <- t(apply(pred_small_pf[[j]],1,
                            function(x)quantile(x,prob=c(0.05, 0.95), na.rm = T)))
  
  
  plot(c(0:y.max),pred_ref[[j]]$ref, type = "n", col = "dark red",
       lty = 1,lwd = 2, ylim = range(pred_quantiles), ylab = " ", xlab =" ", tck = -0.01, 
       xaxp = c(0,y.max,y.max))
  
  polygon(c(rev(0:y.max),0:y.max),c(rev(pred_quantiles[,1]),pred_quantiles[,2]), 
          col ="light blue",border = NA)
  
  lines(c(0:y.max),pred_ref[[j]]$ref,type="l",col="dark blue",lty=1,lwd=2,
        ylab=" ", xlab=" ", tck=-0.01, xaxp=c(0,y.max,y.max))
 
  mtext(paste0("Interval: ", i+1), 3, cex = 0.8, line = 1, font = 2)
  mtext("Number of faults", 2, cex= 0.8, line = 1, font = 2)
  mtext("Density", 2, cex= 0.8, line = 1, font = 2)
  i <- i + 1
}
