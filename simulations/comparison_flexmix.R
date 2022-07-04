rm(list=ls())
graphics.off()
cat("\014")
param_ordering <- function(param, par_dim=2, ncomp=2){
permuted <- array(0,dim(param))
for(i in seq_len(nrow(param))){
  p <- param[i,]
  if (p[1]<=p[3]){
    permuted[i,] <- c(p[c(1:2)], p[c(3:4)], p[-c(1:4)])
  }else{
    permuted[i,] <- c(p[c(3:4)], p[c(1:2)], p[-c(1:4)])
  }
}
return(permuted)
}
# load simulation study file
source('model_selection_simulation_study.R')
library(dmoe)
library(flexmix)
true_par <- cbind(c(-0.1, -0.5), c(0.1, 0.5), c(2,-4))
x_grid <- seq(-1, 1, length.out=100)


compute_mean <- function(beta, x_grid, ncomp=2){
  lin_pred <- cbind(1, x_grid)%*%beta
  mean_val <- apply(lin_pred,1,
        function(eta) sum(exp(eta[1:ncomp])*c(1 , exp(eta[-c(1:ncomp)]))/
                                                sum(1 + exp(eta[-c(1:ncomp)]))))
  return (mean_val)
}

dates<-seq(as.Date("2017-01-01"), as.Date("2017-12-31"), length.out=12)
time_intervals <- c(dates[1]-1, dates) # partition of the time.

n_particles <- 1000


freq_dmoe <- c(0,0,0)
freq_flex <- c(0,0,0)

for(iter in 1:30){
  lps <- bic <-c()
  sim.data <- simulate_data(paramdim=2, dgp="staticpoismix", dates=dates)
  train.data <- sim.data$train
  for(k in 1:3){
  dmoe_model <- dmoe::dmoe(train.data, time_intervals, n_particles, mix_col=c(3),
                           exp_col=c(3), n_comp=k, alpha=.99, v_init=10, static=F)
  lps <- c(lps,dmoe_model$lps)
  flex_model <- flexmix::flexmix(y~x, data = train.data, k=k,
                                 model = FLXMRglm(family = "poisson"),
                                 concomitant = FLXPmultinom(~x))
  bic <- c(bic,BIC(flex_model))
  }
  freq_dmoe[which.max(lps)] <-  freq_dmoe[which.max(lps)] + 1
  freq_flex[which.min(bic)] <-  freq_flex[which.min(bic)] + 1

  print("frequence dmoe:")
  print(freq_dmoe)
  print("frequency flemix:")
  print(freq_flex)
}

# predictive mean function
k=2
start <- Sys.time()
dmoe_model <- dmoe::dmoe(train.data, time_intervals, n_particles, mix_col=c(3),
                         exp_col=c(3), n_comp=k, alpha=.99, v_init=10, static=F)
end <- Sys.time()
print(end-start)
dmoe_model$lps
dmoe_model$ess
flex_model <- flexmix::flexmix(y~x, data = train.data, k=k,
                               model = FLXMRglm(family = "poisson"),
                               concomitant = FLXPmultinom(~x))
dmoe_par <- colSums((dmoe_model$particles)*dmoe_model$weights)
dmoe_par
flex_par <- cbind(parameters(flex_model),
                  parameters(flex_model, which="concomitant")[,-1])

flex_par

flex_mean <- compute_mean(flex_par, x_grid,2)
dmoe_mean <- compute_mean(matrix(dmoe_par, nrow=2), x_grid, k)
true_mean <- compute_mean(true_par, x_grid)



plot(x_grid,true_mean, type = "n",
     lty = 1,lwd = 2, ylim = range(true_mean, dmoe_mean, flex_mean),
     ylab = " ", xlab =" ", tck = -0.01)
lines(x_grid,true_mean,type="l",col="black",lty=3,lwd=1)
lines(x_grid,flex_mean,type="l",col="red",lty=1,lwd=1)
lines(x_grid,dmoe_mean,type="l",col="blue",lty=2,lwd=1)
mtext("x", 1, cex= 1, line = 3, font = 1)
mtext("Predictive mean", 2, cex= 1, line = 3, font = 1)
legend("topright",c("TrueMean", "Flexmix","LinearBayes"),
       col = c("black", "red", "blue"), lty=c(3,1,2))


