rm(list=ls())
graphics.off()
cat("\014")

source("ess_comparison.R")
n_particles <- 2000 # number of particles (recommended number of particles: 10000)
iterations <- 25
alpha <- 0.5
ess_comp_1vars <- ess_comparizon(paramdim=2,n_particles=n_particles,
                                 iterations=iterations, alpha=alpha)
save(ess_comp_1vars,file="results/ess_comp_1vars.R")

load("results/ess_comp_1vars.R")
ess_lb <- ess_comp_1vars$lb
ess_ll <- ess_comp_1vars$ll

J = length(ess_lb)
plot(c(1:J),ess_lb, type = "n",
     lty = 1,lwd = 2, ylim = range(ess_lb, ess_ll),
     ylab = " ", xlab =" ", tck = -0.01, xaxp = c(0,J,J))

lines(c(1:J),ess_ll,type="l",col="red",lty=1,lwd=1)
lines(c(1:J),ess_lb,type="l",col="blue",lty=2,lwd=1)
mtext("Time", 1, cex= 1, line = 3, font = 1)
mtext("ESS per sec", 2, cex= 1, line = 3, font = 1)
legend("topleft",c("LinearBayes", "LocalLinear"),
       col = c( "blue", "red"), lty=c(1,2))
