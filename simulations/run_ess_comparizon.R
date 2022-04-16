dates<-seq(as.Date("2017-01-01"),as.Date("2017-12-31"),length.out=12)

n = 100
dim.design = 1
L <- length(dates) # length of the time intervals
beta1<-cbind(1.11+cumsum(rnorm(2*L,0,0.08)), 2.17+cumsum(rnorm(2*L,0,0.15)))
beta2<-cbind(-0.08+cumsum(rnorm(2*L,0,0.07)),1.94+cumsum(rnorm(2*L,0,0.1)))
beta<-array(NA,c(nrow(beta1),ncol(beta1),2))
beta[,,1] <- beta1
beta[,,2] <- beta2
theta<-cbind(2.63+cumsum(rnorm(2*L,0,0.08)), -4.41+cumsum(rnorm(2*L,0,.17)))


time_intervals<-c(dates[1]-1,dates)
mix_col<-c(3) # column used to define clusters
exp_col<-c(3) # columns used to estimate the mean #TR

n_particles <- 1000 # number of particles (recommended number of particles: 10000)
n_comp <- 2
v_init <-1
alpha <- 0.5


ess_comparizon <- ess_comparizon_study(sim.data,time_intervals, n_particles,
                     mix_col, exp_col, v_init,
                     n_comp, alpha, iterations=10)
