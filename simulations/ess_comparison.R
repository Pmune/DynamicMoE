library(dmoe)
library(LaplacesDemon)
ess_moe <- function(beta, paramdim, dates,n = 100, n_particles=100, v_init=1, alpha=.5){
  time_intervals<-c(dates[1]-1, dates)
  sim.data <- dmoe::dynamic_poisson_mixture(beta[,,1:2], beta[,,3], paramdim-1,dates,n)
  mix_col<-c(3:(2+paramdim-1)) # column used to define clusters
  exp_col<-c(3: (2+paramdim-1)) # columns used to estimate the mean #TR
  ess <- c()
  start_time <- Sys.time()
  output <- dmoe::dmoe(sim.data, time_intervals, n_particles, mix_col, exp_col,
                 n_comp=2, v_init=v_init, alpha=alpha,
                 proposal_method="linearbayes", return_all=F)
  end_time <- Sys.time()
  message(paste0("\n", "computation time for linear Bayes proposal:",
                 end_time-start_time," ", attr(end_time-start_time, "units")))
  message(paste0("LPS: ", output$LPS))
  ess <- cbind(ess, output$ess)
  start_time <- Sys.time()
  output <- dmoe::dmoe(sim.data, time_intervals, n_particles, mix_col, exp_col,
                 n_comp=2, v_init=v_init, alpha=alpha,
                 proposal_method ="locallinear", return_all=F)
  end_time <- Sys.time()
  message(paste0("\n", "computation time for local linear proposal:",
                 end_time-start_time," ", attr(end_time-start_time, "units")))
  ess <- cbind(ess, output$ess)
  message(paste0("LPS: ", output$LPS))
  return (list(lb=ess[,1], ll =ess[,2]))
}

ess_comparizon <- function(paramdim, n_particles, iterations=10, alpha=.5){
    set.seed(600)
    dates<-seq(as.Date("2017-01-01"),as.Date("2017-12-31"),length.out=12)
    time_intervals<-c(dates[1]-1,dates)
    n = 100
    L <- length(dates) # length of the time intervals
    init_par <- cbind(rep(c(0.11, 2.29),paramdim),
                      rep(c(-0.08, 1.94),paramdim),
                      rep(c(2.63,-4.41),paramdim))[1:paramdim,]
    v <- cbind(rep(c(0.08,0.15),paramdim),
               rep(c(0.08,0.1),paramdim),
               rep(c(0.08,0.17),paramdim))[1:paramdim,]
    beta<-array(NA,c(L, paramdim, 3))
    for(d in 1: paramdim){
      beta[,d,1]<-init_par[d,1]+cumsum(rnorm(L,0,v[d,1]))
      beta[,d,2]<-init_par[d,2]+cumsum(rnorm(L,0,v[d,2]))
      beta[,d,3]<-init_par[d,3]+cumsum(rnorm(L,0,v[d,3]))
    }
    ess_ll <- ess_lb <- c()
    for (iter in 1:iterations){
      # sim.data<-static_poisson(train_parameter,dim.design,dates,n)
      skip <- FALSE
      skept_iter <- 0
      result <- tryCatch(
        {
        ess_moe(beta, paramdim, dates, n_particles=n_particles,  alpha=alpha)},
        error = function(e) {
          skip <- TRUE
          skept_iter <- skept_iter + 1
         cat(paste0("Skipping an iteration! ", skept_iter, "are skept","\n"))
         })
      if(skip){
        next
      }else{
        ess_ll <- rbind(ess_ll, result$ll)
        ess_lb <- rbind(ess_lb, result$lb)
        }

    }
   return(list(lb=colMeans(ess_lb), ll=colMeans(ess_ll)))
}
