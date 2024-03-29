# Comparison of the effective sample size rate
# generated by the linear Bayes proposal vs the local linear proposal method.

rm(list=ls())

library(dmoe)
library(LaplacesDemon)
library(tidyverse)

#' Compute effective sample size for the linear Bayes and local linear proposal strategies.
#'
#' @param simul.data: data frame containing the simulated data.
#' @param exp_col: Column indices (in the simul.data) for the variable used as
#'  covariates in the component models.
#' @param mix_col: Column indices (in the simul.data) for the variable used as
#'  covariates in the mixture weights models.
#' @param time_intervals A vector of time intervals partitioning the data into batches.
#' @param n_particles: Number of particles sampled from the posterior distribution.
#' @param discount: Discount factor
#'
#' @return list of the effective sample size rates generated by
#' the linear Bayes proposal, and the local linear proposal.
compute_ess <- function(sim.data, exp_col, mix_col, time_intervals, n_particles, discount){

  ess <- c()
  start_time <- Sys.time()
  # train the particle filter with linear Bayes proposal
  output <- dmoe::dmoe(Data = sim.data,
                       intervals = time_intervals,
                       particle_size = n_particles,
                       exp_col = exp_col,
                       mix_col = mix_col,
                       n_comp = 2,
                       alpha = discount,
                       proposal_method="linearbayes")
  end_time <- Sys.time()
  message(paste0("\n", "computation time for linear Bayes proposal:",
                 end_time-start_time," ", attr(end_time-start_time, "units")))
  ess <- cbind(ess, output$ess)
  start_time <- Sys.time()
  # train the particle filter with local linear proposal
  output <- dmoe::dmoe(Data = sim.data,
                       intervals = time_intervals,
                       particle_size = n_particles,
                       exp_col = exp_col,
                       mix_col = mix_col,
                       n_comp = 2,
                       alpha = discount,
                       proposal_method ="locallinear")
  end_time <- Sys.time()
 # message(paste0("\n", "computation time for local linear proposal:",
 #               end_time-start_time," ", attr(end_time-start_time, "units")))
  ess <- cbind(ess, output$ess)
  return (list(lb=ess[,1], ll =ess[,2]))
}

#' Effective sample size comparison simulation study.
#'
#' A wrapper of the function: compute_ess.
#' The wrapper runs several iterations of the function compute_ess and compute
#' the  average effective sample size rate for each of the linear Bayes and local linear proposal method.
#'
#' @param n_particles Number of particles sampled from the posterior distribution.
#' @param iterations Number of independent iterations the function compute_ess is evaluated.
#' @param paramdim Dimension of the regression coefficients (including the intercept) used for simulating data.
#' @param discount: Discount factor.
#'
#' @return list of mean (averaged over the iterations) effective sample size rates generated by
#' the linear Bayes proposal, and the local linear proposal
ess_comparison_simul <- function(iterations, n_particles,  paramdim = 2, discount = .5){
    set.seed(600)
    dates<-seq(as.Date("2017-01-01"), as.Date("2017-12-31"), length.out=10) #simulate dates used in the data partition.
    time_intervals<-c(dates[1]-1,dates) # time intervals partitioning the data into batches.
    L <- length(dates) # length of the time intervals
    init_par <- cbind(rep(c(1, log(.5)), paramdim),
                      rep(c(-2, log(.5)), paramdim),
                      rep(c(2,-1),paramdim))[1:paramdim,]
    v <- 0.04
    simul_param<-array(NA,c(L, paramdim, 3))
    for(d in 1: paramdim){
      simul_param[,d,1]<-init_par[d,1]+cumsum(rnorm(L,0,v))
      simul_param[,d,2]<-init_par[d,2]+cumsum(rnorm(L,0,v))
      simul_param[,d,3]<-init_par[d,3]+cumsum(rnorm(L,0,v))
    }
    n <- 100 #batch sample size.
    mix_col<-c(3:(2+paramdim-1)) # column used to define clusters
    exp_col<-c(3: (2+paramdim-1)) # columns used to estimate the mean
    ess_ll <- ess_lb <- c()
    for (iter in 1:iterations){
      skip <- FALSE
      skept_iter <- 0
      # simulate data from the dynamic mixture of Poisson data generating process (M3 in the section 5).
      sim.data <- dmoe::dynamic_poisson_mixture(
        simul_param[, ,1:2],
        simul_param[, ,3],
        paramdim-1, dates, n)
      result <- tryCatch(
        {
        compute_ess(
          sim.data = sim.data,
          time_intervals = time_intervals,
          exp_col = exp_col,
          mix_col = mix_col,
          n_particles = n_particles,
          discount=discount)},

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

#--------------------------- Run the simulation -------------------------------------

n_particles <- 1000 # number of particles
iterations <- 25 # number of iterations
ess_comp <- ess_comparison_simul(paramdim = 2, n_particles = n_particles, iterations = iterations)

# files are saved with a tag appended to the file name to avoid overwriting files.
# if you have change the tag please make sure to change the file_tag in the plot files.
file_tag <- "1" # file tag

save(ess_comp,file = paste0("results/ess_comp_", file_tag, ".R"))

#load("results/ess_comp.R")
ess_lb <- ess_comp$lb # effective sample size rate linear Bayes method
ess_ll <- ess_comp$ll # effective sample size rate local linear method

# plot the effective sample size rate. Figure 5.4 in the paper.
essplot <- data.frame(time = seq_len(length(ess_lb)),
           linear_bayes=ess_lb, local_linear=ess_ll)%>%
  select(time, linear_bayes, local_linear) %>%
  gather(key = "Method", value = "value", -time)%>%

  ggplot(aes(x = time, y = value)) +
  geom_line(aes(color = Method, linetype = Method),
            size=1, alpha=1, lineend="round") +
  scale_color_manual(values = c("darkred", "steelblue")) +
  theme_bw() +
  theme(
    plot.title = element_text(size=11, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = c(0.1, .9),
    legend.title = element_blank(),
    strip.background = element_rect(fill = "white", linetype = 0)
  ) +
  xlab("Batch index") +
  ylab("ESS per second") +
  scale_x_continuous(breaks=seq(0,10,by=1))

print(essplot)
png("plots/ess_rate_comparison.png", width=600, height=300)
print(essplot)
dev.off()
