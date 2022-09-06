# Implementation of the simulation study for the data generating processes M1-M3 in section 5.

library(LaplacesDemon)
library(foreach)
library(doParallel)

#' Simulate data from the M1-M3 data generating processes.
#' M1 - static Poisson data generating process.
#' M2 - dynamic Poisson data generating process.
#' M3 - dynamic mixture of Poisson data generating process.
#' @param paramdim Dimension of the regression coefficients (including the intercept)
#' @param dgp Data generating process.
#' @param dates A vector of time points in the date format "yyy-mm-dd".
#' @param n  batch sample size.
#'
#' @details Only three data generating processes are implemented. The `dgp` can
#' take values `dgp = "staticpois"` to simulate data from a Poisson dgp with static parameter,
#' `dgp = "dynamicpois"` to simulate data from a Poisson dgp with dynamic parameter, and
#' `dgp = "dynamicpoismix"` to simulate data from a mixture of Poisson components dgp with dynamic parameter.
#'
#' @return simulated training and test data.

simulate_data <-function(paramdim, dgp, dates, n=100){

  init_par <- cbind(rep(c(1, log(0.5)),paramdim),
                    rep(c(-2, log(0.5)),paramdim),
                    rep(c(2,-1),paramdim))[1:paramdim,]
  v <- 0.04
  L <- 2*length(dates)
  if(dgp=="staticpois"){
    beta <- rbind(init_par[,1])
    sim.data <- dmoe::static_poisson(beta, paramdim-1, dates, n)
  }
  if(dgp=="dynamicpois"){
    beta<-array(NA,c(L, paramdim))
    for(d in 1: paramdim){
      beta[,d]<-init_par[d,1]+cumsum(rnorm(L,0,v))
    }
    sim.data <- dmoe::dynamic_poisson(beta, paramdim-1, dates, n)
  }
  if(dgp=="dynamicpoismix"){
    beta<-array(NA,c(L, paramdim, 3))
    for(d in 1: paramdim){
      beta[,d,1]<-init_par[d,1]+cumsum(rnorm(L,0,v))
      beta[,d,2]<-init_par[d,2]+cumsum(rnorm(L,0,v))
      beta[,d,3]<-init_par[d,3]+cumsum(rnorm(L,0,v))
    }
    sim.data <- dmoe::dynamic_poisson_mixture(beta[,,1:2], beta[,,3],
                                                paramdim-1, dates, n)
  }
  return(sim.data)
}

#' Training models for different combinations of the discount factor and number
#'  of components.
#'
#'  @param sim_data Data frame containing the simulated data.
#'  @param discount_gid A grid on the discount factor.
#'  @param components A vector of different number of components.
#'  @param time_intervals Ordered time points partitioning the time.
#'  @param n_particles Sample size of the particles sample.
#'
#'  @returns Matrix of log predictive scores for all trained models.

train_models <- function(sim_data, discount_grid, components,
                         time_intervals, n_particles = 1000){

  opts <- list(chunkSize=length(discount_grid))
  registerDoParallel(detectCores())

  lps_mat<-foreach(n_comp=components, .combine='cbind',
                   .options.nws=opts, .errorhandling="pass")%:%
    foreach(discount = discount_grid, .combine='c',
            .inorder=FALSE, .errorhandling="pass") %dopar%{

              library(dmoe)
              mix_col <- c(3)
              exp_col <- c(3)
              sim_data <- sim_data
              discount_grid <- discount_grid
              time_intervals <- time_intervals
              n_particles <- n_particles

              # train models with various discount factors
              tryCatch({
                model <- dmoe::dmoe(Data=sim_data, intervals=time_intervals,
                                    particle_size=n_particles, exp_col=exp_col,
                                    mix_col=mix_col, n_comp=n_comp, alpha=discount)
                lps = model$lps

              },
              error = function(e) {
                lps=-Inf
              })

              lps - 2*n_comp

            }

  stopImplicitCluster()
  colnames(lps_mat) <- components
  rownames(lps_mat) <- discount_grid
  return(lps_mat)

}

#' Select the optimal discount and number of components based on the
#' log predictive score (lps).
#' @param LPS.mat Matrix containing lps values.
#'
#' @returns indices of the optimal discount and number of components

select_optimal_hyperpars<- function(LPS.mat){
  # return the index of the best model in terms of the log predictive scores
  selected_alpha_ncomp <- as.vector(which(LPS.mat==max(LPS.mat), arr.ind = T))
  return(selected_alpha_ncomp)
}

#' A wrapper of the simulation study of selecting models using
#' log predictive score (lps). The simulation study represents M1-M3
#' data generating processes in section 5 of the paper.
#'
#' @param dgp Data generating process.
#' @param iterations Number of iterations of the simulation study.
#' @param n_particles Sample size of the particles sample.
#' @returns Frequency of selecting the models under consideration.
#' @returns The difference of the lps for optimal models and their corresponding static models.

#' @details Only three data generating processes are implemented. The `dgp` can
#' take values `dgp = "staticpois"` to simulate data from a Poisson dgp with static parameter,
#' `dgp = "dynamicpois"` to simulate data from a Poisson dgp with dynamic parameter, and
#' `dgp = "dynamicpoismix"` to simulate data from a mixture of Poisson components dgp with dynamic parameter.
#'

simulation_study_wrapper <- function(dgp, iterations, n_particles = 1000){
  # dates used to simulate data
  dates<-seq(as.Date("2017-01-01"), as.Date("2017-12-31"), length.out=10)
  time_intervals <- c(dates[1]-1, dates) # partition of the time.
  discounts <- c(0.3,.4,.5,.6,.7,.8,.9,.95,.99) # discount factors
  components <- c(1,2,3) # number of components in the models.
  paramdim <- 2 # dimension of the regression coefficients
  # object to save the selection frequencies of the trained models
  ModelFrequencyMat<-matrix(0,nrow=length(discounts),ncol=length(components))
  row.names(ModelFrequencyMat) <- discounts
  colnames(ModelFrequencyMat) <- components
  optim_static_diff<-c()
  for(iter in seq_len(iterations)){

    message(paste0("\n", "Starting iteration  ", iter, " of ", iterations, "\n"))
    sim.data <- simulate_data(paramdim=paramdim, dgp=dgp, dates=dates) # simulate data

    models_lps <- train_models(sim_data=sim.data, discount_grid=discounts,
                               components=components,time_intervals=time_intervals,
                               n_particles= n_particles)

    message(paste0("\n", "LPS matrix: ", "\n"))
    print(models_lps)

    selected_model <- select_optimal_hyperpars(models_lps)
    optim_lps_diff <- models_lps[selected_model[1], selected_model[2]] -
      models_lps[length(discounts), selected_model[2]]

    message(paste0("\n", "Optimum LPS diff: ", "\n"))
    print(optim_lps_diff)

    optim_static_diff <- c(optim_static_diff, optim_lps_diff )

    # update the model selection frequencies.
    ModelFrequencyMat[selected_model[1], selected_model[2]] <- ModelFrequencyMat[
      selected_model[1], selected_model[2]] + 1

    message(paste0("\n", "Frequency matrix: ", "\n"))
    print(ModelFrequencyMat)

  }
  return(list(freq=ModelFrequencyMat, LPS.dif=optim_static_diff))
}
