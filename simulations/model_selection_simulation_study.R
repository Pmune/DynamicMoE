library(LaplacesDemon)
library(dmoe)

#' Simulate data from a give data generating process.
#' @param paramdim Dimension of the regression coefficients (including the intercept)
#' @param dgp Data generating process.
#' @param dates A list of dates.
#' @param n  sample size
#'
#' @details Only three data generating processes are implemented. The `dgp` can
#' take values `dgp = "staticpois"` to simulate data from a Poisson dgp with static parameter,
#' `dgp = "dynamicpois"` to simulate data from a Poisson dgp with dynamic parameter, and
#' `dgp = "dynamicpoismix"` to simulate data from a mixture of Poisson components dgp with dynamic parameter.
#'
#' @return simulated training and test data.

simulate_data <-function(paramdim, dgp, dates, n=100){
  # init_par <- cbind(rep(c(0.11, 2.29),paramdim),
  #                   rep(c(-0.08, 1.94),paramdim),
  #                   rep(c(.5,-2),paramdim))[1:paramdim,]
  #
  init_par <- cbind(rep(c(-0.1, -0.5),paramdim),
                    rep(c(-0.1, 0.5),paramdim),
                    rep(c(.5,-2),paramdim))[1:paramdim,]
  # v <- cbind(rep(c(0.08,0.15),paramdim),
  #            rep(c(0.08,0.1),paramdim),
  #            rep(c(0.08,0.17),paramdim))[1:paramdim,]
  v <- 0.1
  L <- 2*length(dates)
  if(dgp=="staticpois"){
    beta <- rbind(init_par[,1])
    train.data <- dmoe::static_poisson(beta, paramdim-1, dates, n)
    test.data <- dmoe::static_poisson(beta, paramdim-1, dates, n)
  }
  if(dgp=="dynamicpois"){
    beta<-array(NA,c(L, paramdim))
    for(d in 1: paramdim){
      beta[,d]<-init_par[d,1]+cumsum(rnorm(L,0,v))
    }
    train.data <- dmoe::dynamic_poisson(beta, paramdim-1, dates, n)
    test.data <- dmoe::dynamic_poisson(beta, paramdim-1, dates, n)
  }
  if(dgp=="dynamicpoismix"){
    beta<-array(NA,c(L, paramdim, 3))
    for(d in 1: paramdim){
      beta[,d,1]<-init_par[d,1]+cumsum(rnorm(L,0,v))
      beta[,d,2]<-init_par[d,2]+cumsum(rnorm(L,0,v))
      beta[,d,3]<-init_par[d,3]+cumsum(rnorm(L,0,v))
    }
    train.data <- dmoe::dynamic_poisson_mixture(beta[,,1:2], beta[,,3],
                                                paramdim-1, dates, n)
    test.data <- dmoe::dynamic_poisson_mixture(beta[,,1:2], beta[,,3],
                                               paramdim-1, dates, n)
  }

  if(dgp=="staticpoismix"){
    beta <- array(NA,c(L, paramdim, 3))
    for(k in 1:3){
      beta[,,k] <-  matrix(rep(init_par[,k], each=L), ncol=paramdim)
    }
    train.data <- dmoe::dynamic_poisson_mixture(beta[,,1:2], beta[,,3],
                                                paramdim-1, dates, n)
    test.data <- dmoe::dynamic_poisson_mixture(beta[,,1:2], beta[,,3],
                                               paramdim-1, dates, n)
  }
  return(list(train = train.data, test = test.data))
}

#' Training models for different combinations of the discount factor and number
#'  of components.
#'
#'  @param train.data Training dataset.
#'  @param discounts A vector of different discount factors.
#'  @param components A vector of different number of components.
#'  @param time_intervals Ordered time points partitioning the time
#'  @param n_particles Number of particles.
#'
#'  @returns Matrix of log predictive scores for all trained models.

train_models <- function(train.data, discounts, components, time_intervals,
                         n_particles){
  LPS.mat<-matrix(0,nrow=length(discounts),ncol=length(components)) # object fot storing the LPS values

  # fit models with different discount factors and number of components
  for(j in components){ # loop over the components
    for(a in seq_len(length(discounts))){ # loop over the discount factors
      output <- dmoe::dmoe(train.data, time_intervals, n_particles, mix_col=c(3),
                         exp_col=c(3), n_comp=j, alpha=discounts[a])
      LPS.mat[a,j] <- output$lps
    }
  }
  return(LPS.mat)
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

#' Compare the optimal model (based on lps) and the corresponding static model.
#' @param test.data Test dataset.
#'  @param optim_ncomp Optimal number of components.
#'  @param optim_alpha Optimal discount factor.
#'  @param time_intervals Ordered time points partitioning the time
#'  @param n_particles Number of particles.
#'
#'  @returns The difference of the lps for optimal models and their corresponding static models.

compare_static_optimal_models <- function(test.data, optim_ncomp, optim_alpha,
                                          time_intervals, n_particles){
  mix_col <- exp_col <- c(3)
  optimal_model <- dmoe::dmoe(test.data, time_intervals, n_particles, mix_col,
                              exp_col, n_comp=optim_ncomp, alpha=optim_alpha)
  static_model <- dmoe::dmoe(test.data, time_intervals, n_particles, mix_col,
                             exp_col, n_comp=optim_ncomp, alpha=0.99)

   return(optimal_model$lps - static_model$lps)

}

#' A wrapper of the simulation study of selecting models using
#' log predictive score (lps). The simulation study is presented
#'  in section 4.2 of the paper.
#'
#' @param dgp Data generating process.
#' @param iterations Number of iterations of the simulation study.
#'
#'  @returns Frequency of selecting the models under consideration.
#'  @returns The difference of the lps for optimal models and their corresponding static models.

#' @details Only three data generating processes are implemented. The `dgp` can
#' take values `dgp = "staticpois"` to simulate data from a Poisson dgp with static parameter,
#' `dgp = "dynamicpois"` to simulate data from a Poisson dgp with dynamic parameter, and
#' `dgp = "dynamicpoismix"` to simulate data from a mixture of Poisson components dgp with dynamic parameter.
#'

simulation_study_wrapper <- function(dgp, iterations){
  set.seed(700)
  # dates used to simulate data
  dates<-seq(as.Date("2017-01-01"), as.Date("2017-12-31"), length.out=12)
  time_intervals <- c(dates[1]-1, dates) # partition of the time.
  discounts <- c(.4,.5,.6,.7,.8,.9,.95,.99) # discount factors (alpha values)
  components <- c(1,2,3) # number of components in the models.
  n_particles <- 2000 # number of particles approximating the posterior density
  paramdim <- 2 # dimension of the regression coefficients
  # object to save the selection frequencies of the trained models
  ModelFrequencyMat<-matrix(0,nrow=length(discounts),ncol=length(components))
  row.names(ModelFrequencyMat) <- discounts
  colnames(ModelFrequencyMat) <- components
  LPS.dif<-c()
  for (iter in 1:iterations){
    skip <- FALSE
    sim.data <- simulate_data(paramdim, dgp, dates) # simulate data
    result <- tryCatch({
      #compute lps for each combination of the discount and number of components
      lps.mat <- train_models(sim.data$train, discounts, components,
                              time_intervals, n_particles)
      # select the best model
      optim_index <- select_optimal_hyperpars(lps.mat)
      # update the model selection frequencies.
      ModelFrequencyMat[optim_index[1], optim_index[2]] <- ModelFrequencyMat[
        optim_index[1], optim_index[2]] + 1

      # optimal discount and number of components.
      optim_ncomp <- components[optim_index[2]]
      optim_alpha <- discounts[optim_index[1]]
      # compare optimal and static models
      lps_opt <- compare_static_optimal_models(sim.data$test, optim_ncomp,
                                                optim_alpha,time_intervals,
                                                n_particles)
      LPS.dif <- c(LPS.dif, lps_opt) # save the comparison results

      cat(paste0("\n", "Current iteration = ", iter, "\n"))
      cat(paste0("\n", "lps diff: ", "\n", LPS.dif))
      cat(paste0("\n", "Frequency matrix: ", "\n"))
      print(ModelFrequencyMat)

        },

      error = function(e) {
        skip <- TRUE
        cat(paste0("\n","Skipping an iteration due to error: ","\n", e))
        })

    if(skip){
      next
    }
  }
 return(list(freq=ModelFrequencyMat,LPS.dif=LPS.dif))
}
