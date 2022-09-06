# Training autoregressive timeseries models.

library(foreach)
library(doParallel)

#' Train autoregressive timeseries mixture of experts (MOE) models.
#'
#'
#' The train_autoregressive_moe function implements the training of
#' autoregressive timeseries models for the data generating processes M4 and M5 in section 5.
#' M4 - the autoregressive timeseries model with static parameters and
#' M5 - autoregressive timeseries model with dynamic parameter.
#'
#' @details The function iterates over N datasets and a vector of discount factor values.
#' For each data set an optimal discount factor is selected based on the log predictive score (lps)
#'
#' @param sim_data: A data frame of the data generated from either M4 or M5.
#' @param discount_grid: A grid on (0, 1) for the discount factor parameter.
#' @param batch_size: The size of the data batches which partition the data.
#' @param n_particles: Number of particles sampled from the posterior distribution.
#'
#' @return list of:
#' - optimal discounts: a vector the optimal discount values ,
#' - lps: log predictive score of the selected models
#' - lps_data: A vector of the lps values per batch in the selected model.

train_autoregressive_models <- function(sim_data, discount_grid, batch_size, n_particles = 1000){
  datasets <- seq_len(ncol(sim_data))
  opts <- list(chunkSize=length(discount_grid))
  # function for aggregating results in the foreach iterations.
  compare_model <- function(model1, model2) if (which.max(c(model1$lps, model2$lps)) == 1) model1 else model2
  registerDoParallel(detectCores())

  simul_output<-foreach(d=datasets,.combine='rbind',
                  .options.nws=opts, .errorhandling="pass")%:%
    foreach(discount = discount_grid,.combine='compare_model',
            .inorder=FALSE, .errorhandling="pass") %dopar%{
              library(dmoe)
              preprocess <- function(series){
                x <- log(1 + series[-length(series)])
                y <- series[-1]
                t <- seq_len(length(y))
                return(data.frame(time=t, y=y, x=x))
              }
              n_particles = n_particles
              sim_data = sim_data
              batch_size = batch_size

              batch_intervals <- seq(1, nrow(sim_data) + 1, by = batch_size)
              mix_col <- c(3)
              exp_col <- c(3)
              n_comp <- 2
      # train models with various discount factors
              model <-  tryCatch({
                dmoe(Data=preprocess(sim_data[,d]), intervals=batch_intervals,
                     particle_size=n_particles, mix_col=mix_col,
                     exp_col=exp_col, n_comp=n_comp, alpha=discount)
                },
                error = function(e) {
                    list(lps=-Inf, lps_list=NA)
                  })
              list(data_id = d,
                   discount = discount,
                   lps = model$lps,
                   lps_list = model$lps_list)
              }

  stopImplicitCluster()

  return(list(optimal_discount=unlist(simul_output[,2]),
              lps=unlist(simul_output[,3]),
            lps_data=as.data.frame(simul_output[,4],
                   col.names = as.character(seq_len(length(simul_output[,4]))))
            ))
}
