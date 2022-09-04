
#' particle filter algorithm for the dynamic generalized Poisson model
#'
#' This function implements the particle filter algorithm for the dynamic generalized Poisson regression model.


#' @param train_data: a data frame containing the training data set.
#' @param intervals: a vector of subsequent time points partitioning the training period.
#' @param n_particles: number of parameters to be sampled (posterior sample size)
#' @param train_col: a vector of column names that are used as covariates. These names should be consistent with the column names in the train_data
#' @param time_ind: column index of the up dates.
#' @param response_ind: column index of the response variable
#' @param alpha: discount factor; it controls the smoothness of the parameters accross intervals.
#' @param r: number of replicates each sampled parameter appears in the posterior sample.
#'    This is used in the resampling step of the particle filter (see file stratified_rejection_sampling.R)
#' @param return_allweights: a logic value that indicates if all weights of the sampled posterior parameters should be returned

#' @return:
#' particles_mean, particles_disp: a sample from the posterior
#' ESS: Effective sample size
#' weights: importance weights
#' lps: log predictive score

#' @export
dgp <- function(train_data, intervals, n_particles, train_col,
                      time_ind = 1, response_ind = 2, alpha = 0.5,
                      r = 2,  return_allweights = F) {
  interval_length <- length(intervals) - 1 # number of time periods
  pb <- txtProgressBar(min = 0, max = interval_length, initial = 0, char = "_", style = 3)
  # Variable declaration and initialization
  design_dim <- length(train_col) + 1 # dimension of the design matrix
  particles_mean <- list() # save sampled particles
  particles_disp <- list()
  weight_state <- rep(1 / n_particles, n_particles) # initialize weights
  covmat_disp <- diag(design_dim)# set the initial covariance
  covmat_mean <- 10*diag(design_dim)
  prev_theta <- rmvnorm(n_particles, rep(0, design_dim), covmat_disp) # sample theta at time point j=0
  mean_theta <- colMeans(prev_theta)
  prev_beta <- rmvnorm(n_particles, rep(0, design_dim), covmat_mean) # sample Beta at time point j=0
  mean_beta <- colMeans(prev_beta)
  allweights <- list()


  # propagate paricles at time point j>1
  ess <- c() # effective sample size
  lps_score <- 0
  for (j in 2:(interval_length + 1)) {

    # generating risk set.
    dj <- train_data[which(train_data[, time_ind] <= intervals[j] & train_data[, time_ind] > intervals[j - 1]), ] # interval data
    yj <- dj[, response_ind]
    # effect covariates
    xj <- as.matrix(cbind(rep(1, nrow(dj)), dj[, c(train_col)])) # design matrix
    # propose particle at time j
    prop_meancov <- gp_linearbayes(yj, xj, mean_beta, mean_theta, covmat_mean, covmat_disp)
    theta <- rmvnorm(n_particles, prop_meancov$disp_effect, prop_meancov$cov_disp)
    beta <- rmvnorm(n_particles, prop_meancov$mean_effect, prop_meancov$cov_mean)
    # compute weights at time j
    new_weight_state <- apply(as.matrix(1:n_particles), 1, function(m) {
      gp_computeweight(
        yj, xj, weight_state, beta[m, ], theta[m, ], prev_beta, prev_theta,
        covmat_mean , covmat_disp, prop_meancov$mean_effect,
        prop_meancov$cov_mean , prop_meancov$disp_effect, prop_meancov$cov_disp
      )
    }) # denominator
    new_weight_state[which(is.na(new_weight_state))] <- log(1e-300)
    new_weight_state[is.infinite(new_weight_state)] <- log(1e-300)

    # normalize weights
    new_weight_state <- normalize(new_weight_state)
    newsample <- reject_strat_resample(new_weight_state, n_particles / r)
    index <- newsample$index # indexes of sampled particles
    weight_state <- newsample$weights
    prev_theta <- particles_disp[[j - 1]] <- theta[index, ]
    prev_beta <- particles_mean[[j - 1]] <- beta[index, ]
    covmat_disp <- cov(theta[index, ]) * (1 - alpha)/ alpha
    covmat_mean <- cov(beta[index, ]) * (1 - alpha)/ alpha
    mean_theta <- colSums(theta[index, ] * weight_state)
    mean_beta <- colSums(beta[index, ] * weight_state)
    allweights[[j - 1]] <- weight_state


    # computation of the effective sample size: ess and the predictive score
    ess[j - 1] <- 1 / sum(weight_state^2)
    if (j > floor(interval_length/2) + 1) {
      lps_score <- lps_score + lps_gpoisson(xj, yj, particles_mean[[j - 2]],
                                    particles_disp[[j - 2]], allweights[[j - 2]])
    }
    setTxtProgressBar(pb, value = j)
  }
  if (return_allweights) {
    return(list(
      particles_mean = particles_mean, particles_disp = particles_disp, ess = ess,
      lps = lps_score, weights = allweights
    ))
  }
  else {
    return(list(
      particles_mean = particles_mean, particles_disp = particles_disp, ess = ess,
      lps = lps_score, weights = weight_state
    ))
  }
}

