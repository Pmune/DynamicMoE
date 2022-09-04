#' density function of the generalized Poisson
#'
#' dgenpois computes the density of a given count value (y).


#' @param y : the response variable.
#' @param mu : is the mean parameter for the generalized poisson distribution
#' @param lambda : is the scale parameter describing the dispersion around the mean.
#' @param log_scale: is a logic indicating whether the density value returned should be in log scale or not


#' @returns density value

dgenpois <- function(y, mu, lambda, log_scale = FALSE) {
  logpdf <- y * (log(mu) - log(1 + lambda * mu)) + (y - 1) * log(1 + lambda * y) -
    log(factorial(y)) - mu * (1 + lambda * y) / (1 + lambda * mu) # density function of the generalized poisson distribution
  logpdf[is.na(logpdf)] <- -Inf # avoid na values
  if (log_scale == TRUE) {
    return(logpdf)
  }
  else {
    return(exp(logpdf))
  }
}

#' generalized Poisson likelihood function

#' dgp_likelihood computes the likelihood value of a response variable y given the covariate x.


#' @param y: a vector of the response values
#' @param x: a matrix of the covariate values. the matrix should have the rows = length(y) and columns= number of covariates+1
#' @param beta: regression coefficients on the mean function of the generalized poisson model
#' @param theta: the regression coefficients on the scale (overdispersion) parameter of the generalized poisson model

#' @returns  log_like: vector of the loglikelihood values for y

gp_likelihood <- function(y, x, beta, theta) {
  log_like <- 0
  for (obs in seq_len(length(y))) {
    mu <- exp(x[obs, ] %*% beta) # mean function of the generalized poisson
    lambda <- exp(x[obs, ] %*% theta) # overdispersion function of the generalized poisson
    log_like <- log_like + y[obs] * (log(mu) - log(1 + lambda * mu)) + (y[obs] - 1) * log(1 + lambda * y[obs]) -
      log(factorial(y[obs])) - mu * (1 + lambda * y[obs]) / (1 + lambda * mu) # log likelihood function
  }
  return(log_like) # returns the log likelihood value
}

#' Compute LPS for generized Poisson model
#'
#' the function lps_gpoisson computes the log predictive score: a measure used to compare models


#' @param y : a vector of the response values.
#' @param x: a matrix of the covariate values. the matrix should have the rows = length(y) and columns= number of covariates+1
#' @param particles_beta:  a posterior sampled of the regression cofficients on the mean function
#' @param particles_theta: a posterior sample of the regression coefficients on the overdispersion function
#' @param importance_weight: a vector of the weights of sampled parameters.

#' @returns  lps value

lps_gpoisson <- function(x, y, particles_beta, particles_theta, importance_weight) {
  m <- dim(particles_beta)[1]
  log_like <- as.numeric(apply(as.matrix(1:m), 1,
                               function(ind) gp_likelihood(y, x,
                                                           particles_beta[ind, ],
                                                           particles_theta[ind, ])))
  lps <- log(sum(exp(log_like) * importance_weight))
  return(lps)
}

#' generalized Poisson predictive density function
#' the following function computes the predictive distribution for the generalized Poisson model

# input:
#' @param test_data: test data
#' @param beta: regression parameter on the mean function
#' @param theta: regression parameter on the overdispersion parameter.
#' @param importance_weight: weight of the parameter.
#' @param y_max: maximum value of the response variable.
#' @param training_col: a vector of the names of variables used in the training.


#' @returns empirical predictive distribution
#' @export
gp_pred_dens <- function(test_data, beta, theta, importance_weight, y_max, training_col) {
  d_ftd <- c()
  for (m in seq_len(nrow(test_data))) {
    x <- cbind(1, test_data[m, training_col])
    lambda <- exp(theta %*% t(x))
    mu <- exp(beta %*% t(x))
    # compute the predictive density
    pmf <- apply(as.matrix(c(0:y_max)), 1, function(i) sum(dgenpois(i, mu, lambda) * importance_weight))
    pmf <- pmf / sum(pmf) # normalize the density
    d_ftd <- rbind(d_ftd, pmf)
  }
  d_ftd <- colMeans(d_ftd)
  return(d_ftd) # return the predictive distribution
}
