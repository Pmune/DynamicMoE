#' function for computing importance weights for the dynaminc mixture of experts model.
#'
#' This function computes the weights of the sampled parameters.
#' The weights are functions of the likelihood and the prior.


#' @param y: a vector of the response values observed in the batch j.
#' @param x: a matrix of the covariates observed in the batch j.
#' @param prev_weights: weights of the parameters in the previous intervals j-1.
#' @param param: A vector of regression coefficeints (one particle sample).
#' @param prev_particles: Posterior sample for the previous interval j-1.
#' @param prior_mean: Prior mean for the regression coefficients.
#' @param prior_var: Prior covariance matrix for regression coefficients.
#' @param proposal_mean: Proposal mean vector for the regression coefficients.
#' @param proposal_cov: Proposal covariance matrix for the regression coefficients.
#' @param n_comp: Number of mixture components.
#' @returns  updated importance weights.
importance_weight_multicomp<-function(y, x, z, param, prev_weights,
                                      prev_particles, prior_mean, prior_var,
                                      proposal_mean, proposal_cov, n_comp){

  dens_proposal <- dmvnorm(param, proposal_mean, proposal_cov, log=TRUE)
  dens_prior<-log(sum(dmvnorm(prev_particles, param, prior_var)*prev_weights))
  log_weight <- mixture_log_likelihood(y, x, z, param, n_comp) + dens_prior -
                                 dens_proposal
  return(log_weight)
}

#' function for computing importance weights for dynamic Poisson model.
#'
#' This function computes the weights of the sampled parameters.
#' The weights are functions of the likelihood and the prior.


#' @param y: a vector of the response values observed in the batch j.
#' @param x: a matrix of the covariates observed in the batch j.
#' @param prev_weights: weights of the parameters in the previous intervals j-1.
#' @param param: A vector of regression coefficeints (one particle sample).
#' @param prev_particles: Posterior sample for the previous interval j-1.
#' @param prior_mean: Prior mean for the regression coefficients.
#' @param prior_var: Prior covariance matrix for regression coefficients.
#' @param proposal_mean: Proposal mean vector for the regression coefficients.
#' @param proposal_cov: Proposal covariance matrix for the regression coefficients.
#' @returns  updated importance weights.

importance_weight_1comp<-function(y, x, param, prev_weights, prev_particles,
                                  prior_mean, prior_var, proposal_mean,
                                  proposal_cov) {

  dens_proposal <- dmvnorm(param, proposal_mean, proposal_cov, log=TRUE)
  dens_prior<-log(sum(dmvnorm(prev_particles, param, prior_var)*prev_weights))
  log_weight <- poisson_log_likelihood(y, x, param) + dens_prior - dens_proposal
  return(log_weight)
  }


#' function for computing importance weights for dynamic generalized Poisson model.
#'
#' This function computes the weights of the sampled parameters.
#' The weights are functions of the likelihood and the prior.


#' @param y: a vector of the response values observed in the batch j.
#' @param x: a matrix of the covariates observed in the batch j.
#' @param prev_weights: weights of the parameters in the previous intervals j-1.
#' @param beta: a vector of proposed regression parameters on the mean function
#' @param theta: a vector of proposed regression parameters on the overdispersion (scale) function
#' @param prev_beta: a vector of regression parameters on the mean function in the previous interval j-1.
#' @param prev_theta: a vector of regression parameters on the overdispersion (scale) function in the previous interval j-1.
#' @param w_beta: transition matrix defining the evolution of beta from prev_beta.
#' @param w_theta: transition matrix defining the evolution of theta from prev_theta.
#' @param proposalmean_mean: a vector of the mean value of the proposal distribution of the regression parameters on the mean function.
#' @param proposalmean_disp: a vector of the mean value of the proposal distribution of the regression parameters on the overdispersion(scale) function.
#' @param proposalcov_mean: covariance matrix of the proposal distribution of the regression parameters on the mean function.
#' @param proposalcov_disp: covariance matrix of the proposal distribution of the regression parameters on the overdispersion(scale) function.

#' @returns  updated importance weights.


gp_computeweight <- function(y, x, prev_weights, beta, theta, prev_beta, prev_theta, w_beta, w_theta,
                             proposalmean_mean, proposalcov_mean, proposalmean_disp, proposalcov_disp) {
  dim_par_beta <- dim(w_beta)[2]
  dim_par_theta <- dim(w_theta)[2]
  dens_proposal <- dmvnorm(beta, proposalmean_mean, proposalcov_mean, log = TRUE) +
    dmvnorm(theta, proposalmean_disp, proposalcov_disp, log = TRUE)
  dens_prior <- as.numeric(dmvnorm(prev_beta, beta, w_beta + diag(1E-7, dim_par_beta), log = TRUE)) +
    as.numeric(dmvnorm(prev_theta, theta, w_theta + diag(1E-7, dim_par_theta), log = TRUE))
  p_prior <- log(sum(exp(dens_prior) * prev_weights))
  log_weight <- gp_likelihood(y, x, beta, theta) + p_prior - dens_proposal
  log_weight[is.na(log_weight)] <- -Inf
  return(log_weight) # return a vector of the importance weights
}
