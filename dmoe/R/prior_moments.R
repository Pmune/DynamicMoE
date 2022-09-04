#' Compute the prior moments of the linear predictor from the degenerate joint
#' prior distribution of the linear predictor and the regression coefficients.
#'
#' @param x  A matrix of dimension h x p containing the covariates data,
#' where h is the number of observations (equal to the length of y) and
#' p represent the dimension the covariates in the component models.
#' Note that the the first column is a column  of 1s which is added in order
#' to allow for the intercept in the model.
#' @param z  A matrix of dimension h x (q+1) containing the covariates data,
#' where h is the number of observations (equal to the length of y) and
#' q represent the dimension the covariates in the mixture weights model.
#' Note that the the first column is a column  of 1s which is added in order
#' to allow for the intercept in the model.
#' @param param The prior mean of the regression coefficients in both
#' the component and mixture weight models.
#' @param cov_mat The prior covariance matrix of the regression coefficients
#' in both the component and mixture weight models.
#' @param n_comp Number of mixture components.
#'
#' @return The prior mean and variance of the linear predictor.
#' Also, the covariance of the linear predictor and the regression coefficients is returned.
linear_pred_prior_moments <-function(x, z, param, cov_mat, n_comp){

  x_mat <- design_matrix(x,z,n_comp) # design matrix
  linear_pred <- param %*% x_mat
  covar_linear_pred <- cov_mat %*% x_mat
  var_linear_pred <- t(x_mat) %*% covar_linear_pred
  var_linear_pred <- 0.5*(var_linear_pred + t(var_linear_pred)) + diag(1E-5, nrow(var_linear_pred))
  return(list(linear_pred=linear_pred, covar_linear_pred=covar_linear_pred, var_linear_pred=var_linear_pred ))
}


#' Compute the prior moments of the regression coefficients.
#'
#' @param coeff An MxP matrix containing sampled regression coefficients.
#'  M is the size of the sample and P is the dimension of the regression coefficients
#' @param weights A vector of length M containing the weights of each sampled regression coefficient.
#' @return The prior mean and variance of the regression coefficients.

reg_coef_prior_moment <- function(coeff, weights){
  mean_coeff <- colSums(coeff * weights)

  if(length(weights)<30000){
    error_vec <- coeff - matrix(rep(mean_coeff, each=nrow(coeff)), ncol= length(mean_coeff))
    var_coeff <- t(error_vec) %*% diag(weights) %*%  error_vec
  }else{ # to avoid running out of memory in the case where weights is a long vector.
    var_coeff <- var(coeff)
  }
  var_coeff <- 0.5*(var_coeff + t(var_coeff)) + diag(1E-5, nrow(var_coeff))

  return(list(mean = mean_coeff, var = var_coeff))
}
