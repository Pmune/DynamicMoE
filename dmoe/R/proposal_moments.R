#' Proposal function's moments for Single component models
#'
#' Computes the updated first and second moments for the proposal density of single component models.
#' @param y A vector of response variable values.
#' @param DesingMatrix  A matrix of dimension h x (p+1) containing the covariates data, where h is the number of observations (equal to the length of y) and p represent the dimension the covariates.
#' Note that the the first column is a column  of 1s which is added in order to allow for the intercept in the model.
#' @param effect.param The prior mean of the regression parameter.
#' @param cov_param The prior covariance matrix of the regression parameter.
#'
#' @return Updated mean and covariance matrix of the regression parameter.

single_component_proposal <- function(y, DesignMatrix, effect.param, cov_param){
  for(h in seq_len(length(y))){
    cov.link <- cov_param %*% DesignMatrix[h,] # covariance value from the link function
    P<-as.vector(DesignMatrix[h,] %*% cov.link) # prior variance
    m<-as.numeric(log(1+P*y[h])-log(1 + P * exp(DesignMatrix[h,] %*% effect.param)))
    effect.param <- effect.param + (cov.link/P) * m
    cov_param <- cov_param - ((cov.link %*% t(cov.link))) * as.numeric(y[h]/(1 + P * y[h]))

  }
  return(list(Means=effect.param,cov=cov_param))
}

#' Proposal function's moments for multiple component models.
#'
#' Computes the updated first and second moments for the proposal density of multiple component models.
#' @param y A vector of response variable values.
#' @param x  A matrix of dimension h x (p+1) containing the covariates data, where h is the number of observations (equal to the length of y) and p represent the dimension the covariates in the component models.
#' Note that the the first column is a column  of 1s which is added in order to allow for the intercept in the model.
#' @param z  A matrix of dimension h x (q+1) containing the covariates data, where h is the number of observations (equal to the length of y) and q represent the dimension the covariates in the mixture weights model.
#' Note that the the first column is a column  of 1s which is added in order to allow for the intercept in the model.
#' @param theta The prior mean of the regression parameter in both component and mixture weight models.
#' @param cov_theta The prior covariance matrix of the regression parameter in both component and mixture weight models.
#' @param n_comp Number of mixture components.
#'
#' @return Updated mean and covariance matrix of the regression parameter.

multiple_component_proposal <-function(y, x, z, theta, cov_theta, n_comp){
  for(h in seq_len(length(y))){
    prior <- prior_moments(x[h,], z[h,], theta, cov_theta, n_comp)
    # compute the inverse of the prior variance of linear predictors
    prior_var_linear_pred_inv <- robust_inverse(prior$var_linear_pred)
    # compute the posterior moments of the linear predictor
    posterior_linear_pred <- posterior_linear_moments(response = y[h],
                                               linear_pred = prior$linear_pred,
                                               prior_var_inv = prior_var_linear_pred_inv,
                                               n_comp = n_comp)
    # update the proposal moments
    theta <- theta + c(prior$covar_linear_pred %*% prior_var_linear_pred_inv %*%
                         t(posterior_linear_pred$location - prior$linear_pred))
    cov_theta <- cov_theta - ((prior$covar_linear_pred %*% prior_var_linear_pred_inv) %*%
                                (diag(length(prior$linear_pred)) -
                                   posterior_linear_pred$variance %*% prior_var_linear_pred_inv) %*%
                                t(prior$covar_linear_pred))

    cov_theta <- 0.5*(cov_theta + t(cov_theta)) + diag(1E-5, nrow(cov_theta)) # ensure that the covariance is symetric
  }

  return(list(mean_value = theta, cov_value = cov_theta ))
}
