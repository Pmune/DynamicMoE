
#' Proposal for Single component models using linear Bayes
#'
#' Computes the updated first and second moments for the proposal density of single component models using linear Bayes approximations.
#' @param y A vector of response variable values.
#' @param x  A matrix of dimension h x (p+1) containing the covariates data, where h is the number of observations (equal to the length of y) and p represent the dimension the covariates.
#' Note that the the first column is a column  of 1s which is added in order to allow for the intercept in the model.
#' @param param The prior mean of the regression parameter.
#' @param cov_param The prior covariance matrix of the regression parameter.
#'
#' @return Updated mean and covariance matrix of the regression parameter.

linear_bayes_proposal_1comp <- function(y, x, param, cov_param){
  for(h in seq_len(length(y))){
    cov_link <- cov_param %*% x[h,] # covariance value from the link function
    P<-as.vector(x[h,] %*% cov_link) # prior variance
    m<-as.numeric(log(1 + P * y[h]) - log(1 + P * exp(x[h,] %*% param)))
    param <- param + (cov_link/P) * m
    cov_param <- cov_param - ((cov_link %*% t(cov_link))) * as.numeric(y[h]/(1 + P * y[h]))

  }
  return(list(mean = param, cov = cov_param))
}

#' Proposal for Single component models using local linearization of the target online posterior density.
#'
#' Computes the updated first and second moments for the proposal density for
#'  single component models using local linearization of the target online posterior density.
#' @param y A vector of response variable values.
#' @param x  A matrix of dimension h x (p+1) containing the covariates data, where h is the number of observations (equal to the length of y) and p represent the dimension the covariates in the component models.
#' Note that the the first column is a column  of 1s which is added in order to allow for the intercept in the model.
#' @param param The prior mean of the regression parameter in both component and mixture weight models.
#' @param cov_param The prior covariance matrix of the regression parameter in both component and mixture weight models.
#' @param n_comp Number of mixture components.
#'
#' @return Updated mean and covariance matrix of the regression parameter.

local_linear_proposal_1comp <- function(y, x, param, cov_param){

  for(h in seq_len(length(y))){

    linear_pred <- x[h,]%*%param
    gradient <- gradient_reg_coeff(x[h,], y[h] - exp(linear_pred))
    hessian_matrix <- hessian_reg_coeff(x[h,], exp(-linear_pred))
    # compute the posterior moments of the linear predictor
    updated_moments <- update_linear_moments(param, robust_inverse(cov_param),
                                             gradient, hessian_matrix)
    param <- updated_moments$location
    cov_param <- (0.5 * (updated_moments$variance +
                           t(updated_moments$variance)) +
                    diag(1E-5, nrow(cov_param)))

  }
  return (list(mean = param, cov = cov_param ))
}



#' Proposal function based on the linear Bayes approximation of the target
#' online posterior density.
#'
#' Computes the updated first and second moments for the proposal density for
#'  mixture models using linear Bayes approximation of the target online posterior density.
#' @param y A vector of response variable values.
#' @param x  A matrix of dimension h x (p+1) containing the covariates data, where h is the number of observations (equal to the length of y) and p represent the dimension the covariates in the component models.
#' Note that the the first column is a column  of 1s which is added in order to allow for the intercept in the model.
#' @param z  A matrix of dimension h x (q+1) containing the covariates data, where h is the number of observations (equal to the length of y) and q represent the dimension the covariates in the mixture weights model.
#' Note that the the first column is a column  of 1s which is added in order to allow for the intercept in the model.
#' @param param The prior mean of the regression parameter in both component and mixture weight models.
#' @param cov_param The prior covariance matrix of the regression parameter in both component and mixture weight models.
#' @param n_comp Number of mixture components.
#'
#' @return Updated mean and covariance matrix of the regression parameter.

linear_bayes_proposal_multicomp <-function( y, x, z, param, cov_param, n_comp){
    expected_allocation <- c()
    for(h in seq_len(length(y))){
      prior <- linear_pred_prior_moments(x[h,], z[h,], param, cov_param, n_comp)
      # compute the inverse of the prior variance of linear predictors
      inv_var_linear_pred <- robust_inverse(prior$var_linear_pred)
      expected_allocation <- colSums(rbind(expected_allocation,
                                   compute_allocation_prob(y[h], prior$linear_pred, n_comp)))/h
      expected_allocation <- expected_allocation/sum(expected_allocation)

      # update the component linear predictors
      gradient<- gradient_linear_pred(expected_allocation, y[h],
                                      prior$linear_pred, n_comp)
      hessian_matrix <- hessian_linear_pred(expected_allocation, y[h],
                                            prior$linear_pred, n_comp)
      # compute the posterior moments of the linear predictor
      posterior_linear_pred <- update_linear_moments(prior$linear_pred,
                                                     inv_var_linear_pred,
                                                     gradient,
                                                     hessian_matrix)
      # update the proposal moments
      param <- param + c(prior$covar_linear_pred %*% inv_var_linear_pred %*%
                           t(posterior_linear_pred$location - prior$linear_pred))
      cov_param <- cov_param - ((prior$covar_linear_pred %*% inv_var_linear_pred) %*%
                                  (diag(length(prior$linear_pred)) -
                                     posterior_linear_pred$variance %*% inv_var_linear_pred) %*%
                                  t(prior$covar_linear_pred))

      cov_param <- 0.5*(cov_param + t(cov_param)) + diag(1E-5, nrow(cov_param)) # ensure that the covariance is symetric
    }

  return(list(mean = param, cov = cov_param ))
}

#' Proposal function based on the local linearization of the target online posterior density.
#'
#' Computes the updated first and second moments for the proposal density for
#'  mixture models using local linearization of the target online posterior density.
#' @param y A vector of response variable values.
#' @param x  A matrix of dimension h x (p+1) containing the covariates data, where h is the number of observations (equal to the length of y) and p represent the dimension the covariates in the component models.
#' Note that the the first column is a column  of 1s which is added in order to allow for the intercept in the model.
#' @param z  A matrix of dimension h x (q+1) containing the covariates data, where h is the number of observations (equal to the length of y) and q represent the dimension the covariates in the mixture weights model.
#' Note that the the first column is a column  of 1s which is added in order to allow for the intercept in the model.
#' @param param The prior mean of the regression parameter in both component and mixture weight models.
#' @param cov_param The prior covariance matrix of the regression parameter in both component and mixture weight models.
#' @param n_comp Number of mixture components.
#'
#' @return Updated mean and covariance matrix of the regression parameter.

local_linear_proposal_multicomp <- function(y, x, z, param, cov_param, n_comp){

    for(h in seq_len(length(y))){
      prior <- linear_pred_prior_moments(x[h,], z[h,], param, cov_param, n_comp)
      # compute the inverse of the prior variance of linear predictors
      expected_allocation <- compute_allocation_prob(y[h],
                                                     prior$linear_pred, n_comp)
      # update the component linear predictors

      gradient <- gradient_reg_coeff(design_matrix(x[h,], z[h,], n_comp),
                                     t(gradient_linear_pred(expected_allocation,
                                                          y[h],
                                                          prior$linear_pred,
                                                          n_comp)))

      hessian_matrix <- hessian_reg_coeff(design_matrix(x[h,], z[h,], n_comp),
                                          hessian_linear_pred(expected_allocation,
                                                              y[h],
                                                              prior$linear_pred,
                                                              n_comp))
      # compute the posterior moments of the linear predictor
      updated_moments <- update_linear_moments(param, robust_inverse(cov_param),
                                               gradient, hessian_matrix)
      param <- updated_moments$location
      cov_param <- 0.5 * (updated_moments$variance + t(updated_moments$variance)) +
        diag(1E-5, nrow(cov_param))

    }

  return (list(mean= param, cov = cov_param ))
}


