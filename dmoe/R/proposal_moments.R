
#' Compute the updated first and second moments for the proposal density of
#'  single component models using linear Bayes approximations.
#'
#' @param y A vector of response variable values.
#' @param x  A matrix of dimension h x (p+1) containing the covariates data,
#' where h is the number of observations (equal to the length of y) and
#' p represent the dimension the covariates.
#' Note that the the first column is a column  of 1s which is added
#' in order to allow for the intercept in the model.
#' @param reg_coef A vector of the regression coefficients.
#' @param cov_reg_coef The covariance matrix of the regression coefficient.
#'
#' @return Updated mean and covariance matrix of the regression parameter.

linear_bayes_proposal_1comp <- function(y, x, reg_coef, cov_reg_coef){
  for(h in seq_len(length(y))){
    cov_link <- cov_reg_coef %*% x[h,] # covariance value from the link function
    P<-as.vector(x[h,] %*% cov_link) # prior variance
    m<-as.numeric(log(1 + P * y[h]) - log(1 + P * exp(x[h,] %*% reg_coef)))
    reg_coef <- reg_coef + (cov_link/P) * m
    cov_reg_coef <- cov_reg_coef - ((cov_link %*% t(cov_link))) * as.numeric(y[h]/(1 + P * y[h]))

  }
  return(list(mean = reg_coef, cov = cov_reg_coef))
}



#' Compute the updated first and second moments for the proposal density for
#'  single component models using local linearization of the target online posterior density.
#'
#' @param y A vector of response variable values.
#' @param x  A matrix of dimension h x (p+1) containing the covariates data,
#' where h is the number of observations (equal to the length of y) and
#' p represent the dimension the covariates in the component models.
#' Note that the the first column is a column  of 1s which is added
#' in order to allow for the intercept in the model.
#' @param reg_coef A vector of the regression coefficients.
#' @param cov_reg_coef The covariance matrix of the regression coefficients.

#' @return Updated mean and covariance matrix of the regression parameter.

local_linear_proposal_1comp <- function( y, x, reg_coef, cov_reg_coef){

  for(h in seq_len(length(y))){

    linear_pred <- x[h,]%*%reg_coef
    gradient <- gradient_reg_coeff(x[h,], y[h] - exp(linear_pred))
    hessian_matrix <- hessian_reg_coeff(x[h,], exp(-linear_pred))
    # compute the posterior moments of the linear predictor
    updated_moments <- update_linear_moments(reg_coef, robust_inverse(cov_reg_coef),
                                             gradient, hessian_matrix)
    reg_coef <- updated_moments$location
    cov_reg_coef <- (0.5 * (updated_moments$variance +
                           t(updated_moments$variance)) +
                    diag(1E-5, nrow(cov_reg_coef)))

  }
  return (list(mean = reg_coef, cov = cov_reg_coef ))
}



#' Proposal function based on the linear Bayes approximation of the target
#' online posterior density.
#'
#' Computes the updated first and second moments for the proposal density for
#'  mixture models using linear Bayes approximation of the target online posterior density.
#' @param y A vector of response variable values.
#' @param x  A matrix of dimension h x (p+1) containing the covariates data,
#'  where h is the number of observations (equal to the length of y) and
#'  p represent the dimension the covariates in the component models.
#' Note that the the first column is a column  of 1s which is added
#' in order to allow for the intercept in the model.
#' @param z  A matrix of dimension h x (q+1) containing the covariates data,
#' where h is the number of observations (equal to the length of y) and
#' q represent the dimension the covariates in the mixture weights model.
#' Note that the the first column is a column  of 1s which is added
#' in order to allow for the intercept in the model.
#' @param reg_coef A vector of the regression coefficients in both
#' component and mixture weight models.
#' @param cov_reg_coef The prior covariance matrix of the regression coefficients
#' in both component and mixture weight models.
#' @param n_comp Number of mixture components.
#'
#' @return Updated mean and covariance matrix of the regression coefficients.

linear_bayes_proposal_multicomp <-function( y, x, z, reg_coef, cov_reg_coef, n_comp){
    expected_allocation <- c()
    for(h in seq_len(length(y))){
      prior <- linear_pred_prior_moments(x[h,], z[h,], reg_coef, cov_reg_coef, n_comp)
      # compute the inverse of the prior variance of linear predictors
      inv_var_linear_pred <- robust_inverse(prior$var_linear_pred)
      expected_allocation <- compute_allocation_prob(y[h], prior$linear_pred, n_comp)
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
      reg_coef <- reg_coef + c(prior$covar_linear_pred %*% inv_var_linear_pred %*%
                           t(posterior_linear_pred$location - prior$linear_pred))
      cov_reg_coef <- cov_reg_coef - ((prior$covar_linear_pred %*% inv_var_linear_pred) %*%
                                  (diag(length(prior$linear_pred)) -
                                     posterior_linear_pred$variance %*% inv_var_linear_pred) %*%
                                  t(prior$covar_linear_pred))

      cov_reg_coef <- 0.5*(cov_reg_coef + t(cov_reg_coef)) + diag(1E-5, nrow(cov_reg_coef)) # ensure that the covariance is symetric
    }

  return(list(mean = reg_coef, cov = cov_reg_coef ))
}


#' Proposal function based on the local linearization of the target online posterior density.
#'
#' Computes the updated first and second moments for the proposal density for
#'  mixture models using local linearization of the target online posterior density.
#' @param y A vector of response variable values.
#' @param x  A matrix of dimension h x (p+1) containing the covariates data,
#' where h is the number of observations (equal to the length of y) and
#' p represent the dimension the covariates in the component models.
#' Note that the the first column is a column  of 1s which is added
#' in order to allow for the intercept in the model.
#' @param z  A matrix of dimension h x (q+1) containing the covariates data,
#' where h is the number of observations (equal to the length of y) and
#' q represent the dimension the covariates in the mixture weights model.
#' Note that the the first column is a column  of 1s which is added
#' in order to allow for the intercept in the model.
#' @param reg_coef A vector of the regression coefficients in both
#' component and mixture weight models.
#' @param cov_reg_coef The prior covariance matrix of the regression coefficients
#' in both component and mixture weight models.
#' @param n_comp Number of mixture components.
#'
#' @return Updated mean and covariance matrix of the regression coefficients.


local_linear_proposal_multicomp <- function(y, x, z, reg_coef, cov_reg_coef, n_comp){

    for(h in seq_len(length(y))){
      prior <- linear_pred_prior_moments(x[h,], z[h,], reg_coef, cov_reg_coef, n_comp)
      # compute the inverse of the prior variance of linear predictors

      inv_var_linear_pred <- robust_inverse(prior$var_linear_pred)
      expected_allocation <- compute_allocation_prob(y[h], prior$linear_pred, n_comp)
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
      updated_moments <- update_linear_moments(reg_coef, robust_inverse(cov_reg_coef),
                                               gradient, hessian_matrix)
      reg_coef <- updated_moments$location
      cov_reg_coef <- 0.5 * (updated_moments$variance + t(updated_moments$variance)) +
        diag(1E-5, nrow(cov_reg_coef))

    }

  return (list(mean = reg_coef, cov = cov_reg_coef ))
}


