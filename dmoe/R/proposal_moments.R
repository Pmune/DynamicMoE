
#' Compute the updated first and second moments for the proposal density of
#'  single component models using linear Bayes approximations.
#'
#' @param y A vector of response variable values.
#' @param x  A matrix of dimension h x p containing the covariates data,
#' where h is the number of observations (equal to the length of y) and
#' p represent the dimension the covariates.
#' @param reg_coef A vector of the regression coefficients.
#' @param cov_reg_coef The covariance matrix of the regression coefficient.
#'
#' @return Updated mean and covariance matrix of the regression parameter.

linear_bayes_proposal_1comp <- function(y, x, reg_coef, cov_reg_coef){
  x <- cbind(1, x) # add extra column for the intercept to the covariate matrix
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
#' @param x  A matrix of dimension h x p containing the covariates data,
#' where h is the number of observations (equal to the length of y) and
#' p represent the dimension the covariates in the component models.
#' @param reg_coef A vector of the regression coefficients.
#' @param cov_reg_coef The covariance matrix of the regression coefficients.

#' @return Updated mean and covariance matrix of the regression parameter.

local_linear_proposal_1comp <- function( y, x, reg_coef, cov_reg_coef){
  x <- cbind(1, x) # add extra column for the intercept to the covariate matrix
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
#' @param x  A matrix of dimension h x p containing the covariates data,
#'  where h is the number of observations (equal to the length of y) and
#'  p represent the dimension the covariates in the component models.
#' Note that the the first column is a column  of 1s which is added
#' in order to allow for the intercept in the model.
#' @param z  A matrix of dimension h x q containing the covariates data,
#' where h is the number of observations (equal to the length of y) and
#' q represent the dimension the covariates in the mixture weights model.
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
#' @param x  A matrix of dimension h x p containing the covariates data,
#' where h is the number of observations (equal to the length of y) and
#' p represent the dimension the covariates in the component models.
#' @param z  A matrix of dimension h x q containing the covariates data,
#' where h is the number of observations (equal to the length of y) and
#' q represent the dimension the covariates in the mixture weights model.
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


#' proposal for the generalized Poisson component model
#'
#' this function computes the moments (mean and covariance) of the proposal distribution
# '
#' @param y: a vecotor of the response values
#' @param design_matrix: a matrix of the covariate values augmented by a column of 1s
#' @param mean_effect: current value of the mean function of the generalized poisson model
#' @param disp_effect: current value of the overdispersion (scale) function of the generalized poisson model
#' @param u_mean: transition matrix of the regression parameters on the mean function
#' @param u_disp: transition matrix of the regression parameters on the overdispersion function.

#' @returns : updated mean_effect, disp_effect, cov_mean and cov_disp.


gp_linearbayes <- function(y, design_matrix, mean_effect, disp_effect, u_mean, u_disp) {
  for (r in seq_len(length(y))) { # loop over obs in the interval
    # prior settings
    # mean level
    cov_val_mean <- u_mean %*% design_matrix[r, ] # covariance value from the mean predictor
    pvar_mean <- as.numeric(design_matrix[r, ] %*% cov_val_mean) # prior variance mean predictor
    prior_mean_pred <- design_matrix[r, ] %*% mean_effect # expected prior of the effect of the mean predictor
    # dispersion level
    cov_val_disp <- u_disp %*% design_matrix[r, ] # covariance value from the dipersion predictor
    pvar_disp <- as.vector(design_matrix[r, ] %*% cov_val_disp) # prior variance for dispersion predictor
    prior_disp_pred <- design_matrix[r, ] %*% disp_effect # expected prior of the effect of the mean predictor

    # updating hyperprior values
    a <- 1 + pvar_disp^-1
    b <- pvar_disp^-1 * exp(-prior_disp_pred)

    # implement the laplace approximation
    mode_mean <- ifelse(y[r] == 0, exp(prior_mean_pred), y[r])
    mode_disp <- ifelse(y[r] == 0, exp(prior_disp_pred), (sqrt(((a - 2) * y[r] - b)^2 +
                                                                 4 * (a - 1) * b * y[r]) + (a - 2) * y[r] - b) / (2 * b * y[r]))
    mode_disp <- ifelse(mode_disp == 0, a - 1, mode_disp) # adjust for singularities in the mode value for the dispersion parameter
    var_mean <- sqrt((1 + mode_disp * y[r])^2 / mode_mean)
    var_disp <- sqrt((1 + mode_disp * y[r])^2 / (mode_disp * (b * (1 + mode_disp * y[r])^2 + y[r])))

    # update the proposal moments
    mean_effect <- mean_effect + (cov_val_mean / pvar_mean) * as.numeric(log(mode_mean) - prior_mean_pred)
    disp_effect <- disp_effect + (cov_val_disp / pvar_disp) * as.numeric(log(mode_disp) - prior_disp_pred)
    u_mean <- u_mean + cov_val_mean %*% t(cov_val_mean) * as.numeric((var_mean - pvar_mean) / pvar_mean^2)
    u_disp <- u_disp + cov_val_disp %*% t(cov_val_disp) * as.numeric((var_disp - pvar_disp) / pvar_disp^2)
  }
  return(list(mean_effect = mean_effect, disp_effect = disp_effect, cov_mean = u_mean, cov_disp = u_disp))
}

