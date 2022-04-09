#' robust_inverse
#'
#' Compute the inverse of a covariance matrix.
#' @param mat the covariance matrix
#' @return The inverse of the covariance matrix.

robust_inverse <- function(mat){

      chol_transfm <- chol(mat)
      inv_mat <- solve(chol_transfm)%*%solve(t(chol_transfm)) + diag(1E-5, nrow(mat))
  return(0.5*(inv_mat + t(inv_mat)))
}

#' prior_moments
#'
#' Compute the prior moments of the linear predictor.
#'
#' Function for computing the prior moments of the degenerate prior distribution of the linear predictor parameter.
#'
#' @param x  A matrix of dimension h x (p+1) containing the covariates data, where h is the number of observations (equal to the length of y) and p represent the dimension the covariates in the component models.
#' Note that the the first column is a column  of 1s which is added in order to allow for the intercept in the model.
#' @param z  A matrix of dimension h x (q+1) containing the covariates data, where h is the number of observations (equal to the length of y) and q represent the dimension the covariates in the mixture weights model.
#' Note that the the first column is a column  of 1s which is added in order to allow for the intercept in the model.
#' @param param The prior mean of the regression coefficients in both the component and mixture weight models.
#' @param cov_mat The prior covariance matrix of the regression coefficients in both the component and mixture weight models.
#' @param n_comp Number of mixture components.
#'
#' @return The prior mean and variance of the linear predictor. Also, the covariance of the linear predictor and the regression coefficients is returned.
prior_moments <-function(x, z, param, cov_mat, n_comp){

  x_mat <- design_matrix(x,z,n_comp) # design matrix
  linear_pred <- param %*% x_mat
  covar_linear_pred <- cov_mat %*% x_mat
  var_linear_pred <- t(x_mat) %*% covar_linear_pred
  var_linear_pred <- 0.5*(var_linear_pred + t(var_linear_pred)) + diag(1E-5, nrow(var_linear_pred))
  return(list(linear_pred=linear_pred, covar_linear_pred=covar_linear_pred, var_linear_pred=var_linear_pred ))
}

#' Computing the component allocation probabilities
#' @param y A vector of response variable values.
#' @param linear_pred A vector of the linear predictors in both components and mixture weights models.
#' @param n_comp Number of mixture components.
#'
#' @return  The components allocation probabilities.
compute_allocation_prob <- function(y,linear_pred,n_comp){
  lambda <- exp(linear_pred[seq_len(n_comp)])
  log_mixtureweight <- dpois(y,lambda,log=TRUE)+c(0,linear_pred[-seq_len(n_comp)])
  alloc_prob <- normalize(log_mixtureweight)
  # avoid some allocation prob to be completely zero
  if (max(alloc_prob)>.95){
    alloc_prob[-which.max(alloc_prob)] <- alloc_prob[-which.max(alloc_prob)] + .05/(n_comp-1)
    alloc_prob[which.max(alloc_prob)] <- alloc_prob[which.max(alloc_prob)] - .05

  }
  return(alloc_prob)
}


#' Likelihood gradients
#'
#' Computing the first derivative of the mixture model's likelihood with respect to the linear predictors.
#'
#' @param y A vector of response variable values.
#' @param linear_pred A vector of the linear predictors in both components and mixture weights models.
#' @param n_comp Number of mixture components.
#' @return  A vector of the likelihood gradients.
likelihood_gradients <- function(y, linear_pred, n_comp){

  # gradient w.r.t to tw.r.t the component model's linear predictors.
  gradient_components <- diag(y - exp(linear_pred[1:n_comp]))
  # gradient w.r.t mixture weights linear predictors.
  linear_pred_weight <- linear_pred[(n_comp+1):length(linear_pred)] # linear predictors for the mixture weights
  weight_normalization_const <- 1 + sum(exp(linear_pred_weight))
  gradient_mix_weights <- matrix(0,n_comp-1, n_comp-1)

  for(k in seq_len(n_comp-1)){
    gradient_mix_weights[k,] <- -exp(linear_pred_weight)/weight_normalization_const
    gradient_mix_weights[k,k] <- 1+gradient_mix_weights[k,k]
  }

  likelihood_gradients_matrix <- cbind(gradient_components, rbind(0,gradient_mix_weights))

  return (likelihood_gradients_matrix)
}

#' Likelihood hessians
#'
#' Computing the second derivative of the mixture model's likelhood with respect to the linear predictors.
#' @param linear_pred A vector of the linear predictors in both components and mixture weights models.
#' @param component An index between 1 and n_comp indicating with component the derivatives are computed for.
#' @param n_comp Number of mixture components.
#' @return  The likelihood hessian matrix for given component index.

likelihood_hessians <- function(linear_pred, component, n_comp){
  likelihood_hessian <- matrix(0,length(linear_pred), length(linear_pred))
  if (component>1){
    # second derivatives w.r.t the component model's linear predictors
    likelihood_hessian[component, component] <- -exp(linear_pred[component])
    # second derivatives w.r.t the mixture weights linear predictors.
    hessian_mix_weights <- matrix(0, n_comp-1, n_comp-1)
    linear_pred_weight <- linear_pred[(n_comp+1):length(linear_pred)] # linear predictors for the mixture weights
    normalization_constant <- 1 + sum(exp(linear_pred_weight))
    for (k in seq_len(n_comp-1)){
      hessian_mix_weights[k,] <- exp(linear_pred_weight + linear_pred_weight[k])
      hessian_mix_weights[k, k] <- -exp(linear_pred_weight[k])*(normalization_constant - exp(linear_pred_weight[k]))
    }

    likelihood_hessian[(n_comp + 1):length(linear_pred), (n_comp + 1):length(linear_pred)] <- hessian_mix_weights/normalization_constant^2

  }else{
    likelihood_hessian[component, component] <- -exp(linear_pred[component])
  }
  return(likelihood_hessian)
}


#' Posterior gradients.
#'
#' Computing the posterior gradients of the linear predictors.
#'
#' @param allocation_prob A vector of component allocation probabilities.
#' @param y A vector of response variable values.
#' @param linear_pred A vector of the linear predictors in both components and mixture weights models.
#' @param n_comp Number of mixture components.
#' @return  A vector of posterior gradients.
gradient_posterior <- function (allocation_prob, y, linear_pred, n_comp){

  gradient_vec <- allocation_prob %*% likelihood_gradients(y = y,
                                                           linear_pred = linear_pred,
                                                           n_comp = n_comp)
  return(gradient_vec)
}

#' Posterior hessian.
#'
#' Computing the posterior gradients with respect to the linear predictors.
#'
#' @param allocation_prob A vector of component allocation probabilities.
#' @param posterior_gradients A vector of posterior gradients of the linear predictors.
#' @param prior_var_inv Inverse of the covariance matrix of the prior for the regression coefficients.
#' @param y A vector of response variable values.
#' @param linear_pred A vector of the linear predictors in both components and mixture weights models.
#' @param n_comp Number of mixture components.
#' @return  The posterior hessian matrix.

hessian_posterior <- function (allocation_prob, posterior_gradients, prior_var_inv,
                               y, linear_pred, n_comp){

  likelihood_gradient_matrix <- likelihood_gradients(y = y,
                                                      linear_pred = linear_pred,
                                                      n_comp = n_comp)

  hessian_matrix <- matrix(0, length(linear_pred), length(linear_pred))

  for (k in seq_len(n_comp)){

    hessian_matrix <- hessian_matrix +
      allocation_prob[k] * likelihood_hessians(linear_pred=linear_pred,
                                                 component=k,
                                                  n_comp=n_comp)
  }
  hessian_matrix <- prior_var_inv - hessian_matrix

  return(0.5*(hessian_matrix + t(hessian_matrix)))
}

#' Local linear approximation of the posterior.
#'
#' Computes the first and second moments of the linear approximation of the posterior using Newton's method.
#'
#' @param y A vector of response variable values.
#' @param linear_pred A vector of the linear predictors in both components and mixture weights models.
#' @param prior_var_inv Inverse of the covariance matrix of the prior for the regression coefficients.
#' @param n_comp Number of mixture components.
#' @return  The the first and second moments of the linear approximation of the posterior.

posterior_linear_moments<- function(y, linear_pred, prior_var_inv, n_comp){
  expected_allocation <- compute_allocation_prob(y,linear_pred,n_comp)
  # update the component linear predictors
  posterior_grad_matrix<- gradient_posterior(allocation_prob =expected_allocation,
                                                     linear_pred = linear_pred,
                                                     y = y,
                                                     n_comp = n_comp)
  posterior_hessian_matrix <- hessian_posterior(allocation_prob = expected_allocation,
                                                posterior_gradients = posterior_grad_matrix,
                                                prior_var_inv = prior_var_inv,
                                                y = y,
                                                linear_pred = linear_pred,
                                                n_comp = n_comp
                                                )
  posterior_var <- robust_inverse(posterior_hessian_matrix)
  posterior_location <- linear_pred + posterior_grad_matrix %*% posterior_var
  return(list(location=posterior_location,variance = posterior_var))
}

