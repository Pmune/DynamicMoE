#' robust_inverse
#'
#' Compute the inverse of a covariance matrix.
#' @param mat the covariance matrix
#' @return The inverse of the covariance matrix.

robust_inverse <- function(mat){
  mat = 0.5 * (mat + mat) # input matrix need to be symmetric
  chol_transfm <- chol(mat)
  inv_mat <- solve(chol_transfm)%*%solve(t(chol_transfm)) + diag(1E-5, nrow(mat))
  return(0.5*(inv_mat + t(inv_mat)))
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
    alloc_prob <- normalize(log(alloc_prob))
  }

  return(alloc_prob)
}


#' Likelihood gradients
#'
#' Computing the first derivative of the mixture model's likelihood with respect to the linear predictors.
#'
#' @param y A  response value.
#' @param linear_pred A vector of the linear predictors in both components and mixture weights models.
#' @param n_comp Number of mixture components.
#' @return  A matrix of the likelihood gradients.
likelihood_gradients <- function(y, linear_pred, n_comp){

  # gradient w.r.t to tw.r.t the component model's linear predictors.
  gradient_components <- diag(y - exp(linear_pred[1:n_comp]))
  # gradient w.r.t mixture weights linear predictors.
  linear_pred_weight <- linear_pred[(n_comp+1):length(linear_pred)] # linear predictors for the mixture weights
  weight_normalization_const <- 1 + sum(exp(linear_pred_weight))
  gradient_mix_weights <- matrix(0,n_comp-1, n_comp-1)

  for(k in seq_len(n_comp-1)){
    gradient_mix_weights[k,] <- -exp(linear_pred_weight) / weight_normalization_const
    gradient_mix_weights[k,k] <- 1 + gradient_mix_weights[k,k]
  }

  likelihood_gradients_matrix <- cbind(gradient_components, rbind(0,gradient_mix_weights))

  return (likelihood_gradients_matrix)
}

#' Likelihood hessian matrix
#'
#' Computing the second derivatives of the mixture model's likelihood with respect to the linear predictors of a given component.
#' @param linear_pred A vector of the linear predictors in both components and mixture weights models.
#' @param component An index between 1 and n_comp indicating which component the derivatives are computed for.
#' @param n_comp Number of mixture components.
#' @return  The likelihood hessian matrix for a given component index.

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
    likelihood_hessian[component, component] <- - exp(linear_pred[component])
  }
  return(likelihood_hessian)
}


#' Posterior gradient with respect to linear predictors.
#'
#' Computing the posterior gradients of the linear predictors.
#'
#' @param allocation_prob A vector of component allocation probabilities.
#' @param y  A response value.
#' @param linear_pred A vector of the linear predictors in both components and mixture weights models.
#' @param n_comp Number of mixture components.
#' @return  A vector of posterior gradients.
gradient_linear_pred <- function (allocation_prob, y, linear_pred, n_comp){

  gradient_vec <- allocation_prob %*% likelihood_gradients(y,linear_pred,
                                                           n_comp)
  return(gradient_vec)
}

#' Hessian matrix with respect to linear predictors.
#'
#' Computing the posterior gradients with respect to the linear predictors.
#'
#' @param allocation_prob A vector of component allocation probabilities.
#' @param y A response value.
#' @param y A vector of response variable values.
#' @param linear_pred A vector of the linear predictors in both components and mixture weights models.
#' @param n_comp Number of mixture components.
#' @return  The hessian matrix with respect to the linear predictors.

hessian_linear_pred <- function (allocation_prob, y, linear_pred, n_comp){
  hessian_matrix <- matrix(0, length(linear_pred), length(linear_pred))

  for (k in seq_len(n_comp)){

    hessian_matrix <- hessian_matrix + (allocation_prob[k] *
                                          likelihood_hessians(linear_pred, k,
                                                              n_comp))
  }
  hessian_matrix <- - hessian_matrix

  return(0.5*(hessian_matrix + t(hessian_matrix)))
}


gradient_reg_coeff <- function(x_mat, grad_linear_pred){

  return (as.vector(x_mat %*% grad_linear_pred))
}


hessian_reg_coeff <- function(x_mat, hess_linear_pred){
  return (x_mat %*% hess_linear_pred %*% t(x_mat))
}

#' Update the first and second moments of regression coefficients given
#' their hessian  and the gradient.
#'
#' @param  prior_mean Prior mean of the regression coefficients.
#' #' @param inv_prior_var Inverse of the covariance matrix of the regression coefficients.
#' @param gradient The gradient of the regression coefficients.
#' @param hessian_matrix The Hessian matrix of the regression coefficients.
#'
#' @return  The updated location and covariance matrix of the regression coefficients.

update_linear_moments<- function(prior_mean, inv_prior_var,
                                 gradient, hessian_matrix){
  updated_var <- robust_inverse(hessian_matrix + inv_prior_var)
  updated_location <- prior_mean + gradient %*% updated_var

   return(list(location = as.vector(updated_location), variance = updated_var))
  }


