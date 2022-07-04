#' function for preparing the design matrix for mixture of experts models
#' @param x A vector of covariates for components models.
#' The vector doesn't include the intercept.
#' @param z A vector of covariates for mixture weights models.
#' The vector doesn't include the intercept.
#' @param n_comp Number of mixture components.
#'
#' @returns The design matrix.
#'
design_matrix <- function(x, z, n_comp){
  x_tild <- kronecker(diag(n_comp),c(1,x)) # expanded covariate for component models
  z_tild <- kronecker(diag(n_comp-1),c(1,z)) # expanded covariate for mixture models
  d_x <- dim(x_tild)
  d_z <- dim(z_tild)
  d <- d_x+d_z
  design_mat <- matrix(0,nrow=d[1],ncol = d[2])
  design_mat[seq_len(d_x[1]),seq_len(d_x[2])] <- x_tild
  design_mat[-seq_len(d_x[1]),-seq_len(d_x[2])] <- z_tild
  return(design_mat)
}


#' Log likelihood for the mixture of experts model.
#'
#' @param y A vector of response variable values.
#' @param x  A matrix of dimension n x (p+1) containing the covariates data,
#' where n is the number of observations (equal to the length of y) and
#' p represent the dimension the covariates in the component models.
#' Note that the the first column is a column  of 1s which is added in order
#' to allow for the intercept in the model.
#' @param z  A matrix of dimension n x (q+1) containing the covariates data,
#' where n is the number of observations (equal to the length of y) and
#' q represent the dimension the covariates in the mixture weights model.
#' Note that the the first column is a column  of 1s which is added
#' in order to allow for the intercept in the model
#' @param reg_coef Vector of the regression coefficients in both the component
#' and mixture weight models.
#' @param n_comp Number of mixture components.
#' @export

mixture_log_likelihood<-function(y, x, z, reg_coef, n_comp){
  log_likelihood<-0
  for(i in seq_len(length(y))){
    x_mat <- design_matrix(x[i,],z[i,],n_comp) # design matrix
    eta <- reg_coef %*% x_mat
    lambda <- exp(eta[seq_len(n_comp)]) # expected effect
    psi<-c(1,exp(eta[-seq_len(n_comp)])) # un-normalized mixing weights
    marginal_densities<-dpois(y[i],lambda)*psi/sum(psi)
    log_likelihood<-log_likelihood+log(sum(marginal_densities))
  }

  return(log_likelihood)
}


#' Log likelihood for Poisson regression model.
#' @param y A vector of response variable values.
#' @param x  A matrix of dimension n x (p+1) containing the covariates data,
#' where n is the number of observations (equal to the length of y) and
#' p represent the dimension the covariates. Note that the the first column
#' is a column  of 1s which is added in order to allow for the intercept.
#' @param reg_coef Vector of the regression coefficients
#' @export

poisson_log_likelihood<-function(y, X, reg_coef)
{
  log_likelihood<-0
  for(i in seq_len(length(y)))
  {
    lambda<-exp(X[i,]%*%reg_coef) # expected effect
    log_likelihood<-log_likelihood + dpois(y[i],lambda,log=TRUE)
  }
  return(log_likelihood)
}
