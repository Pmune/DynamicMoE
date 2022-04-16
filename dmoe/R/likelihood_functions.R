# function for preparing the design matrix
design_matrix <- function(x, z, n_comp){
  x_tild <- kronecker(diag(n_comp),c(1,x)) # expanded covariate for component models
  z_tild <- kronecker(diag(n_comp-1),c(1,z)) # expanded covariate for mixture models
  d_x <- dim(x_tild) # dimension of x_tild
  d_z <- dim(z_tild) #dimension of z_tild
  d <- d_x+d_z
  design_mat <- matrix(0,nrow=d[1],ncol = d[2])
  design_mat[seq_len(d_x[1]),seq_len(d_x[2])] <- x_tild
  design_mat[-seq_len(d_x[1]),-seq_len(d_x[2])] <- z_tild
  return(design_mat)
}

# marginal likelihood for the mixture model
mixture_log_likelihood<-function(y, x, z, param, n_comp){
  log_likelihood<-0
  for(i in seq_len(length(y))){
    x_mat <- design_matrix(x[i,],z[i,],n_comp) # design matrix
    eta <- param %*% x_mat
    lambda <- exp(eta[seq_len(n_comp)]) # expected effect
    psi<-c(1,exp(eta[-seq_len(n_comp)])) # un-normalized mixing weights
    marginal_densities<-dpois(y[i],lambda)*psi/sum(psi)
    log_likelihood<-log_likelihood+log(sum(marginal_densities))
  }

  return(log_likelihood)
}

# likelihood for dynamic poisson model
poisson_log_likelihood<-function(y,X,param)
{
  log_likelihood<-0
  for(i in seq_len(length(y)))
  {
    lambda<-exp(X[i,]%*%param) # expected effect
    log_likelihood<-log_likelihood + dpois(y[i],lambda,log=TRUE)
  }
  return(log_likelihood)
}
