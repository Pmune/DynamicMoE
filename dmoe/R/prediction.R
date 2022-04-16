
dmixpois<-function(y,lambda,Mix.wgt){

  dens<-diag(dpois(y,lambda)%*%t(Mix.wgt))
  return(dens)
}

#' Predictive density for dynamic mixture of experts.
#' @export
predictions_mutlicomp<-function(test_data,particles,importance_weight,y_max,exp_col,mix_col,n_comp){
  d.ftd<-c()
  for(i in 1:nrow(test_data)){
    x <- as.numeric(test_data[i,exp_col])
    z <- as.numeric(test_data[i,mix_col])
    x_mat <- design_matrix(x,z,n_comp) # design matrix
    eta <- particles %*% x_mat
    lambda <- exp(eta[,seq_len(n_comp)]) # expected effect
    psi<-cbind(1,exp(eta[,-seq_len(n_comp)])) # un-normalized mixing weights
    Mix.wgt<-psi/rowSums(psi)
    # compute the predictive density
    pmf<-apply(as.matrix(c(0:y_max)),1,function(y)sum(dmixpois(y,lambda,Mix.wgt)*importance_weight))
    pmf<- pmf/sum(pmf) # normalize the density
    d.ftd<-rbind(d.ftd,pmf)
  }
  d.ftd<-colMeans(d.ftd)
  return(d.ftd)
}

#' Predictive density for dynamic Poisson regression.
#' @expor
predictions_1Comp<-function(test_data,Particles,ImportanceWeight,y.max,exp.col){

  d.ftd<-c()

  for(m in 1:nrow(test_data)){
    X<-cbind(1,test_data[m,exp.col])
    lambda<-exp(Particles%*%t(X))
    # compute the predictive density
    pmf<-apply(as.matrix(c(0:y.max)),1,function(y)sum(dpois(y,lambda)*ImportanceWeight))
    pmf<- pmf/sum(pmf) # normalize the density
    d.ftd<-rbind(d.ftd,pmf)

  }

  d.ftd<-colMeans(d.ftd)
  return(d.ftd)
}
