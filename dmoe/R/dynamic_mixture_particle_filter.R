
compute_reg_coef_moment <- function(theta, weights){
  mean_theta <- colSums(theta * weights)
  error_vec <- theta - matrix(rep(mean_theta, each=nrow(theta)), ncol= length(mean_theta))
  if(length(weights)<30000){
  var_theta <- t(error_vec) %*% diag(weights) %*%  error_vec
  }else{
    var_theta <- var(theta)
  }
  var_theta <- 0.5*(var_theta + t(var_theta)) + diag(1E-5, nrow(var_theta))

  return(list(mean = mean_theta, var = var_theta))
}
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
mixture_likelihood<-function(y, x, z, param, n_comp){
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

#  function for computing the importance weights
compute_importance_weight<-function(y, x, z, prev_weights, param,
                                    prev_particles, prior_var, proposal_mean,
                                    proposal_cov, n_comp){

  #dim_par<-nrow(prior_var)
  dens_proposal<-dmvnorm(param, proposal_mean, proposal_cov, log=TRUE)
  dens_prior<-dmvnorm(prev_particles, param, prior_var,log=TRUE)
  p_prior<-log(sum(exp(dens_prior)*prev_weights))
  log_weight<-mixture_likelihood(y, x, z, param, n_comp) + p_prior - dens_proposal
  return(log_weight)
}

# particle filter algorithm
FAPF<-function(Data, intervals, particle_size, mix_col=NULL, exp_col=NULL,
               n_comp=2, time_ind=1, response_ind=2, alpha=0.45, R=2,
               return_allweights=F,h=1){

  interval_length<-length(intervals)-1 # number of time periods
  pb <- txtProgressBar(min=0,max=interval_length,initial=0,char="_",style = 3)
  # set.seed(600) # set seeds

  # check if required package is installed
  #if (!require("pacman")) install.packages("pacman",quietly=TRUE)
  pacman::p_load(mvtnorm, MCMCpack,LaplacesDemon) # load required packages

  # Variable declaration and initialization

  design_dim <- n_comp*(length(exp_col)+1) + (n_comp-1)*(length(mix_col)+1)# dimension of the design matrix
  particles<-list() # save sampled particles
  importance_weights<-rep(1/particle_size, particle_size) # initialize weights
  #cov_mat<-1*diag(design_dim)# set the initial covariance
  theta_prev<-rmvnorm(particle_size,rep(0,design_dim),1*diag(design_dim)) # sample theta at time point j=0
  theta_moments <- compute_reg_coef_moment(theta_prev, importance_weights)
  mean_theta <- theta_moments$mean
  allWeights<-list()

  # propagate paricles at time point j>1
  ess<-c() # effective sample size
  n_unique<-c()
  score<-0
  for(j in 2:(interval_length+1)){

    # generating risk set.
    Dj<-Data[which(Data[,time_ind]<=intervals[j]&Data[,time_ind]>intervals[j-1]),] #interval data
    yj <- Dj[,response_ind]
    # effect covariates
    xj <- as.matrix(Dj[,c(exp_col)]) # covariate in the component models
    zj <- as.matrix(Dj[,c(mix_col)]) # covariate in the mixture weight models
    # proposal distribution moments

    prop_moments <- tryCatch ({
      multiple_component_proposal(yj, xj, zj, theta_moments$mean,
                                     theta_moments$var*(alpha^-1 - 1), n_comp)
      },
      error = function(err){
        print("something went wrong in the proposal density.Proposing from the prior")
        print(err)
        return(multiple_component_proposal(yj, xj, zj, theta_moments$mean,
                                theta_moments$var, n_comp))
      })

    # propose new particles
    theta<-mvrnorm(particle_size, prop_moments$mean_value, prop_moments$cov_value)

    # compute weights at time j
    newimportance_weights<-apply(as.matrix(1:particle_size),1,function(m)
      compute_importance_weight(yj, xj, zj,
                                importance_weights, theta[m,],
                                theta_prev, theta_moments$var*(alpha**-1 -1),
                                prop_moments$mean_value,
                                prop_moments$cov_value,
                                n_comp))#denominator
    newimportance_weights[which(is.na(newimportance_weights))]<--log(1e-300)
    newimportance_weights[is.infinite(newimportance_weights)]<-log(1e-300)

    # normalize weights
    newimportance_weights<-normalize(newimportance_weights)
    # indexes of the resampled particles
    newsample<-reject_strat_resample(newimportance_weights,particle_size/R)
    index<-newsample$index # indexes of sampled particles
    importance_weights<-newsample$weights

    theta_moments <- compute_reg_coef_moment(theta[index,], importance_weights)
    allWeights[[j-1]]<- importance_weights
    # number of unique particles
    n_unique[j-1]<-length(unique(index))
    # compute the effective sample size
    ess[j-1]<-1/sum(importance_weights^2)
    # compute the log predictive score
    if(j > floor((interval_length-1)/2)){
      log_pred_dens<-apply(as.matrix(1:nrow(theta_prev)),1,function(m)
        mixture_likelihood(yj, xj, zj, theta_prev[m,], n_comp))
      score<- score +log(sum(as.numeric(exp(log_pred_dens))*allWeights[[j-2]]))
    }
    # reset the previous paramater vector
    theta_prev <- theta[index,]
    # save the resampled particles
    particles[[j-1]] <-  theta[index,]
    setTxtProgressBar(pb, value=j)
  }

  if(return_allweights){
    return(list(Particles=particles,ESS=ess,LPS=2*score/(interval_length-1),
                weights=allWeights,n_unique=n_unique, CovMat=theta_moments$var))
  }else{
    return(list(Particles=particles,ESS=ess,LPS=2*score/(interval_length-1),
                weights=newimportance_weights,n_unique=n_unique,CovMat=theta_moments$var))
  }
}

