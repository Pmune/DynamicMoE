
particles_update_multicomp <- function(y, x, z, n_comp, sample_size, particles,
                                       importance_weights, alpha,
                                       proposal_method="linearbayes"){
  prior_moments <- reg_coef_prior_moment(particles, importance_weights)
  if (proposal_method == "linearbayes"){
    proposal_density = linear_bayes_proposal_multicomp
  } else{
    proposal_density = local_linear_proposal_multicomp
  }
  # proposal density moments
  prop_moments <- tryCatch ({
    proposal_density(y, x, z, prior_moments$mean,
                     prior_moments$var*(alpha**-1), n_comp)
  },
  error = function(err){
    print("something went wrong in the proposal density. Recomputing again the proposal density moments")
    print(err)
    return(proposal_density(y, x, z, prior_moments$mean,

                            prior_moments$var*(alpha**-1 -1), n_comp))
  })

  new_particles<-rmvnorm(sample_size, prop_moments$mean, prop_moments$cov)

  newimportance_weights<-apply(as.matrix(1:sample_size),1,function(m)
    importance_weight_multicomp(y, x, z, new_particles[m,], importance_weights,
                                particles,prior_moments$mean,
                                prior_moments$var*(alpha**-1-1),
                                prop_moments$mean,
                                prop_moments$cov, n_comp))
  newimportance_weights[which(is.na(newimportance_weights))]<--log(1e-300)
  newimportance_weights[is.infinite(newimportance_weights)]<-log(1e-300)

  log_likelihood <-apply(as.matrix(1:nrow(particles)),1,function(m)
    mixture_log_likelihood(y, x, z, particles[m,], n_comp))
  log_pred_dens =  log(sum(as.numeric(1e-300 + exp(log_likelihood))*importance_weights))

  return(list(particles = new_particles, weights = newimportance_weights,
              log_pred_dens = log_pred_dens ))
}


particles_update_1comp <- function(y, x, sample_size, particles, importance_weights, alpha,
                                   proposal_method="linearbayes"){

  prior_moments <- reg_coef_prior_moment(particles, importance_weights)

  if (proposal_method == "linearbayes"){
    proposal_density = linear_bayes_proposal_1comp
  } else{
    proposal_density = local_linear_proposal_1comp
  }
  # proposal density moments
  prop_moments <- tryCatch ({
    proposal_density(y, x, prior_moments$mean, prior_moments$var*(alpha**-1))
  },
  error = function(err){
    print("something went wrong in the proposal density. Recomputing again the proposal density moments")
    print(err)
    return(proposal_density(y, x, prior_moments$mean,

                            prior_moments$var*(alpha**-1 -1)))
  })

  new_particles<-rmvnorm(sample_size, prop_moments$mean, prop_moments$cov)

  newimportance_weights<-apply(as.matrix(1:sample_size),1,function(m)
    importance_weight_1comp(y, x, new_particles[m,], importance_weights,
                            particles,prior_moments$mean,
                            prior_moments$var*(alpha**-1-1),
                            prop_moments$mean, prop_moments$cov))
  newimportance_weights[which(is.na(newimportance_weights))]<--log(1e-300)
  newimportance_weights[is.infinite(newimportance_weights)]<-log(1e-300)

  log_likelihood<-apply(as.matrix(1:nrow(particles)),1,function(m)
    poisson_log_likelihood(y, x, particles[m,]))
  log_pred_dens =  log(sum(as.numeric(1e-300 + exp(log_likelihood))*importance_weights))

  return(list(particles = new_particles, weights = newimportance_weights,
              log_pred_dens = log_pred_dens ))
}
