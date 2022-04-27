#  function for computing the importance weights
<<<<<<< HEAD
importance_weight_multicomp<-function(y, x, z, param, prev_weights,
                                      prev_particles, prior_mean, prior_var,
                                      proposal_mean, proposal_cov, n_comp){

  dens_proposal <- dmvnorm(param, proposal_mean, proposal_cov, log=TRUE)
  #dens_prior <- dmvnorm(param, prior_mean, prior_var, log=TRUE)
  dens_prior<-log(sum(dmvnorm(prev_particles, param, prior_var)*prev_weights))
  log_weight <- mixture_log_likelihood(y, x, z, param, n_comp) + dens_prior -
                                 dens_proposal
  return(log_weight)
}


importance_weight_1comp<-function(y, x, param, prev_weights, prev_particles,
                                  prior_mean, prior_var, proposal_mean,
                                  proposal_cov){

  dens_proposal <- dmvnorm(param, proposal_mean, proposal_cov, log=TRUE)
  #dens_prior <- dmvnorm(param, prior_mean, prior_var, log=TRUE)
  dens_prior<-log(sum(dmvnorm(prev_particles, param, prior_var)*prev_weights))
  log_weight <- poisson_log_likelihood(y, x, param) + dens_prior - dens_proposal
=======
importance_weight_multicomp<-function(y, x, z, prev_weights, param,
                                    prev_particles, prior_var, proposal_mean,
                                    proposal_cov, n_comp){

  dens_proposal <- dmvnorm(param, proposal_mean, proposal_cov, log=TRUE)
  dens_prior <- dmvnorm(prev_particles, param, prior_var)
  log_weight <- mixture_log_likelihood(y, x, z, param, n_comp) +
    log(sum(dens_prior * prev_weights)) - dens_proposal
  return(log_weight)
}

# function for computing the importance weights
importance_weight_1comp<-function(y, x, prev_weights, param, prev_particles,
                             prior_var, proposal_mean, proposal_cov){

  dens_proposal <- dmvnorm(param, proposal_mean, proposal_cov, log=TRUE)
  dens_prior <- dmvnorm(prev_particles, param, prior_var)
  log_weight <- poisson_log_likelihood(y, x, param) +
    log(sum(dens_prior * prev_weights)) - dens_proposal
>>>>>>> b448a75bab05897aa5955b3f5b730ddb72c5a5b1
  return(log_weight)
}
