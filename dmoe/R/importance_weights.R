#  function for computing the importance weights
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
  return(log_weight)
}
