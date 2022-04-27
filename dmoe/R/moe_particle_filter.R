

#' particle filter algorithm for dynamic mixture of experts
#' @export
#' @import LaplacesDemon mvtnorm

dmoe<-function(Data, intervals, particle_size, mix_col=NULL, exp_col=NULL,
               n_comp=1, v_init = 1, time_ind=1, response_ind=2, alpha=0.45, R=2,
               proposal_method="linearbayes", mpf=T, return_all=F){

  #pacman::p_load(mvtnorm, LaplacesDemon) # load required packages
  interval_length<-length(intervals)-1 # number of time periods
  pb <- txtProgressBar(min=0,max=interval_length,initial=0,char="_",style = 3)
  # dimension of the design matrix
  design_dim <- n_comp*(length(exp_col)+1) + (n_comp-1)*(length(mix_col)+1)
  particles_sample<-list() # save sampled particles
  ess<-c() # effective sample size
  log_pred_score<-0 # log predictive score

  # Initialization step

  importance_weights<-rep(1/particle_size, particle_size) # initialize weights
  # sample particles (regression coefficients) at time point j=0
  particles<-rmvnorm(particle_size,rep(0,design_dim), v_init*diag(design_dim))
  allWeights<-list()

  # propagate paricles at time point j>1
  for(j in 2:(interval_length+1)){
    start_time <- Sys.time()
    D <-Data[which(Data[,time_ind] <= intervals[j] &
                     Data[,time_ind] > intervals[j-1]),] #interval data
    y <- D[,response_ind]
    x <- as.matrix(D[,c(exp_col)]) # covariate in the component models
    if (n_comp > 1){
      z <- as.matrix(D[,c(mix_col)]) # covariate in the mixture weight models
      if(mpf){
        particles_update = mpf_update_multicomp(y, x, z, n_comp,
                                                particle_size,
                                                particles,
                                                importance_weights,
                                                alpha,
                                                proposal_method)

      } else{
        particles_update = particles_update_multicomp(y, x, z, n_comp,
                                                      particle_size,
                                                      particles,
                                                      importance_weights,
                                                      alpha,
                                                      proposal_method)
      }

    } else{

      if(mpf){
        particles_update = mpf_update_1comp(y, cbind(1,x),
                                            particle_size,
                                            particles,
                                            importance_weights,
                                            alpha,
                                            proposal_method)

      } else{
        particles_update = particles_update_1comp(y, cbind(1,x),
                                                  particle_size,
                                                  particles,
                                                  importance_weights, alpha,
                                                  proposal_method)
      }

    }

    newimportance_weights<-normalize(particles_update$weights)

    # resampling step

    newsample<-reject_strat_resample(newimportance_weights,particle_size/R)
    index<-newsample$index # indexes of sampled particles
    importance_weights<-newsample$weights
    end_time <- Sys.time()
    run_time <- end_time - start_time
    run_time <- ifelse(attr(run_time, "units") == "secs", as.numeric(run_time),
                       ifelse(attr(run_time, "units") == "mins",
                              as.numeric(run_time) * 60,
                              as.numeric(run_time) * 3600))
    # compute the effective sample size
    cat(paste0("ESS for Interval ", j, " = ", 1/sum(importance_weights^2)))
    ess[j-1]<-(1/sum(importance_weights^2))/run_time
    # compute the log predictive score
    if(j > floor((interval_length-1)/2)){
      log_pred_score <- log_pred_score + particles_update$log_pred_dens
    }
    # reset the previous paramater vector
    particles <- particles_update$particles[index,]
    # save the sampled particles
    if(return_all){
      particles_sample[[j-1]] <- particles
      allWeights[[j-1]] <- importance_weights
    }
    setTxtProgressBar(pb, value=j)
  }

  if(return_all){
    return(list(particles=particles_sample, ess=ess,lps=2*log_pred_score/(interval_length-1),
                weights=allWeights))
  }else{
    return(list(particles=particles,ess=ess,lps=2*log_pred_score/(interval_length-1),
                weights=importance_weights))
  }
}
