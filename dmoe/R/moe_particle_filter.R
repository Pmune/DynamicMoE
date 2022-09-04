

#' particle filter algorithm for dynamic mixture of experts
#' @export
#' @import LaplacesDemon mvtnorm

dmoe<-function(Data, intervals, particle_size, exp_col, mix_col=NULL, n_comp=1,
               v_init = 1, time_ind=1, response_ind=2, alpha=0.45,
               duplicate_num=2, proposal_method="linearbayes", return_all=F){

  interval_length<-length(intervals)-1 # number of time periods
  pb <- txtProgressBar(min=0,max=interval_length,initial=0,char="_",style = 3)
  # dimension of the design matrix
  design_dim <- n_comp*(length(exp_col)+1) + (n_comp-1)*(length(mix_col)+1)
  particles_sample<-list() # save sampled particles
  ess_persec <- c() # effective sample size per second
  log_pred_score <-0 # log predictive score
  log_pred_score_list <- c() # list of all batch log predictive scores
  # Initialization step

  importance_weights<-rep(1/particle_size, particle_size) # initialize weights
  # sample particles (regression coefficients) at time point j=0
  particles<-rmvnorm(particle_size,rep(0,design_dim), v_init*diag(design_dim))
  allWeights<-list()

  # propagate particles at time point j>1
  for(j in 2:(interval_length+1)){
    start_time <- Sys.time()
    D <-Data[which(Data[,time_ind] <= intervals[j] &
                     Data[,time_ind] > intervals[j-1]),] #interval data
    y <- D[,response_ind]
    x <- as.matrix(D[,c(exp_col)]) # covariate in the component models

    if (n_comp > 1) z <- as.matrix(D[,c(mix_col)]) # covariate in the mixture weight model

    particles_update <- online_update(y, x, n_comp, particle_size, particles,
                                         importance_weights, alpha, z)

    newimportance_weights<-normalize(particles_update$weights)
    newsample <- reject_strat_resample(newimportance_weights,
                                       particle_size/duplicate_num)
    index <-newsample$index # indexes of sampled particles
    importance_weights<-newsample$weights

    # compute the log predictive score
    if(j > floor(interval_length/2) + 1){
      log_pred_score <- log_pred_score + particles_update$log_pred_dens
      log_pred_score_list <- rbind(log_pred_score_list,
                                   particles_update$log_pred_dens)
    }
    # reset the previous parameter vector
    particles <- particles_update$particles[index,]

    # save the sampled particles
    if(return_all){
      particles_sample[[j-1]] <- particles
      allWeights[[j-1]] <- importance_weights
    }
    # compute the effective sample size per second
    end_time <- Sys.time()
    run_time <- end_time - start_time
    run_time <- ifelse(attr(run_time, "units") == "secs", as.numeric(run_time),
                       ifelse(attr(run_time, "units") == "mins",
                              as.numeric(run_time) * 60,
                              as.numeric(run_time) * 3600))
    # compute the effective sample size per second
    ess_persec[j-1]<-(1/sum(importance_weights^2))/run_time
    setTxtProgressBar(pb, value=j)
  }

  if(return_all){
    return(list(particles=particles_sample, ess=ess_persec,lps=log_pred_score,
                weights=allWeights, lps_list=log_pred_score_list,
                model=list(discount=alpha, components=n_comp,intervals=intervals,
                           exp_col= exp_col, mix_col=mix_col,
                           response=response_ind, proposal = proposal_method)))
  }else{
    return(list(particles=particles,ess=ess_persec,lps=log_pred_score,
                weights=importance_weights, lps_list=log_pred_score_list,
           model=list(discount=alpha, components=n_comp,intervals=intervals,
                      exp_col=exp_col, mix_col=mix_col, response=response_ind,
                      proposal = proposal_method)))
  }
}
