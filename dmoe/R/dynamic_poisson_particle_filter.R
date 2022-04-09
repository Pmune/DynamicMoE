# likelihood for dynamic poisson model
likelihood<-function(y,X,par)
{
  log.likelihood<-0
  for(obs in 1:length(y))
  {
    lambda<-exp(X[obs,]%*%par) # expected effect
    log.likelihood<-log.likelihood+dpois(y[obs],lambda,log=TRUE)
  }
  return(log.likelihood)
}

# function for computing the importance weights
computeWeight.Pois<-function(yj,Xj,weights,param,Prevparticles,
                             W,proposalMean,proposalCov){
  dim.par<-dim(W)[2]
  dens.proposal<-dmvnorm(param,proposalMean,proposalCov,log=TRUE)
  dens.prior<-dmvnorm(Prevparticles,param,W+diag(1E-7,dim.par),log=TRUE)
  p.prior<-log(sum(exp(dens.prior)*weights))
  wt<-likelihood(yj,Xj,param)+ p.prior-dens.proposal
  return(wt)
}

# particle filter for the dynamic poisson model
FAPF1Comp<-function(Data,intervals,M,train.col,time_ind=1,response_ind=2,alpha=0.45,R=5,h=1,returnAllweights=F){

  interval.length<-length(intervals)-1 # number of time periods
  pb <- txtProgressBar(min=0,max=interval.length,initial=0,char="_",style = 3)
  # set.seed(600) # set seeds

  # check if required package is installed
  #if (!require("pacman")) install.packages("pacman",quietly=TRUE)
  pacman::p_load(mvtnorm, MCMCpack,LaplacesDemon) # load required packages

  # Variable declaration and initialization

  design.dim<-length(train.col)+1# dimension of the design matrix
  Particles<-list() # save sampled particles
  weight.state<-rep(1/M,M) # initialize weights
  covMat<-1*diag(design.dim)# set the initial covariance
  Theta.prev<-rmvnorm(M,rep(0,design.dim),covMat) # sample theta at time point j=0
  Mean.Theta<-colMeans(Theta.prev)
  allWeights<-list()
  # Designing the initial  distribution.
  Dj<-Data[which(Data[,time_ind]<=intervals[2]&Data[,time_ind]>intervals[1]),] #interval data
  yj<-Dj[,response_ind]
  # effect covariates
  Xj<-as.matrix(cbind(rep(1,nrow(Dj)),Dj[,c(train.col)])) # design matrix
  PropMeanCov<-single_component_proposal(yj,Xj,colMeans(Theta.prev),covMat)
  Theta<-mvrnorm(M,PropMeanCov$Means,PropMeanCov$cov)

  #update the mean vectors and covariance matrices using data at t=1
  Init.weights<-apply(as.matrix(1:M),1,function(obs)likelihood(yj,Xj,Theta[obs,]))

  Init.weights<-normalize(Init.weights)

  # indexes of the resampled particles
  NewSample<-reject_strat_resample(Init.weights,M/R)
  index<-NewSample$index # indexes of sampled particles
  NewState.weights<-NewSample$weights
  covMat<-cov(Theta[index,])/alpha
  Mean.Theta<-colSums(Theta[index,]*NewState.weights)
  Prev.theta<-Theta[index,] # resampled particles

  # propagate paricles at time point t>1
  ESS<-c() # effective sample size
  N.unique<-c()
  score<-0
  for(t in 2:(interval.length+1)){

    # generating risk set.
    Dj<-Data[which(Data[,time_ind]<=intervals[t]&Data[,time_ind]>intervals[t-1]),] #interval data
    #Dj<-Dj[order(Dj[,response_ind], decreasing = TRUE),]
    yj<-Dj[,response_ind]
    # effect covariates
    Xj<-as.matrix(cbind(rep(1,nrow(Dj)),Dj[,c(train.col)])) # design matrix
    # propose particle at time t
    PropMeanCov<-single_component_proposal(yj,Xj, Mean.Theta,covMat)
    Theta<-mvrnorm(M,PropMeanCov$Means,PropMeanCov$cov*h)

    # compute weights at time j

    weight.state<-apply(as.matrix(1:M),1,function(j)
      computeWeight.Pois(yj,Xj,NewState.weights,Theta[j,],Prev.theta,covMat*(1-alpha),
                         PropMeanCov$Means,PropMeanCov$cov*h))#denominator

    weight.state[which(is.na(weight.state))]<--log(1e-300)
    weight.state[is.infinite(weight.state)]<-log(1e-300)

    # normalize weights
    weight.state<-normalize(weight.state)
    NewSample<-reject_strat_resample(weight.state,M/R)
    index<-NewSample$index # indexes of sampled particles
    NewState.weights<-NewSample$weights
    Prev.theta<-Particles[[t-1]]<-Theta[index,]
    covMat<-cov(Theta[index,])/alpha
    Mean.Theta<-colSums(Theta[index,]*NewState.weights)
    allWeights[[t-1]]<-NewState.weights
    #N.unique[t-1]<-length(unique(NewState.weights))
    N.unique[t-1]<-length(unique(index))
    ESS[t-1]<-1/sum(weight.state^2)
    if(t>floor((interval.length-1)/2)){
      logL<-apply(as.matrix(1:dim(Particles[[t-2]])[1]),1,function(m)likelihood(yj,Xj,Particles[[t-2]][m,]))
      score<- score+log(sum(as.numeric(exp(logL))*allWeights[[t-2]]))

    }
    setTxtProgressBar(pb, value=t)
  }
  if(returnAllweights){
    return(list(Particles=Particles,ESS=ESS,LPS=2*score/(interval.length-1),weights=allWeights,N.unique=N.unique,CovMat=PropMeanCov$cov*alpha))
  }else{
    return(list(Particles=Particles,ESS=ESS,LPS=2*score/(interval.length-1),weights=NewState.weights,N.unique=N.unique,CovMat=PropMeanCov$cov*alpha))

  }
}

