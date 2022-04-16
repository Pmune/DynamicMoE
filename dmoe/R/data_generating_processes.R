
#' Simulate data from single component dynamic Poisson data generating process
#' @export
dynamic_poisson<-function(parameter,dim.design,dates,n){
  log.mu<-c()
  sim.data<-c()
  for(m in 1:length(dates)){
    x.train<-matrix(runif(n*dim.design,-1,1),nrow=n,ncol=dim.design,byrow=T)
    log.mu<-cbind(1,x.train)%*%(parameter[m,])
    y<-c()
    for(n in 1:nrow(x.train)){
      y[n]<-rpois(1,exp(log.mu[n]))
    }
    sim.data<-rbind(sim.data,data.frame(Date=rep(dates[m],length(y)),
                                        y=y,x=x.train))
  }
  return(sim.data)
}


#' Simulate data from dynamic Poisson mixture of experts data generating process
#' @export
dynamic_poisson_mixture<-function(beta,theta,dim.design,dates,n){
  sim.data<-c()
  for(m in 1:length(dates)){
    y<-c()
    x.train<-matrix(runif(n*dim.design,-1,1),nrow=n,
                    ncol=dim.design,byrow=T)
    lambda<-exp(cbind(1,x.train)%*%beta[m,,])
    psi<-cbind(1,x.train)%*%(theta[m,])
    for(n in 1:nrow(psi)){
      prob<-c(1,psi[n,])

      ind<-rcat(1,normalize(prob))

      y[n]<-rpois(1,lambda[n,(ind)])

    }

    sim.data<-rbind(sim.data,data.frame(Date=rep(dates[m],length(y)),
                                        y=y,x=x.train))
  }
  return(sim.data)
}

#' Simulate data from static Poisson data generating process
#' @export
static_poisson<-function(parameter,dim.design,dates,n){
  log.mu<-c()
  sim.data<-c()
  for(m in 1:length(dates)){
    x.train<-matrix(runif(n*dim.design,-1,1),nrow=n,ncol=dim.design,byrow=T)
    log.mu<-cbind(1,x.train)%*%t(parameter)
    y<-c()
    for(n in 1:nrow(x.train)){
      y[n]<-rpois(1,exp(log.mu[n]))
    }
    sim.data<-rbind(sim.data,data.frame(Date=rep(dates[m],length(y)),
                                        y=y,x=x.train))
  }
  return(sim.data)
}

