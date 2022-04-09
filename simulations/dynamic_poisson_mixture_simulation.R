rm(list=ls())
graphics.off()
cat("\014")


library(dmoe)
source('simulation_study_file.R') # simulation study file

# -------------------function for simulating data------------------
dynamic_poisson_mixture<-function(beta,theta,dim.design,dates,n){
  sim.data<-c()
  for(m in 1:length(dates)){
    y<-c()
    x.train<-matrix(runif(n*dim.design,-1,1),nrow=n,ncol=dim.design,byrow=T)
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

# ------- run the simulation ------

results<- simulation_study("dynamic_poisson_mixture", dynamic_poisson_mixture,max_iter=50)
save(results,file="results/dynamic_poisson_mixture_results.R")

