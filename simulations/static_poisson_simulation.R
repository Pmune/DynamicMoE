rm(list=ls())
graphics.off()
cat("\014")


library(dmoe)
# simulation study file
source('simulation_study_file.R')

# -------------------function for simulating data------------------

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

# ------- run the simulation ------

results<- simulation_study("static_poisson",static_poisson, max_iter=50)
save(results,file="results/static_poison_results.R")

