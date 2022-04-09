MPFpredictions1Comp<-function(Test.data,Particles,ImportanceWeight,y.max,exp.col){
  
  d.ftd<-c()
  
  for(m in 1:nrow(Test.data)){
    X<-cbind(1,Test.data[m,exp.col])
    lambda<-exp(Particles%*%t(X))
    # compute the predictive density
    pmf<-apply(as.matrix(c(0:y.max)),1,function(y)sum(dpois(y,lambda)*ImportanceWeight))
    pmf<- pmf/sum(pmf) # normalize the density
    d.ftd<-rbind(d.ftd,pmf)
    
  }
  
  d.ftd<-colMeans(d.ftd)
  return(d.ftd)
}
