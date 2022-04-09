library(hesim)
simulation_study <- function(data_generating,simulation_function, max_iter){
  # initializing the simulation
  dates<-seq(as.Date("2017-01-01"),as.Date("2017-12-31"),length.out=12)
  interval<-dates
  interval<-c(interval[1]-1,interval)
  mix.col<-c(3) # column used to define clusters
  exp.col<-c(3) # columns used to estimate the mean #TR
  L<-length(dates) # length of the time intervals
  N.particles<-10000 # number of particles (recommended number of particles: 10000)
  R=2
  h<-1
  n=100
  dim.design=1
  discounts<-c(.4,.5,.6,.7,.8,.9,.95,.99)
  components<-c(1,2,3)
  LPS.Meanmat<-matrix(0,nrow=length(discounts),ncol=length(components))
  row.names(LPS.Meanmat)<-discounts
  colnames(LPS.Meanmat)<-components
  ModelFrequencyMat<-matrix(0,nrow=length(discounts),ncol=length(components))
  row.names(ModelFrequencyMat)<-discounts
  colnames(ModelFrequencyMat)<-components
  # simulate the parameter
  if(data_generating=="dynamic_poisson"){
    parameter<-cbind(-0.11+cumsum(rnorm(2*L,0,0.17)),2.29+cumsum(rnorm(2*L,0,0.2)))
    train_parameter <- parameter[1:L,]
    test_parameter <- parameter[(L+1):(2*L),]
  } else if(data_generating=="static_poisson"){
    parameter<-cbind(0.11,2.29)
    train_parameter <- parameter
    test_parameter <- parameter
  } else{
    beta1<-cbind(1.11+cumsum(rnorm(2*L,0,0.08)), 2.17+cumsum(rnorm(2*L,0,0.15)))
    beta2<-cbind(-0.08+cumsum(rnorm(2*L,0,0.07)),1.94+cumsum(rnorm(2*L,0,0.1)))
    beta<-array(NA,c(nrow(beta1),ncol(beta1),2))
    beta[,,1]<-beta1
    beta[,,2]<-beta2
    theta<-cbind(2.63+cumsum(rnorm(2*L,0,0.08)), -4.41+cumsum(rnorm(2*L,0,.17)))

  }
  KL.ratio<-c()
  LPS.dif<-c()
  system.time(for(iter in 1:max_iter){
    print(iter)
    LPS.mat<-matrix(0,nrow=length(discounts),ncol=length(components)) # object fot storing the LPS values
    if(data_generating!="dynamic_poisson_mixture"){
      sim.data<-simulation_function(train_parameter,dim.design,dates,n) # simulate dataset
      Test.data<-simulation_function(test_parameter,dim.design,dates,n) # simulate dataset
    }else{
      sim.data<-simulation_function(beta[1:L,,],theta[1:L,],dim.design,dates,n) # simulate dataset
      Test.data<-simulation_function(beta[(L+1):(2*L),,],theta[(L+1):(2*L),],dim.design,dates,n) # simulate dataset
    }
    # fit models with different discount factors and number of components
    for(j in components){ # loop over the components

      for(a in 1:length(discounts)){ # loop over the discount factors
        if(j==1){ # model with one component
          output<-FAPF1Comp(sim.data,interval,M=N.particles,exp.col,
                            time_ind=1,response_ind=2,alpha=discounts[a],R,h)
           #LPS.mat[a,j]<-output$LPS-j+1
          LPS.mat[a,j]<-ifelse(is.finite(output$LPS),output$LPS,log(1E-300))
        }
        else{ # models with more than one components
          output<-FAPF(sim.data,interval,N.particles,mix.col,exp.col,
                       n_comp =j,time_ind=1,response_ind=2,alpha=discounts[a],R,F,h)
          #LPS.mat[a,j]<-output$LPS-j+1
          LPS.mat[a,j]<-ifelse(is.finite(output$LPS),output$LPS,log(1E-300))

        }
      }
    }
    LPS.Meanmat<-LPS.Meanmat+LPS.mat # sum of LPS for different datasets
    #LPS.mat<-LPS.Meanmat/iter
    selectedModel<-as.vector(which(LPS.mat==max(LPS.mat), arr.ind = T)) # return the index of the best model
    ModelFrequencyMat[selectedModel[1],selectedModel[2]]<-ModelFrequencyMat[selectedModel[1],selectedModel[2]]+1 # increment the frequency of the selected model
    print(ModelFrequencyMat)

    # run the model selected based on the highest LPS
    if(selectedModel[2]==1){ # model with one component
      output<-FAPF1Comp(Test.data,interval,N.particles,exp.col,time_ind=1,
                        response_ind=2,alpha=discounts[selectedModel[1]],R,h)

      LPS.selected<-ifelse(is.finite(output$LPS),output$LPS,log(1E-300))
    }
    else{
      output<-FAPF(Test.data,interval,N.particles,mix.col,exp.col,
                   n_comp = selectedModel[2],time_ind=1,response_ind=2,
                   alpha=discounts[selectedModel[1]],R,F,h)
      LPS.selected<-ifelse(is.finite(output$LPS),output$LPS,log(1E-300))
    }
    # run the model selected based on the highest LPS and alpha=1
    if(selectedModel[2]==1){ # model with one component
      output<-FAPF1Comp(Test.data,interval,M=N.particles,exp.col,time_ind=1,
                        response_ind=2,alpha=.99,R,h)
      LPS.static<-ifelse(is.finite(output$LPS),output$LPS,log(1E-300))
    }
    else{
      output<-FAPF(Test.data,interval,N.particles,mix.col,exp.col,
                   n_comp=selectedModel[2],time_ind=1,response_ind=2,
                   alpha=.99,R,F,h)
      LPS.static<-ifelse(is.finite(output$LPS),output$LPS,log(1E-300))
    }
    LPS.dif[iter]<-LPS.selected-LPS.static
    print(LPS.dif)
  })

  return(list(freq=ModelFrequencyMat,LPS.dif=LPS.dif,LPS.Meanmat=LPS.Meanmat/max_iter))
}
