

library(dmoe)

interval_predictions <- function(Data, Test.data, time.intervals, post_output,
                                 pred_intervals, exp.col, mix.col, n_comp, y.max){
  observed_density <- c()
  density_fitted <- list()
  
  for (j in pred_intervals){
    if(j< (length(time.intervals)-1)){
      test.data <- Data[which(Data$Date <= time.intervals[j+2] & Data$Date > time.intervals[j+1]),]
    }else{
      test.data <- Test.data
    }
    observedFrequency<-apply(as.matrix(seq(0,y.max,1)),1,
                             function(x)length(which(test.data[,2]==x)))
    observed_density <- cbind(observed_density, observedFrequency/sum(observedFrequency))
    
    interval_density_fitted <- c()
    model_names <- names(post_output)
    for (m in seq_len(length(post_output))){
      
      particles<-post_output [[m]]$Particles[[j]]
      importanceWeight<-post_output [[m]]$weights[[j]]
      
      if (!model_names[m]== "poisson"){
      predictions<-MPFpredictions(test.data,particles,importanceWeight,
                                  y.max,exp.col,mix.col,n_comp)
      }else{
        predictions<-MPFpredictions1Comp(test.data, particles, importanceWeight, y.max, exp.col)
      }
      interval_density_fitted<-cbind(interval_density_fitted, predictions)
    }
    names(interval_density_fitted) = model_names
    density_fitted[[j]] <- interval_density_fitted
  }
  return(list(observed = observed_density, predicted = density_fitted))
}