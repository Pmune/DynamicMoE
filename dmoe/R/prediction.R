
dmixpois<-function(y,lambda,Mix.wgt){

  dens<-diag(dpois(y,lambda)%*%t(Mix.wgt))
  return(dens)
}


#' Predictive density for dynamic mixture of experts.
#' @export

pred_density<-function(model, test_data, y_max=NULL, batch_ind=NULL){
  pred_dens_vec <- c()
  n_comp <- model$model$components
  y_max <- ifelse(is.null(y_max), max(test_data[, model$model$response]), y_max)
  y_grid <- c(0 : y_max) # grid on the response values
  # selecting the posterior to use: is batch_ind is specified the posterior for the batch is used
  # else the last batch posterior is used.

  if (!is.null(batch_ind)){
    importance_weight <- model$weights[[batch_ind]]
    particles <- model$particles[[batch_ind]]
  }else{
    # in case the batch_ind is not provided, use the last batch's posterior
    num_batches <- length(model$model$intervals[-1])
    print(num_batches)
    # if the posterior for all batches were saved in the model select the last one.
    if(length(model$weights) == num_batches){
      importance_weight <- model$weights[[num_batches]]
      particles <- model$particles[[num_batches]]
    }else{
      # if all posteriors were not save in the model
      importance_weight <- model$weights
      particles <- model$particles
    }
  }

  # predictions
  for(i in seq_len(nrow(test_data))){
    x <- as.numeric(test_data[i, model$model$exp_col])
    if(n_comp > 1){
      z <- as.numeric(test_data[i, model$model$mix_col])
      x_mat <- design_matrix(x,z,n_comp)
      eta <- particles %*% x_mat
      lambda <- exp(eta[,seq_len(n_comp)]) # expected effect
      psi<-cbind(1,exp(eta[,-seq_len(n_comp)])) # unnormalized mixture weights
      mix_wgt<-psi/rowSums(psi) # mixture weights
      # compute the predictive density
      pred_dens<-apply(as.matrix(y_grid),1,function(y)sum(dmixpois(y,lambda, mix_wgt)*importance_weight))
    } else{
      x_mat <- c(1, x)
      eta <- particles %*% x_mat
      lambda <- exp(eta[,seq_len(n_comp)]) # expected effect
      pred_dens <- apply(as.matrix(y_grid),1,function(y)sum(dpois(y,lambda)*importance_weight))
    }
    pred_dens <- pred_dens/sum(pred_dens) # normalize the density
    pred_dens_vec<-rbind(pred_dens_vec, pred_dens)
  }

  return(data.frame(y = y_grid, dens = colMeans(pred_dens_vec)))
}
