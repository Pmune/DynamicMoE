rm(list=ls())
graphics.off()
cat("\014")
setwd("C:/Users/parfait/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/codes/training_models/")
setwd("C:/Users/eparfmu/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/codes/training_models/")

# required paths

mixture_model_path <- 'C:/Users/parfait/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/codes/dynamic_mixture/'
mixture_model_path <- 'C:/Users/eparfmu/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/codes/dynamic_mixture/'
poisson_model_path <- 'C:/Users/parfait/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/codes/dynamic_poisson/'
poisson_model_path <- 'C:/Users/eparfmu/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/codes/dynamic_poisson/'
common_file_path <- 'C:/Users/parfait/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/codes/common/'
common_file_path <- 'C:/Users/eparfmu/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/codes/common/'
data_file_path <- 'C:/Users/parfait/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/codes/training_models/'
data_file_path <- 'C:/Users/eparfmu/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/codes/training_models/'

# dynamic mixture model files
source(paste0(mixture_model_path,"/dynamic_mixture_particle_filter.R"))
source(paste0(mixture_model_path,"/dynamic_mixture_proposal.R"))
source(paste0(mixture_model_path,"/dynamic_mixture_prediction.R"))

# files for  poisson models
source(paste0(poisson_model_path,"/dynamic_poisson_particle_filter.R"))
source(paste0(poisson_model_path,"/dynamic_linear_bayes_poisson_proposal.R"))
source(paste0(poisson_model_path,"/dynamic_poisson_prediction.R"))

# common source files
source(paste0(common_file_path,"/stratified_rejection_resampling.R"))

# load the data 
source(paste0(data_file_path,"/data_cleaner.R"))

# ------set the model parameters  ---------

set.seed(600)
y.max<-11
mix.col<-c(3)
L<-length(time.intervals)-1 # length of the time intervals
N.particles<-15000 # number of particles (recommended number of particles: 10000)
R<-2
max_comp <- 4 #maximun nubr of components
exp.col<-list(Model1=c(9),Model2=c(4,9),Model3=c(4,9,3),Model4=c(4,9,3,8),Model5=c(4,9,3,5,8),Model6=c(4,9,3,5,8,7)) # columns used to estimate the mean #TR
LPS.dynamic<-matrix(0,length(exp.col),max_comp)
LPS.static<-matrix(0,length(exp.col),max_comp)
system.time(
  for(m in 1:length(exp.col)){
    set.seed(600)
    output<-FAPF1Comp(Data,time.intervals,M=N.particles,c(exp.col[[m]]),
                      time_ind=1,response_ind=2,alpha=.99,R,h=1) # run the static model
    LPS.static[m,1]<-output$LPS
    output<-FAPF1Comp(Data,time.intervals,M=N.particles,c(exp.col[[m]]),
                      time_ind=1,response_ind=2,alpha=.5,R,h=1)# run the dynamic model
    LPS.dynamic[m,1]<-output$LPS
    # 2. more than two components
    for(K in 2:max_comp){
      print("running static model")
      output<-FAPF(Data,time.intervals,N.particles,mix.col,c(exp.col[[m]]),
                   n_comp=K,time_ind=1,response_ind=2,alpha=.99,R,F,h=1)
      LPS.static[m,K]<-output$LPS
      
      print("running dynamic model")
      output<-FAPF(Data,time.intervals,N.particles,mix.col,c(exp.col[[m]]),
                   n_comp=K,time_ind=1,response_ind=2,alpha=.5,R,F,h=1)
      LPS.dynamic[m,K]<-output$LPS
    }
    print(LPS.dynamic)
    print(LPS.static)
    
  }
)
model_comparison_lps <- list(dynamic=LPS.dynamic, static=LPS.static)
save(model_comparison_lps,file="model_comparison_lps_new_2.R")
