rm(list=ls())
graphics.off()
cat("\014")
setwd("C:/Users/parfait/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/codes/training_models/")
setwd("C:/Users/eparfmu/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/codes/training_models/")

# ----- required paths -----------------------------

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

set.seed(500)
y.max<-8
exp.col<-c("ChangedModule","FileComplexity")
mix.col<-c("Commits")
L<-length(time.intervals)-1 # length of the time intervals
N.particles<-1000 # number of particles (recommended number of particles: 10000)
alpha<-0.5
R<-2
output_efficiency <- list()
# dynamic model with 2 components
for (iter in seq_len(100)){
system.time(output<-FAPF(Data, time.intervals, N.particles, mix.col, exp.col,
                         n_comp=2, time_ind=1, response_ind=2, alpha, R, T, h=1))

  output_efficiency [[iter]]<- output

}


save(output_efficiency,file="posterior_for_model_efficiency_new_1.R")

