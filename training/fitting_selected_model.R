rm(list=ls())
graphics.off()
cat("\014")

library(dmoe)

# load the data
source("data_cleaner.R")

# ------set the model parameters  ---------
exp_col<-c("ChangedModule","FileComplexity")
mix_col<-c("Commits")

n_particles<-10000 #number of particles (recommended number of particles: 10000)
alpha<-0.5
n_comp <-2
density.fitted<-c()
output_pois <- dmoe(Data, time.intervals, n_particles, mix_col, exp_col,
               n_comp=1, alpha=alpha, return_all=TRUE)

message(paste0("ESS for Poisson model: ",mean(output_pois$ess)))
message(paste0("LPS for Poisson model: ",output_pois$lps))

# dynamic model with 2 components
output_mix <- dmoe(Data, time.intervals, n_particles, mix_col, exp_col,
                    n_comp=n_comp, alpha=alpha, return_all=TRUE)

message(paste0("ESS for 2 components Poisson mixture model: ",mean(output_mix$ess)))
message(paste0("LPS for 2 components Poisson mixture model: ",output_mix$lps))

# save output
model_posterior <- list(poisson=output_pois, poisson_mix=output_mix)
save(model_posterior,file="outputs/selected_model_posterior.R")

