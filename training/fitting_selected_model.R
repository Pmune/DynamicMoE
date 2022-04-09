rm(list=ls())
graphics.off()
cat("\014")

library(dmoe)

# load the data
source("data_cleaner.R")

# ------set the model parameters  ---------

set.seed(600)
exp.col<-c("ChangedModule","FileComplexity")
mix.col<-c("Commits")

N.particles<-100000 #number of particles (recommended number of particles: 10000)
alpha<-0.5
R<-2
n_comp=2
density.fitted<-c()
system.time(output_pois<-FAPF1Comp(Data,time.intervals,M=N.particles,exp.col,
                                   time_ind=1,response_ind=2,alpha,R,h=1)) # run the model with one component
message(paste0("ESS for Poisson model: ( ",mean(output_pois$ESS), " )"))
message(paste0("LPS for Poisson model: ( ",output_pois$LPS, " )"))

# dynamic model with 2 components
system.time(output_mix<-FAPF(Data,time.intervals,N.particles,mix.col,exp.col,
                             n_comp,time_ind=1,response_ind=2,alpha,R,F,h=1))

message(paste0("ESS for 2 components Poisson mixture model: ( ",mean(output_mix$ESS), " )"))
message(paste0("LPS for 2 components Poisson mixture model: ( ",output_mix$LPS, " )"))

# save output
model_posterior <- list(poisson=output_pois, poisson_mix=output_mix)
save(model_posterior,file="outputs/selected_model_posterior.R")

