rm(list=ls())
graphics.off()
cat("\014")

# load simulation study file
source('model_selection_simulation_study.R')

# ------- run the simulation study ------
dynamic_poisson_dgp<- simulation_study_wrapper(dgp = "dynamicpois",iterations=50)
save(dynamic_poisson_dgp,file="results/dynamic_poisson_dgp_results.R")
