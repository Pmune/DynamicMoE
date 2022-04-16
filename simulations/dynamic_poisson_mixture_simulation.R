rm(list=ls())
graphics.off()
cat("\014")

# load simulation study file
source('model_selection_simulation_study.R')

# ------- run the simulation study ------
dynamic_poissonmix_dgp<- simulation_study_wrapper(dgp = "dynamicpoismix", iterations=50)
save(dynamic_poissonmix_dgp,file="results/dynamic_poissonmix_dgp_results.R")

