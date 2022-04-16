rm(list=ls())
graphics.off()
cat("\014")

# load simulation study file
source('model_selection_simulation_study.R')

# ------- run the simulation study ------
static_poisson_dgp<- simulation_study_wrapper(dgp = "staticpois", iterations=50)
save(static_poisson_dgp, file="results/static_poisson_dgp_results.R")

