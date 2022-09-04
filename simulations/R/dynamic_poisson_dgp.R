rm(list=ls())
graphics.off()
cat("\014")

# load simulation study file
source('R/model_selection_simulation_study.R')

# ------- run the simulation study ------
start_time = Sys.time()
dynamic_poisson_dgp<- simulation_study_wrapper(dgp = "dynamicpois",iterations=5)
end_time = Sys.time()
print(end_time - start_time)
save(dynamic_poisson_dgp,file=paste0("results/dynamic_poisson_dgp_results_", Sys.Date(), ".R"))
