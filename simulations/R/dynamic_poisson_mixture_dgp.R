rm(list=ls())
graphics.off()
cat("\014")

# load simulation study file
source('R/model_selection_simulation_study.R')

# ------- run the simulation study ------

start_time = Sys.time()
dynamic_poissonmix_dgp<- simulation_study_wrapper(dgp = "dynamicpoismix", iterations=5)
end_time = Sys.time()
print(end_time - start_time)
save(dynamic_poissonmix_dgp, file=paste0("results/dynamic_poissonmix_dgp_results_", Sys.Date(), ".R"))
