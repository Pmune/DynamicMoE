rm(list=ls())
graphics.off()
cat("\014")

# load simulation study file
source('R/model_selection_simulation_study.R')

# ------- run the simulation study ------
start_time = Sys.time()
static_poisson_dgp<- simulation_study_wrapper(dgp = "staticpois", iterations=5)
end_time = Sys.time()
print(end_time - start_time)
save(static_poisson_dgp, file=paste0("results/static_poisson_dgp_results_", Sys.Date(), ".R"))
