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

# files are saved with a tag appended to the file name to avoid overwriting files.
# if you have change the tag please make sure to change the file_tag in the plot files.
file_tag <- "1" # file tag

save(dynamic_poisson_dgp,file=paste0("results/dynamic_poisson_dgp_results_", file_tag, ".R"))
