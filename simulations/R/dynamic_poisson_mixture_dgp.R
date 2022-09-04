# run the simulation experiment with the dynamic mixture of Poisson data generating process.
# load simulation study file
source('R/model_selection_simulation_study.R')

# ------- run the simulation study ------

message("----------Experiment with dynamic mixture of Poisson DGP starts---------")
start_time = Sys.time()
dynamic_poissonmix_dgp<- simulation_study_wrapper(dgp = "dynamicpoismix", iterations=10)
end_time = Sys.time()
print(end_time - start_time)

# files are saved with a tag appended to the file name to avoid overwriting files.
# if you have change the tag please make sure to change the file_tag in the plot files.
file_tag <- "1" # file tag

save(dynamic_poissonmix_dgp, file=paste0("results/dynamic_poissonmix_dgp_results_", file_tag, ".R"))
message("---------Experiment with dynamic mixture of Poisson DGP ends------------")
