# Training and computing the log predictive score (lps) for
# the autoregressive timeseries model with static parameters (M5) in section 5.
rm(list=ls())
graphics.off()
cat("\014")

source("R/train_autoregressive_timeseries_models.R")

data_path <- "data/simdata_dyn.csv"
sim_data <- read.csv(data_path, header = FALSE)[,1:5]
discount_grid <- c(0.1,0.2,0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99)

# files are saved with a tag appended to the file name to avoid overwriting files.
# if you have change the tag please make sure to change the file_tag in the plot files.
file_tag <- "1" # file tag

# Run the simulation study. In the paper we use batch_size = 10, 25, and 50.
start_time <- Sys.time()
for (batch_size in c(10, 25, 50)){
simul_results <- train_autoregressive_models(sim_data, discount_grid, batch_size, n_particles)

dynamic_dgp_batch_lps <- simul_results$lps_data
dynamic_dgp_lps_summary <- simul_results$lps
dynamic_dgp_optimal_discounts <- simul_results$optimal_discount

# Save the results
write.csv(
    dynamic_dgp_batch_lps,
    paste0("results/dynamic_dgp_batch_lps_batchsize", batch_size,"_", file_tag,".csv"),
    row.names=FALSE
    )

write.csv(
    dynamic_dgp_lps_summary,
    paste0("results/dynamic_dgp_lps_summary_batchsize", batch_size,"_", file_tag,".csv"),
    row.names=FALSE
    )

write.csv(
    dynamic_dgp_optimal_discounts,
    paste0("results/dynamic_dgp_optimal_discounts_batchsize",batch_size,"_",file_tag,".csv"),
    row.names=FALSE
    )
}
print(end_time - start_time)
end_time <- Sys.time()
