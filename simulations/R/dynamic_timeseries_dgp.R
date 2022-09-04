# Training and computing the log predictive score (lps) for
# the autoregressive timeseries model with static parameters (M5) in section 5.
rm(list=ls())
graphics.off()
cat("\014")

source("R/train_autoregressive_timeseries_models.R")

data_path <- "data/simdata_dyn.csv"
sim_data <- read.csv(data_path, header = FALSE)
discount_grid <- c(0.1,0.2,0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99)

# In the paper we use batch_size = 10, 25, and 50.
# Please change the batch_size value bellow accordingly.
batch_size <- 10 # batch size used to partition the data.
n_particles <- 1000
start_time <- Sys.time()
# run the simulation
simul_results <- train_autoregressive_models(sim_data, discount_grid, batch_size, n_particles)
end_time <- Sys.time()
dynamic_dgp_batch_lps <- simul_results$lps_data
dynamic_dgp_lps_summary <- simul_results$lps
dynamic_dgp_optimal_discounts <- simul_results$optimal_discount
print(end_time - start_time)
# Save the results
write.csv(
    dynamic_dgp_batch_lps,
    paste0("results/dynamic_dgp_batch_lps_batchsize", batch_size,"_",Sys.Date(),".csv"),
    row.names=FALSE
    )

write.csv(
    dynamic_dgp_lps_summary,
    paste0("results/dynamic_dgp_lps_summary_batchsize", batch_size,"_",Sys.Date(),".csv"),
    row.names=FALSE
    )

write.csv(
    dynamic_dgp_optimal_discounts,
    paste0("results/dynamic_dgp_optimal_discounts_batchsize",batch_size,"_",Sys.Date(),".csv"),
    row.names=FALSE
    )
