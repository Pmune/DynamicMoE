# Training and computing the log predictive score (lps) for
# the autoregressive timeseries model with static parameters (M4) in section 5.
rm(list=ls())
graphics.off()
cat("\014")

source("R/train_autoregressive_timeseries_models.R")

data_path <- "data/simdata_static.csv" # 50 datasets generated from M4 saved in the csv file.
sim_data <- read.csv(data_path, header = FALSE)
discount_grid <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8, 0.9, 0.99)

# In the paper we use batch_size = 10, 25, and 50.
# Please change the batch_size value bellow accordingly.
batch_size <- 10 # batch size used to partition the data.
start_time <- Sys.time()
# run the simulation
simul_results <- train_autoregressive_models(sim_data, discount_grid, batch_size)
end_time <- Sys.time()
static_dgp_batch_lps = simul_results$lps_data
static_dgp_lps_summary <- simul_results$lps
static_dgp_optimal_discounts <- simul_results$optimal_discount
print(end_time - start_time)
# Save the results
write.csv(
    static_dgp_batch_lps,
    paste0("results/static_dgp_batch_lps_batchsize",batch_size,"_", Sys.Date(), ".csv"),
    row.names=FALSE
    )
write.csv(
    static_dgp_lps_summary,
    paste0("results/static_dgp_lps_summary_batchsize", batch_size,"_", Sys.Date(), ".csv"),
    row.names=FALSE
    )
write.csv(
    static_dgp_optimal_discounts,
    paste0("results/static_dgp_optimal_discounts_batchsize",batch_size,"_",Sys.Date(),".csv"),
    row.names=FALSE
    )
