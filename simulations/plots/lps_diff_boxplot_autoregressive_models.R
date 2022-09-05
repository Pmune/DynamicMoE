# plot lps difference of the dynamic models vs static models fitted
# to the data generating processes M4 and M5.
# M4 - the autoregressive timeseries model with static parameters and
# M5 - autoregressive timeseries model with dynamic parameter.
# Figure 5.3 in the paper.
rm(list=ls())
graphics.off()
cat("\014")

# Libraries
library("tidyverse")
library("hrbrthemes")
library("viridis")
require(gridExtra)
hrbrthemes::import_roboto_condensed()

# Load dynamic MoE model LPS

file_tag <- "1"
dmoe_dyn_10 = read.table(paste0("results/dynamic_dgp_batch_lps_batchsize10_", file_tag, ".csv"),
                       sep = ",", header = TRUE)
dmoe_dyn_25 = read.table(paste0("results/dynamic_dgp_batch_lps_batchsize25_", file_tag, ".csv"),
                         sep = ",", header = TRUE)
dmoe_dyn_50 = read.table(paste0("results/dynamic_dgp_batch_lps_batchsize50_", file_tag, ".csv"),
                         sep = ",", header = TRUE)
dmoe_static_10 = read.table(paste0("results/static_dgp_batch_lps_batchsize10_",file_tag, ".csv"),
                            sep = ",", header = TRUE)
dmoe_static_25 = read.table(paste0("results/static_dgp_batch_lps_batchsize25_", file_tag, ".csv"),
                            sep = ",", header = TRUE)
dmoe_static_50 = read.table(paste0("results/static_dgp_batch_lps_batchsize50_", file_tag, ".csv"),
                            sep = ",", header = TRUE)

# Load MoE model LPS
moe_dyn_10 = read.table("results/LPS_dyn_batchsize_10.csv",
                        sep = ",", header = FALSE)
moe_dyn_25 = read.table("results/LPS_dyn_batchsize_25.csv",
                        sep = ",", header = FALSE)
moe_dyn_50 = read.table("results/LPS_dyn_batchsize_50.csv",
                        sep = ",", header = FALSE)
moe_static_10 = read.table("results/LPS_static_batchsize_10.csv",
                           sep = ",", header = FALSE)
moe_static_25 = read.table("results/LPS_static_batchsize_25.csv",
                           sep = ",", header = FALSE)
moe_static_50 = read.table("results/LPS_static_batchsize_50.csv",
                           sep = ",", header = FALSE)

# create lps diff datasets
batch10 <- data.frame(
  batch =c( rep("Batch size = 10",100)),
  dgp=c( rep("M4",50), rep("M5",50) ),
  lps=c( colSums(dmoe_static_10) - colSums(moe_static_10) ,
         colSums(dmoe_dyn_10) - colSums(moe_dyn_10 ))
)

batch25 <- data.frame(
  batch =c( rep("Batch size = 25",100)),
  dgp=c( rep("M4",50), rep("M5",50) ),
  lps=c( colSums(dmoe_static_25) - colSums(moe_static_25) ,
         colSums(dmoe_dyn_25) - colSums(moe_dyn_25 ) )
)

batch50 <- data.frame(
  batch =c( rep("Batch size = 50",100)),
  dgp = c( rep("M4",50), rep("M5",50) ),
  lps = c( colSums(dmoe_static_50) - colSums(moe_static_50) ,
           colSums(dmoe_dyn_50) - colSums(moe_dyn_50) )
)

lps_data <- rbind(batch10,batch25,batch50 )

# Plot
lps_boxplot <- lps_data %>%
  ggplot(aes(x=reorder(dgp,lps), y=lps, fill=dgp)) +
  geom_boxplot() +
  facet_wrap(~batch ,scales="free_x") +
  scale_fill_viridis(discrete = TRUE, alpha=0.4) +
  theme_bw()+
  theme(
    legend.position="none",
    plot.title = element_text(size=11, hjust = 0.5),
    strip.background = element_rect(fill = "white", linetype = 0))+
  xlab(" ") +
  scale_y_continuous(breaks=seq(-100, 400,100)) +
  ylab("LPS diff")

print(lps_boxplot)
png("plots/lps_diff_boxplot_autoregressive_models.png", width=600, height=300)
print(lps_boxplot)
dev.off()
