library(dplyr)

# Obtain Results Dataset ---------------------------------------------------------
results_dataset <- read.csv("Transfer_Method/summary_GENERAL_AVG_Dataset.csv")

# Get mean and sd lists
mean_lists <- list(list(results_dataset$COR_BayI_Ave), list(results_dataset$NRMSE_BayI_Ave), 
                list(results_dataset$PM_Conv_Ave_10), list(results_dataset$PM_Conv_Ave_20), 
                list(results_dataset$COR_Trans_Ave), list(results_dataset$NRMSE_Trans_Ave),
                list(results_dataset$PM_Transf_Ave_10), list(results_dataset$PM_Transf_Ave_20)
              )
sd_lists <- list(list(results_dataset$COR_BayNI_SD), list(results_dataset$NRMSE_BayI_SD), 
              list(results_dataset$PM_Conv_SD_10), list(results_dataset$PM_Conv_SD_20), 
              list(results_dataset$COR_Trans_SD), list(results_dataset$NRMSE_Trans_SD), 
              list(results_dataset$PM_Transf_SD_10), list(results_dataset$PM_Transf_SD_20)
            )
# Declarate mean and sd values vectors
mean_values <- c()
sd_values <- c()

# Add values to mean and sd values vectors
for (i in 1:length(unique(results_dataset$Dataset))) {
  for (j in 1:length(mean_lists)) {
    mean_values <- c(mean_values, mean_lists[[j]][[1]][1:5])
    mean_lists[[j]][[1]] <- mean_lists[[j]][[1]][-c(1:5)]
    
    sd_values <- c(sd_values, sd_lists[[j]][[1]][1:5])
    sd_lists[[j]][[1]] <- sd_lists[[j]][[1]][-c(1:5)]
  }  
}

# Create results table
results_metrics_table <- data.frame(
  Dataset = rep(unique(results_dataset$Dataset), each = 40),
  Method = rep(c("GBLUP", "Transfer"), each = 20, times = length(unique(results_dataset$Dataset))),
  Metric = rep(c("COR", "NRMSE", "PM_10", "PM_20"), each = 5, times = length(unique(results_dataset$Dataset))),
  Testing_Proportion = rep(c(0.15, 0.30, 0.50, 0.75, 0.80), times = length(unique(results_dataset$Dataset))),
  Mean = mean_values,
  SD = sd_values
)

dim(results_metrics_table)
# Write on a csv file
write.csv(results_metrics_table, file.path("Transfer_Method", "appendix_datasets.csv"))



# Create average results table ---------------------------------------------------
results_avg_metrics_table <- results_metrics_table %>%
  group_by(Method, Metric, Testing_Proportion) %>%
  summarise(
    Mean = mean(Mean),
    SD = mean(SD),
    .groups = 'drop' # To drop the grouping structure after summarization
  )
# Write on a csv file
write.csv(results_avg_metrics_table, file.path("Transfer_Method", "appendix_average.csv"))

# Obtain each table basing on dataset --------------------------------------------
datasets_list <- c(
  "EYT_1", "EYT_2", "EYT_3", "Groundnut", "Indica", # NOT RESOLVED JAPONICA
  "Maize", "Wheat_1", "Wheat_2", "Wheat_3",
  "Wheat_4", "Wheat_5", "Wheat_6"
)

for (i in 1:length(datasets_list)) {
  results_dataset_metrics_table <- data.frame()
  
  dataset_name <- datasets_list[i]
  results_dataset_metrics_table <- results_metrics_table[results_metrics_table$Dataset == dataset_name, ][-1]
  row.names(results_dataset_metrics_table) <- NULL
  
  write.csv(results_dataset_metrics_table, file.path("Transfer_Method", dataset_name, paste0("appendix_summary_", dataset_name, ".csv")))
}