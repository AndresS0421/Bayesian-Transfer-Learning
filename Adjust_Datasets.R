datasets_list <- c(
  "EYT_1", "EYT_2", "EYT_3", "Groundnut", "Indica",
  "Maize", "Wheat_1", "Wheat_2", "Wheat_3",
  "Wheat_4", "Wheat_5", "Wheat_6"
)

for (dataset_name in datasets_list) {
  # Print Dataset Name
  cat("Dataset: ", dataset_name, "\n")
  
  # Get Dataset Info
  dataset_path <- file.path("Transfer_Method", dataset_name, paste0("summary_ALL_", dataset_name, ".csv"))
  dataset_file <- read.csv(dataset_path)
  cat("Dimensions: ", dim(dataset_file), "\n")
  
  # Adjust Dataset Info To Avoid Duplicated Values
  dataset_file <- unique(dataset_file[, -1])
  rownames(dataset_file) <- NULL
  cat("New Dimensions: ", dim(dataset_file), "\n", "\n")
  
  # Write File
  write.csv(dataset_file, dataset_path)
}