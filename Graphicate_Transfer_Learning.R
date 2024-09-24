#rm(list = ls(all = true))

library(dplyr)
library(ggplot2)

source("utils.r")

results_trait <- read.csv("Transfer_Method/summary_GENERAL_AVG_Trait.csv")
results_dataset <- read.csv("Transfer_Method/summary_GENERAL_AVG_Dataset.csv")
head(results_trait) 
dim(results_trait)
head(results_dataset)
dim(results_dataset)

data_sets <- results_dataset$Dataset
data_sets <- unique(data_sets)


# delete all the previous plots
#unlink(plots_dir, recursive = true)
# create the directory where the plots are goind to be stored in
plots_dir<- paste("Results_Graphics",sep="_")
dir.create(plots_dir)

# Results Format To Graph --------------------------------------------------------
results_trait_long <- data.frame(
  Dataset = rep(results_trait$Dataset, 2),
  Trait = rep(results_trait$Trait, 2),
  Testing_Proportion = rep(results_trait$Testing_Proportion, 2),
  Method = rep(c("Conventional", "Transfer"), each = nrow(results_trait)),
  NRMSE = c(results_trait$NRMSE_BayI_Ave, results_trait$NRMSE_Trans_Ave),
  NRMSE_SD = c(results_trait$NRMSE_BayI_SD, results_trait$NRMSE_Trans_SD),
  COR = c(results_trait$COR_BayI_Ave, results_trait$COR_Trans_Ave),
  COR_SD = c(results_trait$COR_BayNI_SD, results_trait$COR_Trans_SD),
  PM_10 = c(results_trait$PM_Conv_Ave_10, results_trait$PM_Transf_Ave_10),
  PM_10_SD = c(results_trait$PM_Conv_SD_10, results_trait$PM_Transf_SD_10),
  PM_20 = c(results_trait$PM_Conv_Ave_20, results_trait$PM_Transf_Ave_20),
  PM_20_SD = c(results_trait$PM_Conv_SD_20, results_trait$PM_Transf_SD_20),
  PM_30 = c(results_trait$PM_Conv_Ave_30, results_trait$PM_Transf_Ave_30),
  PM_30_SD = c(results_trait$PM_Conv_SD_30, results_trait$PM_Transf_SD_30)
)

results_dataset_long <- data.frame(
  Dataset = rep(results_dataset$Dataset, 2),
  Trait = rep("Overall", nrow(results_dataset) * 2),
  Testing_Proportion = rep(results_dataset$Testing_Proportion, 2),
  Method = rep(c("Conventional", "Transfer"), each = nrow(results_dataset)),
  NRMSE = c(results_dataset$NRMSE_BayI_Ave, results_dataset$NRMSE_Trans_Ave),
  NRMSE_SD = c(results_dataset$NRMSE_BayI_SD, results_dataset$NRMSE_Trans_SD),
  COR = c(results_dataset$COR_BayI_Ave, results_dataset$COR_Trans_Ave),
  COR_SD = c(results_dataset$COR_BayNI_SD, results_dataset$COR_Trans_SD),
  PM_10 = c(results_dataset$PM_Conv_Ave_10, results_dataset$PM_Transf_Ave_10),
  PM_10_SD = c(results_dataset$PM_Conv_SD_10, results_dataset$PM_Transf_SD_10),
  PM_20 = c(results_dataset$PM_Conv_Ave_20, results_dataset$PM_Transf_Ave_20),
  PM_20_SD = c(results_dataset$PM_Conv_SD_20, results_dataset$PM_Transf_SD_20),
  PM_30 = c(results_dataset$PM_Conv_Ave_30, results_dataset$PM_Transf_Ave_30),
  PM_30_SD = c(results_dataset$PM_Conv_SD_30, results_dataset$PM_Transf_SD_30)
)

results_to_graph_types <- c("_By_Trait", "_By_Dataset")


# Graphs -------------------------------------------------------------------------
for (data_set in data_sets) {
  data_set = data_set
  results_to_graph <- list(droplevels(results_trait_long[results_trait_long$Dataset == data_set, ]), droplevels(results_dataset_long[results_dataset_long$Dataset == data_set, ]))
  
  cat(data_set, "\n")

  ### Create Dir -----------------------------------------------------------------
  dir.create(file.path(plots_dir, data_set))  
  ### Graph ALL Data -------------------------------------------------------------
  results_long_colnames <- colnames(results_trait_long)
  for (col_data_index in seq(5, length(results_long_colnames), by = 2)) {
    col_data <- results_long_colnames[col_data_index]
    col_data_sd <- results_long_colnames[col_data_index + 1]
    cat("Column: ", col_data, "\n")
    
    for (results_index in 1:length(results_to_graph)) {
      results <- results_to_graph[[results_index]]
      results_type <- results_to_graph_types[results_index]
      
      plot <- ggplot(
        results,
        aes(x = Method, y = !!sym(col_data), fill = factor(Testing_Proportion))
      ) +
        geom_bar(stat = "identity", position = "dodge") +  # Cambiado a "dodge" para dividir las barras
        facet_wrap(~Trait) +
        geom_errorbar(
          aes(
            ymin = !!sym(col_data) - !!sym(col_data_sd), 
            ymax = !!sym(col_data) + !!sym(col_data_sd)
          ),
          width = 0.2,
          position = position_dodge(0.9),
          colour = "black"
        )+
        labs(fill = "Proportions")  # Cambia la etiqueta de la leyenda
      plot <- white_theme(plot)
      Plot <- plot + theme(
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14, face = "bold")
      ) +
        theme_bw() +
        theme(text = element_text(size = 20))
      Plot <- vertical_x(Plot, angle = 90)
      
      ### Create Graph File ----------------------------------------------------------
      save_plot(Plot, file = file.path(plots_dir, data_set, paste0(data_set, paste0("_", col_data), results_type, ".png")))
    }
  }
}
