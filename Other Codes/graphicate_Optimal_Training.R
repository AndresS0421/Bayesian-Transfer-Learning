#rm(list = ls(all = true))

library(dplyr)
library(ggplot2)

source("utils.r")

results1 <- read.csv("archivo_reorganizado_cor.csv")
head(results1) 
dim(results1)

pos_Methods_to_use=which(results1$Method %in% c("Conv", "Opt", "TrainSel"))
results3=results1[pos_Methods_to_use,]


data_sets<-unique(results3$Dataset)
data_sets

# delete all the previous plots
#unlink(plots_dir, recursive = true)
# create the directory where the plots are goind to be stored in
plots_dir<- paste("Optimal_Training_Graphics",sep="_")
dir.create(plots_dir)
for (data_set in data_sets) {
  
  data_set=data_set
  results5=results3[results3$Dataset==data_set, ]
  results5=droplevels(results5)

  cat(data_set, "\n")
  plot <- ggplot(
    results5,
    aes(x = Method, y = Cor, fill = factor(Sample_Size))
  ) +
    geom_bar(stat = "identity", position = "dodge") +  # Cambiado a "dodge" para dividir las barras
    facet_wrap(~Trait) +
    geom_errorbar(
      aes(
        ymin = Cor, ymax = Cor),
      width = 0.2,
      position = position_dodge(0.9),
      colour = "black"
    )+
    labs(fill = "Sample Size")  # Cambia la etiqueta de la leyenda
  plot <- white_theme(plot)
  Plot <- plot + theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold")
  ) +
    theme_bw() +
    theme(text = element_text(size = 18))
  Plot <- vertical_x(Plot, angle = 90)
  
  save_plot(Plot, file = file.path(plots_dir, paste0(data_set,'_COR', ".png")))
  
}

results1 <- read.csv("archivo_reorganizado_nrmse.csv")
results3=results1[pos_Methods_to_use,]

for (data_set in data_sets) {
  
  data_set <- data_set
  results5 <- results3[results3$Dataset == data_set, ]
  results5 <- droplevels(results5)
  
  cat(data_set, "\n")
  
  plot <- ggplot(
    results5,
    aes(x = Method, y = NRMSE, fill = factor(Sample_Size))
  ) +
    geom_bar(stat = "identity", position = "dodge") +  # Cambiado a "dodge" para dividir las barras
    facet_wrap(~Trait) +
    geom_errorbar(
      aes(
        ymin = NRMSE, ymax = NRMSE),
      width = 0.2,
      position = position_dodge(0.9),
      colour = "black"
    )+
    labs(fill = "Sample Size")  # Cambia la etiqueta de la leyenda
  
  plot <- white_theme(plot)
  Plot <- plot + theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold")
  ) +
    theme_bw() +
    theme(text = element_text(size = 18))
  Plot <- vertical_x(Plot, angle = 90)
  
  save_plot(Plot, file = file.path(plots_dir, paste0(data_set, '_NRMSE', ".png")))
}

