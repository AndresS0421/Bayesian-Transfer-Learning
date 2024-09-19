mse <- function(actual, predicted) mean((actual - predicted)^2, na.rm = TRUE)

rrmse <- function(actual, predicted) {
  sqrt(mse(actual, predicted)) / mean(actual, na.rm = TRUE)
}

mkdir <- function(directory) {
  if (!file.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
}

white_theme <- function(Plot) {
  Plot <- Plot + theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14, face = "bold")
    ) +
    theme_bw() +
    theme(text = element_text(size = 18))

  return(Plot)
}

vertical_x <- function(Plot, angle = 45) {
  Plot <- Plot + theme(
    axis.text.x = element_text(angle = angle, hjust = 1, vjust = 0.5)
  )

  return(Plot)
}

save_plot <- function(Plot, file, width = 1000, height = 850, res = 110) {
  png(file = file, width = width, height = height, res = res)
  print(Plot)
  dev.off()
}
