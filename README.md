
# Bayesian Transfer Learning For Genetic Codes

In this R project, we proved that it's possible to use Transfer Learning in order to improve the models, making some comparisons between the Conventional Bayesian Model, and Bayesian Transfer Model, using different metrics such as Correlation, Normalized Medium Square Error (NRMSE), PM_10, PM_20, PM_30, and so on.




## Installation

### dplyr
    install.packages("dplyr")

### lme4
    install.packages("lme4")

### ggplot2
    install.packages("ggplot2")

### devtools
    install.packages("devtools")

### BGLR
    devtools::install_github("gdlc/BGLR-R")

### Keras
    devtools::install_github("rstudio/keras")

### Tensorflow
    devtools::install_github("rstudio/tensorflow")

### SKM
    devtools::install_github("brandon-mosqueda/SKM")

### Crossdes
    install.packages("crossdes")



## How To Execute The Model

Run only the file Bayesian_Transfer_Learning_With_Covariate.R

### Select the dataset
#### Array of datasets, where index goes from 1 - 13
    datasets_list <- c(
                    "EYT_1", "EYT_2", "EYT_3", "Groundnut", "Indica", # NOT RESOLVED JAPONICA
                    "Japonica", "Maize", "Wheat_1", "Wheat_2", "Wheat_3",
                    "Wheat_4", "Wheat_5", "Wheat_6"
                  )
#### Select index on line 19
    datasets_list_index <- 5 #### CHANGE THIS VALUE TO CHANGE DATASET FILE

## License

[MIT](https://choosealicense.com/licenses/mit/)