library(dplyr)
library(SKM)

# Multitrait --------------------------------------------------

testing_to_NA <- function(Pheno, fold) {
  for (trait in setdiff(names(fold), "num")) {
    testing_indices <- fold[[trait]]$testing

    Pheno[[trait]][testing_indices] <- NA
  }

  return(Pheno)
}

env_multivariate_summary <- function(observed, predicted) {
  Summary <- data.frame()

  for (trait in names(observed)) {
    summary <- numeric_summary(
        observed[[trait]],
        predicted[[trait]]$predicted
      ) %>%
      as.data.frame.list() %>%
      mutate(
        Trait = trait,
        nrmse = nrmse(
          observed[[trait]],
          predicted[[trait]]$predicted,
          type = "mean"
        )
      ) %>%
      relocate(Trait, 1)

    Summary <- rbind(Summary, summary)
  }

  return(Summary)
}

print.CVMultivariate <- function(folds) {
  folds_num <- length(folds)
  traits <- setdiff(names(folds[[1]]), "num")
  traits_num <- length(traits)

  for (i in seq(folds_num)) {
    cat("*** Fold", i, "***\n")
    folds[[i]]$num <- NULL
    lines <- c(
      attr(folds[[i]][[1]]$training, "lines"),
      attr(folds[[i]][[1]]$testing, "lines")
    )
    lines_num <- length(lines)

    Matrix <- matrix(1, lines_num, traits_num)
    rownames(Matrix) <- lines
    colnames(Matrix) <- traits

    fold <- folds[[i]]

    for (trait in traits) {
      Matrix[attr(fold[[trait]]$testing, "lines"), trait] <- 0
    }

    print(Matrix)
    cat("\n")
  }
}

multitrait_global_env_summary <- function(Results) {
  traits <- unique(Results$Trait)

  Summary <- data.frame()

  for (trait in traits) {
    Data <- Results %>%
      filter(Trait == trait) %>%
      droplevels()

    Summary <- Summary %>%
      rbind(global_env_summary(Data))
  }

  return(Summary)
}

# Unitrait --------------------------------------------------

env_summary <- function(observed, predicted) {
  Summary <- SKM::numeric_summary(observed, predicted$predicted) %>%
    as.data.frame.list() %>%
    mutate(
      nrmse = SKM::nrmse(
        observed,
        predicted$predicted,
        type = "mean"
      )
    )

  return(Summary)
}

prepare_results_by_fold <- function(Results,
                                    dataset,
                                    algorithm,
                                    model,
                                    data_used,
                                    testing_proportion,
                                    analysis_type,
                                    trait = NULL) {
  # For multitrait analysis that contains trait information in results itself
  if (is.null(trait)) {
    trait <- Results$Trait
  }

  Results <- Results %>%
    mutate(
      Dataset = dataset,
      Trait = trait,
      Algorithm = algorithm,
      Model = model,
      DataUsed = data_used,
      TestingProportion = testing_proportion,
      AnalysisType = analysis_type
    ) %>%
    rename(APC = "pearson") %>%
    rename_at(
      c("mse", "rmse", "nrmse", "maape", "mae"),
      toupper
    )

  return(Results)
}

global_env_summary <- function(Results) {
  Global <- Results %>%
    summarise_if(
      is.numeric,
      list(
        MEAN = ~ mean(., na.rm = TRUE),
        SE = ~ sd(., na.rm = TRUE) / sqrt(n())
      ),
      .groups = "keep"
    ) %>%
    rename_with(~gsub("_MEAN", "", .x)) %>%
    select(
      MSE, MSE_SE,
      RMSE, RMSE_SE,
      NRMSE, NRMSE_SE,
      MAAPE, MAAPE_SE,
      MAE, MAE_SE,
      APC, APC_SE
    )

  Info <- Results %>%
    slice(1) %>%
    select(
      Dataset,
      Trait,
      Model,
      DataUsed,
      Algorithm,
      TestingProportion,
      AnalysisType
    )

  return(data.frame(Info, Global))
}

# Utils --------------------------------------------------

divide <- function(n, divisor) {
  times <- rep(floor(n / divisor), divisor)
  module <- n %% divisor

  if (module > 0) {
    increase_one_indices <- seq(1, module)
    times[increase_one_indices] <- times[1] + 1
  }

  return(times)
}

env_testing_indices <- function(fold) {
  result <- c()
  for (env in setdiff(names(fold), "num")) {
    result <- c(result, fold[[env]]$testing)
  }

  return(result)
}

# Weighted functions -----------------------------------------------------------

weight_markers <- function(y, testing_indices, Line, Markers) {
  y_bin <- rep(0, length(y))
  y_bin[testing_indices] <- 1

  GRM <- Line %*% Markers
  X <- data.frame(y = y_bin, GRM)

  boruta_output <- Boruta::Boruta(y ~ ., data = X, doTrace = 0)
  boruta_signif <- Boruta::getSelectedAttributes(
    boruta_output,
    withTentative = FALSE
  )

  SummaryFeatures <- Boruta::attStats(boruta_output)
  Pos_IF1 <- which(SummaryFeatures$decision == "Confirmed")
  Pos_IF2 <- which(SummaryFeatures$decision == "Tentative")
  Pos_IF3 <- which(SummaryFeatures$decision == "Rejected")

  Pos_IF <- c(Pos_IF1, Pos_IF2, Pos_IF3)
  Sel_IF <- SummaryFeatures[Pos_IF, 1]

  names(Sel_IF) <- rownames(SummaryFeatures)
  Sel_IF <- Sel_IF + abs(min(Sel_IF)) + 1
  importance_sorted <- sort(Sel_IF, decreasing = TRUE)

  P <- length(importance_sorted)
  Inv_IF <- 1 / (importance_sorted)
  weight_Feature <- Inv_IF * P / sum(Inv_IF)

  Weights <- matrix(
    rep(weight_Feature, nrow(Markers)),
    nrow = nrow(Markers),
    ncol = P,
    byrow = TRUE
  )

  top_vars <- names(importance_sorted)
  WeightedMarkers <- SKM::remove_no_variance_cols(
    Markers[, top_vars] * Weights
  )

  return(WeightedMarkers)
}
