library(dplyr)

results <- read.csv("Transfer_Method/summary_GENERAL_ALL_V2.csv")
head(results)
dim(results)

averages_by_dataset <- results %>%
  group_by(Dataset) %>%
  summarise(
    COR_BayI_Ave = mean(COR_BayI_Ave, na.rm = TRUE),
    COR_BayNI_SD = mean(COR_BayNI_SD, na.rm = TRUE),
    COR_Trans_Ave = mean(COR_Trans_Ave, na.rm = TRUE),
    COR_Trans_SD = mean(COR_Trans_SD, na.rm = TRUE),
    MSE_BayI_Ave = mean(MSE_BayI_Ave, na.rm = TRUE),
    MSE_BayNI_SD = mean(MSE_BayNI_SD, na.rm = TRUE),
    MSE_Trans_Ave = mean(MSE_Trans_Ave, na.rm = TRUE),
    MSE_Trans_SD = mean(MSE_Trans_SD, na.rm = TRUE),
    NRMSE_BayI_Ave = mean(NRMSE_BayI_Ave, na.rm = TRUE),
    NRMSE_BayI_SD = mean(NRMSE_BayI_SD, na.rm = TRUE),
    NRMSE_Trans_Ave = mean(NRMSE_Trans_Ave, na.rm = TRUE),
    NRMSE_Trans_SD = mean(NRMSE_Trans_SD, na.rm = TRUE),
    PM_Conv_Ave_10 = mean(PM_Conv_Ave_10, na.rm = TRUE),
    PM_Conv_SD_10 = mean(PM_Conv_SD_10, na.rm = TRUE),
    PM_Transf_Ave_10 = mean(PM_Transf_Ave_10, na.rm = TRUE),
    PM_Transf_SD_10 = mean(PM_Transf_SD_10, na.rm = TRUE),
    PM_Conv_Ave_20 = mean(PM_Conv_Ave_20, na.rm = TRUE),
    PM_Conv_SD_20 = mean(PM_Conv_SD_20, na.rm = TRUE),
    PM_Transf_Ave_20 = mean(PM_Transf_Ave_20, na.rm = TRUE),
    PM_Transf_SD_20 = mean(PM_Transf_SD_20, na.rm = TRUE),
    PM_Conv_Ave_30 = mean(PM_Conv_Ave_30, na.rm = TRUE),
    PM_Conv_SD_30 = mean(PM_Conv_SD_30, na.rm = TRUE),
    PM_Transf_Ave_30 = mean(PM_Transf_Ave_30, na.rm = TRUE),
    PM_Transf_SD_30 = mean(PM_Transf_SD_30, na.rm = TRUE)
  )

averages_by_trait <- results %>%
  group_by(Dataset, Testing_Proportion, Trait) %>%
  summarise(
    COR_BayI_Ave = mean(COR_BayI_Ave, na.rm = TRUE),
    COR_BayNI_SD = mean(COR_BayNI_SD, na.rm = TRUE),
    COR_Trans_Ave = mean(COR_Trans_Ave, na.rm = TRUE),
    COR_Trans_SD = mean(COR_Trans_SD, na.rm = TRUE),
    MSE_BayI_Ave = mean(MSE_BayI_Ave, na.rm = TRUE),
    MSE_BayNI_SD = mean(MSE_BayNI_SD, na.rm = TRUE),
    MSE_Trans_Ave = mean(MSE_Trans_Ave, na.rm = TRUE),
    MSE_Trans_SD = mean(MSE_Trans_SD, na.rm = TRUE),
    NRMSE_BayI_Ave = mean(NRMSE_BayI_Ave, na.rm = TRUE),
    NRMSE_BayI_SD = mean(NRMSE_BayI_SD, na.rm = TRUE),
    NRMSE_Trans_Ave = mean(NRMSE_Trans_Ave, na.rm = TRUE),
    NRMSE_Trans_SD = mean(NRMSE_Trans_SD, na.rm = TRUE),
    PM_Conv_Ave_10 = mean(PM_Conv_Ave_10, na.rm = TRUE),
    PM_Conv_SD_10 = mean(PM_Conv_SD_10, na.rm = TRUE),
    PM_Transf_Ave_10 = mean(PM_Transf_Ave_10, na.rm = TRUE),
    PM_Transf_SD_10 = mean(PM_Transf_SD_10, na.rm = TRUE),
    PM_Conv_Ave_20 = mean(PM_Conv_Ave_20, na.rm = TRUE),
    PM_Conv_SD_20 = mean(PM_Conv_SD_20, na.rm = TRUE),
    PM_Transf_Ave_20 = mean(PM_Transf_Ave_20, na.rm = TRUE),
    PM_Transf_SD_20 = mean(PM_Transf_SD_20, na.rm = TRUE),
    PM_Conv_Ave_30 = mean(PM_Conv_Ave_30, na.rm = TRUE),
    PM_Conv_SD_30 = mean(PM_Conv_SD_30, na.rm = TRUE),
    PM_Transf_Ave_30 = mean(PM_Transf_Ave_30, na.rm = TRUE),
    PM_Transf_SD_30 = mean(PM_Transf_SD_30, na.rm = TRUE)
  )

write.csv(averages_by_dataset, file.path("Transfer_Method/summary_GENERAL_AVG_Dataset.csv"))
write.csv(averages_by_trait, file.path("Transfer_Method/summary_GENERAL_AVG_Trait.csv"))