rm(list = ls())

library(SKM)
#library(foreach)
#library(doParallel)
#library(plyr)
#library(tidyr)
#library(dplyr)
#library(reshape2)
library(dplyr)
library(BGLR)


iterations_number <- 3000
burn_in <- 2500
# Files parameters ---------------------------------------------------------------
datasets_list <- c(
                    "EYT_1", "EYT_2", "EYT_3", "Groundnut", "Indica", # NOT RESOLVED JAPONICA
                    "Japonica", "Maize", "Wheat_1", "Wheat_2", "Wheat_3",
                    "Wheat_4", "Wheat_5", "Wheat_6"
                  )
dataset_folder <- "Dataset_Files"
summary_results_file <- "summary_GENERAL"
results_folder <- "Transfer_Method"
datasets_list_index <- 5 #### CHANGE THIS VALUE TO CHANGE DATASET FILE
dataset_file_name <- datasets_list[datasets_list_index]
dataset_file <- paste(dataset_folder, dataset_file_name, sep = "/")

# Global parameters --------------------------------------------------------------
testing_proportions_list <- c(0.15, 0.3, 0.5, 0.75, 0.8)

# GENERAL FILES PATHS ------------------------------------------------------------
summary_results_file_path <- paste(results_folder, paste(summary_results_file, ".csv", sep = ""), sep = "/")
summary_all_results_file_path <- paste(results_folder, paste(summary_results_file, "_ALL.csv", sep = ""), sep = "/")
predictions_final_file_path <- paste(results_folder, "predictions.csv", sep = "/")
# GENERAL DATA FRAMES ------------------------------------------------------------
general_summary_results <- data.frame()
general_summary_all_results <- data.frame()
general_predictions_results <- data.frame()
# Get general files data ---------------------------------------------------------
if (file.exists(summary_results_file_path)) {
  general_summary_results <- read.csv(summary_results_file_path) %>% select(-X)
  general_summary_all_results <- read.csv(summary_all_results_file_path) %>% select(-X)
  general_predictions_results <- read.csv(predictions_final_file_path) %>% select(-X)
}

# Directory to save the results --------------------------------------------------
results_dir <- file.path(
  results_folder,
  dataset_file_name)
mkdir(results_dir)

# Metric to compute percentage of matching ---------------------------------------
best_lines_match <- function(Data, proportion = 0.1) {
  best_lines_number <- floor(nrow(Data) * proportion)
  
  best_lines <- Data %>%
    arrange(desc(Observed)) %>%
    slice(seq(best_lines_number)) %>%
    pull(Line) %>%
    as.character()
  
  predicted_lines <- Data %>%
    arrange(desc(Predicted)) %>%
    slice(seq(best_lines_number)) %>%
    pull(Line) %>%
    as.character()
  
  percentage_matching <- sum(predicted_lines %in% best_lines) /
    best_lines_number *
    100
  
  return(percentage_matching)
}

# Results variables --------------------------------------------------------------
Predictions_Final=data.frame()
Summary_all=data.frame()
Tb_A_all=data.frame()

# Loading data set ---------------------------------------------------------------
Pheno = data.frame()
Geno = data.frame()
Markers <- data.frame()
load(sprintf("%s.RData", dataset_file), verbose = TRUE)
ls()

# Traits to be evaluated data structure initialization ---------------------------
Traits_to_evaluate <- data.frame()
# Wheat addaptions ---------------------------------------------------------------
if (grepl("Wheat", dataset_file_name)) {
  print("Wheat")
  Pheno = dat_ls$Pheno
  Pheno$Line <- Pheno$GID
  Geno = dat_ls$Geno
  Markers=dat_ls$Markers
  head(Geno)
  ### GET ONLY GY COLUMN DATA
  Traits_to_evaluate <- colnames(Pheno)[3]
} else {
  ### GET ALL COLUMNS TO BE EVALUATED
  Traits_to_evaluate=colnames(Pheno)[c(4:length(colnames(Pheno)))]
}

# Selecting the traits to be evaluated -------------------------------------------
Envs_to_evaluate=unique(Pheno$Env)

# TESTING PROPORTION
for (t in 1:length(testing_proportions_list)) {
  testing_proportion_value <- testing_proportions_list[t]
  
  ### SHOW LOOP STEP
  print(paste("TESTING PROPORTION - ", testing_proportion_value))
  
  # Data preparation -------------------------------------------------------------
  Predictions_Envs=data.frame()
  Summary_Envs=data.frame()
  
  #
  for (e in 2:length(Envs_to_evaluate)){
    #e=2
    Env_e=Envs_to_evaluate[e]
    Pheno_Proxy=Pheno[Pheno$Env==Env_e,]
    Pheno_Proxy <- Pheno_Proxy %>% arrange(Env, Line)
    dim(Pheno_Proxy)
    head(Pheno_Proxy)
    Pheno_Proxy=droplevels(Pheno_Proxy)
    unique(Pheno_Proxy$Env)
    Pheno_Proxy=as.data.frame(Pheno_Proxy)
    
    ####################
    Env_1=Envs_to_evaluate[1]
    Pheno_goal=Pheno[Pheno$Env==Env_1,]
    Pheno_goal <- Pheno_goal %>% arrange(Env, Line)
    dim(Pheno_goal)
    Pheno_goal=droplevels(Pheno_goal)
    unique(Pheno_goal$Env)
    Pheno_goal=as.data.frame(Pheno_goal)
    
    ### Training testing partitions
    folds <- cv_random(length(Pheno_goal$Line),folds_number=10, testing_proportion = testing_proportion_value) 
    
    ### Sorting lines in Geno
    geno_lines <- sort(rownames(Geno))
    Geno <- Geno[geno_lines, geno_lines]
    Predictions_Traits=data.frame()
    Summary_Traits=data.frame()
    for (t in 1:length(Traits_to_evaluate)){
      #  t=1
      Trait_t=Traits_to_evaluate[t]
      ### Response variable in Obregon
      y <- Pheno_Proxy[,Trait_t]
      y_f <- y
      
      for(i in seq_along(folds)) {
        #i=1
        
        #########Training with the whole Obregon data set####################
        ZL <- data.frame()
        Lines_proxy <- data.frame()
        Pos_proxy <- data.frame()
        Markers_Proxy <- data.frame()
        Markers_Proxy_Scale <- data.frame()
        ETA1 <- list()
        
        if (grepl("Wheat", dataset_file_name) || grepl("Japonica", dataset_file_name)) {
          ZL=model.matrix(~0+Line,data=Pheno_Proxy)
          Lines_proxy=unique(Pheno_Proxy$Line)
          Pos_proxy=which(Lines_proxy %in% rownames(Markers))
          Markers_Proxy=Markers[Pos_proxy,]
          dim(Markers_Proxy)
          Markers_Proxy_Scale=scale(Markers_Proxy)
          
          ##########ETA1 con Env y Lines ==P1##############
          ETA1=list(Line=list(model="BRR",X=Markers_Proxy_Scale))
        } else {
          ##########Design matrix of lines
          ZL=model.matrix(~0+Line,data=Pheno_Proxy)
          Geno = data.matrix(Geno) 
          
          ############Singular Value Descomposition of Geno
          X=svd(Geno)
          U=X$u
          d=X$d
          Q_var=quantile(d,probs=0.0)
          Q_var
          Pos_Q=which(d>Q_var)
          U_red=U[,Pos_Q]
          d_red=d[Pos_Q]
          D=diag(sqrt(d_red))
          xD1=U_red%*%D
          xD=ZL%*%xD1
          ##########ETA1 con Env y Lines ==P1##############
          ETA1=list(Line=list(model="BRR",X=xD))
        }
        
        ##############Training with the whole obregon data set#############################
        model_f <- BGLR::BGLR(
          y = y_f,
          ETA = ETA1,
          response_type = "gaussian",
          nIter = iterations_number,
          burnIn = burn_in,
          verbose = FALSE
        )
        ###########Betas learned with the whole obregon data set###############3
        Beta_Line=model_f$ETA$Line$b
        fold_i=folds[[i]]
        
        testing_indices <-fold_i$testing
        ###############Testing data set######
        Tst_final=testing_indices
        yy <- Pheno_goal[,Trait_t]
        y_ff=yy
        y_ff[Tst_final] <- NA
        
        Lines_goal <- data.frame()
        Pos_goal <- data.frame()
        Geno_goal <- data.frame()
        Geno1 <- data.frame()
        K_G <- data.frame()
        if (grepl("Wheat", dataset_file_name) || grepl("Japonica", dataset_file_name)) {
          ZL=model.matrix(~0+Line,data=Pheno_goal)
          Lines_goal=unique(Pheno_goal$Line)
          Pos_goal=which(Lines_goal %in% rownames(Geno))
          Geno_goal=Geno[Pos_goal,Pos_goal]
          
          Geno1=data.matrix(Geno_goal)
          K_G=ZL%*%Geno1%*%t(ZL)
        } else {
          ZL=model.matrix(~0+Line,data=Pheno_goal)
          Geno=data.matrix(Geno)
          K_G=ZL%*%Geno%*%t(ZL)
        }
        
        #######ETA2 2 sin interaccion===PWI######
        ETA2=list(Line=list(model='RKHS',K=K_G)) 
        
        ##############Training the regression model with India data#############################
        model_ff<-BGLR::BGLR(
          y = y_ff,
          ETA = ETA2,
          response_type = "gaussian",
          nIter = iterations_number,
          burnIn = burn_in,
          verbose = FALSE
        )
        
        PredictedA=model_ff$yHat[Tst_final]
        ObservedA=yy[Tst_final]
        
        ################Training with transfer learning################
        yy_s=yy
        y_fff=yy_s
        y_fff[Tst_final] <- NA
        
        Markers_goal <- data.frame()
        Markers_goal_Scale <- data.frame()
        XBeta_Cov <- data.frame()
        if (grepl("Wheat", dataset_file_name) || grepl("Japonica", dataset_file_name)) {
          Lines_goal=unique(Pheno_goal$Line)
          Pos_goal=which(Lines_goal %in% rownames(Markers))
          Markers_goal=Markers[Pos_goal,]
          dim(Markers_goal)
          Markers_goal_Scale=scale(Markers_goal)
          XBeta_Cov=scale(Markers_goal_Scale%*%Beta_Line)
        } else {
          XBeta_Cov=scale(xD%*%Beta_Line) 
        }
        ETA3=list(Line=list(model='RKHS',K=K_G),Cov=list(model='BRR',X=XBeta_Cov))
        
        ##############Training the regression model#############################
        #####Under the inner-cross-validation####################################
        model_fff<-BGLR::BGLR(
          y = y_fff,
          ETA = ETA3,
          response_type = "gaussian",
          nIter = iterations_number,
          burnIn = burn_in,
          verbose = FALSE
        )
        
        Predicted_s=model_fff$yHat[Tst_final]
        ObservedA=yy[Tst_final]
        Predicted_TransfA=Predicted_s
        
        #############Metrics computation
        
        #U_Pred=data.frame(Line=Pheno$Line[fold$testing],Pred=g_Pred_testing1, Observed=g_True_testing1)
        Data_Conv=data.frame(Line=Pheno_goal$Line[Tst_final],Observed=ObservedA, Predicted=PredictedA)
        
        #####Metrics hole testing
        COR_BayI=cor(ObservedA,PredictedA)
        MSE_BayI=mse(ObservedA,PredictedA)
        NRMSE_BayI=nrmse(ObservedA,PredictedA, type="mean")
        
        
        ###########Percentage of metrics####
        
        PM_Conv_10=best_lines_match(Data=Data_Conv,proportion = 0.1) 
        PM_Conv_10
        PM_Conv_20=best_lines_match(Data=Data_Conv,proportion = 0.2) 
        PM_Conv_20
        PM_Conv_30=best_lines_match(Data=Data_Conv,proportion = 0.3) 
        PM_Conv_30
        
        Data_Trans=data.frame(Line=Pheno_goal$Line[Tst_final],Observed=ObservedA, Predicted=Predicted_TransfA)
        
        #####Metrics hole testing
        COR_BayI_T=cor(ObservedA,Predicted_TransfA)
        MSE_BayI_T=mse(ObservedA,Predicted_TransfA)
        NRMSE_BayI_T=nrmse(ObservedA,Predicted_TransfA, type="mean")
        
        ###########Percentage of metrics####
        
        PM_Transf_10=best_lines_match(Data=Data_Trans,proportion = 0.1) 
        PM_Transf_10
        PM_Transf_20=best_lines_match(Data=Data_Trans,proportion = 0.2) 
        PM_Transf_20
        PM_Transf_30=best_lines_match(Data=Data_Trans,proportion = 0.3) 
        PM_Transf_30
        
        Summary=data.frame(
          Dataset =dataset_file_name,
          Testing_Proportion = testing_proportion_value,
          Trait = Trait_t,
          Env_Proxy =Env_e,
          Env_goal =Env_1,
          Fold=i,
          COR_BayI=COR_BayI,
          COR_Trans=COR_BayI_T,
          MSE_BayI=MSE_BayI,
          MSE_Trans=MSE_BayI_T,
          
          NRMSE_BayI=NRMSE_BayI,
          NRMSE_Trans=NRMSE_BayI_T,
          
          PM_Conv_10=PM_Conv_10,
          PM_Transf_10=PM_Transf_10,
          PM_Conv_20=PM_Conv_20,
          PM_Transf_20=PM_Transf_20,
          PM_Conv_30=PM_Conv_30,
          PM_Transf_30=PM_Transf_30
          
        )
        Summary
        Summary_all=rbind(Summary_all,Summary)
        
        Predictions_i=data.frame(
          Dataset =dataset_file_name,
          Testing_Proportion = testing_proportion_value,
          Trait =Trait_t,
          Env =Pheno_goal$Env[Tst_final],
          Fold=i,
          Line = Pheno_goal$Line[Tst_final],
          Observed =ObservedA,
          PredictedI =PredictedA,
          Predicted_Transf =Predicted_TransfA
        )
        Predictions_Final <-rbind(Predictions_Final,
                                  Predictions_i
        )
      }
      Summary_all
      head(Predictions_Final)
      
      Tb_A =Summary_all%>%group_by(Dataset, Testing_Proportion = testing_proportion_value,
                                  Trait,Env_Proxy, Env_goal) %>%summarise(COR_BayI_Ave=mean(COR_BayI, na.rm=TRUE),COR_BayNI_SD=sd(COR_BayI, na.rm=TRUE),
                                                                          COR_Trans_Ave=mean(COR_Trans, na.rm=TRUE),COR_Trans_SD=sd(COR_Trans, na.rm=TRUE),
                                                                          
                                                                          MSE_BayI_Ave=mean( MSE_BayI, na.rm=TRUE), MSE_BayNI_SD=sd( MSE_BayI, na.rm=TRUE),
                                                                          MSE_Trans_Ave=mean( MSE_Trans, na.rm=TRUE), MSE_Trans_SD=sd( MSE_Trans, na.rm=TRUE),
                                                                          
                                                                          NRMSE_BayI_Ave=mean(NRMSE_BayI, na.rm=TRUE), NRMSE_BayI_SD=sd(NRMSE_BayI, na.rm=TRUE),
                                                                          NRMSE_Trans_Ave=mean(NRMSE_Trans, na.rm=TRUE), NRMSE_Trans_SD=sd(NRMSE_Trans, na.rm=TRUE),
                                                                          
                                                                          PM_Conv_Ave_10=mean(PM_Conv_10, na.rm=TRUE), PM_Conv_SD_10=sd(PM_Conv_10, na.rm=TRUE),
                                                                          PM_Transf_Ave_10=mean(PM_Transf_10, na.rm=TRUE), PM_Transf_SD_10=sd(PM_Transf_10, na.rm=TRUE),
                                                                          
                                                                          PM_Conv_Ave_20=mean(PM_Conv_20, na.rm=TRUE), PM_Conv_SD_20=sd(PM_Conv_20, na.rm=TRUE),
                                                                          PM_Transf_Ave_20=mean(PM_Transf_20, na.rm=TRUE), PM_Transf_SD_20=sd(PM_Transf_20, na.rm=TRUE),
                                                                          
                                                                          PM_Conv_Ave_30=mean(PM_Conv_30, na.rm=TRUE), PM_Conv_SD_30=sd(PM_Conv_30, na.rm=TRUE),
                                                                          PM_Transf_Ave_30=mean(PM_Transf_30, na.rm=TRUE), PM_Transf_SD_30=sd(PM_Transf_30, na.rm=TRUE)
                                  )
      Tb_A=data.frame(Tb_A)
    }
    Summary_Traits=rbind(Summary_Traits,Tb_A)
    Summary_Traits
    Predictions_Traits=rbind(Predictions_Traits,Predictions_Final)
    head(Predictions_Traits)
    
    
    Summary_Envs=rbind(Summary_Envs, Summary_Traits)
    Predictions_Envs=rbind(Predictions_Envs,Predictions_Traits)
  }
  Tb_A_all <- rbind(Tb_A_all, Tb_A)
}

# Write each dataset results file ------------------------------------------------
write.csv(Predictions_Final, paste(results_dir, paste("predictions_", dataset_file_name, ".csv", sep = ""), sep = "/"))
write.csv(Tb_A_all  , paste(results_dir, paste("summary_ALL_", dataset_file_name, ".csv", sep = ""), sep = "/"))
write.csv(Summary_all , paste(results_dir, paste("summary_", dataset_file_name, ".csv", sep = ""), sep = "/"))
# Add to general summary results data frame
general_summary_results <- rbind(general_summary_results, Summary_all)
general_summary_all_results <- rbind(general_summary_all_results, Tb_A_all)
general_predictions_results <- rbind(general_predictions_results, Predictions_Final)

# Write the general results file -------------------------------------------------
write.csv(general_summary_results, file.path(summary_results_file_path))
write.csv(general_summary_all_results, file.path(summary_all_results_file_path))
write.csv(general_predictions_results, file.path(predictions_final_file_path))