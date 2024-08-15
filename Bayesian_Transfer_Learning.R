rm(list = ls())

library(SKM)
#library(foreach)
#library(doParallel)
#library(plyr)
#library(tidyr)
#library(dplyr)
#library(reshape2)
library(dplyr)
source("utils.R")
source("cross_validators_env.R")

# FILES --------------------------------------------------------------------------
dataset_files_list <- c("TPE_1_2021_2022", "TPE_2_2021_2022", "TPE_3_2021_2022", 
                        "TPE_1_2022_2023", "TPE_2_2022_2023", "TPE_3_2022_2023"
                      )
dataset_folder <- "Dataset_Files"
sumary_results_file <- "summary_GENERAL"
results_folder <- "Transfer_Method"

# GENERAL FILES PATHS ------------------------------------------------------------
sumary_results_file_path <- paste(results_folder, sumary_results_file, "_Test.csv", sep = "/")
sumary_all_results_file_path <- paste(results_folder, sumary_results_file, "_ALL_Test.csv", sep = "/")
predictions_final_file_path <- paste(results_folder, "predictions_Test.csv", sep = "/")
# GENERAL DATA FRAMES ------------------------------------------------------------
general_sumary_results <- data.frame()
general_sumary_all_results <- data.frame()
general_predictions_results <- data.frame()

# PARAMETERS ---------------------------------------------------------------------
testing_proportions_list <- c(0.15, 0.3, 0.5, 0.75, 0.85)
models_list <- c("M1", "M2", "M3")
etas_list_1 <- c(
              'list(Env=list(model="FIXED",X=ZE),Line=list(model="FIXED",X=xD))', 
              'list(Line=list(model="FIXED",X=xD))'   
          )
etas_list_2 <- c(
              'list(Env=list(model="RKHS",K=K_E),Line=list(model="RKHS",K=K_G),GE=list(model="RKHS",K=K_GE))',
              'list(Env=list(model="RKHS",K=K_E),Line=list(model="RKHS",K=K_G))'
          )

# Metric to compute percentage of matching ---------------------------------------
best_lines_match <- function(Data, proportion = 0.1) {
  best_lines_number <- floor(nrow(Data) * proportion)
  
  best_lines <- Data %>%
    arrange(desc(Observed1)) %>%
    slice(seq(best_lines_number)) %>%
    pull(Line) %>%
    as.character()
  
  predicted_lines <- Data %>%
    arrange(desc(Predicted1)) %>%
    slice(seq(best_lines_number)) %>%
    pull(Line) %>%
    as.character()
  
  percentage_matching <- sum(predicted_lines %in% best_lines) /
    best_lines_number *
    100
  
  return(percentage_matching)
}


# GLOBAL PARAMS ------------------------------------------------------------------
data_used <- "All" # "NoObregon" is the other option
cores_num <- 10
folds_num <- 10
algorithm <- "BRR"
model_No <- ""
iterations_number <- 3000
burn_in <- 2500
base_results_dir <- "results"

# RUN ON EACH DATASET FILE -------------------------------------------------------
for (dataset_file_i in 1:length(dataset_files_list)) {
  # Select dataset file ----------------------------------------------------------
  dataset_file <- dataset_files_list[dataset_file_i]
  
  # RESULTS VARIABLE -------------------------------------------------------------
  Sumary_all <- data.frame()
  Predictions_Final=data.frame()
  
  # Directory to save the results ------------------------------------------------
  results_dir <- file.path(
    results_folder,
    dataset_file,
    model_No)
  mkdir(results_dir)
  
  # Data preparation -------------------------------------------------------------
  # OBREGON #
  load(sprintf("%s.RData", dataset_file), verbose = TRUE)
  Pheno
  dim(Pheno)
  Pheno_Obregon=Pheno[Pheno$TPE=="Obregon",]
  Pheno_Obregon <- Pheno_Obregon %>% arrange(Env, Line)
  dim(Pheno_Obregon)
  head(Pheno_Obregon)
  Pheno_Obregon=droplevels(Pheno_Obregon)
  unique(Pheno_Obregon$Env)
  Pheno_Obregon=as.data.frame(Pheno_Obregon)
  Trait=colnames(Pheno)[4]
  Trait
  # INDIA #
  Pheno_India=Pheno[Pheno$TPE==1,]
  Pheno_India <- Pheno_India %>% arrange(Env, Line)
  dim(Pheno_India)
  Pheno_India=droplevels(Pheno_India)
  unique(Pheno_India$Env)
  Pheno_India=as.data.frame(Pheno_India)
  Env <- model.matrix(~ 0 + Env, data = Pheno_Obregon)
  Line <- model.matrix(~ 0 + Line, data = Pheno_Obregon)
  # GENOMA #
  geno_lines <- sort(rownames(Geno))
  Geno <- Geno[geno_lines, geno_lines] %>% cholesky()
  # EXPECTED RESULTS #
  y <- as.numeric(Pheno_Obregon[,4])
  y_f <- y
  
  
  # TESTING PROPORTION
  for (t in 1:length(testing_proportions_list)) {
    testing_proportion <- testing_proportions_list[t]
    # MODELS
    for (m in 1:length(models_list)) {
      model_No <- models_list[m]
      # ETAs 1
      for (e1 in 1:length(etas_list_1)) {
        eta_1 <- etas_list_1[e1]
        #ETAs 2
        for (e2 in 1:length(etas_list_2)) {
          eta_2 <- etas_list_2[e2]
          
          ###### SHOW LOOP STEP ######
          print(paste("DATASET FILE - ", dataset_file))
          print(paste("TESTING PROPORTION - ", testing_proportion))
          print(paste("MODEL - ", model_No))
          print(paste("ETA 1 - ", e1))
          print(paste("ETA 2 - ", e2))
          ###### MODEL FOLDS (INDIA) ######
          folds <- folds_by_model(
            model = model_No,
            Pheno = Pheno_India,
            folds_num = folds_num,
            testing_proportion = testing_proportion
          )
          
          ######### MODEL TRAINING ##########
          for(i in 1:folds_num) {
            #i=1
            #y_f[Tst_final] <- NA
            
            ZL=model.matrix(~0+Line,data=Pheno_Obregon)
            Geno=data.matrix(Geno)
            K_G=ZL%*%Geno%*%t(ZL)
            ZE=model.matrix(~0+Env,data=Pheno_Obregon)
            K_E=ZE%*%t(ZE)
            K_GE=K_G*K_E
            
            X=svd(Geno)
            U=X$u
            d=X$d
            Q_var=quantile(d,probs=0.02)
            Q_var
            Pos_Q=which(d>Q_var)
            U_red=U[,Pos_Q]
            d_red=d[Pos_Q]
            D=diag(sqrt(d_red))
            xD1=U_red%*%D
            xD=ZL%*%xD1
            dim(xD)
            XGE=model.matrix(~0+xD:Env,data=Pheno_Obregon)
            
            ##########ETA1 con Env y Lines == P1##############
            #ETA=list(Env=list(model="FIXED",X=ZE),Line=list(model="FIXED",X=xD))
            ##########ETA1 con Lines == P2 ##############
            #ETA=list(Line=list(model="FIXED",X=xD))
            ETA = eval(parse(text = eta_1))
            
            ##############Training the regression model#############################
            #####Under the inner-cross-validation####################################
            model_f<-BGLR::BGLR(
              y = y_f,
              ETA = ETA,
              response_type = "gaussian",
              nIter = iterations_number,
              burnIn = burn_in,
              verbose = FALSE
            )
            
            #Beta_Env=model_f$ETA$Env$b
            Beta_Line=model_f$ETA$Line$b
            #Beta_LGE=model_f$ETA$GE$b
            #MU=model_f$mu
            
            Fold_testing=c(folds[[i]]$Delhi$testing,folds[[i]]$Karnal$testing,folds[[i]]$Ludhiana$testing,folds[[i]]$Nowshera$testing)
            
            Tst_final=Fold_testing
            yy <- as.numeric(Pheno_India[,4])
            y_ff=yy
            y_ff[Tst_final] <- NA
            
            ZL=model.matrix(~0+Line,data=Pheno_India)
            Geno=data.matrix(Geno)
            K_G=ZL%*%Geno%*%t(ZL)
            ZE=model.matrix(~0+Env,data=Pheno_India)
            K_E=ZE%*%t(ZE)
            K_GE=K_G*K_E
            
            #######ETA2 1 con interaccion==PI######
            #ETA2=list(Env=list(model='RKHS',K=K_E),Line=list(model='RKHS',K=K_G),GE=list(model='RKHS',K=K_GE))
            #######ETA2 2 sin interaccion===PWI######
            #ETA2=list(Env=list(model='RKHS',K=K_E),Line=list(model='RKHS',K=K_G))
            ETA2 = eval(parse(text = eta_2))
            
            ##############Training the regression model#############################
            #####Under the inner-cross-validation####################################
            model_ff<-BGLR::BGLR(
              y = y_ff,
              ETA = ETA2,
              response_type = "gaussian",
              nIter = iterations_number,
              burnIn = burn_in,
              verbose = FALSE
            )
            
            PredictedA = model_ff$yHat[Tst_final]
            ObservedA = yy[Tst_final]
            
            #U_Pred=data.frame(Line=Pheno$Line[fold$testing],Pred=g_Pred_testing1, Observed=g_True_testing1)
            Data_Conv=data.frame(Line=Pheno_India$Line[Tst_final],Observed1=ObservedA, Predicted1=PredictedA)
            #Data_Conv$Line=as.factor(Data_Conv$Line)
            Summary_Conv <- aggregate(cbind(Observed1, Predicted1) ~ Line, data = Data_Conv, FUN = mean, na.rm = TRUE)
            
            Observed=Summary_Conv$Observed1
            Predicted=Summary_Conv$Predicted1
            #####Metrics hole testing
            COR_BayI=cor(Observed,Predicted)
            MSE_BayI=mse(Observed,Predicted)
            NRMSE_BayI=nrmse(Observed,Predicted)
            
            ######Metrics top 20% testing####
            Top_20_value=quantile(Observed,probs=0.8)
            Pos_Q_tst=which(Observed>Top_20_value)
            COR_BayI20=cor(Observed[Pos_Q_tst],Predicted[Pos_Q_tst])
            MSE_BayI20=mse(Observed[Pos_Q_tst],Predicted[Pos_Q_tst])
            NRMSE_BayI20=nrmse(Observed[Pos_Q_tst],Predicted[Pos_Q_tst])
            
            ###########Percentage of metrics####
            
            PM_Conv=best_lines_match(Data=Summary_Conv,proportion = 0.2) 
            PM_Conv
            
            ################
            xD=ZL%*%xD1
            yy_s=yy-xD%*%Beta_Line
            y_fff=yy_s
            y_fff[Tst_final] <- NA
            ##############Training the regression model#############################
            #####Under the inner-cross-validation####################################
            model_fff<-BGLR::BGLR(
              y = y_fff,
              ETA = ETA2,
              response_type = "gaussian",
              nIter = iterations_number,
              burnIn = burn_in,
              verbose = FALSE
            )
            
            Predicted_s=model_fff$yHat[Tst_final]
            ObservedA=yy[Tst_final]
            Pred_Bs=xD%*%Beta_Line
            Predicted_TransfA=Pred_Bs[Tst_final]+Predicted_s
            
            Data_Trans=data.frame(Line=Pheno_India$Line[Tst_final], Observed1=ObservedA, Predicted1=Predicted_TransfA)
            #Data_Conv$Line=as.factor(Data_Conv$Line)
            Summary_Trans <- aggregate(cbind(Observed1, Predicted1) ~ Line, data =Data_Trans, FUN = mean, na.rm = TRUE)
            
            Observed=Summary_Trans$Observed1
            Predicted_Transf=Summary_Trans$Predicted1
            
            #####Metrics hole testing
            COR_BayI_T=cor(Observed,Predicted_Transf)
            MSE_BayI_T=mse(Observed,Predicted_Transf)
            NRMSE_BayI_T=nrmse(Observed,Predicted_Transf)
            
            ######Metrics top 20% testing####
            Top_20_value=quantile(Observed,probs=0.8)
            Pos_Q_tst=which(Observed>Top_20_value)
            COR_BayI20_T=cor(Observed[Pos_Q_tst],Predicted_Transf[Pos_Q_tst])
            MSE_BayI20_T=mse(Observed[Pos_Q_tst],Predicted_Transf[Pos_Q_tst])
            NRMSE_BayI20_T=nrmse(Observed[Pos_Q_tst],Predicted_Transf[Pos_Q_tst])
            
            ###########Percentage of metrics####
            
            PM_Transf=best_lines_match(Data=Summary_Trans,proportion = 0.2) 
            PM_Transf
            
            Sumary=data.frame(
              Dataset = dataset_file,
              Trait = Trait,
              Fold = i,
              Testing_Proportion = testing_proportion,
              Model = model_No,
              ETA_1 = e1,
              ETA_2 = e2,
              COR_BayI = COR_BayI,
              COR_Trans = COR_BayI_T,
              MSE_BayI = MSE_BayI,
              MSE_Trans = MSE_BayI_T,
              
              NRMSE_BayI = NRMSE_BayI,
              NRMSE_Trans = NRMSE_BayI_T,
              
              COR_BayI20 = COR_BayI20,
              COR_Trans20 = COR_BayI20_T,
              
              MSE_BayI20 = MSE_BayI20,
              MSE_Trans20 = MSE_BayI20_T,
              
              NRMSE_BayI20 = NRMSE_BayI20,
              NRMSE_Trans20 = NRMSE_BayI20_T,
              
              PM_Conv = PM_Conv,
              PM_Transf = PM_Transf
            )
            Sumary_all = rbind(Sumary_all, Sumary)
            
            Predictions_i = data.frame(
              Dataset = dataset_file,
              Trait = Trait,
              Env = Pheno_India$Env[Tst_final],
              Fold = i,
              Line = Pheno_India$Line[Tst_final],
              #Observed = Observed,
              Observed = ObservedA,
              PredictedI = PredictedA,
              Predicted_Transf = Predicted_TransfA
            )
            Predictions_Final <- rbind(Predictions_Final,
                                       Predictions_i
            )
          }
        }
      }
    }
  }
  
  Tb_A = Sumary_all%>%group_by(Dataset,
                               Trait) %>%summarise(COR_BayI_Ave=mean(COR_BayI, na.rm=TRUE),COR_BayNI_SD=sd(COR_BayI, na.rm=TRUE),
                                                   COR_Trans_Ave=mean(COR_Trans, na.rm=TRUE),COR_Trans_SD=sd(COR_Trans, na.rm=TRUE),
                                                   
                                                   MSE_BayI_Ave=mean( MSE_BayI, na.rm=TRUE), MSE_BayNI_SD=sd( MSE_BayI, na.rm=TRUE),
                                                   MSE_Trans_Ave=mean( MSE_Trans, na.rm=TRUE), MSE_Trans_SD=sd( MSE_Trans, na.rm=TRUE),
                                                   
                                                   NRMSE_BayI_Ave=mean(NRMSE_BayI, na.rm=TRUE), NRMSE_BayI_SD=sd(NRMSE_BayI, na.rm=TRUE),
                                                   NRMSE_Trans_Ave=mean(NRMSE_Trans, na.rm=TRUE), NRMSE_Trans_SD=sd(NRMSE_Trans, na.rm=TRUE),
                                                   
                                                   COR_BayI20_Ave=mean(COR_BayI20, na.rm=TRUE),COR_BayI20_SD=sd(COR_BayI20, na.rm=TRUE),
                                                   COR_Trans20_Ave=mean(COR_Trans20, na.rm=TRUE),COR_Trans20_SD=sd(COR_Trans20, na.rm=TRUE),
                                                   
                                                   MSE_BayI20_Ave=mean( MSE_BayI20, na.rm=TRUE), MSE_BayI20_SD=sd( MSE_BayI20, na.rm=TRUE),
                                                   MSE_Trans20_Ave=mean( MSE_Trans20, na.rm=TRUE), MSE_Trans20_SD=sd( MSE_Trans20, na.rm=TRUE),
                                                   
                                                   NRMSE_BayI20_Ave=mean(NRMSE_BayI20, na.rm=TRUE), NRMSE_BayI20_SD=sd(NRMSE_BayI20, na.rm=TRUE),
                                                   NRMSE_Trans20_Ave=mean(NRMSE_Trans20, na.rm=TRUE), NRMSE_Trans20_SD=sd(NRMSE_Trans20, na.rm=TRUE),
                                                   
                                                   PM_Conv_Ave=mean(PM_Conv, na.rm=TRUE), PM_Conv_SD=sd(PM_Conv, na.rm=TRUE),
                                                   PM_Transf_Ave=mean(PM_Transf, na.rm=TRUE), PM_Transf_SD=sd(PM_Transf, na.rm=TRUE)
                               )
  
  Tb_A=data.frame(Tb_A)
  
  # Write each dataset results file ----------------------------------------------
  write.csv(Predictions_Final, file.path(results_dir, paste("predictions_", dataset_file, "_Test.csv", sep = "")))
  write.csv(Tb_A  , file.path(results_dir, paste("summary_ALL_", dataset_file, "_Test.csv", sep = "")))
  write.csv(Sumary_all , file.path(results_dir, paste("summary_", dataset_file, "_Test.csv", sep = "")))
  # Add to general sumary results data frame
  general_sumary_results <- rbind(general_sumary_results, Sumary_all)
  general_sumary_all_results <- rbind(general_sumary_all_results, Tb_A)
  general_predictions_results <- rbind(general_predictions_results, Predictions_Final)
}

# Write the general results file -------------------------------------------------
write.csv(general_sumary_results, file.path(sumary_results_file_path))
write.csv(general_sumary_all_results, file.path(sumary_all_results_file_path))
write.csv(general_predictions_results, file.path(predictions_final_file_path))