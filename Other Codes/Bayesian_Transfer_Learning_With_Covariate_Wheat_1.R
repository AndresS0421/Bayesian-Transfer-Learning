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

#####Loading data set#############
#####Loading data set#############
load("Dataset_Files/Wheat_3.RData")
ls()
Pheno=dat_ls$Pheno
colnames(Pheno)=c("Env","Line","GY")
Geno=dat_ls$Geno
Markers=dat_ls$Markers

head(Pheno)
######Selecting the traits to be evaluated###############
Traits_to_evaluate=colnames(Pheno)[c(3)]
Envs_to_evaluate=unique(Pheno$Env)

iterations_number <- 3000
burn_in <- 2500

##############Metric to compute percentage of matching###
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
dataset_file="Wheat_3"

###########Directory to save the results#####
results_dir <- file.path(
  "Transfer_Method",
  dataset_file)
mkdir(results_dir)

# Data preparation -------------------------------------------------------------
Predictions_Envs=data.frame()
Sumary_Envs=data.frame()
for (e in 2:length(Envs_to_evaluate)){
#e=2
 Env_e=Envs_to_evaluate[e]
Pheno_Proxy=Pheno[Pheno$Env==Env_e,]
Pheno_Proxy <- Pheno_Proxy %>%
  arrange(Env, Line)
dim(Pheno_Proxy)
head(Pheno_Proxy)
Pheno_Proxy=droplevels(Pheno_Proxy)
unique(Pheno_Proxy$Env)
Pheno_Proxy=as.data.frame(Pheno_Proxy)

####################
Env_1=Envs_to_evaluate[1]
Pheno_goal=Pheno[Pheno$Env==Env_1,]
Pheno_goal <- Pheno_goal%>%
  arrange(Env, Line)
dim(Pheno_goal)
Pheno_goal=droplevels(Pheno_goal)
unique(Pheno_goal$Env)
Pheno_goal=as.data.frame(Pheno_goal)

##############Training testing partitions

folds <- cv_random(length(Pheno_goal$Line),folds_number=10, testing_proportion =0.20)  
##########Sorting lines in Geno
geno_lines <- sort(rownames(Geno))
Geno <- Geno[geno_lines, geno_lines]
Predictions_Traits=data.frame()
Sumary_Traits=data.frame()
for (t in 1:length(Traits_to_evaluate)){
#  t=1
  Trait_t=Traits_to_evaluate[t]
#########Response variable in Obregon
y <- Pheno_Proxy[,Trait_t]
y_f <- y

Predictions_Final=data.frame()
Sumary_all=data.frame()
for(i in seq_along(folds )) {
#i=1

#########Training with the whole Obregon data set####################
##########Design matrix of lines
#### START CHANGES 1
ZL=model.matrix(~0+Line,data=Pheno_Proxy)
Lines_proxy=unique(Pheno_Proxy$Line)
Pos_proxy=which(Lines_proxy %in% rownames(Markers))
Markers_Proxy=Markers[Pos_proxy,]
dim(Markers_Proxy)
Markers_Proxy_Scale=scale(Markers_Proxy)

##########ETA1 con Env y Lines ==P1##############
ETA1=list(Line=list(model="BRR",X=Markers_Proxy_Scale))
##########ETA1 con Lines==P2 ##############
#ETA1=list(Line=list(model="FIXED",X=xD))
#### END CHANGES 1

##############Training with the whole oregon data set#############################
model_f<-BGLR::BGLR(
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

#### START CHANGES 2
ZL=model.matrix(~0+Line,data=Pheno_goal)
Lines_goal=unique(Pheno_goal$Line)
Pos_goal=which(Lines_goal %in% rownames(Geno))
Geno_goal=Geno[Pos_goal,Pos_goal]

Geno1=data.matrix(Geno_goal)
K_G=ZL%*%Geno1%*%t(ZL)
#### END CHANGES 2

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
#yy_s=yy-xD%*%Beta_Line
#### START CHANGES 3
yy_s=yy
y_fff=yy_s
y_fff[Tst_final] <- NA

Lines_goal=unique(Pheno_goal$Line)
Pos_goal=which(Lines_goal %in% rownames(Markers))
Markers_goal=Markers[Pos_goal,]
dim(Markers_goal)
Markers_goal_Scale=scale(Markers_goal)
XBeta_Cov=scale(Markers_goal_Scale%*%Beta_Line)
#XBeta_Cov=scale(c(y,y[1:166]))
ETA3=list(Line=list(model='RKHS',K=K_G),Cov=list(model='BRR',X=XBeta_Cov))
#### END CHANGES 3

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

Sumary=data.frame(
  Dataset =dataset_file,
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
Sumary
Sumary_all=rbind(Sumary_all,Sumary)

Predictions_i=data.frame(
  Dataset =dataset_file,
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
Sumary_all
head(Predictions_Final)

Tb_A =Sumary_all%>%group_by(Dataset,
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
Sumary_Traits=rbind(Sumary_Traits,Tb_A)
Sumary_Traits
Predictions_Traits=rbind(Predictions_Traits,Predictions_Final)
head(Predictions_Traits)


Sumary_Envs=rbind(Sumary_Envs, Sumary_Traits)
#Sumary_Envs
Predictions_Envs=rbind(Predictions_Envs,Predictions_Traits)
#head(Predictions_Envs)
}
Sumary_Envs
write.csv(Predictions_Envs, file.path(results_dir, "predictions_Wheat_1.csv"))
write.csv(Sumary_Envs , file.path(results_dir, "summary_Wheat_1.csv"))
