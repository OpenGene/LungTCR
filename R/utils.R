#!/users/miniconda3/envs/r_env/bin/Rscript
# ref: https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
# ref: https://www.machinelearningplus.com/machine-learning/feature-selection/
# ref: https://www.datacareer.ch/blog/ridge-and-lasso-in-r/
# ref : https://remiller1450.github.io/s230f19/caret1.html
library(GenomicRanges)
library(glmnet)
library(caret)
library(pROC)
library(doSNOW)
library(tidyverse)
library(gam)
library(randomForest)
library(gbm)

set.seed(77) 

control <- 'Benign' 
cancer <- 'Malignant'


rebalance_table = function(input_table, fold){
  set.seed(77)
  type_counts_table <- table(input_table$Type)
  control_counts <- type_counts_table[control][[1]]
  cancer_counts <- type_counts_table[cancer][[1]]
  selected_cancer_counts <- as.integer(min(control_counts * fold, cancer_counts))
  print('Type for input table')
  print(table(input_table$Type))
  control_table <- input_table[which(input_table$Type == control), ]
  cancer_table <- sample_n(input_table[which(input_table$Type == cancer), ], selected_cancer_counts)
  # input_table after select samples
  selected_input_table <- rbind(cancer_table, control_table)
  print('Type for selected table')
  print(table(selected_input_table$Type))
  selected_input_table
}

qc_filter = function(input_table, stat_cols){
  stat_cols <- intersect(stat_cols, colnames(input_table))
  print('records before filter:')
  print(dim(input_table)[1])
  input_table <- input_table[which(input_table$Type == cancer | input_table$Type == control), ] # should use which here
  input_table <- input_table[which(input_table$"Raw_Yield(G)" >= 2 & input_table$"Raw_Yield(G)" < 15), ]
  input_table <- input_table[which(input_table$V_primer_percent >= 85 & input_table$J_primer_percent >=85 & input_table$merge_rate >= 0.7 & input_table$clone_reads_ratio >= 0.7), ]
  print('records after filter:')
  print(dim(input_table)[1])
  input_table <- input_table %>% mutate_at(stat_cols, ~replace_na(., 0))
  input_table <- input_table %>% mutate_at(stat_cols, ~replace(., . == '',0)) # replace black values with 0
  # input_table <- input_table %>% replace(is.na(.), 0)
  input_table <- as.data.frame(input_table)
  input_table$Type <- relevel(as.factor(input_table$Type),  control, levels= c(control, cancer))
  input_table
}

qc_filter2 = function(input_table, stat_cols){
  # only filter for type
  stat_cols <- intersect(stat_cols, colnames(input_table))
  print('records before filter:')
  print(dim(input_table)[1])
  input_table <- input_table[which(input_table$Type == cancer | input_table$Type == control), ] # should use which here
  print('records after filter:')
  print(dim(input_table)[1])
  input_table <- input_table %>% mutate_at(stat_cols, ~replace_na(., 0))
  input_table <- input_table %>% mutate_at(stat_cols, ~replace(., . == '',0)) # replace black values with 0
  # input_table <- input_table %>% replace(is.na(.), 0)
  input_table <- as.data.frame(input_table)
  input_table$Type <- relevel(as.factor(input_table$Type),  control, levels= c(control, cancer))
  input_table
}


get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

feature_selection_multiple_methods = function(input_table, stat_cols, output_dir, threads){
  cv_number <- 5
  repeated_num <- 3
  cl <- makeCluster(threads, type = "SOCK")
  stat_cols <- intersect(stat_cols, colnames(input_table))
  # print(stat_cols)
  # print(summary(input_table[,stat_cols]))
  # # Register cluster so that caret will know to train in parallel.
  registerDoSNOW(cl)
  lines <- c()
  feature_importance <- c()
  output_file <- paste0(output_dir, "/feature_selection.csv")
  lines <- cbind(lines, paste(sep=",", collapse=",",c('method', 'feature_number', 'features')))
  
  # recursive feature elimination
  # define the control using a random forest selection function
  subsets <- c(1:6) * 5
  # subsets <- c(1:3) * 5
  rfFuncs$summary <- twoClassSummary
  rfe_control <- rfeControl(functions=rfFuncs, method="repeatedcv", number=cv_number, saveDetails=TRUE, repeats=repeated_num, verbose=TRUE, rerank=TRUE)
  # rfe_control <- rfeControl(functions=nbFuncs, method="repeatedcv", number=cv_number, saveDetails=TRUE, repeats=repeated_num, verbose=TRUE, rerank=TRUE)
  # function options : linear regression (in the object lmFuncs), random forests (rfFuncs), naive Bayes (nbFuncs), bagged trees (treebagFuncs) and functions that can be used with caret’s train function (caretFuncs)
  # run the RFE algorithm
  rfe_model <- rfe(input_table[,stat_cols], input_table$Type, rfeControl=rfe_control, preProcess=c("center", "scale","corr", "nzv"), sizes=subsets, metric = "ROC")

  # print(varImp(object=rfe_model))
  defined_val <- 10
  defined_predictors <- unique(rfe_model$variables[rfe_model$variables$Variables==defined_val,]$var)
  lines <- cbind(lines, paste(sep=',', collapse=',', c('rfe_best', as.character(length(predictors(rfe_model))), paste(sep=';', collapse=';',predictors(rfe_model)))))
  lines <- cbind(lines, paste(sep=',', collapse=',', c('rfe_defined_predictors', as.character(length(defined_predictors)), paste(sep=';',  collapse=';',defined_predictors))))

  ### glmnet feature selection
  ctrl <- trainControl(verboseIter = TRUE, classProbs = TRUE, 
                      summaryFunction = twoClassSummary, method = "repeatedcv",
                      repeats = repeated_num,
                      number=cv_number)
  # glmnet_model <- train(input_table[,stat_cols], input_table$Type, method = "glmnet", metric = "Kappa", trControl = ctrl, preProcess=c("center", "scale","corr", "nzv"))
  glmnet_model <- train(input_table[,stat_cols], input_table$Type, method = "gbm", metric = "ROC", trControl = ctrl, preProcess=c("center", "scale","corr", "nzv"))
  print(glmnet_model)
  # print(varImp(object=glmnet_model))
  feature_importance <- cbind(feature_importance, as.matrix(varImp(object=glmnet_model)$importance))
  lines <- cbind(lines, paste(sep=',', collapse=',', c('glmnet', as.character(length(predictors(glmnet_model))), paste(sep=';', collapse=';',predictors(glmnet_model)))))

  # random forest
  randomForest_model <- train(input_table[,stat_cols], input_table$Type, method = "rf", metric = "ROC", trControl = ctrl,preProcess=c("center", "scale","corr", "nzv"))
  print(randomForest_model)
  print(varImp(object=randomForest_model))
  feature_importance <- cbind(feature_importance, as.matrix(varImp(object=randomForest_model)$importance))
  lines <- cbind(lines, paste(sep=',', collapse=',', c('randomForest', as.character(length(predictors(randomForest_model))), paste(sep=';', collapse=';',predictors(randomForest_model)))))
  stopCluster(cl)


  writeLines(lines, output_file, sep='\n')
  varimp_df <- as.data.frame(varImp(object=rfe_model))
  varimp_df$Overall <- varimp_df$Overall * (100 / max(varimp_df$Overall)) # normalize max score to 100
  feature_importance_all <- merge(as.data.frame(feature_importance), varimp_df, by='row.names', all=TRUE)
  # feature_importance_all <- merge(as.data.frame(feature_importance), as.data.frame(varImp(object=rfe_model)), by='row.names', all=TRUE)
  feature_importance_all <- merge(as.data.frame(feature_importance), varimp_df, by='row.names', all=TRUE)
  colnames(feature_importance_all) <- c('feature', 'glmnet', 'randomForest', 'rfe')
  feature_cols <- c('glmnet', 'randomForest', 'rfe')
  feature_importance_all <- feature_importance_all %>% mutate_at(feature_cols, ~replace_na(., 0))
  feature_importance_all$total <- feature_importance_all$glmnet + feature_importance_all$glmnet + feature_importance_all$rfe
  feature_importance_all <- feature_importance_all[order(feature_importance_all$total, decreasing = TRUE),]
  write.csv(feature_importance_all, file.path(output_dir, paste0("feature_importance.csv")),row.names=FALSE, quote=FALSE) 
}


get_seeds_lst = function (resample_num, model_num){
  set.seed(77)
  candidates <- model_num * 100
  seeds <- vector(mode="list", length=resample_num + 1)
  for (i in 1:resample_num) seeds[[i]] <- sample.int(candidates, model_num)
  seeds[[resample_num + 1]] <- sample.int(candidates, 1)
  seeds
}

feature_selection_multiple_methods_seeds = function(input_table, stat_cols, output_dir, threads){
  cv_number <- 5
  repeated_num <- 3
  cl <- makeCluster(threads, type = "SOCK")
  stat_cols <- intersect(stat_cols, colnames(input_table))

  registerDoSNOW(cl)
  lines <- c()
  feature_importance <- c()
  output_file <- paste0(output_dir, "/feature_selection.csv")
  lines <- cbind(lines, paste(sep=",", collapse=",",c('method', 'feature_number', 'features')))
  
  # recursive feature elimination
  # define the control using a random forest selection function
  subsets <- c(1:6) * 5
  rfe_seeds_lst <- get_seeds_lst(cv_number * repeated_num, length(subsets) + 1)
  # subsets <- c(1:3) * 5
  rfFuncs$summary <- twoClassSummary
  rfe_control <- rfeControl(functions=rfFuncs, method="repeatedcv", number=cv_number, saveDetails=TRUE, repeats=repeated_num, verbose=TRUE, rerank=TRUE, seeds=rfe_seeds_lst)

  rfe_model <- rfe(input_table[,stat_cols], input_table$Type, rfeControl=rfe_control, preProcess=c("center", "scale","corr", "nzv"), sizes=subsets, metric = "ROC")

  # print(varImp(object=rfe_model))
  defined_val <- 10
  defined_predictors <- unique(rfe_model$variables[rfe_model$variables$Variables==defined_val,]$var)
  lines <- cbind(lines, paste(sep=',', collapse=',', c('rfe_best', as.character(length(predictors(rfe_model))), paste(sep=';', collapse=';',predictors(rfe_model)))))
  lines <- cbind(lines, paste(sep=',', collapse=',', c('rfe_defined_predictors', as.character(length(defined_predictors)), paste(sep=';',  collapse=';',defined_predictors))))

  ### glmnet feature selection
  seeds_lst <- get_seeds_lst(cv_number * repeated_num, 100)
  ctrl <- trainControl(verboseIter = TRUE, classProbs = TRUE, 
                      summaryFunction = twoClassSummary, method = "repeatedcv",
                      repeats = repeated_num,
                      number=cv_number,
                      seeds=seeds_lst
                      )
  # glmnet_model <- train(input_table[,stat_cols], input_table$Type, method = "glmnet", metric = "Kappa", trControl = ctrl, preProcess=c("center", "scale","corr", "nzv"))
  glmnet_model <- train(input_table[,stat_cols], input_table$Type, method = "gbm", metric = "ROC", trControl = ctrl, preProcess=c("center", "scale","corr", "nzv"))
  print(glmnet_model)
  # print(varImp(object=glmnet_model))
  feature_importance <- cbind(feature_importance, as.matrix(varImp(object=glmnet_model)$importance))
  lines <- cbind(lines, paste(sep=',', collapse=',', c('glmnet', as.character(length(predictors(glmnet_model))), paste(sep=';', collapse=';',predictors(glmnet_model)))))

  # random forest
  randomForest_model <- train(input_table[,stat_cols], input_table$Type, method = "rf", metric = "ROC", trControl = ctrl,preProcess=c("center", "scale","corr", "nzv"))
  print(randomForest_model)
  print(varImp(object=randomForest_model))
  feature_importance <- cbind(feature_importance, as.matrix(varImp(object=randomForest_model)$importance))
  lines <- cbind(lines, paste(sep=',', collapse=',', c('randomForest', as.character(length(predictors(randomForest_model))), paste(sep=';', collapse=';',predictors(randomForest_model)))))
  stopCluster(cl)


  writeLines(lines, output_file, sep='\n')
  varimp_df <- as.data.frame(varImp(object=rfe_model))
  varimp_df$Overall <- varimp_df$Overall * (100 / max(varimp_df$Overall)) # normalize max score to 100
  feature_importance_all <- merge(as.data.frame(feature_importance), varimp_df, by='row.names', all=TRUE)
  # feature_importance_all <- merge(as.data.frame(feature_importance), as.data.frame(varImp(object=rfe_model)), by='row.names', all=TRUE)
  feature_importance_all <- merge(as.data.frame(feature_importance), varimp_df, by='row.names', all=TRUE)
  colnames(feature_importance_all) <- c('feature', 'glmnet', 'randomForest', 'rfe')
  feature_cols <- c('glmnet', 'randomForest', 'rfe')
  feature_importance_all <- feature_importance_all %>% mutate_at(feature_cols, ~replace_na(., 0))
  feature_importance_all$total <- feature_importance_all$glmnet + feature_importance_all$glmnet + feature_importance_all$rfe
  feature_importance_all <- feature_importance_all[order(feature_importance_all$total, decreasing = TRUE),]
  write.csv(feature_importance_all, file.path(output_dir, paste0("feature_importance.csv")),row.names=FALSE, quote=FALSE) 
}

