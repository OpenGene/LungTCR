#!/users/miniconda3/envs/r_env/bin/Rscript
# ref: https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
# ref: https://www.machinelearningplus.com/machine-learning/feature-selection/
# ref: https://www.datacareer.ch/blog/ridge-and-lasso-in-r/
# ref : https://remiller1450.github.io/s230f19/caret1.html
library(GenomicRanges)
library(caret)
library(pROC)
library(doSNOW)
library(tidyverse)
library(gbm)
library(randomForest)
library(kernlab)
script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "script", "utils.R"), chdir = FALSE)
set.seed(77) 

args <- commandArgs(trailingOnly = TRUE)

train_file <- args[1]
print(train_file)
test_file <- args[2]
print(test_file)
val_file <- args[3]
print(val_file)
output_dir <- args[4]
threads <- as.integer(args[5])
feature_file <- args[6]
resample_tag <- args[7]

control <- "Benign"
cancer <- "Malignant"
lines <- c()
cl <- makeCluster(threads, type = "SOCK")
registerDoSNOW(cl)



qc_filter2 = function(input_table){
  print('records before filter:')
  print(dim(input_table)[1])
  input_table <- input_table[which(input_table$Type == cancer | input_table$Type == control), ] # should use which here
  print('records after filter:')
  print(dim(input_table)[1])
  input_table <- input_table %>% mutate_at(stat_cols, ~replace_na(., 0))
  input_table <- input_table %>% mutate_at(stat_cols, ~replace(., . == '',0)) # replace black values with 0
  # input_table <- input_table %>% replace(is.na(.), 0)
  input_table <- as.data.frame(input_table)
  input_table$Type <- relevel(as.factor(input_table$Type), cancer, levels= c(control,cancer))
  input_table
}

get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}


model_train = function(X, y, method){
# Register cluster so that caret will know to train in parallel.
  model <- caret::train(X, y, method = method,preProcess = c("corr", "nzv", "center", "scale"),metric = "ROC",trControl = ctrl)
  print(model)
  print(varImp(model))
  print(confusionMatrix(model, mode='everything', positive=cancer))
  model
}



model_test_val = function(model, X){
  test_pred <- predict(model, X, type='prob')
  test_pred
}


get_roc_auc = function(test_pred, y){
  roc_auc <- roc(predictor=test_pred[[cancer]], response=y, levels= c(control,cancer), direction = "<", quiet = FALSE)
  roc_auc
}



model_test_val_total = function(models_list, data, tag){
  sprintf("%d Samples for Test Data", dim(data)[1])

  # svg(file.path(output_dir, paste0(tag, ".svg")))
  png(file.path(output_dir, paste0(tag, ".png")))
  color_palettes <- c('blue', 'red', 'green', 'yellow', 'violet')
  pred_val_matrix <- c()
  if ("Anno" %in% names(data)){
    pred_val_df <- data %>% select("Sample_ID", "Type", "Name", "Anno")
  }
  else{
    pred_val_df <- data %>% select("Sample_ID", "Name", "Type")
  }
  # 
  i <- 0
  total_metrics <- c()
  for (model_name in names(models_list)){
    print(model_name)
    i <- i + 1
    model <- models_list[[model_name]]
    # models_list$model_name
    pred_val <- model_test_val(model, data[,stat_cols])
    pred_val_df$Score <- pred_val[[cancer]]
    names(pred_val_df)[dim(pred_val_df)[2]] = model_name # rename the last col to model name, can not directly assign virable as col name in R ...
    roc_auc <- get_roc_auc(pred_val, data$Type)
    print(roc_auc)
    model_auc_metrics <- as.matrix(roc_auc$auc)
    total_metrics <- cbind(total_metrics, model_auc_metrics)
    if (i == 1) {
      plot(roc_auc, col = color_palettes[[i]])
      }
    else
    {
      lines(roc_auc, col = color_palettes[[i]])
    }
  }
  total_metrics <- format(round(total_metrics, 2))
  notes <- c(paste(names(models_list), total_metrics))
  legend("bottomright", legend = notes, col = color_palettes, lty = 1)
  dev.off()
  colnames(total_metrics) <- c('randomForest', 'gbm', 'glmnet', 'svmLinear', 'svmRadial')
  write.csv(total_metrics,file.path(output_dir, paste0(tag, "_data_model_metrics.csv")),row.names=TRUE, quote=FALSE)
  write.csv(pred_val_df,file.path(output_dir, paste0(tag, "_pred_val.csv")),row.names=FALSE, quote=FALSE)
}

model_test_group = function(model, X, y){
  test_pred <- predict(model, X)
  confusionMatrix(reference = y, data = test_pred, mode='everything', positive=cancer)
}


model_test_group_total = function(models_list, data, tag)
{
  # for test data
  sprintf("The Predicted Confusion matrix for %s samples", tag)
  sprintf("%d Samples for Test Data", dim(data)[1])
  total_metrics <- c()
  for (model_name in names(models_list)){
    print(model_name)
    model <- models_list[[model_name]]
    # models_list$model_name
    model_cm <- model_test_group(model, data[,stat_cols], data$Type)
    print(model_cm)
    model_cm_metrics <- as.matrix(model_cm$byClass)
    total_metrics <- cbind(total_metrics, model_cm_metrics)
  }
  colnames(total_metrics) <- c('randomForest', 'gbm', 'glmnet', 'svmLinear', 'svmRadial')
  write.csv(total_metrics,file.path(output_dir, paste0(tag, "_data_model_metrics_group.csv")),row.names=TRUE, quote=FALSE)
}

get_seeds_lst = function (resample_num, model_num){
  set.seed(77)
  candidates <- model_num * 100
  seeds <- vector(mode="list", length=resample_num + 1)
  for (i in 1:resample_num) seeds[[i]] <- sample.int(candidates, model_num)
  seeds[[resample_num + 1]] <- sample.int(candidates, 1)
  seeds
}

trainning_data <- read_csv(train_file)
test_data <- read_csv(test_file)
feature_table <- read_csv(feature_file)
val_data <- read_csv(val_file)
features <- as.vector(unlist(feature_table['feature']))
# cols <- get_selected_features(features, 10)
stat_cols <- features


stat_cols <- intersect(stat_cols, colnames(trainning_data))
print('Features for prediction')
print(stat_cols)
# filter
trainning_data <- qc_filter2(trainning_data)
print(resample_tag)
if (resample_tag == "undersample"){
  print('Undersample')
  trainning_data <- rebalance_table(trainning_data, 2)
}else{
  print('Do not undersample')
}


# get balanced healthy and cancer samples



test_data <- qc_filter2(test_data)
val_data <- qc_filter2(val_data)

sample_type_counts <- c()
sample_type_counts <- cbind(sample_type_counts, table(trainning_data$Type))
sample_type_counts <- cbind(sample_type_counts, table(test_data$Type))
sample_type_counts <- cbind(sample_type_counts, table(val_data$Type))
colnames(sample_type_counts) <- c('train', 'test', 'validation')
write.csv(sample_type_counts,file.path(output_dir, paste0("sample_type.csv")),row.names=TRUE, quote=FALSE)


cv_number <- 5
repeated_num <- 5

# # Register cluster so that caret will know to train in parallel.
seeds_lst <- get_seeds_lst(cv_number * repeated_num, 100) # for repeatability
ctrl <- trainControl(method = "repeatedcv",
                     number = cv_number,
                     repeats = repeated_num,
                     verboseIter = FALSE,
                     savePredictions=TRUE,
                     classProbs=TRUE,
                     summaryFunction = twoClassSummary,
                     seeds=seeds_lst
                     )
# model_gbm <- readRDS(model_file)

# train different models
model_randomForest <- model_train(trainning_data[,stat_cols], trainning_data$Type, 'rf')
model_gbm <- model_train(trainning_data[,stat_cols], trainning_data$Type, 'gbm')
model_glmnet <- model_train(trainning_data[,stat_cols], trainning_data$Type, 'glmnet')
model_svmLinear <- model_train(trainning_data[,stat_cols], trainning_data$Type, 'svmLinear')
model_svmRadial <- model_train(trainning_data[,stat_cols], trainning_data$Type, 'svmRadial')

models_list <- list('rf'=model_randomForest, 'gbm'=model_gbm, 'glmnet'=model_glmnet, 'svmLinear'=model_svmLinear, 'svmRadial'=model_svmRadial)
saveRDS(models_list, file.path(output_dir, "models_list.rds"))
# models_list <- readRDS(file.path(output_dir, "models_list.rds"))


# for trainning data
print('The Predicted Confusion matrix for Test samples, total sample:')
sprintf("%d Samples for Training Data", dim(trainning_data)[1])
feature_importance <- c()
for (model_name in names(models_list)){
  print(model_name)
  model <- models_list[[model_name]]
  # models_list$model_name
  print(confusionMatrix(model, mode='everything', positive=cancer))
  print(varImp(model))
  importance <- as.matrix(varImp(model)$importance)
  feature_importance <- cbind(feature_importance, importance[,1])
}
colnames(feature_importance) <- c('randomForest', 'gbm', 'glmnet', 'svmLinear', 'svmRadial')
feature_importance <- as.data.frame(feature_importance)
feature_importance$total <- feature_importance$randomForest + feature_importance$gbm + feature_importance$glmnet + feature_importance$svmLinear + feature_importance$svmRadial
feature_importance <- feature_importance[order(feature_importance$total, decreasing = TRUE),]
write.csv(feature_importance,file.path(output_dir, paste0("train_data_feature_importance.csv")),row.names=TRUE, quote=FALSE)

# prediction output is probability
model_test_val_total(models_list, test_data, 'test')
# prediction output is class
model_test_group_total(models_list, test_data, 'test')


# prediction output is probability
model_test_val_total(models_list, val_data, 'val')
# prediction output is class
model_test_group_total(models_list, val_data, 'val')
stopCluster(cl)
