#!/users/miniconda3/envs/r_env/bin/Rscript
# ref: https://machinelearningmastery.com/feature-selection-with-the-caret-r-package/
# ref: https://www.machinelearningplus.com/machine-learning/feature-selection/
# ref: https://www.datacareer.ch/blog/ridge-and-lasso-in-r/
# ref : https://remiller1450.github.io/s230f19/caret1.html
library(GenomicRanges)
library(glmnet)
library(caret)
library(pROC)
library(Boruta)

library(tidyverse)
library(gam)
library(randomForest)
library(gbm)

script_dir <- dirname(sys.frame(1)$ofile)
source(file.path(script_dir, "script", "utils.R"), chdir = FALSE)
args <- commandArgs(trailingOnly = TRUE)
set.seed(77)

bruta_feature_selection = function(input_table, cols, output_dir){
    cols <- intersect(cols, colnames(input_table))
    output_file <- file.path(output_dir, paste0("boruta_feature_importance.csv"))
    boruta_res <- Boruta(input_table[,cols], input_table[['Type']], doTrace = 2)
    print(boruta_res)
    # set.seed(77) # if you want to the exactly same results in a script, you have to reset the seed
    # boruta_res2 <- Boruta(input_table[,stat_cols], input_table[['Type']], doTrace = 2)
    boruta_res_df <- as.data.frame(boruta_res$finalDecision)
    colnames(boruta_res_df) = 'finalDecision'
    boruta_res_df <- tibble::rownames_to_column(boruta_res_df, "feature")
    boruta_res_df <- boruta_res_df[(boruta_res_df$finalDecision == 'Confirmed' | boruta_res_df$finalDecision == 'Tentative'), ]
    boruta_res_df <- boruta_res_df %>% arrange(factor(finalDecision, levels = c('Confirmed', 'Tentative', 'Rejected')))
    write.csv(boruta_res_df, output_file,row.names=FALSE, quote=FALSE)
}


input_file <- args[1]   #input training data 
output_dir <- args[2]




sample_info_cols <- c('Sample_ID','Name','Sample_Type','Disease','Type','Sex','Age', 'sample_id')
qc_cols <- c('Raw_Yield(G)','Raw_Reads_Num(M)','Raw_Q30(%)','Raw_Q20(%)','Raw_GC(%)','Clean_Yield(G)','Clean_Reads_Num(M)','Clean_Q30(%)','Clean_Q20(%)','Clean_GC(%)','Effective(%)','Duplication_Rate(%)','Insert_Size_Peak(bp)','V_primer_percent','J_primer_percent','merge_rate','clean_rate','clone_reads_ratio','Clonereads','convert_clonotype')

stat_cols <- c('observedDiversity_mean','d50Index_mean','shannonWienerIndex_mean','normalizedShannonWienerIndex_mean','inverseSimpsonIndex_mean','group_freq_2','group_freq_3','group_freq_4','group_freq_5','group_freq_6','diversity')

VJ_cols <- c("TRBJ1-1","TRBJ1-2","TRBJ1-3","TRBJ1-4","TRBJ1-5","TRBJ1-6","TRBJ2-1","TRBJ2-2","TRBJ2-2P","TRBJ2-3","TRBJ2-4","TRBJ2-5","TRBJ2-6","TRBJ2-7","TRBV1","TRBV10-1","TRBV10-2","TRBV10-3","TRBV11-1","TRBV11-2","TRBV11-3","TRBV12-1","TRBV12-2","TRBV12-3","TRBV12-4","TRBV12-5","TRBV13","TRBV14","TRBV15","TRBV16","TRBV17","TRBV18","TRBV19","TRBV2","TRBV20-1","TRBV21-1","TRBV23-1","TRBV24-1","TRBV25-1","TRBV26","TRBV27","TRBV28","TRBV29-1","TRBV3-1","TRBV3-2","TRBV30","TRBV4-1","TRBV4-2","TRBV4-3","TRBV5-1","TRBV5-3","TRBV5-4","TRBV5-5","TRBV5-6","TRBV5-7","TRBV5-8","TRBV6-1","TRBV6-2","TRBV6-3","TRBV6-4","TRBV6-5","TRBV6-6","TRBV6-7","TRBV6-8","TRBV6-9","TRBV7-1","TRBV7-2","TRBV7-3","TRBV7-4","TRBV7-5","TRBV7-6","TRBV7-7","TRBV7-8","TRBV7-9","TRBV9")

anno_cols <- c("type_CMV","type_Dengue virus","type_EBV","type_HBV","type_HIV","type_HPV","type_Hepatitis C virus","type_Influenza A virus","type_Mycobacterium Tuberculosis","type_SARS-CoV1","type_SARS-CoV2","type_hRSV","type_lymphocytic choriomeningitis virus")
enriched_cols <- c("EGFR_type","KRAS_type","LUNG_CANCER_GDNA_type","LUNG_CANCER_TISSUE_type","EGFR_score","KRAS_score","LUNG_CANCER_GDNA_score","LUNG_CANCER_TISSUE_score", "REF_score", 'REF_type', 'Convergence_type','Convergence_freq','Convergence_score', 'KRAS_T_type','KRAS_T_score','EGFR_T_type','EGFR_T_score','BT_type','BT_score')
aa_len_cols <- c('func_raito','A_len','B_len','C_len','D_len','E_len','A_score','R_score','N_score','D_score','C_score','E_score','Q_score','G_score','H_score','I_score','L_score','K_score','M_score','F_score','P_score','S_score','T_score','W_score','Y_score','V_score')
clinical_info_cols <- c("somking","quit_somking","lung_cancer_history","family_history","COPD","left","right","up","down","middle","size(mm)","SN","pSN","GGN", "LN_num", "glitch")
n_insertion_cols <- c("vd_na","vd_zero_ratio","vd_zero_w_ratio","vd_mean","vd_w_mean","dj_na","dj_zero_ratio","dj_zero_w_ratio","dj_mean","dj_w_mean","vj_na","vj_zero_ratio","vj_zero_w_ratio","vj_mean","vj_w_mean")
input_table <- read_csv(input_file)
stat_cols <- c(stat_cols, VJ_cols, enriched_cols, aa_len_cols)

total_cols <- c(stat_cols, clinical_info_cols)

print('features for build model')
print(total_cols)
input_table <- qc_filter(input_table, total_cols)
# get balanced healthy and cancer samples
input_table <- rebalance_table(input_table, 3)

multiple_feature_group_lst <- list(tcr=stat_cols, clinical=clinical_info_cols, all=total_cols)
# multiple_feature_group_lst <- list(all=total_cols, clinical=clinical_info_cols, tcr=stat_cols)
for (key in names(multiple_feature_group_lst)){
  print(paste("Group Name:",key))
  print("Input features")
  print(multiple_feature_group_lst[[key]])
  group_output_dir <- paste0(output_dir, "/", key)
  dir.create(group_output_dir)
  print(paste("Outputdir:",group_output_dir))
  bruta_feature_selection(input_table, multiple_feature_group_lst[[key]], group_output_dir)
}