## This script is to build and test the machine learning LASSO model with the prepared data
## There are 380 genes with consistent direction across 3 studies
## There are 261 samples in total : HP (n = 118) + IPF (n = 143)
## ML dataset is 261 samples x 380 genes for the models
## The platform corrected dataset is filtered using top 50 genes from mRMR method

# Set working directory
setwd("C:/IIT KGP/Method4/ML_v4")

# Install libraries
install.packages("glmnet")

# Import libraries
library(tidyverse)
library(glmnet)
library(ggplot2)
library(reshape2)
library(pROC)


###################################################
########### 70-30 split : Split 1 #################
###################################################


# Load the data
train_data_split1 <- read.csv("C:/IIT KGP/Method4/ML_v4/70-30 split/Train_Dataset_70.csv", row.names = 1)
test_data_split1 <- read.csv("C:/IIT KGP/Method4/ML_v4/70-30 split/Test_Dataset_30.csv", row.names = 1)

mRMR_genes <- c('BORCS6', 'KLK7', 'ZNF443', 'OLIG1', 'SDHAF1', 'GBP4', 'RNF208',
                'GPR25', 'C19orf47', 'DCANP1', 'AGER', 'TMEM102', 'MDP1', 'ALPP',
                'NKG7', 'ARGLU1', 'TIGD5', 'APOM', 'CSAG3', 'MROH6', 'ADAMTS1',
                'DPM3', 'MYRF', 'TMIGD2', 'SLC27A5', 'IRF4', 'ANXA2R', 'LGALS4',
                'SULT1E1', 'KDM5D', 'TNNC1', 'INPP5J', 'ZNF668', 'IGSF22', 'ITGAV',
                'KIF12', 'EREG', 'ATP6V1F', 'DDX3Y', 'COL4A1', 'RGS9BP', 'MTLN', 'CXCL8',
                'CLDN18', 'COL5A3', 'C20orf204', 'DNASE1L2', 'CLIC3', 'DDIT3', 'SLC5A2')

train <- train_data_split1[ ,colnames(train_data_split1) %in% mRMR_genes] #Filter the training set with only mRMR genes
train$disease_label <- train_data_split1$disease_label #Add the disease label

X_train <- train[,1:50] #182 NGS samples, 50 genes
y_train <- factor(train[,51], levels = c("IPF", "HP")) #182 NGS samples info separated
y_train_bin <- ifelse(y_train == "IPF", 0, 1)


test <- test_data_split1[, colnames(test_data_split1) %in% mRMR_genes]
test$disease_label <- test_data_split1$disease_label

X_test <- test[, 1:50] 
y_test <- factor(test[, 51], levels = c("IPF", "HP"))
y_test_bin  <- ifelse(y_test == "IPF", 0, 1)

# Scaling data
X_train_scaled <- scale(X_train)
X_test_scaled <- scale(X_test, center = attr(X_train_scaled, "scaled:center"),
                       scale = attr(X_train_scaled, "scaled:scale")) 
# X_test is scaled based on information from X_train to prevent data leakage


# ------------------------ #
# LASSO (glmnet)
# ------------------------ #

auc_list <- c()
acc_list <- c()
gene_counts <- 2:length(mRMR_genes)
selected_gene_sets <- list()

for (k in gene_counts) {
  genes_k <- mRMR_genes[1:k]
  
  X_train_sub <- X_train_scaled[, genes_k, drop = FALSE]
  X_test_sub  <- X_test_scaled[, genes_k, drop = FALSE]
  
  # LASSO with cross-validation
  cv_fit <- cv.glmnet(as.matrix(X_train_sub), y_train_bin, 
                      alpha = 1, family = "binomial",
                      type.measure = "auc", nfolds = 10)
  
  # Predict probabilities on test data
  pred_prob <- predict(cv_fit, newx = as.matrix(X_test_sub), 
                       s = "lambda.min", type = "response")
  
  # AUC
  roc_obj <- roc(y_test_bin, as.vector(pred_prob))
  auc_list[k] <- auc(roc_obj)
  
  # Accuracy
  pred_label <- ifelse(pred_prob > 0.5, 1, 0)
  acc <- mean(pred_label == y_test_bin)
  acc_list[k] <- acc
  
  # Save gene list
  selected_gene_sets[[k]] <- genes_k
}


# Combine results
performance_df <- data.frame(
  NumGenes = gene_counts,
  AUC = na.omit(auc_list),
  Accuracy = na.omit(acc_list)
)


# Optional: Plot performance

png(filename = "LASSO_AUC_and_Accuracy_plot.png",
    height = 6, width = 8,
    units = 'in', res = 300)
# Base plot: AUC
plot(performance_df$NumGenes, performance_df$AUC, type = "b", pch = 19, col = "blue",
     xlab = "Number of Genes", ylab = "AUC", ylim = c(0, 1),
     main = "LASSO: AUC and Accuracy vs. Number of Genes")
grid()

# Add second y-axis on right for Accuracy
par(new = TRUE)
plot(performance_df$NumGenes, performance_df$Accuracy, type = "b", pch = 17, col = "red",
     axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1))
axis(side = 4, col = "red", col.axis = "red")
mtext("Accuracy", side = 4, line = 3, col = "red")

# Add legend
legend("bottomright", legend = c("AUC", "Accuracy"),
       col = c("blue", "red"), pch = c(19, 17))

dev.off()

########################################################
#####  Based on the principle of minimum features, #####
#####  maximum performance, 8 genes were selected  #####
#####  having AUC of 0.92 and Accuracy of 0.89     #####
########################################################


# Extract the 8 genes
lasso_genes <- selected_gene_sets[[8]]



################# ROC-AUC on train data ##########################
auc_results_train <- data.frame(Gene = character(), AUC = numeric(), CI_lower = numeric(),
                          CI_upper = numeric(), stringsAsFactors = FALSE)

for (gene in lasso_genes) {
  roc_obj <- roc(y_train_bin, X_train[[gene]])
  auc_val <- as.numeric(auc(roc_obj))
  ci_vals <- ci.auc(roc_obj)
  auc_results_train <- rbind(auc_results_train, data.frame(Gene = gene, AUC = auc_val, CI_lower = ci_vals[1],
                                               CI_upper = ci_vals[3]))
}

# View sorted by AUC
auc_results_train <- auc_results_train[order(-auc_results_train$AUC), ]
print(auc_results_train)


png("Predictor_genes_roc_train_data.png", width = 10, height = 8, units = "in", res = 300)

par(mfrow = c(2, 4))  # 2 row, 4 columns, 8 plots

for (gene in head(auc_results_train$Gene, 8)) {
  roc_obj <- roc(y_train, X_train[[gene]])
  plot(roc_obj, main = paste(gene), col = "red",
       xlab = "False Positive Rate",    # relabel X-axis
       ylab = "True Positive Rate")     # relabel Y-axis
  legend("bottomright", legend = paste("AUC =", round(auc(roc_obj), 3)))
}

dev.off()

# Individual ROC plots for training data
for (gene in auc_results_train$Gene) {
  roc_obj <- roc(y_train, X_train[[gene]])
  
  auc_val <- as.numeric(auc(roc_obj))
  ci_vals <- ci.auc(roc_obj)
  ci_lower <- round(ci_vals[1], 3)
  ci_upper <- round(ci_vals[3], 3)
  
  png(filename = paste0("Train_ROC_", gene, ".png"),
      width = 8, height = 8, units = "in", res = 300)
  
  plot(roc_obj,
       xlab = "False Positive Rate",    # relabel X-axis
       ylab = "True Positive Rate",     # relabel Y-axis,
       main = gene, 
       col = "red",
       cex.main = 1.0,                  # bigger title font
       cex.lab = 1.4,                   # bigger axis label font
       cex.axis = 1.2)                  # bigger axis number font
  
  # Add bold AUC + 95% CI text
  
  text(0.6, 0.1, paste0("AUC = ", round(auc_val,3), 
                        " \n(95% CI: ", ci_lower, "-", ci_upper, ")"),
       cex = 1.2, font = 2, col = "black")
  
  dev.off()
}



################# ROC-AUC on test data ##########################
auc_results_test <- data.frame(Gene = character(), AUC = numeric(), CI_lower = numeric(),
                          CI_upper = numeric(), stringsAsFactors = FALSE)

for (gene in lasso_genes) {
  roc_obj <- roc(y_test_bin, X_test[[gene]])
  auc_val <- as.numeric(auc(roc_obj))
  ci_vals <- ci.auc(roc_obj)
  auc_results_test <- rbind(auc_results_test, data.frame(Gene = gene, AUC = auc_val, CI_lower = ci_vals[1],
                                               CI_upper = ci_vals[3]))
}

# View sorted by AUC
auc_results_test <- auc_results_test[order(-auc_results_test$AUC), ]
print(auc_results_test)


png("Predictor_genes_roc_test_data.png", width = 10, height = 8, units = "in", res = 300)

par(mfrow = c(2, 4))  # 2 row, 4 columns, 8 plots

for (gene in head(auc_results_test$Gene, 8)) {
  roc_obj <- roc(y_test, X_test[[gene]])
  plot(roc_obj, main = paste(gene), col = "red",
       xlab = "False Positive Rate",    # relabel X-axis
       ylab = "True Positive Rate")     # relabel Y-axis
  legend("bottomright", legend = paste("AUC =", round(auc(roc_obj), 3)))
}

dev.off()

# Individual ROC plots for test data
for (gene in auc_results_test$Gene) {
  roc_obj <- roc(y_test, X_test[[gene]])
  
  auc_val <- as.numeric(auc(roc_obj))
  ci_vals <- ci.auc(roc_obj)
  ci_lower <- round(ci_vals[1], 3)
  ci_upper <- round(ci_vals[3], 3)
  
  png(filename = paste0("Test_ROC_", gene, ".png"),
      width = 8, height = 8, units = "in", res = 300)
  
  plot(roc_obj,
       xlab = "False Positive Rate",    # relabel X-axis
       ylab = "True Positive Rate",     # relabel Y-axis,
       main = gene, 
       col = "red",
       cex.main = 1.0,                  # bigger title font
       cex.lab = 1.4,                   # bigger axis label font
       cex.axis = 1.2)                  # bigger axis number font
  
  # Add bold AUC + 95% CI text
  
  text(0.6, 0.1, paste0("AUC = ", round(auc_val,3), 
                        " \n(95% CI: ", ci_lower, "-", ci_upper, ")"),
       cex = 1.2, font = 2, col = "black")
  
  dev.off()
}


##################################################
############## Gene combinations #################
##################################################

# Store results
combo_results <- data.frame(
  Combination = character(),
  NumGenes = integer(),
  AUC = numeric(),
  Accuracy = numeric(),
  stringsAsFactors = FALSE
)

# Function to evaluate a gene combination
evaluate_combination <- function(gene_combo) {
  X_train_sub <- X_train_scaled[, gene_combo, drop = FALSE]
  X_test_sub  <- X_test_scaled[, gene_combo, drop = FALSE]
  
  # LASSO with cross-validation
  cv_fit <- cv.glmnet(as.matrix(X_train_sub), y_train_bin,
                      alpha = 1, family = "binomial",
                      type.measure = "auc", nfolds = 10)
  
  # Predict probabilities on test data
  pred_prob <- predict(cv_fit, newx = as.matrix(X_test_sub),
                       s = "lambda.min", type = "response")
  
  # ROC-AUC
  roc_obj <- roc(y_test_bin, as.vector(pred_prob))
  auc_val <- as.numeric(auc(roc_obj))
  
  # Accuracy
  pred_label <- ifelse(pred_prob > 0.5, 1, 0)
  acc <- mean(pred_label == y_test_bin)
  
  return(list(auc = auc_val, acc = acc, roc = roc_obj, genes = gene_combo))
}

# Store top combinations with ROC for later plotting
top_combos_list <- list()

# Evaluate combinations
for (k in 2:5) {
  combos <- combn(lasso_genes, k, simplify = FALSE)
  
  for (combo in combos) {
    combo_name <- paste(combo, collapse = ", ")  # consistent combo name
    
    res <- evaluate_combination(combo)
    
    # Save combo info
    combo_results <- rbind(combo_results, data.frame(
      Combination = combo_name,
      NumGenes = k,
      AUC = res$auc,
      Accuracy = res$acc,
      stringsAsFactors = FALSE
    ))
    
    # Save ROC with exact matching name
    top_combos_list[[combo_name]] <- res$roc
  }
}

# Sort by AUC descending
combo_results <- combo_results[order(-combo_results$AUC), ]

# View results
print(head(combo_results))

# Save to CSV if needed
write.csv(combo_results, "LASSO_combinations_AUC_results.csv", row.names = FALSE)


# ======================================
# Save ROC plots for top 10 combinations
# ======================================

# Create output folder (optional)
output_dir <- "Combination_Top10_ROC_Plots"
if (!dir.exists(output_dir)) dir.create(output_dir)

top_10 <- head(combo_results, 10)

#Plot
for (i in 1:nrow(top_10)) {
  combo_name <- top_10$Combination[i]
  roc_obj <- top_combos_list[[combo_name]]
  
  auc_val <- round(top_10$AUC[i], 3)
  
  # Calculate 95% CI
  ci_vals <- ci.auc(roc_obj)
  ci_lower <- round(ci_vals[1], 3)
  ci_upper <- round(ci_vals[3], 3)
  
  # Clean file name (remove illegal characters)
  file_safe_name <- gsub("[^a-zA-Z0-9_]", "_", combo_name)
  file_name <- file.path(output_dir, paste0("ROC_", i, "_", file_safe_name, ".png"))
  
  # Save PNG
  png(filename = file_name, width = 8, height = 8, units = "in", res = 300)  # 300 dpi
  
  # Plot ROC
  plot(roc_obj,
       col = "red", lwd = 3,            # thick red ROC line
       xlab = "False Positive Rate",    # relabel X-axis
       ylab = "True Positive Rate",     # relabel Y-axis
       main = paste0("Combination ", i, ":\n", combo_name),
       cex.main = 1.0,                  # bigger title font
       cex.lab = 1.4,                   # bigger axis label font
       cex.axis = 1.2)                  # bigger axis number font
  
  # Add bold AUC + 95% CI text
  text(0.6, 0.1, paste0("AUC = ", auc_val, 
                        " (95% CI: ", ci_lower, "-", ci_upper, ")"),
       cex = 1.3, font = 2, col = "black")
  
  dev.off()
}









###################################################
########### 60-40 split : Split 2 #################
###################################################


# Load the data
train_data_split2 <- read.csv("Train_Dataset_60.csv", row.names = 1)
test_data_split2 <- read.csv("Test_Dataset_40.csv", row.names = 1)

mRMR_genes <- c('RNF208', 'KRT6A', 'SLC27A5', 'BORCS6', 'TDO2', 'TMEM102', 'ARGLU1', 
                'SDHAF1', 'TMIGD2', 'PSPN', 'ZNF443', 'OLIG1', 'DPM3', 'GBP4', 'MROH6', 
                'INPP5J', 'ZNF668', 'C19orf47', 'DNASE1L2', 'ANXA2R', 'AGER', 'ZGLP1', 
                'TIGD5', 'COL5A3', 'DCANP1', 'NOXO1', 'TMEM132D', 'NKG7', 'PPFIA3', 
                'SLC16A11', 'IGSF22', 'APOM', 'IRF4', 'C20orf204', 'MDP1', 'HLA.B',  #HLA-B is HLA.B in R
                'COL4A1', 'ATP6V1F', 'MTLN', 'LGALS4', 'KDM5D', 'GPR25', 'ITGAV', 
                'IER5L', 'KIF12', 'SULT1E1', 'EXOC3L1', 'NKX6.2', 'C22orf39', 'MYRF') #Same for NKX6-2

train <- train_data_split2[ ,colnames(train_data_split2) %in% mRMR_genes] #Filter the training set with only mRMR genes
train$disease_label <- train_data_split2$disease_label #Add the disease label

X_train <- train[,1:50] #182 NGS samples, 50 genes
y_train <- factor(train[,51], levels = c("IPF", "HP")) #182 NGS samples info separated
y_train_bin <- ifelse(y_train == "IPF", 0, 1)


test <- test_data_split2[, colnames(test_data_split2) %in% mRMR_genes]
test$disease_label <- test_data_split2$disease_label

X_test <- test[, 1:50] 
y_test <- factor(test[, 51], levels = c("IPF", "HP"))
y_test_bin  <- ifelse(y_test == "IPF", 0, 1)

# Scaling data
X_train_scaled <- scale(X_train)
X_test_scaled <- scale(X_test, center = attr(X_train_scaled, "scaled:center"),
                       scale = attr(X_train_scaled, "scaled:scale")) 
# X_test is scaled based on information from X_train to prevent data leakage


# ------------------------ #
# LASSO (glmnet)
# ------------------------ #

auc_list <- c()
acc_list <- c()
gene_counts <- 2:length(mRMR_genes)
selected_gene_sets <- list()

for (k in gene_counts) {
  genes_k <- mRMR_genes[1:k]
  
  X_train_sub <- X_train_scaled[, genes_k, drop = FALSE]
  X_test_sub  <- X_test_scaled[, genes_k, drop = FALSE]
  
  # LASSO with cross-validation
  cv_fit <- cv.glmnet(as.matrix(X_train_sub), y_train_bin, 
                      alpha = 1, family = "binomial",
                      type.measure = "auc", nfolds = 10)
  
  # Predict probabilities on test data
  pred_prob <- predict(cv_fit, newx = as.matrix(X_test_sub), 
                       s = "lambda.min", type = "response")
  
  # AUC
  roc_obj <- roc(y_test_bin, as.vector(pred_prob))
  auc_list[k] <- auc(roc_obj)
  
  # Accuracy
  pred_label <- ifelse(pred_prob > 0.5, 1, 0)
  acc <- mean(pred_label == y_test_bin)
  acc_list[k] <- acc
  
  # Save gene list
  selected_gene_sets[[k]] <- genes_k
}


# Combine results
performance_df <- data.frame(
  NumGenes = gene_counts,
  AUC = na.omit(auc_list),
  Accuracy = na.omit(acc_list)
)


# Optional: Plot performance

png(filename = "LASSO_AUC_and_Accuracy_plot_2.png",
    height = 6, width = 8,
    units = 'in', res = 300)
# Base plot: AUC
plot(performance_df$NumGenes, performance_df$AUC, type = "b", pch = 19, col = "blue",
     xlab = "Number of Genes", ylab = "AUC", ylim = c(0, 1),
     main = "LASSO: AUC and Accuracy vs. Number of Genes")
grid()

# Add second y-axis on right for Accuracy
par(new = TRUE)
plot(performance_df$NumGenes, performance_df$Accuracy, type = "b", pch = 17, col = "red",
     axes = FALSE, xlab = "", ylab = "", ylim = c(0, 1))
axis(side = 4, col = "red", col.axis = "red")
mtext("Accuracy", side = 4, line = 3, col = "red")

# Add legend
legend("bottomright", legend = c("AUC", "Accuracy"),
       col = c("blue", "red"), pch = c(19, 17))

dev.off()

########################################################
#####  Based on the principle of minimum features, #####
#####  maximum performance, 30 genes were selected #####
#####  having AUC of 0.86 and Accuracy of 0.81     #####
########################################################


# Extract the 8 genes
lasso_genes <- selected_gene_sets[[29]]



# ROC-AUC on training data
auc_results_train_2 <- data.frame(Gene = character(), AUC = numeric(), CI_lower = numeric(),
                                CI_upper = numeric(), stringsAsFactors = FALSE)

for (gene in lasso_genes) {
  roc_obj <- roc(y_train_bin, X_train[[gene]])
  auc_val <- as.numeric(auc(roc_obj))
  ci_vals <- ci.auc(roc_obj)
  auc_results_train_2 <- rbind(auc_results_train_2, data.frame(Gene = gene, AUC = auc_val, CI_lower = ci_vals[1],
                                                           CI_upper = ci_vals[3]))
}

# View sorted by AUC
auc_results_train_2 <- auc_results_train_2[order(-auc_results_train_2$AUC), ]
print(auc_results_train_2)


png("Predictor_genes_roc_train_data_2.png", width = 12, height = 10, units = "in", res = 300)

par(mfrow = c(5, 6))  # 2 row, 4 columns, 8 plots

for (gene in head(auc_results_train_2$Gene, 30)) {
  roc_obj <- roc(y_train, X_train[[gene]])
  plot(roc_obj, main = paste(gene), col = "red")
  legend("bottomright", legend = paste("AUC =", round(auc(roc_obj), 3)))
}

dev.off()



################# ROC-AUC on test data ##########################
auc_results_test_2 <- data.frame(Gene = character(), AUC = numeric(), CI_lower = numeric(),
                               CI_upper = numeric(), stringsAsFactors = FALSE)

for (gene in lasso_genes) {
  roc_obj <- roc(y_test_bin, X_test[[gene]])
  auc_val <- as.numeric(auc(roc_obj))
  ci_vals <- ci.auc(roc_obj)
  auc_results_test_2 <- rbind(auc_results_test_2, data.frame(Gene = gene, AUC = auc_val, CI_lower = ci_vals[1],
                                                         CI_upper = ci_vals[3]))
}

# View sorted by AUC
auc_results_test_2 <- auc_results_test_2[order(-auc_results_test_2$AUC), ]
print(auc_results_test_2)


png("Predictor_genes_roc_test_data_2.png", width = 12, height = 10, units = "in", res = 300)

par(mfrow = c(5, 6))  # 2 row, 4 columns, 7 plots

for (gene in head(auc_results_test_2$Gene, 30)) {
  roc_obj <- roc(y_test, X_test[[gene]])
  plot(roc_obj, main = paste(gene), col = "red")
  legend("bottomright", legend = paste("AUC =", round(auc(roc_obj), 3)))
}

dev.off()

