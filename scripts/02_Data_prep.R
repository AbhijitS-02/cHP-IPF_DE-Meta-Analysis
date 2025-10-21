## This script is to prepare the data for the machine learning models
## There are 380 genes with consistent direction across three studies

## There are 261 samples in total : HP (n = 118) + IPF (n = 143)
## ML dataset should be 261 samples x 380 genes for LASSO model


# Import libraries
library(tidyverse)
library(here)
library(ggplot2)
library(reshape2)
library(sva) # sva library has the Combat() function required for batch correction

### Load the expression/count matrices for all the datasets

#Load normalized gene count data of GSE184316 
## vst() normalization done using DESeq2
count.GSE184316 <- read.csv(here("data", "GSE184316", "Normalized_count_HP_IPF_GSE184316.csv"),
                            header = TRUE) %>%
  column_to_rownames(., var = "X") # Convert the first column "X" (gene symbols) to row names


# Load normalized gene count data of GSE150910 (includes both biopsy and explant samples)
count.GSE150910 <- read.csv(here("data", "GSE150910", "Normalized_count_HP_IPF_GSE150910.csv"),
                            header = TRUE) %>%
  column_to_rownames(., var = "X") # Convert the first column "X" (gene symbols) to row names



### Load and prepare the sample information data for all the datasets

## Sample info GSE184316: HP (n = 36) and IPF (n = 40)
sample_info.GSE184316 <- read.csv(here("data", "GSE184316", "Sample_information_HP_IPF_GSE184316.csv"),
                                  header = TRUE)

sample_info.GSE184316 <- sample_info.GSE184316[, c(2,3)] #remove the first column named "X"
colnames(sample_info.GSE184316)[2] = "disease" #rename condition column to disease


## Sample info GSE150910: Explant (HP: n = 56, IPF: n = 67) and Biopsy (HP: n = 26, IPF: n = 36)
sample_info.GSE150910 <- read.csv(here("data", "GSE150910", "Sample_information_HP_IPF_GSE150910.csv"),
                                  header = TRUE)

sample_info.GSE150910 <- sample_info.GSE150910[, c(2,3)] #remove the first column named "X"
colnames(sample_info.GSE150910)[2] = "disease" #rename condition column to disease

sample_info.GSE150910$disease[sample_info.GSE150910$disease == "CHP"] <- "HP"  #Rename CHP as HP for proper labelling



## Precautionary check to see whether the sample names and count are same in counts matrix and sample info

#GSE184316
all(sample_info.GSE184316$sample %in% colnames(count.GSE184316))#if TRUE, all the samples are present in both matrices

#GSE150910
all(sample_info.GSE150910$sample %in% colnames(count.GSE150910))#if TRUE, all the samples are present in both matrices



#Load the FINM results with the 380 genes having consistent direction in the 2 NGS datasets
common_genes <- read.csv(here("results", "Common_genes_same_direction_expr.csv"), header = TRUE)



## Filter the counts datasets with only the 380 genes

#GSE184316
count.GSE184316 <- count.GSE184316[rownames(count.GSE184316) %in% common_genes$Gene_symbol, ]
write.csv(count.GSE184316, here("results", "counts_GSE184316_HPvIPF_380genes.csv")) #optional

#GSE150910
count.GSE150910 <- count.GSE150910[rownames(count.GSE150910) %in% common_genes$Gene_symbol, ]
write.csv(count.GSE150910, here("results", "counts_GSE150910_HPvIPF_380genes.csv")) #optional



# In order to merge the 2 datasets, the rownames, i.e., gene symbols, must be in order for cbind() to work
# Reorder count.GSE150910 to match the rownames of count.GSE184316
count.GSE150910 <- count.GSE150910[rownames(count.GSE184316), ] 


## Optional - to reconfirm if they are in same order
all(rownames(count.GSE150910) == rownames(count.GSE184316)) #TRUE, if same order


### Merge the 3 count/expr datasets column wise, matching their corresponding gene symbols (rownames)
counts_data_combined <- cbind(count.GSE184316, count.GSE150910)


# Merge the sample information of the 2 NGS datasets
sample_info_combined <- rbind(sample_info.GSE184316, sample_info.GSE150910)


## Match the order of sample names in coluumn of expr_data_combined with the sample column of sample_info_combined
## This is done to ensure proper sample labelling required for machine learning dataset

counts_data_combined <- counts_data_combined[, sample_info_combined$sample]

all(colnames(counts_data_combined) == sample_info_combined$sample) #TRUE, if order matches


### Prepare the data for machine learning

# Step 1) Transpose the combined expression data as samples are in rows and genes are columns/features
counts_data_combined <- as.data.frame(t(counts_data_combined))

# Step 2) Add the disease column at the end as label (this will be separated again later in Python script)
data_final <- cbind(counts_data_combined, sample_info_combined$disease)

sum(is.na(data_final)) # Should return 0

colnames(data_final)[381] <- "disease_label"

write.csv(data_final, here("data", "Data_final_v4.csv"))

## HP (n = 118) + IPF (n = 143) = 261 samples x 380 genes




### PCA plot to check platform effects


# Only get the numerical values, excluding the last column i.e, disease labels
expr_t <- data_final[,1:380]

# Run PCA
pca_res <- prcomp(expr_t, scale. = TRUE)

# Get PC scores
pca_df <- as.data.frame(pca_res$x)

# calculate the explained variance of each principal component
explained_var <- round(100 * (pca_res$sdev^2) / sum(pca_res$sdev^2), 1)

# Labels to display the percentage of variance explained by each principal component
x_lab <- paste0("PC1 (", explained_var[1], "%)")
y_lab <- paste0("PC2 (", explained_var[2], "%)")

# Define the requied platform vector with information about the platform of the samples
platform_vector <- c(
  rep("illumina", 76),
  rep("iontorrent", 185))

pca_df$Platform <- platform_vector  # Add platform label
pca_df$disease <- data_final$disease_label #add disease label

# Plot PC1 vs PC2 of platform
pca_platform <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Platform)) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, 
               aes(fill = Platform), alpha = 0.25, color = NA) +  # <-- map fill
  geom_point(size = 3, alpha = 0.9) +
  theme_bw() +
  labs(title = "PCA of Expression Data by Platform",
       x = x_lab, y = y_lab) +
  theme(legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.background = element_blank())

print(pca_platform)
ggsave(filename = here("results", "PCA_(before_combat)_Platform.png"), pca_platform, dpi = 300, height = 6, width = 8)

# Plot PC1 vs PC2 of disease
pca_disease <- ggplot(pca_df, aes(x = PC1, y = PC2, color = disease)) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, 
               aes(fill = disease), alpha = 0.25, color = NA) + # 95% CI zone
  geom_point(size = 3, alpha = 0.8) +
  theme_bw() +
  labs(title = "PCA of Expression Data by Disease",
       x = x_lab, y = y_lab) +
  theme(legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.background = element_blank())

print(pca_disease)
ggsave(filename = here("results", "PCA_(before_combat)_Disease.png"), pca_disease, dpi = 300, height = 6, width = 8)





### Platform effect correction


# Create pheno data
pheno <- data.frame(disease = data_final$disease_label,
                    platform = platform_vector,
                    row.names = rownames(expr_t))

# Input matrix must be genes × samples (as original)
expr_matrix <- t(expr_t) #we are transposing expr_t as it was subset from data_final that has samples as rows and genes as columns

pheno$disease <- factor(pheno$disease, levels = c("IPF", "HP"))
# define design matrix to preserve the disease variable
modcombat <- model.matrix(~ disease, data = pheno)

# Apply ComBat
combat_expr <- ComBat(dat = as.matrix(expr_matrix),
                      batch = pheno$platform, # Platform to remove platform effect
                      mod = modcombat, # Design matrix to preserve biological disease variable
                      par.prior = TRUE,  # Empirical Bayes (default)
                      prior.plots = FALSE)


#Get the expression data after platform correction
combat_expr <- as.data.frame(t(combat_expr)) #convert genes to columns and samples as rows
combat_expr$disease_label <- data_final$disease_label #add disease label

write.csv(combat_expr, here("data", "Data_Final_v4_platformcorrected.csv")) #export



#subset the numerical columns to calculate principal components
combat_t <- combat_expr[,1:380]
# Run PCA
pca_res2 <- prcomp(combat_t, scale. = TRUE)

pca_df2 <- as.data.frame(pca_res2$x)
pca_df2$Platform <- platform_vector #Add platform
pca_df2$disease <- combat_expr$disease_label # add disease label

# Calculate variance explained by each principal component
explained_var <- round(100 * (pca_res2$sdev^2) / sum(pca_res2$sdev^2), 1)

# Display the explained variance in x and y labels
x_lab <- paste0("PC1 (", explained_var[1], "%)")
y_lab <- paste0("PC2 (", explained_var[2], "%)")

# PLot PC1 vs PC2 after platform correction - Platform
pca_platform <- ggplot(pca_df2, aes(x = PC1, y = PC2, fill = Platform)) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, alpha = 0.25) + # 95% CI zone
  geom_point(aes(color = Platform), size = 3, alpha = 0.9) +
  theme_bw() +
  labs(title = "PCA After ComBat: Platform Bias Correction",
       x = x_lab, y = y_lab) +
  theme(legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.background = element_blank())

print(pca_platform)
ggsave(filename = here("results", "PCA_(after_combat)_Platform.png"), pca_platform, dpi = 300, height = 6, width = 8)

# Plot PC1 vs PC2 after platform correction - separation based on disease
pca_disease <- ggplot(pca_df2, aes(x = PC1, y = PC2, fill = disease)) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.95, alpha = 0.25) +
  geom_point(aes(color = disease), size = 3, alpha = 0.9) +
  theme_bw() +
  labs(title = "PCA After ComBat: Disease based separation",
       x = x_lab, y = y_lab) +
  theme(legend.position = "right",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.background = element_blank())

print(pca_disease)
ggsave(filename = here("results", "PCA_(after_combat)_Disease.png"), pca_disease, dpi = 300, height = 6, width = 8)

# combat_expr is then used for machine learning model building
# END
