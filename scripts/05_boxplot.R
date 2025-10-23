library(ggplot2)
library(ggpubr)
library(here)

#Create folder to store boxplot figures
boxplot_dir <- here("results", "Boxplots")
if (!dir.exists(boxplot_dir)) dir.create(boxplot_dir)


#' Create and save a boxplot for a given gene and dataset
#'
#' This function generates a boxplot with jittered points and Wilconxon t-test 
#' significance for a specified gene in a given dataset, and saves the plot 
#' as a PNG file.
#'
#' @param gene Character. Name of the gene to plot.
#' @param dataset Data frame. Expression data with a column for disease labels
#'               (train and test datasets can be taken directly from 04_lasso.R)
#' @param dataset_name Character. Used in the output filename (e.g., "training_data").
#' @param boxplot_dir Character. Path to the folder where plots will be saved.
#'
#' @return Saves a PNG file of the boxplot and prints the plot to the device.
#'
#' @examples
#' plot_gene_boxplot("ZNF443", train, "training_data", "results/Boxplots")
plot_gene_boxplot <- function(gene, dataset, dataset_name, boxplot_dir) {
  expr_df <- as.data.frame(dataset[, gene])
  expr_df$group <- dataset$disease_label
  colnames(expr_df)[1] <- gene
  expr_df$group <- factor(expr_df$group, levels = c("HP", "IPF"))
  
  p <- ggplot(expr_df, aes(x = group, y = .data[[gene]], color = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.4, fill = NA, size = 1.2) +
    geom_jitter(width = 0.15, size = 2) +
    scale_color_manual(values = c("steelblue", "firebrick")) +
    labs(y = paste0(gene, " expression"), x = NULL) +
    theme_classic(base_size = 16) +
    theme(legend.position = "top") +
    stat_compare_means(
      comparisons = list(c("HP", "IPF")),
      method = "t.test",
      label = "p.format",
      label.y = max(expr_df[[gene]]) * 0.2
    )
  
  print(p)
  
  ggsave(file.path(boxplot_dir, paste0(gene, "_exp_boxplot_", dataset_name, ".png")),
         p, dpi = 300, height = 6, width = 8)
}

# List of genes
genes <- c("ZNF443", "BORCS6", "RNF208", "SDHAF1")

# Loop over genes for train and test datasets
for (gene in genes) {
  plot_gene_boxplot(gene, train, "training_data", boxplot_dir)
  plot_gene_boxplot(gene, test, "test_data", boxplot_dir)
}