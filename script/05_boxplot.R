library(ggpubr)

##### 70-30 split

top_predictor <- "ZNF443"



train_expr <- as.data.frame(train[, top_predictor])
train_expr$group <- train$disease_label
colnames(train_expr)[1] <- "ZNF443"

test_expr <- as.data.frame(test[, top_predictor])
test_expr$group <- test$disease_label
colnames(test_expr)[1] <- "ZNF443"


######################### Expr in train data ###############################

train_expr$group <- factor(train_expr$group, levels = c("HP", "IPF"))

# PLot the bar plots with significance levels
box_train <- ggplot(train_expr, aes(x = group, y = ZNF443, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, fill = NA, size = 1.2) +   # Box outline
  geom_jitter(width = 0.15, size = 2) +                                    # Dots
  scale_color_manual(values = c("steelblue", "firebrick")) +               # Custom colors
  labs(y = "ZNF443 expression", x = NULL) +
  theme_classic(base_size = 16) +
  theme(legend.position = "top")
# Add significance annotation (t-test)
box_train <- box_train + stat_compare_means(
  comparisons = list(c("HP", "IPF")),
  method = "t.test",
  label = "p.format",
  label.y = max(train_expr$ZNF443) + 0.2
)

print(box_train)
ggsave("./results/70-30_split/ZNF443_Gene exp barplot_training_data.png", box_train, dpi = 300, height = 6, width = 8)




######################### Expr in test data ###############################

test_expr$group <- factor(test_expr$group, levels = c("HP", "IPF"))

# PLot the bar plots with significance levels
box_test <- ggplot(test_expr, aes(x = group, y = ZNF443, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, fill = NA, size = 1.2) +   # Box outline
  geom_jitter(width = 0.15, size = 2) +                                    # Dots
  scale_color_manual(values = c("steelblue", "firebrick")) +               # Custom colors
  labs(y = "ZNF443 expression", x = NULL) +
  theme_classic(base_size = 16) +
  theme(legend.position = "top")
# Add significance annotation (t-test)
box_test <- box_test + stat_compare_means(
  comparisons = list(c("HP", "IPF")),
  method = "t.test",
  label = "p.format",
  label.y = max(test_expr$ZNF443) + 0.2
)

print(box_test)
ggsave("./results/70-30_split/ZNF443_Gene exp barplot_test_data.png", box_test, dpi = 300, height = 6, width = 8)







top_predictor <- "BORCS6"



train_expr <- as.data.frame(train[, top_predictor])
train_expr$group <- train$disease_label
colnames(train_expr)[1] <- "BORCS6"

test_expr <- as.data.frame(test[, top_predictor])
test_expr$group <- test$disease_label
colnames(test_expr)[1] <- "BORCS6"


######################### Expr in train data ###############################

train_expr$group <- factor(train_expr$group, levels = c("HP", "IPF"))

# PLot the bar plots with significance levels
box_train <- ggplot(train_expr, aes(x = group, y = BORCS6, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, fill = NA, size = 1.2) +   # Box outline
  geom_jitter(width = 0.15, size = 2) +                                    # Dots
  scale_color_manual(values = c("steelblue", "firebrick")) +               # Custom colors
  labs(y = "BORCS6 expression", x = NULL) +
  theme_classic(base_size = 16) +
  theme(legend.position = "top")
# Add significance annotation (t-test)
box_train <- box_train + stat_compare_means(
  comparisons = list(c("HP", "IPF")),
  method = "t.test",
  label = "p.format",
  label.y = max(train_expr$BORCS6) + 0.2
)

print(box_train)
ggsave("./results/70-30_split/BORCS6_Gene exp barplot_training_data.png", box_train, dpi = 300, height = 6, width = 8)




######################### Expr in test data ###############################

test_expr$group <- factor(test_expr$group, levels = c("HP", "IPF"))

# PLot the bar plots with significance levels
box_test <- ggplot(test_expr, aes(x = group, y = BORCS6, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, fill = NA, size = 1.2) +   # Box outline
  geom_jitter(width = 0.15, size = 2) +                                    # Dots
  scale_color_manual(values = c("steelblue", "firebrick")) +               # Custom colors
  labs(y = "BORCS6 expression", x = NULL) +
  theme_classic(base_size = 16) +
  theme(legend.position = "top")
# Add significance annotation (t-test)
box_test <- box_test + stat_compare_means(
  comparisons = list(c("HP", "IPF")),
  method = "t.test",
  label = "p.format",
  label.y = max(test_expr$BORCS6) + 0.2
)

print(box_test)
ggsave("./results/70-30_split/BORCS6 exp barplot_test_data.png", box_test, dpi = 300, height = 6, width = 8)







top_predictor <- "RNF208"



train_expr <- as.data.frame(train[, top_predictor])
train_expr$group <- train$disease_label
colnames(train_expr)[1] <- "RNF208"

test_expr <- as.data.frame(test[, top_predictor])
test_expr$group <- test$disease_label
colnames(test_expr)[1] <- "RNF208"


######################### Expr in train data ###############################

train_expr$group <- factor(train_expr$group, levels = c("HP", "IPF"))

# PLot the bar plots with significance levels
box_train <- ggplot(train_expr, aes(x = group, y = RNF208, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, fill = NA, size = 1.2) +   # Box outline
  geom_jitter(width = 0.15, size = 2) +                                    # Dots
  scale_color_manual(values = c("steelblue", "firebrick")) +               # Custom colors
  labs(y = "RNF208 expression", x = NULL) +
  theme_classic(base_size = 16) +
  theme(legend.position = "top")
# Add significance annotation (t-test)
box_train <- box_train + stat_compare_means(
  comparisons = list(c("HP", "IPF")),
  method = "t.test",
  label = "p.format",
  label.y = max(train_expr$RNF208) + 0.2
)

print(box_train)
ggsave("./results/70-30_split/RNF208_Gene exp barplot_training_data.png", box_train, dpi = 300, height = 6, width = 8)




######################### Expr in test data ###############################

test_expr$group <- factor(test_expr$group, levels = c("HP", "IPF"))

# PLot the bar plots with significance levels
box_test <- ggplot(test_expr, aes(x = group, y = RNF208, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, fill = NA, size = 1.2) +   # Box outline
  geom_jitter(width = 0.15, size = 2) +                                    # Dots
  scale_color_manual(values = c("steelblue", "firebrick")) +               # Custom colors
  labs(y = "RNF208 expression", x = NULL) +
  theme_classic(base_size = 16) +
  theme(legend.position = "top")
# Add significance annotation (t-test)
box_test <- box_test + stat_compare_means(
  comparisons = list(c("HP", "IPF")),
  method = "t.test",
  label = "p.format",
  label.y = max(test_expr$RNF208) + 0.2
)

print(box_test)
ggsave("./results/70-30_split/RNF208 exp barplot_test_data.png", box_test, dpi = 300, height = 6, width = 8)





top_predictor <- "SDHAF1"



train_expr <- as.data.frame(train[, top_predictor])
train_expr$group <- train$disease_label
colnames(train_expr)[1] <- "SDHAF1"

test_expr <- as.data.frame(test[, top_predictor])
test_expr$group <- test$disease_label
colnames(test_expr)[1] <- "SDHAF1"


######################### Expr in train data ###############################

train_expr$group <- factor(train_expr$group, levels = c("HP", "IPF"))

# PLot the bar plots with significance levels
box_train <- ggplot(train_expr, aes(x = group, y = SDHAF1, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, fill = NA, size = 1.2) +   # Box outline
  geom_jitter(width = 0.15, size = 2) +                                    # Dots
  scale_color_manual(values = c("steelblue", "firebrick")) +               # Custom colors
  labs(y = "SDHAF1 expression", x = NULL) +
  theme_classic(base_size = 16) +
  theme(legend.position = "top")
# Add significance annotation (t-test)
box_train <- box_train + stat_compare_means(
  comparisons = list(c("HP", "IPF")),
  method = "t.test",
  label = "p.format",
  label.y = max(train_expr$SDHAF1) + 0.2
)

print(box_train)
ggsave("./results/70-30_split/SDHAF1_Gene exp barplot_training_data.png", box_train, dpi = 300, height = 6, width = 8)




######################### Expr in test data ###############################

test_expr$group <- factor(test_expr$group, levels = c("HP", "IPF"))

# PLot the bar plots with significance levels
box_test <- ggplot(test_expr, aes(x = group, y = SDHAF1, color = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, fill = NA, size = 1.2) +   # Box outline
  geom_jitter(width = 0.15, size = 2) +                                    # Dots
  scale_color_manual(values = c("steelblue", "firebrick")) +               # Custom colors
  labs(y = "SDHAF1 expression", x = NULL) +
  theme_classic(base_size = 16) +
  theme(legend.position = "top")
# Add significance annotation (t-test)
box_test <- box_test + stat_compare_means(
  comparisons = list(c("HP", "IPF")),
  method = "t.test",
  label = "p.format",
  label.y = max(test_expr$SDHAF1) + 0.2
)

print(box_test)
ggsave("./results/70-30_split/SDHAF1 exp barplot_test_data.png", box_test, dpi = 300, height = 6, width = 8)