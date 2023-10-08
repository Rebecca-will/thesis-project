library(limma)




Clara_T3$cluster <- ifelse(Clara_T3$cluster == "#1F77B4", 1,
                        ifelse(Clara_T3$cluster == "#2CA02C", 2,
                               ifelse(Clara_T3$cluster == '#FF7F0E',3,
                                      NA)))

Clara_T3$cluster <- as.numeric(Clara_T3$cluster)
cluster_labels <- Clara_T3$cluster
cluster_labels <- as.factor(cluster_labels)


colnames(Clara_T3)
expression_matrix <- subset(Clara_T3, select = c(3:22,24))
name<- colnames(expression_matrix)
cluster_1 <- expression_matrix$cluster

expression_matrix <- t(expression_matrix)

design<- model.matrix(~cluster_labels-1)

levels(cluster_labels)


fit <- lmFit(expression_matrix, design)
fit <- eBayes(fit)

logFC <- fit$coefficients[,'cluster_labels2']
pvalue <- fit$p.value[,'cluster_labels2']


DE_results <- data.frame(Gene=name, LogFC=logFC, Pvalue = pvalue, cluster = cluster, row.names = row.names(expression_matrix))

DE_results$DES <- -log10(DE_results$Pvalue) * abs(DE_results$LogFC)
DE_results$cluster <- factor(DE_results$cluster)
DE_results$DES <- factor(DE_results$DES)

# Create a volcano plot
volcano_plot <- ggplot(DE_results, aes(x = LogFC, y = -log10(Pvalue), color = cluster, size = DES)) +
  geom_point() +
  scale_color_manual(values = c("#1F77B4", "#FF7F0E", "#2CA02C")) +
  labs(x = "Log Fold Change", y = "-log10(P-value)", color = "cluster", size = "DES") +
  ggtitle("Volcano Plot of Differential Expression") +
  theme_bw()

# Print the volcano plot
print(volcano_plot)


str(DE_results)

