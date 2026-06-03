# Figure 6: Visualization of the ODK-means induced topology properties on the 
# Old Faithful dataset using calculated alpha* significance levels.

library(ggplot2)
source("ODKm_pub.R")
source("alpha_star.R")

set.seed(123)
X <- scale(faithful)
res <- odkm(data = X, K = 2, restart = 50, alpha = 0.1)

U_matrix <- res$U 
cluster_vector <- apply(U_matrix, 1, which.max)

color_cluster1 <- "#440154" 
color_cluster2 <- "#21908C" 

df_plot <- data.frame(
  eruptions  = X[, 1],
  waiting    = X[, 2],
  Unit_ID    = 1:nrow(X),
  Cluster    = as.factor(cluster_vector),
  Is_Outlier = FALSE,
  Alpha_Star = "Inlier"
)

outlier_labels <- res$out

if (length(outlier_labels) > 0) {
  df_plot$Is_Outlier[df_plot$Unit_ID %in% outlier_labels] <- TRUE
  
  alpha_vector <- numeric(length(outlier_labels))
  for (i in 1:length(outlier_labels)) {
    alpha_vector[i] <- evaluate_alpha(data = X, number = outlier_labels[i], U = U_matrix, pos = FALSE)
  }
  
  alpha_rounded <- round(alpha_vector, 4)
  
  for (i in 1:length(outlier_labels)) {
    unit <- outlier_labels[i]
    df_plot$Alpha_Star[df_plot$Unit_ID == unit] <- as.character(alpha_rounded[i])
  }
}

unique_alphas <- if (length(outlier_labels) > 0) sort(unique(alpha_rounded)) else c()
alpha_labels  <- as.character(unique_alphas)
legend_elements <- c("Cluster 1", "Cluster 2", alpha_labels)
cluster_colors <- c("Cluster 1" = color_cluster1, "Cluster 2" = color_cluster2)

num_alphas <- length(alpha_labels)
if (num_alphas > 0) {
  outlier_colors <- colorRampPalette(c("#FFA500", "#D32F2F", "#7A0010"))(num_alphas)
  names(outlier_colors) <- alpha_labels
} else {
  outlier_colors <- character(0)
}

legend_colors <- c(cluster_colors, outlier_colors)
legend_shapes <- c("Cluster 1" = 16, "Cluster 2" = 16, setNames(rep(4, length(alpha_labels)), alpha_labels))
legend_sizes  <- c("Cluster 1" = 2,  "Cluster 2" = 2,  setNames(rep(4, length(alpha_labels)), alpha_labels))
legend_alphas <- c("Cluster 1" = 0.6, "Cluster 2" = 0.6, setNames(rep(1.0, length(alpha_labels)), alpha_labels))

df_plot$Display_Group <- "Inlier"
if (length(outlier_labels) > 0) {
  df_plot$Display_Group[df_plot$Is_Outlier] <- as.character(df_plot$Alpha_Star[df_plot$Is_Outlier])
}
df_plot$Display_Group[!df_plot$Is_Outlier] <- paste("Cluster", df_plot$Cluster[!df_plot$Is_Outlier])
df_plot$Display_Group <- factor(df_plot$Display_Group, levels = legend_elements)

p <- ggplot(df_plot, aes(x = eruptions, y = waiting)) +
  geom_point(aes(shape = Display_Group, color = Display_Group, size = Display_Group, alpha = Display_Group), stroke = 1.5) +
  geom_text(data = subset(df_plot, Is_Outlier), aes(label = Unit_ID), vjust = 2, color = "black", size = 3, fontface = "bold") +
  coord_fixed(ratio = 1) +
  scale_color_manual(name = expression("Clusters & " * alpha^"*"), values = legend_colors,  breaks = legend_elements) +
  scale_shape_manual(name = expression("Clusters & " * alpha^"*"), values = legend_shapes,  breaks = legend_elements) +
  scale_size_manual(name = expression("Clusters & " * alpha^"*"), values = legend_sizes,   breaks = legend_elements) +
  scale_alpha_manual(name = expression("Clusters & " * alpha^"*"), values = legend_alphas,  breaks = legend_elements) +
  labs(
    subtitle = expression("ODK-means ("* alpha * "=0.01)"),
    x = "eruptions",
    y = "waiting"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = "grey90"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    legend.key.size = unit(0.9, "lines")
  )

print(p)
ggsave("plot2_odkm_alphastar_refined.jpg", plot = p, width = 1920, height = 1080, units = "px", dpi = 300, scale = 1.75)
