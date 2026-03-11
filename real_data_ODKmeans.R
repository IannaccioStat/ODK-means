################################################################################
# SCRIPT: Real Data Application - Old Faithful Geyser Dataset
# PROJECT: ODK-means: A Simultaneous Approach to Clustering and Outliers Detection
# AUTHOR: Tiziano Iannaccio
# DESCRIPTION: 
# This script executes the comparative analysis presented in Section 8 of the 
# manuscript. It evaluates ODK-means against Trimmed K-means and LOF using 
# the 'faithful' dataset. The analysis highlights:
#  1. Sensitivity to alpha vs. the search for the 'Top N' isolated points.
#  2. The relationship between the probabilistic threshold alpha and the 
#     empirical trimming proportion p.
#  3. The topological difference between distance-based and isolation-based 
#     outlier identification.
################################################################################

# NOTICE THAT COLORS MAY CHANGE DUE TO LABEL SWITCH

# ==============================================================================
# IMPORT ALL THE NEEDED PACKAGES
# ==============================================================================
library(ggplot2)
library(gridExtra)
library(trimcluster)
library(Rlof)
library(tclust)
library(FNN)

# ==============================================================================
# IMPORT ODKmeans
# ==============================================================================
source("ODKm.R")

# ==============================================================================
# INITIAL SETUP
# ==============================================================================
data(faithful)
X <- scale(as.matrix(faithful))
set.seed(123)

color_cluster1 <- "#440154" 
color_cluster2 <- "#21908C" 
color_outlier  <- "red"
color_centroid <- "orange"

panel_theme <- theme_minimal() + 
  theme(plot.margin = margin(5, 5, 5, 5), legend.position = "right")

# This ensures "Cluster" is always the top legend
legend_order <- guides(
  color = guide_legend(order = 1),
  shape = guide_legend(order = 2)
)

# ==============================================================================
# 1. ODKM (STANDARD) vs ODKM (TOP 5 ISOLATED)
# ==============================================================================
res_alpha05 <- odkm(data = X, K = 2, restart = 20, alpha = 0.05)
cluster_assignments <- apply(res_alpha05$U, 1, which.max)
df_odkm <- as.data.frame(X)
df_odkm$Cluster <- as.factor(cluster_assignments)
df_odkm$ID <- 1:nrow(X)

p1a <- ggplot(df_odkm, aes(x = eruptions, y = waiting)) +
  geom_point(aes(color = Cluster), size = 2, alpha = 0.6) +
  geom_point(data = df_odkm[1:nrow(X) %in% res_alpha05$out, ], 
             aes(shape = "Outlier"), color = color_outlier, size = 4, stroke = 1.5) +
  geom_text(data = df_odkm[1:nrow(X) %in% res_alpha05$out, ], 
            aes(label = ID), vjust = 2, color = color_outlier, size = 3, fontface = "bold") +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("1" = color_cluster1, "2" = color_cluster2)) +
  scale_shape_manual(values = c("Outlier" = 4)) +
  labs(title = expression(paste("ODK-means (", alpha, " = 0.05)")), 
       color = "Cluster", shape = "Reference") + panel_theme + legend_order

od_info_full <- outdet(1:nrow(X), as.matrix(dist(X)), alpha = 0.05)
top5_iso_idx <- as.numeric(names(tail(od_info_full$sums, 5)))

p1b <- ggplot(df_odkm, aes(x = eruptions, y = waiting)) +
  geom_point(aes(color = Cluster), size = 2, alpha = 0.6) +
  geom_point(data = df_odkm[1:nrow(X) %in% top5_iso_idx, ], 
             aes(shape = "Outlier"), color = color_outlier, size = 4, stroke = 1.5) +
  geom_text(data = df_odkm[1:nrow(X) %in% top5_iso_idx, ], 
            aes(label = ID), vjust = 2, color = color_outlier, size = 3, fontface = "bold") +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("1" = color_cluster1, "2" = color_cluster2)) +
  scale_shape_manual(values = c("Outlier" = 4)) +
  labs(title = "ODK-means (Top 5 Isolated)", 
       color = "Cluster", shape = "Reference") + panel_theme + legend_order

comp1 <- grid.arrange(p1a, p1b, ncol = 2)
ggsave("comp1_odkm_std_vs_iso.jpg", comp1, width = 12, height = 6, dpi = 300)

# ==============================================================================
# 2. ODKM (ALPHA STAR)
# ==============================================================================
alpha_triplet <- 0.0776 
res_triplet <- odkm(data = X, K = 2, restart = 20, alpha = alpha_triplet)
df_triplet <- df_odkm
df_triplet$Cluster <- as.factor(apply(res_triplet$U, 1, which.max))

p2 <- ggplot(df_triplet, aes(x = eruptions, y = waiting)) +
  geom_point(aes(color = Cluster), size = 2, alpha = 0.6) +
  geom_point(data = df_triplet[1:nrow(X) %in% res_triplet$out, ], 
             aes(shape = "Outlier"), color = color_outlier, size = 4, stroke = 1.5) +
  geom_text(data = df_triplet[1:nrow(X) %in% res_triplet$out, ], 
            aes(label = ID), vjust = 2, color = color_outlier, size = 3, fontface = "bold") +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("1" = color_cluster1, "2" = color_cluster2)) +
  scale_shape_manual(values = c("Outlier" = 4)) +
  labs(title = bquote("ODK-means (" * alpha == .(round(alpha_triplet, 4)) * ")"), 
       color = "Cluster", shape = "Reference") + theme_minimal() + legend_order

print(p2)
ggsave("plot2_odkm_alphastar.jpg", p2, width = 1920, height = 1080, units = "px", dpi = 300, scale = 1.75)

# ==============================================================================
# 3. TKM (EVALUATION OF P)
# ==============================================================================
jpeg("plot3_tclust_diagnostic.jpg", width = 1920, height = 1080, res = 300)
ctrl <- ctlcurves(x = X, k = 1:3, alpha = seq(0, 0.5, by = 0.01))
plot(ctrl, xlab = "p (Trimming Proportion)", main = "TClust Diagnostic Curves")
dev.off()
plot(ctrl, xlab = "p (Trimming Proportion)", main = "TClust Diagnostic Curves")

# ==============================================================================
# 4. TKM (P=0.018) vs ODKM (TOP 5 DISTANT TO CENTROID)
# ==============================================================================
p_val <- 0.018
t_km <- trimkmeans(X, k = 2, trim = p_val)
df_tkm <- as.data.frame(X)
df_tkm$Cluster <- as.factor(t_km$classification)
df_tkm$ID <- 1:nrow(X)

# --- FIXED CENTROID LOGIC FOR TRIMKMEANS ---
# The package trimcluster uses $means, not $centers
tkm_cents <- as.data.frame(t_km$means) 
colnames(tkm_cents) <- c("eruptions", "waiting")

p4a <- ggplot(df_tkm, aes(x = eruptions, y = waiting)) +
  geom_point(aes(color = Cluster), size = 2, alpha = 0.5) +
  geom_point(data = subset(df_tkm, Cluster == 3), aes(shape = "Outlier"), color = color_outlier, size = 4, stroke = 1.5) +
  geom_point(data = tkm_cents, aes(x = eruptions, y = waiting, shape = "Centroid"), color = color_centroid, size = 5, stroke = 2) +
  geom_text(data = subset(df_tkm, Cluster == 3), aes(label = ID), vjust = 2, color = color_outlier, size = 2.5, fontface = "bold") +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("1" = color_cluster1, "2" = color_cluster2, "3" = "grey80")) +
  scale_shape_manual(values = c("Outlier" = 4, "Centroid" = 8)) +
  labs(title = paste0("Trimmed K-means (p = ", p_val, ")"), color = "Cluster", shape = "Reference") + panel_theme + legend_order

dist_to_c <- sqrt(rowSums((X - res_alpha05$c[cluster_assignments, ])^2))
top5_dist_idx <- order(dist_to_c, decreasing = TRUE)[1:5]
centroids_df <- as.data.frame(res_alpha05$c); colnames(centroids_df) <- c("eruptions", "waiting")

p4b <- ggplot(df_odkm, aes(x = eruptions, y = waiting)) +
  geom_point(aes(color = Cluster), size = 2, alpha = 0.6) +
  geom_point(data = df_odkm[1:nrow(X) %in% top5_dist_idx, ], aes(shape = "Outlier"), color = color_outlier, size = 4, stroke = 1.5) +
  geom_point(data = centroids_df, aes(x = eruptions, y = waiting, shape = "Centroid"), color = color_centroid, size = 5, stroke = 2) +
  geom_text(data = df_odkm[1:nrow(X) %in% top5_dist_idx, ], aes(label = ID), vjust = 2, color = color_outlier, size = 3, fontface = "bold") +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("1" = color_cluster1, "2" = color_cluster2)) +
  scale_shape_manual(values = c("Outlier" = 4, "Centroid" = 8)) +
  labs(title = "ODK-means (Top 5 Distance)", color = "Cluster", shape = "Reference") + panel_theme + legend_order

comp4 <- grid.arrange(p4a, p4b, ncol = 2)
ggsave("comp4_tkm_vs_odkdist.jpg", comp4, width = 12, height = 6, dpi = 300)

# ==============================================================================
# 5. LOF (EVALUATION OF LOF)
# ==============================================================================
scores <- lof(X, k = 5)
sorted_scores <- sort(scores, decreasing = TRUE)
lof_top5_val <- sorted_scores[5]
prop_top5 <- 5 / nrow(X)
df_scree <- data.frame(rank = 1:nrow(X), score = sorted_scores)

p5 <- ggplot(df_scree, aes(x = rank, y = score)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 1.5, alpha = 0.6) +
  # Top 5 Line and Annotation
  geom_hline(yintercept = lof_top5_val, color = "darkgreen", linetype = "dashed") +
  annotate("text", x = 25, y = lof_top5_val, 
           label = paste0("p = ", round(prop_top5, 2), ", LOF = ", round(lof_top5_val, 2)),
           vjust = -1, color = "darkgreen", fontface = "bold") +
  # LOF 1.38 Line and Annotation
  geom_hline(yintercept = 1.38, color = "red", linetype = "dotted") +
  annotate("text", x = 25, y = 1.38, label = "LOF = 1.38",
           vjust = -1, color = "red", fontface = "bold") +
  labs(title = "LOF Score Scree Plot", x = "Observation Rank", y = "LOF Score") +
  coord_cartesian(xlim = c(1, 50)) + 
  theme_minimal()

print(p5)
ggsave(filename = "plot5_lof_scree.jpg", plot = p5, width = 1920, height = 1080, units = "px", dpi = 300, scale = 1.75)

# ==============================================================================
# 6. LOF (TOP 5) vs LOF (1.38)
# ==============================================================================
km_std <- kmeans(X, centers = 2, nstart = 20)
df_lof <- as.data.frame(X); df_lof$Cluster <- as.factor(km_std$cluster)

p6a <- ggplot(df_lof, aes(x = eruptions, y = waiting)) +
  geom_point(aes(color = Cluster), size = 2, alpha = 0.5) +
  geom_point(data = df_lof[scores >= lof_top5_val, ], aes(shape = "Outlier"), color = color_outlier, size = 4, stroke = 1.5) +
  geom_text(data = df_lof[scores >= lof_top5_val, ], aes(label = which(scores >= lof_top5_val)), vjust = 2, color = color_outlier, size = 2.5, fontface = "bold") +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("1" = color_cluster1, "2" = color_cluster2)) +
  scale_shape_manual(values = c("Outlier" = 4)) +
  labs(title = "LOF (Top 5 Isolated)", color = "Cluster", shape = "Reference") + panel_theme + legend_order

p6b <- ggplot(df_lof, aes(x = eruptions, y = waiting)) +
  geom_point(aes(color = Cluster), size = 2, alpha = 0.5) +
  geom_point(data = df_lof[scores >= 1.38, ], aes(shape = "Outlier"), color = color_outlier, size = 4, stroke = 1.5) +
  geom_text(data = df_lof[scores >= 1.38, ], aes(label = which(scores >= 1.38)), vjust = 2, color = color_outlier, size = 2.5, fontface = "bold") +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("1" = color_cluster1, "2" = color_cluster2)) +
  scale_shape_manual(values = c("Outlier" = 4)) +
  labs(title = "LOF (Threshold 1.38)", color = "Cluster", shape = "Reference") + panel_theme + legend_order

comp6 <- grid.arrange(p6a, p6b, ncol = 2)
print(comp6)
ggsave("comp6_lof_comparison.jpg", comp6, width = 12, height = 6, dpi = 300)
