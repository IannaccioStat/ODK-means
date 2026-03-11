################################################################################
# SCRIPT: Figure 4 - Conceptual Map: Cell-wise vs. Row-wise Outliers
# PROJECT: ODK-means: A Simultaneous Approach to Clustering and Outliers Detection
# AUTHOR: Tiziano Iannaccio
# DESCRIPTION: 
# This script generates a conceptual visualization to distinguish between 
# different types of anomalies. It simulates a concentrated central mass of 
# regular observations and injects a "Cell-wise" outlier (atypical in only one 
# dimension) and a "Row-wise" outlier (atypical in the multivariate space).
################################################################################

# ==============================================================================
# 1. SETUP AND PACKAGE LOAD
# ==============================================================================
library(ggplot2)
library(dplyr)
library(mvtnorm) 

# ==============================================================================
# 2. GENERATE AND CONSOLIDATE THE DATASET
# ==============================================================================
set.seed(42) # Fixed seed for publication-consistent plot

# --- Core Data (Central Mass) ---
n_regular <- 500
# REDUCED SIGMA: lower values (0.2) make the data much more concentrated
# We also reduced the covariance (0.05) to keep the "ball" shape tight
tight_sigma <- matrix(c(0.2, 0.05, 0.05, 0.2), 2)

df_core <- as.data.frame(rmvt(n_regular, delta = c(0, 0), sigma = tight_sigma, df = 5))
colnames(df_core) <- c("Feature_1", "Feature_2")
df_core$Anomalous <- "Regular"

# --- Outlier 1 (Univariate / Cell-wise) ---
outlier_1 <- data.frame(
  Feature_1 = runif(1, min = -0.2, max = 0.2), # Kept tight on X
  Feature_2 = 6,                               # High on Y
  Anomalous = "Cell-wise Outlier"
)

# --- Outlier 2 (Multivariate / Row-wise) ---
outlier_2 <- data.frame(
  Feature_1 = 6,                               # High on X
  Feature_2 = 5.8,                             # High on Y
  Anomalous = "Row-wise Outlier"
)

df_combined <- bind_rows(df_core, outlier_1, outlier_2)

# ==============================================================================
# 3. CREATE THE DATA VISUALIZATION (GGPLOT2)
# ==============================================================================



p_anomalies <- ggplot(df_combined, aes(x = Feature_1, y = Feature_2, color = Anomalous, shape = Anomalous)) +
  # Regular points
  geom_point(data = subset(df_combined, Anomalous == "Regular"), alpha = 0.5, size = 1.8) +
  # Outliers: bold stroke for thickness
  geom_point(data = subset(df_combined, Anomalous != "Regular"), size = 6, stroke = 1.8) +
  
  scale_color_manual(values = c("Regular" = "deepskyblue4", 
                                "Cell-wise Outlier" = "#D55E00", 
                                "Row-wise Outlier" = "#CC79A7")) +
  
  # Note: shape 16 is a solid circle for the regular points
  scale_shape_manual(values = c("Regular" = 16, "Cell-wise Outlier" = 4, "Row-wise Outlier" = 8)) +
  
  # Use identical titles for both scales to merge legends
  labs(
    title = "Conceptual Map: Cell-wise vs. Row-wise Outliers",
    x = "Feature 1",
    y = "Feature 2",
    color = "Observation Type",
    shape = "Observation Type"
  ) +
  theme_minimal() + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "right") +
  
  # Merge legends and control aesthetics in the legend box
  guides(
    color = guide_legend(override.aes = list(alpha = 1, size = 5, stroke = 1.5)),
    shape = guide_legend(override.aes = list(size = 5, stroke = 1.5))
  ) +
  expand_limits(x = 0, y = 0)

# Display the plot in the console/IDE
print(p_anomalies)

# Save the high-resolution JPG with 12:9 (4:3) equivalent proportion
# 1920x1440 would be 4:3; 1920x1080 is 16:9. 
# Keeping your requested 1920x1080 for high-definition widescreen.
ggsave("Figure4_ConceptualMap.jpg", p_anomalies, width = 1920, height = 1080, units = "px", dpi = 300)
