################################################################################
# SCRIPT: Figure 5 - Tethra Dataset Generation 
# PROJECT: ODK-means: A Simultaneous Approach to Clustering and Outliers Detection
# AUTHOR: Tiziano Iannaccio
# DESCRIPTION: 
# This script generates and visualizes the synthetic 3D "Tethra" dataset. 
# It illustrates the geometric framework of the ODK-means model by plotting 
# clusters within local hyper-rectangles and a global bounding domain.
################################################################################

# Reproducibility
set.seed(42) 

# Data Loading
source("ODKm.R")
data <- tethraset(200, 0.05, 0.02, 4.5, 1)
cubes <- data$cubes

library(scatterplot3d)

# --- FUNCTIONS ---

# Function to draw small cube (wireframe)
draw_cube_scatterplot3d <- function(cube_matrix, s3d, color = "red") {
  x_min <- cube_matrix[1, 1]; x_max <- cube_matrix[1, 2]
  y_min <- cube_matrix[2, 1]; y_max <- cube_matrix[2, 2]
  z_min <- cube_matrix[3, 1]; z_max <- cube_matrix[3, 2]
  
  vertices <- matrix(c(
    x_min, y_min, z_min, x_max, y_min, z_min, x_max, y_max, z_min, x_min, y_max, z_min,
    x_min, y_min, z_max, x_max, y_min, z_max, x_max, y_max, z_max, x_min, y_max, z_max 
  ), ncol = 3, byrow = TRUE)
  
  edges <- list(c(1, 2), c(2, 3), c(3, 4), c(4, 1), c(5, 6), c(6, 7), c(7, 8), c(8, 5), c(1, 5), c(2, 6), c(3, 7), c(4, 8))
  
  for (edge in edges) {
    s3d$points3d(vertices[edge, ], type = "l", col = color, lwd = 2)
  }
}

# Function to draw big bounding cube
draw_big_cube_scatterplot3d <- function(s, r, s3d, color = "gray") {
  x_min <- -3 * r; x_max <- s + 3 * r
  y_min <- -3 * r; y_max <- s + 3 * r
  z_min <- -3 * r; z_max <- s + 3 * r
  
  vertices <- matrix(c(
    x_min, y_min, z_min, x_max, y_min, z_min, x_max, y_max, z_min, x_min, y_max, z_min,
    x_min, y_min, z_max, x_max, y_min, z_max, x_max, y_max, z_max, x_min, y_max, z_max 
  ), ncol = 3, byrow = TRUE)
  
  edges <- list(c(1, 2), c(2, 3), c(3, 4), c(4, 1), c(5, 6), c(6, 7), c(7, 8), c(8, 5), c(1, 5), c(2, 6), c(3, 7), c(4, 8))
  
  for (edge in edges) {
    s3d$points3d(vertices[edge, ], type = "l", col = color, lwd = 3, lty = 2)
  }
}

# --- PLOTTING & EXPORT ---

# Saving high-resolution JPG with 12:9 proportion
# 1200x900 pixels maintains the 12:9 ratio at 150-300 DPI quality
jpeg("Figure5_Tethra.jpg", width = 1200, height = 900, quality = 100, res = 150)

radius <- 1
side <- 4.5

x_min <- -3 * radius; x_max <- side + 3 * radius
y_min <- -3 * radius; y_max <- side + 3 * radius
z_min <- -3 * radius; z_max <- side + 3 * radius

dummy_points <- matrix(c(x_min, y_min, z_min, x_max, y_max, z_max), ncol = 3, byrow = TRUE)

s3d <- scatterplot3d(dummy_points[,1], dummy_points[,2], dummy_points[,3],
                     xlim = c(x_min, x_max), ylim = c(y_min, y_max), zlim = c(z_min, z_max),
                     color = "white", pch = NA, 
                     xlab = expression(x[1]), ylab = expression(x[2]), zlab = expression(x[3]),
                     main = "", angle = 45)

# Layer 1: Big Cube
draw_big_cube_scatterplot3d(4.5, 1, s3d, color = "gray")

# Layer 2: Data Points
items <- data$data
col <- rep("blue", nrow(items))
col[data$outliers] <- "red"
s3d$points3d(items[,1], items[,2], items[,3], col = col, pch = 16, cex = 0.7)

# Layer 3: Small Cubes
draw_cube_scatterplot3d(cubes[,,1], s3d, color = "black")
draw_cube_scatterplot3d(cubes[,,2], s3d, color = "black") 
draw_cube_scatterplot3d(cubes[,,3], s3d, color = "black")
draw_cube_scatterplot3d(cubes[,,4], s3d, color = "black")

# Layer 4: Opaque Legend
legend("topright", 
       legend = c("Regular", "Outliers", "Small cubes","Big cube"),
       col = c("blue", "red", "black", "gray"),
       pch = c(16, 16, NA, NA),
       lty = c(NA, NA, 1, 2),
       lwd = c(NA, NA, 2, 1),
       cex = 0.8,
       bg = "white")

dev.off() # Closes the file and saves to working directory
