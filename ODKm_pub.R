################################################################################
# MODULE: ODK-means Algorithm & Synthetic Data Generation (Tethraset)
# PROJECT: ODK-means: A Simultaneous Approach to Clustering and Outliers Detection
# AUTHOR: Tiziano Iannaccio
# DESCRIPTION: 
# This library implements the ODK-means clustering algorithm, which uses a 
# probabilistic "isolation" metric to detect outliers during the optimization
# process. It also includes the 'Tethraset' generator for 3D benchmarking.
################################################################################

# ------------------------------------------------------------------------------
# FUNCTION: odkm
# DESCRIPTION: The main iterative algorithm. It alternates between cluster 
# assignment, outlier detection via isolation thresholds, and centroid updates
# using weighted observations (where outliers are penalized).
# ------------------------------------------------------------------------------
odkm <- function(data, K, margin = 0, restart = 10, mult = 1, alpha = 0.05) {
  data <- scale(data) 
  nr <- nrow(data) 
  epsilon <- .0000001 
  c_final <- NA; f_final <- Inf; lab_final <- NA
  di <- as.matrix(dist(data))
  
  for (s in 1:restart) {
    idx <- sample(1:nr, K)
    c <- data[idx, ]
    convergence <- FALSE
    f0 <- Inf
    W <- numeric(nr)
    restart_flag = FALSE
    
    while (!convergence) {
      lab <- new_assign(data, c) 
      if (length(unique(lab)) != K) {
        restart_flag = TRUE
        break
      }
      for (k in 1:K) {
        u_k <- which(lab == k)
        W[u_k] <- outdet(u_k, di, margin, mult, alpha)$w
        if (length(u_k) == 1) {
          c[k, ] <- data[u_k, ]
        } else {
          c[k, ] <- colSums(W[u_k] * data[u_k, ])
        }
      }
      f <- 0
      for (i in 1:nr) { f = f + sum((data[i, ] - c[lab[i], ])^2) }
      convergence <- ifelse(f0 - f < epsilon, TRUE, FALSE)
      f0 <- f
    }
    
    if (restart_flag) { s = s - 1; next }
    
    if (f0 < f_final) {
      f_final <- f0; c_final <- c; lab_final <- lab
    }
  }
  
  U <- matrix(0, nr, K)
  for (i in 1:nr) { U[i, lab_final[i]] = 1 }
  outliers <- numeric(0)
  for (k in 1:K) {
    lab <- which(U[, k] == 1)
    outliers <- c(outliers, outdet(lab, di, margin, mult, alpha)$out)
  }
  if (identical(numeric(0), outliers)) { outliers = NULL }
  return(list("U" = U, "c" = c_final, "f" = f_final, "out" = sort(outliers)))
}

# ------------------------------------------------------------------------------
# FUNCTION: new_assign
# DESCRIPTION: Performs the E-step (Assignment) of the algorithm. It assigns 
# each data point to the cluster of the nearest centroid using Manhattan distance.
# ------------------------------------------------------------------------------
new_assign <- function(data, c) {
  k <- nrow(c); n <- nrow(data); j <- ncol(data)
  dif <- matrix(NA, k, n)
  for (idx in 1:k) {
    dif[idx, ] <- rowSums(abs(matrix(c[idx, ], n, j, TRUE) - data))
  }
  clusters <- numeric(n)
  for (i in 1:n) { clusters[i] <- which.min(dif[, i]) }
  return(clusters)
}

# ------------------------------------------------------------------------------
# FUNCTION: outdet
# DESCRIPTION: The core outlier detection engine. It calculates pseudo-isolation 
# sums for a subset of data and applies a probabilistic threshold (Cantelli) 
# based on the significance level alpha to identify anomalous units.
# ------------------------------------------------------------------------------
outdet <- function(labels, distance, margin = 0, mult = 1, alpha = 0.05) {
  dist_sub <- distance[labels, labels]
  l <- length(labels)
  if (l == 1) return(list("k" = 1, "sums" = numeric(0), "out" = numeric(0), "w" = 1))
  if (l == 2) return(list("k" = 1, "sums" = numeric(0), "out" = numeric(0), "w" = c(0.5, 0.5)))
  
  for (i in 1:l) { dist_sub[i, ] <- sort(dist_sub[i, ]) }
  k <- ifelse(l > 2, floor(log(l)), 1)
  h <- ceiling(k / 10)
  
  if (l < (2 * k + 1)) {
    sums_k <- sort(rowSums(dist_sub))
  } else {
    sums_k <- sort(rowSums(dist_sub[, h:(k + 1)]))
  }
  
  tvec <- numeric(0); told <- Inf; deltat <- Inf; epsilon <- .00001
  K_val <- sqrt((1 / alpha) - 1)
  
  while (deltat > 0 + epsilon) {
    negatives <- sums_k[which(sums_k <= told)]
    m <- mean(negatives); s <- sd(negatives)
    tnew <- m + K_val * s
    tvec <- c(tvec, tnew); deltat <- told - tnew; told <- tnew
  }
  
  outl <- which(sums_k > told)
  w <- penalty(sums_k, outl, mult)
  return(list("k" = k, "sums" = sums_k, "out" = as.numeric(names(sums_k[outl])), "w" = w, "tvec" = tvec))
}

# ------------------------------------------------------------------------------
# FUNCTION: penalty
# DESCRIPTION: Calculates the weights for data points. Outliers identified by 
# outdet receive a weight proportional to their isolation, effectively reducing 
# their influence on the movement of centroids.
# ------------------------------------------------------------------------------
penalty <- function(sums, outliers, mult) {
  pen <- rep(1, length(sums))
  non_outliers <- if (length(outliers) > 0) sums[-outliers] else sums
  m <- mean(non_outliers); se <- sd(non_outliers)
  if (length(outliers) > 0) {
    pen[outliers] <- mult * (se / (sums[outliers] - m))^2
  } 
  return(pen / sum(pen))
}

# ------------------------------------------------------------------------------
# FUNCTION: metrics
# DESCRIPTION: Computes the F1-Score and Confusion Matrix for outlier detection 
# performance by comparing estimated outliers against the ground truth.
# ------------------------------------------------------------------------------
metrics <- function(n, true, est) {
  true_m <- est_m <- matrix(c(0, 1), n, 2, byrow = T)
  est_m[est, ] <- matrix(c(1, 0), length(est), 2, byrow = T)
  true_m[true, ] <- matrix(c(1, 0), length(true), 2, byrow = T)
  conf_mat <- t(est_m) %*% true_m
  precision <- conf_mat[1, 1] / (conf_mat[1, 1] + conf_mat[1, 2])
  sensitivity <- conf_mat[1, 1] / (conf_mat[1, 1] + conf_mat[2, 1])
  F1 <- (2 * precision * sensitivity) / (precision + sensitivity)
  return(list("measure" = F1, "conf" = conf_mat))
}

# ------------------------------------------------------------------------------
# FUNCTION: tethraset
# DESCRIPTION: Generates the "Tethra" synthetic dataset: four 3D clusters 
# arranged in a tetrahedron with added Internal and External outliers.
# ------------------------------------------------------------------------------
tethraset <- function(nk, perc_out, var_fluc, side, r) {
  if (perc_out > 0.05) perc_out <- 0.05
  Z <- matrix(rnorm(nk * 3), nk, 3)
  X <- t(apply(Z, 1, function(x) { x / (sqrt(sum(x * x))) })) * r
  
  centers <- matrix(c(0,0,0, side,0,0, side/2,side,0, 
                      side/2,side/2,sqrt(side^2-2*(side/2)^2)), 4, 3, T)
  
  data <- rbind(X + matrix(rnorm(nk*3, 0, var_fluc), nk, 3),
                X + matrix(rnorm(nk*3, 0, var_fluc), nk, 3) + matrix(centers[2,], nk, 3, T),
                X + matrix(rnorm(nk*3, 0, var_fluc), nk, 3) + matrix(centers[3,], nk, 3, T),
                X + matrix(rnorm(nk*3, 0, var_fluc), nk, 3) + matrix(centers[4,], nk, 3, T))
  
  cubes <- array(NA, c(3, 2, 4))
  for (c in 1:4) { cubes[, , c] <- smallcube(var_fluc + 1.5 * r, centers[c, ]) }
  
  if (perc_out > 0) {
    toberemoved <- floor(perc_out * (4 * nk))
    data <- data[-sample(1:nrow(data), toberemoved), ]
    m <- floor(toberemoved * 0.3)
    inl <- centers[sample(1:4, m, T), ] + matrix(rnorm(3 * m, 0, var_fluc), m, 3)
    outl <- bigcube(matrix(0, 0, 3), toberemoved - m, side, r, cubes)
    data <- rbind(data, inl, outl)
    return(list("data" = data, "outliers" = (nrow(data) - toberemoved + 1):nrow(data), "cubes" = cubes))
  }
  return(list("data" = data, "outliers" = NULL, "cubes" = cubes))
}

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS: smallcube & bigcube
# DESCRIPTION: Geometric utilities for Tethraset. 'smallcube' defines cluster 
# boundaries; 'bigcube' generates External outliers using oversampling.
# ------------------------------------------------------------------------------
smallcube <- function(dev, coord) {
  mat <- matrix(NA, 3, 2)
  for (i in 1:3) { mat[i, 1] <- coord[i] - dev; mat[i, 2] <- coord[i] + dev }
  return(mat)
}

bigcube <- function(ag, tbg, s, r, cube) {
  over <- 3; total <- tbg * over
  ng <- matrix(runif(3 * total, -3 * r, s + 3 * r), total, 3)
  in_any <- rep(FALSE, total)
  for (c in 1:4) {
    in_any <- in_any | (ng[,1] > cube[1,1,c] & ng[,1] < cube[1,2,c] &
                          ng[,2] > cube[2,1,c] & ng[,2] < cube[2,2,c] &
                          ng[,3] > cube[3,1,c] & ng[,3] < cube[3,2,c])
  }
  valid <- ng[!in_any, ]
  if (nrow(valid) >= tbg) return(valid[1:tbg, ])
  return(rbind(valid, bigcube(matrix(0, 0, 3), tbg - nrow(valid), s, r, cube)))
}