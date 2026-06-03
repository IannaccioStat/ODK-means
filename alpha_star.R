# ==============================================================================
# FUNCTION: evaluate_alpha
# DESCRIPTION: 
# Rigorous and streamlined calculation of the refined significance level alpha*.
# It evaluates the naive alpha values from the target unit's rank (m) up to the 
# most isolated unit (n_c) and returns the maximum value as the true alpha*.
# ==============================================================================
evaluate_alpha <- function(data, number, U, pos = TRUE, h = NULL, l = NULL, epsilon = 0.001) {
  library(FNN)
  n_total <- nrow(data)
  
  # --- CASE 1: UNIT INDEX PROVIDED (Optimized Path) ---
  if (!pos) {
    target_unit <- number
    cluster_idx <- which.max(U[target_unit, ])
    cluster_members <- which(U[, cluster_idx] == 1)
    
    # Subset data to the relevant cluster immediately
    cluster_data <- data[cluster_members, , drop = FALSE]
    n_c <- nrow(cluster_data)
    
    # Local k-NN parameter setup
    k_log <- floor(log(n_c))
    lk <- if (is.null(l)) k_log else l
    hk <- if (is.null(h)) ceiling(k_log / 10) else h
    lk <- min(lk, n_c - 1)
    
    knn_dist <- get.knn(cluster_data, k = lk)$nn.dist
    y_local <- if (hk == lk) knn_dist[, hk] else rowSums(knn_dist[, hk:lk])
    
    # Find the target unit's rank within this local cluster
    internal_idx <- which(cluster_members == target_unit)
    y_target <- y_local[internal_idx]
    sorted_y <- sort(y_local)
    m_rank <- which(sorted_y == y_target)[1]
    
  } else {
    # --- CASE 2: GLOBAL RANK PROVIDED (Full Map Path) ---
    global_y <- numeric(n_total)
    K_clusters <- ncol(U)
    
    for (k in 1:K_clusters) {
      members <- which(U[, k] == 1)
      n_c <- length(members)
      if (n_c < 2) { global_y[members] <- 0; next }
      
      # Local parameters for this specific cluster
      k_log <- floor(log(n_c))
      lk <- if (is.null(l)) k_log else l
      hk <- if (is.null(h)) ceiling(k_log / 10) else h
      lk <- min(lk, n_c - 1)
      
      knn_dist <- get.knn(data[members, , drop=FALSE], k = lk)$nn.dist
      global_y[members] <- if (hk == lk) knn_dist[, hk] else rowSums(knn_dist[, hk:lk])
    }
    
    # Find the rank globally
    sorted_global <- sort(global_y)
    y_target <- sorted_global[number]
    
    # Identify which cluster this rank belongs to
    target_unit_idx <- which(global_y == y_target)[1]
    target_unit <- target_unit_idx 
    cluster_idx <- which.max(U[target_unit_idx, ])
    cluster_members <- which(U[, cluster_idx] == 1)
    
    # Isolate the cluster's scores for the search framework
    sorted_y <- sort(global_y[cluster_members])
    m_rank <- which(sorted_y == y_target)[1]
  }
  
  n_c <- length(sorted_y)
  if (m_rank < 2) return(1)
  
  # Scan only from the target unit's rank (m) to the most isolated edge (n_c)
  alphas_to_check <- numeric(n_c - m_rank + 1)
  idx_counter <- 1
  
  for (i in m_rank:n_c) {
    S_i <- sorted_y[1:i]
    sd_i <- sd(S_i)
    
    if (sd_i > 0) {
      alphas_to_check[idx_counter] <- 1 / (((sorted_y[i] - mean(S_i)) / sd_i)^2 + 1)
    } else {
      alphas_to_check[idx_counter] <- 0
    }
    idx_counter <- idx_counter + 1
  }
  
  # The true alpha_star is simply the maximum value in this outer sequence
  alpha_star <- max(alphas_to_check)
  
  return(alpha_star)
}
