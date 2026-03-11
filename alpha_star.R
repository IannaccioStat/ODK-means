##########################################################################################
# FUNCTION: evaluate_alpha
# DESCRIPTION: 
# Optimized calculation of refined significance level alpha*. 
# If a unit index is provided, it only processes that unit's cluster.
# If a rank is provided, it processes all units intra-cluster to find the target.
################################################################################

evaluate_alpha <- function(data, number, U, pos = TRUE, h = NULL, l = NULL) {
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
    
    # Local k-NN setup
    k_log <- floor(log(n_c))
    if (is.null(l)) l_val <- k_log + 1 else l_val <- l
    if (is.null(h)) h_val <- ceiling(k_log / 10) else h_val <- h
    l_val <- min(l_val, n_c - 1)
    
    knn_dist <- get.knn(cluster_data, k = l_val)$nn.dist
    y_local <- if (h_val == l_val) knn_dist[, h_val] else rowSums(knn_dist[, h_val:l_val])
    
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
      lk <- if (is.null(l)) k_log + 1 else l
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
    cluster_idx <- which.max(U[target_unit_idx, ])
    cluster_members <- which(U[, cluster_idx] == 1)
    
    # Isolate the cluster's scores for the backward search
    sorted_y <- sort(global_y[cluster_members])
    m_rank <- which(sorted_y == y_target)[1]
  }
  
  # --- SHARED BACKWARD SEARCH (The Refinement) ---
  if (m_rank < 2) return(1)
  
  y_m <- sorted_y[m_rank]
  mu_m <- mean(sorted_y[1:m_rank])
  sd_m <- sd(sorted_y[1:m_rank])
  alpha_star <- 1 / (((y_m - mu_m) / sd_m)^2 + 1)
  
  # Backward Loop to handle coalescence
  search_idx <- m_rank - 1
  is_coalesced <- FALSE
  
  while (search_idx >= 2) {
    S_j <- sorted_y[1:search_idx]
    y_j <- sorted_y[search_idx]
    mj <- mean(S_j); sj <- sd(S_j)
    if (sj == 0) break
    
    alpha_j <- 1 / (((y_j - mj) / sj)^2 + 1)
    
    if (alpha_j >= alpha_star) {
      alpha_star <- alpha_j
      is_coalesced <- TRUE
    } else {
      if (is_coalesced) break 
    }
    search_idx <- search_idx - 1
  }
  
  return(alpha_star)
}
