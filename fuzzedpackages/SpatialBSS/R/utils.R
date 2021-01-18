white_data <- function(x) {
  n <- nrow(x)
  
  mu <- colMeans(x)
  x_0 <- sweep(x, MARGIN = 2, STATS = mu, FUN = '-')
  
  s <- crossprod(x_0) / (n - 1)
  
  s_evd <- eigen(s, symmetric = TRUE)
  s_inv_sqrt <- s_evd$vectors %*% tcrossprod(diag(1 / sqrt(s_evd$values)), s_evd$vectors)
  s_sqrt <- s_evd$vectors %*% tcrossprod(diag(sqrt(s_evd$values)), s_evd$vectors)
  
  x_w <- tcrossprod(x_0, s_inv_sqrt)
  colnames(x_w) <- colnames(x_0)
  
  return(list(mu = mu, x_0 = x_0, x_w = x_w, s_inv_sqrt = s_inv_sqrt, s_sqrt = s_sqrt))
}

spatial_kernel_matrix <- function(coords, kernel_type = c('ring', 'ball', 'gauss'), kernel_parameters) {
  kernel_type <- match.arg(kernel_type)
  
  if (kernel_type == 'ball') {
    kernel_list <- lapply(kernel_parameters, 
                          function(h) if (h >= 0) k_mat_ball(coords = coords, h = h)
                                      else stop('Radius must be zero or positive.'))    
  } else if (kernel_type == 'ring' && (length(kernel_parameters) %% 2 == 0)) {
    kernel_list <- lapply(seq(from = 1, to = length(kernel_parameters), by = 2), 
                          function(idx) if (kernel_parameters[idx] >= kernel_parameters[idx + 1]) stop('Inner radius must be smaller than outer radius.') 
                                        else k_mat_ring(coords = coords, h1 = kernel_parameters[idx], h2 = kernel_parameters[idx + 1]))  
  } else if (kernel_type == 'gauss') {
    kernel_list <- lapply(kernel_parameters, 
                          function(h) if (h >= 0) k_mat_exp(coords = coords, h = h)
                                      else stop('Parameter must be zero or positive.'))
  } else {
    stop('Invalid input. Note that the length(kernel_parameters) must be an even number for the ring kernel.')
  }
  
  return(kernel_list)
}

local_covariance_matrix <- function(x, kernel_list, lcov = c('lcov', 'ldiff'), whitening = TRUE) {
  lcov <- match.arg(lcov)
  
  if (whitening) {
    x <- white_data(x)$x_w
  }
  
  cov_sp_list <- switch(lcov,
    'lcov'  = lapply(kernel_list, function(k_mat) sp_lcov_sparse(x = x, k = k_mat)),
    'ldiff' = lapply(kernel_list, function(k_mat) sp_ldiff_sparse(x = x, k = k_mat))
  )
  
  cov_sp_list <- lapply(cov_sp_list, function(x) (x + t(x)) / 2)
  
  attr(cov_sp_list, 'lcov') <- lcov
  
  return(cov_sp_list)
}

predict_idw <- function(vals, coords, p, n_grid) {
  coords_pred <- as.matrix(expand.grid(seq(from = floor(min(coords[, 1])), to = ceiling(max(coords[, 1])), length.out = n_grid), 
                                       seq(from = floor(min(coords[, 2])), to = ceiling(max(coords[, 2])), length.out = n_grid)))
  colnames(coords_pred) <- colnames(coords)
  vals_pred <- idw(coords_pred = coords_pred, coords_vals = coords, vals = vals, p = p)
  colnames(vals_pred) <- paste0(colnames(vals), '.pred')
  return(list(vals_pred_idw = vals_pred, coords_pred_idw = coords_pred))
}

diag_scatters <- function(cov_list, ordered, ...) {
  decr <- if (attr(cov_list, 'lcov') == 'ldiff') FALSE else TRUE
  k <- length(cov_list)
  if (k == 1) {
    cov_evd <- eigen(cov_list[[1]], symmetric = TRUE)
    u <- cov_evd$vectors
    d <- diag(cov_evd$values)
  } else {
    jade <- JADE::frjd(do.call(rbind, cov_list), ...)
    u <- jade$V
    d <- jade$D
  }
  
  p <- ncol(d)
  if (ordered) {
    diags_mat <- matrix(0, nrow = k, ncol = p)
    for (idx in 1:k) {
      diags_mat[idx, ] <- diag(d[(1:p) + (idx - 1) * p, ])
    }
    diag_order <- order(colSums(diags_mat ^ 2), decreasing = decr)
    u <- u[, diag_order]
    for (idx in 1:k) {
      d[(1:p) + (idx - 1) * p, ] <- d[(1:p) + (idx - 1) * p, ][diag_order, diag_order]
    }
  } 
  
  return(list(u = u, d = d))
}
