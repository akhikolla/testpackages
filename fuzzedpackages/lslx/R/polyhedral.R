## functions for package lslx to implement post-selection inference
## written by Po-Hsien Huang psyphh@gmail.com
## some functions are written based on codes in R package selectiveInference

## \code{compute_tnorm_quantity()} computes quantity for calculating truncated normal probability 
compute_tnorm_quantity <- function(i, a_ph, b_ph, 
                                   debiased_coefficient, coefficient_acov,
                                   is_pen, is_active, is_selected) {
  theta_i <- debiased_coefficient[[i]]
  sigma_i <- sqrt(coefficient_acov[i, i])
  if (is_pen[i]) {
    if (is_selected[i]) {
      coefficient_acov[is.na(coefficient_acov)] <- 0
      a_ph <- a_ph[is_selected, is_active, drop = FALSE]
      b_ph <- b_ph[is_selected, 1, drop = FALSE]
      n_is_active <- sum(is_active)
      cumsum_is_active <- cumsum(is_active)
      c_i <- (coefficient_acov[is_active, i, drop = FALSE]) / (sigma_i^2)
      c_i_expand <- diag(0, n_is_active)
      c_i_expand[, cumsum_is_active[i]] <- c_i
      z_i <-
        (diag(n_is_active) - c_i_expand) %*% matrix(debiased_coefficient[is_active])
      w_i <- c(a_ph %*% c_i)
      r_i <- c(b_ph - (a_ph %*% z_i)) / w_i
      left_i <- ifelse(any(w_i < 0), yes = max(r_i[w_i < 0]), no = - Inf) 
      right_i <- ifelse(any(w_i > 0), yes = min(r_i[w_i > 0]), no = Inf) 
    } else {
      left_i <- NA
      right_i <- NA
    }
  } else {
    left_i <- - Inf
    right_i <- Inf
  }
  tnorm_quantity <- 
    list(theta = theta_i, sigma = sigma_i, left = left_i, right = right_i)
  return(tnorm_quantity)
}

## \code{compute_tnorm_p_value()} computes p-value based on truncated normal
compute_tnorm_p_value <- function(theta, mu, sigma, left, right) {
  if (is.na(left) | is.na(right) | is.na(sigma)) {
    p_value <- NA_real_
  } else {
    if ((left == -Inf) & (right == Inf)) {
      p_value <- 2 * stats::pnorm(- abs((theta - mu) / sigma))
    } else {
      if (theta != 0) {
        if (theta > 0) {
          p_value <- 2 * compute_tnorm_prob(-theta, -mu, sigma, -right, -left)
        } else {
          p_value <- 2 * compute_tnorm_prob(theta, mu, sigma, left, right)
        }
      } else {
        p_value <- NA_real_
      }
    }
  }
  return(p_value)
}

## \code{compute_tnorm_interval()} computes confidence interval based on truncated normal
compute_tnorm_interval <- function(theta, sigma, left, right, alpha_level, 
                                   grid_range, grid_length, depth_max = 3) {
  if (is.na(left) | is.na(right) | is.na(sigma)) {
    lower <- NA_real_
    upper <- NA_real_
  } else {
    if ((left == -Inf) & (right == Inf)) {
      lower <- theta + stats::qnorm(alpha_level / 2) * sigma
      upper <- theta + stats::qnorm(1 - alpha_level / 2) * sigma
    } else {
      grid_point <- 
        seq(from = grid_range[1] * sigma, 
            to = grid_range[2] * sigma, 
            length.out = grid_length)
      grid_prob <- 
        sapply(X = grid_point, FUN = function(mu_i) {
          grid_prob_i <- compute_tnorm_prob(theta, mu_i, sigma, left, right) 
          return(grid_prob_i)
        })
      grid_right <- which(grid_prob <= 1 - alpha_level / 2)
      grid_left <- which(grid_prob >= alpha_level / 2)
      if (length(grid_left) == 0) {
        lower <- - Inf
        upper <- grid_point[1]
      } else if (length(grid_right) == 0) {
        lower <- grid_point[grid_length]
        upper <- Inf
      } else {
        if (grid_right[1] == 1) {
          lower <- -Inf
        } else {
          i_depth <- 1
          lower_left <- grid_point[grid_right[1] - 1]
          lower_right <- grid_point[grid_right[1]]
          while (i_depth < depth_max) {
            grid_point_lower <- 
              seq(from = lower_left, 
                  to = lower_right, 
                  length.out = grid_length)
            grid_prob_lower <- 
              sapply(X = grid_point_lower, FUN = function(mu_i) {
                grid_prob_lower_i <- 
                  compute_tnorm_prob(theta, mu_i, sigma, left, right) 
                return(grid_prob_lower_i)
              })
            grid_lower_right <- which(grid_prob_lower <= 1 - alpha_level / 2)
            lower_left <- grid_point_lower[grid_lower_right[1] - 1]
            lower_right <- grid_point_lower[grid_lower_right[1]]
            i_depth <- i_depth + 1
          }
          lower <- lower_left
        }
        if (grid_left[length(grid_left)] == grid_length) {
          upper <- Inf
        } else {
          i_depth <- 1
          upper_left <- grid_point[grid_left[length(grid_left)]]
          upper_right <- grid_point[grid_left[length(grid_left)] + 1]
          while (i_depth < depth_max) {
            grid_point_upper <- 
              seq(from = upper_left, 
                  to = upper_right, 
                  length.out = grid_length)
            grid_prob_upper <- 
              sapply(X = grid_point_upper, FUN = function(mu_i) {
                grid_prob_upper_i <- 
                  compute_tnorm_prob(theta, mu_i, sigma, left, right) 
                return(grid_prob_upper_i)
              })
            grid_upper_left <- which(grid_prob_upper >= alpha_level / 2)
            upper_left <- grid_point_upper[grid_upper_left[length(grid_upper_left)]]
            upper_right <- grid_point_upper[grid_upper_left[length(grid_upper_left)] + 1]
            i_depth <- i_depth + 1
          }
          upper <- upper_right
        }
      }
    }
  }
  return(c(lower = lower, upper = upper))
}

## \code{compute_tnorm_prob()} computes truncated normal probability
compute_tnorm_prob <- function(theta, mu, sigma, left, right) {
  z_center <- (theta - mu) / sigma
  z_left <- (left - mu) / sigma
  z_right <- (right - mu) / sigma
  if (z_center <= z_left) {
    tnorm_prob <- 0
  } else if (z_center >= z_right) {
    tnorm_prob <- 1
  } else {
    tnorm_prob <- (stats::pnorm(z_center) - stats::pnorm(z_left)) / (stats::pnorm(z_right) - stats::pnorm(z_left)) 
  }
  if (is.na(tnorm_prob)) {
    if (z_left > -Inf) {
      term_left <- compute_bryc(z_left) * exp(-(z_left^2 - z_center^2) / 2)
    } else {
      term_left <- exp(z_center^2)
    }
    if (z_right < Inf) {
      term_right <- compute_bryc(z_right) * exp(-(z_right^2 - z_center^2) / 2)
    } else {
      term_right <- 0
    }
    tnorm_prob <- (term_left - compute_bryc(z_center)) / (term_left - term_right)
  }
  return(tnorm_prob)
}

## \code{compute_bryc()} computes tailed probability
## Bryc, Wlodzimierz. "A uniform approximation to the right normal tail integral." Applied mathematics and computation 127.2 (2002): 365-374.
compute_bryc <- function(x) {
  y <- (x^2 + 5.575192695 * x + 12.7743632) / 
    (sqrt(2 * pi) * x^3 + 14.38718147 * x^2 + 31.53531977 * x + 2 * 12.77436324)
  return(y)
}