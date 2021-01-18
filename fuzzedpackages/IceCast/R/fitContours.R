# -------------------------------------
# main function
# -------------------------------------
#' Sets up MCMC to fit the parameters of the contour Model in R, then runs
#' the sampler in C++
#' @param r number indicating which region in the \code{reg_info} list is being
#'          considered
#' @param n_iter number of iterations to run the MCMC, must be a multiple of \code{w}
#' @param y_obs output of \code{y_obs} function. This is a list of matrices,
#'              one per region, giving the observed \eqn{y} values. Each row
#'              corresponds to the lines and each column corresponds to
#'              a training year
#' @param reg_info  a \code{reg_info} list (see documentation for \code{reg_info})
#' @param dists symmetric matrix of the same dimension as the number of
#'                lines being used, specifying distances among indices.
#'                Defaults to \code{NULL}, which means atrix will be computed by
#'                the \code{dist_mat} function
#' @param sigma_min minimum value for all \eqn{\sigma} parameters.
#'                 Typically close to but not exactly zero (defaults to 0.01).
#'                 Not used if \code{sigma0_lb} is set to \code{NULL}
#' @param sigma0_lb vector of the same length as the number of lines which
#'                 specifies the lower bound of the uniform prior for each sigma
#'                 value. Defaults to \code{NULL}, meaning \code{sigma0_lb} is
#'                 set to be a vector with all values set to \code{sigmaMin}
#' @param sigma0_ub vector of the same length as the number of lines which
#'                 specifies the upper bound of the uniform prior for each sigma
#'                 value. Defaults to \code{NULL}
#' @param xU_prop_sd_def Standard deviation for proposals for xU when xU can take on
#'                    an infinite set of values
#' @param mu_ini vector of the same length as the number of lines which
#'               specifies the values from which each element of \eqn{\mu} will
#'               be initialized in the MCMC. Defaults to \code{NULL}, meaning
#'               \eqn{\mu} will be initialized with the mean of the observed y's
#' @param mu0 vector of the same length as the number of lines which specifies
#'            the prior mean for \eqn{\mu}. Defaults to \code{NULL}, meaning
#'            each element in \code{mu0} will be set to be in the middle of its
#'            corresponding line
#' @param lambda0 matrix of the same dimension as the number of lines which
#'                specifices the prior covariance matrix for \eqn{\mu}. Defaults
#'                to \code{NULL}, which gives a diagonal matrix with diagonal
#'                elements corresponding to the variance that would be required
#'                for 80% of the data to fall between the minimum and maximum
#'                values of the corresponding line if the data were normally
#'                distributed.
#' @param sigma_ini vector of the same length as the number of lines which
#'                  specifies the values from which each element in \eqn{\Sigma}
#'                  will be initialized from. Defaults to \code{NULL}, meaning
#'                  each element of \eqn{\Sigma} will be initialized with
#'                  the observed standard deviation of its corresponding y's,
#'                  bounded by \code{sigma0_lb} and \code{sigma0_ub}.
#' @param sigma_prop_cov covariance matrix of the same length as the number of
#'                     lines that is used in sampling \eqn{\Sigma} values.
#'                     Defaults to \code{NULL}, meaning a diagonal matrix
#'                     is used. The elements on the diagonal of this matrix are
#'                     generally set to have value \code{sigma_ini}/20 unless
#'                     the corresponding observed y's have zero variance, in
#'                     which case these values are set to 0.1.
#' @param sigma_sp integer specifying how many elements in the \eqn{\Sigma}
#'                 matrix should be sampled together in the MCMC. Defaults to
#'                 25.
#' @param rho_ini double between 0 and 1 from which the value of \code{rho} will
#'                be initialized. Defaults to 0.5
#' @param rho0_lb double between 0 and 1 which gives the lower bound of the
#'               uniform prior for \code{rho}. Detauls to 0.
#' @param rho0_ub double between 0 and 1 which gives the upper bound of the
#'               uniform prior for \code{rho}. Defaults to 1.
#' @param rho_prop_sd standard deviation for the normal proposal distribution used
#'                  when proposing value for \code{rho} in the sampler. Defaults
#'                  to 0.01
#' @param w integer specifying how many samples of the parameters will be
#'          maintained. Samples from every w-th iteration is stored.
#' @return List that gives the values of the MCMC chain for
#'         \code{xU}, \code{mu}, \code{sigma} and \code{rho} along with
#'         indicators of acceptance on each iteration: \code{xURate},
#'         \code{sigmaRate}, and \code{rhoRate}. Background information
#'         is also outputted including the upper and lower bounds for
#'         unobserved x's (\code{xU_lb}, \code{xU_ub}), vectors giving
#'         the first and last indices of each grouping in sampling \eqn{\Sigma}
#'         (\code{sigma_ind_1},\code{sigma_ind_2}), the distance matrix
#'         (\code{dists}), and the integer specifying how many samples of the
#'         parameters will be maintained \code{w}
#' @importFrom stats cov qnorm sd
#' @examples \dontrun{
#' y_obs <- y_obs(maps = obs_maps, reg_info)
#' res <- fit_cont_pars(r = 3, n_iter = 1000, y_obs, reg_info)
#' }
fit_cont_pars <- function(r, n_iter, y_obs, reg_info,
                          dists = NULL, sigma_min = .01, sigma0_lb = NULL,
                          sigma0_ub = NULL, xU_prop_sd_def = .03, mu_ini = NULL,
                          mu0 = NULL, lambda0 = NULL, sigma_ini = NULL,
                          sigma_prop_cov = NULL, sigma_sp = 25, rho_ini = 0.5,
                          rho0_lb = 0, rho0_ub = 0.99, rho_prop_sd = 0.01,
                          w = 20) {
  #check that spacing of thinning works with number of iterations
  stopifnot(n_iter%%w == 0)

  n_lines <- nrow(y_obs[[r]])

  #information on hidden x's and bounds
  b_info <- bound_info(y_obs = y_obs[[r]], dist = reg_info$dist[[r]],
                       loop = reg_info$loop[r])
  b_ind <- which(b_info$xUnobs == 1, arr.ind = T)
  xU_lb <- b_info$lb[b_ind]; xU_ub <- b_info$ub[b_ind]
  xU_vecs <- b_ind[,1]; xU_years <- b_ind[,2]
  all_finite <- which(is.finite(xU_lb) & is.finite(xU_ub))
  xU_prop_sd <- rep(xU_prop_sd_def, length(xU_years))
  xU_prop_sd[all_finite] <- (xU_ub[all_finite] - xU_lb[all_finite])/150

  #information on observed covariance
  emp_cov_y <- cov(t(y_obs[[r]]))
  zero_var_ind <- which(diag(emp_cov_y) == 0)

  #information on max y lengths
  y_min <- unlist(lapply(reg_info$dist[[r]], min))
  y_max <- unlist(lapply(reg_info$dist[[r]], max))
  y_mid <- (y_min + y_max)/2

  #distance matrix
  if (is.null(dists)) {
    dists <- dist_mat(n_lines)
  }

  #priors
  if (is.null(sigma0_lb)) {
    sigma0_lb <- rep(sigma_min, n_lines)
  }
  if (is.null(sigma0_ub)) {
      sigma0_ub <- ((y_max - y_min)/2)/qnorm(.995)
  }
  if (is.null(mu0)) {
    mu0 <- y_mid
  }
  if (is.null(lambda0)) {
    lambda0 <- diag((((y_max - y_min)/2)/qnorm(.9))^2)
  }

  #Initial values
  if (is.null(mu_ini)) {
    mu_ini <- apply(y_obs[[r]], 1, mean)
  }
  if (is.null(sigma_ini)) {
    #generally start sigma at observed SD, but move slightly away from
    #values right at the boundaries
    sigma_ini <- sqrt(diag(emp_cov_y))
    high_sigma <- which(sigma_ini >= .95*sigma0_ub)
    high_sigma_adj <- .2*(sigma0_ub[high_sigma] - sigma0_lb[high_sigma])
    sigma_ini[high_sigma] <- sigma0_ub[high_sigma] - high_sigma_adj
    low_sigma <- which(sigma_ini <= .05*sigma0_lb)
    low_sigma_adg <- .2*(sigma0_ub[low_sigma] - sigma0_lb[low_sigma])
    sigma_ini[low_sigma] <- sigma0_lb[low_sigma] + low_sigma_adg
    #if no observed variability, just start at lower bound
    sigma_ini[which(sqrt(diag(emp_cov_y)) == 0)] <- sigma0_lb[which(sqrt(diag(emp_cov_y)) == 0)]
  }

  #MCMC parameters
  #sigma groupings for block updates
  cz <- contig_zero(diag(emp_cov_y))
  sigma_ind_1 <- sigma_ind_2 <- c()
  for (k in 1:nrow(cz)) {
    if (cz[k,"isAllZeros"]) {
      sigma_ind_1 <- c(sigma_ind_1, cz[k, "first"])
      sigma_ind_2 <- c(sigma_ind_2, cz[k, "last"])
    } else {
      nCurr <- cz[k, "last"] - cz[k, "first"] + 1
      if (nCurr > sigma_sp) {
        sigma_ind_1_curr <-  seq(cz[k, "first"], cz[k, "last"], sigma_sp)
        sigma_ind_1 <-  c(sigma_ind_1, sigma_ind_1_curr)
        sigma_ind_2 <- c(sigma_ind_2, sigma_ind_1_curr[2:(length(sigma_ind_1_curr))] - 1,
                       cz[k, "last"])
      } else {
        sigma_ind_1 <- c(sigma_ind_1, cz[k, "first"])
        sigma_ind_2 <- c(sigma_ind_2, cz[k, "last"])
      }
    }
  }
  if (is.null(sigma_prop_cov)) {
    temp <- sigma_ini/20
    temp[zero_var_ind] <- .1
    sigma_prop_cov <- diag(temp)
  }

  #Run MCMC
  res <- RunMCMC(n_iter = n_iter, dists = dists,
                 x = y_obs[[r]], xU_vecs = xU_vecs - 1,
                 xU_years = xU_years - 1, xU_prop_sd = xU_prop_sd,
                 xU_lb = xU_lb, xU_ub = xU_ub,
                 mu = mu_ini, mu0 = mu0, lambda0 = lambda0,
                 sigma = sigma_ini, sigma_ind_1 = sigma_ind_1 - 1,
                 sigma_ind_2 = sigma_ind_2 - 1, sigma_prop_cov = sigma_prop_cov,
                 rho = rho_ini, rho0_lb = rho0_lb, rho0_ub = rho0_ub,
                 rho_prop_sd = rho_prop_sd,
                 sigma0_lb = sigma0_lb, sigma0_ub = sigma0_ub, w = w)
  return(res)
}

# -------------------------------------------------------
# Related to setting up MCMC: priors, boundary info, etc.
# -------------------------------------------------------
#' Identify the indices of sequences of repeated zero values and indices
#' of sequences of non-zero values
#' @title Find indices of sequences of contiguous zeros
#' @param y vector to consider
#' @return Matrix of three columns where each row gives the first and last index
#'         of a sequence of numbers. The third column is a boolean. If
#'         \code{TRUE}, the indices are for a sequence zeros. If \code{FALSE},
#'         the indices are for a sequence non-zero values.
contig_zero <- function(y) {
  zero_Ind <- (y == 0)
  #index of switches between contiguous sets of zeros
  yZ <- which(zero_Ind == 1)
  y_switch <- which(diff(yZ) != 1) #jump index
  ind_Z <- matrix(ncol = 2, byrow = T,
                  data = sort(c(yZ[1], yZ[length(yZ)],
                                yZ[y_switch], yZ[y_switch + 1])))
  ind_Z <- cbind(ind_Z, rep(TRUE, nrow(ind_Z)))

  #index of switches between contiguous sets of non-zeros
  yNZ <- which(zero_Ind == 0)
  y_switch <- which(diff(yNZ) != 1) #jump index
  indNZ <- matrix(ncol = 2, byrow = T,
                  data = sort(c(yNZ[1], yNZ[length(yNZ)],
                                yNZ[y_switch], yNZ[y_switch + 1])))
  indNZ <- cbind(indNZ, rep(FALSE, nrow(indNZ)))

  #put together and order
  ind <- rbind(ind_Z, indNZ)
  ind <- matrix(data = ind[order(ind[,1]),], ncol = 3)
  colnames(ind) <- c("first", "last", "isAllZeros")

  return(ind)
}

#' Determine which y values are on the boundaries and what the corresponding
#' bounds of those y values are
#' @title Get boundary info
#' @param y_obs matrix of observed distances (dimension: number of lines
#'             by number of years)
#' @param dist a list of the lengths for the corresponding \code{lines}
#' @param loop boolean which if true \code{TRUE} indicates that the \code{lines}
#'             extend outward from a single point forming a circle and if
#'             \code{FALSE} indicates that the lines are mapped along a fixed
#'             contour such as a land boundary
#' @return list of 3 matrices each of dimension number of lines by number of
#'         years giving the lower bound for hidden x values, the upper bound
#'         for hidden x values, and an indicator of whether the
#'         value is bounded at all. The values in the list are named
#'         \code{ub}, \code{lb}, and \code{xUnObs} respectively.
#'
bound_info <- function(y_obs, dist, loop) {
  #check if length of y's match length of distances
  stopifnot(nrow(y_obs) == length(dist))

  #if only 1 year of data in a vector form, convert it to a matrix
  if (is.vector(y_obs)) {
    y_obs <- matrix(data = y_obs, nrow = length(y_obs), ncol = 1)
  }
  n_lines <- nrow(y_obs); n_years <- ncol(y_obs)
  #matrix indicating whether an x value is not observed and needs to be imputed
  x_unobs <- matrix(nrow = n_lines, ncol = n_years, data = 0)
  if (loop) { #if loop, all y-values with value 0 are equal to x
    lb <- matrix(nrow = n_lines, ncol = n_years, data = 0)
  } else {  #else, zero-length y values can come from negative length x-values
    lb <- matrix(nrow = n_lines, ncol = n_years, data = -Inf)
  }
  ub <- matrix(nrow = n_lines, ncol = n_years, data = Inf)

  #if y is on a boundary, identify the lower and upper bounds of corresponding x
  for (i in 1:n_lines) {
    for (j in 1:length(dist[[i]])) {
      if ((!loop) || (j != 1)) { #if loop and touches first bound, ub unaffected
        #determine if on (or very near) boundary
        on_b_ind <- which(abs(y_obs[i,] - dist[[i]][j]) <= 1e-5)
        if (j == 1) {
          #Add unobserved indices for the special case where length zero
          #observed (no ice), but a non-zero length would be required if ice
          #were to  grow (occurs when line starts out by going through land)
          on_b_ind <- c(on_b_ind, which(y_obs[i,] < dist[[i]][j]))
        }
        if (length(on_b_ind) > 0) {
          x_unobs[i, on_b_ind] <- 1 #record that value will need to be imputed
          if (j == 1) {
            ub[i, on_b_ind] <- 0
          } else if (j == length(dist[[i]])) { #touches last length
            lb[i, on_b_ind] <- dist[[i]][j]
          } else if (j%%2 == 0) { #touches even-valued boundary
            lb[i, on_b_ind] <- dist[[i]][[j]]
            ub[i, on_b_ind] <- dist[[i]][[j]] +
                             (dist[[i]][[j + 1]] - dist[[i]][[j]])/2
          } else { #touches odd-valued boundary
            lb[i, on_b_ind] <- dist[[i]][[j]] -
                             (dist[[i]][[j]] - dist[[i]][[j - 1]])/2
            ub[i, on_b_ind] <- dist[[i]][[j]]
          }
        }
      }
    }
  }

  #For observed x-values, the lower and upper bound are just the observed value
  #(no censoring to account for)
  lb[x_unobs == 0] <- y_obs[x_unobs == 0]
  ub[x_unobs == 0] <- y_obs[x_unobs == 0]

  return(list("lb" = lb, "ub" = ub, "xUnobs" = x_unobs))
}



#' Creates a matrix specifying how difference indices are among an
#' ordered set of indices. For lines with indices \eqn{i} and \eqn{j},
#' the 'distance' computed is \eqn{|i -j|}. Results are used to define covariance.
#' @title Compute 'distances' among $n$ lines
#' @param n_lines number of lines in matrix
#' @return Matrix of `distances' among the indices
dist_mat <- function(n_lines) {
  dists <- matrix(nrow = n_lines, ncol = n_lines)
  #compute distance as sequences of locations along a path
  for (i in 1:n_lines) {
    for (j in 1:i) {
      dists[i, j] <- dists[j, i] <- abs(i - j)
    }
  }
  return(dists)
}

# ---------------------------------------------
# related to finding y and info from it
# ----------------------------------------------

#' Identify the years on which to train accounting for years with missing data
#' @title Find indices on which to train contour model
#' @param maps output of a \code{create_mapping} object
train_ind <- function(maps) {
  train_ind <- trainInd <- 1:length(maps$start_year:maps$end_year)
  na_only <- sapply(maps$obs_list, function(x){apply(x, 1, function(y){all(is.na(y))})})
  missing <- which(apply(na_only, 1, function(x){all(x)}))
  if (length(missing) > 0) {
    missing_ind <- which(train_ind == missing)
    if (length(missing_ind) > 0) {
      train_ind <- train_ind[-missing_ind]
    }
  }
  return(train_ind)
}


#' Compute y values from the output of the \code{create_mapping} object
#' @title Compute \eqn{y}
#' @param maps output of the \code{create_mapping} function
#' @param reg_info  a \code{reg_info} list (see documentation for \code{reg_info})
#' @return List of matrices, one per region, giving the observed \eqn{y} values.
#'         Each row corresponds to the lines in \eqn{L} and each column corresponds
#'         to  a training year
#' @export
#' @examples
#' y_obs <- y_obs(maps = obs_maps, reg_info)
y_obs <- function(maps, reg_info) {
  train_ind <- train_ind(maps)
  n_reg <- length(reg_info$regions)
  n_train <- length(train_ind)
  y_obs <- list()
  for (r in 1:n_reg) {
    n_lines_r <- nrow(reg_info$start_coords[[r]])
    y_obs_r <- matrix(nrow = n_lines_r, ncol = n_train)
    for (l in 1:n_train) {
      y_obs_r[,l] <- sqrt(maps$obs_list[[r]][train_ind[l],,5]^2
                          + maps$obs_list[[r]][train_ind[l],,6]^2)
    }
    y_obs[[r]] <- y_obs_r
  }
  return(y_obs)
}

#' Determine which regions are completely ice-filled (full) or ice-free (empty)
#' in all years in the training period. Also, make the polygon corresponding to
#' regions that are fully ice-covered.
#' @title Identify fully ice-covered and ice-free regions
#' @param y_obs list of y values outputted from \code{y_obs} function
#' @param reg_info  a \code{reg_info} list (see documentation for \code{reg_info})
#' @importFrom maptools spRbind
#' @importFrom raster aggregate
to_fit <- function(y_obs, reg_info) {
  n_reg <- length(reg_info$regions)
  empty <- full <- rep(FALSE, n_reg)

  #determine which regions are completely ice-filled (full) or ice-free in all
  #years in training period
  for (r in 1:n_reg) {
    if (all((y_obs[[r]] == 0))) {
      empty[r] <- TRUE
    } else if ((all(apply(y_obs[[r]], 1, sd) == 0))) {
      full[r] <- TRUE
    }
  }
  regs_to_fit <- which(!empty & !full)
  reg_full <- which(full)

  #Make a polygon of all regions that are fully filled with sea ice
  full <- NA
  nRegsFull <- length(reg_full)
  if (length(reg_full) > 0) {
    full <- reg_info$regions[[reg_full[[1]]]]
    full@polygons[[1]]@ID <- "regsFull"
    if (length(reg_full) > 1) {
      for (i in 2:nRegsFull) {
        full <- aggregate(spRbind(full, reg_info$regions[[reg_full[[i]]]]))
        full@polygons[[1]]@ID <- "regsFull"
      }
    }
    full <- gDifference(full, land)
  }

  return(list("regs_to_fit" = regs_to_fit, "full" = full))
}

# ----------------------------------------------
# computing parameters from MCMC results
# ----------------------------------------------
#' Compute parameter estimates for the contour model using MCMC output from
#' \code{fit_cont_pars}
#' @title Compute Parameters Estimates
#' @param res_r output of MCMC run from function \code{fit_cont_pars} for
#'              one region
#' @param burn_in number of iterations to discard as burn-in. This is the number
#'               before thinning. Value will be divided by \code{w}.
#' @param w integer specifying the thinning used. Samples from every w-th
#'          iteration are stored.
#' @return List of a list of parameters for each region. Each list contains two
#'         elements, \code{muEst} and \code{sigmaEst}. These which give
#'         estimates for the \code{mu} and \code{sigma} parameters used to
#'         generate contours.
#' @export
#' @examples \dontrun{
#' y_obs <- y_obs(maps = obs_maps, reg_info)
#' res <- fit_cont_pars(r = 3, n_iter = 1000, y_obs, reg_info)
#' calc_pars(res, burn_in = 100, w = res$w)
#' }
calc_pars <- function(res_r, burn_in, w) {
  burn_in <- round(burn_in/w)
  n_samp <- length(res_r$rho)
  mu_est <- apply(res_r$mu[,(burn_in + 1):n_samp], 1, mean)
  n_lines <- length(mu_est)
  sigma_est_vec <- apply(res_r$sigma[,(burn_in + 1):n_samp], 1, mean)
  rho_est <- mean(res_r$rho[(burn_in + 1):n_samp])
  sigma_est <- matrix(nrow = n_lines, ncol = n_lines, data = NA)
  for (i in 1:n_lines) {
    for (j in 1:i) {
      sigma_est[i, j] <- sigma_est[j, i] <- (sigma_est_vec[i]*sigma_est_vec[j]*
                                             rho_est^(res_r$dists[i, j]))
    }
  }
  return(list("mu_est" = mu_est, "sigma_est" = sigma_est))
}
