# ============================= gp_loglik =====================================

gp_loglik <- function(pars, data, m, xm, sum_gp) {
  #
  # Evaluates the GP log-likelihood at pars = (sigma_u,xi) based on data data.
  # For values of xi close to zero the value is approximated using a Taylor
  # series expansion.
  #
  # Args:
  #   pars   : numeric vector (sigma_u, xi).
  #   data   : sample of threshold excesses.
  #   m      : sample size.
  #   xm     : largest sample excess.
  #   sum_gp : sum of sample excesses.
  #
  # Returns:
  #   the value of the log-likelihood.
  #
  if (pars[1] <= 0 | pars[2] <= -pars[1] / xm) {
    return(-Inf)
  }
  sdat <- data / pars[1]
  zz <- 1 + pars[2] * sdat
  if (abs(pars[2]) > 1e-6) {
    val <- -m * log(pars[1]) - (1 + 1 / pars[2]) * sum(log(zz))
  } else {
    i <- 1:4
    g_fun <- function(x) {
      t1 <- x ^ i
      t2 <- (i * x - i - 1)
      sum((-1) ^ i * t1 * t2 * pars[2] ^ i / i / (i + 1))
    }
    val <- -m * log(pars[1]) - sum_gp / pars[1] - sum(sapply(sdat, g_fun))
  }
  return(val)
}

# ============================= gev_loglik =====================================

gev_loglik <- function(pars, data, m, x1 = NULL, xm = NULL, sum_gev) {
  #
  # Evaluates the GEV log-likelihood at pars = (mu, sigma, xi) based on data
  # data.  For values of xi close to zero the value is approximated using
  # Taylor series expansions.
  #
  # Args:
  #   pars    : numeric vector (mu, sigma, xi).
  #   data    : sample of block maxima.
  #   m       : sample size.
  #   x1      : the smallest block maximum. (Not used in gev_loglik.)
  #   xm      : the largest block maximum. (Not used in gev_loglik.)
  #   sum_gev : sum of all block maxima.
  #
  # Returns:
  #   the value of the log-likelihood.
  #
  if (pars[2] <= 0) {
    return(-Inf)
  }
  sdat <- (data - pars[1]) / pars[2]
  zz <- 1 + pars[3] * sdat
  if (any(zz <= 0)) {
    return(-Inf)
  }
  val <- -m * log(pars[2])
  if (abs(pars[3]) > 1e-6) {
    val <- val - (1 + 1 / pars[3]) * sum(log(zz)) - sum(zz ^ (- 1 / pars[3]))
  } else {
    t0 <- (sum_gev - m * pars[1]) / pars[2]
    i <- 1:4
    g_fun <- function(x) {
      t1 <- x ^ i
      t2 <- (i * x - i - 1)
      sum((-1) ^ i * t1 * t2 * pars[3] ^ i / i / (i + 1))
    }
    g2_fun <- function(x) {
      sum((-1) ^ i * x ^ (i + 1) * pars[3] ^ i / (i + 1))
    }
    val <- val - t0 - sum(sapply(sdat, g_fun)) -
       sum(exp(-(data - pars[1]) / pars[2] - sapply(sdat, g2_fun)))
  }
  return(val)
}

# ============================= os_loglik =====================================

os_loglik <- function(pars, data, min_data, max_data = NULL, nos, x1 = NULL,
                       xm = NULL, sum_os, gumbel = FALSE) {
  #
  # Evaluates the GEV log-likelihood at pars = (mu, sigma, xi) based on data
  # data.  For values of xi close to zero the value is approximated using
  # Taylor series expansions.
  #
  # Args:
  #   pars     : numeric vector (mu, sigma, xi).
  #   data     : sample of all order statistics.
  #   min_data : sample of smallest order statistics in the data.
  #   max_data : sample of largest order statistics in the data. (Not used in
  #              os_loglik.)
  #   nos      : total number of order statistics in the data.
  #   x1       : the smallest order statistic. (Not used in os_loglik.)
  #   xm       : the largest order statistic. (Not used in os_loglik.)
  #   sum_os   : sum of all the order statistics.
  #
  # Returns:
  #   the value of the log-likelihood.
  #
  if (gumbel) {
    pars <- c(pars, 0)
  }
  if (pars[2] <= 0) {
    return(-Inf)
  }
  sdat <- (data - pars[1]) / pars[2]
  zz <- 1 + pars[3] * sdat
  if (any(zz <= 0)) {
    return(-Inf)
  }
  smindat <- (min_data - pars[1]) / pars[2]
  zz_min <- 1 + pars[3] * smindat
  val <- -nos * log(pars[2])
  if (abs(pars[3]) > 1e-6) {
    val <- val - (1 + 1 / pars[3]) * sum(log(zz)) -
      sum(zz_min ^ (- 1 / pars[3]))
  } else {
    t0 <- (sum_os - nos * pars[1]) / pars[2]
    i <- 1:4
    g_fun <- function(x) {
      t1 <- x ^ i
      t2 <- (i * x - i - 1)
      sum((-1) ^ i * t1 * t2 * pars[3] ^ i / i / (i + 1))
    }
    g2_fun <- function(x) {
      sum((-1) ^ i * x ^ (i + 1) * pars[3] ^ i / (i + 1))
    }
    val <- val - t0 - sum(sapply(sdat, g_fun)) -
      sum(exp(-(min_data - pars[1]) / pars[2] - sapply(smindat, g2_fun)))
  }
  return(val)
}

# ============================== pp_loglik ====================================

pp_loglik <- function(pars, data, n_exc, thresh, xm = NULL, noy, sum_pp,
                      m = NULL) {
  #
  # Evaluates the PP log-likelihood at pars = (mu,sigma,xi) based on data data.
  # For values of xi close to zero the value is approximated using Taylor
  # series expansions.
  #
  # Args:
  #   pars      : numeric vector (mu, sigma, xi).
  #   data      : sample of threshold exceedances.
  #   n_exc     : number of threshold excesses.
  #   thresh    : threshold.
  #   xm        : the maximum data value. (Not used in pp_loglik.)
  #   noy       : number of years (or periods) of data.
  #   sum_pp    : sum of the threshold exceedances.
  #   m         : the number of raw data values. (Not used in pp_loglik.)
  #
  # Returns:
  #   the value of the log-likelihood.
  #
  if (pars[2] <= 0) {
    return(-Inf)
  }
  udat <- (thresh - pars[1]) / pars[2]
  zz_u <- 1 + pars[3] * udat
  if (zz_u <= 0) {
    return(-Inf)
  }
  sdat <- (data - pars[1]) / pars[2]
  zz <- 1 + pars[3] * sdat
  if (any(zz <= 0)) {
    return(-Inf)
  }
  val <- -n_exc * log(pars[2])
  if (abs(pars[3]) > 1e-6) {
    val <- val - (1 + 1 / pars[3]) * sum(log(zz)) -
      noy * zz_u ^ (- 1 / pars[3])
  } else {
    t0 <- (sum_pp - n_exc * pars[1]) / pars[2]
    i <- 1:4
    g_fun <- function(x) {
      t1 <- x ^ i
      t2 <- (i * x - i - 1)
      sum((-1) ^ i * t1 * t2 * pars[3] ^ i / i / (i + 1))
    }
    g2_fun <- function(x) {
      sum((-1) ^ i * x ^ (i + 1) * pars[3] ^ i / (i + 1))
    }
    val <- val - t0 - sum(sapply(sdat, g_fun)) -
      noy * exp(-(thresh - pars[1]) / pars[2] - g2_fun(udat))
  }
  return(val)
}
