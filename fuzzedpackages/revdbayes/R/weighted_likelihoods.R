# ============================= gp_loglik =====================================

gp_wloglik <- function(pars, data, m, xm, sum_gp, w, sumw) {
  #
  # Evaluates the weighted GP log-likelihood at pars = (sigma_u,xi) based on
  # data data.  For values of xi close to zero the value is approximated using
  # a Taylor series expansion.
  #
  # Args:
  #   pars   : numeric vector (sigma_u, xi).
  #   data   : sample of threshold excesses.
  #   m      : sample size.
  #   xm     : largest sample excess.
  #   sum_gp : sum of sample excesses.
  #   w      : numeric vector of the same length as data.
  #   sumw   : the sum of the values in w.
  #
  # Returns:
  #   the value of the weighted log-likelihood.
  #
  if (pars[1] <= 0 | pars[2] <= -pars[1] / xm) {
    return(-Inf)
  }
  sdat <- data / pars[1]
  zz <- 1 + pars[2] * sdat
  if (abs(pars[2]) > 1e-6) {
    val <- -sumw * log(pars[1]) - (1 + 1 / pars[2]) * sum(w * log(zz))
  } else {
    i <- 1:4
    g_fun <- function(x) {
      t1 <- x ^ i
      t2 <- (i * x - i - 1)
      sum((-1) ^ i * t1 * t2 * pars[2] ^ i / i / (i + 1))
    }
    val <- -sumw * log(pars[1]) - sum(w * data) / pars[1] -
      sum(w * sapply(sdat, g_fun))
  }
  return(val)
}
