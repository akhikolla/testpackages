new_outlier_model <- function(A, sims, mu, sigma, cdf, cms, sms) {
  structure(
    list(
      A     = A,
      sims  = sims,
      mu    = mu,
      sigma = sigma,
      cdf   = cdf,
      cms   = cms,
      sms   = sms
    ),
    class = c("outlier_model", "list")
  )
}

new_novelty <- function(m, dev, pv, cv, a) {
  # m : outlier_model object
  m$dev    <- dev
  m$pval   <- pv
  m$cv     <- cv
  m$alpha  <- a
  class(m) <- c("novelty", class(m))
  return(m)
}

new_outlier <- function(m, cv, a) {
  # m : outlier_model object
  m$cv     <- cv
  m$alpha  <- a
  class(m) <- c("outlier", class(m))
  return(m)
}

new_mixed_outlier <- function(m, type, alpha, dev = NULL) {
  # m : outlier_model object

  out <- if (type == "novelty") {
    new_novelty(m, dev, pval(m, dev), critval(m, alpha), alpha)
  } else {
    new_outlier(m, critval(m, alpha), alpha)
  }

  class(out) <- c("mixed_outlier", class(out))
  return(out)
}

convolute <- function(m1, m2) {
  # m1 and m2 : outlier_model objects
  .sims <- m1$sims + m2$sims
  .cdf  <- stats::ecdf(.sims)
  .mu   <- m1$mu + m2$mu
  .sig  <- m1$sigma + m2$sigma
  new_outlier_model(
    A = list(A1 = m1$A, A2 = m2$A),
    .sims,
    .mu,
    .sig,
    .cdf,
    cms = list(cms1 = m1$cms, cms2 = m2$cms),
    sms = list(sms1 = m1$sms, sms2 = m2$sms)
  )
}
