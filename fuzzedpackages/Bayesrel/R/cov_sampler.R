# this function uses gibbs sampling to estimate the posterior distribution
# of a sample's covariance matrix
# sources: https://en.wikipedia.org/wiki/Normal-inverse-Wishart_distribution,
# Murphy, K. P. (2007). Conjugate bayesian analysis of the gaussian distribution (Tech. Rep.). University of British Columbia.

covSamp <- function(data, n.iter, n.burnin, thin, n.chains, pairwise, callback = function(){}){
  n <- nrow(data)
  p <- ncol(data)

  c_post <- array(0, c(n.chains, n.iter, p, p))
  inds <- which(is.na(data), arr.ind = T)
  dat_imp <- array(0, c(n.chains, n.iter, nrow(inds)))

  for (z in 1:n.chains) {

    if (pairwise) {
      dat_complete <- data
      # initial generation of complete data set with means as substitutes
      dat_complete[inds] <- colMeans(data, na.rm = T)[inds[, 2]]
      # now the missing are being replaced in each iteration with draws from the conditional joints
      for (i in 1:n.iter) {
        cc <- sampleCovParams(dat_complete)
        # ms <- MASS::mvrnorm(1, mun, cc/kn)
        ms <- numeric(p)
        c_post[z, i, , ] <- cc
        # substitute missing values one by one, where each value is drawn conditional on the rest of the data
        cols <- unique(inds[, 2])
        for (ccc in cols) {
          rows <- inds[which(inds[, 2] == ccc), 1]
          mu1 <- ms[ccc]
          mu2 <- ms[-ccc]
          cc11 <- cc[ccc, ccc]
          cc21 <- cc[-ccc, ccc]
          cc12 <- cc[ccc, -ccc]
          cc22 <- cc[-ccc, -ccc]
          ccq <- cc11 - cc12 %*% try(solve(cc22)) %*% cc21
          for (r in rows) {
            muq <- mu1 + cc12 %*% try(solve(cc22)) %*% (as.numeric(dat_complete[r, -ccc]) - mu2)
            dat_complete[r, ccc] <- rnorm(1, muq, sqrt(ccq))
          }
        }
        callback()
        dat_imp[z, i, ] <- dat_complete[inds]
      }

    } else {
      k0 <- 1e-10
      v0 <- p
      t <- diag(p)
      T0 <- diag(k0, nrow = p, ncol = p) # matrix inversion of diagonal matrix
      mu0 <- rep(0, p) # prior means
      kn <- k0 + n
      vn <- v0 + n

      ym <- .colMeans(data, n, p)
      mun <- (k0 * mu0 + n * ym) / (k0 + n)
      S <- cov(sweep(data, 2L, ym, `-`)) * (n - 1)

      Tn <- T0 + S + (k0 * n / (k0 + n)) * (ym - mu0) %*% t(ym - mu0)
      # drawing samples from posterior:
      Tn <- chol(chol2inv(chol(Tn)))
      dfChisq <- vn:(vn-p+1)
      utz <- upper.tri(matrix(0, p, p))
      for (i in 1:n.iter){
        c_post[z, i, , ] <- rinvwishart2(vn, Tn, p, dfChisq, utz) # sample from inverse Wishart
        callback()
      }
    }
  }

  c_post_burned <- c_post[, (n.burnin + 1):n.iter, , , drop = F]
  c_post_out <- c_post_burned[, seq(1, dim(c_post_burned)[2], thin), , , drop = F]

  dat_imp_burned <- dat_imp[, (n.burnin + 1):n.iter, , drop = F]
  dat_out <- dat_imp_burned[, seq(1, dim(dat_imp_burned)[2], thin), , drop = F]


  return(list(cov_mat = c_post_out, dat_mis_samp_cov = coda::mcmc(dat_out)))
}

sampleCovParams <- function(data) {
  n <- nrow(data)
  p <- ncol(data)
  # posterior covariance matrix ---------------------------------------------------
  k0 <- 1e-10
  v0 <- p
  t <- diag(p)
  T0 <- diag(k0, nrow = p, ncol = p) # matrix inversion of diagonal matrix
  mu0 <- rep(0, p) # prior means
  kn <- k0 + n
  vn <- v0 + n

  ym <- .colMeans(data, n, p)
  mun <- (k0 * mu0 + n * ym) / (k0 + n)
  S <- cov(sweep(data, 2L, ym, `-`)) * (n - 1)

  Tn <- T0 + S + (k0 * n / (k0 + n)) * (ym - mu0) %*% t(ym - mu0)
  # drawing samples from posterior:
  Tn <- chol(chol2inv(chol(Tn)))
  dfChisq <- vn:(vn-p+1)
  utz <- upper.tri(matrix(0, p, p))
  cc <- rinvwishart2(vn, Tn, p, dfChisq, utz) # sample from inverse Wishart

  return(cc)
}

# ------- customized covariance matrix sampling with cholesky decomposition -----------
rinvwishart2 <- function(nu, S, k = nrow(S), dfChisq = nu:(nu-k+1), utz = upper.tri(matrix(0, k, k))) {

  # LaplacesDemon::rwishartc
  Z <- matrix(0, k, k)
  x <- rchisq(k, dfChisq)
  x[x == 0] <- 1e-100
  diag(Z) <- sqrt(x)
  if (k > 1) {
    # kseq <- 1:(k - 1)
    # Z[rep(k * kseq, kseq) + unlist(lapply(kseq, seq))] <- rnorm(k * {k - 1} / 2)
    # --end of copied code
    Z[utz] <- rnorm(k * {k - 1} / 2)
  }
  # LaplacesDemon::rinvwishart
  return(chol2inv(Z %*% S))
}
