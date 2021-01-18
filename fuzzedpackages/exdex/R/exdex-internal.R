#' Internal exdex functions
#'
#' Internal exdex functions
#' @details
#' These functions are not intended to be called by the user.
#' @name exdex-internal
#' @keywords internal
NULL

# =============================== To check spm() ============================ #

#' @keywords internal
#' @rdname exdex-internal
spm_check <- function(data, b, sliding = TRUE,
                      bias_adjust = c("BB3", "BB1", "N", "none"),
                      constrain = TRUE, varN = TRUE,
                      which_dj = c("last", "first")) {
  Call <- match.call(expand.dots = TRUE)
  bias_adjust <- match.arg(bias_adjust)
  # We don't check inputs here because this function is only used for
  # testing cases where I have ensured that the input are OK
  which_dj <- match.arg(which_dj)
  # Check that the value of b satisfies the inequality in Proposition 4.1
  # of Berghaus and Bucher (2018)
  k_n <- floor(length(data) / b)
  if (k_n < 1) {
    stop("b is too large: it is larger than length(data)")
  }
  # Assume that the value of b is OK (for variances to be estimated)
  b_ok <- TRUE
  # If we want the last set of disjoint maxima then reverse the data vector
  if (which_dj == "last") {
    pass_data <- rev(data)
  } else {
    pass_data <- data
  }
  # A function that returns N2015 and BB2018 estimates of theta
  # (constrained to (0, 1] if constrain = TRUE))
  spm_estimates <- function(data) {
    # Calculate the block maxima
    if (sliding) {
      temp <- sliding_maxima(data, b)
    } else{
      temp <- disjoint_maxima(data, b)
    }
    # Extract x ~ F (only xs contributing to y are included) and y ~ G
    x <- temp$x
    y <- temp$y
    # Empirical c.d.f. of raw (`daily') values
    Fhat <- stats::ecdf(x)
    # Evaluate Fx at y
    Fhaty <- Fhat(y)
    if (bias_adjust == "N") {
      # `Bias-adjust' the empirical c.d.f. of Y based on the Xs: by subtracting
      # b in numerator and denominator we remove Xs that are in the same block
      # as Y. We use m-b in the denominator rather than the m-b+1 in
      # Northrop (2015)
      m <- length(x)
      Fhaty <- (m * Fhaty - b) / (m - b)
      # In the unlikely event that an element of Fhaty is equal to zero,
      # i.e. a block maximum is less than all the data from outside that
      # block, we force Fhaty to be positive
      Fhaty[Fhaty == 0] <- 1 / (m - b + length(y))
    }
    # Calculate the estimate of theta:
    # theta_N: Northrop (2015) and theta_BB: Berghaus and Bucher (2018)
    theta_N <- -1 / mean(b * log(Fhaty))
    theta_BB <- 1 / (b * mean(1 - Fhaty))
    # Estimate sigma2_dj based on Section 4 of Berghaus and Bucher (2018)
    # We require the disjoint maxima to do this.
    # Only do this if b_ok = TRUE
    if (b_ok) {
      if (sliding) {
        sigma2hat_dj <- spm_sigmahat_dj(data = data, b = b)
      } else {
        sigma2hat_dj <- spm_sigmahat_dj(data = data, b = b)
      }
      # If sliding = TRUE then estimate sigma2hat_sl
      # Otherwise use sigma2hat_dj
      indexN <- ifelse(varN, 2, 1)
      if (sliding) {
        sigma2hat_N <- sigma2hat_dj[indexN] - (3 - 4 * log(2)) / theta_N ^ 2
        sigma2hat_BB <- sigma2hat_dj[1] - (3 - 4 * log(2)) / theta_BB ^ 2
      } else {
        sigma2hat_N <- sigma2hat_dj[indexN + 2]
        sigma2hat_BB <- sigma2hat_dj[1 + 2]
      }
    }
    # Estimate the sampling variances of the estimators
    theta <- c(theta_N, theta_BB)
    if (b_ok) {
      vars <- theta ^ 4 * c(sigma2hat_N, sigma2hat_BB) / k_n
    }
    # Perform BB2018 bias-adjustment if required
    bias_N <- bias_BB <- 0
    if (bias_adjust == "BB3") {
      bias_N <- theta_N / k_n + theta_N ^ 3 * sigma2hat_N / k_n
      theta_N <- theta_N * (1 - 1 / k_n) - theta_N ^ 3 * sigma2hat_N / k_n
      bias_BB <- theta_BB / k_n + theta_BB ^ 3 * sigma2hat_BB / k_n
      theta_BB <- theta_BB * (1 - 1 / k_n) - theta_BB ^ 3 * sigma2hat_BB / k_n
      theta <- c(theta_N, theta_BB)
    } else if (bias_adjust == "BB1") {
      bias_N <- theta_N / k_n
      theta_N <- theta_N * (1 - 1 / k_n)
      bias_BB <- theta_BB / k_n
      theta_BB <- theta_BB * (1 - 1 / k_n)
      theta <- c(theta_N, theta_BB)
    }
    # Save the unconstrained estimates, so that they can be returned
    unconstrained_theta <- theta
    # Constrain to (0, 1] if required
    if (constrain) {
      theta <- pmin(theta, 1)
    }
    se <- sqrt(vars)
    return(list(theta = theta, se = se,
                unconstrained_theta = unconstrained_theta,
                N2015_data = -b * log(Fhaty),
                BB2018_data = b * (1 - Fhaty),
                bias_val = c(bias_N, bias_BB)))
  }
  # End of function spm_estimates() ----------
  #
  # Find the point estimate of theta and the raw data that contribute to it
  res <- spm_estimates(data = pass_data)
  estimator_names <- c("N2015", "BB2018")
  names(res$theta) <- estimator_names
  names(res$se) <- estimator_names
  names(res$unconstrained_theta) <- estimator_names
  names(res$bias_val) <- estimator_names
  res$bias_adjust <- bias_adjust
  res$b <- b
  res$sliding <- sliding
  res$call <- Call
  class(res) <- c("spm", "exdex")
  return(res)
}

#' @keywords internal
#' @rdname exdex-internal
spm_sigmahat_dj <- function(data, b, check = FALSE){
  all_dj_maxima <- all_disjoint_maxima(data, b)
  # The number of blocks and the number of raw observations that contribute
  k_n <- nrow(all_dj_maxima$y)
  m <- nrow(all_dj_maxima$x)
  lenxx <- m - b
  const <- -log(m - b + k_n)
  # block indicator
  block <- rep(1:k_n, each = b)
  # Set up some functions
  #  BB2018_fn <- function(x, y) {
  #    return(mean(1 - ecdf2(x, y)))
  #  }
  BB2018_fn <- function(x, y) {
    return(1 - sum(ecdf2(x, y)) / k_n)
  }
  #  loobBB2018_fn <- function(x, block, xvec, y) {
  #    xx <- xvec[block != x]
  #    return(mean(1 - ecdf2(xx, y)))
  #  }
  loobBB2018_fn <- function(x, block, xvec, y) {
    xx <- xvec[block != x]
    return(1 - sum(ecdf2(xx, y)) / k_n)
  }
  # In the unlikely event that an element of Fhaty is equal to zero,
  # i.e. a block maximum is less than all the data from outside that
  # block, we force Fhaty to be positive
  #  loobN2015_fn <- function(x, block, xvec, y) {
  #    xx <- xvec[block != x]
  #    return(mean(-log0const(ecdf2(xx, y), const)))
  #  }
  loobN2015_fn <- function(x, block, xvec, y) {
    xx <- xvec[block != x]
    return(sum(-log0const(ecdf2(xx, y), const)) / k_n)
  }
  ests_fn <- function(i) {
    # y: disjoint block maxima, x: (only) the raw values that contribute to y
    x <- all_dj_maxima$x[, i]
    y <- all_dj_maxima$y[, i]
    # The ecdf of the data evaluated at the block maxima
    Nhat <- ecdf2(x, y)
    # BB2018
    Zhat <- b * (1 - Nhat)
    That <- mean(Zhat)
    # N2015
    ZhatN <- -b * log(Nhat)
    ThatN <- mean(ZhatN)
    Usum <- b * tapply(x, block, BB2018_fn, y = y)
    UsumN <- b * vapply(1:k_n, loobN2015_fn, 0, block = block, xvec = x, y = y)
    UsumN <-  k_n * ThatN - (k_n - 1) * UsumN
    #
    Bhat <- Zhat + Usum - 2 * That
    BhatN <- ZhatN + UsumN - 2 * ThatN
    # Bhat is mean-centred, but BhatN isn't (quite)
    BhatN <- BhatN - mean(BhatN)
    sigmahat2_dj <- mean(Bhat ^ 2)
    sigmahat2_djN <- mean(BhatN ^ 2)
    return(c(sigmahat2_dj, sigmahat2_djN))
  }
  if (!check) {
    temp <- vapply(1:ncol(all_dj_maxima$y), ests_fn, c(0, 0))
    ests <- rowMeans(temp)
    # For dj maxima add the first values.  This is the case because, if
    # which_dj = "last" we used rev(data) to reverse the data
    ests <- c(ests, temp[, 1])
    return(ests)
  }
  # 4 (effectively 3) other ways to calculate Usum
  # (Usum1 is commented out because it incurs rounding error from
  # Nhat -> Zhat -> Nhat that can be non-negligible)
  x <- all_dj_maxima$x[, 1]
  y <- all_dj_maxima$y[, 1]
  #  Fhat <- stats::ecdf(x)
  #  Nhat <- Fhat(y)
  Nhat <- ecdf2(x, y)
  Zhat <- b * (1 - Nhat)
  That <- mean(Zhat)
  Usum <- b * tapply(x, block, BB2018_fn, y = y)
  #  Uhats <- Fhat(x)
  Uhats <- ecdf2(x, x)
  # Usum1 <- colSums(vapply(Uhats, function(x) x > 1 - Zhat / b, rep(0, k_n)))
  Usum2 <- colSums(vapply(Uhats, function(x) x > Nhat, rep(0, k_n)))
  Usum3 <- colSums(vapply(x, function(x) x > y, rep(0, k_n)))
  # Usum4 is the analogous calculation to UsumN
  Usum4 <- b * vapply(1:k_n, loobBB2018_fn, 0, block = block, xvec = x, y = y)
  Usum4 <-  k_n * That - (k_n - 1) * Usum4
  # Aggregate the first 3
  # Usum1 <- tapply(Usum1, block, sum) / k_n
  Usum2 <- tapply(Usum2, block, sum) / k_n
  Usum3 <- tapply(Usum3, block, sum) / k_n
  return(cbind(Usum, Usum2, Usum3, Usum4))
}

# ============================== spm_R_quick() ============================== #

#' @keywords internal
#' @rdname exdex-internal
spm_R_quick <- function(data, b, bias_adjust = c("BB3", "BB1", "N", "none"),
                        constrain = TRUE, varN = TRUE,
                        which_dj = c("last", "first")) {
  Call <- match.call(expand.dots = TRUE)
  #
  # Check inputs
  #
  if (missing(data) || length(data) == 0L || mode(data) != "numeric") {
    stop("'data' must be a non-empty numeric vector")
  }
  if (any(!is.finite(data))) {
    stop("'data' contains missing or infinite values")
  }
  if (is.matrix(data)) {
    stop("'data' must be a vector")
  }
  data <- as.vector(data)
  if (!is.numeric(b) || length(b) != 1) {
    stop("'b' must be a numeric scalar (specifically, a positive integer)")
  }
  if (b < 1) {
    stop("'b' must be no smaller than 1")
  }
  bias_adjust <- match.arg(bias_adjust)
  if (!is.logical(constrain) || length(constrain) != 1) {
    stop("'constrain' must be a logical scalar")
  }
  if (!is.logical(varN) || length(varN) != 1) {
    stop("'varN' must be a logical scalar")
  }
  which_dj <- match.arg(which_dj)
  # Find the number of (disjoint) blocks
  k_n <- floor(length(data) / b)
  if (k_n < 1) {
    stop("b is too large: it is larger than length(data)")
  }
  #
  # Estimate sigma2_dj based on Section 4 of Berghaus and Bucher (2018)
  # We require the disjoint maxima to do this.  If sliding = TRUE then
  # pass these to spm_sigmahat_dj using the dj_maxima argument
  # Only do this is b_ok = TRUE.
  # Otherwise, just calculate point estimates of theta
  # At this point these estimates have not been bias-adjusted, unless
  # bias_adjust = "N".
  #
  # Find all sets of maxima of disjoint blocks of length b
  all_max <- all_maxima(data, b)
  res <- ests_sigmahat_dj(all_max, b, which_dj, bias_adjust)
  # Sliding maxima
  Fhaty <- ecdf2(all_max$xs, all_max$ys)
  # Avoid over-writing the `disjoint' sample size k_n: it is needed later
  k_n_sl <- length(all_max$ys)
  m <- length(all_max$xs)
  const <- -log(m - b + k_n_sl)
  if (bias_adjust == "N") {
    Fhaty <- (m * Fhaty - b) / (m - b)
  }
  res$theta_sl <- c(-1 / mean(b * log0const(Fhaty, const)),
                    1 / (b * mean(1 - Fhaty)))
  names(res$theta_sl) <- c("N2015", "BB2018")
  #
  # Add the values of the Y-data and the Z-data to the output
  res$data_sl <- cbind(N2015 = -b * log(Fhaty), BB2018 = b * (1 - Fhaty))
  #
  # Estimate the sampling variances of the estimators
  #
  res$sigma2sl <- res$sigma2dj_for_sl - (3 - 4 * log(2)) / res$theta_sl ^ 2
  # res$sigma2sl could contain non-positive values
  # If it does then replace them with NA
  res$sigma2sl[res$sigma2sl <= 0] <- NA
  indexN <- ifelse(varN, 2, 1)
  if (varN) {
    index <- 1:2
  } else {
    index <- c(2, 2)
  }
  res$se_dj <- res$theta_dj ^ 2 * sqrt(res$sigma2dj[index] / k_n)
  res$se_sl <- res$theta_sl ^ 2 * sqrt(res$sigma2sl[index] / k_n)
  #
  # Perform BB2018 bias-adjustment if required
  #
  if (bias_adjust == "BB3") {
    res$bias_dj <- res$theta_dj / k_n + res$theta_dj ^ 3 * res$sigma2dj / k_n
    res$theta_dj <- res$theta_dj - res$bias_dj
    BB3adj_sl <- res$theta_sl / k_n + res$theta_sl ^ 3 * res$sigma2sl / k_n
    BB1adj_sl <- res$theta_sl / k_n
    res$bias_sl <- ifelse(is.na(res$se_sl), BB1adj_sl, BB3adj_sl)
    res$theta_sl <- res$theta_sl - res$bias_sl
    if (is.na(res$se_sl[1])) {
      warning("'bias_adjust' has been changed to ''BB1'' for estimator N2015")
    }
    if (is.na(res$se_sl[2])) {
      warning("'bias_adjust' has been changed to ''BB1'' for estimator BB2018")
    }
  } else if (bias_adjust == "BB1") {
    res$bias_dj <- res$theta_dj / k_n
    res$theta_dj <- res$theta_dj - res$bias_dj
    res$bias_sl <- res$theta_sl / k_n
    res$theta_sl <- res$theta_sl - res$bias_sl
  } else {
    res$bias_dj <- res$bias_sl <- c(N2015 = 0, BB2018 = 0)
  }
  #
  # Save the unconstrained estimates, so that they can be returned
  res$uncon_theta_dj <- res$theta_dj
  res$uncon_theta_sl <- res$theta_sl
  #
  # Constrain to (0, 1] if required
  if (constrain) {
    res$theta_dj <- pmin(res$theta_dj, 1)
    res$theta_sl <- pmin(res$theta_sl, 1)
  }
  #
  res$bias_adjust <- bias_adjust
  res$b <- b
  res$call <- Call
  class(res) <- c("spm", "exdex")
  return(res)
}

#' @keywords internal
#' @rdname exdex-internal
ests_sigmahat_dj <- function(all_max, b, which_dj, bias_adjust){
  # Which of the raw values in x are <= each of the values in y?
  # For each of the block maxima in y calculate the numbers of the raw
  # values in each block that are <= the block maximum
  # k_n is the number of blocks
  k_n <- nrow(all_max$yd)
  # m is the number of raw observations
  m <- nrow(all_max$xd)
  # Value to replace log(0), in the unlikely event that this happens
  const <- -log(m - b + k_n)
  block <- rep(1:k_n, each = b)
  sum_fun <- function(x, y) {
    return(vapply(y, function(y) sum(x <= y), 0))
  }
  # This returns a list with k_n elements.  The ith element of the list
  # (a of length vector k_n) contains the numbers of values in the ith block
  # that are <= each of the block maxima in y
  #
  # Function to calculate Fhaty and UsumN for each set of disjoint block maxima
  UsumN_fn <- function(i) {
    y <- all_max$yd[, i]
    x <- all_max$xd[, i]
    nums_list <- tapply(x, block, sum_fun, y = y)
    # Make this into a matrix
    # Column j contains the numbers of values in the ith block that are <= each
    # of the block maxima in y
    # Row i contains the numbers of values in each block that are <=
    # block maximum i in y
    nums_mat <- do.call(cbind, nums_list)
    # Therefore, the row sums contain the total number of values that are <=
    # each block maximum in y
    # The corresponding proportion is Fhaty in spm(): ecdf of x evaluated at y,
    # in the disjoint (sliding = FALSE) case
    Fhaty <- rowSums(nums_mat) / m
    # For each block, we want an equivalent vector obtained when we delete that
    # block
    fun <- function(i, x) {
      rowSums(x[, -i, drop = FALSE])
    }
    # Column j contains the numbers of values outside of block j that are <= each
    # of the block maxima in y
    # Row i contains the number of values that are outside block 1, ..., k_n
    # and <= block maximum i in y
    # The proportions are Fhat_{-j}(M_{ni}), i, j = 1, ..., k_n
    FhatjMni <- vapply(1:k_n, fun, numeric(k_n), x = nums_mat) / (m - b)
    # Column j enables us to calculate Yhatni(j) and/or Zhatni(j)
    UsumN <- -b * colMeans(log0const(FhatjMni, const))
    Usum <- b * (1 - colMeans(FhatjMni))
    return(list(Nhat = Fhaty, UsumN = UsumN, Usum = Usum))
  }
  fun_value <- list(numeric(k_n), numeric(k_n), numeric(k_n))
  # Use all sets of disjoint maxima to estimate sigmahat2_dj for sliding maxima
  which_vals <- 1:ncol(all_max$yd)
  temp <- vapply(which_vals, UsumN_fn, fun_value)
  Nhat <- do.call(cbind, temp[1, ])
  # BB2018
  Zhat <- b * (1 - Nhat)
  That <- colMeans(Zhat)
  Usum <- do.call(cbind, temp[3, ])
  Usum <-  t(k_n * That - (k_n - 1) * t(Usum))
  Bhat <- t(t(Zhat + Usum) - 2 * That)
  # N2015
  ZhatN <- -b * log(Nhat)
  ThatN <- colMeans(ZhatN)
  UsumN <- do.call(cbind, temp[2, ])
  UsumN <-  t(k_n * ThatN - (k_n - 1) * t(UsumN))
  # Bhat was mean-centred, but BhatN isn't (quite)
  BhatN <- t(t(ZhatN + UsumN) - 2 * ThatN)
  BhatN <- t(t(BhatN) - colMeans(BhatN))
  # Estimate sigma2_dj.
  # First calculate an estimate for each set of disjoint block maxima
  sigmahat2_dj <- apply(Bhat, 2, function(x) sum(x ^ 2) / length(x))
  sigmahat2_djN <- apply(BhatN, 2, function(x) sum(x ^ 2) / length(x))
  # Then, for sliding maxima, take the mean value over all sets of maxima
  sigmahat2_dj_for_sl <- sum(sigmahat2_dj) / length(sigmahat2_dj)
  sigmahat2_dj_for_slN <- sum(sigmahat2_djN) / length(sigmahat2_djN)
  sigma2dj_for_sl <- c(sigmahat2_dj_for_slN, sigmahat2_dj_for_sl)
  # For disjoint maxima pick either the first or last value, based on which_dj
  which_dj <- switch(which_dj, first = 1, last = length(sigmahat2_dj))
  sigma2dj <- c(sigmahat2_djN[which_dj], sigmahat2_dj[which_dj])
  #
  # Point estimates: disjoint maxima. Component i of ThatN (N) and That (BB)
  # contains (the reciprocal of) point estimates of theta based on set of
  # disjoint maxima i.  Use the mean of these estimates as an overall estimate.
  # Perform the Northrop (2015) `bias-adjustment' of Fhaty (Nhat here) if
  # requested and recalculate That and ThatN
  #
  if (bias_adjust == "N") {
    Nhat <- (m * Nhat - b) / (m - b)
    That <- colMeans(b * (1 - Nhat))
    ThatN <- colMeans(-b * log0const(Nhat, const))
  }
  theta_dj <- 1 / c(ThatN[which_dj], That[which_dj])
  names(theta_dj) <- names(sigma2dj) <- names(sigma2dj_for_sl) <-
    c("N2015", "BB2018")
  return(list(sigma2dj = sigma2dj, sigma2dj_for_sl = sigma2dj_for_sl,
              theta_dj = theta_dj,
              data_dj = cbind(N2015 = -b * log(Nhat[, which_dj]),
                              BB2018 = b * (1 - Nhat[, which_dj]))))
}

# ================= Functions used in spm() and spm_R_quick() =============== #

# log(x), but return a constant const for an x = 0

#' @keywords internal
#' @rdname exdex-internal
log0const_slow <- function(x, const) {
  ifelse(x == 0, const, log(x))
}

#' @keywords internal
#' @rdname exdex-internal
log0const <- function(x, const) {
  return(log(x + !x) + const * !x)
}

# =================== Empirical c.d.f. of x, evaluated at y ================= #

#' @keywords internal
#' @rdname exdex-internal
ecdf3 <- function(x, y) {
  return(vapply(y, function(y) mean(x <= y), 0))
}

#' @keywords internal
#' @rdname exdex-internal
ecdf2 <- function(x, y) {
  return(vapply(y, function(y) sum(x <= y) / length(x), 0))
}

#' @keywords internal
#' @rdname exdex-internal
ecdf1 <- function(x, y, lenx) {
  return(vapply(y, function(y) sum(x <= y) / lenx, 0))
}

# ======== Functions to calculate block maxima, used only in testing ======== #

#' @keywords internal
#' @rdname exdex-internal
sliding_maxima <- function(x, b = 1){
  y <- as.numeric(zoo::rollapply(data = zoo::zoo(x), width = b, FUN = max,
                                 na.rm = TRUE))
  return(list(y = y, x = x))
}

#' @keywords internal
#' @rdname exdex-internal
disjoint_maxima <- function(x, b = 1, which_dj = c("first", "last")){
  which_dj <- match.arg(which_dj)
  if (which_dj == "last") {
    x <- rev(x)
  }
  n <- length(x)
  # number of maxima of blocks of length b
  n_max <- floor(n / b)
  # take only the first 1 to n_max*b observations
  x <- x[1:(n_max * b)]
  # block indicator: 1, ..., 1, ..., n_max, ..., n_max
  ind <- rep(1:n_max, each = b)
  # calculate block maxima
  y <- as.numeric(tapply(x, ind, max, na.rm = TRUE))
  if (which_dj == "last") {
    x <- rev(x)
    y <- rev(y)
  }
  return(list(y = y, x = x))
}

#' @keywords internal
#' @rdname exdex-internal
all_disjoint_maxima <- function(x, b = 1,
                                which_dj = c("all", "first", "last")){
  which_dj <- match.arg(which_dj)
  n <- length(x)
  # The number of maxima of blocks of length b
  n_max <- floor(n / b)
  # Set the possible first indices
  first_value <- switch(which_dj,
                        all = 1:(n - n_max * b + 1),
                        first = 1,
                        last = n - n_max * b + 1)
  # block indicator: 1, ..., 1, ..., n_max, ..., n_max
  ind <- rep(1:n_max, each = b)
  get_maxima <- function(first) {
    last <- first + n_max * b - 1
    # take only the first_value to n_max * b + first_value - 1 observations
    xx <- x[first:last]
    # calculate block maxima
    y <- as.numeric(tapply(xx, ind, max, na.rm = TRUE))
    return(c(y, xx))
  }
  temp <- vapply(first_value, FUN = get_maxima, numeric(n_max * (b + 1)))
  y_mat <- temp[1:n_max, , drop = FALSE]
  x_mat <- temp[-(1:n_max), , drop = FALSE]
  return(list(y_mat = y_mat, x_mat = x_mat))
}

#' @keywords internal
#' @rdname exdex-internal
all_maxima <- function(x, b = 1, which_dj = c("all", "first", "last")){
  which_dj <- match.arg(which_dj)
  # First calculate the sliding block maxima.  All the disjoint maxima that
  # we need are contained in s_max, and we need the sliding maxima anyway
  s_max <- sliding_maxima(x = x, b = b)
  # The number of maxima of blocks of length b
  n <- length(x)
  n_max <- floor(n / b)
  # Set the possible first indices
  first_value <- switch(which_dj,
                        all = 1:(n - n_max * b + 1),
                        first = 1,
                        last = n - n_max * b + 1)
  # A function to return block maxima and contributing values starting from
  # the first value first
  get_maxima <- function(first) {
    s_ind <- seq.int(from = first, by = b, length.out = n_max)
    return(c(s_max$y[s_ind], x[first:(first + n_max * b - 1)]))
  }
  temp <- vapply(first_value, FUN = get_maxima, numeric(n_max * (b + 1)))
  yd <- temp[1:n_max, , drop = FALSE]
  xd <- temp[-(1:n_max), , drop = FALSE]
  return(list(ys = s_max$y, xs = s_max$x, yd = yd, xd = xd))
}

# ============== Functions used by kgaps() and confint.kgaps() ============== #

# =============================== kgaps_loglik ================================
#' @keywords internal
#' @rdname exdex-internal
kgaps_loglik <- function(theta, N0, N1, sum_qs, n_kgaps){
  if (theta < 0 || theta > 1) {
    return(-Inf)
  }
  loglik <- 0
  if (N1 > 0) {
    loglik <- loglik + 2 * N1 * log(theta) - sum_qs * theta
  }
  if (N0 > 0) {
    loglik <- loglik + N0 * log(1 - theta)
  }
  return(loglik)
}

# ============================== kgaps_conf_int ===============================
#' @keywords internal
#' @rdname exdex-internal
kgaps_conf_int <- function(theta_mle, ss, conf = 95) {
  cutoff <- stats::qchisq(conf / 100, df = 1)
  theta_list <- c(list(theta = theta_mle), ss)
  max_loglik <- do.call(kgaps_loglik, theta_list)
  ob_fn <- function(theta) {
    theta_list$theta <- theta
    loglik <- do.call(kgaps_loglik, theta_list)
    return(2 * (max_loglik - loglik) - cutoff)
  }
  ci_low <- 0
  ci_up <- 1
  if (ss$N1 > 0) {
    ci_low <- stats::uniroot(ob_fn, c(0, theta_mle))$root
  }
  if (ss$N0 > 0) {
    ci_up <- stats::uniroot(ob_fn, c(theta_mle, 1))$root
  }
  return(c(ci_low, ci_up))
}

# ============================== kgaps_quad_solve =============================
#' @keywords internal
#' @rdname exdex-internal
kgaps_quad_solve <- function(N0, N1, sum_qs) {
  aa <- sum_qs
  bb <- -(N0 + 2 * N1 + sum_qs)
  cc <- 2 * N1
  qq <- -(bb - sqrt(bb ^ 2 - 4 * aa * cc)) / 2
  theta_mle <- cc / qq
  return(theta_mle)
}

# ========================= Function used by iwls() ========================= #

# ================================== iwls_fun =================================
#' @keywords internal
#' @rdname exdex-internal
iwls_fun <- function(n_wls, N, S_1_sort, exp_qs, ws, nx) {
  #
  # This function implements the algorithm on page 46 of Suveges (2007).
  # [In step (1) there is a typo in the paper: in x_i the N_C+1 should be N.]
  #
  # Args:
  # n_wls    : A numeric scalar.  The number of the largest 1-gaps to include
  #            in the weighted least squares estimation.
  # N        : A numeric scalar.  The number of threshold excesses.
  # S_1_sort : A numeric N-vector.  Sorted (largest to smallest) scaled 1-gaps.
  #            The scaling multiplies the raw 1-gaps by the sample proportion
  #            of values that exceed u.
  # exp_qs   : A numeric N-vector.  Standard exponential quantiles (order
  #            statistics) for a sample of size N.
  # ws       : A numeric N-vector.  Weights for the least squares fit.
  # nx       : A numeric scalar.  The number of raw observations.
  #
  # Returns: A list with components
  #    theta : A numeric scalar.  The new estimate of theta.
  #    n_wls : A numeric scalar.  The new value of n_wls.
  #
  # Extract the values corresponding to the largest n_wls 1-gaps
  # Extract the largest n_wls scaled 1-gaps (ordered largest to smallest)
  chi_i <- S_1_sort[1:n_wls]
  # Standard exponential quantiles, based on N 1-gaps (largest to smallest)
  x_i <- exp_qs[1:n_wls]
  # Extract the weights for the values in chi_i
  ws <- ws[1:n_wls]
  # Weighted least squares for (chi_i, x_i)
  temp <- stats::lm(chi_i ~ x_i, weights = ws)
  ab <- temp$coefficients
  # Estimate theta
  theta <- min(exp(ab[1] / ab[2]), 1)
  # Update n_wls
  n_wls <- floor(theta * (N - 1))
  return(list(theta = theta, n_wls = n_wls))
}
