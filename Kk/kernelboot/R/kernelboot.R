
#' Smoothed bootstrap
#'
#' Smoothed bootstrap is an extension of standard bootstrap using kernel densities.
#'
#' @param data       vector, matrix, or data.frame. For non-numeric values standard bootstrap
#'                   is applied (see below).
#' @param statistic  a function that is applied to the \code{data}. The first argument of
#'                   the function will always be the original data.
#' @param R          the number of bootstrap replicates.
#' @param bw         the smoothing bandwidth to be used (see \code{\link[stats]{density}}).
#'                   The kernels are scaled such that this is the standard deviation,
#'                   or covariance matrix of the smoothing kernel. By default
#'                   \code{\link[stats]{bw.nrd0}} is used for univariate data,
#'                   and \code{\link{bw.silv}} is used for multivariate data. When using
#'                   \code{kernel = "multivariate"} this parameter should be a
#'                   \emph{covariance matrix} of the smoothing kernel.
#' @param kernel     a character string giving the smoothing kernel to be used.
#'                   This must partially match one of "multivariate", "gaussian",
#'                   "rectangular", "triangular", "epanechnikov", "biweight", "cosine",
#'                   "optcosine", or "none" with default "multivariate", and may be abbreviated.
#'                   Using \code{kernel = "multivariate"} forces multivariate Gaussian kernel
#'                   (or univariate Gaussian for univariate data). Using \code{kernel = "none"}
#'                   forces using standard bootstrap (no kernel smoothing).
#' @param adjust     scalar; the bandwidth used is actually \code{adjust*bw}. This makes
#'                   it easy to specify values like 'half the default' bandwidth.
#' @param weights    vector of importance weights. It should have as many elements
#'                   as there are observations in \code{data}. It defaults to uniform
#'                   weights.
#' @param shrinked   logical; if \code{TRUE} random generation algorithm preserves
#'                   means and variances of the variables. This parameter is ignored for
#'                   "multivariate" kernel.
#' @param ignore     vector of names of columns to be ignored during the smoothing phase of
#'                   bootstrap procedure (their values are not altered using random noise).
#' @param parallel   if \code{TRUE}, parallel computing is used (see \code{\link[future.apply]{future_lapply}}).
#'                   \emph{Warning:} using parallel computing does not necessary have to
#'                   lead to improved performance.
#' @param workers    the number of workers used for parallel computing (see \code{\link[future]{multiprocess}}).
#'
#'
#' @details
#'
#' \emph{Smoothed bootstrap} is an extension of standard bootstrap procedure, where instead
#' of drawing samples with replacement from the empirical distribution, they are drawn
#' from kernel density estimate of the distribution.
#'
#' For smoothed bootstrap, points (in univariate case), or rows (in multivariate case), are drawn with
#' replacement, to obtain samples of size \eqn{n} from the initial dataset of size \eqn{n}, as with
#' standard bootstrap. Next, random noise from kernel density \eqn{K} is added to each of the drawn
#' values. The procedure is repeated \eqn{R} times and \code{statistic} is evaluated on each of the
#' samples.
#'
#' The noise is added \emph{only} to the numeric columns, while non-numeric columns (e.g.
#' \code{character}, \code{factor}, \code{logical}) are not altered. What follows, to the
#' non-numeric columns and columns listed in \code{ignore} parameter standard bootstrap procedure
#' is applied.
#'
#'
#' \strong{Univariate kernel densities}
#'
#' Univariate kernel density estimator is defined as
#'
#' \deqn{
#' \hat{f_h}(x) = \sum_{i=1}^n w_i \, K_h(x-y_i)
#' }{
#' f(x) = sum[i](w[i] * Kh(x-y[i]))
#' }
#'
#' where \eqn{w} is a vector of weights such that all \eqn{w_i \ge 0}{w[i] \ge 0}
#' and \eqn{\sum_i w_i = 1}{sum(w) = 1} (by default uniform \eqn{1/n} weights are used),
#' \eqn{K_h = K(x/h)/h}{Kh = K(x/h)/h} is kernel \eqn{K} parametrized by bandwidth \eqn{h}
#' and \eqn{y} is a vector of data points used for estimating the kernel density.
#'
#' To draw samples from univariate kernel density, the following procedure can be applied (Silverman, 1986):
#'
#' \emph{Step 1} Sample \eqn{i} uniformly with replacement from \eqn{1,\dots,n}.
#'
#' \emph{Step 2} Generate \eqn{\varepsilon}{\epsilon} to have probability density \eqn{K}.
#'
#' \emph{Step 3} Set \eqn{x = y_i + h\varepsilon}{x = y[i] + h\epsilon}.
#'
#' If samples are required to have the same variance as \code{data}
#' (i.e. \code{shrinked = TRUE}), then \emph{Step 3} is modified
#' as following:
#'
#' \emph{Step 3'} \eqn{
#' x = \bar y + (y_i - \bar y + h\varepsilon)/(1 + h^2 \sigma^2_K/\sigma^2_Y)^{1/2}
#' }{
#' x = m + (y[i] - m + h\epsilon)/sqrt(1 + h^2 var(K)/var(y))
#' }
#'
#' where \eqn{\sigma_K^2}{var(K)} is variance of the kernel (fixed to 1 for kernels used in this package).
#'
#' When shrinkage described in \emph{Step 3'} is applied, the smoothed bootstrap density function changes
#' it's form to
#'
#' \deqn{
#' \hat{f}_{h,b}(x) = (1 + r) \; \hat{f_h}(x + r(x - \bar{y}))
#' }{
#' fb(x) = (1+r) f(x + r (x-mean(y)))
#' }
#'
#' where \eqn{r = \left(1 + h^2 \sigma_K^2 / \sigma_y^2 \right)^{1/2}-1}{r = sqrt(1 + h^2 var(K)/var(y)) - 1}.
#'
#' This package offers the following univariate kernels:
#'
#' \tabular{ll}{
#' \emph{Gaussian}     \tab \eqn{\frac{1}{\sqrt{2\pi}} e^{-{u^2}/2}}{1/sqrt(2\pi) exp(-(u^2)/2)} \cr
#' \emph{Rectangular}  \tab \eqn{\frac{1}{2} \ \mathbf{1}_{(|u|\leq1)}}{1/2} \cr
#' \emph{Triangular}   \tab \eqn{(1-|u|) \ \mathbf{1}_{(|u|\leq1)}}{1 - |u|} \cr
#' \emph{Epanchenikov} \tab \eqn{\frac{3}{4}(1-u^2) \ \mathbf{1}_{(|u|\leq1)}}{3/4 (1 - u^2)} \cr
#' \emph{Biweight}     \tab \eqn{\frac{15}{16}(1-u^2)^2 \ \mathbf{1}_{(|u|\leq1)}}{15/16 (1 - u^2)^2} \cr
#' \emph{Cosine}       \tab \eqn{\frac{1}{2} \left(1 + \cos(\pi u)\right) \ \mathbf{1}_{(|u|\leq1)}}{1/2 (1 + cos(\pi u))} \cr
#' \emph{Optcosine}    \tab \eqn{\frac{\pi}{4}\cos\left(\frac{\pi}{2}u\right) \ \mathbf{1}_{(|u|\leq1)}}{\pi/4 cos(\pi/2 u)}
#' }
#'
#' All the kernels are re-scalled so that their standard deviations are equal to 1,
#' so that bandwidth parameter controls their standard deviations.
#'
#' Random generation from Epanchenikov kernel is done using algorithm
#' described by Devroye (1986). For optcosine kernel inverse transform
#' sampling is used. For biweight kernel random values are drawn from
#' \eqn{\mathrm{Beta}(3, 3)}{Beta(3, 3)} distribution and
#' \eqn{\mathrm{Beta}(3.3575, 3.3575)}{Beta(3.3575, 3.3575)}
#' distribution serves as a close approximation of cosine kernel.
#' Random generation for triangular kernel is done by taking difference
#' of two i.i.d. uniform random variates. To sample from rectangular
#' and Gaussian kernels standard random generation algorithms are used
#' (see \code{\link[stats]{runif}} and \code{\link[stats]{rnorm}}).
#'
#'
#' \strong{Product kernel densities}
#'
#' Univariate kernels may easily be extended to multiple dimensions by
#' using product kernel
#'
#' \deqn{
#' \hat{f_H}(\mathbf{x}) = \sum_{i=1}^n w_i \prod_{j=1}^m
#' K_{h_j}(x_i - y_{ij})
#' }{
#' f(x) = sum[i](w[i] * prod[j]( Kh[j](x[i]-y[i,j) ))
#' }
#'
#' where \eqn{w} is a vector of weights such that all \eqn{w_i \ge 0}{w[i] \ge 0}
#' and \eqn{\sum_i w_i = 1}{sum(w) = 1} (by default uniform \eqn{1/n} weights are used),
#' and \eqn{K_{h_j}}{Kh[j]} are univariate kernels \eqn{K} parametrized by bandwidth
#' \eqn{h_j}{h[j]}, where \eqn{\boldsymbol{y}}{y} is a matrix of data points used for
#' estimating the kernel density.
#'
#' Random generation from product kernel is done by drawing with replacement
#' rows of \eqn{y}, and then adding to the sampled values random noise from
#' univariate kernels \eqn{K}, parametrized by corresponding bandwidth parameters
#' \eqn{h_j}{h[j]}.
#'
#'
#' \strong{Multivariate kernel densities}
#'
#' Multivariate kernel density estimator may also be defined in terms of multivariate kernels
#' \eqn{K_H}{KH} (e.g. multivariate normal distribution, as in this package)
#'
#' \deqn{
#' \hat{f_H}(\mathbf{x}) = \sum_{i=1}^n w_i \, K_H( \mathbf{x}-\boldsymbol{y}_i)
#' }{
#' f(x) = sum[i](w[i] * KH(x-y[i]))
#' }
#'
#' where \eqn{w} is a vector of weights such that all \eqn{w_i \ge 0}{w[i] \ge 0}
#' and \eqn{\sum_i w_i = 1}{sum(w) = 1} (by default uniform \eqn{1/n} weights are used),
#' \eqn{K_H}{KH} is kernel \eqn{K} parametrized by bandwidth matrix \eqn{H} and
#' \eqn{\boldsymbol{y}}{y} is a matrix of data points used for estimating the kernel density.
#'
#' \emph{Notice:} When using multivariate normal (Gaussian) distribution as a kernel \eqn{K}, the
#' bandwidth parameter \eqn{H} is a \emph{covariance matrix} as compared to standard deviations
#' used in univariate and product kernels.
#'
#' Random generation from multivariate kernel is done by drawing with replacement
#' rows of \eqn{y}, and then adding to the sampled values random noise from
#' multivariate normal distribution centered at the data points and parametrized
#' by corresponding bandwidth matrix \eqn{H}. For further details see \code{\link{rmvg}}.
#'
#'
#' @references
#' Silverman, B. W. (1986). Density estimation for statistics and data analysis.
#' Chapman and Hall/CRC.
#'
#' @references
#' Scott, D. W. (1992). Multivariate density estimation: theory, practice,
#' and visualization. John Wiley & Sons.
#'
#' @references
#' Efron, B. (1981). Nonparametric estimates of standard error: the jackknife,
#' the bootstrap and other methods. Biometrika, 589-599.
#'
#' @references
#' Hall, P., DiCiccio, T.J. and Romano, J.P. (1989). On smoothing and the bootstrap.
#' The Annals of Statistics, 692-704.
#'
#' @references
#' Silverman, B.W. and Young, G.A. (1987). The bootstrap: To smooth or not to smooth?
#' Biometrika, 469-479.
#'
#' @references
#' Scott, D.W. (1992). Multivariate density estimation: theory, practice,
#' and visualization. John Wiley & Sons.
#'
#' @references
#' Wang, S. (1995). Optimizing the smoothed bootstrap. Annals of the Institute of
#' Statistical Mathematics, 47(1), 65-80.
#'
#' @references
#' Young, G.A. (1990). Alternative smoothed bootstraps. Journal of the Royal
#' Statistical Society. Series B (Methodological), 477-484.
#'
#' @references
#' De Angelis, D. and Young, G.A. (1992). Smoothing the bootstrap.
#' International Statistical Review/Revue Internationale de Statistique, 45-56.
#'
#' @references
#' Polansky, A.M. and Schucany, W. (1997). Kernel smoothing to improve bootstrap
#' confidence intervals. Journal of the Royal Statistical Society: Series B
#' (Statistical Methodology), 59(4), 821-838.
#'
#' @references
#' Devroye, L. (1986). Non-uniform random variate generation. New York: Springer-Verlag.
#'
#' @references
#' Parzen, E. (1962). On estimation of a probability density function and mode.
#' The annals of mathematical statistics, 33(3), 1065-1076.
#'
#' @references
#' Silverman, B.W. and Young, G.A. (1987). The bootstrap: To smooth or not to smooth?
#' Biometrika, 469-479.
#'
#' @references
#' Jones, M.C. (1991). On correcting for variance inflation in kernel density estimation.
#' Computational Statistics & Data Analysis, 11, 3-15.
#'
#'
#' @seealso \code{\link{bw.silv}}, \code{\link[stats]{density}},
#'          \code{\link[stats]{bandwidth}}, \code{\link{kernelboot-class}}
#'
#'
#' @examples
#'
#' set.seed(1)
#'
#' # smooth bootstrap of parameters of linear regression
#'
#' b1 <- kernelboot(mtcars, function(data) coef(lm(mpg ~ drat + wt, data = data)) , R = 250)
#' b1
#' summary(b1)
#'
#' b2 <- kernelboot(mtcars, function(data) coef(lm(mpg ~ drat + wt, data = data)) , R = 250,
#'                  kernel = "epanechnikov")
#' b2
#' summary(b2)
#'
#' # smooth bootstrap of parameters of linear regression
#' # smoothing phase is not applied to "am" and "cyl" variables
#'
#' b3 <- kernelboot(mtcars, function(data) coef(lm(mpg ~ drat + wt + am + cyl, data = data)) , R = 250,
#'                  ignore = c("am", "cyl"))
#' b3
#' summary(b3)
#'
#' # standard bootstrap (without kernel smoothing)
#'
#' b4 <- kernelboot(mtcars, function(data) coef(lm(mpg ~ drat + wt + am + cyl, data = data)) , R = 250,
#'                  ignore = colnames(mtcars))
#' b4
#' summary(b4)
#'
#' # smooth bootstrap for median of univariate data
#'
#' b5 <- kernelboot(mtcars$mpg, function(data) median(data) , R = 250)
#' b5
#' summary(b5)
#'
#'
#' @importFrom stats rnorm bw.SJ bw.bcv bw.nrd bw.nrd0 bw.ucv
#' @importFrom future plan multiprocess
#' @importFrom future.apply future_lapply
#'
#' @export

kernelboot <- function(data, statistic, R = 500L, bw = "default",
                       kernel = c("multivariate", "gaussian", "epanechnikov",
                                  "rectangular", "triangular", "biweight",
                                  "cosine", "optcosine", "none"),
                       weights = NULL, adjust = 1,
                       shrinked = TRUE, ignore = NULL,
                       parallel = FALSE, workers = 1L) {

  call <- match.call()
  kernel <- match.arg(kernel)
  seed <- get0(".Random.seed", envir = .GlobalEnv, ifnotfound = NULL)
  n <- NROW(data)
  m <- NCOL(data)
  vars <- NULL

  if (!(is.simple.vector(data) || is.data.frame(data) || is.matrix(data)))
    stop("unsupported data type")

  if (is.data.frame(data) || is.matrix(data)) {
    num_cols <- numericColumns(data)
    ignr_cols <- colnames(data) %in% ignore
    incl_cols <- num_cols & !ignr_cols
  }

  if (is.character(bw)) {
    method <- bw
    if (is.data.frame(data) || is.matrix(data)) {
      bw <- matrix(0, m, m)
      if (!is.null(colnames(data)))
        rownames(bw) <- colnames(bw) <- colnames(data)
      bw[incl_cols, incl_cols] <- calculate_bandwidth(data[, incl_cols], method, kernel == 'multivariate')
    } else {
      bw <- calculate_bandwidth(data, method, FALSE)
    }
  }

  if (!is.simple.vector(adjust) || length(adjust) > 1L)
    stop("adjust is not a scalar")

  bw <- bw * adjust

  # check for non-numeric, NAs, NaNs, infinite values
  if (!all(is.finite(bw)))
    stop("inappropriate values of bw")

  if (!is.null(weights)) {
    if (!all(is.finite(weights)))
      stop("inappropriate values of weights")
  }

  # equally weighted
  if (length(weights) == 1L)
    weights <- NULL

  # try evaluating statistic() on the original data

  tryCatch(
    orig.stat <- statistic(data),
    error = function(e) {
      message("applying the statistic on the original data resulted in an error")
      stop(e)
    }
  )

  # looping functions - paralell or lapply

  if (parallel && workers > 1L) {

    # using future for parallel computing
    repeatFun <- function(n, FUN, workers) {
      plan(multiprocess, workers = workers)
      future_lapply(1:n, FUN, future.seed = TRUE)
    }

  } else {

    parallel <- FALSE
    repeatFun <- function(n, FUN, workers) lapply(1:n, FUN)

  }

  if (is.data.frame(data) || is.matrix(data)) {

    if (!is.null(colnames(data))) {
      vars <- list(
        smoothed = colnames(data)[incl_cols],
        ignored  = colnames(data)[!incl_cols]
      )
    } else {
      vars <- list(
        smoothed = which(incl_cols),
        ignored  = which(!incl_cols)
      )
    }

    if (kernel == "none" || !any(incl_cols) || is.allzeros(bw)) {

      # standard bootstrap

      kernel.type <- "multivariate"
      kernel <- "none"

      res <- repeatFun(R, function(i) {

        idx <- sample.int(n, n, replace = TRUE, prob = weights)
        boot.data <- data[idx, , drop = FALSE]
        statistic(boot.data)

      }, workers = workers)

    } else {

      # smoothed bootstrap

      data_mtx <- as.matrix(data[, incl_cols])

      if (is.null(weights))
        weights <- rep(1/n, n)

      if (kernel != "multivariate") {

        # product kernel

        kernel.type <- "product"

        if (is.simple.vector(bw)) {
          if (length(bw) == 1L)
            bw <- rep(bw, m)
        } else {
          if (!is.square(bw))
            stop("bw is not a square matrix")
          bw <- diag(bw)
        }

        bw <- bw[incl_cols]

        res <- repeatFun(R, function(i) {

          samp <- cpp_rmvk(n, data_mtx, bw, weights, kernel)
          idx <- attr(samp, "boot_index")
          attr(samp, "boot_index") <- NULL
          boot.data <- data[idx, , drop = FALSE]
          boot.data[, incl_cols] <- samp
          statistic(boot.data)

        }, workers = workers)

      } else {

        # MVN kernel

        kernel.type <- "multivariate"

        # is this check really needed?
        # if (qr(data_mtx)$rank < min(dim(data_mtx)))
        #   warning("data matrix is rank deficient")

        if (is.simple.vector(bw)) {
          if (length(bw) == 1L)
            bw <- diag(bw, nrow = ncol(data))
          else
            bw <- diag(bw)
        }

        if (ncol(bw) != m)
          stop("dimensions of data and bw do not match")

        if (!is.square(bw))
          stop("bw is not a square matrix")

        bw <- as.matrix(bw)
        bw <- bw[incl_cols, incl_cols]
        bw_chol <- chol(bw)
        mm <- sum(incl_cols)

        res <- repeatFun(R, function(i) {

          idx <- sample.int(n, n, replace = TRUE, prob = weights)
          boot.data <- data[idx, , drop = FALSE]
          samp <- matrix(rnorm(n*mm), n, mm) %*% bw_chol
          boot.data[, incl_cols] <- boot.data[, incl_cols] + samp
          statistic(boot.data)

        }, workers = workers)

      }
    }

  } else if (is.simple.vector(data)) {

    # data is a vector

    kernel.type <- "univariate"

    if (kernel == "none" || !is.numeric(data) || is.allzeros(bw)) {

      # standard bootstrap

      kernel <- "none"

      res <- repeatFun(R, function(i) {

        idx <- sample.int(n, n, replace = TRUE, prob = weights)
        boot.data <- data[idx]
        statistic(boot.data)

      }, workers = workers)

    } else {

      # smoothed bootstrap

      if (kernel == "multivariate") {
        kernel <- "gaussian"
        message("data is univariate, switching to 'gaussian' kernel")
      }

      if (!is.simple.vector(bw))
        stop("bw is not a scalar")
      if (length(bw) > 1L) {
        bw <- bw[1L]
        message("bw has length > 1 and only the first element will be used")
      }

      if (is.null(weights))
        weights <- rep(1/n, n)

      res <- repeatFun(R, function(i) {

        samp <- cpp_ruvk(n, data, bw, weights, kernel, shrinked)
        boot.data <- as.vector(samp)
        statistic(boot.data)

      }, workers = workers)

    }

  } else {

    stop("unsupported data type")

  }

  # simplify the results to data.frame
  samples <- do.call(rbind, res)

  structure(list(
    orig.stat     = orig.stat,
    boot.samples  = samples,
    call          = call,
    statistic     = statistic,
    orig.data     = data,
    variables     = vars,
    type          = kernel.type,
    param = list(
      R           = R,
      bw          = bw,
      adjust      = adjust,
      weights     = weights,
      kernel      = kernel,
      shrinked    = shrinked,
      parallel    = parallel,
      random.seed = seed
    )
  ), class = "kernelboot")

}


calculate_bandwidth <- function(data, method, multivariate) {
  method <- tolower(method)
  if (method == "default") {
    if (is.simple.vector(data)) {
      bw <- bw.nrd0(data)
    } else {
      bw <- bw.silv(data)
      if (!multivariate)
        bw <- sqrt(diag(bw))
    }
  } else {
    bw <- switch(method, nrd0 = bw.nrd0(data), nrd = bw.nrd(data),
                 ucv = bw.ucv(data), bcv = bw.bcv(data), sj = ,
                 `sj-ste` = bw.SJ(data, method = "ste"),
                 `sj-dpi` = bw.SJ(data, method = "dpi"),
                 silv = bw.silv(data), scott = bw.scott(data),
                 stop("unknown bandwidth rule"))
  }
  return(bw)
}
