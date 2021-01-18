

#' 'kernelboot' class object
#'
#' @details
#'
#' Object of class \code{"kernelboot"}, is a list with components including
#'
#' \tabular{ll}{
#' \code{orig.stat}          \tab  estimates from \code{statistic} on the original data, \cr
#' \code{boot.samples}       \tab  samples drawn, \cr
#' \code{call}               \tab  function call, \cr
#' \code{statistic}          \tab  actual \code{statistic} function that was used, \cr
#' \code{orig.data}          \tab  original data used for bootstrapping, \cr
#' \code{variables}          \tab  used variables: it is \code{NULL} for univariate data and
#'                                 for multivariate data it contains two lists of \code{smoothed}
#'                                 and \code{ignored} variables (names or column indexes) during
#'                                 the smoothing phase. \cr
#' \code{type}               \tab  type of kernel density that was used ("univariate", "product",
#'                                 "multivariate"), \cr
#' \code{param}              \tab  list of parameters that were used.
#' }
#'
#' \code{param} section contains:
#'
#' \tabular{ll}{
#' \code{R}                  \tab  number of bootstrap iterations, \cr
#' \code{bw}                 \tab  the bandwidth that was used, \cr
#' \code{weights}            \tab  vector of the weights that were applied, \cr
#' \code{kernel}             \tab  name of the kernel that was used ("multivariate",
#'                                 "gaussian", "epanechnikov", "rectangular",
#'                                 "triangular", "biweight", "cosine", "optcosine",
#'                                 "none"), \cr
#' \code{shrinked}           \tab  value of the \code{shrinked} parameter, \cr
#' \code{parallel}           \tab  indicates if parallel computation was used, \cr
#' \code{random.seed}        \tab  random seed used to initialize the random number
#'                                 generator (see \code{\link[base]{.Random.seed}}).
#' }
#'
#' @seealso \code{\link{kernelboot}}
#'
#' @name kernelboot-class
NULL


#' Summarize the result of kernelboot
#'
#' @param object     \code{kernelboot} class object.
#' @param probs      quantiles returned by \code{summary} (see \code{\link{quantile}}).
#' @param \dots      further arguments passed to or from other methods.
#' @param na.rm      a logical value indicating whether \code{NA} values should be
#'                   stripped before the computation proceeds.
#'
#' @importFrom stats sd quantile na.omit
#' @export

summary.kernelboot <- function(object, probs = c(0.025, 0.5, 0.975), ..., na.rm = FALSE) {
  samp <- object$boot.samples
  res <- lapply(1:ncol(samp), function(i) {
    x <- samp[, i]
    if (is.numeric(x)) {
      if (na.rm)
        x <- na.omit(x)
      c(mean = mean(x), sd = sd(x), quantile(x, probs = probs))
    } else {
      warning("skipping non-numeric variable")
      NA
    }
  })
  names(res) <- colnames(samp)
  do.call(rbind, res)
}


#' @export

print.kernelboot <- function(x, ...) {

  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), sep = "")
  cat("\n\n")

  cat(x$param$R, " samples were generated", sep = "")
  if (x$param$kernel == "none") {
    cat(" using standard bootstrap.\n", sep = "")
  } else {
    kernel <- if (x$type == "multivariate") "multivariate gaussian" else x$param$kernel
    cat(" from ", x$type, " kernel density with ", kernel, " kernel.\n", sep = "")
  }

  invisible(x)
}
