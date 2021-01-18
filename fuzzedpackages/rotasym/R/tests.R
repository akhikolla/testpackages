

#' @title Tests of rotational symmetry for hyperspherical data
#'
#' @description Tests for assessing the rotational symmetry of a unit-norm
#' random vector \eqn{\mathbf{X}}{X} in
#' \eqn{S^{p-1}:=\{\mathbf{x}\in R^p:||\mathbf{x}||=1\}}{
#' S^{p-1} := \{x \in R^p : ||x|| = 1\}}, \eqn{p \ge 2}, about a location
#' \eqn{\boldsymbol{\theta}\in S^{p-1}}{\theta \in S^{p-1}}, from a
#' hyperspherical sample \eqn{\mathbf{X}_1,\ldots,\mathbf{X}_n\in S^{p-1}}{
#' X_1, \ldots, X_n \in S^{p-1}}.
#'
#' The vector \eqn{\mathbf{X}}{X} is said to be rotational symmetric about
#' \eqn{\boldsymbol{\theta}}{\theta} if the distributions of
#' \eqn{\mathbf{OX}}{OX} and \eqn{\mathbf{X}}{X} coincide, where
#' \eqn{\mathbf{O}}{O} is any \eqn{p\times p}{p x p} rotation matrix
#' that fixes \eqn{\boldsymbol{\theta}}{\theta}, \emph{i.e.},
#' \eqn{\mathbf{O}\boldsymbol{\theta}=\boldsymbol{\theta}}{O\theta = \theta}.
#'
#' @param data hyperspherical data, a matrix of size \code{c(n, p)} with unit
#' norm rows. Normalized internally if any row does not have unit norm
#' (with a \code{warning} message). \code{NA}s are ignored.
#' @param theta either a unit norm vector of size \code{p} giving the axis of
#' rotational symmetry (for the specified-\eqn{\boldsymbol{\theta}}{\theta}
#' case) or a function that implements an estimator
#' \eqn{\hat{\boldsymbol{\theta}}}{\hat \theta} of
#' \eqn{\boldsymbol{\theta}}{\theta} (for the
#' unspecified-\eqn{\boldsymbol{\theta}}{\theta} case). The default calls
#' the \code{\link{spherical_mean}} function. See examples.
#' @param type a character string (case insensitive) indicating the type of
#' test to conduct:
#' \itemize{
#' \item \code{"sc"}: "scatter" test based on the statistic
#' \eqn{Q_{\boldsymbol{\theta}}^{\mathrm{sc}}}{Q_\theta^{sc}}. Evaluates if the
#' covariance matrix of the multivariate signs is isotropic.
#' \item \code{"loc"}: "location" test based on the statistic
#' \eqn{Q_{\boldsymbol{\theta}}^{\mathrm{loc}}}{Q_\theta^{loc}}. Evaluates if
#' the expectation of the multivariate signs is zero.
#' \item \code{"loc_vMF"}: adapted "location" test, based on the statistic
#' \eqn{Q_{\mathrm{vMF}}^{\mathrm{loc}}}{Q_{vMF}^{loc}}.
#' \item \code{"hyb"}: "hybrid" test based on the statistics
#' \eqn{Q_{\boldsymbol{\theta}}^{\mathrm{sc}}}{Q_\theta^{sc}} and
#' \eqn{Q_{\boldsymbol{\theta}}^{\mathrm{loc}}}{Q_\theta^{loc}}.
#' \item \code{"hyb_vMF"} (default): adapted "hybrid" test based on the
#' statistics \eqn{Q_{\boldsymbol{\theta}}^{\mathrm{sc}}}{Q_\theta^{sc}} and
#' \eqn{Q_{\mathrm{vMF}}^{\mathrm{loc}}}{Q_{vMF}^{loc}}.
#' }
#' See the details below for further explanations of the tests.
#' @param Fisher if \code{TRUE}, then Fisher's method is employed to aggregate
#' the scatter and location tests in the hybrid test, see details below.
#' Otherwise, the hybrid statistic is the sum of the scatter and location
#' statistics. Defaults to \code{FALSE}.
#' @param U \emph{multivariate signs} of \code{data}, a matrix of size
#' \code{c(n, p - 1)}. Computed if \code{NULL} (the default).
#' @param V \emph{cosines} of \code{data}, a vector of size \code{n}. Computed
#' if \code{NULL} (the default).
#' @return An object of the \code{htest} class with the following elements:
#' \itemize{
#' \item \code{statistic}: test statistic.
#' \item \code{parameter}: degrees of freedom of the chi-square distribution
#' appearing in all the null asymptotic distributions.
#' \item \code{p.value}: \eqn{p}-value of the test.
#' \item \code{method}: information on the type of test performed.
#' \item \code{data.name}: name of the value of \code{data}.
#' \item \code{U}: multivariate signs of \code{data}.
#' \item \code{V}: cosines of \code{data}.
#' }
#' @details
#' Descriptions of the tests:
#' \itemize{
#' \item The "scatter" test is locally and asymptotically optimal against
#' \link[=tangent-elliptical]{tangent elliptical} alternatives to rotational
#' symmetry. However, it is not consistent against
#' \link[=tangent-vMF]{tangent von Mises--Fisher} (vMF) alternatives.
#' The asymptotic null distribution of
#' \eqn{Q_{\boldsymbol{\theta}}^{\mathrm{sc}}}{Q_\theta^{sc}}
#' is unaffected if \eqn{\boldsymbol{\theta}}{\theta} is estimated, that is,
#' the asymptotic null distributions of
#' \eqn{Q_{\boldsymbol{\theta}}^{\mathrm{sc}}}{Q_\theta^{sc}} and
#' \eqn{Q_{\hat{\boldsymbol{\theta}}}^{\mathrm{sc}}}{Q_{\hat \theta}^{sc}} are
#' the same.
#' \item The "location" test is locally and asymptotically most powerful
#' against vMF alternatives to rotational symmetry. However, it is not
#' consistent against tangent elliptical alternatives. The asymptotic
#' null distribution of
#' \eqn{Q_{\boldsymbol{\theta}}^{\mathrm{loc}}}{Q_\theta^{loc}}
#' for known \eqn{\boldsymbol{\theta}}{\theta} (the one implemented in
#' \code{test_rotasym}) \emph{does change} if \eqn{\boldsymbol{\theta}}{\theta}
#' is estimated by \eqn{\hat{\boldsymbol{\theta}}}{\hat \theta}. Therefore, if
#' the test is performed with an estimated \eqn{\boldsymbol{\theta}}{\theta}
#' (if \code{theta} is a function)
#' \eqn{Q_{\hat{\boldsymbol{\theta}}}^{\mathrm{loc}}}{Q_{\hat \theta}^{loc}}
#' will not be properly calibrated. \code{test_rotasym} will give a warning in
#' such case.
#' \item The "vMF location" test is a modification of the "location" test
#' designed to make its null asymptotic distribution invariant from the
#' estimation of \eqn{\boldsymbol{\theta}}{\theta} (as the "scatter" test is).
#' The test is optimal against tangent vMF alternatives with a \emph{specific},
#' vMF-based, angular function \code{\link{g_vMF}}. Despite not
#' being optimal against all tangent vMF alternatives, it is
#' consistent for all of them. As the location test,
#' it is not consistent against tangent elliptical alternatives.
#' \item The "hybrid" test combines (see below how) the "scatter" and
#' "location" tests. The test is neither optimal against tangent elliptical nor
#' tangent vMF alternatives, but it is consistent against both. Since it is
#' based on the "location" test, if computed with an estimator
#' \eqn{\hat{\boldsymbol{\theta}}}{\hat \theta}, the test statistic will not
#' be properly calibrated. \code{test_rotasym} will give a warning in such case.
#' \item The "vMF hybrid" test is the analogous of the "hybrid" test but
#' replaces the "location" test by the "vMF location" test.
#' }
#'
#' The combination of the scatter and location tests in the hybrid tests is
#' done in two different ways:
#' \itemize{
#' \item If \code{Fisher = FALSE}, then the scatter and location tests
#' statistics give the hybrid test statistic
#' \deqn{Q^{\mathrm{hyb}}:=Q_{\boldsymbol{\theta}}^{\mathrm{sc}}+
#' Q_{\boldsymbol{\theta}}^{\mathrm{loc}}.}{
#' Q_\theta^{hyb} := Q_\theta^{sc} + Q_\theta^{loc}.}
#' \item If \code{Fisher = TRUE}, then Fisher's method for aggregating
#' independent tests (the two test statistics are independent under rotational
#' symmetry) is considered, resulting the hybrid test statistic:
#' \deqn{Q_{\boldsymbol{\theta}}^{\mathrm{hyb}}
#' :=-2(\log(p_{\mathrm{sc}})+\log(p_{\mathrm{loc}}))}{
#' Q_\theta^{hyb} := -2(log(p_{sc}) + log(p_{loc}))}
#' where \eqn{p_{\mathrm{sc}}}{p_{sc}} and \eqn{p_{\mathrm{loc}}}{p_{loc}} are
#' the \eqn{p}-values of the scatter and location tests, respectively.
#' }
#' The hybrid test statistic \eqn{Q_{\mathrm{vMF}}^{\mathrm{hyb}}}{
#' Q_{vMF}^{hyb}} follows analogously to
#' \eqn{Q_{\boldsymbol{\theta}}^{\mathrm{hyb}}}{Q_\theta^{hyb}} by replacing
#' \eqn{Q_{\boldsymbol{\theta}}^{\mathrm{loc}}}{Q_\theta^{loc}} with
#' \eqn{Q_{\mathrm{vMF}}^{\mathrm{loc}}}{Q_{vMF}^{loc}}.
#'
#' Finally, recall that the tests are designed to test \emph{implications} of
#' rotational symmetry. Therefore, the tests are not consistent against
#' \emph{all} types of alternatives to rotational symmetry.
#' @author Eduardo García-Portugués, Davy Paindaveine, and Thomas Verdebout.
#' @references
#' García-Portugués, E., Paindaveine, D., Verdebout, T. (2020) On optimal tests
#' for rotational symmetry against new classes of hyperspherical distributions.
#' \emph{Journal of the American Statistical Association}, to appear.
#' \url{https://doi.org/10.1080/01621459.2019.1665527}
#' @examples
#' ## Rotational symmetry holds
#'
#' # Sample data from a vMF (rotational symmetric distribution about mu)
#' n <- 200
#' p <- 10
#' theta <- c(1, rep(0, p - 1))
#' set.seed(123456789)
#' data_0 <- r_vMF(n = n, mu = theta, kappa = 1)
#'
#' # theta known
#' test_rotasym(data = data_0, theta = theta, type = "sc")
#' test_rotasym(data = data_0, theta = theta, type = "loc")
#' test_rotasym(data = data_0, theta = theta, type = "loc_vMF")
#' test_rotasym(data = data_0, theta = theta, type = "hyb")
#' test_rotasym(data = data_0, theta = theta, type = "hyb", Fisher = TRUE)
#' test_rotasym(data = data_0, theta = theta, type = "hyb_vMF")
#' test_rotasym(data = data_0, theta = theta, type = "hyb_vMF", Fisher = TRUE)
#'
#' # theta unknown (employs the spherical mean as estimator)
#' test_rotasym(data = data_0, type = "sc")
#' test_rotasym(data = data_0, type = "loc") # Warning
#' test_rotasym(data = data_0, type = "loc_vMF")
#' test_rotasym(data = data_0, type = "hyb") # Warning
#' test_rotasym(data = data_0, type = "hyb", Fisher = TRUE) # Warning
#' test_rotasym(data = data_0, type = "hyb_vMF")
#' test_rotasym(data = data_0, type = "hyb_vMF", Fisher = TRUE)
#'
#' ## Rotational symmetry does not hold
#'
#' # Sample non-rotational symmetric data from a tangent-vMF distribution
#' # The scatter test is blind to these deviations, while the location tests
#' # are optimal
#' n <- 200
#' p <- 10
#' theta <- c(1, rep(0, p - 1))
#' mu <- c(rep(0, p - 2), 1)
#' kappa <- 2
#' set.seed(123456789)
#' r_V <- function(n) {
#'   r_g_vMF(n = n, p = p, kappa = 1)
#' }
#' data_1 <- r_TM(n = n, r_V = r_V, theta = theta, mu = mu, kappa = kappa)
#'
#' # theta known
#' test_rotasym(data = data_1, theta = theta, type = "sc")
#' test_rotasym(data = data_1, theta = theta, type = "loc")
#' test_rotasym(data = data_1, theta = theta, type = "loc_vMF")
#' test_rotasym(data = data_1, theta = theta, type = "hyb")
#' test_rotasym(data = data_1, theta = theta, type = "hyb", Fisher = TRUE)
#' test_rotasym(data = data_1, theta = theta, type = "hyb_vMF")
#' test_rotasym(data = data_1, theta = theta, type = "hyb_vMF", Fisher = TRUE)
#'
#' # theta unknown (employs the spherical mean as estimator)
#' test_rotasym(data = data_1, type = "sc")
#' test_rotasym(data = data_1, type = "loc") # Warning
#' test_rotasym(data = data_1, type = "loc_vMF")
#' test_rotasym(data = data_1, type = "hyb") # Warning
#' test_rotasym(data = data_1, type = "hyb", Fisher = TRUE) # Warning
#' test_rotasym(data = data_1, type = "hyb_vMF")
#' test_rotasym(data = data_1, type = "hyb_vMF", Fisher = TRUE)
#'
#' # Sample non-rotational symmetric data from a tangent-elliptical distribution
#' # The location tests are blind to these deviations, while the
#' # scatter test is optimal
#' n <- 200
#' p <- 10
#' theta <- c(1, rep(0, p - 1))
#' Lambda <- matrix(0.5, nrow = p - 1, ncol = p - 1)
#' diag(Lambda) <- 1
#' set.seed(123456789)
#' r_V <- function(n) {
#'   r_g_vMF(n = n, p = p, kappa = 1)
#' }
#' data_2 <- r_TE(n = n, r_V = r_V, theta = theta, Lambda = Lambda)
#'
#' # theta known
#' test_rotasym(data = data_2, theta = theta, type = "sc")
#' test_rotasym(data = data_2, theta = theta, type = "loc")
#' test_rotasym(data = data_2, theta = theta, type = "loc_vMF")
#' test_rotasym(data = data_2, theta = theta, type = "hyb")
#' test_rotasym(data = data_2, theta = theta, type = "hyb", Fisher = TRUE)
#' test_rotasym(data = data_2, theta = theta, type = "hyb_vMF")
#' test_rotasym(data = data_2, theta = theta, type = "hyb_vMF", Fisher = TRUE)
#'
#' # theta unknown (employs the spherical mean as estimator)
#' test_rotasym(data = data_2, type = "sc")
#' test_rotasym(data = data_2, type = "loc") # Warning
#' test_rotasym(data = data_2, type = "loc_vMF")
#' test_rotasym(data = data_2, type = "hyb") # Warning
#' test_rotasym(data = data_2, type = "hyb", Fisher = TRUE) # Warning
#' test_rotasym(data = data_2, type = "hyb_vMF")
#' test_rotasym(data = data_2, type = "hyb_vMF", Fisher = TRUE)
#'
#' ## Sunspots births data
#'
#' # Load data
#' data("sunspots_births")
#' sunspots_births$X <-
#'   cbind(cos(sunspots_births$phi) * cos(sunspots_births$theta),
#'         cos(sunspots_births$phi) * sin(sunspots_births$theta),
#'         sin(sunspots_births$phi))
#'
#' # Test rotational symmetry for the 23rd cycle, specified theta
#' sunspots_23 <- subset(sunspots_births, cycle == 23)
#' test_rotasym(data = sunspots_23$X, type = "sc", theta = c(0, 0, 1))
#' test_rotasym(data = sunspots_23$X, type = "loc", theta = c(0, 0, 1))
#' test_rotasym(data = sunspots_23$X, type = "hyb", theta = c(0, 0, 1))
#'
#' # Test rotational symmetry for the 23rd cycle, unspecified theta
#' spherical_loc_PCA(sunspots_23$X)
#' test_rotasym(data = sunspots_23$X, type = "sc", theta = spherical_loc_PCA)
#' test_rotasym(data = sunspots_23$X, type = "loc_vMF",
#'              theta = spherical_loc_PCA)
#' test_rotasym(data = sunspots_23$X, type = "hyb_vMF",
#'              theta = spherical_loc_PCA)
#'
#' # Test rotational symmetry for the 22nd cycle, specified theta
#' sunspots_22 <- subset(sunspots_births, cycle == 22)
#' test_rotasym(data = sunspots_22$X, type = "sc", theta = c(0, 0, 1))
#' test_rotasym(data = sunspots_22$X, type = "loc", theta = c(0, 0, 1))
#' test_rotasym(data = sunspots_22$X, type = "hyb", theta = c(0, 0, 1))
#'
#' # Test rotational symmetry for the 22nd cycle, unspecified theta
#' spherical_loc_PCA(sunspots_22$X)
#' test_rotasym(data = sunspots_22$X, type = "sc", theta = spherical_loc_PCA)
#' test_rotasym(data = sunspots_22$X, type = "loc_vMF",
#'              theta = spherical_loc_PCA)
#' test_rotasym(data = sunspots_22$X, type = "hyb_vMF",
#'              theta = spherical_loc_PCA)
#' @seealso \code{\link{tangent-elliptical}}, \code{\link{tangent-vMF}},
#' \code{\link{spherical_mean}}.
#' @export
test_rotasym <- function(data, theta = spherical_mean,
                         type = c("sc", "loc", "loc_vMF", "hyb", "hyb_vMF")[5],
                         Fisher = FALSE, U = NULL, V = NULL) {

  # Data name
  data_name <- deparse(substitute(data))

  # Remove NA's
  data <- data[complete.cases(data), , drop = FALSE]

  # Check if data is spherical
  data <- check_unit_norm(x = data, warnings = TRUE)

  # type to lowercase
  type <- tolower(type)

  # Compute theta if not provided
  if (is.function(theta)) {

    # Warning for tests affected by the estimation of theta
    if (type == "loc" | type == "hyb") {

      warning(paste0("\"loc\" and \"hyb\" tests are not properly ",
                     "calibrated when theta is estimated.\n Consider using ",
                     "\"loc_vMF\" and \"hyb_vMF\" tests instead for accurate ",
                     "results."))

    }

    # Estimate theta
    theta <- theta(data)

  }

  # Check unit norm of theta
  theta <- check_unit_norm(x = theta, warnings = TRUE)

  # Compute U's
  if (is.null(U)) {

    U <- signs(X = data, theta = theta, check_X = FALSE)

  }

  # Compute V's
  if (type == "loc_vmf" & is.null(V)) {

    V <- cosines(X = data, theta = theta, check_X = FALSE)

  }

  # Dimensions
  n <- nrow(U)
  p <- ncol(U) + 1

  # Kind of test
  if (type == "sc") {

    # Test
    U_bar <- colMeans(U)
    S <- (n - 1) / n * cov(U) + tcrossprod(U_bar)
    statistic <- 0.5 * n * (p * p - 1) * (sum(diag(S %*% S)) - 1 / (p - 1))
    df <- 0.5 * (p - 2) * (p + 1)
    p_value <- pchisq(q = statistic, df = df, lower.tail = FALSE)

  } else if (type == "loc") {

    # Test
    statistic <- n * (p - 1) * sum(colMeans(U)^2)
    df <- p - 1
    p_value <- pchisq(q = statistic, df = df, lower.tail = FALSE)

  } else if (type == "loc_vmf") {

    # Quantities for computing the statistic
    V2 <- 1 - V * V
    V_sqrt <- sqrt(V2)
    V_inv_sqrt <- 1 / sqrt(V2)
    V_inv_sqrt[!is.finite(V_inv_sqrt)] <- NA
    D_pg <- (p - 2) * mean(V * V_inv_sqrt, na.rm = TRUE) / ((p - 1) * mean(V))

    # Delta
    Delta <- colSums((1 - D_pg * V_sqrt) * U) / sqrt(n)

    # Gamma
    inv_Gamma <- diag((p - 1) / (1 - 2 * D_pg * mean(V_sqrt) +
                                   D_pg * D_pg * mean(V2)),
                      nrow = p - 1, ncol = p - 1)

    # Test
    statistic <- drop(t(Delta) %*% inv_Gamma %*% Delta)
    df <- p - 1
    p_value <- pchisq(q = statistic, df = df, lower.tail = FALSE)

  } else if (type == "hyb" | type == "hyb_vmf") {

    # Individual tests
    type_loc <- switch(type, "hyb" = "loc", "hyb_vmf" = "loc_vMF")
    test_sc <- test_rotasym(data = data, theta = theta, type = "sc", U = U)
    test_loc <- test_rotasym(data = data, theta = theta, type = type_loc,
                             U = U, V = V)

    # Test considering Fisher's method or addition of test statistics
    if (Fisher) {

      statistic <- -2 * (log(test_sc$p.value) + log(test_loc$p.value))
      df <- 4

    } else {

      statistic <- test_sc$statistic + test_loc$statistic
      df <- test_sc$parameter + test_loc$parameter

    }
    p_value <- pchisq(q = statistic, df = df, lower.tail = FALSE)

  } else {

    stop(paste0("Incorrect type. Must be one of the following: \"sc\", ",
                "\"loc\", \"loc_vMF\", \"hyb\", \"hyb_vMF\""))

  }

  # Result as htest class
  names(statistic) <- switch(type,
                             "sc" = "Q_sc",
                             "loc" = "Q_loc",
                             "loc_vmf" = "Q_loc_vMF",
                             "hyb" = "Q_hyb",
                             "hyb_vmf" = "Q_hyb_vMF")
  names(p_value) <- names(statistic)
  method <- switch(type,
    "sc" = "Scatter test for rotational symmetry",
    "loc" = "Location test for rotational symmetry",
    "loc_vmf" = "Location vMF test for rotational symmetry",
    "hyb" =
      ifelse(Fisher,
             "Hybrid test (Fisher's method) for rotational symmetry",
             "Hybrid test (addition of statistics) for rotational symmetry"),
    "hyb_vmf" =
      ifelse(Fisher,
             "Hybrid vMF test (Fisher's method) for rotational symmetry",
             paste("Hybrid vMF test (addition of statistics) for",
                   "rotational symmetry")))
  names(df) <- "df"
  test <- list(statistic = statistic, parameter = df, p.value = p_value,
               method = method, data.name = data_name, U = U, V = V)
  class(test) <- "htest"
  return(test)

}
