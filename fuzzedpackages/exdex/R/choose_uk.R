# ================================= choose_uk =================================

#' Threshold \eqn{u} and runs parameter \eqn{K} diagnostic for the \eqn{K}-gaps
#' estimator
#'
#' Creates data for a plot to aid the choice of the threshold and
#' run parameter \eqn{K} for the \eqn{K}-gaps estimator (see
#' \code{\link{kgaps}}).  \code{\link{plot.choose_uk}} creates the plot.
#'
#' @param data A numeric vector of raw data.  No missing values are allowed.
#' @param u,k Numeric vectors.  \code{u} is a vector of
#'   extreme value thresholds applied to data.  \code{k} is a vector of values
#'   of the run parameter \eqn{K}, as defined in Suveges and Davison (2010).
#'   See \code{\link{kgaps}} for more details.
#' @details For each combination of threshold in \code{u} and \eqn{K}
#'   in \code{k} the functions \code{\link{kgaps}} and \code{\link{kgaps_imt}}
#'   are called in order to estimate \eqn{\theta} and to perform the
#'   information matrix test of Suveges and Davison (2010).
#' @references Suveges, M. and Davison, A. C. (2010) Model
#'   misspecification in peaks over threshold analysis, \emph{The Annals of
#'   Applied Statistics}, \strong{4}(1), 203-221.
#'   \url{https://doi.org/10.1214/09-AOAS292}
#' @return An object (a list) of class \code{c("choose_uk", "exdex")}
#'   containing
#'   \item{imt }{an object of class \code{c("kgaps_imt", "exdex")} returned
#'     from \code{\link{kgaps_imt}}.}
#'   \item{theta }{a \code{length(u)} by \code{length(k)} matrix.
#'     Element (i,j) of \code{theta} contains an object (a list) of class
#'     \code{c("kgaps", "exdex")}, a result of a call
#'     \code{kgaps(data, u[j], k[i])} to \code{\link{kgaps}}.}
#' @seealso \code{\link{kgaps}} for maximum likelihood estimation of the
#'   extremal index \eqn{\theta} using the \eqn{K}-gaps model.
#' @seealso \code{\link{kgaps_imt}} for the information matrix test under the
#'   \eqn{K}-gaps model
#' @seealso \code{\link{plot.choose_uk}} to produce the diagnostic plot.
#' @examples
#' ### S&P 500 index
#'
#' # Multiple thresholds and run parameters
#' u <- quantile(sp500, probs = seq(0.1, 0.9, by = 0.1))
#' imt_theta <- choose_uk(sp500, u = u, k = 1:5)
#' plot(imt_theta)
#' plot(imt_theta, uprob = TRUE)
#' plot(imt_theta, y = "theta")
#'
#' # One run parameter K, many thresholds u
#' u <- quantile(sp500, probs = seq(0.1, 0.9, by = 0.1))
#' imt_theta <- choose_uk(sp500, u = u, k = 1)
#' plot(imt_theta)
#' plot(imt_theta, y = "theta")
#'
#' # One threshold u, many run parameters K
#' u <- quantile(sp500, probs = 0.9)
#' imt_theta <- choose_uk(sp500, u = u, k = 1:5)
#' plot(imt_theta)
#' plot(imt_theta, y = "theta")
#'
#' ### Newlyn sea surges
#'
#' u <- quantile(newlyn, probs = seq(0.1, 0.9, by = 0.1))
#' imt_theta <- choose_uk(newlyn, u = u, k = 1:5)
#' plot(imt_theta, uprob = TRUE)
#' @export
choose_uk <- function(data, u, k = 1) {
  n_u <- length(u)
  n_k <- length(k)
  theta <- matrix(rep(list(), n_u * n_k), n_k, n_u)
  # Function to set the correct element of the matrix of lists theta
  # i indexes k, j indexes u
  comp <- function(i, j) {
    return((i - 1) * n_u + j)
  }
  for (i in 1:n_k) {
    for (j in 1:n_u) {
      theta[[comp(i, j)]] <- kgaps(data, u[j], k[i])
    }
  }
  imt <- kgaps_imt(data, u, k)
  res <- list(imt = imt, theta = theta)
  class(res) <- c("choose_uk", "exdex")
  return(res)
}

# ============================= plot.choose_uk ===============================

#' Plot Threshold \eqn{u} and runs parameter \eqn{K} diagnostic for the
#' \eqn{K}-gaps estimator
#'
#' \code{plot} method for objects inheriting from class \code{"choose_uk"},
#' returned from \code{\link{choose_uk}}
#'
#' @param x an object of class \code{c("choose_uk", "exdex")}, a result of a
#'   call to \code{\link{choose_uk}}.
#' @param y A character scalar indicating what should be plotted on the
#'   vertical axes of the plot: information matrix test statistics (IMTS)
#'   if \code{y = "imts"} and estimates of \eqn{\theta} if \code{y = "theta"}.
#'   If \code{y = "theta"}, and either \code{x$u} or \code{x$k} have length
#'   one, then 100\code{level}\% confidence intervals are added to the plot.
#' @param level A numeric scalar in (0, 1).  The confidence level used in
#'   calculating confidence intervals for \eqn{\theta}.  Only relevant if
#'   \code{y = "theta"} and either \code{x$u} or \code{x$k} have length one.
#' @param interval_type A character scalar.  The type of confidence interval
#'   to be plotted, if \code{y = "theta"}.  See \code{\link{confint.kgaps}}.
#' @param conf_scale A character scalar.  If \code{interval_type = "norm"} then
#'   \code{conf_scale} determines the scale on which we use approximate
#'   large-sample normality of the estimator to estimate confidence intervals.
#'   See \code{\link{confint.kgaps}}.
#' @param alpha A numeric vector with entries in (0, 1). The size of the test
#'   to be performed.
#' @param constrain A logical scalar.  The argument \code{constrain} to
#'  \code{\link{confint.kgaps}}.
#' @param for_abline Only relevant when \code{y = "imts"} and at one of
#'   \code{u} or \code{k} is scalar. A list of graphical parameters to be
#'   passed to \code{\link{abline}} to indicate the critical value of the
#'   information matrix test (IMT) implied by \code{alpha}.
#' @param digits An integer. Used for formatting the value of the threshold
#'   with \code{\link[base:Round]{signif}} before adding its value to a plot.
#' @param uprob A logical scalar. Should we plot \code{x$u} on the
#'   horizontal axis (\code{uprob = FALSE}) or the approximate sample quantile
#'   to which \code{x$u} corresponds (\code{uprob = FALSE})?
#' @param leg_pos A character scalar.  The position of any legend added to
#'   a plot.  Only relevant when both the arguments \code{u} and \code{k}
#'   in the call to \code{\link{choose_uk}} have length greater than one.
#' @param ... Additional arguments passed to \code{\link[graphics]{matplot}}.
#' @details The type of plot produced depends mainly on \code{y}.
#'
#'   If \code{y = "imts"} then the values of IMTS are plotted against the
#'   thresholds in \code{x$u} (or their corresponding approximate sample
#'   quantile levels if \code{uprob = TRUE}) for each value of \eqn{K}
#'   in \code{x$k}.  Horizontal lines are added to indicate the critical
#'   values of the IMT for the significance levels in \code{alpha}.
#'   We would not reject at the 100\code{alpha}\% level combinations of
#'   threshold and \eqn{K} corresponding to values of the IMTS that fall
#'   below the line.
#'
#'   If \code{y = "theta"} then estimates of \eqn{\theta} are plotted on the
#'   vertical axis.  If both \code{x$u} and \code{x$k$} have length greater
#'   than one then only these estimates are plotted.  If either \code{x$u}
#'   or \code{x$k} have length one then approximate 100\code{level}\%
#'   confidence intervals are added to the plot and the variable,
#'   \code{x$u} or \code{x$k} that has length greater than one is plotted on
#'   the horizontal axis.
#' @return Nothing is returned.
#' @seealso \code{\link{choose_uk}}.
#' @section Examples:
#' See the examples in \code{\link{choose_uk}}.
#' @export
plot.choose_uk <- function(x, y = c("imts", "theta"), level = 0.95,
                           interval_type = c("norm", "lik"),
                           conf_scale = c("theta", "log"), alpha = 0.05,
                           constrain = TRUE,
                           for_abline = list(lty = 2, lwd = 1, col = 1),
                           digits = 3, uprob = FALSE,
                           leg_pos = if (y == "imts") "topright" else "topleft",
                           ...) {
  y <- match.arg(y)
  interval_type <- match.arg(interval_type)
  conf_scale <- match.arg(conf_scale)
  # Extract the values of k and u
  k <- x$imt$k
  u <- x$imt$u
  n_k <- length(k)
  n_u <- length(u)
  # Approximate sample quantiles of threshold u
  u_ps <- rownames(x$imt$imt)
  if (n_k == 1 && n_u == 1) {
    stop("Object contains only 1 threshold and 1 value of K")
  }
  # Function to set the correct element of the matrix of lists theta
  # i indexes k, j indexes u
  comp <- function(i, j) {
    return((i - 1) * n_u + j)
  }
  # My plotting functions: to give defaults but allow the user to override
  my_matplot <- function(x, y, ..., type = "l", lty = my_lty, col = my_col,
                         xlab = my_xlab, ylab = my_ylab, ylim = my_ylim,
                         lwd = my_lwd) {
    graphics::matplot(x = x, y = y, ..., type = type, lty = lty, col = col,
                      xlab = xlab, ylab = ylab, ylim = ylim, lwd = lwd)
  }
  my_title <- function(..., main = my_main) {
    graphics::title(..., main = main)
  }
  # Critical value for the IMT (in case we need it)
  crit <- stats::qchisq(alpha, df = 1, lower.tail = FALSE)
  # One of k or u is scalar
  cond1 <- n_k == 1 && n_u > 1
  cond2 <- n_k > 1 && n_u == 1
  my_lwd <- 1
  if (cond1 || cond2)  {
    my_col <- my_lty <- 1
    max_uk <- max(n_k, n_u)
    ymat <- matrix(NA, ncol = 3, nrow = max_uk)
    if (cond1) {
      if (uprob) {
        xvec <- u_ps
      } else {
        xvec <- u
      }
    } else {
      xvec <- k
    }
    loop_vec <- 1:max_uk
    if (y == "theta") {
      for (ij in loop_vec) {
        if (cond1) {
          kgaps_object <- x$theta[[comp(1, ij)]]
        } else {
          kgaps_object <- x$theta[[comp(ij, 1)]]
        }
        temp <- confint(kgaps_object, level = level,
                        interval_type = interval_type,
                        conf_scale = conf_scale, constrain = constrain)
        ymat[ij, 1] <- kgaps_object$theta
        ymat[ij, 2:3] <- temp
      }
      my_ylab <- "theta"
      my_xlab <- ifelse(cond1, "threshold u", "run parameter K")
      if (uprob & cond1) {
        my_xlab <- "sample quantile level of threshold u"
      }
      my_ylim <- c(0, 1)
      my_matplot(xvec, ymat, ...)
      my_main <- ifelse(cond1, paste0("run parameter K = ", k),
                        paste0("threshold u = ", signif(u, digits = digits),
                               " (", u_ps, "% quantile)"))
      my_title(...)
    } else {
      my_ylab <- "IMTS"
      my_xlab <- ifelse(cond1, "threshold u", "run parameter K")
      if (uprob & cond1) {
        my_xlab <- "sample quantile level of threshold u"
      }
      my_ylim <- c(0, max(x$imt$imt, na.rm = TRUE))
      if (cond1) {
        if (uprob) {
          xvec <- u_ps
        } else {
          xvec <- x$imt$u
        }
        ymat <- x$imt$imt
      } else {
        xvec <- x$imt$k
        ymat <- t(x$imt$imt)
      }
      my_matplot(xvec, ymat, ...)
      my_main <- ifelse(cond1, paste0("run parameter K = ", k),
                        paste0("threshold u = ", signif(u, digits = digits),
                               " (", u_ps, "% quantile)"))
      my_title(...)
      for_abline <- c(for_abline, list(h = crit))
      do.call(graphics::abline, for_abline)
      graphics::mtext(as.character(alpha), 4, at = crit, las = 1, cex = 0.8,
                      adj = 1, padj = 1)
    }
  } else {
    my_lty <- 1
    my_col <- 1:n_k
    if (uprob) {
      xvec <- u_ps
    } else {
      xvec <- x$imt$u
    }
    if (y == "theta") {
      ymat <- x$imt$theta
      my_ylab <- "theta"
      my_ylim <- c(0, 1)
    } else {
      ymat <- x$imt$imt
      my_ylab <- "IMTS"
      my_ylim <- c(0, max(x$imt$imt, na.rm = TRUE))
    }
    if (uprob) {
      my_xlab <- "sample quantile level of threshold u"
    } else {
      my_xlab <- "threshold u"
    }
    my_matplot(xvec, ymat, ...)
    if (y == "imts") {
      for_abline <- c(for_abline, list(h = crit))
      do.call(graphics::abline, for_abline)
      graphics::mtext(as.character(alpha), 4, at = crit, las = 1, cex = 0.8,
                      adj = 1, padj = 1)
    }
    user_args <- list(...)
    if (is.null(user_args$lty)) {
      leg_lty <- my_lty
    } else {
      leg_lty <- user_args$lty
    }
    if (is.null(user_args$col)) {
      leg_col <- my_col
    } else {
      leg_col <- user_args$col
    }
    if (is.null(user_args$lwd)) {
      leg_lwd <- my_lwd
    } else {
      leg_lwd <- user_args$lwd
    }
    graphics::legend(leg_pos, legend = paste0("K = ", k), lty = leg_lty,
                     col = leg_col, lwd = leg_lwd)
  }
  return(invisible())
}
