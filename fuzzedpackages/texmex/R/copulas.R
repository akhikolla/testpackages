#' Compute empirical distribution function
#'
#' Compute the empirical distribution function
#'
#'
#' @usage edf(x, na.last = NA)
#' @param x A numeric vector
#' @param na.last How to treat missing values. See \code{\link{rank}} for
#' details.
#' @return A vector of quantiles relating to the observations in \code{x}.
#' @author Harry Southworth
#' @seealso \code{\link{copula}}
#' @keywords univar
#' @examples
#'
#' plot(winter$NO, edf(winter$NO))
#'
#' @export edf
edf <- function(x, na.last=NA){
    res <- rank(x) / (length(x) + 1)
    oldClass(res) <- "edf"
    invisible(res)
}


#' Calculate the copula of a matrix of variables
#'
#' Returns the copula of several random variables.
#'
#' The result is obtained by applying \code{\link{edf}} to each column of
#' \code{x} in turn.
#'
#' Print and plot methods are available for the copula class.
#'
#' @param x A matrix or data.frame containing numeric variables.
#' @param na.last How to treat missing values. See \code{rank} for details.
#' @param ... further arguments
#' @return A matrix with the same dimensions as \code{x}, each column of which
#' contains the quantiles of each column of \code{x}. This object is of class
#' \code{copula}.
#' @author Harry Southworth
#' @seealso \code{\link{edf}} \code{\link{plot.copula}} \code{\link{ggplot.copula}}
#' @keywords multivariate
#' @examples
#'
#'   D <- liver[liver$dose == "D",]
#'   Dco <- copula(D)
#'   plot(Dco)
#'
#' @export copula
copula <- function(x, na.last=NA, ...) {
    UseMethod("copula")
}

#' @describeIn copula default method
#' @export
copula.default <- function(x, na.last=NA, ...) {
    stop("Can't calculate copula")
}

#' @describeIn copula data frame method
#' @export
copula.data.frame <- function(x, na.last=NA, ...) {
    theCall <- match.call()

    really.numeric <- function(x){
        class(x) %in% c("integer", "numeric")
    }

    wh <- sapply(x, really.numeric)

    if (sum(wh) == 0){
        stop("x contains no numeric columns")
    }

    if (sum(wh) < length(wh)){
        warning(paste("Some variables have been dropped:", paste(colnames(x)[!wh], collapse=", ")))
    }

    result <- copula(as.matrix(x[, wh]), na.last=na.last)
    result$call <- theCall
    result
}

#' @describeIn copula matrix method
#' @export
copula.matrix <- function (x, na.last = NA, ...) {
    theCall <- match.call()

    res <- apply(x, 2, edf)

    res <- list(call=theCall, copula=res)
    oldClass(res) <- "copula"
    res
}


#' @export
print.copula <- function(x, ...){
    print(x$call)
    cat("A copula of", ncol(x$copula), "variables.\n")
    invisible(x)
}

#' Plot copulas
#'
#' @param x A copula object
#' @param jitter. If \code{jitter=TRUE}, the values are jittered
#'     before plotting. Defaults to \code{jitter. = FALSE}.
#' @param jitter.factor How much jittering to use. Defaults to
#'     \code{jitter.factor = 1}.
#' @param ... Other arguments to pass through to \code{plot}.
#' @export
plot.copula <- function(x, jitter. = FALSE, jitter.factor=1, ...){
    x <- x$copula

    thecall <- match.call()
    jitter. <- FALSE
    if (is.element("jitter.", names(thecall))){
    	jitter. <- thecall[["jitter."]]
    }

	if (jitter.){
		x <- apply(x, 2, jitter, factor=jitter.factor)
	}
    pairs(x, ...)
    invisible()
}


#' Fancy plotting for copulas
#'
#' @param jitter If \code{jitter=TRUE}, the values are jittered
#'     before plotting. Defaults to \code{jitter. = FALSE}.
#' @param jitter.factor How much jittering to use. Defaults to
#'     \code{jitter.factor = .05}.
#' @param data A data.frame.
#' @param mapping Not used.
#' @param color Defaults to \code{color = "blue"}.
#' @param alpha Defaults to \code{alpha = 0.7}.
#' @param point.size Defaults to \code{point.size = 1}.
#' @param smooth Defaults to \code{smooth = FALSE}.
#' @param smooth.method Defaults to \code{smooth.method = "auto"} and is passed
#'   to \code{geom_smooth} only when \code{smooth = TRUE}.
#' @param smooth.se Defaults to \code{smooth.se = TRUE} and is used only when
#'   \code{smooth = TRUE}.
#' @param smooth.level Defaults to \code{smooth.level = 0.95} and is used only
#'   when \code{smooth = TRUE}.
#' @param smooth.formula A formula, defaulting to \code{smooth.formula = y ~ x}
#'   to be passed as the \code{formula} argument to \code{geom_smooth}.
#' @param legend.position Passed into \code{theme}, defaults to \code{legend.position="none"}.
#' @param legend.title Passed into \code{theme}. Defaults to \code{legend.title = waiver()}.
#' @param diag Defaults to \code{diag = FALSE} and panels on the diagonal are not
#'   produced.
#' @param lower Defaults to \code{lower = TRUE} and only the lower triangle is plotted.
#' @param ticks Defaults to \code{ticks = TRUE} and ticks and their labels are put
#'   on the axes. Otherwise, no tick or labels are used.
#' @param environment Not used.
#' @param ... Not used.
#' @export
ggplot.copula <-
  function (data, mapping = aes(), color = "blue",
            alpha = 0.7, jitter = FALSE, jitter.factor = 0.05, point.size = 1,
            smooth = FALSE, smooth.method = "auto", smooth.se = TRUE,
            smooth.level = 0.95, smooth.formula = y ~ x, legend.position = "none",
            legend.title = ggplot2::waiver(), diag = FALSE, lower = TRUE,
            ticks = TRUE, ..., environment = parent.frame()) {

    data <- as.data.frame(data$copula)

    lvls <- names(data)

    data$.XXidXX. <- 1:nrow(data)

    ljdata <- data.frame(.XXidXX. = data[, ".XXidXX."], stringsAsFactors = FALSE)

    yy <- tidyr::gather(data, H, xval, -.XXidXX.)
    yy$H <- factor(yy$H, levels = lvls)
    ww <- yy
    names(ww) <- c(".XXidXX.", "V", "yval")

    zz <- dplyr::left_join(yy, ww, by = ".XXidXX.")
    zz <- dplyr::left_join(zz, ljdata, by = ".XXidXX.")

    if (!diag) {
      zz <- zz[zz$H != zz$V, ]
    }
    if (lower){
      zz <- zz[as.numeric(zz$H) <= as.numeric(zz$V), ]
    }
    if (jitter) {
      jw <- jh <- jitter.factor
    }
    else {
      jw <- jh <- 0
    }

    res <- ggplot(zz, aes(xval, yval))

    if (smooth) {
      res <- res + geom_smooth(mapping = aes(xval, yval), inherit.aes = FALSE,
                               method = smooth.method, formula = smooth.formula,
                               se = smooth.se, level = smooth.level)
    }


    res <- res + geom_point(color = color, size = point.size,
                            alpha = alpha,
                            position = position_jitter(height = jh, width = jw))

    res <- res + ggplot2::facet_grid(V ~ H, scales = "free") +
      labs(x = "", y = "")
    if (!ticks) {
      res <- res + theme(axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
                         axis.ticks.x = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
    }
    if (length(unique(ljdata$color)) > 1) {
      res <- res + theme(legend.position = legend.position,
                         legend.title = legend.title)
    }

    res
  }
