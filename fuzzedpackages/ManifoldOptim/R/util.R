printf <- function(msg, ...) { cat(sprintf(msg, ...)) }

`%notin%` <- Negate(`%in%`)

#' Print summary from \code{manifold.optim} results
#'
#' Print results
#' @param x A \code{ManifoldOptim} object output by \code{manifold.optim}.
#' @param ... Not currently used
print.ManifoldOptim <- function(x, ...)
{
	cat("-----", "ManifoldOptim result", "-----\n")
	cat("xopt: <length ", length(x$xopt), ">\n", sep = "")
	cat("fval:", x$fval, "\n")
	cat("normgf:", x$normgf, "\n")
	cat("normgfgf0:", x$normgfgf0, "\n")
	cat("iter:", x$iter, "\n")
	cat("num.obj.eval:", x$num.obj.eval, "\n")
	cat("num.grad.eval:", x$num.grad.eval, "\n")
	cat("nR:", x$nR, "\n")
	cat("nV:", x$nV, "\n")
	cat("nVp:", x$nVp, "\n")
	cat("nH:", x$nH, "\n")
	cat("elapsed:", x$elapsed, "\n")

	if (nchar(x$message) > 0) {
		cat("NOTE:", x$message, "\n")
	}
	# head(x$funSeries)
	# head(x$gradSeries)
	# head(x$timeSeries)

	return(invisible(x))
}

#' Compute the trace of a square matrix
#'
#' @param X A matrix
Trace <- function(X)
{	
	if (!is.matrix(X)) stop("Argument to Trace is not a matrix")
	stopifnot(nrow(X) == ncol(X))
	return(sum(diag(X))) 
}

#' Orthonormalize the columns of a matrix
#'
#' @param u A matrix
orthonorm <- function(u) 
{
	if (is.null(u)) return(NULL)
	if (!(is.matrix(u))) u <- as.matrix(u)
	if (sum(is.nan(u) > 0)) stop("Matrix contains NaNs")
	if (sum(is.na(u) > 0)) stop("Matrix contains NAs")
	dd <- dim(u); n <- dd[1]; p <-dd[2]
	
	if (prod(abs(La.svd(u)$d) > 1e-08) == 0) stop("collinear vectors in orthonorm")
	if (n < p)
	{
		warning("There are too many vectors to orthonormalize in orthonorm.")
		u <- as.matrix(u[, 1:p])
		n <- p
	}
	v <- u
	if (p > 1)
	{
		for (i in 2:p)
		{
			coef.proj <- c(crossprod(u[, i], v[, 1:(i - 1)]))/diag(crossprod(v[, 1:(i - 1)]));
			v[, i] <- u[, i] - matrix(v[, 1:(i - 1)], nrow = n) %*% matrix(coef.proj, nrow = i - 1)
		}
	}
	coef.proj <- 1/sqrt(diag(crossprod(v)))
	return(t(t(v) * coef.proj))
}

