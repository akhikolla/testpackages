
#' Plotting method for vennLasso fitted objects
#' 
#' @param x fitted \code{vennLasso} or \code{cv.vennLasso} model object
#' @param which.subpop which row in the coefficient matrix should be plotting? Each row corresponds to
#' a particular combination of the specified stratifying variables
#' @param xvar What is on the X-axis. "norm" plots against the L1-norm of the coefficients, "lambda" against the log-lambda sequence, and "dev"
#' against the percent deviance explained.
#' @param xlab character value supplied for x-axis label
#' @param ylab character value supplied for y-axis label
#' @param ... other graphical parameters for the plot
#' @rdname plot
#' @importFrom graphics matplot
#' @importFrom stats na.omit
#' @export
#' @examples
#' library(Matrix)
#' 
#' set.seed(123)
#' n.obs <- 200
#' n.vars <- 50
#'
#' true.beta.mat <- array(NA, dim = c(3, n.vars))
#' true.beta.mat[1,] <- c(-0.5, -1, 0, 0, 2, rep(0, n.vars - 5))
#' true.beta.mat[2,] <- c(0.5, 0.5, -0.5, -0.5, 1, -1, rep(0, n.vars - 6))
#' true.beta.mat[3,] <- c(0, 0, 1, 1, -1, rep(0, n.vars - 5))
#' rownames(true.beta.mat) <- c("1,0", "1,1", "0,1")
#' true.beta <- as.vector(t(true.beta.mat))
#'
#' x.sub1 <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
#' x.sub2 <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
#' x.sub3 <- matrix(rnorm(n.obs * n.vars), n.obs, n.vars)
#'
#' x <- as.matrix(rbind(x.sub1, x.sub2, x.sub3))
#'
#' conditions <- as.matrix(cbind(c(rep(1, 2 * n.obs), rep(0, n.obs)),
#'                               c(rep(0, n.obs), rep(1, 2 * n.obs))))
#'
#' y <- rnorm(n.obs * 3, sd = 3) + drop(as.matrix(bdiag(x.sub1, x.sub2, x.sub3)) %*% true.beta)
#'
#' fit <- vennLasso(x = x, y = y, groups = conditions)
#'
#' layout(matrix(1:3, ncol = 3))
#' plot(fit, which.subpop = 1)
#' plot(fit, which.subpop = 2)
#' plot(fit, which.subpop = 3)
#'
plot.vennLasso <- function(x, which.subpop = 1,
                           xvar = c("norm", "lambda", "loglambda", "dev"),
                           xlab = iname, ylab = "Coefficients",
                           ...)
{
    num.condition.combinations <- dim(x$beta)[1]
    if (which.subpop > num.condition.combinations)
    {
        err.txt <- paste0("Row ", which.subpop, " specified, but only ",
                          num.condition.combinations,
                          " condition combinations present in the model.")
        stop(err.txt)
    }
    xvar <- match.arg(xvar)
    nbeta <- as.matrix(x$beta[which.subpop,,])
    switch(xvar,
           "norm" = {
               index <- apply(abs(nbeta), 2, sum)
               iname <- "L1 Norm"
               xlim <- range(index)
           },
           "lambda" = {
               index <- x$lambda
               iname <- "Lambda"
               xlim <- rev(range(index))
           },
           "loglambda" = {
               index <- log(x$lambda)
               iname <- "Log Lambda"
               xlim <- rev(range(index))
           },
           "dev" = {
               index = x$sumSquare
               iname = "Sum of Squares"
               xlim <- range(index)
           }
    )
    rn <- rownames(x$beta)
    if (!is.null(rn))
    {
        main.txt <- rn[which.subpop]
    } else
    {
        main.txt <- ""
    }
    matplot(index, t(nbeta), lty = 1,
            xlab = xlab, ylab = ylab,
            xlim = xlim, main = main.txt,
            type = 'l', ...)
}
