#' log likelihood function for fitted vennLasso objects
#'
#' @param object fitted "vennLasso" model object.
#' @param ... not used
#' @rdname logLik
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
#' logLik(fit)
#'
logLik.vennLasso <- function(object, ...) {
    # taken from ncvreg. Thanks to Patrick Breheny.
    n <- as.numeric(object$nobs)
    df <- object$df

    if (object$family == "gaussian")
    {

        resid.ss <- object$loss
        logL <- -0.5 * n * (log(2 * pi) - log(n) + log(resid.ss)) - 0.5 * n
    } else if (object$family == "binomial")
    {
        logL <- -1 * object$loss
    } else if (object$family == "poisson")
    {
        stop("poisson not complete yet")
        #y <- object$y
        #ind <- y != 0
        #logL <- -object$loss + sum(y[ind] * log(y[ind])) - sum(y) - sum(lfactorial(y))
    } else if (object$family == "coxph")
    {
        logL <- -1e99
    }

    #attr(logL,"df") <- df
    #attr(logL,"nobs") <- n
    #class(logL) <- "logLik"
    logL
}


