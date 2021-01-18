#' @include allClasses.R
NULL

#' Method for discretizing a new input dataset given a discretization scheme of class \code{\link{glmdisc}}.
#' @rdname discretize
#' @name discretize
#' @aliases discretize,glmdisc-method
#' @description This defines the method "discretize" which will discretize a new input dataset given a discretization scheme of S4 class \code{\link{glmdisc}}
#' @examples
#' # Simulation of a discretized logit model
#' set.seed(1)
#' x <- matrix(runif(300), nrow = 100, ncol = 3)
#' cuts <- seq(0, 1, length.out = 4)
#' xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
#' theta <- t(matrix(c(0, 0, 0, 2, 2, 2, -2, -2, -2), ncol = 3, nrow = 3))
#' log_odd <- rowSums(t(sapply(seq_along(xd[, 1]), function(row_id) {
#'   sapply(
#'     seq_along(xd[row_id, ]),
#'     function(element) theta[xd[row_id, element], element]
#'   )
#' })))
#' y <- rbinom(100, 1, 1 / (1 + exp(-log_odd)))
#'
#' sem_disc <- glmdisc(x, y,
#'   iter = 50, m_start = 4, test = FALSE,
#'   validation = FALSE, criterion = "aic"
#' )
#' discretize(sem_disc, data.frame(x))
methods::setMethod("discretize", methods::signature(object = "glmdisc"), function(object, data) {
  discretize_link(object@best.disc[[2]], data, object@parameters$m_start)
})
