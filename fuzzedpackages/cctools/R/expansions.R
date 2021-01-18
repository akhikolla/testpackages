#' Numeric model matrix for continuous convolution
#'
#' Turns ordered variables into integers and expands factors as binary dummy
#' codes. [cont_conv()] additionally adds noise to discrete variables, but this is only
#' useful for estimation. `[cc_prepare()]` can be used to evaluate an already
#' fitted estimate.
#'
#' @param x a vector or data frame with numeric, ordered, or factor columns.
#'
#' @return A numeric matrix containing the expanded variables. It has additional
#'   type `expanded_as_numeric` and `attr(, "i_disc")` cntains the indices of
#'   discrete variables.
#'
#' @examples
#' # dummy data with discrete variables
#' dat <- data.frame(
#'     F1 = factor(rbinom(100, 4, 0.1), 0:4),
#'     Z1 = as.ordered(rbinom(100, 5, 0.5)),
#'     Z2 = as.ordered(rpois(100, 1)),
#'     X1 = rnorm(100),
#'     X2 = rexp(100)
#' )
#'
#' pairs(dat)
#' pairs(expand_as_numeric(dat))  # expanded variables without noise
#' pairs(cont_conv(dat))          # continuously convoluted data
#'
#' @export
expand_as_numeric <- function(x) {
    if (inherits(x, "expanded_as_numeric"))
        return(x)
    if (!inherits(x, "data.frame"))
        x <- as.data.frame(x, stringsAsFactors = FALSE)
    # which variables will be discrete in the output data frame?
    i_disc <- get_i_disc(x)
    # levels and names of expanded variables
    new_names <- expand_names(x)
    new_levels <- expanded_levels(x)

    # ordered -> integer, factors -> dummy coding
    x <- do.call(cbind, lapply(x, cc_prepare_one))
    colnames(x) <- new_names

    # indicate which variables are discrete
    attr(x, "i_disc") <- i_disc
    attr(x, "levels") <- new_levels
    class(x) <- c("numeric", "matrix", "expanded_as_numeric")

    x
}

#' @importFrom stats model.matrix
#' @noRd
cc_prepare_one <- function(x) {
    if (is.numeric(x)) {
        # nothing to do
    } else if (is.ordered(x)) {
        x <- as.numeric(x) - 1
    } else if (is.factor(x)) {
        # expand factors, first column is intercept
        x <- model.matrix(~ x)[, -1, drop = FALSE]
    } else if (is.character(x)) {
        stop("Don't know how to treat character variables; ",
             "use either numeric, ordered, or factor.")
    } else {
        stop("x has unsupported type (", class(x), "); ",
             "use either numeric, ordered, or factor.")
    }
    x
}

get_i_disc <- function(x)
    which(unlist(lapply(x, is_disc)))

is_disc <- function(x) {
    if (is.numeric(x)) {
        return(FALSE)
    } else if (is.ordered(x)) {
        return(TRUE)
    } else if (is.character(x)) {
        stop("Don't know how to treat character variables; ",
             "use either numeric, ordered, or factor.")
    } else {
        return(rep(TRUE, length(levels(x)) - 1))
    }
}

expanded_levels <- function(x) {
    lvls <- lapply(x, function(y) {
        if (is.numeric(y) | is.ordered(y))
            return(list(levels(y)))
        lapply(seq_along(levels(y)), function(l) as.character(0:1))
    })
    do.call(c, lvls)
}


#' Expand a vector like expand_as_numeric
#'
#' Expands each element according to the factor expansions of columns in
#' [expand_as_numeric()].
#'
#' @param y a vector of length 1 or `ncol(x)`.
#' @param x as in [expand_as_numeric()].
#'
#' @return A vector of size `ncol(expand_as_numeric(x))`.
#'
#' @export
expand_vec <- function(y, x) {
    if (length(y) == 1)
        y <- rep(y, ncol(x))
    if (length(y) == ncol(x)) {
        # replicate number of level times y
        y <- sapply(seq_along(y), function(i) {
            if (is.factor(x[, i]) & !is.ordered(x[, i])) {
                rep(y[i], length(levels(x[, i])) - 1)
            } else {
                y[i]
            }
        })
        y <- unlist(y)
    }

    y
}

#' Expands names for expand_as_numeric
#'
#' Expands each element according to the factor expansions of columns in
#' [expand_as_numeric()].
#'
#' @param x as in [expand_as_numeric()].
#'
#' @return A vector of size `ncol(expand_as_numeric(x))`.
#'
#' @export
expand_names <- function(x) {
    nms <- names(x)
    if (length(nms) == 1)
        nms <- rep(nms, ncol(x))
    if (length(nms) == ncol(x)) {
        # replicate number of level times nms
        nms <- sapply(seq_along(nms), function(i) {
            if (is.factor(x[, i]) & !is.ordered(x[, i])) {
                paste(nms[i], levels(x[, i])[-1], sep = ".")
            } else {
                nms[i]
            }
        })
        nms <- unlist(nms)
    }

    nms
}
