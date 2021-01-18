
#' @rdname append
#' @export
append <- function(x, y, ...) {
  UseMethod("append")
}

#' Append a vector to an lvec
#' 
#' @param x \code{\link{lvec}} to append to.
#' @param y vector to append to \code{x}. Is converted to \code{\link{lvec}}
#'   using \code{\link{as_lvec}}.
#' @param clone should \code{x} be cloned first. If not, the input \code{x} is
#'   modified. 
#' @param ... ignored; used to pass additional arguments to other methods.
#'
#' @return
#' Returns an lvec combining both \code{x} and \code{y}. When \code{x} is 
#' \code{NULL} a clone of \code{y} is returned.
#'
#' @rdname append
#' @export
append.lvec <- function(x, y, clone = TRUE, ...) {
  if (!is_lvec(y)) y <- as_lvec(y)
  if (is.null(x) || length(x) == 0) return(clone(y))
  if (!is_lvec(x)) stop("x should be of type lvec (or NULL)")
  lx <- length(x)
  if (clone) x <- clone(x)
  if (lvec_type(x) == "character" && lvec_type(x) == "character" &&
    strlen(x) < strlen(y)) strlen(x) <- strlen(y)
  if (!is.null(rattr(x, "class")) && rattr(x, "class") == "factor") {
    if (is.null(rattr(y, "class")) || rattr(y, "class") != "factor")
      stop("y is not a factor vector")
    if (!isTRUE(all.equal(rattr(x, "levels"), rattr(y, "levels")))) {
      newlevels <- union(rattr(x, "levels"), rattr(y, "levels"))
      rattr(x, "levels") <- newlevels
      y <- elementwise(y, function(d) factor(d, levels = newlevels))
    }
  }
  length(x) <- lx + length(y)
  lset(x, range = c(lx+1, length(x)), values = y)
}

#' @rdname append
#' @export
append.ldat <- function(x, y, clone = TRUE, ...) {
  if (!is_ldat(x) && !is.null(x)) stop("x should be of type ldat (or NULL)")
  if (is.null(x)) return(clone(y))
  stopifnot(length(x) == length(y))
  for (i in seq_along(x)) {
    x[[i]] <- append(x[[i]], y[[i]], clone = clone)
  }
  x
}

