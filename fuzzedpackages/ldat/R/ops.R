
#' Implementation of Ops group generics for lvec
#'
#' @param e1 an object of type \code{\link{lvec}}.
#' @param e2 an object of type \code{\link{lvec}}.
#'
#' @details
#' Math is group generic implementing the following functions: 
#' ‘"+"’, ‘"-"’, ‘"*"’, ‘"/"’, ‘"^"’, ‘"%%"’, ‘"%/%"’ ‘"&"’, 
#' ‘"|"’, ‘"!"’ ‘"=="’, ‘"!="’, ‘"<"’, ‘"<="’, ‘">="’, ‘">"’.
#' For more information see \code{\link{Ops}}. 
#'
#' @return
#' Returns an \code{link{lvec}} of the same length as the input.
#'
#' @import lvec
#' @export
Ops.lvec <- function(e1, e2) {
  # Process input
  if (!missing(e2)) {
    if (!(length(e2) == length(e1)) && !(length(e1) == 1 || length(e2) == 1)) {
      stop(paste0("Vectors do not have the same length and none of them has a", 
        " length of one."))
    }
  }
  # Process the input in chunks
  chunks <- if (length(e1) != 1 || missing(e2)) chunk(e1)  else chunk(e2)
  length <- if (length(e1) != 1 || missing(e2)) length(e1) else length(e2)
  result <- NULL
  for (c in chunks) {
    # Read chunk from first vector
    d1 <- if (length(e1) == 1) as_rvec(e1) else 
      as_rvec(lget(e1, range = c))
    # Call the specific function ('+', '/' etc) for the chunk
    if (!missing(e2)) {
      d2 <- if (length(e2) == 1) as_rvec(e2) else 
        as_rvec(lget(e2, range = c))
      r <- do.call(.Generic, list(d1, d2))
    } else {
      r <- do.call(.Generic, list(d1))
    }
    # Append to result
    r <- as_lvec(r)
    if (is.null(result)) {
      result <- r
      length(result) <- length
    } else {
      if (lvec_type(r) == 'character' && strlen(r) > strlen(result)) {
        warning("Changing maximum string length")
        strlen(result) <- strlen(r)
      }
      lset(result, range = c, r)
    }
  }
  result
}

