make_order <- function(exprs, times, order) {
  stopifnot(is.list(exprs))
  if (order == "block") {
    res <- rep.int(seq_along(exprs), rep.int(times, length(exprs)))
  } else {
    res <- rep.int(seq_along(exprs), times)
    if (order == "random") res <- sample(res)
  }
  structure(res, class = "factor", levels = names(exprs))
}
