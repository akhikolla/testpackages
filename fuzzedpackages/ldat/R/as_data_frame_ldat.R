
#' @export
as.data.frame.ldat <- function(x, ...) {
  res <- lapply(x, as_rvec)
  as.data.frame(res, stringsAsFactors = FALSE)
}

