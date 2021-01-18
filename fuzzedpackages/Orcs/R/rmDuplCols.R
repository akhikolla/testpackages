#' Remove duplicated columns from \code{data.frame}
#' 
#' @description
#' Automatically detect and remove columns from a \code{data.frame}
#' based on duplicated headers.
#' 
#' @param x Input \code{data.frame}.
#' @param keep_first A \code{logical} determining whether the first column of an 
#' otherwise duplicated header should be kept, defaults to \code{TRUE}. 
#' @param ... Currently not in use.
#' 
#' @return 
#' Revised \code{data.frame}.
#' 
#' @author 
#' Florian Detsch
#' 
#' @seealso
#' \code{\link{duplicated}}
#' 
#' @examples
#' ## sample data
#' set.seed(123)
#' dat <- data.frame(matrix(rnorm(28), nc = 7))
#' names(dat) <- c("Col1", "Col1", "Col1", "Col2", "Col3", "Col3", "Col4")
#' 
#' dat
#' rmDuplCols(dat)
#' rmDuplCols(dat, keep_first = FALSE)
#' 
#' @export rmDuplCols
#' @name rmDuplCols
rmDuplCols <- function(x, keep_first = TRUE, ...) {
  
  ## Identify duplicated column names
  dupl.cols.ff <- duplicated(names(x))
  dupl.cols.fl <- duplicated(names(x), fromLast = TRUE)
  
  dupl.cols = if (keep_first) {
    which(dupl.cols.ff)
  } else {
    unique(c(which(dupl.cols.ff), which(dupl.cols.fl)))
  }
  
  # Remove duplicated columns from initial dataframe
  x[, -dupl.cols, drop = FALSE]
}