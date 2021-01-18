#' Create \code{data.frame} from \code{list}
#' 
#' Create a \code{data.frame} from a \code{list} directly, \emph{i.e.} without 
#' being required to explicitly call \code{\link{rbind}} first.
#' 
#' @param x A \code{list} object.
#' @param bind Binding direction. Available options are \code{"rows"} (default) 
#' and \code{"cols"} for \code{\link{rbind}} and \code{\link{cbind}}, 
#' respectively.
#' @param ... Additional arguments passed to \code{\link{data.frame}}.
#' 
#' @return A \code{data.frame} object.
#' 
#' @seealso 
#' \code{\link{data.frame}}, \code{\link{rbind}}, \code{\link{cbind}}.
#' 
#' @examples 
#' lst <- list(letters[1:3], letters[4:6], letters[7:9])
#' 
#' do.call("rbind", lst) # results in matrix
#' list2df(lst)          # results in data.frame created using rbind()
#' list2df(lst, bind = "cols") # same for cbind()
#' 
#' @export list2df
#' @name list2df
list2df <- function(x, bind = c("rows", "cols"), ...) {
  
  bind = bind[1]
  
  mat_x = if (bind == "rows") {
    mat_x = do.call("rbind", x)
  } else if (bind == "cols") {
    do.call("cbind", x)
  } else {
    stop("Binding direction unknown, choose either 'rows' or 'cols'.\n")
  }
  
  data.frame(mat_x, ...)
}