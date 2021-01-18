#' Unlist the outcome of \code{strsplit}
#' 
#' @description
#' Per default, \code{\link{strsplit}} returns a \code{list}, with each entry 
#' holding the vector of splits of the initial string(s). This function is a 
#' simple wrapper that casts \code{\link{unlist}} upon the returned list to 
#' produce a concatenated \code{character} vector consisting of the single split 
#' elements. 
#' 
#' @param x A \code{character} vector with elements to be split. 
#' @param split A \code{character} vector used for splitting, see 
#' \code{\link{strsplit}}. 
#' @param ... Additional arguments passed to \code{\link{strsplit}}.
#' 
#' @author 
#' Florian Detsch
#' 
#' @seealso
#' \code{\link{strsplit}}
#' 
#' @examples
#' ## 1st example
#' x <- "This is a test."
#' unlistStrsplit(x, " ")
#' 
#' ## 2nd example; note that 'split' defaults to 'whitespace'
#' x2 <- "This is a 2nd test."
#' unlistStrsplit(c(x, x2))
#' 
#' @export unlistStrsplit
#' @name unlistStrsplit
unlistStrsplit <- function(x, split, ...) {
  
  if (missing(split))
    split <- " "
  
  ls_split <- strsplit(x, split = split, ...)
  ch_split <- unlist(ls_split)
  
  return(ch_split)
}