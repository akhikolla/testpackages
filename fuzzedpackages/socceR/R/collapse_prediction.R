#' Create a matrix to collapse tournament predictions to ranks
#'
#' Creates a matrix to collapse the rows of a tournamewnt prediction matrix
#'
#' Returns a vector of numeric values. Elements in the input factor that cannot be converted to numeric will produce NA.
#'
#' @param ranks An integer vector of R ordered elements giving the cut offs of the ranks to create 
#' @return Returns a numeric matrix with R rows and T columns that can be multiplied on a square prediction matrix to obtain the collapsed predictions
#' @author Claus Ekstrom \email{ekstrom@@sund.ku.dk}
#' @keywords manip
#' @examples
#'
#' m2 <- matrix(c(.5, .5, 0, 0, .5, .5, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), 4)
#' # Collapse into ranks 1, 2, and 3+4
#' collapse <- collapse_prediction(c(1, 2, 4))
#'
#' collapsed_prediction <- collapse %*% m2
#' collapsed_prediction
#' 
#' @export
collapse_prediction <- function(ranks=c(1, 2, 3, 4, 8, 16, 32)) {
  m <- matrix(0, ncol=max(ranks), nrow=length(ranks))
  start <- 1
  for (i in 1:length(ranks)) {
    m[i, start:ranks[i]] <- 1
    start <- ranks[i]+1
  }
  m
}