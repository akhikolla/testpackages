#' Get Yacas variables
#' 
#' @return Vector of variables defined in yacas
#' 
#' @export
y_ls <- function() {
  x <- yacas("Variables()", addSemi = FALSE, retclass = "character")
  z <- x$LinAlgForm
  z <- setdiff(z, "i") # i = imaginary unit: I in Yacas, i in R
  return(z)
}