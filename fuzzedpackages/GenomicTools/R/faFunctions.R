seqlength <- function(x){
  out <- nchar(x)
  names(out) <- names(x)
  out
}