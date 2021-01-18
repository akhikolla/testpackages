tensorCentering <-
function(x, location = NULL){
  r <- length(dim(x)) - 1
  if(is.null(location)){
    location <- apply(x, 1:r, mean)
  }
  x <- sweep(x, 1:r, location, '-')
  attr(x, "location") <- location
  return(x)
}
