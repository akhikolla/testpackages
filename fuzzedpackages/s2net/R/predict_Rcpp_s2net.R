# predict.Rcpp_s2net S4 method

predict_Rcpp_s2net <- function(object, newX, type = "default"){
  if(class(newX)=="s2Data"){
    # No need to pre-process
    newX = newX$xL
  }
  if(!is.matrix(newX)){
    stop("newX shoud be a matrix or s2Data containing xL")
  }
  switch (type,
    response = {code = 1},
    probs = {code = 2},
    class = {code = 3},
    {code = 0}
  )
  return(object$predict(newX, code))
}
