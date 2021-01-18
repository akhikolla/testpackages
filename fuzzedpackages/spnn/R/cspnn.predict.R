
cspnn.predict <- function(nn,
                          newData){
  
  newData <- as.matrix(newData, ncol = nn$k)
  
  probs <- .cspnn_predict_cpp(xr = nn$xr,
                              newData = as.matrix(newData),
                              sigmaInverse = nn$sigmaInverse)
  
  colnames(probs) <- nn$categories
  
  categories <- sapply(max.col(m = probs, ties.method = "first"), function(x) nn$categories[x])
  
  results <- list(categories = categories, 
                  probabilities = probs)
  
  return(results)
}

