
spnn.predict <- function(nn, 
                         newData){
  
  setCube <- array(NaN, # intialize input cube
                   c(nrow(nn$set),
                     ncol(nn$set) - 1,
                     length(nn$categories)))
  
  for(i in 1:length(nn$categories)){ # split set into an array where each slice is a class
    subset <- nn$set[which(nn$set[,nn$category.column] == nn$categories[i]), -nn$category.column]
    setCube[1:nrow(subset), 1:ncol(subset) ,i] <- as.matrix(subset)
  }
  
  probs <- .spnn_predict_cpp(setCube = setCube,
                             newData = as.matrix(newData),
                             sigmaInverse = nn$sigmaInverse)
  
  colnames(probs) <- nn$categories
  
  categories <- sapply(max.col(m = probs, ties.method = "first"), function(x) nn$categories[x])
  
  results <- list(categories = categories, 
                  probabilities = probs)
  
  return(results)
}
