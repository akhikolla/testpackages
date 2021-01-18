
spnn.learn <- function(set, # training data set
                       nn, # trained probabilistic neural network (optional)
                       sigma, # input covariance matrix (optional)
                       category.column = 1){
  
  if(missing(set)){ stop("input set is missing") }
  if(missing(nn)){ nn <- .spnn.create() }

  if(is.null(nn$set)){
    nn$category.column <- category.column
    nn$set <- set
  }else{
    nn$set <- rbind(nn$set, set)
  }
  
  if(missing(sigma)){ nn$sigma <- cov(nn$set[,-nn$category.column]) }
  
  nn$categories <- levels(factor(nn$set[,nn$category.column]))
  nn$sigmaInverse <- MASS::ginv(nn$sigma)
  nn$k <- length(nn$set[1,]) - 1
  nn$n <- length(nn$set[,1])
  
  return(nn)
}
