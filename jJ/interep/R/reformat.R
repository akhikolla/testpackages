#' This function changes the format of the longitudinal data from wide format to long format
#'
#' @param k the number of repeated measurement.
#' @param y the longitudinal response.
#' @param x a matrix of predictors, consisting of omics and environment factors, as well as their interactions. In
#' the case study, the omics measurements are lipidomics data.
#' @export
reformat <- function(k,y,x){
  n=dim(y)[1]
  response=y
  id=rep(0,n*k)
  y=rep(0,n*k)
  for (i in 1:n) {
    for (j in 1:k) {
      id[(i-1)*k+j]=i
      y[(i-1)*k+j]=response[i,j]
    }
  }

  data=cbind(id=id,y,x[rep(1:nrow(x), times = rep(k,n)), ])

  data=as.data.frame(data)

  y=data[,2]
  x=data[,-c(1,2)]
  x=cbind(data.frame(rep(1,length(data$id))),x)
  x=data.matrix(x)
  return(list("y"=y,"x"=x,"id"=id))
}
