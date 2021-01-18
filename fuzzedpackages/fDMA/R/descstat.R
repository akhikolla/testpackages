
descstat <- function(data)
  {
    ### data - a matrix or xts object 

    if (missing(data)) { stop("please, specify data") }
    if (! is.matrix(data)) 
      { 
        data <- as.matrix(data)
        warning("data should be a matrix, the function tried to convert data to a matrix") 
      }

    x <- describe(data)[,c(3:5,8,9,11,12)]
    x <- cbind(x,(x[,2])^2)
    x <- x[,c(1,2,8,3:7)]
    colnames(x)[3] <- "variance"
    x <- cbind(x,x[,1]/x[,2])
    colnames(x)[9] <- "coeff. of variation"
    
    x <- format(x,scientific=FALSE,digits=4)
    
    return(data.matrix(x))
  }
  