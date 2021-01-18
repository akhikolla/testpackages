### predict.basad.R
### The default predict method for the basad
###
###
### Author: Qingyan Xiang

predict.basad <- function(object, testx = NULL, ...)
{
    
    
    #### Check the testData
	if (missing(testx))
        stop("testing data missing ...")
    if( any( is.na(testx)) )
        stop("NA not permitted in x")
    
    if (is.null(colnames(testx)))
    {
    if ( (  ncol(testx) + 1 ) != ncol( object$x )  ) stop("test data dimension does not match training data, variable names are not supplied...")
    }
    
    testx = cbind( rep(1, nrow(testx)), testx )


    newB <- numeric( ncol( object$x ) )
    newB[object$model.index] <- object$est.B[object$model.index]
    y = testx %*% newB
    
    return(y)
}
