
hit.ratio <- function (y,y.hat,d=NULL)
{

### computes hit ratio for forecast, i.e., 
### checks in which cases the direction of a change given by forecast
### agrees with the change in real data

### y - forecasted time-series,
###     numeric, vector, or one row or one column matrix or xts object

### y.hat - forecast prediction
###     numeric, vector, or one row or one column matrix or xts object


### d - logical, d = FALSE for level time-series,
###     d = TRUE if time-series represent changes,
###     by default d = FALSE

if (missing(y)) { stop("please, specify y") }
if (missing(y.hat)) { stop("please, specify y.hat") }
if (is.xts(y)) { y <- as.matrix(y) } 
if (is.xts(y.hat)) { y.hat <- as.matrix(y.hat) } 
if (! (is.vector(y) || is.numeric(y) || is.matrix(y))) { stop("y must be numeric, vector or matrix") }
if (! (is.vector(y.hat) || is.numeric(y.hat) || is.matrix(y.hat))) { stop("y must be numeric, vector or matrix") }
if (is.matrix(y) && ! ((ncol(y) == 1) || nrow(y) == 1)) { stop("y must be a one column or one row matrix") }
if (is.matrix(y.hat) && ! ((ncol(y.hat) == 1) || nrow(y.hat) == 1)) { stop("y.hat must be a one column or one row matrix") }
if (!(length(y) == length(y.hat))) { stop("y and y.hat must have the same length") }
if (is.null(d)) { d <- FALSE }
if (! is.logical(d)) { stop("d must be logical, i.e., TRUE or FALSE") }

y <- as.vector(y)
y.hat <- as.vector(y.hat)
    
if (d==FALSE)
  {
    test <- length(which(sign((y.hat - (c(NA,y))[1:length(y)])[-1])==sign((c(NA,diff(y)))[-1]))) / (length(y)-1)
  }
else
  {
    test <- length(which(sign(y)==sign(y.hat))) / length(y)
  }

return(round(as.numeric(test),digits=4))

}
