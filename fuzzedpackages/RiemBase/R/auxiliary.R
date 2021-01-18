##  01. aux_stack3d  : from 'riemdata' instance, create a 3d slice of data
##  02. aux_rndivide : random division of integers



# 01. aux_stack3d : from 'riemdata' instance, create a 3d slice of --------
#' @keywords internal
#' @noRd
aux_stack3d <- function(riemdata){
  msize = riemdata$size
  ndata = length(riemdata$data)
  
  matdata = array(0,c(msize[1], msize[2], ndata))
  for (i in 1:ndata){
    matdata[,,i] = (riemdata$data[[i]])
  }
  return(matdata)
}


# 02. aux_rndivide --------------------------------------------------------
#' @keywords internal
#' @noRd
aux_rndivide <- function(n, k){
  kk = as.integer(k)
  randomDraw = rnorm(n) # number of data
  kQuantiles = quantile(randomDraw, 0:kk/kk)
  whichK <- cut(randomDraw, kQuantiles, include.lowest = TRUE)  # Divide randomDraw into kk equally-sized groups
  levels(whichK) <- 1:kk
  
  output = list()
  for (i in 1:kk){
    output[[i]] = which(whichK==i)
  }
  return(output)
}
