#' @references
#' N Benjamin Erichson, Sergey Voronin, Steven L Brunton, and J Nathan Kutz. Randomized matrix decompositions using r. arXiv preprint arXiv:1608.02148, 2016.
#' @useDynLib PINSPlus
#' @import RcppParallel 
#' @importFrom Rcpp evalCpp
rpca.para <- function(A, k=NULL, center=TRUE, scale=TRUE, retx=TRUE, p=10, q=2, rand = TRUE) {
  
  A <- as.matrix(A)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Checks
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (any(is.na(A))) {
    warning("Missing values are omitted: na.omit(A).")
    A <- stats::na.omit(A)
  }   
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Init rpca object
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rpcaObj = list(rotation = NULL,
                 eigvals = NULL,
                 sdev = NULL,
                 var = NULL,
                 center = center,
                 scale = scale,
                 x=NULL)
  
  m <- nrow(A)
  n <- ncol(A)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Set target rank
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(is.null(k)) rand <- FALSE
  if(is.null(k)) k <- min(n,m)
  if(k > min(n,m)) k <- min(n,m)
  if(k<1) stop("Target rank is not valid!")
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Center/Scale data
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(center == TRUE) {
    rpcaObj$center <- colMeans(A)
    A <- sweep(A, MARGIN = 2, STATS = rpcaObj$center, FUN = "-", check.margin = TRUE)
    #A <- H(H(A) - rpcaObj$center)
  } else { rpcaObj$center <- FALSE }
  
  if(scale == TRUE) {
    rpcaObj$scale <- sqrt(colSums(A**2) / (m-1))
    if(is.complex(rpcaObj$scale)) { rpcaObj$scale[Re(rpcaObj$scale) < 1e-8 ] <- 1+0i  
    } else {rpcaObj$scale[rpcaObj$scale < 1e-8] <- 1}
    A <- sweep(A, MARGIN = 2, STATS = rpcaObj$scale, FUN = "/", check.margin = TRUE)
    #A <- H(H(A) / rpcaObj$scale)
    #A[is.nan(A)] <- 0
  } else { rpcaObj$scale <- FALSE }
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Compute randomized svd / eigen decomposition
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(rand == TRUE) {
    svdalg = 'rsvd'
  }else { 
    svdalg = 'svd' 
  }
  
  out <- switch(svdalg,
                svd = svd(A, nu = k, nv = k),
                rsvd = rsvd.para(A, k = k, p = p, q = q),
                stop("Selected SVD algorithm is not supported!")
  )
  
  rpcaObj$eigvals <- switch(svdalg,
                            svd = out$d[1:k]**2 / (m-1),
                            rsvd = out$d**2 / (m-1)
  )
  
  rpcaObj$rotation <- switch(svdalg,
                             svd = out$v,
                             rsvd = out$v
  )
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Explained variance and explained variance ratio
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rpcaObj$sdev <-  sqrt( rpcaObj$eigvals )
  rpcaObj$var <- sum( matrixStats::colVars( Re(A) ) )
  if(is.complex(A)) rpcaObj$var <- Re(rpcaObj$var + sum( apply( Im(A) , 2, stats::var ) ))
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Add row and col names
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rownames(rpcaObj$rotation) <- colnames(A)
  colnames(rpcaObj$rotation) <- paste(rep('PC', k), 1:k, sep = "")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Compute rotated data
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(retx==TRUE) {
    rpcaObj$x <- mat_mul_para(A, rpcaObj$rotation) # slow
    #rpcaObj$x <- H(H(out$u[,1:k]) * out$d[1:k])
    #rpcaObj$x <- sweep(out$u[, 1:k, drop=FALSE], MARGIN = 2, STATS = out$d[1:k], FUN = "*", check.margin = TRUE)
    rownames(rpcaObj$x) <- rownames(A)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Return
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  class(rpcaObj) <- "rpca"
  return( rpcaObj )
  
}#End rPCA



predict.rpca.para <- function( object, newdata, ...)
{
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Predict
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(!is.logical(object$center)) {
    #newdata <- H(H(newdata) - object$center)
    newdata <- sweep(newdata, MARGIN = 2, STATS = object$center, FUN = "-", check.margin = TRUE)
  }
  
  if(!is.logical(object$scale)) {
    #newdata <- H(H(newdata) / object$scale)
    newdata <- sweep(newdata, MARGIN = 2, STATS = object$scale, FUN = "/", check.margin = TRUE)
    newdata[is.nan(newdata)] <- 0
  }
  
  x <- mat_mul_para(as.matrix(newdata), as.matrix(object$rotation) )
  rownames(x) <- rownames(newdata)
  
  return( x )
}


rsvd.para <- function(A, k=NULL, nu=NULL, nv=NULL, p=10, q=2, sdist="normal") {
  #*************************************************************************
  #***        Author: N. Benjamin Erichson <nbe@st-andrews.ac.uk>        ***
  #***                              <2015>                               ***
  #***                       License: BSD 3 clause                       ***
  #*************************************************************************
  
  #Dim of input matrix
  m <- nrow(A)
  n <- ncol(A)
  
  #Flipp matrix, if wide
  if(m < n){
    A <- H(A)
    m <- nrow(A)
    n <- ncol(A)
    flipped <- TRUE
  } else flipped <- FALSE
  
  #Set target rank
  if(is.null(k)) k = n
  if(k > n) k <- n
  if(is.character(k)) stop("Target rank is not valid!")
  if(k < 1) stop("Target rank is not valid!")
  
  #Set oversampling parameter
  l <- round(k) + round(p)
  if(l > n) l <- n
  if(l < 1) stop("Target rank is not valid!")
  
  #Check if array is real or complex
  if(is.complex(A)) {
    isreal <- FALSE
  } else {
    isreal <- TRUE
  }
  
  #Set number of singular vectors
  if(is.null(nu)) nu <- k
  if(is.null(nv)) nv <- k
  if(nu < 0) nu <- 0
  if(nv < 0) nv <- 0
  if(nu > k) nu <- k
  if(nv > k) nv <- k
  if(flipped==TRUE) {
    temp <- nu
    nu <- nv
    nv <- temp
  }
  
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Generate a random sampling matrix O
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  O <- switch(sdist,
              normal = matrix(stats::rnorm(l*n), n, l),
              unif = matrix(stats::runif(l*n), n, l),
              rademacher = matrix(sample(c(-1,1), (l*n), replace = TRUE, prob = c(0.5,0.5)), n, l),
              stop("Selected sampling distribution is not supported!"))
  
  if(isreal==FALSE) {
    O <- O + switch(sdist,
                    normal = 1i * matrix(stats::rnorm(l*n), n, l),
                    unif = 1i * matrix(stats::runif(l*n), n, l),
                    rademacher = 1i * matrix(sample(c(-1,1), (l*n), replace = TRUE, prob = c(0.5,0.5)), n, l),
                    stop("Selected sampling distribution is not supported!"))
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Build sample matrix Y : Y = A * O
  #Note: Y should approximate the range of A
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Y <- mat_mul_para(A, O)
remove(O)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Orthogonalize Y using economic QR decomposition: Y=QR
#If q > 0 perfrom q subspace iterations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if( q > 0 ) {
  for( i in 1:q) {
    Y <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
    Z <- crossprod_help(A , Y )
    Z <- qr.Q( qr(Z, complete = FALSE) , complete = FALSE )
    Y <- mat_mul_para(A, Z)
  }#End for
  remove(Z)
}#End if

Q <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
remove(Y)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Project the data matrix a into a lower dimensional subspace
#B := Q.T * A
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
B <- crossprod_help(Q , A )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Singular Value Decomposition
#Note: B =: U * S * Vt
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rsvdObj <- svd(B, nu=nu, nv=nv) # Compute SVD
rsvdObj$d <- rsvdObj$d[1:k] # Truncate singular values

if(nu != 0) rsvdObj$u <- mat_mul_para(Q, rsvdObj$u) # Recover left singular vectors

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Flip SVD back
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(nu == 0){ rsvdObj$u <- NULL}
if(nv == 0){ rsvdObj$v <- NULL}
if(flipped == TRUE) {
  u_temp <- rsvdObj$u
  rsvdObj$u <- rsvdObj$v
  rsvdObj$v <- u_temp
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Return
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class(rsvdObj) <- "rsvd"
return(rsvdObj) 

} # End rsvd

crossprod_help <- function( A , B ) {
  return( mat_mul_para( c_transpose(A) , B ) )
}

H <- function( X ) {
  return( c_transpose(X) )
  
}