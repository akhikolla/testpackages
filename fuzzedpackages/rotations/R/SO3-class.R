#' `SO3` class for storing rotation data as rotation matrices
#'
#' Creates or tests for objects of class "SO3".
#'
#' Construct a single or sample of rotations in 3-dimensions in 3-by-3 matrix
#' form.  Several possible inputs for \code{x} are possible and they are
#' differentiated based on their class and dimension.
#'
#' For \code{x} an n-by-3 matrix or a vector of length 3, the angle-axis
#' representation of rotations is utilized.  More specifically, each rotation
#' matrix can be interpreted as a rotation of some reference frame about the
#' axis \eqn{U} (of unit length) through the angle \eqn{\theta}.  If a single
#' axis (in matrix or vector format) or matrix of axes are provided for
#' \code{x}, then for each axis and angle the matrix is formed through
#' \deqn{R=\exp[\Phi(U\theta)]}{R=exp[\Phi(U\theta)]} where \eqn{U} is replace
#' by \code{x}.  If axes are provided but \code{theta} is not provided then the
#' length of each axis is taken to be the angle of rotation, theta.
#'
#' For \code{x} an n-by-4 matrix of quaternions or an object of class
#' \code{"Q4"}, this function will return the rotation matrix equivalent of
#' \code{x}.  See \code{\link{Q4}} or the vignette "rotations-intro" for more
#' details on quaternions.
#'
#' For \code{x} an n-by-9 matrix, rows are treated as 3-by-3 matrices; rows that
#' don't form matrices in SO(3) are projected into SO(3) and those that are
#' already in SO(3) are returned untouched.  See \code{\link{project.SO3}} for
#' more on projecting arbitrary matrices into SO(3).  A message is printed if
#' any of the rows are not proper rotations.
#'
#' For \code{x} a \code{"data.frame"}, it is translated into a matrix of the
#' same dimension and the dimensionality of \code{x} is used to determine the
#' data type: angle-axis, quaternion or rotation. As demonstrated below,
#' \code{is.SO3} may return \code{TRUE} for a data frame, but the functions
#' defined for objects of class \code{"SO3"} will not be called until
#' \code{as.SO3} has been used.
#'
#' @name SO3
#' @include preliminary.R
#'
#' @param x object to be coerced or tested; see details for possible forms
#' @param theta vector or single rotation angle; if \code{length(theta)==1} the
#'   same theta is used for all axes
#' @param ... additional arguments.
#'
#' @format \code{id.SO3} is the identity rotation given by the the 3-by-3
#'   identity matrix.
#'
#' @return
#'   \item{as.SO3}{coerces provided data into an SO3 type.}
#'   \item{is.SO3}{returns \code{TRUE} or \code{False} depending on whether its
#'     argument satisfies the conditions to be an rotation matrix.  Namely, has
#'     determinant one and its transpose is its inverse.}
#'
#' @examples
#' # Select one location to focus on
#' Loc698 <- subset(nickel, location == 698)
#'
#' is.SO3(Loc698[,5:13])          #Some of the rows are not rotations due to rounding or entry errors
#'                                #as.SO3 will project matrices not in SO(3) to SO(3)
#'
#' Rs <- as.SO3(Loc698[,5:13])    #Translate the Rs data.frame into an object of class 'SO3'
#'                                #Rows 4, 6 and 13 are not in SO(3) so they are projected to SO(3)
#'
#' mean(Rs)                       #Estimate the central orientation with the average
#' median(Rs)                     #Re-estimate central orientation robustly
#' Qs <- as.Q4(Rs)                #Coerse into "SO3" format, see ?as.SO3 for more
#'
#' #Visualize the location, there appears to be two groups
#' \donttest{
#'   plot(Rs, col = c(1, 2, 3))
#' }
NULL

#' @rdname SO3
#' @export
as.SO3 <- function(x, ...) {
  UseMethod("as.SO3")
}

#' @rdname SO3
#' @export
as.SO3.default <- function(x, theta=NULL,...) {
  
  p<-ncol(x)
  n<-nrow(x)
  
  if(is.null(p)){
    p<-length(x)
    n<-1
  }
  
  if(n==3 && p==3 && is.SO3(x) && is.null(theta)){
    
    #If there are 3 rows and columns and the object is already a rotation matrix, the same rotation is returned
    class(x) <- "SO3"
    return(x)
    
  }else if(p==9){
    
    rots<-is.SO3(x)
    notRots<-which(!rots)
    
    if(all(rots)){
      #If there are 9 columns and the data are already rotation matrices then the SO3 class is appeneded and object returned
      
      class(x) <- "SO3"
      return(x)
      
    }else{
      #If there are 9 columns and some are not rotations,
      #those that aren't rotations are projected to SO(3) and the others are left alone
      
      for(i in notRots){
        x[i,]<-project.SO3(x[i,])
      }
      
      txt<-"Row(s)"
      for(i in notRots){
        txt<-paste(txt,i)
      }
      txt<-paste(txt,"was(were) not in SO(3).")
      message(txt)
      class(x)<-"SO3"
      return(x)
      
    }
  }else if(p==3){
    
    #If there are 3 columns, it's assumed the input R is the matrix of unit axes of rotations and the theta vector are the angles,
    #or the length of the axes is the angle of rotation
    
    U<-matrix(x,n,3)
    
    ulen<-sqrt(rowSums(U^2))
    ntheta<-length(theta)
    
    if(is.null(theta)){
      
      theta<-ulen%%(pi)
      
    }else if(ntheta!=n){
      
      if(ntheta==1){
        
        theta<-rep(theta,n)
        
      }else{
        
        stop("Number of angles must match number of axes")
        
      }
      
    }
    
    R<-matrix(NA,n,9)
    
    for(i in 1:n){
      
      if(ulen[i]!=0)
        U[i,]<-U[i,]/ulen[i]
      
      P <- U[i,] %*% t(U[i,])
      
      R[i,] <- P + (diag(3) - P) * cos(theta[i]) + eskew(U[i,]) * sin(theta[i])
    }
    class(R) <- "SO3"
    return(R)
    
  }else if(p==4){
    
    #If there are 4 columns, it's assumed the input is an n-by-4 matrix with rows corresponding to quaternions
    R<-as.Q4(x)
    return(as.SO3(R))
    
  }
  
  stop("Unknown data type.  Please see ?SO3 for more details.")
  
}

#' @rdname SO3
#' @export
as.SO3.Q4 <- function(x, ...) {
  q <- formatQ4(x)
  
  if (any((rowSums(q^2) - 1) > 1.0e-9)) {
    warning("Unit quaternions required. Input was normalized.")
    nonq <- which((rowSums(q^2) - 1) > 1.0e-9)
    q[nonq, ] <- as.Q4(q[nonq, ] / sqrt(rowSums(q[nonq, ]^2)))
  } else
    class(q) <- "Q4"
  
  theta <- mis.angle(q)
  u <- mis.axis(q)
  as.SO3.default(u, theta)
}

#' @rdname SO3
#' @export
as.SO3.SO3 <- function(x, ...) {
  x
}

#' @rdname SO3
#' @export
as.SO3.data.frame <- function(x, ...) {
  n <- nrow(x)
  p <- ncol(x)
  R <- as.matrix(x, n, p)
  as.SO3.default(R)
}

#' @rdname SO3
#' @export
is.SO3 <- function(x) {
  Rlen <- length(x)
  
  if (Rlen %% 9 != 0)
    return(FALSE)
  
  R <- x
  
  if (class(x)[1] == 'data.frame')
    R <- data.matrix(R)
  
  R <- matrix(R, ncol = 9)
  apply(R, 1, function(R) {
    R <- matrix(R, 3, 3)
    if (any(is.na(R)))
      return(FALSE)
    if (abs(det(R) - 1) > 1.0e-9)
      return(FALSE)
    all(abs(t(R) %*% R - diag(1, 3)) < 1.0e-4)
  })
}

#' @rdname SO3
#' @export
id.SO3 <- as.SO3(c(1, 0, 0), 0)