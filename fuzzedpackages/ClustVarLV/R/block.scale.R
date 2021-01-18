#' @title Scaling of a three-way array
#'
#' @description Centering and scaling of a three-way array
#'
#' @param X : a three-way array
#' @param xcenter : centering of X. By default X will be centered for both mode 2 and 3 (xcenter=TRUE), otherwise xcenter=FALSE
#' @param xscale : scaling parameter applied to X. By default no scaling (xscale=0) \cr
#'             0 : no scaling only centering - the default \cr
#'             1 : scaling with standard deviation of  (mode 2 x mode 3) elements \cr
#'             2 : global scaling (each block i.e. each mode 2 slice will have the same inertia ) \cr
#'             3 : global scaling (each block i.e. each mode 3 slice will have the same inertia )
#'
#' @return Xscaled     : the scaled three-way array
#' 
#' @export


block.scale<-function(X,xcenter=TRUE,xscale=0) {
  p <- dim(X)[2] # number of elements  (mode 2)
  n <- dim(X)[1] # number of elements (mode 1) # in sensory ususally products
  q <- dim(X)[3] # number of elements (mode 3) # in QDA assessors
  if (xcenter) {
    meanX<-plyr::aaply(.data=X,.margins=c(2,3),.fun=mean) # centering each descriptor for each assessor
    X<-plyr::aaply(X,1,function(s) {s-meanX})
  }
  if (xscale==1) {
    # each column mode2 x mode 3 will be standardized
    sdX<-plyr::aaply(.data=X,.margins=c(2,3),.fun=sd)
    X<-plyr::aaply(X,1,function(s) {s/sdX})
  }
  else {
    if (xscale==2) {
      # tables associated with mode 2 are set at the same footing that is to say the same inertia
      sst<-function(x) {return(x^2)}
      sstX<-plyr::aaply(.data=X,.margins=c(2,3),.fun=sst)
      localInertia<-apply(sstX,1,sum) #same inertia for all consumers
      globalInertia<-sum(localInertia)

      for (ip in 1:p) {
        X[,ip,]<-X[,ip,]*(sqrt(globalInertia)/sqrt(p*localInertia[ip]))
      }
    }
    else if (xscale==3) {
      # tables associated with mode 3 are set at the same footing that is to say the same inertia
      sst<-function(x) {return(x^2)}
      sstX<-plyr::aaply(.data=X,.margins=c(2,3),.fun=sst)
      localInertia<-apply(sstX,2,sum)
      globalInertia<-sum(localInertia)

      for (iq in 1:q) {
        X[,,iq]<-X[,,iq]*globalInertia/(q*localInertia[iq])
      }

    }
  }
  Xscaled=X
  return(Xscaled)
}
