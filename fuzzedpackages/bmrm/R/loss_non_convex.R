



#' @describeIn binaryClassificationLoss Hinge Loss for Linear Support Vector Machine (SVM)
#' @export
hingeLoss <- function(x,y,loss.weights=1) {
  if (!is.logical(y)) stop("y must be logical")
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  loss.weights <- rep(loss.weights,length.out=length(y))
  loss.weights <- loss.weights/sum(loss.weights)
  
  f <- function(w) {
    w <- cbind(matrix(numeric(),ncol(x),0),w)
    f <- x %*% w
    y[is.na(y)] <- f[is.na(y)]>0
    loss <- loss.weights * pmax(f*ifelse(y,-1,+1)+1,0)
    grad <- loss.weights * (loss>0) * ifelse(y,-1,+1)
    lvalue(w) <- colSums(loss)
    gradient(w) <- crossprod(x,grad)
    class(w) <- c("hingeLoss","binaryClassificationLoss")
    return(w)
  }
  is.convex(f) <- all(!is.na(y))
  return(f)
}



#' Soft Margin Vector Loss function for multiclass SVM
#' 
#' @param x instance matrix, where x(t,) defines the features of instance t
#' @param y target vector where y(t) is an integer encoding target of x(t,). If it contains NAs, the return function is
#'        a non-convex loss for transductive multiclass-SVM.
#' @param l loss matrix. l(t,p(t)) must be the loss for predicting target p(t) instead of y(t) 
#'        for instance t. By default, the parameter is set to character value "0/1" so that the loss is set to a 0/1 loss matrix.
#' @return a function taking one argument w and computing the loss value and the gradient at point w
#' @export
#' @references Teo et al.
#'   A Scalable Modular Convex Solver for Regularized Risk Minimization.
#'   KDD 2007
#' @examples
#'   # -- Build a 2D dataset from iris, and add an intercept
#'   x <- cbind(intercept=100,data.matrix(iris[c(1,2)]))
#'   y <- iris$Species
#'   
#'   # -- build the multiclass SVM model
#'   w <- nrbm(softMarginVectorLoss(x,y))
#'   table(predict(w,x),y)
#'   
#'   # -- Plot the dataset, the decision boundaries, the convergence curve, and the predictions
#'   gx <- seq(min(x[,2]),max(x[,2]),length=200) # positions of the probes on x-axis
#'   gy <- seq(min(x[,3]),max(x[,3]),length=200) # positions of the probes on y-axis
#'   Y <- outer(gx,gy,function(a,b) {predict(w,cbind(100,a,b))})
#'   image(gx,gy,unclass(Y),asp=1,main="dataset & decision boundaries",
#'         xlab=colnames(x)[2],ylab=colnames(x)[3])
#'   points(x[,-1],pch=19+as.integer(y))
softMarginVectorLoss <- function(x,y,l=1 - table(seq_along(y),y)) {
  if (!is.matrix(x)) stop('x must be a numeric matrix')
  if (!is.factor(y)) stop('y must be a factor')
  if (nrow(x) != length(y)) stop('dimensions of x and y mismatch')
  if (!identical(nrow(x),nrow(l))) stop('dimensions of x and l mismatch')
  if (any(levels(y)!=colnames(l))) stop('colnames(l) must match with levels(y)')
  
  f <- function(w) {
    W <- matrix(w,ncol(x),ncol(l),dimnames=list(colnames(x),levels(y)))
    fp <- x %*% W
    y[is.na(y)] <- levels(y)[max.col(fp[is.na(y),],ties.method="first")]
    fy <- rowSums(x * t(W[,y]))
    lp <- fp - fy + l
    p <- max.col(lp,ties.method='first')
    lp <- lp[cbind(1:length(p),p)]
    
    # compute gradient
    gy <- gp <- matrix(0,length(y),ncol(W))
    gp[cbind(seq_along(y),p)] <- 1
    gy[cbind(seq_along(y),y)] <- 1
    grad <- gp - gy
    
    w <- as.vector(W)
    attr(w,"model.dim") <- dim(W)
    attr(w,"model.dimnames") <- dimnames(W)
    lvalue(w) <- sum(lp)
    gradient(w) <- as.vector(crossprod(x,grad))
    class(w) <- "softMarginVectorLoss"
    return(w)
  }
  is.convex(f) <- all(!is.na(y))
  return(f)
}





#' @export
predict.softMarginVectorLoss <- function(object,x,...) {
  W <- array(object,attr(object,"model.dim"),attr(object,"model.dimnames"))
  f <- x %*% W
  y <- max.col(f,ties.method="first")
  y <- factor(colnames(W)[y],colnames(W))
  attr(y,"decision.value") <- f
  y
}
