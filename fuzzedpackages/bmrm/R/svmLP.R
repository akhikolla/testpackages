

#' Linearly Programmed SVM
#' 
#' Linearly Programmed L1-loss Linear Support Vector Machine with L1 regularization  
#' 
#' svmLP solves a linear program implementing a linear SVM 
#' with L1 regularization and L1 loss. It solves:
#' min_w LAMBDA*|w| + sum_i(e_i);
#' s.t. y_i * <w.x_i> >= 1-e_i; e_i >= 0
#' where |w| is the L1-norm of w
#'
#' svmMulticlassLP solves a linear program implementing multiclass-SVM 
#' with L1 regularization and L1 loss. It solves:
#' min_w LAMBDA*|w| + sum_i(e_i);
#' s.t. <w.x_i> - <w.x_j> >= 1-e_i; e_i >= 0
#' where |w| is the L1-norm of w
#'
#' @name lpSVM
#' @param x a numeric data matrix to predict
#' @param y a response factor for each row of x. It must be a 2 levels factor for svmLP, or a >=2 levels factor for svmMulticlassLP
#' @param loss.weights numeric vector of loss weights to incure for each instance of x. 
#'        Vector length should match length(y), but values are cycled if not of identical size.
#' @param LAMBDA control the regularization strength in the optimization process. This is the value used as coefficient of the regularization term.
#' @param object an object of class svmLP or svmMLP
#' @param ... unused, present to satisfy the generic predict() prototype
#' @return the optimized weights matrix, with class svmLP
#' @author Julien Prados
#' 
#' @examples
#'   x <- cbind(100,data.matrix(iris[1:4]))
#'   y <- iris$Species
#'   w <- svmMulticlassLP(x,y)
#'   table(predict(w,x),y)
#'
#'   w <- svmLP(x,y=="setosa")
#'   table(predict(w,x),y)
NULL




#' @describeIn lpSVM linear programm solving binary-SVM with L1-regularization and L1-norm
#' @export
#' @import lpSolve
svmLP <- function(x,y,LAMBDA=1,loss.weights=1) {
  if (!is.logical(y)) stop("y must be logical")
  if (nrow(x)!=length(y)) stop("length(y) must match nrow(x)")
  loss.weights <- rep_len(loss.weights,nrow(x))
  y.num <- ifelse(y,+1,-1)
  opt <- lp(direction = "min",
            objective.in = c(rep(LAMBDA,2L*ncol(x)),loss.weights),
            const.mat = cbind(y.num*x,-y.num*x,diag(nrow(x))),
            const.dir = ">=",
            const.rhs = 1
  )
  u <- opt$solution[seq(1L,length.out=ncol(x))] 
  v <- opt$solution[seq(ncol(x)+1L,length.out=ncol(x))]
  w <- u-v
  class(w) <- "svmLP"
  w
}



#' @rdname lpSVM
#' @return predict() return predictions for row of x, with an attribute "decision.value"
#' @export
predict.svmLP <- function(object,x,...) {
  f <- x %*% object
  y <- f>0
  attr(y,"decision.values") <- f
  return(y)
}






#' @describeIn lpSVM linear programm solving multiclass-SVM with L1-regularization and L1-norm
#' @export
#' @import lpSolve
svmMulticlassLP <- function(x,y,LAMBDA=1,loss.weights=1) {
  y <- as.factor(y)
  if (nlevels(y)<2) stop("nlevels(y) must be >=2")
  if (nrow(x)!=length(y)) stop("length(y) must match nrow(x)")
  loss.weights <- rep_len(loss.weights,nrow(x))

  opt <- lp(direction = "min",
            objective.in = c(rep(LAMBDA,2L*ncol(x)*nlevels(y)),loss.weights),
            const.mat = local({
              L <- merge(data.frame(i=seq_along(y),y=y),data.frame(Y=factor(levels(y),levels(y))))
              L <- L[L$y != L$Y,]
              
              m <- matrix(rep(seq_len(nlevels(y)),each=ncol(x)),nlevels(y),ncol(x)*nlevels(y),byrow=TRUE)
              m <- (m == row(m)) + 0
              m <- m[L$y,] - m[L$Y,]
              m <- m * as.vector(x[L$i,])
              
              n <- matrix(0,nrow(L),nrow(x))
              n[cbind(seq_len(nrow(L)),L$i)] <- 1
              
              cbind(m,-m,n)
            }),
            const.dir = ">=",
            const.rhs = 1
  )
  u <- matrix(opt$solution[seq(1L,length.out=ncol(x)*nlevels(y))],ncol(x)) 
  v <- matrix(opt$solution[seq(ncol(x)*nlevels(y)+1L,length.out=ncol(x)*nlevels(y))],ncol(x))
  w <- u - v
  colnames(w) <- levels(y)
  class(w) <- "svmMLP"
  w
}

    


#' @rdname lpSVM
#' @return predict() return predictions for row of x, with an attribute "decision.value"
#' @export
predict.svmMLP <- function(object,x,...) {
  f <- x %*% object
  y <- colnames(f)[max.col(f,ties.method = "first")]
  attr(y,"decision.values") <- f
  return(y)
}













