


#' Convex and non-convex risk minimization with L2 regularization and limited memory
#' 
#' Use algorithm of Do and Artieres, JMLR 2012 to find w minimizing: 
#' f(w) = 0.5*LAMBDA*l2norm(w) + riskFun(w)
#' where riskFun is either a convex or a non-convex risk function.
#' @name nrbm
#' @param riskFun the risk function to use in the optimization (e.g.: hingeLoss, softMarginVectorLoss). 
#'   The function must evaluate the loss value and its gradient for a given point vector (w).
#'   The function must return the given point vector w, with attributes "lvalue" and "gradient" set.
#' @param LAMBDA control the regularization strength in the optimization process. 
#'   This is the value used as coefficient of the regularization term.
#' @param MAX_ITER the maximum number of iteration to perform. 
#'   The function stop with a warning message if the number of iteration exceed this value
#' @param EPSILON_TOL a numeric value between 0 and 1 controling stoping criteria: 
#'   the optimization end when the ratio between the optimization gap and the
#'   objective value is below this threshold
#' @param w0 initial weight vector where optimization start
#' @param maxCP mximal number of cutting plane to use to limit memory footprint
#' @param convexRisk a length 1 logical telling if the risk function riskFun is convex. 
#'    If TRUE, use CRBM algorithm; if FALSE use NRBM algorithm from Do and Artieres, JMLR 2012
#' @param LowRankQP.method a single character value defining the method used by LowRankQP (should be either "LU" or "CHOL")
#' @param line.search a logical, when TRUE use line search to speed up convergence
#' @return the optimal weight vector (w)
#' @references Do and Artieres
#'   Regularized Bundle Methods for Convex and Non-Convex Risks
#'   JMLR 2012
#' @export
#' @import LowRankQP
#' @importFrom utils head
#' @examples
#'   # -- Create a 2D dataset with the first 2 features of iris, with binary labels
#'   x <- data.matrix(iris[1:2])
#'
#'   # -- Add a constant dimension to the dataset to learn the intercept
#'   x <- cbind(intercept=1000,x)
#'
#'   # -- train scalar prediction models with maxMarginLoss and fbetaLoss
#'   models <- list(
#'     svm_L1 = nrbmL1(hingeLoss(x,iris$Species=="setosa"),LAMBDA=1),
#'     svm_L2 = nrbm(hingeLoss(x,iris$Species=="setosa"),LAMBDA=1),
#'     f1_L1 = nrbmL1(fbetaLoss(x,iris$Species=="setosa"),LAMBDA=1),
#'     tsvm_L2 = nrbm(hingeLoss(x,
#'                    ifelse(iris$Species=="versicolor",NA,iris$Species=="setosa")),
#'                    LAMBDA=1)
#'   )
#'
#'   # -- Plot the dataset and the predictions
#'   plot(x[,-1],pch=ifelse(iris$Species=="setosa",1,2),main="dataset & hyperplanes")
#'   legend('bottomright',legend=names(models),col=seq_along(models),lty=1,cex=0.75,lwd=3)
#'   for(i in seq_along(models)) {
#'     w <- models[[i]]
#'     if (w[3]!=0) abline(-w[1]*1000/w[3],-w[2]/w[3],col=i,lwd=3)
#'   }
#'
#'
#'   # -- fit a least absolute deviation linear model on a synthetic dataset
#'   # -- containing 196 meaningful features and 4 noisy features. Then
#'   # -- check if the model has detected the noise
#'   set.seed(123)
#'   X <- matrix(rnorm(4000*200), 4000, 200)
#'   beta <- c(rep(1,ncol(X)-4),0,0,0,0)
#'   Y <- X%*%beta + rnorm(nrow(X))
#'   w <- nrbm(ladRegressionLoss(X/100,Y/100),maxCP=50)
#'   barplot(as.vector(w))
NULL


#' @describeIn nrbm original L2-regularized version of nrbm
#' @export
nrbm <- function(riskFun,LAMBDA=1,MAX_ITER=1000L,EPSILON_TOL=0.01,w0=0,maxCP=50L,convexRisk=is.convex(riskFun),LowRankQP.method="LU",line.search=!convexRisk) {
  # check parameters
  if (maxCP<3) stop("maxCP should be >=3")
  cat(ifelse(convexRisk%in%TRUE,"Run nrbm with convex loss\n","Run nrbm with non-convex loss\n"))
  
  # intialize first point estimation
  neval <- 1
  w <- riskFun(w0)
  at <- as.vector(gradient(w))
  bt <- as.vector(lvalue(w)) - crossprod(as.vector(w),at)
  f <- as.vector(LAMBDA*0.5*crossprod(as.vector(w)) + lvalue(w))
  st <- 0
  
  # initialize aggregated working plane and working set
  A <- rbind(at,at);b <- c(bt,bt);s <- c(st,st)
  inactivity.score <- c(NA_real_,NA_real_)
  ub <- f;ub.w <- w;ub.t <- 2
  for (i in 1:MAX_ITER) {
    
     # optimize underestimator
     H <- matrix(0,1L+nrow(A),1L+nrow(A))
     H[-1,-1] <- tcrossprod(A)
     alpha <- LowRankQP(H,c(0,-LAMBDA*b),matrix(1,1L,nrow(H)),1,rep(1,nrow(H)),method=LowRankQP.method)$alpha[-1L]
     
     # compute the optimum point and corresponding objective value
     w <- as.vector(-crossprod(A,alpha) / LAMBDA)
     lb <- as.vector(LAMBDA*0.5*crossprod(as.vector(w)) + max(0,A %*% as.vector(w) + b))
     
     if (line.search) {
       w <- wolfe.linesearch(riskFun,x0=ub.w,s0=w-as.vector(ub.w),f.adjust=function(w) {
         lvalue(w) <- lvalue(w) + as.vector(LAMBDA*0.5*crossprod(w))
         gradient(w) <- as.vector(gradient(w)) + LAMBDA*as.vector(w)
         w
       })
       neval <- neval + attr(w,"neval")
     } else {
       neval <- neval + 1
     }
     
     # estimate loss at the new underestimator optimum
     w <- riskFun(w)
     f <- as.vector(LAMBDA*0.5*crossprod(as.vector(w)) + lvalue(w))

    # update inactivity score
    inactivity.score <- inactivity.score + pmax(1-alpha,0)

    # deduce parameters of the new cutting plane
    at <- as.vector(gradient(w))
    bt <- lvalue(w) - crossprod(as.vector(w),at)

    if (!(convexRisk %in% TRUE)) {
      # solve possible conflicts with the new cutting plane
      if (f<ub) { # descent step
        st <- 0
        s <- s + as.vector(0.5*LAMBDA*crossprod(as.vector(ub.w)-as.vector(w)))
        b <- pmin(b,lvalue(w) - (A %*% as.vector(w)) - s)
      } else { # null step
        st <- 0.5*LAMBDA*crossprod(as.vector(w)-as.vector(ub.w))
        if (lvalue(ub.w) < st + crossprod(at,as.vector(ub.w)) + bt) {
          U <- lvalue(ub.w) - crossprod(at,as.vector(ub.w)) - st
          L <- ub - crossprod(at,as.vector(w)) - 0.5*LAMBDA*crossprod(as.vector(w))
          if (L<=U) {
            bt <- L
          } else {
            at <- -LAMBDA*as.vector(ub.w)
            bt <- ub - 0.5*LAMBDA*crossprod(as.vector(w)) - crossprod(at,as.vector(w))
          }
        }
      }
    }

    # update aggregated cutting plane
    A[1,] <- alpha %*% A
    b[1] <- alpha %*% b
    
    # add the new cutting plane to the working set
    if (nrow(A)<maxCP) {
      A <- rbind(A,at);b <- c(b,bt);s <- c(s,st)
      inactivity.score <- c(inactivity.score,0)
      t <- length(b)
    } else {
      t <- which.max(inactivity.score)
      A[t,] <- at;b[t] <- bt;s[t] <- st
      inactivity.score[t] <- 0
    }
    
    # if new best found
    if (f<ub) {
      inactivity.score[ub.t] <- 0
      ub <- f;ub.w <- w;ub.t <- t
      inactivity.score[ub.t] <- NA
    }
    
    # test end of the convergence
    cat(sprintf("%d:ncall=%d gap=%g obj=%g reg=%g risk=%g w=[%g,%g]\n",i,neval,ub-lb,ub,LAMBDA*0.5*crossprod(as.vector(ub.w)),lvalue(ub.w),min(ub.w),max(ub.w)))
    if (ub-lb < ub*EPSILON_TOL) break
  }
  if (i >= MAX_ITER) warning('max # of itertion exceeded')
  attr(ub.w,"obj") <- ub
  return(ub.w)
}







#' @describeIn nrbm L1-regularized version of nrbm that can only handle convex risk
#' @export
nrbmL1 <- function(riskFun,LAMBDA=1,MAX_ITER=300L,EPSILON_TOL=0.01,w0=0,maxCP=+Inf,line.search=FALSE) {
  # check parameters
  if (maxCP<3) stop("maxCP should be >=3")

  # intialize first point estimation
  neval <- 1
  w <- riskFun(w0)
  at <- as.vector(gradient(w))
  bt <- as.vector(lvalue(w)) - crossprod(as.vector(w),at)
  f <- LAMBDA*sum(abs(w)) + lvalue(w)

  # initialize aggregated working plane and working set
  A <- rbind(at,at);b <- c(bt,bt);
  inactivity.score <- c(NA_real_,NA_real_)
  ub <- f;ub.w <- w;ub.t <- 2
  for (i in 1:MAX_ITER) {
     # optimize underestimator
     opt <- lp(compute.sens = 1, direction = "max",
               objective.in = c(-1,rep_len(-LAMBDA,2L*ncol(A))),
               const.mat = cbind(-1,A,-A), const.dir = rep("<=",nrow(A)),const.rhs = -b
     )
     if (opt$status!=0) warning("issue in the LP solver:",opt$status)
     opt$W <- matrix(opt$solution[-1L],ncol=2L)
     alpha <- opt$duals[1:nrow(A)]
     
     # compute the optimum point and corresponding objective value    
     w <- opt$W[,1L] - opt$W[,2L]
     lb <- -opt$objval
     
     if (line.search) {
       w <- wolfe.linesearch(riskFun,x0=ub.w,s0=w-as.vector(ub.w),f.adjust=function(w){
         lvalue(w) <- lvalue(w) + LAMBDA*sum(abs(w))
         gradient(w) <- as.vector(gradient(w) + LAMBDA*sign(w))
         w
       })
       neval <- neval + attr(w,"neval")
     } else {
       neval <- neval + 1
     }
     
     # estimate loss at the new underestimator optimum
     w <- riskFun(w)
     f <- LAMBDA*sum(abs(w)) + lvalue(w)
             
    # update inactivity score
    inactivity.score <- inactivity.score + pmax(1-alpha,0)
    
    # deduce parameters of the new cutting plane
    at <- as.vector(gradient(w))
    bt <- lvalue(w) - crossprod(as.vector(w),at)

    # update aggregated cutting plane
    A[1,] <- alpha %*% A
    b[1] <- alpha %*% b
    
    # add the new cutting plane to the working set
    if (nrow(A)<maxCP) {
      A <- rbind(A,at);b <- c(b,bt);
      inactivity.score <- c(inactivity.score,0)
      t <- length(b)
    } else {
      t <- which.max(inactivity.score)
      A[t,] <- at;b[t] <- bt;
      inactivity.score[t] <- 0
    }
    
    # if new best found
    if (f<ub) {
      inactivity.score[ub.t] <- 0
      ub <- f;ub.w <- w;ub.t <- t
      inactivity.score[ub.t] <- NA
    }
    
    # test end of the convergence
    cat(sprintf("%d:ncall=%d gap=%g obj=%g reg=%g risk=%g w=[%g,%g]\n",i,neval,ub-lb,ub,LAMBDA*sum(abs(ub.w)),lvalue(ub.w),min(ub.w),max(ub.w)))
    if (ub-lb < ub*EPSILON_TOL) break
  }
  if (i >= MAX_ITER) warning('max # of itertion exceeded')
  attr(ub.w,"obj") <- ub
  return(ub.w)
}

