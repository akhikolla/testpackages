


#' Split a dataset for Cross Validation taking into account class balance
#' 
#' @param y the class labels of each sample of the dataset
#' @param num.cv number of cross validation required
#' @return a factor of num.cv levels that assign to each sample a test fold
#' @export
balanced.cv.fold <- function(y,num.cv=10) {
  fold <- factor(rep_len(sample(seq_len(num.cv)),length.out=length(y)))
  i <- sample(seq_along(y))
  o <- order(y[i])
  fold[i[o]] <- fold
  fold
}

#' Compute loss.weights so that total losses of each class is balanced
#' 
#' @param y a object coerced to factor that represent the class labels of each sample of the dataset
#' @return a numeric vector of the same length as y
#' @export
balanced.loss.weights <- function(y) {
  #if (is.logical(y)) return(ifelse(y,mean(!y),mean(y)))
  y <- as.factor(y)
  cw <- 1/(tabulate(y)/length(y))
  cw <- cw / sum(cw)
  cw[y]
}


#' Rank linear weight of a linear model
#' 
#' @param w a numeric vector of linear weights
#' @return a data.frame with a rank for each feature as well as z-score, p-value, and false discovery rate.
#' @export
#' @importFrom stats sd pnorm p.adjust
#' @importFrom methods as
rank.linear.weights <- function(w) {
  w <- as(w,"vector")
  R <- data.frame(stringsAsFactors = FALSE,
    feature.name = if (is.null(names(w))) seq_along(w) else names(w),
    w = w,
    rk = pmin(rank(+w,ties.method="first"),rank(-w,ties.method="first"))
  )
  R$z <- as.vector((w-mean(w))/sd(w))
  R$pval <- as.vector(2*pnorm(abs(R$z),lower.tail=FALSE))
  R$fdr <- as.vector(p.adjust(R$pval,method="fdr"))
  R
}


#' Compute statistics for ROC curve plotting
#' 
#' @param f decision value for each instance
#' @param y a logical that specify binary labels
#' @return a data.frame() that compute for each threshold value 'f' roc curve statistics: TP, FP, TN, FN, FPR, TPR, sensitivity, specificity, precision, recall, accuracy
#' @author Julien Prados, adapted from Bob Horton code
#' @export
#' @examples
#'   x <- cbind(data.matrix(iris[1:4]))
#'   w <- nrbmL1(rocLoss(x,iris$Species=="versicolor"),LAMBDA=0.01)
#'   plot(roc.stat(x %*% w,iris$Species=="versicolor"))
#'   lines(roc.stat(-x[,2],iris$Species=="versicolor"),col="blue")
roc.stat <- function(f,y) {
  if (!is.logical(y)) stop("y must be a logical vector")
  if (length(f)!=length(y)) stop("f and y must have same length")
  o <- order(f, decreasing=TRUE)
  roc <- data.frame(
    f = c(-Inf,f[o]),
    TP = c(0,cumsum(y[o])),
    FP = c(0,cumsum(!y[o]))
  )
  roc <- roc[rev(!duplicated(rev(roc$f))),]
  roc$TN <- sum(!y) - roc$FP
  roc$FN <- sum(y) - roc$TP
  roc$FPR <- roc$FP/sum(!y)
  roc$sensitivity <- roc$recall <- roc$TPR <- roc$TP/sum(y)
  roc$accuracy <- (roc$TP+roc$TN)/length(y)
  roc$specificity <- roc$TNR <- roc$TN / (roc$TN + roc$FP)
  roc$precision <- roc$TP/(roc$TP+roc$FP)
  
  dx <- diff(roc$FPR)
  dy <- diff(roc$TPR)
  attr(roc,"AUC") <- sum(dx*roc$TPR[seq_along(dx)] + dx*dy/2)
  class(roc) <- c("roc.stat",class(roc))
  return(roc)
}

#' Generic method overlad to print object of class roc.stat
#' 
#' @param x a roc.stat object return by the function roc.stat
#' @param ... additional parameters
#' @export
print.roc.stat <- function(x,...) {
  NextMethod()
  cat(sprintf("AUC: %.2f\n",attr(x,"AUC")))
}

#' @export
#' @importFrom graphics plot
plot.roc.stat <- function(x,y,type="o",xlab="sensitivity",ylab="specificity",...) {
  plot(x$sensitivity,x$specificity,xlab=xlab,ylab=ylab,type=type,...)
}

#' @export
#' @importFrom graphics lines
lines.roc.stat <- function(x,type="o",...) {
  lines(x$sensitivity,x$specificity,type=type,...)
}


#' Perform multiple hierachical clustering on random subsets of a dataset
#' 
#' @param x the numeric matrix containing the data to cluster (one instance per row)
#' @param seeds a vector of random seed to use.
#' @param row.rate,col.rate numeric value in [0,1] to specify the proportion of instance 
#'        (resp. feature) to subset at each random iteration.
#' @param max.cluster upper bound on the number of expected cluster (can by +Inf).
#' @param ret.height a logical to specify whether the average merging height should be returned.
#' @param hc.method a clustering method of arity 1, taking as input a random subset of the 
#'        input matrix x and returning an hclust object
#' @param ... additional arguments are passed to the hc.method
#' @return a list of 3 square matrices N,H,K of size nrow(x): N is the number of 
#'         time each pair of instance as been seen in the random subsets; H is the
#'         corresponding sum of heights for the pairs; K is the sum of the number of split
#'         possible that still preserve the two samples into the same cluster.
#' @author Julien Prados
#' @importFrom stats as.dist dist hclust prcomp
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
iterative.hclust <- function(x,seeds=1:100,
  row.rate=0.3,col.rate=0.1,max.cluster=10L,
  ret.height=FALSE,
  hc.method=function(x,PCs=1:6,...) {hclust(dist(prcomp(x,rank.=max(PCs))$x[,PCs,drop=FALSE]),...)},
  ...
) {
  ret <- list()
  ret$K <- ret$N <- matrix(0L,nrow(x),nrow(x))
  if (ret.height) ret$H <- matrix(0,nrow(x),nrow(x))

  pb <- txtProgressBar(0,length(seeds),style=3)
  for(k in seq_along(seeds)) {
    set.seed(seeds[k])
    i <- sort(sample(nrow(x),nrow(x)*row.rate))
    j <- sample(ncol(x),ncol(x)*col.rate)
    hc <- hc.method(x[i,j],...)
    
    ri <- rep(seq_along(hc$order),rev(seq_along(hc$order)))
    ci <- seq_along(ri) + cumsum(seq_along(hc$order)-1L)[ri]
    ci <- (ci-1L)%%length(hc$order) + 1L
    A <- hclust_fca(hc,ri,ci)

    ij <- cbind(i[ci],i[ri])
    if (ret.height) ret$H[ij] <- ret$H[ij] + hc$height[A]
    ret$N[ij] <- ret$N[ij] + 1L
    A <- nrow(hc$merge)+1L-A
    A[A>max.cluster] <- max.cluster
    ret$K[ij] <- ret$K[ij] + A
    setTxtProgressBar(pb,k)
  }
  ret$K <- as.dist(ret$K)
  ret$N <- as.dist(ret$N)
  if (ret.height) ret$H <- as.dist(ret$H)
  ret
}




#' Compute Bhattacharyya coefficient needed for Hellinger distance
#' 
#' @param x a numeric matrix
#' @return a square matrix containing Bhattacharyya coefficient for each pair of row in x
#' @author Julien Prados
#' @export
bhattacharyya.coefficient <- function(x) {
  pmin(pmax(crossprod(sqrt(t(x))),0),1)
}

#' Compute Hellinger distance
#' 
#' @param x a numeric matrix
#' @return an object of class "dist" with Hellinger distance of each pair of row in x
#' @author Julien Prados
#' @export
hellinger.dist <- function(x) {
  as.dist(sqrt(1-bhattacharyya.coefficient(x)))
}





