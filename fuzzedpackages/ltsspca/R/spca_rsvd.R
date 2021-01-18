#########################################
#' @title Sparse Principal Component Analysis via Regularized Singular Value Decompsition (sPCA-rSVD)
#' @description the function that computes sPCA_rSVD
#' @param x the input data matrix
#' @param k the maximal number of PC's to seach for in the initial stage
#' @param method threshold method used in the algorithm; If \code{method = "hard"} (defauls), the hard threshold function is used;
#'     if \code{method = "soft"}, the soft threshold function is used; if \code{method = "scad"}, the scad threshold function is used
#' @param center if \code{center = TRUE} the data will be centered by the columnwise means; default is \code{center = FALSE}
#' @param scale if \code{scale = TRUE} the data will be scaled by the columnwise standard deviations; default is \code{scaled = FALSE}
#' @param  l.search a list of length kmax which contains the search grids chosen by the user; default is NULL
#' @param ls.min the smallest grid step when searching for the sparsity of each PC; default is 1
#' @return an object of class "sPCA_rSVD" is returned \cr
#' \item{loadings}{the sparse loading matrix estimated with sPCA_rSVD}
#' \item{scores}{the estimated score matrix}
#' \item{eigenvalues}{the estimated eigenvalues}
#' \item{spca.it}{the list that contains the results of sPCA_rSVD when searching for the individual PCs}
#' \item{ls}{the list that contains the final search grid for each PC direction}
#' @references Shen, H. and Huang, J. (2008), ``Sparse principal component anlysis via regularized low rank matrix decomposition'', \emph{Journal of Multivariate Analysis}, 99, 1015--1034.
#' @references Shen, D., Shen, H., and Marron, J. (2013). ``Consistency of sparse PCA in high dimensional low sample size context'', \emph{Journal of Multivariate Analysis}, 115, 315--333.
#' @examples
#' \dontrun{
#' nonrobM <- sPCA_rSVD(x = x, k = 2, center =  T, scale = F)
#' }
sPCA_rSVD <- function(x, k, method="hard", center = FALSE, scale=FALSE,  l.search = NULL, ls.min = 1){

  x <- as.matrix(x)
  n <- dim(x)[1]
  p <- dim(x)[2]
  v <- matrix(NA,p,k)
  Pb <- matrix(NA,p,k)

  xs <- scale(x,center=center,scale=scale)
  xe <- xs

  ls <- list()
  spca.it <- list()

  l.search0 <- l.search

  for (iter in 1:k){

    if (is.null(l.search0)){

      l.max <- p
      l.min <- 1

      l.s <- Inf

      while (l.s > ls.min){

        l.search <- unique(round(seq(l.max,l.min,length.out = 10)))
        spca.it[[iter]] <- sPCA.iterate(xe, method=method, k=1, l.search = l.search)

        ix.bic <- which.min(spca.it[[iter]]$bic)
        l.max <- l.search[max((ix.bic-1),1)]
        l.min <- l.search[min((ix.bic+1),length(l.search))]
        l.s <- (l.max - l.min + 1)/10
      }

      l.search <- seq(l.max,l.min,by=-ls.min)
      spca.it[[iter]] <- sPCA.iterate(xe,method=method,k=1, l.search = l.search)

    }else{
      l.search <- l.search0[[iter]]
      spca.it[[iter]] <- sPCA.iterate(xe,method=method,k=1, l.search = l.search)
    }

    ls[[iter]] <- l.search
    v[,iter] <- spca.it[[iter]]$v
    xe <-  Deflate.PCA(xe, v[,iter])
  }

  scores <-  xs %*% v
  eigenvalues <- apply(scores,2,var)

  return(list(loadings = v, scores=scores, eigenvalues= eigenvalues, ls = ls, spca.it = spca.it))
}


## sPCA.iterate provides the main function to search for the first PC for the deflated matrix in sPCA_rSVD
sPCA.iterate <- function(x,method="hard",k=1,  l.search = l.search){

  x <- as.matrix(x)
  n <- dim(x)[1]
  p <- dim(x)[2]

  svd.x <- svd(x)

  u0 <-  as.matrix(svd.x$u[,1:k] ,ncol=1)
  v0 <- as.matrix(svd.x$v[,1:k]* svd.x$d[1:k],ncol=1)

  v.err <- sum((x - u0%*%t(v0))^2)

  bic <- rep(NA,length(l.search))
  s.init <- sum(apply((x - u0 %*%t(v0))^2,2,sum))

  for(l in 1:(length(l.search))){

    s <- Inf;   res <- Inf
    v <- v0

    while (abs(res) > 1e-6){
      s0 <- s
      tl <- sort(abs(v),decreasing = T)[l.search[l]]
      v <- apply(v,2,threshold,tl,method=method)
      ua <- x %*% v
      u <- ua %*% diag(apply(ua^2,2,sum)^(-1/2),k,k)
      x.h <- u %*% t(v)
      s <- sum(apply((x - x.h)^2,2,sum))
      res <- s0 - s
    }
    bic[l] <- (s/(v.err)) + (sum(apply(v^2,1,sum)!=0))*log(n*p)/(n*p)
  }

  l.opt <- l.search[which.min(bic)]

  s <- Inf;   res <- Inf
  v <- v0

  while (abs(res) > 1e-6){
    s0 <- s
    tl <- sort(abs(v),decreasing = TRUE)[l.opt]
    v <- apply(v,2,threshold,tl,method=method)
    ua <- x %*% v
    u <- ua %*% diag(apply(ua^2,2,sum)^(-1/2),k,k)
    x.h <- u %*% t(v)
    s <- sum(apply((x - x.h)^2,2,sum))
    res <- s0 - s
  }

  v <-  v %*% diag(apply(v^2,2,sum)^(-1/2),k,k)

  return(list(v=v,l.opt=l.opt,bic=bic))
}


## threshold function whcih corresponds to different penalty function
threshold <- function(x,lambda,method =c("hard","soft","scad"), a=3.7){
  return(switch(method,
                "hard" =  x * ((abs(x) - lambda)>= 0),
                "soft" = sign(x) * (pmax(abs(x) - lambda,0)),
                "scad" = (sign(x) * (pmax(abs(x) - lambda,0)))*((abs(x) - 2*lambda) <= 0)
                + ((((a - 1)*x - sign(x)*a*lambda))/(a-2))*((abs(x) - 2*lambda) > 0)*((abs(x) - a*lambda) <= 0)
                + x * (abs(x) - a * lambda > 0)))
}



