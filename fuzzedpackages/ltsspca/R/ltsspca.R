#'@title  Sparse Principal Component Analysis Based on Least Trimmed Squaers (LTS-SPCA)
#' @description the function that computes the initial LTS-SPCA
#' @param x the input data matrix
#' @param kmax the maximal number of PCs searched by the intial LTS-SPCA
#' @param alpha the robust parameter which takes value between 0 to 0.5, default is 0.5
#' @param mu.choice the center estimate fixed by the user; by default, the center will be estimated automatically by the algorithm
#' @param  l.search a list of length kmax which contains the search grids chosen by the user; default is NULL
#' @param ls.min the smallest grid step when searching for the sparsity of each PC; default is 1
#' @param tol convergence criterion
#' @param N1  the number controls the updates for a without updating b in the concentration step for LTS-PCA
#' @param N2 the number controls outer loop in the concentration step for LTS-PCA
#' @param N2bis the number controls the outer loop for the selected b for both LTS-PCA and LTS-SPCA
#' @param Npc the number controls the inner loop for both LTS-PCA and LTS-SPCA
#' @return the object of class "ltsspca" is returned \cr
#' \item{loadings}{the initially estimated loading matrix by LTS-SPCA}
#' \item{mu}{the center estimates associated with each PC}
#' \item{spca.it}{the list that contains the results of LTS-SPCA when searching for the individual PCs}
#' \item{ls}{the list that contains the final search grid for each PC direction}
#' @author Yixin Wang
#' @references Wang, Y., Van Aelst, S. (2019), `` Sparse Principal Component Based On Least Trimmed Squares'', \emph{Technometrics, accepted}.
#' @examples
#' library(mvtnorm)
#' dataM <- dataSim(n = 200, p = 20, bLength = 4, a = c(0.9, 0.5, 0),
#'                 SD = c(10, 5, 2), eps = 0, seed = 123)
#' x <- dataM$data
#' ltsspcaMI <- ltsspca(x = x, kmax = 5, alpha = 0.5)
#' ltsspcaMR <- ltsspcaRw(x = x, obj = ltsspcaMI, k = 2, alpha = 0.5)
#' matplot(ltsspcaMR$loadings,type="b",ylab="Loadings")

ltsspca <- function(x, kmax, alpha=0.5, mu.choice=NULL,
                     l.search = NULL, ls.min = 1,tol=1e-6, N1=3, N2=2, N2bis=10, Npc=10){

  x <- as.matrix(x)

  n <- dim(x)[1]
  p <- dim(x)[2]

  a <- matrix(NA,nrow = n, ncol = kmax)
  b <- matrix(NA,nrow = p, ncol = kmax)

  xe <- x
  spca.it <- list()
  ls <- list()
  mu <- matrix(NA, nrow=p, ncol=kmax)

  l.search0 <- l.search

  for (j in 1:kmax){

    if (is.null(l.search0)){

      l.max <- p
      l.min <- 1

      l.s <- Inf

      while (l.s > ls.min){

        l.search <- unique(round(seq(l.max,l.min,length.out = 10)))
        spca.it[[j]] <- lts.spca.iterate(x = xe, q = 1, l.search = l.search, N1 = N1, N2 = N2,
                                         N2bis = N2bis, Npc = Npc, tol=tol, alpha = alpha, mu.choice = mu.choice)

        ix.bic <- which.min(spca.it[[j]]$bic)
        l.max <- l.search[max((ix.bic-1),1)]
        l.min <- l.search[min((ix.bic+1),length(l.search))]
        l.s <- (l.max - l.min + 1)/10
      }

      l.search <- seq(l.max,l.min,by=-ls.min)

      spca.it[[j]] <- lts.spca.iterate(x = xe, q = 1, l.search = l.search, N1 = N1, N2 = N2,
                                       N2bis = N2bis, Npc = Npc, tol=tol, alpha = alpha, mu.choice = mu.choice)
    }else{
      l.search <- l.search0[[j]]
      spca.it[[j]] <- lts.spca.iterate(x = xe, q = 1, l.search = l.search, N1 = N1, N2 = N2,
                                       N2bis = N2bis, Npc = Npc, tol=tol, alpha = alpha, mu.choice = mu.choice)
    }


    ls[[j]] <- l.search

    b[,j] <- spca.it[[j]]$spca.opt$b
    a[,j] <- spca.it[[j]]$spca.opt$a
    mu[,j] <- spca.it[[j]]$spca.opt$mu

    xe <- Deflate.PCA(scale(xe,center=spca.it[[j]]$spca.opt$mu,scale = F),b[,j])

  }

  return(list(loadings = b, mu = mu, spca.it = spca.it, ls = ls))
}

## lts.spca.iterate provides the main function to search for the first PC for the deflated matrix in ltsSpca
lts.spca.iterate <- function(x=x, q=q,alpha=alpha, l.search = l.search, N1=N1, N2=N2, N2bis=N2bis, Npc=Npc ,tol=tol,
                             mu.choice=mu.choice){

  x <- as.matrix(x)
  n <- dim(x)[1]
  p <- dim(x)[2]

  h <- n - floor(n*alpha)



  if(is.null(mu.choice)){
    mui <- apply(x,2,median)
  }else{
    mui <- mu.choice
  }

  lts.x <- ltspca(x = x, q = q,  alpha=alpha, N1=N1, N2=N2, N2bis=N2bis, Npc=Npc ,tol=tol)

  lts_mu <- t(lts.x$mu)
  lts_b <-  qr.Q(qr(lts.x$b))

  ## update the center given by LTS-PCA
  s.init <- scale(x,center = lts_mu,scale = F) %*% lts_b
  cov.s.init <- robustbase::covMcd(s.init, alpha = 1-alpha)
  lts_mu_new <- lts_mu + lts_b %*% cov.s.init$center

  # Improvement of the best deterministic candidate
  bic <- rep(NA,length(l.search))

  for(l in 1:length(l.search)){
    spca.0 <- findsparsePC(x = x, alpha = alpha, mu = lts_mu_new, b= lts_b, l=l.search[l], N2bis = N2bis, Npc = Npc, tol = tol)
    bic[l] <-  spca.0$s/lts.x$best.obj + (l.search[l])*(log(h*p))/(h*p)
  }

  l.opt <- l.search[which.min(bic)]
  spca.opt.0 <- findsparsePC(x=x,alpha = alpha, mu=lts_mu_new, b=lts_b, l=l.opt, N2bis = N2bis, Npc = Npc, tol = tol)

  # Results: storage of the optimal value
  res <- list(spca.opt = spca.opt.0, l.opt = l.opt, bic = bic)
  return(res)
}

