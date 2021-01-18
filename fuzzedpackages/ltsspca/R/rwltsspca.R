#' @title  Reweighted LTS-SPCA
#' @description the function that computes the reweighted LTS-SPCA
#' @param x the input data matrix
#' @param obj initial LTS-SPCA object given by ltsspca function
#' @param k dimension of the PC subspace; by default is NULL then k takes the value of kmax in the initial LTS-SPCA
#' @param alpha the robust parameter which takes value between 0 to 0.5, default is 0.5
#' @param co.sd cutoff value for score outlier weight, default is 0.25
#' @return the object of class "ltsspcaRw" is returned \cr
#' \item{loadings}{the sparse loading matrix estimated with reweighted LTS-SPCA}
#' \item{scores}{the estimated score matrix}
#' \item{eigenvalues}{the estimated eigenvalues}
#' \item{mu}{the center estimate}
#' \item{rw.obj}{the list that contains the results of sPCA_rSVD on the reduced data}
#' \item{od}{the orthonal distances with respect to the initially estimated PC subspace with all the noisy variables removed}
#' \item{co.od}{the cutoff value for the orthogonal distances}
#' \item{ws.od}{if the observation is outlying in the orthgonal complement of the initially estimated PC subspace \code{ws.od}=0; otherwise \code{ws.od}=1}
#' \item{sc.wt}{the score outlier weight, which is compared with 0.25 (by default) to flag score outliers}
#' \item{co.sd}{the cutoff value for score outlier weight, default is 0.25}
#' \item{ws.sd}{if the observation is outlying with the PC subspace \code{ws.sd}=0; otherwise \code{ws.sd}=1}
#' \item{sc.out}{the retruned object when computing the score outlier weights}
ltsspcaRw <- function(x, obj, k = NULL, alpha = 0.5, co.sd = 0.25){

  if (is.null(k)){
    ks <- dim(obj$loadings)[2]
  }else{
    ks <- k
  }

  mu <-  apply(as.matrix(obj$mu[,1:ks],ncol=ks),1,sum)
  ind <- which(apply(obj$loadings[,1:ks]^2,1,sum) != 0)
  h <- dim(x)[1]-floor(((dim(x)[1])*alpha))
  xc <- scale(x[,ind],center = mu[ind],scale = F)


  ##detecting orthogonal outliers with noisy variables removed
  OD <- sqrt(apply((xc - xc%*%obj$loadings[ind,1:ks]%*%t(obj$loadings[ind,1:ks]))^2,1,sum))
  co.od <- coOD(od = OD, h =(dim(x)[1]-floor(((dim(x)[1])*alpha))))$co.od
  ws.od <- OD <= co.od

  ##detecting orthogonal outliers within the PC subspace
  scores <- scale(x[,ind],center = mu[ind],scale = F) %*% obj$loadings[ind,1:ks]
  iout <- IdOUT(scores,ws.od=ws.od,outbound=co.sd)
  ws.sd <- as.logical(iout$wfinal01)

  ##apply sPCA_rSVD on the reduced data set
  rw.spca.x <- sPCA_rSVD(x[(ws.od & ws.sd),ind],center = T, scale = F , k=ks)

  loadings <- matrix(0, nrow = dim(x)[2], ncol = ks)
  loadings[ind,] <-  rw.spca.x$loadings
  mu[ind] <- apply(x[ws.od&ws.sd,ind],2,mean)
  scores <-  scale(x[,ind], center = mu[ind], scale = F) %*% loadings[ind,]
  eigenvalues <- apply(scores[(ws.od & ws.sd),],2,var)

  ix <- sort(eigenvalues,decreasing = T,index=T)$ix

  return(list(loadings = loadings[,ix], scores = scores[,ix], eigenvalues = eigenvalues[ix], mu = mu, rw.obj = rw.spca.x,
              od = OD, co.od = co.od,  ws.od = ws.od, wt.sc = iout$wfinal, co.sd = co.sd, ws.sd = ws.sd, sc.out = iout))
}

