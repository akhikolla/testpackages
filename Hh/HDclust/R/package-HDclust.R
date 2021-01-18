#' @title Clustering high dimensional data with Hidden Markov Model on Variable Blocks
#' @name HDclust-package
#' @aliases HDclust
#' @description Clustering of high dimensional data with Hidden Markov Model on Variable Blocks (HMM-VB)
#' fitted via Baum-Welch algorithm. Clustering is performed by the Modal Baum-Welch
#' algorithm (MBW), which finds modes of the density function.
#' @details 
#' For a quick introduction to \pkg{HDclust} see the vignette \href{../doc/HDclust.html}{\code{vignette("HDclust")}}.
#' @docType package
#' @author{ Lin Lin, Yevhen Tupikov, Lixiang Zhang and Jia Li.
#' 
#' Maintainer: Jia Li \email{jiali@@psu.edu}
#' }
#' @references 
#' Lin Lin and Jia Li, "Clustering with hidden Markov model on variable blocks," \strong{Journal of Machine Learning Research}, 18(110):1-49, 2017.
#' @seealso \code{\link{hmmvbTrain}}, \code{\link{hmmvbClust}}
#' @examples
#' data("sim3")
#' set.seed(12345)
#' Vb <- vb(2, dim=40, bdim=c(10,30), numst=c(3,5), varorder=list(c(1:10),c(11:40)))
#' hmmvb <- hmmvbTrain(sim3[,1:40], VbStructure=Vb)
#' clust <- hmmvbClust(sim3[,1:40], model=hmmvb)
#' show(clust)
NULL
