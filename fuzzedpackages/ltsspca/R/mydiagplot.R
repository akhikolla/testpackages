#' @title Make diagnostic plot using the estimated PC subspace
#' @param x the input data matrix
#' @param obj the returned output from rwltsspca
#' @param k dimension of the PC subspace
#' @param alpha the robust parameter which takes value between 0 to 0.5, default is 0.5
#' @param co.sd cutoff value for score outlier weight, default is 0.25
#' @return the diagnostics of outliers \cr
#' \item{od}{the orthgonal distances with respect to the k-dimensional PC subspace}
#' \item{ws.od}{if the observation is outlying in the orthgonal complement of the PC subspace \code{ws.od}=0; otherwise \code{ws.sd}=1}
#' \item{co.od}{the cutoff value for orthogonal distances}
#' \item{sc.wt}{the score outlier weight, which is compared with 0.25 (by default) to flag score outliers}
#' \item{ws.sd}{if the observation is outlying with the PC subspace \code{ws.sd}=0; otherwise \code{ws.sd}=1}
#' \item{co.sd}{the cutoff value for score outlier weight, default is 0.25}
#' \item{sc.out}{the retruned object when computing the score outlier weights}
mydiagPlot <- function(x, obj, k, alpha = 0.5, co.sd = 0.25){

  ws.od <- obj$ws.od
  ws.sd <- obj$ws.sd

  xs  <- scale(x, center = obj$mu,scale = F)

  od.rw <-  sqrt(apply((xs - (xs %*% obj$loadings[,1:k]) %*% t(obj$loadings[,1:k]))^2,1,sum))
  co.od.rw <- coOD(od = od.rw,  h =(dim(x)[1]-floor(((dim(x)[1])*alpha))))$co.od
  ws.od.rw <- od.rw <= co.od.rw

  iout.new <- IdOUT(obj$scores[,1:k],ws.od=ws.od.rw,outbound=co.sd)
  ws.sd.rw <- as.logical(iout.new$wfinal01)

  ws00 <- which(ws.sd.rw &  ws.od.rw)
  ws10 <- which((!ws.sd.rw ) & ws.od.rw)
  ws01 <-  which( ws.sd.rw  & (! ws.od.rw))
  ws11 <- which((! ws.sd.rw ) & (!ws.od.rw))

  min.od <- min(od.rw)/co.od.rw
  max.od <- max(od.rw)/co.od.rw
  s.od <- (max.od - max.od) / length(ws00)

  opar <- par(mar=c(5,5,5,1),no.readonly = TRUE)
  on.exit(par(opar))

  plot(iout.new$wfinal[ws00],od.rw[ws00]/co.od.rw,
                 main = "Diagnostic plot", xlab="Score outlier weight", ylab="Scaled OD", cex.lab = 2, cex.main=2,cex=1.5,
                 pch=16, col=1,ylim=c((min.od-s.od),(max.od+s.od)),xlim=c(1,0)) +
    points(iout.new$wfinal[ws10],od.rw[ws10]/co.od.rw,pch = 18, col= 3,cex=1.5) +
    points(iout.new$wfinal[ws01],od.rw[ws01]/co.od.rw,pch = 15, col= 4,cex=1.2) +
    points(iout.new$wfinal[ws11],od.rw[ws11]/co.od.rw,pch = 17, col= 2,cex=1.5) +
    abline(h=1,v=co.sd)

  return(list(od = od.rw, ws.od = ws.od.rw, sc.wt = iout.new$wfinal, ws.sd = ws.sd.rw, co.od = co.od.rw, co.sd = 0.25, sc.out = iout.new))
}
