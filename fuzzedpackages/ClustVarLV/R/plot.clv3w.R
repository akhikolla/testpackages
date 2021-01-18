#'@title Graphical representation of the CLV3W hierarchical clustering stages
#'
#'@description
#'This function plots either the CLV3W dendrogram or the variations of the consolidated CLV3W criterion.
#'
#' @param x : an object of class \code{clv3w}
#' @param type : What to plot. \cr
#'      "dendrogram" : the dendrogram of the hierchical clustering algorithm, \cr
#'      "delta"      : a barplot showing the variation of the clustering criterium after consolidation.
#' @param cex : Character expansion for labels.
#' @param \dots Additional arguments passed on to the real \code{print}.
#'
#' @seealso CLV3W
#'
#' @export


plot.clv3w =  function (x, type="dendrogram",cex=0.8,...)
{
  resclv3w <- x
  if (!inherits(resclv3w, "clv3w"))    stop("non convenient object")

  if(!inherits(resclv3w, "clv3wHCA")) stop("plot_clv3w() only applied on hierarchical analysis")

  if (type=="dendrogram") {
    plot(as.dendrogram(resclv3w$hclust), type ="rectangle",  main="CLV3W Dendrogram", axes=TRUE, cex=cex)
  }

  if (type=="delta") {
    appel      <- as.list(resclv3w$call)
    moddendoinertie <- eval(appel$moddendoinertie)
    p          <- dim(resclv3w$param$X)[[2]]
    gmax       <- resclv3w$param$gmax
    results<-resclv3w$tabres

    delta <- results[,4]
    delta[(length(delta)-gmax+3):length(delta)] <- diff(results[(length(delta)-gmax+2):length(delta),7])


    dev.new() ## tester pour le faire en delta cumsum ou en delta prendre en compte le gmax
    if (moddendoinertie==FALSE){ ## le faire en delta cumsum
      graphics::barplot(cumsum(delta)[(length(delta)-gmax+2):length(delta)],col=4,xlab="Nb clusters after aggregation", ylab="cum.delta", main="Evolution of the aggregation criterium",axisnames=TRUE,names.arg=paste(gmax:2,"->",(gmax-1):1),las=2,cex.names=cex,cex.main = 0.8) }
    else if (moddendoinertie==TRUE) {
      graphics::barplot(delta[(length(delta)-gmax+2):length(delta)],col=4,xlab="Nb clusters after aggregation", ylab="delta", main="Evolution of the aggregation criterium",axisnames=TRUE,names.arg=paste(gmax:2,"->",(gmax-1):1),las=2,cex.names=cex,cex.main = 0.8) }

  }


}
