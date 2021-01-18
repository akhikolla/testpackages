# weights associated with a clv3w class
#
#
# To get the weights (mode 3) associated with each cluster.
#
# resclv3w : result of CLV3W(), CLV3W_kmeans()
# K : the number of groups chosen (already defined if CLV3W_kmeans is used)
# graph : boolean, if TRUE, the barplot associated with the scores is displayed (default : graph=FALSE)
# cex.lab : magnification to be used for labels (1 by default)
#
#
get_weight.clv3w <- function(resclv3w, K=NULL,graph=FALSE,cex.lab=1) {
  if (!inherits(resclv3w, "clv3w") )   stop("non convenient object")
  appel      <- as.list(resclv3w$call)
  if(is.null(eval.parent(appel$K))) {
    if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
    if ((K > resclv3w$param$gmax+1) | (K<1))
      stop(paste("Consolidation not performed for ",K,"clusters"))
    weight     <- resclv3w[[K]]$weight
  } else {
    K<-eval.parent(appel$K)
    weight     <- resclv3w$weight
  }
  if (graph) {
    dev.new()
    barplot(t(weight), ylab = 'Weights',beside=TRUE,names.arg=rownames(weight),cex.names=cex.lab,legend.text=paste('Cluster',1:K),las=2)
    title( main = 'Weights associated with the CLV3W partition ' )

  }
  return(weight=weight)
}
