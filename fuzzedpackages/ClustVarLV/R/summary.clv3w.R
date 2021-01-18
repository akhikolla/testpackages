#' @title Summary and description of the clusters of (mode 2) variables associated with CLV3W or CLV3W_kmeans
#'
#' @description This function provides the list of the variables within each group and complementary informations.
#' Users will be asked to specify the number of clusters,
#'
#' @param object : result of CLV3W() or CLV3W_kmeans()
#' @param K : the number of clusters (unless if CLV3W_kmeans was used)
#' @param \dots Additional arguments passed on to the real \code{summary}.
#'
#'
#' @details The ouputs include :
#' \itemize{
#' \item the size of the groups, \cr
#' \item the proportion of the variance within each group explained by its latent variable, \cr
#'  \item the proportion of the whole dataset accounted by the group latent variables \cr
#' \item the latent components (mode 1) associated to the various groups, \cr
#' \item the weights (mode 3) associated to the various groups, \cr
#' \item the list of the variables within each group. For each cluster, the loading (mode 2) of the variable is given together with the correlation of the block component with its group latent component
#'   and the correlation with the next neighbouring group latent component are given.  \cr
#'  \item the matrix of correlation between the latent variables.}
#'
#' @export
#' @importFrom stats sd cov var cor
#'
summary.clv3w <-
function(object,K=NULL,...) {

  resclv3w<-object
  if (!inherits(resclv3w, "clv3w") )   stop("non convenient object")
  appel      <- as.list(resclv3w$call)
  X          <- resclv3w$param$X

  n <- dim(X)[[1]]
  p <- dim(X)[[2]]
  q <- dim(X)[[3]]

  NN         <-  eval.parent(appel$NN)
# group membership of the variables
  if(is.null(eval.parent(appel$K))) {
    if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
    if ((K >0) & (K<resclv3w$param$gmax+1)) {
      clusters <- resclv3w[[K]]$clusters[2,]
      comp     <- resclv3w[[K]]$comp
      loading  <- resclv3w[[K]]$loading
      weight   <- round(resclv3w[[K]]$weight,3)
      criterion<- resclv3w[[K]]$criterion
      pk <- table(clusters)
    }
    else stop("invalid number of clusters")
  } else {
    clusters<-resclv3w$clusters[2,]
    K<-eval.parent(appel$K)
    comp     <- resclv3w$comp
    loading  <- resclv3w$loading
    weight   <- round(resclv3w$weight,3)
    criterion<- resclv3w$criterion
    pk <- table(clusters)
  }

  # initialization : computation of the block component associated with each variable (mode 2)
  # computation of the inertia associated with each block
  correlation.j <- matrix(0,nrow=p,ncol=3)
  rownames(correlation.j) <- dimnames(X)[[2]]
  if (NN) {
    colnames(correlation.j) <- c("loading", "cor in.group", "cor next.group") } else {
      colnames(correlation.j) <- c("loading", "|cor| in.group", "|cor| next.group")}


  Tblock      <- matrix(0,nrow=n,ncol=p)
  inertia     <- sapply(1:K,function(k){sum(X[,which(clusters==k),]^2)})
  inertia.hat <- inertia-criterion
  if (K>1) {
    for (j in 1:p) {
      Tblock <- sapply(1:K,function(k) {apply(sweep(X[,j,],2,STATS=weight[,k],FUN="*"),1,sum)})
      cor.j <- sapply(1:K,function(k){if (var(Tblock[,k])!=0) cor(comp[,k],Tblock[,k]) else 0})
      if (NN)
        correlation.j[j,] <- c(loading[j,clusters[j]],cor.j[clusters[j]],max(cor.j[-clusters[j]]))
      else
        correlation.j[j,] <- c(loading[j,clusters[j]],abs(cor.j[clusters[j]]),max(abs(cor.j[-clusters[j]])))
    }
    correlation <- lapply(1:K,function(k){round(correlation.j[which(clusters==k),,drop=FALSE],3);
    })
    names(correlation) <- paste("Cluster",1:K)
  }
  else {
    for (j in 1:p) {
      Tblock <- apply(sweep(X[,j,],2,STATS=weight[,1],FUN="*"),1,sum)
      if (NN)
        correlation.j[j,] <- c(loading[j,1],cor(comp[,1],Tblock),0)
      else
        correlation.j[j,] <- c(loading[j,1],abs(cor(comp[,1],Tblock)),0)
    }
    correlation <- list(round(correlation.j,3));
    names(correlation) <- "Cluster 1"
  }

  inertia.perc       <-  100* (inertia -criterion)/inertia
  inertia.whole.perc <-  100* (sum(inertia) -sum(criterion))/sum(inertia)
#correlation between latent components
  correl.comp <- matrix(cor(comp),nrow=K,ncol=K)
  dimnames(correl.comp) <- list(paste("Comp",1:K,sep="."),paste("Comp",1:K,sep="."))



  sumclv3w <- list(size=pk,prop_explained_per_cluster=round(inertia.perc,3), prop_explained_total=round(inertia.whole.perc,3), comp=round(comp,3),weight=weight, groups=correlation,cormatrix=round(correl.comp,3))
  print(sumclv3w)
}
