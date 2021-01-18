#' Print the CLV3W results
#'
#' @param x an object of class \code{clv3w}
#' @param \dots Additional arguments passed on to the real \code{print}.
#'
#' @seealso CLV3W, CLV3W_kmeans
#'
#' @export
#'
print.clv3w =  function (x, ...)
{
  if (!inherits(x, "clv3w"))
    stop("non convenient object")
  resclv3w   <- x
  appel      <- as.list(resclv3w$call)
    n <- dim(resclv3w$param$X)[[1]]
    p <- dim(resclv3w$param$X)[[2]]
    q <- dim(resclv3w$param$X)[[3]]



  NN         <-  eval.parent(appel$NN)

  cat("\n")
  cat(paste("number of observations in mode 1 : n=", n), sep = " ")
  cat("\n")
  cat(paste("number of variables in mode 2 : p=", p), sep = " ")
  cat("\n")
  cat(paste("number of parameters in mode 3 : q=", q), sep = " ")
  cat("\n")

  if (NN) cat("Cluster Analysis of mode 2 associated with a one-rank PARAFAC model. Non negativity is set on the loadings of mode 2 variables.")
  else
    cat("Cluster Analysis of mode 2 associated with a one-rank PARAFAC model. ")
  cat("\n")
  if (inherits(resclv3w,"clv3wHCA")) {
    gmax <- resclv3w$param$gmax
    cat(paste("Consolidation for K in c(",gmax,":2)",sep = ""))
    cat("\n")
    cat("\n")
    cat("$tabres: results of the hierarchical clustering")
    cat("\n")
    cat("$partitionK or [[K]]: partition into K clusters")
    cat("\n")
    cat("    [[K]]$clusters: cluster's membership (1st line: before and 2nd line: after consolidation)")
    cat("\n")
    cat("    [[K]]$comp: latent components of the clusters (after consolidation),matrix of size (n x K)")
    cat("\n")
    cat("    [[K]]$loading: loadings  associated with the second mode by cluster (after consolidation),matrix of size (p x K)")
    cat("\n")
    cat("    [[K]]$weight: weights associated with the third mode by cluster (after consolidation),matrix of size (q x K)")
    cat("\n")
    cat("    [[K]]$criterion: loss criterion giving the residual amount between the sub-array and its reconstitution ")
    cat("\n")
    cat("                              obtained by the one rank PARAFAC model (after consolidation),vector of size K")
    cat("\n")
  } else {
    cat(paste("number of clusters: ", eval.parent(appel$K)), sep = " ")
    cat("\n")
    cat("\n")
    cat("$clusters: cluster's membership (1st line: intial partition and 2nd line: partition at convergence")
    cat("\n")
    cat("$comp: latent components of the clusters;matrix of size (n x K)")
    cat("\n")
    cat("$loading: loadings  associated with the second mode by cluster,matrix of size (p x K)")
    cat("\n")
    cat("$weight: weights associated with the third mode by cluster,matrix of size (q x K)")
    cat("\n")
    cat("$criterion: loss criterion giving the residual amount between the sub-array and its reconstitution ")
    cat("\n")
    cat("             obtained by the one rank PARAFAC model (after consolidation),vector of size K")
    cat("\n")
  }
}
