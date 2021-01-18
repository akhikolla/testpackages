#' @title clusters memberships for a partition into K clusters. 
#' 
#' @description To get the clusters memberships of the variables. 
#' 
#' @param resclv : result of CLV(), CLV_kmeans(), LCLV(), CLV3W() or CLV3W_km()
#' @param K : the number of clusters chosen (already defined if CLV_kmeans or CLV3W_kmeans is used)
#' @param type : "vector" (by default) for output given as a vector of integers between 1 and K (with 0 for "kplusone" strategy), \cr
#'               "matrix", the output given as a binary matrix of size p x n.
#'
#' @return \item{partition}{the group's membership for the variables, in a vector or matrix form. \cr
#'                      For CLV3W object, a vector of memberships with mode 2) }
#'   
#' @examples data(apples_sh)
#' resclvX <- CLV(X = apples_sh$senso, method = "directional", sX = TRUE)
#' parti4G<-get_partition(resclvX, K = 4) 
#' 
#' @export
#' 
get_partition <-
  function(resclv, K=NULL,type="vector")
  {
    
    
    if (!inherits(resclv, c("clv","lclv","clv3w"))) 
      stop("non convenient objects")
    
    X<-resclv$param$X
    libel<-colnames(X)
    
    if(inherits(resclv,c("clv","lclv"))) {  
     if(is.null(resclv$param$K)) { 
      if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
       partition<-as.vector(resclv[[K]]$clusters[2,])
     } else {
       partition<-as.vector(resclv$clusters[2,])
      K<-resclv$param$K
     }
    names(partition)<-libel
    }
    
    if(inherits(resclv,"clv3w")) {  
     appel      <- as.list(resclv$call)
     if(is.null(eval.parent(appel$K))) { 
      if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
      if ((K > resclv$param$gmax+1) | (K<1))
        stop(paste("Consolidation not performed for ",K,"clusters"))
      partition     <- as.vector(resclv[[K]]$clusters[2,])
     } else {
       K<-eval.parent(appel$K)
       partition     <- as.vector(resclv$clusters[2,])
     } 
     names(partition)<-libel
    }
    
   
    if (type=="vector") return(partition=partition)
    if (type=="matrix") { 
      partition<-as.data.frame(partition)
      names(partition)<-"G"
      tab<-tabdisj(partition)
      return(partition=tab)
    }
  }