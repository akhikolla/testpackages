#' @title Weights of the external variables, or additional mode, on the latent component in each cluster.
#' 
#' @description To get the weigths associated with each cluster.
#'  For CLV(_kmeans) or LCLV, applies only when external variables (Xr, Xu or both) are involved.
#'  For CLV3W(_kmeans), the weights are associated with the third mode of the 3-way array.
#' 
#' @param resclv : result of CLV(), CLV_kmeans() or LCLV(), CLV3W() or CLV3W_kmeans()
#' @param K : the number of clusters chosen (already defined if CLV_kmeans or CLV3W_kmeans is used)
#' @param graph : boolean, if TRUE, the barplot associated with the scores is displayed (default : graph=FALSE)
#' @param cex.lab : magnification to be used for labels (1 by default)
#' 
#' @return \item{weight}{Weights in each cluster (associated with mode 3 for CLV3W object) \cr
#'       For each cluster, the vector of weights is set to length 1 \cr
#'       output provided as matrix with K columns (K: number of clusters) \cr
#'       In the special case of LCLV, two matrices of weights are defined : \cr
#'       weight_v : weights of the external Xr variables, \cr
#'       weight_u : weights of the external Xu variables.
#'      }
#'       
#' @export
#' 
get_weight <-  function(resclv, K=NULL,graph=FALSE,cex.lab=1)
  {
    
   if (!inherits(resclv, c("clv","lclv","clv3w")))    stop("non convenient objects")

   # for object clv
   if(inherits(resclv,c("clv","lclv"))) {  
     if(is.null(resclv$param$K)) { 
       if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
     }
     weight=get_weight.clv(resclv,K,graph=graph,cex.lab=cex.lab) 
    }
   
   # for object clv3w
   if(inherits(resclv,"clv3w"))   {  
     appel      <- as.list(resclv$call)
     if(is.null(eval.parent(appel$K))) { 
       if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
     }
     weight=get_weight.clv3w(resclv,K,graph=graph,cex.lab=cex.lab)
   }
   
   return(weight=weight)
   
  }