#' @title latent components associated with each cluster
#' 
#' @description To get the latent components associated with each cluster.
#' 
#' @param resclv : result of CLV(), CLV_kmeans(), LCLV(), CLV3W() or CLV3W_km()
#' @param K : the number of clusters chosen (already defined if CLV_kmeans or CLV3W_kmeans is used)
#' @param graph : boolean, if TRUE, the barplot associated with the scores is displayed (default : graph=FALSE)
#' @param cex.lab : magnification to be used for labels (1 by default)
#'
#' @return \item{comp}{the group latent components (centered) \cr
#'      For results of CLV(_kmeans), the latent components returned have their own norm \cr
#'      For results of CLV3W(_kmeans), the latent component associated with mode 1 (centered, but not standardized)  \cr
#'      For results of LCLV, two types of latent components are available : \cr
#'      compt : The latent components of the clusters defined according to the Xr variables, \cr
#'      compc : The latent components of the clusters defined according to the Xu variables
#'      }
#' 
#'         
#' @examples data(apples_sh)
#' resclvX <- CLV(X = apples_sh$senso, method = "directional", sX = TRUE)
#' comp4G<-get_comp(resclvX, K = 4) 
#' 
#' @export
#' 
get_comp <-
  function(resclv, K=NULL,graph=FALSE,cex.lab=1)
  {
    
   if (!inherits(resclv, c("clv","lclv","clv3w"))) 
      stop("non convenient objects")
    
   
   if(inherits(resclv,"clv")) {
    if(is.null(resclv$param$K)) { 
      if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
      comp<-resclv[[K]]$comp
    } else {
      comp<-resclv$comp
      K<-resclv$param$K
    }
     rownames(comp)=rownames(resclv$param$X)
     if(is.null(rownames(comp))) rownames(comp)=1:nrow(resclv$param$X)
     if (graph) {
       dev.new()
       if (nrow(comp)<10) {
         barplot(t(comp), ylab = 'Scores',beside=TRUE,names.arg=rownames(comp),cex.names=cex.lab,
                 las=2, col=1:K,legend.text=paste('Cluster',1:K), main ="Scores associated with the cluster's latent components")
       } else {
         n=nrow(resclv$param$X)
         if (n<20) {
           NC=3
         } else {
           if (n<50){
             NC=2
         } else {
             NC=1
         }}
         NL=ceiling(K/NC)
         par(mfrow=c(NL,NC))
         for(k in 1:K) {
           barplot(t(comp[,k]), ylab = 'Scores',beside=TRUE,names.arg=rownames(comp),cex.names=cex.lab,
                   las=2,  main = paste("Scores associated with the cluster",k) )
         }
       }
     }
  return(comp=comp)
   }
   
   if(inherits(resclv,"clv3w")) {
     appel      <- as.list(resclv$call)
     if(is.null(eval.parent(appel$K))) { 
       if (is.null(K)) {
         K<- as.numeric(readline("Please, give the number of groups : "))
       }
       if ((K > resclv$param$gmax+1) | (K<1))
         stop(paste("Consolidation not performed for ",K,"clusters"))
       comp     <- resclv[[K]]$comp
     } else {
       K<-eval.parent(appel$K)
       comp     <- resclv$comp
     }   
     if (graph) {
        dev.new()
        barplot(t(comp), ylab = 'Scores',beside=TRUE,names.arg=rownames(comp),cex.names=cex.lab,legend.text=paste('Cluster',1:K),las=2)
        title( main = 'Scores associated with the CLV3W partition ' )
     }
     return(comp=comp)
   }
   
   
   if(inherits(resclv,"lclv")) {
     if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
     compt<-resclv[[K]]$compt
     compc<-resclv[[K]]$compc
     return(list(compt=compt,compc=compc))
   }
     
     
    
      
  }