#' @title Loadings of the variables on the latent component, in each cluster.
#' 
#' @description To get the variables loadings for the latent component in each cluster. \cr
#' For CLV(_kmeans), the loadings are of particular interest when method="directional" or when strategy="sparse LV" \cr
#' For CLV3W(_kmeans), the loadings are given for the variables associated with mode 2 of the 3-way array.
#' 
#' @param resclv : result of CLV(), CLV_kmeans(), LCLV(), CLV3W() or CLV3W_kmeans()
#' @param K : the number of clusters chosen (already defined if CLV_kmeans or CLV3W_kmeans is used)
#' @param type : outputs in the form of a "list" (one element by cluster, by default) or a "vector" (available only if "sparselv" strategy is used)
#' @param graph : boolean, if TRUE, the barplot associated with the scores is displayed (default : graph=FALSE)
#' @param cex.lab : magnification to be used for labels (1 by default)
#' 
#' @return \item{loading}{the loadings of the variables on each cluster's latent component}
#'       
#' @export
#' 
get_loading <-
  function(resclv, K=NULL, type="list",graph=FALSE, cex.lab=1)
  {
    if (!inherits(resclv, c("clv","clv3w")))  stop("non convenient objects")
    
    # for object clv
    if(inherits(resclv,"clv")) {
      if (is.null(resclv$param$K)) {              
        if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
      }
      
      if (is.null(resclv$param$nmax)&(resclv$param$strategy=="sparselv"))   {
        # if CLV_kmeans and strategy="sparselv"
        K=resclv$param$K
        loading=resclv$sloading
        if (type=="list"){
          if (graph) {
            dev.new()
            varnam<-colnames(resclv$param$X)
            if (length(varnam)<10) {
              # create a matrix p var x K groups od the loadings with NA if the var does nor belong to the group
              matload=matrix(NA,length(varnam),K)
              for (k in 1:K) {
                loadk=loading[[k]]
                namk=dimnames(loadk)[[1]]
                for (j in 1:length(namk)) {matload[which(varnam==namk[j]),k]<-loadk[j]}
              }
              barplot(t(matload), ylab = 'Loadings',beside=TRUE,names.arg=varnam,cex.names=cex.lab,
                      las=2, col=1:K, legend.text=paste('Cluster',1:K),main= "Loadings associated with the CLV clusters")
            } else {
              p=ncol(resclv$param$X)
              if ((p/K)<20) {
                NC=3
              } else {
                if ((p/K)<50){
                  NC=2
                } else {
                  NC=1
                }}
              NL=ceiling(K/NC)
              par(mfrow=c(NL,NC))
              for(k in 1:K) {
                namk=dimnames(loading[[k]])[[1]]
                barplot(loading[[k]],beside=TRUE, ylab = 'Loadings',names.arg=namk,cex.names=cex.lab,
                        las=2, col=k, main= paste("Loadings associated with the cluster",k))
              }
            }
          }
          return(loading=resclv$sloading)
        }
        if (type=="vector") {
            varnam<-colnames(resclv$clusters)
            K=resclv$param$K
            tab=rep(0,length(varnam))
            for (k in 1:K) {
              sploadk=resclv$sloading[[k]]
              namk=dimnames(sploadk)[[1]]
              for (j in 1:length(namk)) {tab[which(varnam==namk[j])]<-sploadk[j]}
            }
            return(loading=tab)
        }
       } else {   
       # other cases with CLV or CLV_kmeans
         if (type=="vector") { 
           print("available only if 'sparselv' strategy is used. 'list' is used instead.")
           type="list"
         }
         loading=list()
         X<-resclv$param$X
         method=resclv$param$method
         if (!is.null(resclv$param$K)) {
           compnorm=scale(resclv$comp)
           clusters<-as.vector(resclv$clusters[2,])
           K=resclv$param$K
         } else {
           compnorm=scale(resclv[[K]]$comp)
           clusters<-as.vector(resclv[[K]]$clusters[2,])
         }
         if (method==1) {
         for (k in 1:K) {
           Xk=X[,which(clusters==k)]
           valmq=FALSE
           if (sum(is.na(Xk))>0)  valmq=TRUE
           if(!valmq) cova=cov(Xk,compnorm[,k])
           if(valmq) cova=t(covamiss(Xk,compnorm[,k],method))
           loading[[k]]=cova/sqrt(sum(cova^2))
           rownames(loading[[k]])=colnames(Xk)
          }
         }
         if (method==2) {
           for (k in 1:K) {
             Xk=X[,which(clusters==k)]
             pk=sum(clusters==k)
             loading[[k]]=as.matrix(rep(1/pk,pk))
             rownames(loading[[k]])=colnames(Xk)
           }
         }
         if (graph) {
           dev.new()
           varnam<-colnames(X)
           if (length(varnam)<10) {
             # create a matrix p var x K groups od the loadings with NA if the var does nor belong to the group
             matload=matrix(NA,length(varnam),K)
             for (k in 1:K) {
               loadk=loading[[k]]
               namk=dimnames(loadk)[[1]]
               for (j in 1:length(namk)) {matload[which(varnam==namk[j]),k]<-loadk[j]}
             }
             barplot(t(matload), ylab = 'Loadings',beside=TRUE,names.arg=varnam,cex.names=cex.lab,
                   las=2, col=1:K, legend.text=paste('Cluster',1:K),main= "Loadings associated with the CLV clusters")
           } else {
             p=ncol(X)
             if ((p/K)<20) {
               NC=3
             } else {
             if ((p/K)<50){
                NC=2
             } else {
               NC=1
             }}
             NL=ceiling(K/NC)
             par(mfrow=c(NL,NC))
             for(k in 1:K) {
               namk=dimnames(loading[[k]])[[1]]
               barplot(loading[[k]],beside=TRUE, ylab = 'Loadings',names.arg=namk,cex.names=cex.lab,
                       las=2, col=k, main= paste("Loadings associated with the cluster",k))
             }
           }
         }
         return(loading=loading)
       }
    }
    
    # for object clv3w 
    if(inherits(resclv,"clv3w")) {
      appel      <- as.list(resclv$call)
      if(is.null(eval.parent(appel$K))) { 
        if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
        if ((K > resclv$param$gmax+1) | (K<1))
          stop(paste("Consolidation not performed for ",K,"clusters"))
        loading     <- resclv[[K]]$loading
        ########################################mise sous forme d'une liste ???????
      } else {
        K<-eval.parent(appel$K)
        loading     <- resclv$loading
        ########################################mise sous forme d'une liste ???????
      }   
      if (graph) {
        dev.new()
        barplot(t(loading), ylab = 'Loadings',beside=TRUE,names.arg=rownames(loading),cex.names=cex.lab,legend.text=paste('Cluster',k),las=2)
        title( main = 'Loadings associated with the CLV3W partition ' )
      }
      return(loading=loading)
    }
  }