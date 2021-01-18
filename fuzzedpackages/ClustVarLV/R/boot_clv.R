#'@title Boostrapping for assessing the stability of a CLV result
#'
#'@description
#'Bootstrapping on the samples is performed. 
#'Each boostrapped data matrix is submitted to CLV in order to get partitions from 1 to nmax clusters.
#'For each number of clusters, K, the adjusted Rand Index between actual and the bootstrapped partitions are computed and used in order to assess the stability of the solution into K clusters.
#'Parallel computing is performed in order to save time.
#' 
#' @param object : result of CLV()
#' @param B : the number of bootstrap to be run (100 by default)
#' @param nmax : maximal size of the partitions to be considered (if NULL, the value of nmax used for the object is used)
#' 
#' @return \item{matARI}{a matrix of the Adjusted Rand Index of size (B x nmax).}
#'
#' @seealso CLV
#' 
#' @export


 boot_clv =  function (object, B=100,nmax=NULL) 
 {
  
    resclv<-object   
    if (!inherits(resclv, "clv"))   stop("non convenient objects")
    
    X<-resclv$param$X
    if(is.null(nmax)) nmax=resclv$param$nmax
    method<-resclv$param$method
    n<-resclv$param$n
    p<-resclv$param$p
    sX<-resclv$param$sX
  
    # boot_fx ------------------------------------------
    boot_fx <- function(b) {
      resari=c()
      bsamp=sample(1:n,replace=TRUE)
      if (p<=1500) {
        resboot=CLV(X[bsamp,],method=method,sX=sX,nmax=nmax)
        for (k in 2:nmax) {
          resari=c(resari,ARI(get_partition(resclv,k),get_partition(resboot,k)))
        }
      } else {
        for (k in 2:nmax) {
          resboot=CLV_kmeans(X[bsamp,],method=method,sX=sX,clust=k,nstart=50)
          resari=c(resari,ARI(get_partition(resclv,k),get_partition(resboot)))
        } 
      }
      #     
      return(resari)
    }  
    # ARI -------------------------------------------------
    ARI<-function(P1,P2) {
      # compute the Adjusted Rand Index between two partitions
      tab <- table(as.vector(P1), as.vector(P2))
      totm1=apply(tab,1,sum)
      totm2=apply(tab,2,sum)
      tot=sum(tab)
      index=sum(choose(tab,2))
      expect=sum(choose(totm1,2))*sum(choose(totm2,2))/sum(choose(tot,2))
      maxindex=(sum(choose(totm1,2))+sum(choose(totm2,2)))/2
      ari=(index-expect)/(maxindex-expect)
      return(ari)
    } 
  
  numCores <- parallel::detectCores()
  #  numCores
  cl <- parallel::makeCluster(numCores)
  doParallel::registerDoParallel(cl)
  b<-NULL
     resari<- foreach(b=1:B,.combine="rbind",.packages="ClustVarLV") %dopar%
      {
           boot_fx(b)
      }
   parallel::stopCluster(cl)   
  
   resari=cbind(rep(NA,B),resari)
   mari=apply(resari,2,mean)
   binf=apply(resari,2,mean)-apply(resari,2,sd)
   bsup=apply(resari,2,mean)+apply(resari,2,sd)
   maxi=ceiling(max(bsup,na.rm=TRUE)*10)/10
   mini=floor(min(binf,na.rm=TRUE)*10)/10
  
#   boxplot(resari,ylim=c(mini,maxi),xlab="nb clusters",ylab="Adjusted Rand Index")
    plot(1:nmax,apply(resari,2,mean),ylim=c(mini,maxi),type="b",xlab="nb clusters",ylab="Adjusted Rand Index (mean+/- 1 sd)")
    for (k in 1:nmax) arrows(x0=k,y0=binf[k],y1=bsup[k],code=0,lty=2)
   
   return(matARI=resari)  
 }