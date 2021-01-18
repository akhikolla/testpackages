# sparse loadings in each cluster when using the "sparselv" strategy
# Applies only on CLV_kmeans output with strategy="sparselv".
# resclv : result of CLV_kmeans() 
# type : presented in the form of a "list" (one element by cluster, by default) or a "vector"
#  
get_sparseload <-
  function(resclv,type="list")
  {
    
   if (!inherits(resclv, c("clv")))   stop("non convenient objects")
   if (!(is.null(resclv$param$nmax)&(resclv$param$strategy=="sparselv")))   
                   stop("only for output of strategy \"sparselv\" from CLV-kmeans")
  
      
   if(resclv$param$strategy=="sparselv") {
     if (type=="list") { return(sparse_loadings=resclv$sloading)}
     if (type=="vector") {
          varnam<-colnames(resclv$clusters)
          K=resclv$param$K
          tab=rep(0,length(varnam))
          for (k in 1:K) {
               sploadk=resclv$sloading[[k]]
               namk=dimnames(sploadk)[[1]]
               for (j in 1:length(namk)) {tab[which(varnam==namk[j])]<-sploadk[j]}
          }
          return(sparse_loadings=tab)
      }
   }
    
    
  }