# weigths of the external variables in each cluster
# Applies only when external variables (Xr, Xu or both) are involved.
# 
get_weight.clv <-
  function(resclv, K=NULL,graph=FALSE,cex.lab=1)
  {
    
   if (!inherits(resclv, c("clv","lclv"))) 
      stop("non convenient objects")
    
 #  X<-resclv$param$X
      
   if(inherits(resclv,"clv")) {
    if(is.null(resclv$param$K)) { 
      if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
      if (is.null(resclv[[K]]$loading)) stop("applies only if external variables are taken into account")
      weight<-resclv[[K]]$loading
    } else {
      if (is.null(resclv$loading)) stop("applies only if external variables are taken into account")
      weight<-resclv$loading
      K<-resclv$param$K
    }
    if (graph) {
       dev.new()
       barplot(t(weight), ylab = 'Weights',beside=TRUE,names.arg=rownames(weight),cex.names=cex.lab,legend.text=paste('Cluster',1:K),las=2)
       title( main = 'Weights associated with the external variables' )
     }
    return(weight=weight)
   }
   
   if(inherits(resclv,"lclv")) {
     if (is.null(K)) {K<- as.numeric(readline("Please, give the number of groups : "))}
     weight_u<-resclv[[K]]$loading_u
     weight_v<-resclv[[K]]$loading_v
     if (graph) {
       dev.new()
       barplot(t(weight_u), ylab = 'Weights',beside=TRUE,names.arg=rownames(weight_u),cex.names=cex.lab,legend.text=paste('Cluster',1:K),las=2)
       title( main = 'Weights associated with the external Xu variables' )
       dev.new()
       barplot(t(weight_v), ylab = 'Weights',beside=TRUE,names.arg=rownames(weight_v),cex.names=cex.lab,legend.text=paste('Cluster',1:K),las=2)
       title( main = 'Weights associated with the external Xr variables' )
     }
     return(list(weight_u=weight_u,weight_v=weight_v))
   }
      
  }