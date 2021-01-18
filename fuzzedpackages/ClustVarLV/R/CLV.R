#' Hierarchical clustering of variables with consolidation
#'
#' Hierarchical Cluster Analysis of a set of variables with consolidation.
#' Directional or local groups may be defined. Each group of variables is associated with a latent component.
#' Moreover, the latent component may be constrained using external information collected on the observations or on the variables.
#' 
#' If external variables are used, define either Xr or Xu, but not both.
#' Use the LCLV function when Xr and Xu are simultaneously provided.
#' 
#' @param X : The matrix of variables to be clustered
#' @param Xu : The external variables associated with the columns of X
#' @param Xr : The external variables associated with the rows of X
#' @param method : The criterion to be use in the cluster analysis.  \cr
#'        1 or "directional" : the squared covariance is used as a measure of proximity (directional groups). \cr    
#'        2 or "local"       : the covariance is used as a measure of proximity (local groups)
#' @param sX ,TRUE/FALSE : standardization or not of the columns X (TRUE by default)\cr
#'        (predefined -> cX = TRUE : column-centering of X)
#' @param sXr ,TRUE/FALSE : standardization or not of the columns Xr (FALSE by default)\cr
#'        (predefined -> cXr    = TRUE : column-centering of Xr)
#' @param sXu ,TRUE/FALSE : standardization or not of the columns Xu (FALSE by default)\cr
#'        (predefined -> cXu= FALSE : no centering, Xu considered as a weight matrix)
#' @param nmax : maximum number of partitions for which the consolidation will be done (by default nmax=20)
#' @param maxiter : maximum number of iterations allowed for the consolidation/partitioning algorithm (by default maxiter=20)
#' 
#' 
#' @return \item{tabres}{ Results of the clustering algorithm.
#'         In each line you find the results of one specific step of the hierarchical clustering.
#'         \itemize{
#'                \item {Columns 1 and 2}{ : The numbers of the two groups which are merged}
#'                \item {Column 3}{ : Name of the new cluster}
#'                \item {Column 4}{ : The value of the aggregation criterion for the Hierarchical Ascendant Clustering (HAC)}
#'                \item {Column 5}{ : The value of the clustering criterion for the HAC }
#'                \item {Column 6}{ : The percentage of the explained initial criterion value\cr
#'                (method 1 => \% var. expl. by the latent comp.)}
#'                \item {Column 7}{ : The value of the clustering criterion after consolidation}
#'                \item {Column 8}{ : The percentage of the explained initial criterion value after consolidation}        
#'                \item {Column 9}{ : The number of iterations in the partitioning algorithm. \cr
#'                Remark : A zero in columns 7 to 9 indicates that no consolidation was done}
#'        }}                
#' @return \item{partition K}{ contains a list for each number of clusters of the partition, K=2 to nmax with
#'          \itemize{
#'                \item {clusters}{ :  in line 1, the groups membership before consolidation; in line 2 the groups membership after consolidation} 
#'                \item {comp}{ : The latent components of the clusters (after consolidation)}
#'                \item {loading}{ : if there are external variables Xr or Xu :  The loadings of the external variables (after consolidation)}
#'          }}
#' @seealso CLV_kmeans, LCLV
#' 
#' @references Vigneau E., Qannari E.M. (2003). Clustering of variables around latents components. Comm. Stat, 32(4), 1131-1150.
#' @references Vigneau E., Chen M., Qannari E.M. (2015). ClustVarLV:  An R Package for the clustering of Variables around Latent Variables. The R Journal, 7(2), 134-148
#' 
#' @examples data(apples_sh)
#' #directional groups
#' resclvX <- CLV(X = apples_sh$senso, method = "directional", sX = TRUE)
#' plot(resclvX,type="dendrogram")
#' plot(resclvX,type="delta")
#' #local groups with external variables Xr
#' resclvYX <- CLV(X = apples_sh$pref, Xr = apples_sh$senso, method = "local", sX = FALSE, sXr = TRUE)
#'
#' @export                
#'                 
CLV <- function(X,Xu=NULL,Xr=NULL,method=NULL,sX=TRUE,sXr=FALSE,sXu=FALSE,nmax=20,maxiter=20)
{
  
 if(is.null(method)) stop('parameter method should be =1/"directional" or =2/"local"')
 if (method=="directional") method=1
 if (method=="local") method=2
 
 
 cX=TRUE
 cXr=TRUE
 cXu=FALSE
 

 # verification if some variables have constant values (standard deviation=0)
 who<-which(apply(X,2,sd,na.rm=TRUE)==0)
 if ((length(who)>0)&(sX==TRUE)) {
   listwho<-c(": ")
   for (r in 1:length(who)) {listwho=paste(listwho,colnames(X)[who[r]],",")}
   stop("The variables",listwho," have constant values (standard deviation=0). Please remove these variables from the X matrix.")
 }
 if (length(who)>0) {
   listwho<-c(": ")
   for (r in 1:length(who)) {listwho=paste(listwho,colnames(X)[who[r]],",")}
   warning("The variables",listwho," have constant values (standard deviation=0). Please remove these variables from the X matrix.")
 }
 

 X<- scale(X, center=cX, scale=sX)
 p <- ncol(X) 
 n <- nrow(X)
 
 valmq=FALSE
 # verification if there are NA values
 if (sum(is.na(X))>0)  {
   valmq=TRUE
   tauxNA=sum(is.na(X))/(n*p)
 }
 epsil=0.00001
 if(nmax>p) nmax<-p
 
 if (is.null(Xr)) {
   EXTr<-0
   } else {
      EXTr<-1                                    
      Xr<- scale(Xr, center=cXr, scale=sXr)
      ntilde <- dim(Xr)[1]
      q<-dim(Xr)[2] 
      if (n != ntilde) {stop("X and Xr must have the same number of observations")}
   } 
 if (is.null(Xu)) {
   EXTu<-0
   } else {
   EXTu<-1                   
   Xu<- scale(Xu, center=cXu, scale=sXu) 
   ptilde <- dim(Xu)[1]  
   m<-dim(Xu)[2] 
   if (p != ptilde) {stop("The number of consumers in X and Xu must are not the same") }                      
   if (EXTr==1) {stop("This procedure doesn't allow Xr and Xu to be defined simultaneously. Use LCLV instead")}
   }
 if (valmq & ((EXTr==1)|(EXTu==1))) stop("The matrix X contains missing values. Use a X matrix without missing value for CLV with external data")
 
 ###################################################################################
 # Hierarchical Ascendant algorithm
 ###################################################################################
 # first step : one variable=one cluster
 groupes <- 1:p                     
 fusions <- -groupes
 crit<-crit_init(method,X,EXTr,Xr,EXTu,Xu)
 inertie <- sum(crit)
 sbegin <- sum(crit) 
 ncluster <- p
 #print(paste('initial value of the criterion : ',round(inertie,2)))
 results = matrix(0,p-1,9)
 resultscc <- list()
 critstep <-  matrix(nrow=p, ncol=p)
 deltamin<-matrix(nrow=p,ncol=p)
 hmerge <-  matrix(0,nrow=p-1, ncol=2)
 delta <- matrix(0,nrow=p-1, ncol=1)
 ordr <- c()
 

 if (method == 1 & EXTu == 0) {
   if((EXTr==0)) { 
     if (!valmq)  {
       matcov = t(X)%*%X /(n-1)
       rescpp= critcpp(matcov,crit)   
       critstep = rescpp[[1]]
       deltamin = rescpp[[2]]
       critstep[lower.tri(critstep, diag = T)] = NA
       deltamin[lower.tri(critstep, diag = T)] = NA
     }
     if (valmq)   {
#     matcov = cov(X,use="pairwise.complete.obs")
       pg <- max(groupes)
       for (i in 1:(pg - 1)) {
         for (j in (i + 1):pg) {
           X12=X[,c(i,j)]
           if (EXTu==1) {Xu12<-Xu[c(i,j),]}
           critval<-crit_h(method,X12,EXTr,Xr,EXTu,Xu12,tauxNA=tauxNA)
           critstep[i,j]<-critval
           deltamin[i,j]<-crit[i]+crit[j]-critstep[i,j]
         }
       }
     }
   }
   if((EXTr==1)) {
     matcov = t(X)%*%Xr%*%t(Xr)%*%X/(n-1)
     rescpp= critcpp(matcov,crit)   
     critstep = rescpp[[1]]
     deltamin = rescpp[[2]]
     critstep[lower.tri(critstep, diag = T)] = NA
     deltamin[lower.tri(critstep, diag = T)] = NA
   }
 } else {
   pg <- max(groupes)
   for (i in 1:(pg - 1)) {
     for (j in (i + 1):pg) {
       X12=X[,c(i,j)]
       if (EXTu==1) {Xu12<-Xu[c(i,j),]}
       critval<-crit_h(method,X12,EXTr,Xr,EXTu,Xu12,tauxNA=tauxNA)
       critstep[i,j]<-critval
       deltamin[i,j]<-crit[i]+crit[j]-critstep[i,j]
     }
   }  
 }

 

# next steps
 for (level in 1:(p-1)) { 
   cmerge<-mincpp(deltamin)
   cmerge1 = cmerge[[1]]
   cmerge2 =  cmerge[[2]]
 		ind1 <- which(groupes == cmerge1)
    ind2 <- which(groupes == cmerge2)
    ind<-c(fusions[ind1[1]],fusions[ind2[1]])
    hmerge[level,]<-ind
    listind<-c(ind1,ind2)
    ordr <- order_var(ordr,listind)
    fusions[listind]<- level
    crit[listind[1]]<- critstep[cmerge1,cmerge2];
    autre <- listind[2:length(listind)]
    crit[autre]<-NA
    inertie <- rbind(inertie, (inertie[level] - deltamin[cmerge1,cmerge2]))
	  delta[level]<-deltamin[cmerge1,cmerge2]
	  inertie_inv<- -inertie+sbegin
	  ncluster <- ncluster-1	                     
    results[level,1:6]<- c(ind,level,delta[level],inertie[level+1],100*inertie[level+1]/sbegin)     
    groupes[groupes==cmerge2]<-cmerge1
    groupes[groupes>cmerge2]<-groupes[groupes>cmerge2]-1
    deltamin<-deltamin[-cmerge2,-cmerge2]
    critstep<-critstep[-cmerge2,-cmerge2]
    gr2<-which(groupes==cmerge1)
   
   itercm = c(1:max(groupes))[-cmerge1]
   gr = lapply(itercm,function(x) unique(c(which(groupes == x),gr2)))
   if(method==1 & EXTr==0 & EXTu==0) {
                  critval = lapply(gr,function(x) critval = crit_h(method,X[,x],EXTr,Xr,EXTu,Xu,tauxNA=tauxNA))
   } else {
     if (EXTu==1) critval = lapply(gr,function(x) critval = crit_h(method,X[,x],EXTr,Xr,EXTu,Xu[x,],tauxNA=tauxNA))
     else         critval = lapply(gr,function(x) critval = crit_h(method,X[,x],EXTr,Xr,EXTu,Xu,tauxNA=tauxNA))
   } 
   critgr = lapply(gr,function(x) crit[x[1]])
   
   deltaval = unlist(critgr)+crit[gr2[1]]-unlist(critval)
   
   if (cmerge1 != 1){
     critstep[1:(cmerge1-1),cmerge1]<-unlist(critval)[1:(cmerge1-1)]  
     deltamin[1:(cmerge1-1),cmerge1]<-deltaval[1:(cmerge1-1)]
   }
   if (cmerge1 != max(groupes)){
     critstep[cmerge1,(cmerge1+1):max(groupes)] = unlist(critval)[(cmerge1):(max(groupes)-1)]
     deltamin[cmerge1,(cmerge1+1):max(groupes)] = deltaval[(cmerge1):(max(groupes)-1)]
   }

    #---------------------------------------------
    # consolidation phase
    if ((ncluster <= nmax) & (ncluster > 1) ) {
    cc_consol <- t(t(groupes))
    K <- ncluster
    T<-c()   
    pcritav=0
                 
   for (i in 1:maxiter) {       
        critere <-rep(0,K)
        groupes_tmp <- cc_consol[,i]
        out<-mat_init(X,EXTr,Xr,EXTu,Xu,K)
        comp<-out$comp
        if (EXTr==1)  a<-out$a
        if (EXTu==1)  u<-out$u
       
        for (k in 1:K) {
            ind <- which(groupes_tmp == k)
            if (length(ind) > 0) {
              res<-consol_calcul(method,X,EXTr,Xr,EXTu,Xu,ind)
              critere[k]<-res$critere
              comp[,k]<-res$comp
              if (EXTr==1)  a[,k]<-res$a
              if (EXTu==1)  u[,k]<-res$u
           }
        }    
       pcrit<-sum(critere)/sbegin
 
       groupes_tmp<-consol_affect(method,X,Xr,Xu,EXTr,EXTu,comp,a,u)
       if(sum(is.na(groupes_tmp))>0) warning("a variable has not been allocated to any cluster at step ",i)
 
       if (length(which((cc_consol[, i] == groupes_tmp) == FALSE, arr.ind = TRUE)) == 0)  break
       if((pcrit-pcritav)<epsil) break
       cc_consol = cbind(cc_consol, groupes_tmp)
       pcritav=pcrit
    }
    rownames(cc_consol) <- colnames(X)      
    names(cc_consol) = NULL
                                    
    initgroupes<-cc_consol[,1]
    lastgroupes<-cc_consol[,ncol(cc_consol)]
    if ((EXTu==0)&(EXTr==0)) listcc = list(clusters = rbind(initgroupes,lastgroupes),  comp=comp)
    else if ((EXTu==0)&(EXTr==1)) listcc = list(clusters = rbind(initgroupes,lastgroupes),  comp=comp,loading=a)
    else if ((EXTu==1)&(EXTr==0)) listcc = list(clusters = rbind(initgroupes,lastgroupes),  comp=comp,loading=u)
    results[level,7:9]<- c(sum(critere),(sum(critere)/sbegin)*100, i)  
    resultscc[[K]] <- listcc 
              	 	                                   
    }   # --------------------  end of the consolidation 
 } # end of the current hierarchical step



# if there is only one cluster (K=1)
  results[p-1,7:9]<- c(results[p-1,5],results[p-1,6], 0)  
  group1<-matrix(1,nrow=2,ncol=p)
  colnames(group1) <- colnames(X)
  rownames(group1) <- c("initgroupes","lastgroupes")
  out<-mat_init(X,EXTr,Xr,EXTu,Xu,1)
  comp<-out$comp
  if (EXTr==1)  a<-out$a
  if (EXTu==1)  u<-out$u
  ind<-c(1:p)
  res<-consol_calcul(method,X,EXTr,Xr,EXTu,Xu,ind) 
   comp[,1]<-res$comp
  if (EXTr==1) a=res$a
  if (EXTu==1) u=res$u       
  if ((EXTu==0)&(EXTr==0)) listcc = list(clusters = group1,  comp=comp)
  if ((EXTu==0)&(EXTr==1)) listcc = list(clusters = group1,  comp=comp,loading=a)  
  if ((EXTu==1)&(EXTr==0)) listcc = list(clusters = group1,  comp=comp,loading=u)      
  resultscc[[1]] <- listcc 
     
  # consolidation result when K=p (only required for homogeneity of the results)
  if (nmax == p) {
      groupes<-matrix(1:p,nrow=p,ncol=1)
      rownames(groupes) <- colnames(X)
      critere <-rep(0,p)
      out<-mat_init(X,EXTr,Xr,EXTu,Xu,p)
      comp<-out$comp
      if (EXTr==1)  a<-out$a
      if (EXTu==1)  u<-out$u
      for (k in 1:p) {
          ind <- which(groupes == k)
          res<-consol_calcul(method,X,EXTr,Xr,EXTu,Xu,ind)
          critere[k]<-res$critere
          comp[,k]<-res$comp
          if (EXTr==1)  a[,k]<-res$a
          if (EXTu==1)  u[,k]<-res$u
      }    
      groupes = cbind(groupes,groupes)  
      colnames(groupes) <- c("initgroupes","lastgroupes")
      if ((EXTu==0)&(EXTr==0)) listcc = list(clusters = t(groupes),  comp=comp)
      if ((EXTu==0)&(EXTr==1)) listcc = list(clusters = t(groupes),  comp=comp,loading=a) 
      if ((EXTu==1)&(EXTr==0)) listcc = list(clusters = t(groupes),  comp=comp,loading=u)      
      resultscc[[p]] <- listcc      
   }  
 
                          
 colnames(results)= c("merg1","merg2","new.clust","agg.crit.hac","clust.crit.hac",
                      "%S0expl.hac","clust.crit.cc","%S0expl.cc","iter")
 names(resultscc) = paste("partition",1:min(p,nmax),sep="")
 resultscc$tabres=results
 resultscc$param<-list(X=X,method=method,n = n, p = p,nmax = nmax,EXTu=EXTu,EXTr=EXTr,sX=sX,sXr=sXr,cXu=cXu,sXu=sXu,sbegin=sbegin,strategy="none")
 resultscah=list(labels=colnames(X),inertie=inertie, height=delta, merge=hmerge,order=ordr)    
 mytot<-resultscah  
 class(mytot)="hclust"
# clvclt= c(resultscc, list(mydendC = mytot)) 
 mydendC=as.dendrogram(mytot)
 clvclt= c(resultscc, list(mydendC = mydendC))  

 
 class(clvclt) = "clv"
 return(clvclt) 
}
