#'@title  Imputation of a data matrix based on CLV results
#' 
#' @description
#' For each variable, its missing data will be imputed according to the values of the latent variable of the group in which the variable belong to.
#' 
#' @param x : an object of class \code{clv}
#' @param X0 : the initial data matrix with missing values (NA)
#' @param K : the number of Latent Variables to be considered, each of them being associated with a group of variables.
#'
#' @return X0imput : the imputed data matrix, in the original scale
#' @return Ximput  : the imputed matrix, centered and scaled according to the pretratment parameters chosen in CLV
#' 
#' @details It is adviced to use a larger number of latent variables, on the basis of which the imputation will be done, than the suspected 'true' number of groups of variables
#'  
#' @export
#' 
imput_clv =  function (x,X0,K=NULL) 
{
  if (!inherits(x, "clv")) 
    stop("non convenient object")
  
  # verification if there are NA values in X0
  valmq=FALSE
  if (sum(is.na(X0))>0)  {
    valmq=TRUE
  }
  if (!valmq) {stop("There is no missing (NA) values in the data matrix. Imputation is not required")}
 
  
  if(is.null(x$param$K)) { 
    if (is.null(K)) {K<- as.numeric(readline("Please, give the number of latent variables used for imputation"))}
    if(K>x$param$nmax) {stop("Please, perform again CLV with a value for the 'nmax' parameter greater than the number of LVs chosen for the imputation")}  
  } else {
    K<-x$param$K
  }
  
  X=  x$param$X
  p = x$param$p
  n = x$param$n
  method =  x$param$method  
  sX<-x$param$sX
  strategy<-x$param$strategy

  clust<-get_partition(x,K)
  comp<-get_comp(x,K)
  

  # verification if there are NA values in comp
  # and that there is no empty cluster among the K chosen
  if ((sum(is.na(comp))>0) |(length(table(clust))<K)) {
    if(is.null(x$param$K)) { 
        stop("At least one Latent Variable contains a missing values or K distincts Latent Variables can not be defined. \nTry the imputation with a lower value for parameter K")
    } else {
       stop("At least one Latent Variable contains a missing values or K distincts Latent Variables can not be defined. \nperform again the CLV_kmeans functions with a lower value for parameter K")
    }
  }

  if (length(table(clust))<K) {stop("K distincts Latent Variables can not be defined. Try the imputation with a lower value for parameter K.")}
  
  out<-NULL
  if (strategy=="kplusone") {
    warning("The variables assigned to the noise cluster are set apart. The corresponding columns in the data set will be set to NA.")
    out<-which(clust==0)
  }
  if (strategy=="sparselv") {
    warning("The variables with zero loadings regarding their goup's LV will be imputed as the others. Caution must be taken in term of interpretation.")
  }
  
  
  
  # 1st step : imputation of the data matrix (centered and scaled matrix)
  Ximput<-matrix(NA,nrow=n,ncol=p)
  for (k in 1:K) {
      Xk<-as.matrix(X[,which(clust==k)])
      sinit=apply(Xk,2,sd,na.rm=TRUE)
      lXk<-as.list(data.frame(Xk))
      if(method==1) compl<-lapply(X=lXk,FUN=compl1j, y=comp[,k])
      if(method==2) compl<-lapply(X=lXk,FUN=compl2j, y=comp[,k])
      Xk<-do.call(cbind,compl)
      snow<-apply(Xk,2,sd,na.rm=TRUE)
      Xk<-scale(Xk,center=TRUE,scale=snow/sinit)
      Ximput[,which(clust==k)]<-Xk
   }
  
  
  # 2nd step : reconstruction of the missing values in X0 with initial means and sd)
  m=apply(X0,2,mean,na.rm=TRUE)
  s=apply(X0,2,sd,na.rm=TRUE)
  snow=apply(Ximput,2,sd)
  Ximputback<-scale(Ximput,center=FALSE,scale=snow/s)
  Ximputback<-scale(Ximputback,center=-m,scale=FALSE)
  Ximputback<-as.data.frame(Ximputback)
  indic<-is.na(X0)
  XX0<-X0
  XX0[(is.na(X0))]<-0
  X0imput<-((!indic)*XX0)+ (indic*(Ximputback))
 
  rownames(X0imput)=rownames(X0)
  rownames(Ximput)=rownames(X0)
  colnames(X0imput)=colnames(X0)
  colnames(Ximput)=colnames(X0)
  
  if (length(out)>0) {
    Ximput[,out]<-matrix(NA,nrow=n,ncol=length(out))
    X0imput[,out]<-matrix(NA,nrow=n,ncol=length(out))
    
  }
  
  
  return(list(X0imput=X0imput,Ximput=Ximput))
  
}

compl1j<-function(x,y) {x[which(is.na(x))]<-sign(cor(x,y,use="pairwise.complete.obs"))*y[which(is.na(x))]; return(x)}

compl2j<-function(x,y) {x[which(is.na(x))]<-y[which(is.na(x))]; return(x)}