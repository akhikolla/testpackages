#'Deviance
#'
#'Calculate deviances at provided coefficients
#'
#'@param X Input matrix
#'@param z Response vector
#'@param py1 True prevalence Pr(Y=1)
#'@param coefMat A coefficient matrix whose column corresponds to a set of coefficients
#'@param weights observation weights. Default is 1 for each observation.
#'@return deviances
#'@examples
#'data("simulPU")
#'coef0<-replicate(2,runif(ncol(simulPU$X)+1))
#'deviances(simulPU$X,simulPU$z,py1=simulPU$truePY1,coefMat = coef0)
#'@importFrom Rcpp evalCpp
#'@importFrom methods as
#'@useDynLib PUlasso
#'@export
#'
deviances <-function(X,z,py1,coefMat,weights = NULL)
{
  if(is.null(dim(X))){stop("not a valid X")}
  if(is.null(colnames(X))){colnames(X) <- paste("V",1:ncol(X),sep = "")}
  
  row_ordering= order(z,decreasing = T); col_ordering = 1:ncol(X); group = 1:ncol(X);
  ordering_res = ordering_data(row_ordering,col_ordering, X, z, group, weights)
  X_lu = ordering_res$X_lu; z_lu = ordering_res$z_lu; w_lu = ordering_res$w_lu;
  remove(X,z,ordering_res)
  
  # Normalize weights
  if(!is.null(w_lu)){weiOption<- TRUE; w_lu <- w_lu/sum(w_lu)*length(w_lu)}else{
    weiOption <- FALSE; w_lu <- rep(1,nrow(X_lu))}
  
  is.sparse = FALSE
  if (inherits(X_lu, "sparseMatrix")) {
    is.sparse = TRUE
    X_lu = as(X_lu, "CsparseMatrix")
    X_lu = as(X_lu, "dgCMatrix")
  } else if (inherits(X_lu,"dgeMatrix")){
    X_lu = as.matrix(X_lu)
  }
  if(!(class(X_lu)=="matrix"||class(X_lu)=="dgCMatrix")){stop("X must be a matrix, or a sparse matrix")}
  if(typeof(coefMat)=="double"){coefMat <- as.matrix(coefMat)}
  if(nrow(coefMat)!=(ncol(X_lu)+1)){stop("nrow(coefMat) must be the same as p+1")}
  
  
  if(!is.sparse){
    dev<- deviances_dense_cpp(X_ = X_lu,z_ = z_lu,pi_ = py1,coefMat_ = coefMat, wei_ = w_lu, weiOption_ = weiOption)
  }else{
    dev<- deviances_sparse_cpp(X_ = X_lu,z_ = z_lu,pi_ = py1,coefMat_ = coefMat,wei_ = w_lu,  weiOption_= weiOption)
  }
  return(c(dev))
  
}

