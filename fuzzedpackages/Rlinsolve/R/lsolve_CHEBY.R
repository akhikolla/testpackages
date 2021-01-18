#' Chebyshev Method
#'
#' Chebyshev method - also known as Chebyshev iteration - avoids computation of inner product,
#' enabling distributed-memory computation to be more efficient at the cost of requiring
#' a priori knowledge on the range of spectrum for matrix \code{A}.
#'
#' @param A an \eqn{(m\times n)} dense or sparse matrix. See also \code{\link[Matrix]{sparseMatrix}}.
#' @param B a vector of length \eqn{m} or an \eqn{(m\times k)} matrix (dense or sparse) for solving \eqn{k} systems simultaneously.
#' @param xinit a length-\eqn{n} vector for initial starting point. \code{NA} to start from a random initial point near 0.
#' @param reltol tolerance level for stopping iterations.
#' @param maxiter maximum number of iterations allowed.
#' @param preconditioner an \eqn{(n\times n)} preconditioning matrix; default is an identity matrix.
#' @param adjsym a logical; \code{TRUE} to symmetrize the system by transforming the system into normal equation, \code{FALSE} otherwise.
#' @param verbose a logical; \code{TRUE} to show progress of computation.
#'
#' @return a named list containing \describe{
#' \item{x}{solution; a vector of length \eqn{n} or a matrix of size \eqn{(n\times k)}.}
#' \item{iter}{the number of iterations required.}
#' \item{errors}{a vector of errors for stopping criterion.}
#' }
#'
#' @examples
#' ## Overdetermined System
#' set.seed(100)
#' A = matrix(rnorm(10*5),nrow=10)
#' x = rnorm(5)
#' b = A%*%x
#'
#' out1 = lsolve.sor(A,b,w=0.5)
#' out2 = lsolve.cheby(A,b)
#' matout = cbind(x, out1$x, out2$x);
#' colnames(matout) = c("original x","SOR result", "Chebyshev result")
#' print(matout)
#'
#' @references
#' \insertRef{gutknecht_chebyshev_2002}{Rlinsolve}
#'
#' @rdname krylov_CHEBY
#' @export
lsolve.cheby <- function(A,B,xinit=NA,reltol=1e-5,maxiter=10000,
                      preconditioner=diag(ncol(A)),adjsym=TRUE,verbose=TRUE){
  ###########################################################################
  # Step 0. Initialization
  if (verbose){
    message("* lsolve.cheby : Initialiszed.")
  }
  if (any(is.na(A))||any(is.infinite(A))||any(is.na(B))||any(is.infinite(B))){
    stop("* lsolve.cheby : no NA or Inf values allowed.")
  }
  sparseformats = c("dgCMatrix","dtCMatrix","dsCMatrix")
  if (aux.is.sparse(A)||aux.is.sparse(B)||aux.is.sparse(preconditioner)){
    A = Matrix(A,sparse=TRUE)
    B = Matrix(B,sparse=TRUE)
    preconditioner = Matrix(preconditioner,sparse=TRUE)
    sparseflag = TRUE
  } else {
    A = matrix(A,nrow=nrow(A))
    if (is.vector(B)){
      B = matrix(B)
    } else {
      B = matrix(B,nrow=nrow(B))
    }
    preconditioner = matrix(preconditioner,nrow=nrow(preconditioner))
    sparseflag = FALSE
  }
  # xinit
  if (length(xinit)==1){
    if (is.na(xinit)){
      xinit = matrix(rnorm(ncol(A)))
    } else {
      stop("* lsolve.cheby : please use a valid 'xinit'.")
    }
  } else {
    if (length(xinit)!=ncol(A)){
      stop("* lsolve.cheby : 'xinit' has invalid size.")
    }
    xinit = matrix(xinit)
  }
  ###########################################################################
  # Step 1. Preprocessing
  # 1-1. Neither NA nor Inf allowed.
  if (any(is.infinite(A))||any(is.na(A))||any(is.infinite(B))||any(is.na(B))){
    stop("* lsolve.cheby : no NA, Inf, -Inf values are allowed.")
  }
  # 1-2. Size Argument
  m = nrow(A)
  if (is.vector(B)){
    mB = length(B)
    if (m!=mB){
      stop("* lsolve.cheby : a vector B should have a length of nrow(A).")
    }
  } else {
    mB = nrow(B)
    if (m!=mB){
      stop("* lsolve.cheby : an input matrix B should have the same number of rows from A.")
    }
  }
  if (is.vector(B)){
    B = as.matrix(B)
  }
  # 1-3. Adjusting Case
  if (m > ncol(A)){        ## Case 1. Overdetermined
    B = t(A)%*%B
    A = t(A)%*%A
  } else if (m < ncol(A)){ ## Case 2. Underdetermined
    stop("* lsolve.cheby : underdetermined case is not supported.")
  } else {                 ## Case 3. Square Size
    if (norm(abs(t(A)-A),"f")>1e-10){
      if (verbose){
        message("* lsolve.cheby : A may not be symmetric.")
      }
      if (adjsym){
        B = t(A)%*%B
        A = t(A)%*%A
        if (verbose){
          message("* lsolve.cheby : making it normal equation form via 'adjsym' flag.")
        }
      }
    }
  }
  # 1-4. Preconditioner : only valid for square case
  if (!all.equal(dim(A),dim(preconditioner))){
    stop("* lsolve.cheby : Preconditioner is a size-matching.")
  }
  if (verbose){message("* lsolve.cheby : preprocessing finished ...")}
  ###########################################################################
  # Step 2. Main Computation
  eigs = solve(preconditioner, A)
  eigsdec = eigen(eigs, only.values = TRUE)
  eigmax  = max(eigsdec$values)
  eigmin  = min(eigsdec$values)

  ncolB = ncol(B)
  if (ncolB==1){
    if (!sparseflag){
      vecB = as.vector(B)
      res = linsolve.cheby.single(A,vecB,xinit,reltol,maxiter,preconditioner,eigmax,eigmin)
    } else {
      vecB = B
      res = linsolve.cheby.single.sparse(A,vecB,xinit,reltol,maxiter,preconditioner,eigmax,eigmin)
    }
  } else {
    x      = array(0,c(ncol(A),ncolB))
    iter   = array(0,c(1,ncolB))
    errors = list()
    for (i in 1:ncolB){
      if (!sparseflag){
        vecB = as.vector(B[,i])
        tmpres = linsolve.cheby.single(A,vecB,xinit,reltol,maxiter,preconditioner,eigmax,eigmin)
      } else {
        vecB = Matrix(B[,i],sparse=TRUE)
        tmpres = linsolve.cheby.single.sparse(A,vecB,xinit,reltol,maxiter,preconditioner,eigmax,eigmin)
      }
      x[,i]       = tmpres$x
      iter[i]     = tmpres$iter
      errors[[i]] = tmpres$errors
      if (verbose){
        message(paste("* lsolve.cheby : B's column.",i,"being processed.."))
      }
    }
    res = list("x"=x,"iter"=iter,"errors"=errors)
  }

  if ("flag"%in%names(res)){
    flagval = res$flag;
    if (flagval==0){
      if (verbose){
        message("* lsolve.cheby : convergence was well achieved.")
      }
    } else {
      if (verbose){
        message("* lsolve.cheby : convergence was not achieved within maxiter.")
      }
    }
    res$flag = NULL
  }
  return(res)
  if (verbose){
    message("* lsolve.cheby : computations finished.")
  }
}
