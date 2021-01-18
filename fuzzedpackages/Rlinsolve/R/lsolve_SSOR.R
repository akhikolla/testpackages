#' Symmetric Successive Over-Relaxation method
#'
#' Symmetric Successive Over-Relaxation(SSOR) method is a variant of Gauss-Seidel method for solving a system of linear equations,
#' with a decomposition \eqn{A = D+L+U} where \eqn{D} is a diagonal matrix and
#' \eqn{L} and {U} are strictly lower/upper triangular matrix respectively.
#' For a square matrix \eqn{A}, it is required to be diagonally dominant or symmetric and positive definite like GS method.
#' For an overdetermined system where \code{nrow(A)>ncol(A)},
#' it is automatically transformed to the normal equation. Underdetermined system -
#' \code{nrow(A)<ncol(A)} - is not supported.
#'
#' @param A an \eqn{(m\times n)} dense or sparse matrix. See also \code{\link[Matrix]{sparseMatrix}}.
#' @param B a vector of length \eqn{m} or an \eqn{(m\times k)} matrix (dense or sparse) for solving \eqn{k} systems simultaneously.
#' @param xinit a length-\eqn{n} vector for initial starting point. \code{NA} to start from a random initial point near 0.
#' @param reltol tolerance level for stopping iterations.
#' @param maxiter maximum number of iterations allowed.
#' @param w a weight value in \eqn{(0,2).}; \code{w=1} leads to Gauss-Seidel method.
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
#' out1 = lsolve.ssor(A,b)
#' out2 = lsolve.ssor(A,b,w=0.5)
#' out3 = lsolve.ssor(A,b,w=1.5)
#' matout = cbind(matrix(x),out1$x, out2$x, out3$x);
#' colnames(matout) = c("true x","SSOR w=1", "SSOR w=0.5", "SSOR w=1.5")
#' print(matout)
#'
#' @references
#' \insertRef{demmel_applied_1997}{Rlinsolve}
#' @rdname basic_SSOR
#' @export
lsolve.ssor <- function(A,B,xinit=NA,reltol=1e-5,maxiter=1000,w=1,adjsym=TRUE,verbose=TRUE){
  if (verbose){
    message("* lsolve.ssor : Initialiszed.")
  }
  if (any(is.na(A))||any(is.infinite(A))||any(is.na(B))||any(is.infinite(B))){
    stop("* lsolve.ssor : no NA or Inf values allowed.")
  }
  # Preprocessing : sparsity
  # http://dirk.eddelbuettel.com/tmp/RcppArmadillo-sparseMatrix.pdf
  sparseformats = c("dgCMatrix","dtCMatrix","dsCMatrix")
  if (aux.is.sparse(A)||aux.is.sparse(B)){
    A = Matrix(A,sparse=TRUE)
    B = Matrix(B,sparse=TRUE)
    sparseflag = TRUE
  } else {
    A = matrix(A,nrow=nrow(A))
    if (is.vector(B)){
      B = matrix(B)
    } else {
      B = matrix(B,nrow=nrow(B))
    }
    sparseflag = FALSE
  }
  # xinit
  if (length(xinit)==1){
    if (is.na(xinit)){
      xinit = matrix(rnorm(ncol(A)))
    } else {
      stop("* lsolve.ssor : please use a valid 'xinit'.")
    }
  } else {
    if (length(xinit)!=ncol(A)){
      stop("* lsolve.ssor : 'xinit' has invalid size.")
    }
    xinit = matrix(xinit)
  }
  # Preprocessing : symmetricity warning
  if (nrow(A)==ncol(A)){
    if (norm(abs(t(A)-A),"f")>1e-10){
      if (verbose){
        message("* lsolve.ssor : A may not be symmetric.")
      }
      if (adjsym){
        B = t(A)%*%B
        A = t(A)%*%A
        if (verbose){
          message("* lsolve.ssor : making it normal equation form via 'adjsym' flag.")
        }
      } else {
        stop("* lsolve.ssor : SSOR must be applied to symmetric matrix A.")
      }
    }
  }

  # Preprocessing : SSOR only : w
  if ((w<=0)||(w>=2)){
    stop("* lsolve.ssor : weight value w should be in (0,2).")
  }

  # Preprocessing : no NA or Inf
  if (any(is.infinite(A))||any(is.na(A))||any(is.infinite(B))||any(is.na(B))){
    stop("* lsolve.ssor : no NA, Inf, -Inf values are allowed.")
  }

  # Preprocessing : size argument : A and B
  m = nrow(A)
  if (is.vector(B)){
    mB = length(B)
    if (m!=mB){
      stop("* lsolve.ssor : a vector B should have a length of nrow(A).")
    }
  } else {
    mB = nrow(B)
    if (m!=mB){
      stop("* lsolve.ssor : an input matrix B should have the same number of rows from A.")
    }
  }
  if (is.vector(B)){
    B = as.matrix(B)
  }
  # Preprocessing : size argument : A case
  # Overdetermined  - A'Ax = A'b
  # Underdetermined - not supporting this case.
  n = ncol(A)
  if (m<n){
    stop("* lsolve.ssor : underdetermined case is not supported.")
  } else if (m>n){
    B = (t(A)%*%B)
    A = (t(A)%*%A)
    if (verbose){
      message("* lsolve.ssor : overdetermined case : turning into normal equation.")
    }
  }

  # Preprocessing : aux.is.dd
  if (aux.is.dd(A)==FALSE){
    if (verbose){
      message("* lsolve.ssor : LHS matrix A is not diagonally dominant.")
    }
  }
  # Preprocessing : adjust diagonal entries for A
  if (any(diag(A)==0)){
    cvec     = rnorm(10)
    adjconst = cvec[sample(which(cvec!=0),1)]/(1e+5)
    diag(A)  = diag(A)+adjconst
  }

  # Main Computation
  ncolB = ncol(B)
  if (ncolB==1){
    if (!sparseflag){
      vecB = as.vector(B)
      res = linsolve.ssor.single(A,vecB,xinit,reltol,maxiter,w)
    } else {
      vecB = B
      res = linsolve.ssor.single.sparse(A,vecB,xinit,reltol,maxiter,w)
    }
  } else {
    x      = array(0,c(n,ncolB))
    iter   = array(0,c(1,ncolB))
    errors = list()
    for (i in 1:ncolB){
      if (!sparseflag){
        vecB = as.vector(B[,i])
        tmpres = linsolve.ssor.single(A,vecB,xinit,reltol,maxiter,w)
      } else {
        vecB = Matrix(B[,i],sparse=TRUE)
        tmpres = linsolve.ssor.single.sparse(A,vecB,xinit,reltol,maxiter,w)
      }
      x[,i]       = tmpres$x
      iter[i]     = tmpres$iter
      errors[[i]] = tmpres$errors
      if (verbose){
        message(paste("* lsolve.ssor : B's column.",i,"being processed.."))
      }
    }
    res = list("x"=x,"iter"=iter,"errors"=errors)
  }

  # Return
  if (verbose){
    message("* lsolve.ssor : computations finished.")
  }
  return(res)
}

