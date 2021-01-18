#' Jacobi method
#'
#' Jacobi method is an iterative algorithm for solving a system of linear equations,
#' with a decomposition \eqn{A = D+R} where \eqn{D} is a diagonal matrix.
#' For a square matrix \eqn{A}, it is required to be diagonally dominant. For an overdetermined system where \code{nrow(A)>ncol(A)},
#' it is automatically transformed to the normal equation. Underdetermined system -
#' \code{nrow(A)<ncol(A)} - is not supported.
#'
#' @param A an \eqn{(m\times n)} dense or sparse matrix. See also \code{\link[Matrix]{sparseMatrix}}.
#' @param B a vector of length \eqn{m} or an \eqn{(m\times k)} matrix (dense or sparse) for solving \eqn{k} systems simultaneously.
#' @param xinit a length-\eqn{n} vector for initial starting point. \code{NA} to start from a random initial point near 0.
#' @param reltol tolerance level for stopping iterations.
#' @param maxiter maximum number of iterations allowed.
#' @param weight a real number in \eqn{(0,1]}; 1 for native Jacobi.
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
#' out1 = lsolve.jacobi(A,b,weight=1,verbose=FALSE)   # unweighted
#' out2 = lsolve.jacobi(A,b,verbose=FALSE)            # weight of 0.66
#' out3 = lsolve.jacobi(A,b,weight=0.5,verbose=FALSE) # weight of 0.50
#' print("* lsolve.jacobi : overdetermined case example")
#' print(paste("*   error for unweighted    Jacobi case : ",norm(out1$x-x)))
#' print(paste("*   error for 0.66 weighted Jacobi case : ",norm(out2$x-x)))
#' print(paste("*   error for 0.50 weighted Jacobi case : ",norm(out3$x-x)))
#'
#' @references
#' \insertRef{demmel_applied_1997}{Rlinsolve}
#'
#' @rdname basic_JACOBI
#' @export
lsolve.jacobi <- function(A,B,xinit=NA,reltol=1e-5,
                          maxiter=1000,weight=2/3,adjsym=TRUE,verbose=TRUE){
  if (verbose){
    message("* lsolve.jacobi : Initialiszed.")
  }
  if (any(is.na(A))||any(is.infinite(A))||any(is.na(B))||any(is.infinite(B))){
    stop("* lsolve.jacobi : no NA or Inf values allowed.")
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
      stop("* lsolve.jacobi : please use a valid 'xinit'.")
    }
  } else {
    if (length(xinit)!=ncol(A)){
      stop("* lsolve.jacobi: 'xinit' has invalid size.")
    }
    xinit = matrix(xinit)
  }
  # Preprocessing : symmetricity warning
  if (nrow(A)==ncol(A)){
    if (norm(abs(t(A)-A),"f")>1e-10){
      if (verbose){
        message("* lsolve.jacobi : A may not be symmetric.")
      }
      if (adjsym){
        B = t(A)%*%B
        A = t(A)%*%A
        if (verbose){
          message("* lsolve.jacobi : making it normal equation form via 'adjsym' flag.")
        }
      }
    }
  }
  # Preprocessing : JACOBI ONLY : weight should be (0,1]
  if ((weight<=0)||(weight>1)){
    stop("* lsolve.jacobi : weight should be a positive real number in (0,1].")
  }
  # Preprocessing : no NA or Inf
  if (any(is.infinite(A))||any(is.na(A))||any(is.infinite(B))||any(is.na(B))){
    stop("* lsolve.jacobi : no NA, Inf, -Inf values are allowed.")
  }
  # Preprocessing : size argument : A and B
  m = nrow(A)
  if (is.vector(B)){
    mB = length(B)
    if (m!=mB){
      stop("* lsolve.jacobi : a vector B should have a length of nrow(A).")
    }
  } else {
    mB = nrow(B)
    if (m!=mB){
      stop("* lsolve.jacobi : an input matrix B should have the same number of rows from A.")
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
    stop("* lsolve.jacobi : underdetermined case is not supported.")
  } else if (m>n){
    B = (t(A)%*%B)
    A = (t(A)%*%A)
    if (verbose){
      message("* lsolve.jacobi : overdetermined case : turning into normal equation.")
    }
  }
  # Preprocessing : aux.is.dd
  if (aux.is.dd(A)==FALSE){
    if (verbose){
      message("* lsolve.jacobi : LHS matrix A is not diagonally dominant.")
      message("*               : solution from Jacobi method is not guaranteed.")
    }
  }

  # # Preprocessing : aux.is.dd
  # if (aux.is.dd(A)==FALSE){
  #   if (verbose){
  #     message("* lsolve.jacobi : LHS matrix A is not diagonally dominant.")
  #     message("*               : solution from Jacobi method is not guaranteed.")
  #   }
  # }
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
      res = linsolve.jacobi.single(A,vecB,xinit,reltol,maxiter,weight)
    } else {
      vecB = B
      res = linsolve.jacobi.single.sparse(A,vecB,xinit,reltol,maxiter,weight)
    }
  } else {
    x      = array(0,c(n,ncolB))
    iter   = array(0,c(1,ncolB))
    errors = list()
    for (i in 1:ncolB){
      if (!sparseflag){
        vecB = as.vector(B[,i])
        tmpres = linsolve.jacobi.single(A,vecB,xinit,reltol,maxiter,weight)
      } else {
        vecB = Matrix(B[,i],sparse=TRUE)
        tmpres = linsolve.jacobi.single.sparse(A,vecB,xinit,reltol,maxiter,weight)
      }

      x[,i]       = tmpres$x
      iter[i]     = tmpres$iter
      errors[[i]] = tmpres$errors
      if (verbose){
        message(paste("* lsolve.jacobi : B's column.",i,"being processed.."))
      }
    }
    res = list("x"=x,"iter"=iter,"errors"=errors)
  }

  # Return
  if (verbose){
    message("* lsolve.jacobi : computations finished.")
  }
  return(res)
}
