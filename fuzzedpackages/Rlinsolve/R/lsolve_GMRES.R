#' Generalized Minimal Residual method
#'
#' GMRES is a generic iterative solver for a nonsymmetric system of linear equations. As its name suggests, it approximates
#' the solution using Krylov vectors with minimal residuals.
#'
#' @param A an \eqn{(m\times n)} dense or sparse matrix. See also \code{\link[Matrix]{sparseMatrix}}.
#' @param B a vector of length \eqn{m} or an \eqn{(m\times k)} matrix (dense or sparse) for solving \eqn{k} systems simultaneously.
#' @param xinit a length-\eqn{n} vector for initial starting point. \code{NA} to start from a random initial point near 0.
#' @param reltol tolerance level for stopping iterations.
#' @param maxiter maximum number of iterations allowed.
#' @param preconditioner an \eqn{(n\times n)} preconditioning matrix; default is an identity matrix.
#' @param restart the number of iterations before restart.
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
#' out1 = lsolve.cg(A,b)
#' out3_1 = lsolve.gmres(A,b,restart=2)
#' out3_2 = lsolve.gmres(A,b,restart=3)
#' out3_3 = lsolve.gmres(A,b,restart=4)
#' matout = cbind(matrix(x),out1$x, out3_1$x, out3_2$x, out3_3$x);
#' colnames(matout) = c("true x","CG", "GMRES(2)", "GMRES(3)", "GMRES(4)")
#' print(matout)
#'
#'
#' @references
#' \insertRef{saad_gmres:_1986}{Rlinsolve}
#'
#'
#' @rdname krylov_GMRES
#' @export
lsolve.gmres <- function(A,B,xinit=NA,reltol=1e-5,maxiter=1000,
                         preconditioner=diag(ncol(A)),restart=(ncol(A)-1),verbose=TRUE){
  ###########################################################################
  # Step 0. Initialization
  if (verbose){
    message("* lsolve.gmres : Initialiszed.")
  }
  if (any(is.na(A))||any(is.infinite(A))||any(is.na(B))||any(is.infinite(B))){
    stop("* lsolve.gmres : no NA or Inf values allowed.")
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
      stop("* lsolve.gmres : please use a valid 'xinit'.")
    }
  } else {
    if (length(xinit)!=ncol(A)){
      stop("* lsolve.gmres : 'xinit' has invalid size.")
    }
    xinit = matrix(xinit)
  }
  
  if ((restart<2)||(is.na(restart))||(is.infinite(restart))||(abs(restart-round(restart))>sqrt(.Machine$double.eps))){
    stop("* lsolve.gmres : 'restart' should be a positive integer >= 2.")
  }
  restart = round(restart)
  if (restart>=ncol(A)){
    stop("* lsolve.gmres : take a restart value smaller than ncol(A).")
  }
  ###########################################################################
  # Step 1. Preprocessing
  # 1-1. Neither NA nor Inf allowed.
  if (any(is.infinite(A))||any(is.na(A))||any(is.infinite(B))||any(is.na(B))){
    stop("* lsolve.gmres : no NA, Inf, -Inf values are allowed.")
  }
  # 1-2. Size Argument
  m = nrow(A)
  if (is.vector(B)){
    mB = length(B)
    if (m!=mB){
      stop("* lsolve.gmres : a vector B should have a length of nrow(A).")
    }
  } else {
    mB = nrow(B)
    if (m!=mB){
      stop("* lsolve.gmres : an input matrix B should have the same number of rows from A.")
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
    stop("* lsolve.gmres : underdetermined case is not supported.")
  }
  # 1-4. Preconditioner : only valid for square case
  if (!all.equal(dim(A),dim(preconditioner))){
    stop("* lsolve.gmres : Preconditioner is a size-matching.")
  }
  if (verbose){message("* lsolve.gmres : preprocessing finished ...")}
  ###########################################################################
  # Step 2. Main Computation
  ncolB = ncol(B)
  if (ncolB==1){
    if (!sparseflag){
      vecB = as.vector(B)
      res = linsolve.gmres.single(A,vecB,xinit,reltol,maxiter,preconditioner,restart)
    } else {
      vecB = B
      res = linsolve.gmres.single.sparse(A,vecB,xinit,reltol,maxiter,preconditioner,restart)
    }
  } else {
    x      = array(0,c(ncol(A),ncolB))
    iter   = array(0,c(1,ncolB))
    errors = list()
    for (i in 1:ncolB){
      if (!sparseflag){
        vecB = as.vector(B[,i])
        tmpres = linsolve.gmres.single(A,vecB,xinit,reltol,maxiter,preconditioner,restart)
      } else {
        vecB = Matrix(B[,i],sparse=TRUE)
        tmpres = linsolve.gmres.single.sparse(A,vecB,xinit,reltol,maxiter,preconditioner,restart)
      }
      x[,i]        = tmpres$x
      iter[i]      = tmpres$iter
      errors[[i]]  = tmpres$errors
      if (verbose){
        message(paste("* lsolve.gmres : B's column.",i,"being processed.."))
      }
    }
    res = list("x"=x,"iter"=iter,"errors"=errors)
  }

  ###########################################################################
  # Step 3. Finalize
  if ("flag" %in% names(res)){
    flagval = res$flag
    if (flagval==0){
      if (verbose){
        message("* lsolve.gmres : convergence well achieved.")
      }
    } else if (flagval==1){
      if (verbose){
        message("* lsolve.gmres : convergence not achieved within maxiter.")
      }
    } else {
      if (verbose){
        message("* lsolve.gmres : breakdown.")
      }
    }
    res$flag = NULL
  }
  if (verbose){
    message("* lsolve.gmres : computations finished.")
  }
  return(res)
}
