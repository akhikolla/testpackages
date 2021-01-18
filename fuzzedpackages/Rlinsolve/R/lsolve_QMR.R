#' Quasi Minimal Residual Method
#'
#' Quasia-Minimal Resudial(QMR) method is another remedy of the BiCG which shows
#' rather irregular convergence behavior. It adapts to solve the reduced tridiagonal system
#' in a least squares sense and its convergence is known to be quite smoother than BiCG.
#'
#' @param A an \eqn{(m\times n)} dense or sparse matrix. See also \code{\link[Matrix]{sparseMatrix}}.
#' @param B a vector of length \eqn{m} or an \eqn{(m\times k)} matrix (dense or sparse) for solving \eqn{k} systems simultaneously.
#' @param xinit a length-\eqn{n} vector for initial starting point. \code{NA} to start from a random initial point near 0.
#' @param reltol tolerance level for stopping iterations.
#' @param maxiter maximum number of iterations allowed.
#' @param preconditioner an \eqn{(n\times n)} preconditioning matrix; default is an identity matrix.
#' @param verbose a logical; \code{TRUE} to show progress of computation.
#'
#' @return a named list containing \describe{
#' \item{x}{solution; a vector of length \eqn{n} or a matrix of size \eqn{(n\times k)}.}
#' \item{iter}{the number of iterations required.}
#' \item{errors}{a vector of errors for stopping criterion.}
#' }
#'
#' @examples
#' \dontrun{
#' ## Overdetermined System
#' set.seed(100)
#' A = matrix(rnorm(10*5),nrow=10)
#' x = rnorm(5)
#' b = A%*%x
#'
#' out1 = lsolve.cg(A,b)
#' out2 = lsolve.bicg(A,b)
#' out3 = lsolve.qmr(A,b)
#' matout = cbind(matrix(x),out1$x, out2$x, out3$x);
#' colnames(matout) = c("true x","CG result", "BiCG result", "QMR result")
#' print(matout)
#' }
#'
#' @references
#' \insertRef{freund_qmr:_1991}{Rlinsolve}
#'
#' @rdname krylov_QMR
#' @export
lsolve.qmr <- function(A,B,xinit=NA,reltol=1e-5,maxiter=1000,
                       preconditioner=diag(ncol(A)),verbose=TRUE){
  ###########################################################################
  # Step 0. Initialization
  if (verbose){
    message("* lsolve.qmr : Initialiszed.")
  }
  if (any(is.na(A))||any(is.infinite(A))||any(is.na(B))||any(is.infinite(B))){
    stop("* lsolve.qmr : no NA or Inf values allowed.")
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
      stop("* lsolve.qmr : please use a valid 'xinit'.")
    }
  } else {
    if (length(xinit)!=ncol(A)){
      stop("* lsolve.qmr : 'xinit' has invalid size.")
    }
    xinit = matrix(xinit)
  }
  ###########################################################################
  # Step 1. Preprocessing
  # 1-1. Neither NA nor Inf allowed.
  if (any(is.infinite(A))||any(is.na(A))||any(is.infinite(B))||any(is.na(B))){
    stop("* lsolve.qmr : no NA, Inf, -Inf values are allowed.")
  }
  # 1-2. Size Argument
  m = nrow(A)
  if (is.vector(B)){
    mB = length(B)
    if (m!=mB){
      stop("* lsolve.qmr : a vector B should have a length of nrow(A).")
    }
  } else {
    mB = nrow(B)
    if (m!=mB){
      stop("* lsolve.qmr : an input matrix B should have the same number of rows from A.")
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
    stop("* lsolve.qmr : underdetermined case is not supported.")
  }
  # 1-4. Preconditioner : only valid for square case
  if (!all.equal(dim(A),dim(preconditioner))){
    stop("* lsolve.qmr : Preconditioner is a size-matching.")
  }
  if (verbose){message("* lsolve.qmr : preprocessing finished ...")}
  ###########################################################################
  # Step 2. Main Computation
  ncolB = ncol(B)
  if (ncolB==1){
    if (sparseflag){
      luM = lu(preconditioner)
      M1 = luM@L
      M2 = luM@U
    }
    if (!sparseflag){
      vecB = as.vector(B)
      res = linsolve.qmr.single(A,vecB,xinit,reltol,maxiter,preconditioner)
    } else {
      vecB = B
      res = linsolve.qmr.single.sparse(A,vecB,xinit,reltol,maxiter,preconditioner,M1,M2)
    }
  } else {
    x      = array(0,c(ncol(A),ncolB))
    iter   = array(0,c(1,ncolB))
    errors = list()
    if (sparseflag){
      luM = lu(preconditioner)
      M1 = luM@L
      M2 = luM@U
    }
    for (i in 1:ncolB){
      if (!sparseflag){
        vecB = as.vector(B[,i])
        tmpres = linsolve.qmr.single(A,vecB,xinit,reltol,maxiter,preconditioner)
      } else {
        vecB = Matrix(B[,i],sparse=TRUE)
        tmpres = linsolve.qmr.single.sparse(A,vecB,xinit,reltol,maxiter,preconditioner,M1,M2)
      }
      x[,i]        = tmpres$x
      iter[i]      = tmpres$iter
      errors[[i]]  = tmpres$errors
      if (verbose){
        message(paste("* lsolve.qmr : B's column.",i,"being processed.."))
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
        message("* lsolve.qmr : convergence well achieved.")
      }
    } else if (flagval==1){
      if (verbose){
        message("* lsolve.qmr : convergence not achieved within maxiter.")
      }
    } else if (flagval==-1){
      if (verbose){
        message("* lsolve.qmr : breakdown due to degenerate 'rho' value.")
      }
    } else if (flagval==-2){
      if (verbose){
        message("* lsolve.qmr : breakdown due to degenerate 'beta' value.")
      }
    } else if (flagval==-3){
      if (verbose){
        message("* lsolve.qmr : breakdown due to degenerate 'gamma' value.")
      }
    } else if (flagval==-4){
      if (verbose){
        message("* lsolve.qmr : breakdown due to degenerate 'delta' value.")
      }
    } else if (flagval==-5){
      if (verbose){
        message("* lsolve.qmr : breakdown due to degenerate 'ep' value.")
      }
    } else if (flagval==-6){
      if (verbose){
        message("* lsolve.qmr : breakdown due to degenerate 'xi' value.")
      }
    }
    res$flag = NULL
  }
  if (verbose){
    message("* lsolve.qmr : computations finished.")
  }
  return(res)
}
