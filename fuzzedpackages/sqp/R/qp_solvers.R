#' @title
#' Quadratic optimization solver
#' 
#' @description
#' Dense & Sparse solvers for linearly constrained quadratic optimization problems
#' \insertCite{@cf. @fletcher71; @nocedal99; @powell78; @wilson63}{sqp}.
#' 
#' @details 
#' Sequential quadratic programming relies on
#' iteratively solving linear approximations of the optimality conditions 
#' \insertCite{@cf. @kjeldsen00; @kuhn51}{sqp}.

#'  
#' This is equivalent to minimizing a quadratic approximation of the distance function 
#' under linearised constraint functions. 
#' \code{qp_solver} can be used to solve this quadratic sub-problem.
#' Solving a quadratic problem under linear equalits constraints is equivalent to 
#' solving a system of linear equations.  
#' The inequality constraints are handeled by an active set strategy, where the
#' binding ones are treated as equalities, and the active set is found iteratively
#' \insertCite{@cf. @fletcher71; @nocedal99; @powell78; @wilson63}{sqp}.
#' 
#' @note 
#' Although there is already an implementation for using the
#' SuperLU sparse solver within this package,
#' it is currently disabled due to licensing considerations.
#' 
#' Sparse matrices are converted to dense ones in the solving procedure.
#' 
#' Hopefully, this can be updated in the near future.
#' 
#' @param Q,C_eq,C_ineq
#' \strong{Dense or sparse numeric matrices}:
#' \describe{
#' \item{\bold{Q}      }{\eqn{N \times N}{N x N}-\bold{matrix}:\cr
#'                       Quadratic distance (loss) multiplier for the optimization problem.}
#' \item{\bold{C_eq}   }{\eqn{N_{eq} \times N}{N_eq x N}-\bold{matrix}:\cr
#'                       Equality constraint multiplier for the \eqn{N_{eq}}{N_eq} equality constraints.}
#' \item{\bold{C_ineq} }{\eqn{N_{ineq} \times N}{N_ineq x N}-\bold{matrix}:\cr
#'                       Inequality constraint multiplier for the \eqn{N_{ineq}}{N_ineq} inequality constraints.}
#' }
#'  
#' @param l,t_eq,t_ineq
#' \strong{Numeric vectors}:
#' \describe{
#' \item{\bold{l}      }{\bold{Vector} of size \eqn{N}{N}:\cr
#'                       Linear distance (loss) multiplier for the optimization problem.}
#' \item{\bold{t_eq}   }{\bold{Vector} of size \eqn{N_{eq}}{N_eq}:\cr
#'                       Targets for equality constraints.}
#' \item{\bold{t_ineq} }{\bold{Vector} of size \eqn{N_{ineq}}{N_ineq}:
#'                       \bold{upper} bound for inequality constraints.}
#' }
#' 
#' @param x
#' \strong{Numeric vector}                   
#' of size \code{N}:\cr
#' Initial values for optimization parameters.
#' Slack variables are only used for constraints violated by this 
#' \code{x}, unless \code{all_slack} is \code{TRUE}.
#' 
#' @param penalty
#' \strong{Numeric value}:\cr
#' Penalty multiplier for slack variables in distance function.
#' 
#' @param tol
#' \strong{Numeric value}:\cr
#' Tolerance for assessing convgergence criteria & constraints.
#' 
#' @param max_iter
#' \strong{Integer value}:\cr
#' Tolerance for assessing convgergence criteria & constraints.
#' 
#' @param fast
#' \strong{Boolean}:\cr
#' Whether to use faster (but lower quality)
#' solver (cf. \href{http://arma.sourceforge.net/docs.html#solve}{Armadillo documentation}:\cr
#' fast mode: disable determining solution quality via rcond, disable iterative refinement, disable equilibration.
#' 
#' @param all_slack 
#' \strong{Boolean}:\cr
#' Whether to use slack variables for all constraints instead of \cr
#' only for the ones violated by the initial values
#' 
#' @param debug 
#' \strong{Boolean}:\cr
#' Whether to print debugging status messages.
#' 
#' @param solver
# \strong{Integer value}:\cr
#' Solver identification used for optimization in the \bold{dense} matrix case. Not yet used.
#' 
#' @return 
#' A \strong{named list} with values
#' \describe{
#' \item{\bold{x}}{Final values for optimization parameters}
#' \item{\bold{lagrange_eq}, \bold{lagrange_ineq}}{Lagrange multipliers for equality and inequality constraints}
#' \item{\bold{slack_eq_positive}, \bold{slack_eq_negative}}{Positive and negative slack variables for equality constraints}
#' \item{\bold{slack_ineq}}{Slack variables for inequalits constraints}
#' \item{\bold{lagrange_slack_eq_positive}, \bold{lagrange_slack_eq_negative}, \bold{lagrange_slack_ineq}}{Lagrange multipliers for positivity of slack variables}
#' }
#' 
#' 
#' @examples
#' 
#' set.seed(1)
#' n <- 5
#' 
#' x_init <- cbind(runif(n))
#' 
#' w <- runif(n)
#' 
#' Q <- 3*diag(n) # minimize sum(3*x^2 + 3*x)
#' l <- cbind(rep(3,n)) # minimize sum(3*x^2 + 3*x)
#' 
#' C_eq <- rbind(1,w) # constraints: sum(x) == 1, sum(w*x) == 5
#' C_ineq <- rbind(diag(n),-diag(n)) # constraints: all(x >= -4) & all(x <= 4)
#' 
#' t_eq <- rbind(1,5) # constraints: sum(x) == 1, sum(w*x) == 5 
#' t_ineq <- cbind(rep(c(4,4),each=n)) # constraints: all(x >= -4) & all(x <= 4) 
#' 
#' output <- qp_solver(Q = Q,
#'                C_eq = C_eq,
#'                C_ineq = C_ineq,
#'                l=l,
#'                t_eq = t_eq,
#'                t_ineq = t_ineq,
#'                x = x_init,
#'                tol = 1e-15)
#' 
#' 
#' sum(output$x) # constraints: sum(x) == 1
#' sum(w*output$x) # constraints: sum(w*x) == 5
#' 
#' all(output$x >= -4) # constraints: all(x >= -4)
#' all(output$x <= 4) # constraints: all(x <= 4)

#' 
#' @references \insertAllCited{}
#' 
#' @export

qp_solver <- function(                      
  Q,        # Quadratic Multiplier
  C_eq      = NULL,  # Equality  Constraint Multiplier
  C_ineq    = NULL,  # Inquality Constraint Multiplier
  l         = NULL,  # Linear Multiplier
  t_eq      = NULL,  # Equality   Constraint RHS
  t_ineq    = NULL,  # Inequality Constraint RHS (upper bound)
  x         = NULL,  # Optimization variable (initial value)
  penalty   = 1e+10, # Penalty parameter for slack variables
  tol       = 1e-7,
  max_iter  = 500,
  fast      = FALSE,
  all_slack = FALSE,
  debug     = FALSE,
  solver    = 0)
{
  if(is.null(x))
  {
    x <- rep(0,nrow(Q))
  }
  
  if(is.null(nrow(x)))
  {
    x <- cbind(x)
  }
  
  if(!is.numeric(x))
    stop("Argument 'x' must be numeric!\n")
  
  if(nrow(x)!=ncol(Q))
    stop(paste0("Number of elements in x (",nrow(x),") != number of rows / columns in Q (",ncol(Q),")"))
  
  if(!is.null(l))
  {
    if(is.null(nrow(l)))
    {
      l <- cbind(l)
    }
    if(nrow(l)!=ncol(Q))
      stop(paste0("Number of elements in l (",nrow(l),") != number of rows / columns in Q (",ncol(Q),")"))
  }
  
  if(!is.null(C_eq))
  {
    if(!is.null(t_eq))
    {
      if(is.null(nrow(t_eq)))
      {
        t_eq <- cbind(t_eq)
      }
      
      if(nrow(t_eq)!=nrow(C_eq))
        stop(paste0("Number of elements in t_eq (",nrow(t_eq),") != number of rows in C_eq (",nrow(C_eq),")"))
      
    } else
    {
      stop("Argument 'C_eq' specified, but not 't_eq'!\n")
    }
  } else
  {
    if(!is.null(nrow(t_eq)))
    {
      stop("Argument 't_eq' specified, but not 'C_eq'!\n")
    }
  }
  
  if(!is.null(C_ineq))
  {
    if(!is.null(t_ineq))
    {
      if(is.null(nrow(t_ineq)))
      {
        t_ineq <- cbind(t_ineq)
      }
      
      if(nrow(t_ineq)!=nrow(C_ineq))
        stop(paste0("Number of elements in t_ineq (",nrow(t_ineq),") != number of rows in C_ineq (",nrow(C_ineq),")"))
      
    } else
    {
      stop("Argument 'C_ineq' specified, but not 't_ineq'!\n")
    }
  } else
  {
    if(!is.null(nrow(t_ineq)))
    {
      stop("Argument 't_ineq' specified, but not 'C_ineq'!\n")
    }
  }
  
  sparse <- c(Q=(class(Q) == "dsCMatrix"),
              C_eq=(class(C_eq) == "dsCMatrix"),
              C_ineq=(class(C_ineq) == "dsCMatrix"))
  
  dim_eq    <- nrow(C_eq)
  dim_ineq  <- nrow(C_ineq)
  dim_Q     <- nrow(Q)
  
  
  if(sum(sparse)>0)
  {
    if(sum(sparse)<3)
    {
      warn_text <- paste0(paste0(names(sparse)[which(sparse)], collapse = " & "),
                          if(sum(sparse)>1) " are sparse matrices, but " else " is sparse matrix, but ",
                          paste0(names(sparse)[which(!sparse)], collapse = " & "),
                          if(sum(sparse)<2) " are " else " is ",
                          "not.")

      if(!sparse["Q"])
        Q <- Matrix::Matrix(Q, sparse = TRUE)
      if(!sparse["C_eq"])
        C_eq <- Matrix::Matrix(C_eq, sparse = TRUE)
      if(!sparse["C_ineq"])
        C_ineq <- Matrix::Matrix(C_ineq, sparse = TRUE)
      
      sparse[] <- TRUE
    }
    
    return(.solvers_slacked_sparse(x,
                                   Q,
                                   C_eq,
                                   C_ineq,
                                   l,
                                   t_eq,
                                   t_ineq,
                                   penalty,
                                   tol,
                                   max_iter,
                                   dim_eq,
                                   dim_ineq,
                                   dim_Q,
                                   all_slack,
                                   debug))
  } else
  {
    return(.solvers_slacked_dense(x,
                                  Q,
                                  C_eq,
                                  C_ineq,
                                  l,
                                  t_eq,
                                  t_ineq,
                                  penalty,
                                  tol,
                                  max_iter,
                                  dim_eq,
                                  dim_ineq,
                                  dim_Q,
                                  0,
                                  fast,
                                  all_slack,
                                  debug))
  }
}