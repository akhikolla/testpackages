##' @importFrom Rcpp evalCpp
##' @importFrom methods is
##' @useDynLib lbfgsb3c, .registration=TRUE
"lbfgsb3c"


##' Interfacing wrapper for the Nocedal - Morales LBFGSB3 (Fortran) limited memory BFGS solver.
##'
##' @param par A parameter vector which gives the initial guesses to
##'     the parameters that will minimize \code{fn}. This can be
##'     named, for example, we could use par=c(b1=1, b2=2.345,
##'     b3=0.123)
##' @param fn A function that evaluates the objective function to be
##'     minimized.  This can be a R function or a Rcpp function
##'     pointer.
##' @param gr If present, a function that evaluates the gradient
##'     vector for the objective function at the given parameters
##'     computing the elements of the sum of squares function at the
##'     set of parameters \code{start}. This can be a R function or a
##'     Rcpp function pointer.
##' @param lower Lower bounds on the parameters. If a single number,
##'     this will be applied to all parameters. Default -Inf.
##' @param upper Upper bounds on the parameters. If a single number,
##'     this will be applied to all parameters. Default Inf.
##' @param control An optional list of control settings. See below in
##'     details.
##' @param ... Any data needed for computation of the objective
##'     function and gradient.
##' @param rho An Environment to use for function evaluation.  If
##'     present the arguments in ... are ignored.  Otherwise the
##'     ... are converted to an environment for evaluation.
##' @details
##'
##' See the notes below for a general appreciation of this package.
##'
##' The control list can contain:
##' \itemize{
##' \item{trace} If positive, tracing information on the progress of the optimization is produced. Higher values may produce more tracing information: for method "L-BFGS-B" there are six levels of tracing. (To understand exactly what these do see the source code: higher levels give more detail.)
##' \item{factr} controls the convergence of the "L-BFGS-B" method. Convergence occurs when the reduction in the objective is within this factor of the machine tolerance. Default is 1e7, that is a tolerance of about 1e-8.
##' \item{pgtol} helps control the convergence of the "L-BFGS-B" method. It is a tolerance on the projected gradient in the current search direction. This defaults to zero, when the check is suppressed.
##' \item{abstol} helps control the convergence of the "L-BFGS-B" method. It is an absolute tolerance difference in x values. This defaults to zero, when the check is suppressed.
##' \item{reltol} helps control the convergence of the "L-BFGS-B" method. It is an relative tolerance difference in x values. This defaults to zero, when the check is suppressed.
##' \item{lmm} is an integer giving the number of BFGS updates retained in the "L-BFGS-B" method, It defaults to 5.
##' \item{maxit} maximum number of iterations.
##' \item{iprint} Provided only for compatibility with older codes. This control is no longer active.)
##' \item{info} a boolean to indicate if more optimization information is captured and output in a $info list
##' }
##'
##' @return
##'   A list of the following items
##' \itemize{
##' \item{par} The best set of parameters found.
##' \item{value} The value of fn corresponding to par.
##' \item{counts} A two-element integer vector giving the number of calls to fn and gr respectively. This excludes any calls to fn to compute a finite-difference approximation to the gradient.
##' \item{convergence} An integer code. 0 indicates successful completion
##' }
##' @seealso Packages \code{\link{optim}} and \code{optimx}.
##' @keywords nonlinear parameter optimization
##' @author Matthew Fidler (move to C and add more options for adjustments),
##'     John C Nash <nashjc@uottawa.ca> (of the wrapper and edits to Fortran code to allow R output)
##'     Ciyou Zhu, Richard Byrd, Jorge Nocedal, Jose Luis Morales (original Fortran packages)
##'
##' @references
##'     Morales, J. L.; Nocedal, J. (2011). "Remark on 'algorithm 778: L-BFGS-B:
##'           Fortran subroutines for large-scale bound constrained optimization' ".
##'           ACM Transactions on Mathematical Software 38: 1.
##'
##'     Byrd, R. H.; Lu, P.; Nocedal, J.; Zhu, C. (1995). "A Limited Memory Algorithm
##'           for Bound Constrained Optimization". SIAM J. Sci. Comput. 16 (5): 1190-1208.
##'
##'     Zhu, C.; Byrd, Richard H.; Lu, Peihuang; Nocedal, Jorge (1997). "L-BFGS-B:
##'           Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained
##'           optimization". ACM Transactions on Mathematical Software 23 (4): 550-560.
##'
##' @note
##'   This package is a wrapper to the Fortran code released by Nocedal and Morales.
##'   This poses several difficulties for an R package. While the \code{.Fortran()}
##'   tool exists for the interfacing, we must be very careful to align the arguments
##'   with those of the Fortran subroutine, especially in type and storage.
##'
##'   A more annoying task for interfacing the Fortran code is that Fortran WRITE or
##'   PRINT statements must all be replaced with calls to special R-friendly output
##'   routines. Unfortunately, the Fortran is full of output statements. Worse, we may
##'   wish to be able to suppress such output, and there are thus many modifications
##'   to be made. This means that an update of the original code cannot be simply
##'   plugged into the R package \code{src} directory.
##'
##'   Finally, and likely because L-BFGS-B has a long history, the Fortran code is far
##'   from well-structured. For example, the number of function and gradient evaluations
##'   used is returned as the 34'th element of an integer vector. There does not appear
##'   to be an easy way to stop the program after some maximum number of such evaluations
##'   have been performed.
##'
##'   On the other hand, the version of L-BFGS-B in \code{optim()} is a \code{C} translation
##'   of a now-lost Fortran code. It does not implement the improvements Nocedal and
##'   Morales published in 2011. Hence, despite its deficiencies, this wrapper has been
##'   prepared.
##'
##' In addition to the above reasons for the original lbfgsb3 package,
##' this additional package allows C calling of L-BFGS-B 3.0 by a
##' program as well as adjustments to the tolerances that were not
##' present in the original CRAN package.  Also adjustments were made
##' to have outputs conform with R's optim routine.
##' @examples
##' # Rosenbrock's banana function
##' n=3; p=100
##'
##' fr = function(x)
##' {
##'     f=1.0
##'     for(i in 2:n) {
##'         f=f+p*(x[i]-x[i-1]**2)**2+(1.0-x[i])**2
##'     }
##'     f
##' }
##'
##' grr = function(x)
##' {
##'     g = double(n)
##'     g[1]=-4.0*p*(x[2]-x[1]**2)*x[1]
##'     if(n>2) {
##'         for(i in 2:(n-1)) {
##'             g[i]=2.0*p*(x[i]-x[i-1]**2)-4.0*p*(x[i+1]-x[i]**2)*x[i]-2.0*(1.0-x[i])
##'         }
##'     }
##'     g[n]=2.0*p*(x[n]-x[n-1]**2)-2.0*(1.0-x[n])
##'     g
##' }
##' x = c(a=1.02, b=1.02, c=1.02)
##' (opc <- lbfgsb3c(x,fr, grr))
##' (op <- lbfgsb3(x,fr, grr, control=list(trace=1)))
##' (opx <- lbfgsb3x(x,fr, grr))
##' (opf <- lbfgsb3f(x,fr, grr))
##' @export
lbfgsb3c <- function(par, fn, gr=NULL, lower = -Inf, upper = Inf,
                     control=list(), ..., rho=NULL){
# control defaults -- idea from spg
    ctrl <- list(trace= 0L,
                 maxit=100L,
                 iprint= -1L,
                 lmm=5,
                 factr=1e7,
                 pgtol=0,
                 reltol=0,
                 abstol=0,
                 info=FALSE);
    callstak <- sys.calls() # get the call stack
    lcs <- length(callstak)
    fstr <- as.character(callstak[lcs])
    fstr <- strsplit(fstr, "(", fixed=TRUE)[[1]][1]
    if (ctrl$trace > 0) { cat("Using function ",fstr,"\n") }
    if ( (fstr == "lbfgsb3") || (fstr == "lbfgsb3f") ) { ctrl$info <- TRUE }
    # This emits more information from lbfgsb3 Fortran code.
    namc <- names(control)
    if (!all(namc %in% names(ctrl)))
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    if (missing(rho)){
        rho <- as.environment(list(...));
    }
    if (is.null(gr)){
        gr <- function(x, ...){
            numDeriv::grad(fn, x, ...);
        }
    }
    if (is(fn, "function") & is (gr, "function")){
##        cat("USING fnR, grR\n")
        fnR <- function(x, rho){
            do.call(fn, c(list(x), as.list(rho)));
        }
        grR <- function(x, rho){
            do.call(gr, c(list(x), as.list(rho)));
        }
        return(lbfgsb3cpp(par, fnR, grR, lower, upper, ctrl, rho));
    } else {
        return(lbfgsb3cpp(par, fn, gr, lower, upper, ctrl, rho));
    }
} # end of lbfgsb3()

##'@rdname lbfgsb3c
##'@export
lbfgsb3 <- lbfgsb3c

##'@rdname lbfgsb3c
##'@export
lbfgsb3f <- lbfgsb3c

##'@rdname lbfgsb3c
##'@export
lbfgsb3x <- lbfgsb3c
