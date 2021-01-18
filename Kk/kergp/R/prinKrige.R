## DO NOT ROXYGENISE THAT!!!
##' Principal Kriging Functions.
##'
##' The Principal Kriging Functions (PKF) are the eigenvectors of a
##' symmetric positive matrix \eqn{\mathbf{B}}{B} named the
##' \emph{Bending Energy Matrix} with dimension \eqn{n \times n}{n *
##' n} and rank {n - p}. The PKF are given in the \emph{ascending}
##' order of the eigenvalues \eqn{e_i}{e[i]} \deqn{e_1 = e_2 = \dots =
##' e_p = 0 < e_{p+1} \leq e_{p+2} \leq \dots \leq e_n.}{ e[1] = e[2]
##' = ... = e[p] = 0 < e[p + 1] <= e[p + 2] <= ... <= e[n].}  The
##' \eqn{p} first PKF generate the same space as do the \eqn{p}
##' columns of the trend matrix \eqn{\mathbf{F}}{F}, say
##' \eqn{\text{colspan}(F)}{colspan(F)}. The following \eqn{n-p} PKFs
##' generate a supplementary of the subspace
##' \eqn{\text{colspan}(F)}{colspan(F)}, and they have a decreasing
##' influence on the response. So the \eqn{p +1}-th PKF can give a
##' hint on a possible deterministic trend functions that could be
##' added to the \eqn{p} existing ones.
##'
##' The matrix \eqn{\mathbf{B}}{B} is such that \eqn{\mathbf{B}
##' \mathbf{F} = \mmathbf{0}}{B F = 0}, so the colums of
##' \eqn{\mathb{F}}{F} can be thought of as the eigenvectors that are
##' associated with the null eigenvalues \eqn{e_1}{e[1]},
##' \eqn{\dots}{...}, \eqn{e_p}{e[p]}.
##' 
##' @title Principal Kriging Functions
##' 
##' @param object An object with class \code{"gp"}.
##'
##' @return A list
##' \itemize{
##'
##' \item{\code{values}}{The eigenvalues of the energy bending matrix
##' in \emph{ascending} order. The first \eqn{p} values must be very
##' close to zero, but will not be zero since they are provided by
##' numeric linear algebra.}
##'
##' \item{\code{vectors}}{A matrix \eqn{\mathbf{U}}{U} with its
##' columns \eqn{\mathbd{u}_i}{U[ , i]} equal to the eigenvectors of
##' the energy bending matrix, in correspondence with the
##' eigenvalues \eqn{e_i}{e[i]}.}
##'
##' \item{\code{B}}{The Energy Bending Matrix
##' \eqn{\mathbf{B}}{B}. Remind that the eigenvectors are used here in
##' the ascending order of the eigenvalues, which is the reverse of
##' the usual order.}
##'
##' }
##'
##' @note When an eigenvalue \eqn{e_i}{e[i]} is such that \eqn{e_{i-1}
##' < e_i < e_{i+1}}{e[i-1] < e[i] < e[i+1]} (which can happen only
##' for \eqn{i > p}), the corresponding PKF is unique up to a change
##' of sign. However a run of \eqn{r > 1} identical eigenvalues is associated
##' with a \eqn{r}-dimensional eigenspace and the corresponding PKFs
##' have no meaning when they are considered individually.
##' 
##' @references
##'
##' Sahu S.K. and Mardia K.V. (2003).  A Bayesian kriged Kalman model
##' for short-term forecasting of air pollution
##' levels. \emph{Appl. Statist.} 54 (1), pp. 223-244.
##' 
##' @author Yves Deville
##'
##' @examples
##' library(kergp)
##' set.seed(314159)
##' n <- 100
##' x <- sort(runif(n))
##' y <- 2 + 4 * x  + 2 * x^2 + 3 * sin(6 * pi * x ) + 1.0 * rnorm(n)
##' nNew <- 60; xNew <- sort(runif(nNew))
##' df <- data.frame(x = x, y = y)
##'
##' ##-------------------------------------------------------------------------
##' ## use a Matern 3/2 covariance and a mispecified trend. We should guess
##' ## that it lacks a mainily linear and slightly quadratic part.
##' ##-------------------------------------------------------------------------
##'
##' myKern <- k1Matern3_2
##' inputNames(myKern) <- "x"
##' mygp <- gp(formula = y ~ sin(6 * pi * x),
##'            data = df, inputs = "x",
##'            parCovLower = c(0.01, 0.01), parCovUpper = c(10, 100),
##'            cov = myKern, estim = TRUE, noise = TRUE)
##' PK <- prinKrige(mygp)
##'
##' ## the third PKF suggests a possible linear trend term, and the
##' ## fourth may suggest a possible quadratic linear trend
##' 
##' matplot(x, PK$vectors[ , 1:4], type = "l", lwd = 2)
##' 
prinKrige <- function(object) {
   
    n <- nrow(object$F)
    
    B <- forwardsolve(t(object$RStar), t(object$FStar))
    B <- -crossprod(B)
    diag(B) <- diag(B) + 1.0
    B <- t(backsolve(t(object$L), B))
    B <- backsolve(t(object$L), B)

    eB <- eigen(B)
    e <- rev(eB$values)
    U <- eB$vectors[ , rev(1:n)]
    list(values = e, vectors = U, EBM = B)
    
}
