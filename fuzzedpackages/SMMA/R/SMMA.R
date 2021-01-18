#
#     Description of this R script:
#     R interface/wrapper for the Rcpp function pga in the SMMA package.
#
#     Intended for use with R.
#     Copyright (C) 2017 Adam Lund
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>
#

#' @name SMMA
#' @aliases pga
#' @title Soft Maximin Estimation for Large Scale Array Data with Known Groups
#'
#' @description  Efficient design matrix free procedure for solving a soft maximin problem for
#' large scale array-tensor structured models, see  \cite{Lund et al., 2020}.
#' Currently Lasso and SCAD penalized estimation is implemented.
#'
#' @usage  softmaximin(X,
#'             Y,
#'             zeta,
#'             penalty = c("lasso", "scad"),
#'             alg = c("npg", "fista"),
#'             nlambda = 30,
#'             lambda.min.ratio = 1e-04,
#'             lambda = NULL,
#'             penalty.factor = NULL,
#'             reltol = 1e-05,
#'             maxiter = 15000,
#'             steps = 1,
#'             btmax = 100,
#'             c = 0.0001,
#'             tau = 2,
#'             M = 4,
#'             nu = 1,
#'             Lmin = 0,
#'             log = TRUE)
#'
#' @param X list containing the Kronecker components (1, 2 or 3) of the Kronecker design matrix.
#'  These are  matrices of sizes \eqn{n_i \times p_i}.
#' @param Y  array of size \eqn{n_1 \times\cdots\times n_d \times G} containing the response values.
#' @param zeta strictly positive float controlling  the softmaximin approximation accuracy.
#' @param penalty string specifying the penalty type. Possible values are \code{"lasso", "scad"}.
#' @param alg string specifying the optimization algorithm. Possible values are \code{"npg", "fista"}.
#' @param nlambda positive integer giving the number of \code{lambda} values. Used when lambda is not specified.
#' @param lambda.min.ratio strictly positive float giving the smallest value for \code{lambda}, as a fraction of
#' \eqn{\lambda_{max}}; the (data dependent) smallest value for which all coefficients are zero.
#' Used when lambda is not specified.
#' @param lambda sequence of strictly positive floats used  as penalty parameters.
#' @param penalty.factor array of size \eqn{p_1 \times \cdots \times p_d} of positive floats. Is multiplied
#' with each element in \code{lambda} to allow differential penalization on the coefficients.
#' @param reltol strictly positive float giving the convergence tolerance for the inner loop.
#' @param maxiter positive integer giving the maximum number of  iterations allowed for each \code{lambda}
#' value, when  summing over all outer iterations for said \code{lambda}.
#' @param steps strictly positive integer giving the number of steps used in the multi-step adaptive lasso algorithm for non-convex penalties.
#' Automatically set to 1 when \code{penalty = "lasso"}.
#' @param btmax strictly positive integer giving the maximum number of backtracking steps allowed in each iteration. Default is \code{btmax = 100}.
#' @param c strictly positive float used in the NPG algorithm. Default is \code{c = 0.0001}.
#' @param tau strictly positive float used to control the stepsize for NPG. Default is \code{tau = 2}.
#' @param M positive integer giving the look back for the NPG. Default is \code{M = 4}.
#' @param nu strictly positive float used to control the stepsize. A  value less that 1 will decrease
#' the stepsize and a value larger than one will increase it. Default is \code{nu = 1}.
#' @param Lmin non-negative float used by the NPG algorithm to control the stepsize. For the default  \code{Lmin = 0}
#' the maximum step size is the same as for the FISTA algorithm.
#' @param log logical variable indicating whether to use log-loss.  TRUE is default and yields the loss below.

#' @details Following \cite{Lund et al., 2020}  this package solves the optimization problem for a linear
#' model for heterogeneous \eqn{d}-dimensional array data (\eqn{d=1,2,3}) organized in \eqn{G} known groups,
#' and with identical tensor structured design matrix \eqn{X} across all groups.  Specifically \eqn{n = \prod_i^d n_i} is the
#' number of observations in each group, \eqn{Y_g}  the  \eqn{n_1\times \cdots \times n_d} response array
#' for group \eqn{g \in \{1,\ldots,G\}}, and \eqn{X}  a \eqn{n\times p} design matrix, with tensor structure
#'  \deqn{X = \bigotimes_{i=1}^d X_i.}
#' For \eqn{d =1,2,3}, \eqn{X_1,\ldots, X_d} are the marginal \eqn{n_i\times p_i} design matrices (Kronecker components).
#' Using the GLAM framework  the model equation for group \eqn{g\in \{1,\ldots,G\}} is expressed as
#'  \deqn{Y_g = \rho(X_d,\rho(X_{d-1},\ldots,\rho(X_1,B_g))) + E_g,}
#' where \eqn{\rho} is the so called rotated \eqn{H}-transfrom (see  \cite{Currie et al., 2006}),
#' \eqn{B_g} for each \eqn{g} is a (random) \eqn{p_1\times\cdots\times p_d} parameter array
#'  and \eqn{E_g}  is  a \eqn{n_1\times \cdots \times n_d} error array.
#'
#' This package solves the penalized soft maximin problem from \cite{Lund et al., 2020}, given by
#' \deqn{\min_{\beta}\frac{1}{\zeta}\log\bigg(\sum_{g=1}^G \exp(-\zeta \hat V_g(\beta))\bigg) + \lambda  \Vert\beta\Vert_1, \quad \zeta > 0,\lambda \geq 0}
#'  for the setup described above. Note that
#' \deqn{\hat V_g(\beta):=\frac{1}{n}(2\beta^\top X^\top vec(Y_g)-\beta^\top X^\top X\beta),}
#'  is the empirical explained variance from \cite{Meinshausen and B{u}hlmann, 2015}.  See \cite{Lund et al., 2020} for more details and references.
#'
#'  For \eqn{d=1,2,3}, using only the marginal matrices \eqn{X_1,X_2,\ldots} (for \eqn{d=1} there is only one marginal), the function \code{softmaximin}
#' solves the soft maximin problem for a sequence of penalty parameters \eqn{\lambda_{max}>\ldots >\lambda_{min}>0}.
#'
#' Two optimization algorithms  are implemented, a non-monotone
#' proximal gradient (NPG) algorithm and a fast iterative soft thresholding algorithm (FISTA).
#' We note that this package also solves the problem above with the penalty given by the SCAD
#'  penalty, using the multiple step adaptive lasso procedure to loop over the proximal algorithm.
#'
#'
#' @return An object with S3 Class "SMMA".
#' \item{spec}{A string indicating the array dimension (1, 2 or 3) and the penalty.}
#' \item{coef}{A \eqn{p_1\cdots p_d \times} \code{nlambda} matrix containing the estimates of
#' the model coefficients (\code{beta}) for each \code{lambda}-value.}
#' \item{lambda}{A vector containing the sequence of penalty values used in the estimation procedure.}
#' \item{Obj}{A matrix containing the objective values for each iteration and each model.}
#' \item{df}{The number of nonzero coefficients for each value of \code{lambda}.}
#' \item{dimcoef}{A vector giving the dimension of the model coefficient array \eqn{\beta}.}
#' \item{dimobs}{A vector giving the dimension of the observation (response) array \code{Y}.}
#' \item{Iter}{A list with 4 items:
#' \code{bt_iter}  is total number of backtracking steps performed,
#' \code{bt_enter} is the number of times the backtracking is initiated,
#' and \code{iter_mat} is a vector containing the  number of  iterations for each \code{lambda} value
#' and  \code{iter} is total number of iterations i.e. \code{sum(Iter)}.}
#'
#' @author  Adam Lund
#'
#' Maintainer: Adam Lund, \email{adam.lund@@math.ku.dk}
#'
#' @references
#' Lund, A., S. W. Mogensen and N. R. Hansen (2020). Soft Maximin Estimation for Heterogeneous Array Data.
#' \emph{Preprint}.
#'
#' Meinshausen, N and P. B{u}hlmann (2015). Maximin effects in inhomogeneous large-scale data.
#' \emph{The Annals of Statistics}. 43, 4, 1801-1830. url = {https://doi.org/10.1214/15-AOS1325}.
#'
#' Currie, I. D., M. Durban, and P. H. C. Eilers (2006). Generalized linear
#' array models with applications to multidimensional smoothing.
#' \emph{Journal of the Royal Statistical Society. Series B}. 68, 259-280. url = {http://dx.doi.org/10.1111/j.1467-9868.2006.00543.x}.
#'
#'
#' @keywords package
#'
#' @examples
#'
#' ##size of example
#' n1 <- 65; n2 <- 26; n3 <- 13; p1 <- 13; p2 <- 5; p3 <- 4
#'
#' ##marginal design matrices (Kronecker components)
#' X1 <- matrix(rnorm(n1 * p1), n1, p1)
#' X2 <- matrix(rnorm(n2 * p2), n2, p2)
#' X3 <- matrix(rnorm(n3 * p3), n3, p3)
#' X <- list(X1, X2, X3)
#'
#' component <- rbinom(p1 * p2 * p3, 1, 0.1)
#' Beta1 <- array(rnorm(p1 * p2 * p3, 0, 0.1) + component, c(p1 , p2, p3))
#' mu1 <- RH(X3, RH(X2, RH(X1, Beta1)))
#' Y1 <- array(rnorm(n1 * n2 * n3), dim = c(n1, n2, n3)) + mu1
#' Beta2 <- array(rnorm(p1 * p2 * p3, 0, 0.1) + component, c(p1 , p2, p3))
#' mu2 <- RH(X3, RH(X2, RH(X1, Beta2)))
#' Y2 <- array(rnorm(n1 * n2 * n3), dim = c(n1, n2, n3)) + mu2
#' Beta3 <- array(rnorm(p1 * p2 * p3, 0, 0.1) + component, c(p1 , p2, p3))
#' mu3 <- RH(X3, RH(X2, RH(X1, Beta3)))
#' Y3 <- array(rnorm(n1 * n2 * n3), dim = c(n1, n2, n3)) + mu3
#' Beta4 <- array(rnorm(p1 * p2 * p3, 0, 0.1) + component, c(p1 , p2, p3))
#' mu4 <- RH(X3, RH(X2, RH(X1, Beta4)))
#' Y4 <- array(rnorm(n1 * n2 * n3), dim = c(n1, n2, n3)) + mu4
#' Beta5 <- array(rnorm(p1 * p2 * p3, 0, 0.1) + component, c(p1 , p2, p3))
#' mu5 <- RH(X3, RH(X2, RH(X1, Beta5)))
#' Y5 <- array(rnorm(n1 * n2 * n3), dim = c(n1, n2, n3)) + mu5
#'
#' Y <- array(NA, c(dim(Y1), 5))
#' Y[,,, 1] <- Y1; Y[,,, 2] <- Y2; Y[,,, 3] <- Y3; Y[,,, 4] <- Y4; Y[,,, 5] <- Y5;
#'
#' fit <- softmaximin(X, Y, zeta = 10, penalty = "lasso", alg = "npg")
#' Betafit <- fit$coef
#'
#' modelno <- 15
#' m <- min(Betafit[ , modelno], c(component))
#' M <- max(Betafit[ , modelno], c(component))
#' plot(c(component), type="l", ylim = c(m, M))
#' lines(Betafit[ , modelno], col = "red")
#'

#' @export
#' @useDynLib SMMA, .registration = TRUE
#' @importFrom Rcpp evalCpp
softmaximin <-function(X,
               Y,
               zeta,
               penalty = c("lasso", "scad"),
               alg = c("npg", "fista"),
               nlambda = 30,
               lambda.min.ratio = 1e-04,
               lambda = NULL,
               penalty.factor = NULL,
               reltol = 1e-05,
               maxiter = 15000,
               steps = 1,
               btmax = 100,
               c = 0.0001,
               tau = 2,
               M = 4,
               nu = 1,
               Lmin = 0,
               log = TRUE) {

dimglam <- length(X)

if (dimglam > 3){

stop(paste("the dimension of the model must be 1, 2 or 3!"))

}else if (dimglam == 1){

X[[2]] <- matrix(1, 1, 1)
X[[3]] <- matrix(1, 1, 1)

}else if (dimglam == 2){X[[3]] <- matrix(1, 1, 1)}

X1 <- X[[1]]
X2 <- X[[2]]
X3 <- X[[3]]

dimX <- rbind(dim(X1), dim(X2), dim(X3))

n1 <- dimX[1, 1]
n2 <- dimX[2, 1]
n3 <- dimX[3, 1]
p1 <- dimX[1, 2]
p2 <- dimX[2, 2]
p3 <- dimX[3, 2]
n <- prod(dimX[,1])
p <- prod(dimX[,2])
G <- dim(Y)[length(dim(Y))]

Z <- array(NA, c(n1, n2 * n3, G))

if(dimglam == 1){

  for(i in 1:G){Z[, , i] <- matrix(Y[,  i], n1, n2 * n3)}

}else if(dimglam == 2){

  for(i in 1:G){Z[, , i] <- matrix(Y[, ,  i], n1, n2 * n3)}

}else{
for(i in 1:G){Z[, , i] <- matrix(Y[, , , i], n1, n2 * n3)}
}

if(sum(alg == c("npg", "fista")) != 1){stop(paste("algorithm must be correctly specified"))}

if(alg == "npg"){alg <- 1}else{alg <- 0}

if(log == TRUE){ll <- 1}else{ll <- 0}

if(c <= 0){stop(paste("c must be strictly positive"))}

if(Lmin < 0){stop(paste("Lmin must be positive"))}

if(zeta <= 0){stop(paste("zeta must be strictly positive"))}

if(sum(penalty == c("lasso", "scad")) != 1){stop(paste("penalty must be correctly specified"))}

if(penalty == "lasso"){steps <- 1}

if(is.null(lambda)){

makelamb <- 1
lambda <- rep(NA, nlambda)

}else{
nlambda<-length(lambda)
makelamb <- 0

}

if(is.null(penalty.factor)){

penalty.factor <- matrix(1, p1, p2 * p3)

}else if(prod(dim(penalty.factor)) != p){

stop(
paste("number of elements in penalty.factor (", length(penalty.factor),") is not equal to the number of coefficients (", p,")", sep = "")
)

}else {

if(min(penalty.factor) < 0){stop(paste("penalty.factor must be positive"))}

penalty.factor <- matrix(penalty.factor, p1, p2 * p3)

}

res <- pga(X1, X2, X3,
           Z,
           penalty,
           zeta,
           c,
           lambda, nlambda, makelamb, lambda.min.ratio,
           penalty.factor,
           reltol,
           maxiter,
           steps,
           btmax,
           M,
           tau,
           nu,
           alg,
           ll,
           Lmin)

if(res$Stops[2] == 1){

warning(paste("program exit due to maximum number of inner iterations (",maxiter,") reached for model no ",res$endmodelno + 1,""))

}

if(res$Stop[3] == 1){

warning(paste("program exit due to maximum number of backtraking steps reached for model no ",res$endmodelno + 1,""))

}

endmodelno <- res$endmodelno #converged models since c++ is zero indexed
# Iter <- res$Iter
#
# maxiterpossible <- sum(Iter > 0)
# maxiterreached <- sum(Iter >= (maxiter - 1))
#
# if(maxiterreached > 0){
#
# warning(
# paste("maximum number of inner iterations (",maxiter,") reached ",maxiterreached," time(s) out of ",maxiterpossible," possible")
# )
#
#}

out <- list()

class(out) <- "SMMA"

out$spec <- paste("", dimglam,"-dimensional ", penalty," penalized model")
out$coef <- res$Beta[ , 1:endmodelno]
out$lambda <- res$lambda[1:endmodelno]
out$df <- res$df[1:endmodelno]
out$dimcoef <- c(p1, p2, p3)[1:dimglam]
out$dimobs <- c(n1, n2, n3)[1:dimglam]
out$Obj <- res$Obj
out$endno <- res$endmodelno
out$L <- res$L
out$L1 <- res$L1
out$sumsqdiff <- res$Sumsqdiff##remove
out$Delta <- res$Delta##remove
out$deltamax <- res$deltamax
out$BT <- res$BT

Iter <- list()
Iter$bt_enter <- res$btenter
Iter$bt_iter <- res$btiter
Iter$iter_mat <- res$Iter[1:endmodelno, ]
Iter$iter <- sum(Iter$iter_mat, na.rm = TRUE)

out$Iter <- Iter

return(out)

}

