#
# Description of this R script:
# R interface/wrapper for the Rcpp function brgdpg in the glamlasso package.
#
# Intended for use with R.
# Copyright (C) 2017 Adam Lund
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.If not, see <http://www.gnu.org/licenses/>
#

#' @name glamlassoRR
#' 
#' @title Penalized reduced rank regression  in a GLAM
#' 
#' @description Efficient design matrix free procedure for fitting large scale penalized reduced rank
#'  regressions in a 3-dimensional generalized linear array model. To obtain a factorization of the parameter array, 
#'  the \code{glamlassoRR} function performes a block relaxation scheme within the gdpg algorithm, see \cite{Lund et al., 2017}. 
#'
#' @usage  glamlassoRR(X, 
#'             Y, 
#'             Z = NULL,
#'             family = "gaussian",
#'             penalty = "lasso",
#'             intercept = FALSE,
#'             weights = NULL,
#'             thetainit = NULL,
#'             alphainit = NULL,
#'             nlambda = 100,
#'             lambdaminratio = 1e-04,
#'             lambda = NULL,
#'             penaltyfactor = NULL,
#'             penaltyfactoralpha = NULL,
#'             reltolinner = 1e-07,
#'             reltolouter = 1e-04,
#'             reltolalt = 1e-04,
#'             maxiter = 15000,
#'             steps = 1,
#'             maxiterinner = 3000,
#'             maxiterouter = 25,
#'             maxalt = 10,
#'             btinnermax = 100,
#'             btoutermax  = 100,
#'             iwls = "exact",
#'             nu = 1)
#'
#' @param X A list containing the 3 tensor components of the tensor design matrix. These are  matrices of sizes \eqn{n_i   \times p_i}.
#' @param Y The response values, an array of size \eqn{n_1 \times n_2\times n_3}. For option 
#' \code{family = "binomial"} this array must contain the proportion of successes and the 
#' number of trials is then specified as \code{weights} (see below).
#' @param Z The non tensor structrured part of the design matrix. A matrix of size \eqn{n_1 n_2 n_3\times q}. 
#' Is set to \code{NULL} as default. 
#' @param family A string specifying the model family (essentially the response distribution). Possible values 
#' are \code{"gaussian", "binomial", "poisson", "gamma"}.
#' @param penalty A string specifying the penalty. Possible values are \code{"lasso", "scad"}.
#' @param intercept  Logical variable indicating if the model includes an intercept.  When \code{intercept = TRUE} the first 
#' coulmn in the non-tensor design component \code{Z} is all 1s. Default is \code{FALSE}. 
#' @param weights Observation weights, an array of size \eqn{n_1 \times \cdots \times n_d}. For option 
#' \code{family = "binomial"} this array must contain the number of trials and must be provided.
#' @param thetainit A list (length 2) containing the initial parameter values for each of the parameter factors. 
#' Default is NULL in which case all parameters are initialized at 0.01.
#' @param alphainit A \eqn{q\times 1} vector containing the initial parameter values for the non-tensor parameter. 
#' Default is NULL in which case all parameters are initialized at 0.
#' @param nlambda The number of \code{lambda} values.
#' @param lambdaminratio The smallest value for \code{lambda}, given as a fraction of 
#' \eqn{\lambda_{max}}; the (data derived) smallest value for which all coefficients are zero.
#' @param lambda The sequence of penalty parameters for the regularization path.
#' @param penaltyfactor A list of length two containing an array of size \eqn{p_1 \times  p_2} and a \eqn{p_3 \times  1} vector.  
#' Multiplied  with each element in \code{lambda} to allow differential shrinkage on the (tensor) coefficients blocks. 
#' @param penaltyfactoralpha A \eqn{q \times 1} vector multiplied with each element in \code{lambda} to allow differential shrinkage on the non-tensor coefficients. 
#' @param reltolinner The convergence tolerance for the inner loop
#' @param reltolouter The convergence tolerance for the outer loop.
#' @param reltolalt The convergence tolerance for the alternation loop over the two parameter blocks.
#' @param maxiter The maximum number of inner iterations allowed for each \code{lambda}
#' value, when  summing over all outer iterations for said \code{lambda}.
#' @param steps The number of steps used in the multi-step adaptive lasso algorithm for non-convex penalties. Automatically set to 1 when \code{penalty = "lasso"}.
#' @param maxiterinner The maximum number of inner iterations allowed for each outer iteration.
#' @param maxiterouter The maximum number of outer iterations allowed for each lambda.
#' @param maxalt The maximum number of  alternations over parameter blocks.
#' @param btinnermax Maximum number of backtracking steps allowed in each inner iteration. Default is \code{btinnermax = 100}.
#' @param btoutermax Maximum number of backtracking steps allowed in each outer iteration. Default is \code{btoutermax = 100}.
#' @param iwls A string indicating whether to use the exact iwls weight matrix or use a tensor structured approximation to it.
#' @param nu A number between 0 and 1 that controls the step size \eqn{\delta} in the proximal algorithm (inner loop) by 
#' scaling the upper bound \eqn{\hat{L}_h} on the Lipschitz constant \eqn{L_h} (see \cite{Lund et al., 2017}). 
#' For \code{nu = 1} backtracking never occurs and the proximal step size is always \eqn{\delta = 1 / \hat{L}_h}. 
#' For \code{nu = 0} backtracking always occurs and the proximal step size is initially \eqn{\delta = 1}. 
#' For \code{0 < nu < 1} the proximal step size is initially \eqn{\delta = 1/(\nu\hat{L}_h)} and backtracking 
#' is only employed if the objective function does not decrease. A \code{nu} close  to 0 gives large step 
#' sizes and presumably more backtracking in the inner loop. The default is \code{nu = 1} and the option is only 
#' used if \code{iwls = "exact"}.
#' 
#' @details Given the setting from \code{\link{glamlasso}} we  place a reduced rank
#'  restriction on the \eqn{p_1\times p_2\times p _3} parameter array \eqn{\Theta} given by
#' \deqn{\Theta=(\Theta_{i,j,k})_{i,j,k} = (\gamma_{k}\beta_{i,j})_{i,j,k}, \ \ \ \gamma_k,\beta_{i,j}\in \mathcal{R}.}  
#'  The  \code{glamlassoRR} function  solves the PMLE problem by combining a block relaxation scheme with the gdpg algorithm. This scheme alternates between  optimizing over the first 
#' parameter block \eqn{\beta=(\beta_{i,j})_{i,j}} and  the second block \eqn{\gamma=(\gamma_k)_k} while fixing the second resp. first block. We note that the 
#' individual parameter blocks are only identified up to a multiplicative constant.
#' 
#' @return An object with S3 Class "glamlasso". 
#' \item{spec}{A string indicating the model family and the penalty.}  
#' \item{coef12}{A \eqn{p_1 p_2 \times} \code{nlambda} matrix containing the estimates of 
#' the first model coefficient factor  (\eqn{\beta}) for each \code{lambda}-value.}
#' \item{coef3}{A \eqn{p_3 \times} \code{nlambda} matrix containing the estimates of 
#' the second model coefficient factor  (\eqn{\gamma}) for each \code{lambda}-value.}
#' \item{alpha}{A \eqn{q \times} \code{nlambda} matrix containing the estimates of 
#' the parameters for the non tensor structured part of the model (\code{alpha}) for each \code{lambda}-value.}. 
#' If \code{intercept = TRUE} the first row contains the intercept estimate for each \code{lambda}-value.}.
#' \item{lambda}{A vector containing the sequence of penalty values used in the estimation procedure.}
#' \item{df}{The number of nonzero coefficients for each value of \code{lambda}.}	
#' \item{dimcoef}{A vector giving the dimension of the model coefficient array \eqn{\beta}.}
#' \item{dimobs}{A vector giving the dimension of the observation (response) array \code{Y}.}
#' \item{Iter}{A list with 4 items:  
#' \code{bt_iter_inner}  is total number of backtracking steps performed in the inner loop,
#' \code{bt_enter_inner} is the number of times the backtracking is initiated in the inner loop,
#' \code{bt_iter_outer} is total number of backtracking steps performed in the outer loop,
#' and \code{iter_mat} is a \code{nlambda} \eqn{\times} \code{maxiterouter} matrix containing the  number of 
#' inner iterations for each \code{lambda} value and each outer iteration and  \code{iter} is total number of iterations i.e. \code{sum(Iter)}.}  
#'
#' @author Adam Lund
#' 
#' Maintainer: Adam Lund, \email{adam.lund@@math.ku.dk}
#' 
#' @references 
#' Lund, A. and N. R. Hansen (2017). Sparse Network  Estimation for  Dynamical Spatio-temporal Array Models. 
#'  \emph{ArXiv}. 
#'  
#' @examples 
#' \dontrun{
#' ##size of example 
#' n1 <- 65; n2 <- 26; n3 <- 13; p1 <- 12; p2 <- 6; p3 <- 4
#' 
#' ##marginal design matrices (tensor components)
#' X1 <- matrix(rnorm(n1 * p1), n1, p1) 
#' X2 <- matrix(rnorm(n2 * p2), n2, p2) 
#' X3 <- matrix(rnorm(n3 * p3), n3, p3) 
#' X <- list(X1, X2, X3)
  
#' Beta12 <- matrix(rnorm(p1 * p2), p1, p2) * matrix(rbinom(p1 * p2, 1, 0.5), p1, p2)
#' Beta3 <- matrix(rnorm(p3) * rbinom(p3, 1, 0.5), p3, 1)
#' Beta <- outer(Beta12, c(Beta3))
#' Mu <- RH(X3, RH(X2, RH(X1, Beta)))
#' Y <- array(rnorm(n, Mu), dim = c(n1, n2, n3))  
#' 
#' system.time(fit <- glamlassoRR(X, Y))
#' 
#' modelno  <- length(fit$lambda)
#' par(mfrow = c(1, 3))
#' plot(c(Beta), type = "h")
#' points(c(Beta))
#' lines(c(outer(fit$coef12[, modelno], c(fit$coef3[, modelno]))), col = "red", type = "h")
#' plot(c(Beta12), ylim = range(Beta12, fit$coef12[, modelno]), type = "h")
#' points(c(Beta12))
#' lines(fit$coef12[, modelno], col = "red", type = "h")
#' plot(c(Beta3), ylim = range(Beta3, fit$coef3[, modelno]), type = "h")
#' points(c(Beta3))
#' lines(fit$coef3[, modelno], col = "red", type = "h")
#' 
#' ###with non tensor design component Z
#' q <- 5
#' alpha <- matrix(rnorm(q)) * rbinom(q, 1, 0.5)
#' Z <- matrix(rnorm(n1 * n2 * n3 * q), n1 * n2 * n3, q) 
#' Y <- array(rnorm(n1 * n2 * n3, Mu + array(Z %*% alpha, c(n1, n2, n3))), c(n1, n2, n3))
#' system.time(fit <- glamlassoRR(X, Y, Z))
#' 
#' modelno <- length(fit$lambda)
#' par(mfrow = c(2, 2))
#' plot(c(Beta), type = "h")
#' points(c(Beta))
#' lines(c(outer(fit$coef12[, modelno], c(fit$coef3[, modelno]))), col = "red", type = "h")
#' plot(c(Beta12), ylim = range(Beta12,fit$coef12[, modelno]), type = "h")
#' points(c(Beta12))
#' lines(fit$coef12[, modelno], col = "red", type = "h")
#' plot(c(Beta3), ylim = range(Beta3, fit$coef3[, modelno]), type = "h")
#' points(c(Beta3))
#' lines(fit$coef3[, modelno], col = "red", type = "h")
#' plot(c(alpha), ylim = range(alpha, fit$alpha[, modelno]), type = "h")
#' points(c(alpha))
#' lines(fit$alpha[, modelno], col = "red", type = "h")
#' 
#' ################ poisson example
#' Beta12 <- matrix(rnorm(p1 * p2, 0, 0.5), p1, p2) * matrix(rbinom(p1 * p2, 1, 0.1), p1, p2)
#' Beta3 <-  matrix(rnorm(p3, 0, 0.5) * rbinom(p3, 1, 0.5), p3, 1)
#' Beta <- outer(Beta12, c(Beta3))
#' Mu <- RH(X3, RH(X2, RH(X1, Beta)))
#' Y <- array(rpois(n1 * n2 * n3, exp(Mu)), dim = c(n1, n2, n3))
#' system.time(fit <- glamlassoRR(X, Y, family = "poisson"))
#' 
#' modelno <- length(fit$lambda)
#' par(mfrow = c(1, 3))
#' plot(c(Beta), type = "h")
#' points(c(Beta))
#' lines(c(outer(fit$coef12[, modelno], c(fit$coef3[, modelno]))), col = "red", type = "h")
#' plot(c(Beta12), ylim = range(Beta12, fit$coef12[, modelno]), type = "h")
#' points(c(Beta12))
#' lines(fit$coef12[, modelno], col = "red", type = "h")
#' plot(c(Beta3), ylim = range(Beta3, fit$coef3[, modelno]), type = "h")
#' points(c(Beta3))
#' lines(fit$coef3[, modelno], col = "red", type = "h")
#' 
#'}

glamlassoRR <- function(X,
                        Y, 
                        Z = NULL, 
                        family = "gaussian",
                        penalty = "lasso",
                        intercept = FALSE,
                        weights = NULL,
                        thetainit =  NULL,
                        alphainit = NULL,
                        nlambda = 100,
                        lambdaminratio = 0.0001,
                        lambda = NULL,
                        penaltyfactor = NULL,
                        penaltyfactoralpha = NULL,
                        reltolinner = 1e-07,
                        reltolouter = 1e-04,
                        reltolalt = 1e-04,
                        maxiter = 15000,
                        steps = 1,
                        maxiterinner = 3000,
                        maxiterouter = 25,
                        maxalt = 10,
                        btinnermax = 100,
                        btoutermax  = 100,
                        iwls = "exact",
                        nu = 1) { 

##get dimensions of problem
dimglam <- length(X)

if (dimglam < 3 || dimglam > 3){

stop(paste("the dimension of the GLAM must be  3!"))

}#else if(dimglam == 2){X[[3]] <- matrix(1, 1, 1)} does the rcpp code work for d=2? or is it hardcoded to d=3?

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

##restrictions/initial par.
if(is.null(thetainit)){
  
thetainit <- list()
thetainit[[1]] <- matrix(0.01, p1, p2)
thetainit[[2]] <- matrix(0.01, p3, 1)

}else{

if(length(thetainit) == 2){

if(max(abs(c(dim(thetainit[[1]]), dim(thetainit[[2]])[1])) - c(p1, p2, p3)) != 0){

stop("dimensions of parameter array (",c(dim(thetainit[[1]]), dim(thetainit[[2]])[1]),") incorrect")

}

}else{stop(paste("initial parameter values(parameter restrictions) must be (correctly) specified")) }

}

Yrot <- aperm(Y, c(3, 1, 2))
Yrot <- matrix(Yrot, dim(Yrot)[1], dim(Yrot)[2] * dim(Yrot)[3])

####check on weights 
if(family == "binomial" & is.null(weights)){stop(paste("for binomial model number of trials (weights) must be specified"))}

if(is.null(weights)){
  
weightedgaussian <- 0  
weights <- matrix(1, n1, n2 * n3)
weightsrot <- matrix(1, n3, n1 * n2) 
  
}else{
  
if(min(weights) < 0){stop(paste("only positive weights allowed"))}    
  
weightedgaussian <- 1
weightsrot <- aperm(weights, c(3, 1, 2))
weights <- matrix(weights, n1, n2 * n3)
weightsrot <- matrix(weightsrot, n3, n1 * n2) #right??
  
}

## nontensor component Z and intercept
if(is.null(Z) == TRUE & intercept == FALSE){
  
nonten <- 0
q <- 1
Z <- matrix(0, 1, q)
Zrot <- Z
  
}else if(is.null(Z) == TRUE & intercept == TRUE){
  
nonten <- 1
q <- 1
Z <- matrix(1, n, q)
Zrot <- Z
  
}else if(is.null(Z) ==  FALSE & intercept == FALSE){
  
nonten <- 1
q <- dim(Z)[2]
Zrot <- Z[c(aperm(array(1:(n1 * n2 * n3), c(n1, n2, n3)), c(3, 1, 2))), ]
  
}else{
  
nonten <- 1  
Z <- cbind(matrix(1, n, 1), Z)
q <- dim(Z)[2]
Zrot <- Z[c(aperm(array(1:(n1 * n2 * n3), c(n1, n2, n3)), c(3, 1, 2))), ]
  
}

if(is.null(alphainit)){alphainit <- matrix(0, q, 1)}

####reshape Y into matrix
Y <- matrix(Y, dim(Y)[1], dim(Y)[2] * dim(Y)[3])

####check on family
if(sum(family == c("gaussian", "binomial", "poisson", "gamma")) != 1){stop(paste("family must be correctly specified"))}

####check on penalty
if(sum(penalty == c("lasso", "scad")) != 1){stop(paste("penalty must be correctly specified"))}

if(penalty == "lasso"){steps <- 1}

####check on lambda 
if(is.null(lambda)){

makelamb <- 1
lambda <- rep(NA, nlambda)

}else{

if(length(lambda) != nlambda){

#warning(paste("number of elements in lambda (", length(lambda),") is not equal to default nlambda (", nlambda,")", sep = ""))
nlambda <- length(lambda)

}

makelamb <- 0

}

####check on penaltyfactor
if(length(thetainit) == 2){

if(is.null(penaltyfactor[1]) & is.null(penaltyfactor[2])){
penaltyfactor[[1]] <- matrix(1, p1, p2)
penaltyfactor[[2]] <- matrix(1, p3, 1)

}else if(is.null(penaltyfactor[1])) {

penaltyfactor[[1]] <- matrix(1, p1, p2)

}else if(is.null(penaltyfactor[2])){

penaltyfactor[[2]] <- matrix(1, p3, 1)

}

if(length(penaltyfactor) != length(thetainit)){

stop(paste("parameter list length (",length(thetainit),") not equal to penaltyfactor length (",length(penaltyfactor),")"))

}

if(max(abs(c(dim(penaltyfactor[[1]]), dim(penaltyfactor[[2]])[1])) - c(p1, p2, p3)) != 0){

stop("dimensions of penalty array incorrect")

}}else{stop()}


##alpha
if(is.null(penaltyfactoralpha)){
  
penaltyfactoralpha <- matrix(1, q, 1)
  
}else if(prod(dim(penaltyfactoralpha)) != q){
  
stop(paste("number of elements in penaltyfactoralpha (", length(penaltyfactoralpha),") is not equal to the number of coefficients (", q,")", sep = ""))
  
}else {if(min(penaltyfactoralpha) < 0){stop(paste("penaltyfactoralpha must be positive"))}} 

####check on iwls
if(is.null(iwls)){ 
  
iwls <- "exact"
  
}else{if(sum(iwls == c("exact", "identity", "kron1", "kron2")) != 1){stop(paste("iwls must be correctly specified"))}}

####check on nu
if(nu < 0 || nu > 1){stop(paste("nu must be between 0 and 1"))}

##run br-gdpg algorithm
res <- gdpg(X1, X2, X3,
            Z, Zrot,##non ten design
            Y, Yrot,
            matrix(0, n1, n2 * n3), #V
            weights, weightsrot,
            matrix(outer(thetainit[[1]], c(thetainit[[2]])), p1, p2 * p3), thetainit[[1]], thetainit[[2]],
            alphainit, ## non ten par
            family,  
            penalty,
            nonten, #### non ten ind
            iwls, 
            nu,
            lambda, makelamb, nlambda, lambdaminratio,
            matrix(outer(penaltyfactor[[1]], c(penaltyfactor[[2]])), p1, p2 * p3),  penaltyfactor[[1]], penaltyfactor[[2]],
            penaltyfactoralpha, ##non ten penaltyfactor
            reltolinner, 
            reltolouter, 
            reltolalt,
            maxiter,
            steps,
            maxiterinner,
            maxiterouter,
            maxalt,
            btinnermax,
            btoutermax,
            weightedgaussian, 0, 1, n1, n2, n3)

endmodelno <- res$endmodelno + 1
Iter <- res$Iter

out <- list()

class(out) <- "glamlasso"
out$spec <- paste("", dimglam,"-dimensional ", penalty," penalized reduced rank", family," GLAM")
out$coef <- res$Beta[ , 1:endmodelno]
out$coef12 <- res$Beta12[ , 1:endmodelno]
out$coef3 <- res$Beta3[ , 1:endmodelno]
out$alpha <- res$alpha[, 1:endmodelno]
out$lambda <- res$lambda[1:endmodelno]
out$df <- res$df[1:endmodelno]
out$dimcoef <- c(p1, p2, p3)[1:dimglam]
out$dimobs <- c(n1, n2, n3)[1:dimglam]

Iter <- list()
Iter$bt_enter_inner <- res$btenterprox
Iter$bt_iter_inner <- res$btiterprox
Iter$bt_iter_outer <- res$btiternewt
Iter$iter_mat <- res$Iter[1:endmodelno, ]
Iter$iter <- sum(Iter$iter_mat, na.rm = TRUE)

out$Iter <- Iter
out$Lambdas <- res$Lambdas
return(out)

}