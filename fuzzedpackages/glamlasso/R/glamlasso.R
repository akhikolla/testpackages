#
#     Description of this R script:
#     R interface/wrapper for the Rcpp function gdpg in the glamlasso package.
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

#' @name glamlasso
#' 
#' @title Penalization in Large Scale Generalized Linear Array Models
#' 
#' @description  Efficient design matrix free procedure for fitting large scale penalized  2 or 3-dimensional
#' generalized linear array models (GLAM). It is also possible to fit an additional non-tensor structured component 
#' - e.g an intercept - however this can reduce the computational efficiency of the procedure substanstially. 
#' Currently the LASSO penalty and the SCAD penalty are both implemented. Furthermore,
#' the Gaussian model with identity link,  the Binomial model with logit link, the Poisson model
#' with log link and the Gamma model with log link is currently implemented. The underlying algorithm combines gradient descent and proximal gradient (gdpg algorithm), see  \cite{Lund et al., 2017}. 
#'  
#' @usage  glamlasso(X, 
#'           Y, 
#'           Z = NULL,
#'           family = "gaussian",
#'           penalty = "lasso",
#'           intercept = FALSE,
#'           weights = NULL,
#'           thetainit = NULL,
#'           alphainit = NULL,
#'           nlambda = 100,
#'           lambdaminratio = 1e-04,
#'           lambda = NULL,
#'           penaltyfactor = NULL,
#'           penaltyfactoralpha = NULL,
#'           reltolinner = 1e-07,
#'           reltolouter = 1e-04,
#'           maxiter = 15000,
#'           steps = 1,
#'           maxiterinner = 3000,
#'           maxiterouter = 25,
#'           btinnermax = 100,
#'           btoutermax = 100,
#'           iwls = "exact",
#'           nu = 1)
#'            
#' @param X A list containing the tensor components (2 or 3) of the tensor design matrix.
#'  These are  matrices of sizes \eqn{n_i   \times p_i}.
#' @param Y The response values, an array of size \eqn{n_1 \times\cdots\times n_d}. For option 
#' \code{family = "binomial"} this array must contain the proportion of successes and the 
#' number of trials is then specified as \code{weights} (see below).
#' @param Z The non tensor structrured part of the design matrix. A matrix of size \eqn{n_1 \cdots n_d\times q}. 
#' Is set to \code{NULL} as default. 
#' @param family A string specifying the model family (essentially the response distribution). Possible values 
#' are \code{"gaussian", "binomial", "poisson", "gamma"}.
#' @param penalty A string specifying the penalty. Possible values 
#' are \code{"lasso", "scad"}.
#' @param intercept  Logical variable indicating if the model includes an intercept. 
#' When \code{intercept = TRUE} the first coulmn in the non-tensor design component \code{Z} is all 1s.
#'  Default is \code{FALSE}. 
#' @param weights Observation weights, an array of size \eqn{n_1 \times \cdots \times n_d}. For option 
#' \code{family = "binomial"} this array must contain the number of trials and must be provided.
#' @param thetainit The initial parameter values. Default is NULL in which case all parameters are initialized at zero.
#' @param alphainit A \eqn{q\times 1} vector containing the initial parameter values for the non-tensor parameter. 
#'  Default is NULL in which case all parameters are initialized at 0.
#' @param nlambda The number of \code{lambda} values.
#' @param lambdaminratio The smallest value for \code{lambda}, given as a fraction of 
#' \eqn{\lambda_{max}}; the (data derived) smallest value for which all coefficients are zero.
#' @param lambda The sequence of penalty parameters for the regularization path.
#' @param penaltyfactor An array of size \eqn{p_1 \times \cdots \times p_d}. Is multiplied 
#' with each element in \code{lambda} to allow differential shrinkage on the coefficients. 
#' @param penaltyfactoralpha A \eqn{q \times 1} vector multiplied with each element in \code{lambda} to allow differential 
#' shrinkage on the non-tensor coefficients. 
#' @param reltolinner The convergence tolerance for the inner loop
#' @param reltolouter The convergence tolerance for the outer loop.
#' @param maxiter The maximum number of inner iterations allowed for each \code{lambda}
#' value, when  summing over all outer iterations for said \code{lambda}.
#' @param steps The number of steps used in the multi-step adaptive lasso algorithm for non-convex penalties. Automatically set to 1 when \code{penalty = "lasso"}.
#' @param maxiterinner The maximum number of inner iterations allowed for each outer iteration.
#' @param maxiterouter The maximum number of outer iterations allowed for each lambda.
#' @param btinnermax Maximum number of backtracking steps allowed in each inner iteration. Default is \code{btinnermax = 100}.
#' @param btoutermax Maximum number of backtracking steps allowed in each outer iteration. Default is \code{btoutermax = 100}.
#' @param iwls A string indicating whether to use the exact iwls weight matrix or use a kronecker structured approximation to it.
#' @param nu A number between 0 and 1 that controls the step size \eqn{\delta} in the proximal algorithm (inner loop) by 
#' scaling the upper bound \eqn{\hat{L}_h} on the Lipschitz constant \eqn{L_h} (see \cite{Lund et al., 2017}). 
#' For \code{nu = 1} backtracking never occurs and the proximal step size is always \eqn{\delta = 1 / \hat{L}_h}. 
#' For \code{nu = 0} backtracking always occurs and the proximal step size is initially \eqn{\delta = 1}. 
#' For \code{0 < nu < 1} the proximal step size is initially \eqn{\delta = 1/(\nu\hat{L}_h)} and backtracking 
#' is only employed if the objective function does not decrease. A \code{nu} close  to 0 gives large step 
#' sizes and presumably more backtracking in the inner loop. The default is \code{nu = 1} and the option is only 
#' used if \code{iwls = "exact"}.
#' 
#' @details Consider a (two component) generalized linear model (GLM) 
#' \deqn{g(m) = X\theta + Z\alpha =: \eta.}
#' Here \eqn{g} is a link function, \eqn{m} is a \eqn{n\times 1} vector containing the mean of the 
#' response variable  \eqn{Y}, \eqn{Z} is a \eqn{n\times q} matrix and  \eqn{X} a \eqn{n\times p} matrix  with tensor structure 
#'  \deqn{X =  X_d\otimes\ldots\otimes X_1,}
#'   where \eqn{X_1,\ldots,X_d} are the marginal \eqn{n_i\times p_i} design matrices (tensor factors) such that
#'  \eqn{p = p_1\cdots p_d} and  \eqn{n=n_1\cdots n_d}. Then \eqn{\theta} is the \eqn{p\times 1} parameter associated with the tensor component
#'   \eqn{X} and \eqn{\alpha} the \eqn{q\times 1} parameter associated with the non-tensor component \eqn{Z}, e.g. the intercept. 
#'  
#' The related log-likelihood is a function of  \eqn{\tilde\theta:=(\theta,\alpha)} through the linear predictor \eqn{\eta} i.e. \eqn{\tilde\theta \mapsto l(\eta(\tilde\theta))}.
#' In the usual exponential family framework this can be expressed as
#' \deqn{l(\eta(\tilde\theta)) = \sum_{i = 1}^n a_i \frac{y_i \vartheta(\eta_i(\tilde\theta)) - b(\vartheta(\eta_i(\tilde\theta)))}{\psi}+c(y_i,\psi)} 
#' where \eqn{\vartheta}, the canonical parameter map,  is linked to the  linear predictor via the identity
#'  \eqn{\eta(\tilde\theta) = g(b'(\vartheta))} with \eqn{b} the cumulant function. Here \eqn{a_i \ge 0, i = 1,\ldots,n} are observation weights and
#'  \eqn{\psi} is the dispersion parameter.
#' 
#' By ignoring the non-tensor component \eqn{Z} (assume \eqn{\alpha = 0}) we can use the generalized linear array model (GLAM) framework to write the model equation as
#'  \deqn{g(M) = \rho(X_d,\rho(X_{d-1},\ldots,\rho(X_1,\Theta))),}
#' where \eqn{\rho} is the so called rotated \eqn{H}-transform and \eqn{M} and \eqn{\Theta} 
#' are the  array versions of \eqn{m} and \eqn{\theta} respectively. See \cite{Currie et al., 2006} for more details.
#'         
#' For \eqn{d = 3} or \eqn{d = 2}, using only the marginal matrices \eqn{X_1,X_2,\ldots}, the function \code{glamlasso} solves the penalized estimation problem 
#' \deqn{\min_{\theta} -l(\eta(\theta)) + \lambda J (\theta),} 
#' for \eqn{J} either the LASSO or SCAD penalty function,  in the GLAM setup for a sequence of penalty parameters \eqn{\lambda>0}. The underlying algorithm is based on an outer 
#' gradient descent loop and an inner proximal gradient based loop. We note that if \eqn{J} is not 
#' convex, as with the SCAD penalty, we use the multiple step adaptive lasso procedure to loop over the inner proximal algorithm, see \cite{Lund et al., 2017} for more details.
#'   
#' Furthermore, the function \code{glamlasso} also solves the penalized estimation problem for a model that  includes  a non-tensor component \eqn{Z}, e.g. an intercept. However, 
#' not without incurring a  potentially substantial computational cost. Especially it is not advisable to inlcude a very large non-tensor component in the model (large \eqn{q}) 
#' and even  adding an intecept to the model (\eqn{q=1}) will result in a reduction of computational efficiency.     
#'   
#' @return An object with S3 Class "glamlasso". 
#' \item{spec}{A string indicating the GLAM dimension (\eqn{d = 2, 3}), the model family and the penalty.}  
#' \item{beta}{A \eqn{p_1\cdots p_d \times} \code{nlambda} matrix containing the estimates of 
#' the parameters for the tensor structured part of the model (\code{beta}) for each \code{lambda}-value.}
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
#' @author  Adam Lund
#' 
#' Maintainer: Adam Lund, \email{adam.lund@@math.ku.dk}
#' 
#' @references 
#' Lund, A., M. Vincent, and N. R. Hansen (2017). Penalized estimation in 
#' large-scale generalized linear array models. 
#' \emph{Journal of Computational and Graphical Statistics}, 26, 3, 709-724.  url = {https://doi.org/10.1080/10618600.2017.1279548}.
#' 
#' Currie, I. D., M. Durban, and P. H. C. Eilers (2006). Generalized linear
#' array models with applications to multidimensional smoothing. 
#' \emph{Journal of the Royal Statistical Society. Series B}. 68, 259-280. url = {http://dx.doi.org/10.1111/j.1467-9868.2006.00543.x}.
#' 
#' @keywords package 
#'
#' @examples 
#' ##size of example 
#' n1 <- 65; n2 <- 26; n3 <- 13; p1 <- 12; p2 <- 6; p3 <- 4
#' 
#' ##marginal design matrices (tensor components)
#' X1 <- matrix(rnorm(n1 * p1), n1, p1) 
#' X2 <- matrix(rnorm(n2 * p2), n2, p2) 
#' X3 <- matrix(rnorm(n3 * p3), n3, p3) 
#' X <- list(X1, X2, X3)
#' 
#' ##############gaussian example 
#' Beta <- array(rnorm(p1 * p2 * p3) * rbinom(p1 * p2 * p3, 1, 0.1), c(p1 , p2, p3))
#' Mu <- RH(X3, RH(X2, RH(X1, Beta)))
#' Y <- array(rnorm(n1 * n2 * n3, Mu), c(n1, n2, n3))
#' 
#' system.time(fit <- glamlasso(X, Y))
#' 
#' modelno <- length(fit$lambda)
#' plot(c(Beta), type = "h", ylim = range(Beta, fit$coef[, modelno]))
#' points(c(Beta))
#' lines(fit$coef[ , modelno], col = "red", type = "h")
#' 
#' \dontrun{
#' ###with non tensor design component Z
#' q <- 5
#' alpha <- matrix(rnorm(q)) * rbinom(q, 1, 0.5)
#' Z <- matrix(rnorm(n1 * n2 * n3 * q), n1 * n2 *n3, q) 
#' Y <- array(rnorm(n1 * n2 * n3, Mu + array(Z %*% alpha, c(n1, n2, n3))), c(n1, n2, n3))
#' system.time(fit <- glamlasso(X, Y, Z))
#' 
#' modelno <- length(fit$lambda)
#' par(mfrow = c(1, 2))
#' plot(c(Beta), type = "l", ylim = range(Beta, fit$coef[, modelno]))
#' points(c(Beta))
#' lines(fit$coef[ , modelno], col = "red")
#' plot(c(alpha), type = "h", ylim = range(Beta, fit$alpha[, modelno]))
#' points(c(alpha))
#' lines(fit$alpha[ , modelno], col = "red", type = "h")
#' 
#' ################ poisson example
#' Beta <- array(rnorm(p1 * p2 * p3, 0, 0.1) * rbinom(p1 * p2 * p3, 1, 0.1), c(p1 , p2, p3))
#' Mu <- RH(X3, RH(X2, RH(X1, Beta)))
#' Y <- array(rpois(n1 * n2 * n3, exp(Mu)), dim = c(n1, n2, n3))
#' system.time(fit <- glamlasso(X, Y, family = "poisson", nu = 0.1))
#' 
#' modelno <- length(fit$lambda)
#' plot(c(Beta), type = "h", ylim = range(Beta, fit$coef[, modelno]))
#' points(c(Beta))
#' lines(fit$coef[ , modelno], col = "red", type = "h")
#' 
#' }

glamlasso <-function(X,
                     Y, 
                     Z = NULL, 
                     family = "gaussian",
                     penalty = "lasso",
                     intercept = FALSE,
                     weights = NULL,
                     thetainit = NULL,
                     alphainit = NULL,
                     nlambda = 100,
                     lambdaminratio = 0.0001,
                     lambda = NULL,
                     penaltyfactor = NULL,
                     penaltyfactoralpha = NULL,
                     reltolinner = 1e-07,
                     reltolouter = 1e-04,
                     maxiter = 15000,
                     steps = 1,
                     maxiterinner = 3000,
                     maxiterouter = 25,
                     btinnermax = 100,
                     btoutermax = 100,
                     iwls = "exact",
                     nu = 1) {

##get dimensions of problem
dimglam <- length(X)

if (dimglam < 2 || dimglam > 3){
  
  stop(paste("the dimension of the GLAM must be 2 or 3!"))
  
}else if(dimglam == 2){X[[3]] <- matrix(1, 1, 1)} 

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

## nontensor component Z and intercept
if(is.null(Z) == TRUE & intercept == FALSE){

nonten <- 0
q <- 1
Z <- matrix(0, 1, q)

}else if(is.null(Z) == TRUE & intercept == TRUE){
  
nonten <- 1
q <- 1
Z <- matrix(1, n, q)

}else if(is.null(Z) ==  FALSE & intercept == FALSE){
  
nonten <- 1
q <- dim(Z)[2]

}else{
  
nonten <- 1  
Z <- cbind(matrix(1, n, 1), Z)
q <- dim(Z)[2]

}

####reshape Y into matrix
Y <- matrix(Y, n1, n2 * n3)

####check on family
if(sum(family == c("gaussian", "binomial", "poisson", "gamma")) != 1){

stop(paste("family must be correctly specified"))

}
 
####check on penalty
if(sum(penalty == c("lasso", "scad")) != 1){stop(paste("penalty must be correctly specified"))}

if(penalty == "lasso"){steps <- 1}

####check on weights 
if(family == "binomial" & is.null(weights)){stop(paste("for binomial model number of trials (weights) must be specified"))}

if(is.null(weights)){
  
weightedgaussian <- 0  
weights <- matrix(1, n1, n2 * n3)

}else{

if(min(weights) < 0){stop(paste("only positive weights allowed"))}    

weightedgaussian <- 1
weights <- matrix(weights, n1, n2 * n3)

}


####check on initial values 
if(is.null(thetainit)){thetainit <- matrix(0, p1, p2 * p3)}else{thetainit <- matrix(thetainit, p1, p2 * p3)}

if(is.null(alphainit)){alphainit <- matrix(0, q, 1)}

####check on lambda 
if(is.null(lambda)){
  
  makelamb <- 1
  lambda <- rep(NA, nlambda)
  
}else{
  
  if(length(lambda) != nlambda){
    
 #   warning(paste("number of elements in lambda (", length(lambda),") is not equal to default nlambda (", nlambda,")", sep = ""))
    nlambda <- length(lambda)
  }
  
  makelamb <- 0
  
}

####check on penaltyfactor
if(is.null(penaltyfactor)){
  
penaltyfactor <- matrix(1, p1, p2 * p3)
  
}else if(prod(dim(penaltyfactor)) != p){
  
stop(
paste("number of elements in penaltyfactor (", length(penaltyfactor),") is not equal to the number of coefficients (", p,")", sep = "")
)
  
}else {
  
if(min(penaltyfactor) < 0){stop(paste("penaltyfactor must be positive"))}    
  
penaltyfactor <- matrix(penaltyfactor, p1, p2 * p3)
  
} 

##alpha
if(is.null(penaltyfactoralpha)){
  
  penaltyfactoralpha <- matrix(1, q, 1)
  
}else if(prod(dim(penaltyfactoralpha)) != q){
  
  stop(
    paste("number of elements in penaltyfactoralpha (", length(penaltyfactoralpha),") is not equal to the number of coefficients (", q,")", sep = "")
  )
  
}else {
  
if(min(penaltyfactoralpha) < 0){stop(paste("penaltyfactoralpha must be positive"))}    
  
} 


####check on iwls
if(is.null(iwls)){ 
  
iwls <- "exact"

}else{

if(sum(iwls == c("exact", "identity", "kron1", "kron2")) != 1){stop(paste("iwls must be correctly specified"))}

}

####check on nu
if(nu < 0 || nu > 1){stop(paste("nu must be between 0 and 1"))}

##run gdpg algorithm
res <- gdpg(X1, X2, X3,
            Z, matrix(0, n, q),##non ten design
            Y, matrix(0, n3, n1 * n2),
            matrix(0, n1, n2 * n3), #maybe V could be used in source only when S==1!!!!!
            weights, matrix(0, n3, n1 * n2),
            thetainit, matrix(0, 1, 1), matrix(0, 1, 1),
            alphainit, ## non ten par
            family,  
            penalty,
            nonten, #### non ten ind
            iwls, 
            nu,
            lambda, makelamb, nlambda, lambdaminratio,
            penaltyfactor, matrix(0, 1, 1), matrix(0, 1, 1),
            penaltyfactoralpha, ##non tenpenaltyfactor
            reltolinner, 
            reltolouter, 
            0,
            maxiter,
            steps,
            maxiterinner,
            maxiterouter,
            1,
            btinnermax,
            btoutermax,
            weightedgaussian, 0, 0, n1, n2, n3)

####checks
if(res$STOPmaxiter == 1){

warning(paste("program exit due to maximum number of iterations (",maxiter,") reached for model no ",res$endmodelno,""))

}

if(res$STOPprox == 1){

warning(paste("program exit due to maximum number of backtraking steps in the inner loop reached for model no ",res$endmodelno,""))

}

if(res$STOPnewt == 1){

warning(paste("program exit due to max number of backtraking steps in the outer loop reached for model no ",res$endmodelno,""))

}

endmodelno <- res$endmodelno + 1
Iter <- res$Iter

if(family != "gaussian" || (weightedgaussian == 1 & iwls == "kron")){
  
maxiterouterreached <- ifelse(Iter[1:endmodelno, maxiterouter] > 0, 1, 0) * (1:endmodelno)
maxiterouterreached <- maxiterouterreached[maxiterouterreached > 0]

if(length(maxiterouterreached)){

models <- paste("", maxiterouterreached)
warning(paste("maximum number of outer iterations (", maxiterouter,") reached for model(s):"), models)

}

}

maxiterinnerpossible <- sum(Iter > 0)
maxiterinnerreached <- sum(Iter >= (maxiterinner - 1)) 

if(maxiterinnerreached > 0){

warning(
paste("maximum number of inner iterations (",maxiterinner,") reached ",maxiterinnerreached," time(s) out of ",maxiterinnerpossible," possible")
)

}

out <- list()

class(out) <- "glamlasso"

out$spec <- paste("", dimglam,"-dimensional ", penalty," penalized ", family," GLAM") 
out$coef <- res$Beta[ , 1:endmodelno]
out$alpha <- res$alpha[ , 1:endmodelno]
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

return(out)

}