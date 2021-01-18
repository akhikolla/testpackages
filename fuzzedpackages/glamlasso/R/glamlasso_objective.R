#
#     Description of this R script:
#     R interface for the function objective in the glamlasso package.
#
#     Intended for use with R.
#     Copyright (C) 2015 Adam Lund
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
 
#' @name objective
#' 
#' @aliases glamlasso_objective 
#' 
#' @title Compute objective values 
#' 
#' @description Computes the objective values of the penalized log-likelihood problem
#'              for the models implemented in the package glamlasso.
#'              
#' @usage objective(Y, 
#'           Weights, 
#'           X, 
#'           Beta, 
#'           lambda,
#'           penalty.factor, 
#'           family,
#'           penalty)
#' 
#' @param Y The response values, an array of size \eqn{n_1 \times \cdots \times n_d}. 
#' @param Weights Observation weights, an array of size \eqn{n_1 \times \cdots \times n_d}.
#' @param X A list containing the tensor components of the tensor design matrix, each of size
#'  \eqn{n_i \times p_i}.
#' @param Beta A coefficient matrix of size \eqn{p_1\cdots p_d \times }\code{nlambda}.
#' @param lambda The sequence of penalty parameters for the regularization path.
#' @param penalty.factor An array of size  \eqn{p_1 \times \cdots \times p_d}. Is multiplied with each 
#' element in  \code{lambda} to allow differential shrinkage on the coefficients. 
#' @param family A string specifying the model family (essentially the response distribution).
#' @param penalty A string specifying the penalty.
#' 
#' @return
#' A vector of length \code{length(lambda)} containing the objective values for each \code{lambda} value. 
#' 
#' @examples
#' \dontrun{
#' n1 <- 65; n2 <- 26; n3 <- 13; p1 <- 13; p2 <- 5; p3 <- 4
#' X1 <- matrix(rnorm(n1 * p1), n1, p1) 
#' X2 <- matrix(rnorm(n2 * p2), n2, p2) 
#' X3 <- matrix(rnorm(n3 * p3), n3, p3) 
#' Beta <- array(rnorm(p1 * p2 * p3) * rbinom(p1 * p2 * p3, 1, 0.1), c(p1 , p2, p3))
#' mu <- RH(X3, RH(X2, RH(X1, Beta)))
#' Y <- array(rnorm(n1 * n2 * n3, mu), dim = c(n1, n2, n3))
#' fit <- glamlasso(list(X1, X2, X3), Y, family = "gaussian", penalty = "lasso", iwls = "exact")
#' objfit <- objective(Y, NULL, list(X1, X2, X3), fit$coef, fit$lambda, NULL, fit$family)
#' plot(objfit, type = "l")
#' }

objective <-function(Y, Weights, X, Beta, lambda, penalty.factor, family, penalty) {
  
##get dimensions of problem
  
  dimglam <- length(X)
  
  if (dimglam < 2 || dimglam > 3){
    
    stop(paste("the dimension of the GLAM must be 2 or 3!"))
    
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

##check if lambda is specified
if(is.null(lambda)){stop(paste("no lambda sequence is specified"))}

##check if penalty.factor is specified; if not set to all ones
if(is.null(penalty.factor)){
  
penalty.factor <- matrix(1, p1, p2 * p3)

}else if(prod(dim(penalty.factor)) != p){
  
stop(
paste("number of elements in penalty.factor (", length(penalty.factor),") is not equal to the number of coefficients (", p,")", sep = "")
)

}else {penalty.factor <- matrix(penalty.factor, p1, p2 * p3)}  

##check if weights/nof trials is specfied for the binomial model
if(family == "binomial" & is.null(Weights)){

stop(paste("for binomial model number of trials (Weights) must be specified"))

}

if(is.null(Weights)){
  
Weights <- matrix(1, n1, n2 * n3)

}else{Weights <- matrix(Weights, n1, n2 * n3)}

Y <- matrix(Y, n1, n2 * n3)
BetaArr <- array(Beta, c(p1, p2 * p3, length(lambda)))

##get objective values
out <- getobj(Y, Weights,  X1, X2, X3, BetaArr, lambda, penalty.factor, family, penalty)$Obj

return(out)

}