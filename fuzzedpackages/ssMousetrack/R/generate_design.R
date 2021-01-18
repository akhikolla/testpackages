#' Generate the design of a mouse-tracking experiment
#'
#' @details The function generates a dataframe containing the experimental design of a mouse-tracking study. The design is of the order (sbj,trial,variable1,...variableQ), where variable1,...,variableQ are Q categorical variables each with K_1,...,K_Q levels. The levels are codified using hundreds.
#' This is an internal function, generally not to be called by the user.
#' @param I (integer) number of individuals
#' @param J (integer) number of trials
#' @param K (list of integers) list of length Q of the number of levels for each categorical variable
#' @param Z.type (list of characters) list of length Q of the methods (symmetric or random) to generate the matrix (see \code{\link{generate_Z}})
#' @return a dataframe of the order (sbj,trial,variable1,...variableQ)
#' @export
#' @examples
#' 
#' ## Generate a design with Q = 2 categorical variables:
#' ## the first variable has K = 4 levels generated via symmetric method
#' ## the second variable has K = 3 levels generated via random method.
#' X <- generate_design(I = 10, J = 12, K = c(4,3), Z.type = c("symmetric","random"))
#' print(X)

generate_design <- function(I=10,J=12,K=c(4),Z.type=c("symmetric")){
  if(I<1 | J<1)
    stop("Positive integer should be provided for I,J")
  if(length(K)<1)
    stop("At least an integer should be provided for K")
  if(length(Z.type)<1)
    stop("At least a method should be provided to generate Z")
  

  Q <- length(K)
  Z.data <- mapply(FUN = function(q){Z <- generate_Z(I = I,J = J,K = K[q],type = Z.type[q])%*%matrix(seq(100,100*K[q],length.out = K[q]),K[q],1)},1:Q)
  X.data <- cbind(data.frame("sbj" = rep(1:I,each=J),"trial" = rep(1:J,I)),data.frame(apply(Z.data,2,as.factor)));
  names(X.data)[3:dim(X.data)[2]] <-paste("Z",seq(1,Q),sep="") #adjust for names

  return(X.data)
}
