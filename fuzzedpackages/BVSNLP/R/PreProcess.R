#' Preprocessing the design matrix, preparing it for variable selection
#' procedure
#'
#' @description This function preprocesses the design matrix by removing
#' columns that contain \code{NA}'s or are all zero. It also standardizes
#' non-binary columns to have mean zero and variance one. The user has the
#' choice of log transforming continuous covariates before scaling them.
#' @param X The \code{n} times \code{p} design matrix. The columns should
#' represent genes and rows represent the observations. The column names are
#' used as gene names so they should not be left as \code{NULL}. Note that the
#' input matrix \code{X} should NOT contain vector of \code{1}'s representing
#' the intercept.
#' @param logT A boolean variable determining if log transform should be done
#' on continuous columns before scaling them. Note that those columns should
#' not contain any zeros or negative values.
#' @author Amir Nikooienejad
#' @return It returns a list having the following objects:
#' \item{X}{The filtered design matrix which can be used in variable selection
#' procedure. Binary columns are moved to the end of the design matrix.}
#' \item{gnames}{Gene names read from the column names of the filtered design
#' matrix.}
#' @examples
#' ### Constructing a synthetic design matrix for the purpose of preprocessing
#' ### imposing columns with different scales
#' n <- 40
#' p1 <- 50
#' p2 <- 150
#' p <- p1 + p2
#' X1 <- matrix(rnorm(n*p1, 1, 2), ncol = p1)
#' X2 <- matrix(rnorm(n*p2), ncol = p2)
#' X <- cbind(X1, X2)
#'
#' ### putting NA elements in the matrix
#' X[3,85] <- NA
#' X[25,85] <- NA
#' X[35,43] <- NA
#' X[15,128] <- NA
#' colnames(X) <- paste("gene_",c(1:p),sep="")
#'
#' ### Running the function. Note the intercept column that is added as the
#' ### first column in the "logistic" family
#' Xout <- PreProcess(X)
#' dim(Xout$X)[2] == (p + 1) ## 1 is added because intercept column is included
#' ## This is FALSE because of the removal of columns with NA elements
PreProcess <- function(X, logT = FALSE){
  Xin <- X

  ex0 <- which(apply(Xin,2,function(x) all(x==0))) 
  if (length(ex0)) Xin <- Xin[,-ex0]
  
  ex1 <- which(colSums(is.na(Xin)) != 0)
  if (length(ex1)) Xin <- Xin[, -ex1]

  XX2 <- Xin
  bincols <- which(apply(Xin, 2, function(x) {all(x %in% 0:1)}))
  if (length(bincols)){
    if(logT) XX2[,-bincols] <- log(XX2[,-bincols])
    XX2[,-bincols] <- scale(XX2[,-bincols])
  } else {
    if(logT) Xin <- log(Xin)
    XX2 <- scale(Xin)
  }

  ex2 <- which(colSums(is.na(XX2)) != 0)
  if (length(ex2)) XX2 <- XX2[, -ex2]

  gene_names <- colnames(XX2)
  return(list(X = XX2, gnames = gene_names))
}
