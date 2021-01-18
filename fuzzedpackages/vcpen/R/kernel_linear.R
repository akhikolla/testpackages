#' Variance Component Linear Kernel Matrix 
#'
#' Variance component Linear kernel matrix from genotype dosage 
#'
#' @param dose data.frame or matrix with 
#' @param method type of kernel; currently only linear kernel implemented
#' @return square symmetric kernel matrix for subject similarity by genotype dosage
#' @examples
#' data(vcexample)
#' Kern1 <- kernel_linear(dose[,which(doseinfo[,1]==1)], method="linear")
#' Kern1[1:5,1:5]
#' 
#' @author JP Sinnwell, DJ Schaid
#' @seealso \code{\link{vcpen}}
#' @name kernel_linear
NULL
#> NULL
#' @rdname kernel_linear
#' @export

kernel_linear <- function(dose, method="linear"){
  ## linear kernel matrix based matrix of SNP doses of minor alleles
  n.snp <- ncol(dose)
  ## kernel matrix based on centered dose (consistent with how genetic relationship
  ## matrix is calcualted in Gemma)
  dose.mean <- apply(dose, 2, mean, na.rm=TRUE)
  dose.sd  <- sqrt(apply(dose, 2, var, na.rm=TRUE))
  zdose <- t((t(dose) - dose.mean)/dose.sd)
  kmat <- zdose %*% t(zdose) / n.snp
  return(kmat)
}




