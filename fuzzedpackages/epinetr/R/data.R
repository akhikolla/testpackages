#' Example genotype matrix.
#'
#' An example genotype matrix for 100 SNPs across 500 individuals.
#'
#' @format \code{geno100snp} is a matrix with 500 rows and 200
#' columns, used in examples when constructing a population from
#' genotypes using 500 individuals and 100 SNPs.
#'
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#'
#' @seealso \code{\link{Population}}, \code{\link{map100snp}}
"geno100snp"

#' Example map.
#'
#' An example map \code{data.frame} for 100 SNPs across 22
#' chromosomes.
#'
#' @format \code{map100snp} is a \code{data.frame} with 100 rows and
#' 3 variables:
#' \describe{
#'   \item{V1}{SNP ID}
#'   \item{V2}{chromosome ID}
#'   \item{v3}{position on chromosome, in base pairs}
#' }
#'
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#'
#' @seealso \code{\link{Population}}, \code{\link{geno100snp}}
"map100snp"

#' Example incidence matrix.
#'
#' An example incidence matrix for 19 pairwise interactions across
#' 20 QTLs.
#'
#' @format \code{rincmat100snp} is a matrix with 20 rows and 19
#' columns, used in examples when constructing an epistatic network
#' using a user-supplied incidence matrix.
#'
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#'
#' @seealso \code{\link{attachEpiNet}}
"rincmat100snp"
