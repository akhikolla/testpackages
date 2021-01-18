
#' Individual level dataset
#'
#' A simulated individual level dataset for PMR.
#' 
#' @format   A list contains the following objects:
#' \describe{
#'   \item{zx}{the standardized genotype matrix for 465 individuals and 50 cis-SNPs in eQTL data.}
#'   \item{zy}{the standardized genotype matrix for 2000 individuals and 50 cis-SNPs in GWAS data.}
#'   \item{x}{the standarized gene expression vector for 465 individuals in eQL data.}
#'   \item{y}{the standarized complex trait vector for 2000 individuals in GWAS data.}
#' }
#' 
#' @export
"Exampleindividual"

#' Summary level dataset 
#'
#' A simulated summary level dataset for PMR
#' 
#' @format   A list contains the following objects:
#' \describe{
#'   \item{betax}{the cis-SNP effect size vector for one specific gene in eQTL data.}
#'   \item{betay}{the cis-SNP effect size vector for one specific gene in GWAS data.}
#'   \item{Sigma1}{the LD matrix in eQTL data.}
#'   \item{Sigma2}{the LD matrix in GWAS data.}
#'   \item{n1}{the sample size of eQTL data.}
#'   \item{n2}{the sample size of GWAS data.}
#' }
#' 
#' @export
"Examplesummary"