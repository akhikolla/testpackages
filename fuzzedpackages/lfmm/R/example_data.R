#' Genetic and phenotypic data for Arabidopsis thaliana
#'
#' A dataset containing SNP frequency and simulated phenotypic data for 170 plant accessions.
#' The variables are as follows:
#'
#' \itemize{
#'   \item genotype: binary (0 or 1) SNP frequency for 170 individuals (26943 SNPs).
#'   \item phenotype: simulated phenotypic data for 170 individuals.
#'  \item causal.set: set of indices for causal SNPs.
#'   \item chrpos: genetic map including chromosome position of each SNP.
#' }
#'
#' @docType data
#' @keywords datasets
#' @details Reference: Atwell et al (2010). 
#' Genome-wide association study of 107 phenotypes in Arabidopsis thaliana inbred lines. 
#' Nature 465, 627–631.
#' @name example.data
#' @usage data(example.data)
#' @format A list with 4 arguments: genotype, phenotype, causal.set, chrpos
NULL
