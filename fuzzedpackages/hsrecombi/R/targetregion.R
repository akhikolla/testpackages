#' @title Description of the targetregion data set
#' @name targetregion
#' @docType data
#' @description The data set contains sire haplotypes, assignment of progeny to
#'   sire, progeny genotypes and physical map information in a target region
#' \describe{
#'   The raw data can be downloaded at the source given below. Then,
#'   executing the following R code leads to the data provided in
#'   \code{targetregion.RData}.
#'   \item{\code{hapSire}}{matrix of sire haplotypes of each sire; 2 lines per
#'    sire; 1. column contains sireID}
#'   \item{\code{daughterSire}}{vector of sire ID for each progeny}
#'   \item{\code{genotype.chr}}{matrix of progeny genotypes}
#'   \item{\code{map.chr}}{SNP marker map in target region}
#' }
#' @source The data are available from the RADAR repository
#'   \url{https://dx.doi.org/10.22000/280}
#' @examples
#' \donttest{
#' ## list of haplotypes of sires for each chromosome
#' load('sire_haplotypes.RData')
#' ## assign progeny to sire
#' daughterSire <- read.table('assign_to_family.txt')[, 1]
#' ## progeny genotypes
#' X <- as.matrix(read.table('XFam-ARS.txt'))
#' ## physical and approximated genetic map
#' map <- read.table('map50K_ARS_reordered.txt', header = T)
#' ## select target region
#' chr <- 1
#' window <- 301:600
#' ## map information of target region
#' map.chr <- map[map$Chr == chr, ][window, 1:5]
#' ## matrix of sire haplotypes in target region
#' hapSire <- rlist::list.rbind(haps[[chr]])
#' sireID <- 1:length(unique(daughterSire))
#' hapSire <- cbind(rep(sireID, each = 2), hapSire[, window])
#' ## matrix of progeny genotypes
#' genotype.chr <- X[, map.chr$SNP]
#' }
#' @importFrom rlist list.rbind
NULL
