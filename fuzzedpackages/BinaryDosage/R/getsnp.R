#' Read SNP data from a binary dosage file
#'
#' Routine to read the dosage and genetic probabilities about
#' a SNP from a binary dosage file
#'
#' @param bdinfo Information about a binary dosage file return
#' from getbdinfo
#' @param snp The SNP to read the information about. This may
#' be the SNP ID or the index of the SNP in the snps dataset in
#' the bdinfo list
#' @param dosageonly Indicator to return the dosages only or the
#' dosages allowing with the genetic probabilities.
#' Default value is TRUE
#'
#' @return
#' A list with either the dosages or the dosages and the genetic
#' probabilities.
#' @export
#'
#' @examples
#' # Get the information about the file
#' vcf1abdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
#' bdinfo <- getbdinfo(bdfiles = vcf1abdfile)
#'
#' # Read the first SNP
#' getsnp(bdinfo, 1, FALSE)
getsnp <- function(bdinfo, snp, dosageonly = TRUE) {
  if (missing(bdinfo) == TRUE)
    stop("bdinfo missing")
  if (is.na(match("genetic-info", class(bdinfo))) == TRUE)
      stop("bdinfo is not of class genetic-info")
  if (is.na(match("bdose-info", class(bdinfo$additionalinfo))) == TRUE)
    stop("bdinfo does not contain information about a binary dosage file")

  if (missing(snp) == TRUE)
    stop("No SNP specified")
  if (length(snp) != 1)
    stop("snp must be of length one")
  if (is.character(snp) == TRUE) {
    x <- match(snp, bdinfo$snps$snpid)
    if (is.na(x))
      stop("Cannot find SNP in bdinfo")
    snp <- x
  } else if (is.numeric(snp) == TRUE) {
    if (floor(snp) != snp)
      stop("snp must be a character or integer value")
    snp <- as.integer(floor(snp))
  }
  if (is.integer(snp) == FALSE)
    stop("snp must be a character or integer value")
  if (snp < 1 | snp > nrow(bdinfo$snps))
    stop("snp value out or range")

  dosage <- numeric(nrow(bdinfo$samples))
  p0 <- numeric(nrow(bdinfo$samples))
  p1 <- numeric(nrow(bdinfo$samples))
  p2 <- numeric(nrow(bdinfo$samples))
  us <- numeric(2*nrow(bdinfo$samples))
  dosage[1:nrow(bdinfo$samples)] <- NA
  p0[1:nrow(bdinfo$samples)] <- NA
  p1[1:nrow(bdinfo$samples)] <- NA
  p2[1:nrow(bdinfo$samples)] <- NA
  ReadBinaryDosageData(bdinfo, snp, dosage, p0, p1, p2, us)
  if (dosageonly == TRUE)
    return(list(dosage = dosage))
  return(list(dosage = dosage,
              p0 = p0,
              p1 = p1,
              p2 = p2))
}
