###########################################################
#                Binary Dosage                            #
###########################################################

#' Apply a function to each SNP in a
#' binary dosage file
#'
#' A routine that reads in the SNP data serially
#' from a binary dosage file and applies a user
#' specified function to the data.
#'
#' @param bdinfo List with information about the
#' binary dosage file returned from getbdinfo
#' @param func A user supplied function to apply
#' to the data for each snp. The function must be
#' provide with the following parameters, dosage,
#' p0, p1, and p2, where dosage is the dosage values
#' for each subject and p0, p1, and p2 are the
#' probabilities that a subject has zero, one,
#' and two copies of the alternate allele,
#' respectively.
#' @param ... Additional parameters needed by the
#' user supplied function
#'
#' @return A list with length equal to the number
#' of SNPs in the binary dosage file. Each element
#' of the list is the value returned by the user
#' supplied function
#' @export
#' @family Iterating functions
#' @examples
#' # Get information about a binary dosage file
#'
#' vcf1abdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
#' bdinfo <- getbdinfo(bdfiles = vcf1abdfile)
#'
#' # Apply the getaaf, get alternate allele frequency, function
#' # to all the SNPs in the binary dosage file
#'
#' aaf <- bdapply(bdinfo = bdinfo,
#'                func = BinaryDosage:::getaaf)
bdapply <- function(bdinfo, func, ...) {
  if (missing(bdinfo) == TRUE)
    stop("No binary dosage file information specified")
  if (is.na(match("genetic-info", class(bdinfo))) == TRUE)
    stop("bdinfo does not contain information about a binary dosage file")
  if (is.na(match("bdose-info", class(bdinfo$additionalinfo))) == TRUE)
    stop("bdinfo does not contain information about a binary dosage file")

  if (missing(func) == TRUE)
    stop("No function specified")
  if (is.na(match("function", class(func))) == TRUE)
    stop("func is not a function")

  retval <- vector("list", nrow(bdinfo$snps))
  dosage <- numeric(nrow(bdinfo$samples))
  p0 <- numeric(nrow(bdinfo$samples))
  p1 <- numeric(nrow(bdinfo$samples))
  p2 <- numeric(nrow(bdinfo$samples))
  us <- integer(2 * nrow(bdinfo$samples))
  for (i in 1:nrow(bdinfo$snps)) {
    dosage[1:nrow(bdinfo$samples)] <- NA
    p0[1:nrow(bdinfo$samples)] <- NA
    p1[1:nrow(bdinfo$samples)] <- NA
    p2[1:nrow(bdinfo$samples)] <- NA
    ReadBinaryDosageData(bdinfo, i, dosage, p0, p1, p2, us)
    retval[[i]] <- func(dosage, p0, p1, p2, ...)
  }
  return (retval)
}

###########################################################
#                        VCF                              #
###########################################################

#' Apply a function to each SNP in a vcf file
#'
#' A routine that reads in the SNP data serially
#' from a vcf file and applies a user specified
#' function to the data.
#'
#' @param vcfinfo List with information about the
#' vcf file returned from getvcfinfo
#' @param func A user supplied function to apply
#' to the data for each snp. The function must be
#' provide with the following parameters, dosage,
#' p0, p1, and p2, where dosage is the dosage values
#' for each subject and p0, p1, and p2 are the
#' probabilities that a subject has zero, one,
#' and two copies of the alternate allele,
#' respectively.
#' @param ... Additional parameters needed by the
#' user supplied function
#'
#' @return A list with length equal to the number
#' of SNPs in the vcf file. Each element
#' of the list is the value returned by the user
#' supplied function
#' @export
#' @family Iterating functions
#' @examples
#' # Get information about a vcf file
#'
#' vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
#' vcfinfo <- getvcfinfo(vcffiles = vcf1afile)
#'
#' # Apply the getaaf, get alternate allele frequency, function
#' # to all the SNPs in the vcf file
#'
#' aaf <- vcfapply(vcfinfo = vcfinfo,
#'                 func = BinaryDosage:::getaaf)
vcfapply <- function(vcfinfo, func, ...) {
  if (missing(vcfinfo) == TRUE)
    stop("No vcf file information specified")
  if (is.na(match("genetic-info", class(vcfinfo))) == TRUE)
    stop("vcfinfo does not appear to contain information about a vcf file")
  if (is.na(match("vcf-info", class(vcfinfo$additionalinfo))) == TRUE)
    stop("vcfinfo does not appear to contain information about a vcf file")

  if (missing(func) == TRUE)
    stop("No function specified")
  if (is.na(match("function", class(func))) == TRUE)
    stop("func is not a function")

  retval <- vector("list", nrow(vcfinfo$snps))
  if (vcfinfo$additionalinfo$gzipped == FALSE)
    con <- file(vcfinfo$filename, "r")
  else
    con <- gzfile(vcfinfo$filename, "r")
  line <- readLines(con, n = vcfinfo$additionalinfo$headerlines)

  dosage <- numeric(nrow(vcfinfo$samples))
  p0 <- numeric(nrow(vcfinfo$samples))
  p1 <- numeric(nrow(vcfinfo$samples))
  p2 <- numeric(nrow(vcfinfo$samples))

  for (i in 1:nrow(vcfinfo$snps)) {
    line <- readLines(con, n = 1)
    x <- unlist(strsplit(line, "\t"))
    y <- unlist(strsplit(x[10:length(x)], ":"))
    if (length(vcfinfo$additionalinfo$datacolumns$dosage) == 1) {
      dosagecol <- vcfinfo$additionalinfo$datacolumns$dosage
      gpcol <- vcfinfo$additionalinfo$datacolumns$genotypeprob
      numcolumns <- vcfinfo$additionalinfo$datacolumns$numcolumns
    } else {
      dosagecol <- vcfinfo$additionalinfo$datacolumns$dosage[i]
      gpcol <- vcfinfo$additionalinfo$datacolumns$genotypeprob[i]
      numcolumns <- vcfinfo$additionalinfo$datacolumns$numcolumns[i]
    }
    if(is.na(dosagecol) == FALSE) {
      dosage[1:nrow(vcfinfo$samples)] <- as.numeric(y[seq(dosagecol, length(y) - numcolumns + dosagecol, numcolumns)])
    }
    if(is.na(gpcol) == FALSE) {
      gpstring <- y[seq(gpcol, length(y) - numcolumns + gpcol, numcolumns)]
      z <- unlist(strsplit(gpstring, ","))
      p0[1:nrow(vcfinfo$samples)] <- as.numeric(z[seq(1, length(z) - 2, 3)])
      p1[1:nrow(vcfinfo$samples)] <- as.numeric(z[seq(2, length(z) - 1, 3)])
      p2[1:nrow(vcfinfo$samples)] <- as.numeric(z[seq(3, length(z), 3)])
    } else {
      p0[1:nrow(vcfinfo$samples)] <- NA
      p1[1:nrow(vcfinfo$samples)] <- NA
      p2[1:nrow(vcfinfo$samples)] <- NA
    }
    if(is.na(dosagecol) == FALSE) {
      dosage[1:nrow(vcfinfo$samples)] <- as.numeric(y[seq(dosagecol, length(y) - numcolumns + dosagecol, numcolumns)])
    } else {
      dosage[1:nrow(vcfinfo$samples)] <- p1 + p2 + p2
    }
    retval[[i]] <- func(dosage = dosage,
                        p0 = p0,
                        p1 = p1,
                        p2 = p2
                        , ...)
  }
  close(con)
  return (retval)
}

###########################################################
#                 GEN (Impute2)                           #
###########################################################

#' Apply a function to each SNP in a
#' gen, impute2, file
#'
#' A routine that reads in the SNP data serially
#' from a gen file and applies a user
#' specified function to the data.
#'
#' @param geninfo List with information about the
#' gen, impute2, file returned from [getgeninfo]
#' @param func A user supplied function to apply
#' to the data for each snp. The function must be
#' provide with the following parameters, dosage,
#' p0, p1, and p2, where dosage is the dosage values
#' for each subject and p0, p1, and p2 are the
#' probabilities that a subject has zero, one,
#' and two copies of the alternate allele,
#' respectively.
#' @param ... Additional parameters needed by the
#' user supplied function
#'
#' @return A list with length equal to the number
#' of SNPs in the vcf file. Each element
#' of the list is the value returned by the user
#' supplied function
#' @export
#'
#' @family Iterating functions
#' @examples
#' # Get information about a gen, impute2, file
#'
#' gen1afile <- system.file("extdata", "set1a.imp", package = "BinaryDosage")
#' geninfo <- getgeninfo(genfiles = gen1afile,
#'                       snpcolumns = c(1L, 3L, 2L, 4L, 5L),
#'                       header = TRUE)
#'
# Apply the getaaf, get alternate allele frequency, function
# to all the SNPs in the vcf file
#'
#' aaf <- genapply(geninfo = geninfo,
#'                 func = BinaryDosage:::getaaf)
genapply <- function(geninfo, func, ...) {
  if (missing(geninfo) == TRUE)
    stop("No gen file information specified")
  if (is.na(match("genetic-info", class(geninfo))) == TRUE)
    stop("geninfo does not appear to contain information about a gen file")
  if (is.na(match("gen-info", class(geninfo$additionalinfo))) == TRUE)
    stop("geninfo does not appear to contain information about a gen file")

  if (missing(func) == TRUE)
    stop("No function specified")
  if (is.na(match("function", class(func))) == TRUE)
    stop("func is not a function")

  retval <- vector("list", nrow(geninfo$snps))
  if (geninfo$additionalinfo$gzipped == FALSE)
    con <- file(geninfo$filename, "r")
  else
    con <- gzfile(geninfo$filename, "r")
  line <- readLines(con, n = geninfo$additionalinfo$headersize)

  dosage <- numeric(nrow(geninfo$samples))
  p0 <- numeric(nrow(geninfo$samples))
  p0[1:length(p0)] <- NA
  p1 <- p0
  p2 <- p0
  for (i in 1:nrow(geninfo$snps)) {
    line <- readLines(con, n = 1)
    x <- unlist(strsplit(line, geninfo$additionalinfo$sep))
    y <- x[geninfo$additionalinfo$startcolumn:length(x)]
    if (geninfo$additionalinfo$format == 1) {
      dosage[1:length(dosage)] <- as.numeric(y)
    } else if (geninfo$additionalinfo$format == 2) {
      p0[1:length(p0)] <- as.numeric(y[seq(1, length(y) - 1, 2)])
      p1[1:length(p1)] <- as.numeric(y[seq(2, length(y), 2)])
      p2[1:length(p2)] <- ifelse(p0 + p1 < 1, 1 - p0 - p1, 0.)
      dosage[1:length(dosage)] <- ifelse(p1 + p2 + p2 < 2, p1 + p2 + p2, 2.)
    } else {
      p0[1:length(p0)] <- as.numeric(y[seq(1, length(y) - 2, 3)])
      p1[1:length(p1)] <- as.numeric(y[seq(2, length(y) - 1, 3)])
      p2[1:length(p2)] <- as.numeric(y[seq(3, length(y), 3)])
      dosage[1:length(dosage)] <- ifelse(p1 + p2 + p2 < 2, p1 + p2 + p2, 2.)
    }

    retval[[i]] <- func(dosage = dosage,
                        p0 = p0,
                        p1 = p1,
                        p2 = p2
                        , ...)
  }
  close(con)
  return (retval)
}
