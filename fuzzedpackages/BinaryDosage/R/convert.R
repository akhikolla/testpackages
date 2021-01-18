#' @useDynLib BinaryDosage, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom digest digest
#' @importFrom prodlim row.match
#' @importFrom utils read.table
NULL

validatebdinput <- function(bdfiles,
                            format,
                            subformat,
                            snpidformat,
                            bdoptions) {
  if (is.character(bdfiles) == FALSE)
    stop("bdfiles must be a vector of characters")

  if (is.numeric(format) == FALSE && is.integer(format) == FALSE)
    stop("format must be an integer value")
  if (length(format) != 1)
    stop("format must be an integer vector of length 1")
  if (is.numeric(format) == TRUE) {
    if (floor(format) != format)
      stop("format must be an integer")
    format <- floor(format)
  }
  if (format < 1 || format > 4)
    stop("format must be an integer value from 1 to 4")

  if (is.numeric(subformat) == FALSE && is.integer(subformat) == FALSE)
    stop("subformat must be an integer value")
  if (length(subformat) != 1)
    stop("subformat must be an integer vector of length 1")
  if (is.numeric(subformat) == TRUE) {
    if (floor(subformat) != subformat)
      stop("subformat must be an integer")
    subformat <- floor(subformat)
  }
  if (subformat < 0 || subformat > 4)
    stop("subformat must be an integer value from 0 to 4")
  if (format < 3 && subformat > 2)
    stop("subformat must be an integer value from 0 to 2 for formats 1 and 2")

  if (format == 4 & length(bdfiles) != 1)
    stop("Only one output file name is needed when using format 4")
  if (format < 4 & length(bdfiles) != 3)
    stop("Three output file names are required when using formats 1, 2, and 3")
  if (is.na(match("", bdfiles)) == FALSE)
    stop("Output file names cannot be blank")

  if (is.numeric(snpidformat) == FALSE && is.integer(snpidformat) == FALSE)
    stop("snpidformat must be an integer value")
  if (length(snpidformat) != 1)
    stop("snpidformat must be an integer vector of length 1")
  if (is.numeric(snpidformat) == TRUE) {
    if (floor(snpidformat) != snpidformat)
      stop("snpidformat must be an integer")
    snpidformat = floor(snpidformat)
  }
  if (snpidformat < -1 | snpidformat > 3)
    stop("snpidformat must be and integer from -1 to 3")

  if (is.character(bdoptions) == FALSE)
    stop("bdoptions must be a character array")
  if (length(bdoptions) > 0 & format != 4)
    stop("bdoptions can only be used with format 4")
  if (length(bdoptions) > 0) {
    if (any(is.na(match(bdoptions, c("aaf", "maf", "rsq")))) == TRUE)
      stop("Only valid bdoptions are aaf, maf, and rsq")
  }
  return(list(format = format,
              subformat = subformat))
}

###########################################################
#                  VCF to Binary Dosage                   #
###########################################################

#' Convert a VCF file to a binary dosage file
#'
#' Routine to read information from a VCF file and create
#' a binary dosage file. The function is designed to use
#' files return from the Michigan Imputation Server but will
#' run on other VCF files if they contain dosage and genetic
#' probabilities. Note: This routine can take a long time to
#' run if the VCF file is large.
#'
#' @param vcffiles A vector of file names.
#' The first is the name of the vcf file. The
#' second is name of the file that contains information
#' about the imputation of the SNPs. This file is produced
#' by minimac 3 and 4.
#' @param gz Indicator if VCF file is compressed using gzip.
#' Default value is FALSE.
#' @param bdfiles Vector of names of the output files.
#' The binary dosage file name is first. The family and
#' map files follow. For format 4, no family and map file
#' names are needed.
#' @param format The format of the output binary dosage file.
#' Allowed values are 1, 2, 3, and 4. The default value is 4.
#' Using the default value is recommended.
#' @param subformat The subformat of the format of the output
#' binary dosage file. A value of 1 or 3 indicates that only the
#' dosage value is saved. A value of 2 or 4 indicates
#' the dosage and genetic probabilities will be output. Values
#' of 3 or 4 are only allowed with formats 3 and 4. If a value
#' of zero if provided, and genetic probabilities are in the vcf
#' file, subformat 2 will be used for formats 1 and 2, and
#' subformat 4 will be used for formats 3 and 4. If the vcf file
#' does not contain genetic probabilities, subformat 1 will be
#' used for formats 1 and 2, and subformat 3 will be used for
#' formats 3 and 4. The default value is 0.
#' @param snpidformat The format that the SNP ID will be saved as.
#' -1 SNP ID not written
#' 0 - same as in the VCF file
#' 1 - chromosome:location
#' 2 - chromosome:location:reference_allele:alternate_allele
#' If snpidformat is 1 and the VCF file uses format 2, an error is
#' generated. Default value is 0.
#' @param bdoptions Character array containing any of the following
#' value, "aaf", "maf", "rsq". The presence of any of these
#' values indicates that the specified values should be
#' calculates and stored in the binary dosage file. These values only
#' apply to format 4.
#'
#' @return
#' None
#' @export
#'
#' @examples
#' # Find the vcf file names
#' vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
#' vcf1ainfo <- system.file("extdata", "set1a.info", package = "BinaryDosage")
#' bdfiles <- tempfile()
#' # Convert the file
#' vcftobd(vcffiles = c(vcf1afile, vcf1ainfo), bdfiles = bdfiles)
#' # Verify the file was written correctly
#' bdinfo <- getbdinfo(bdfiles)
vcftobd <- function(vcffiles,
                    gz = FALSE,
                    bdfiles,
                    format = 4L,
                    subformat = 0L,
                    snpidformat = 0,
                    bdoptions = character(0)) {
  if (missing(vcffiles) == TRUE)
    stop("No VCF file specified")

  if (missing(bdfiles) == TRUE)
    stop("No output files specified")

  validation <- validatebdinput(bdfiles = bdfiles,
                                format = format,
                                subformat = subformat,
                                snpidformat = snpidformat,
                                bdoptions = bdoptions)
  format <- validation$format
  subformat <- validation$subformat

  if (snpidformat == -1)
    readsnpformat = 0
  else
    readsnpformat = snpidformat
  vcfinfo <- getvcfinfo(vcffiles = vcffiles,
                        gz = gz,
                        index = FALSE,
                        snpidformat = readsnpformat)
  if (snpidformat == -1)
    vcfinfo$snpidformat = 1
  else
    vcfinfo$snpidformat = 0

  if (subformat == 0) {
    if (anyNA(vcfinfo$additionalinfo$datacolumns$genotypeprob) == TRUE)
      subformat <- 1
    else
      subformat <- 2
  }
  WriteBinaryDosageHeader(format = format,
                          subformat = subformat,
                          filename = bdfiles,
                          genefileinfo = vcfinfo,
                          bdoptions = bdoptions)
  headerinfo <- ReadBinaryDosageHeader(filename = bdfiles)
  bdwriteinfo <- AllocateBinaryDosageWriteMemory(headerinfo = headerinfo)
  vcfapply(vcfinfo = vcfinfo,
           func = WriteBinaryDosageData,
           writeinfo = bdwriteinfo)
  WriteBinaryDosageIndices(writeinfo = bdwriteinfo)
  bdinfo <- getbdinfo(bdfiles = bdfiles)
  if (is.na(match("aaf", bdoptions)) == FALSE)
    updateaaf(bdinfo)
  if (is.na(match("maf", bdoptions)) == FALSE)
    updatemaf(bdinfo)
  if (is.na(match("rsq", bdoptions)) == FALSE)
    updatersq(bdinfo)
  ##return (0)
}

###########################################################
#                  Gen to Binary Dosage                   #
###########################################################

#' Convert a gen file to a binary dosage file
#'
#' Routine to read information from a gen file and create
#' a binary dosage file. Note: This routine can take a long
#' time to run if the gen file is large.
#'
#' @param genfiles A vector of file names.
#' The first is the name of the gen file. The
#' second is name of the sample file that contains
#' the subject information.
#' @param snpcolumns Column numbers containing chromosome,
#' snpid, location, reference allele, alternate allele,
#' respectively. This must be an integer vector. All
#' values must be positive except for the chromosome.
#' The value for the chromosome may be -1 or -0.
#' -1 indicates that the chromosome value is passed to
#' the routine using the chromosome parameter.
#' 0 indicates that the chromosome value is in the snpid
#' and that the snpid has the format chromosome:other_data.
#' Default value is c(1L, 2L, 3L, 4L, 5L).
#' @param startcolumn Column number of first column with
#' genetic probabilities or dosages. Must
#' be an integer value. Default value is 6L.
#' @param impformat Number of genetic data values per
#' subject. 1 indicates dosage only, 2 indicates P(g=0)
#' and P(g=1) only, 3 indicates P(g=0), P(g=1), and
#' P(g=2). Default value is 3L.
#' @param chromosome Chromosome value to use if the
#' first value of the snpcolumns is equal to 0.
#' Default value is character().
#' @param header Indicators if the gen and sample files
#' have headers. If the gen file does not have a
#' header. A sample file must be included.
#' Default value is c(FALSE, TRUE).
#' @param gz Indicator if file is compressed using gzip.
#' Default value is FALSE.
#' @param sep Separator used in the gen file. Default
#' value is `"\t"`
#' @param bdfiles Vector of names of the output files.
#' The binary dosage file name is first. The family and
#' map files follow. For format 4, no family and map file
#' names are needed.
#' @param format The format of the output binary dosage file.
#' Allowed values are 1, 2, 3, and 4. The default value is 4.
#' Using the default value is recommended.
#' @param subformat The subformat of the format of the output
#' binary dosage file. A value of 1 or 3 indicates that only the
#' dosage value is saved. A value of 2 or 4 indicates
#' the dosage and genetic probabilities will be output. Values
#' of 3 or 4 are only allowed with formats 3 and 4. If a value
#' of zero if provided, and genetic probabilities are in the vcf
#' file, subformat 2 will be used for formats 1 and 2, and
#' subformat 4 will be used for formats 3 and 4. If the vcf file
#' does not contain genetic probabilities, subformat 1 will be
#' used for formats 1 and 2, and subformat 3 will be used for
#' formats 3 and 4. The default value is 0.
#' @param snpidformat The format that the SNP ID will be saved as.
#' -1 - SNP ID not written.
#' 0 - same as in the VCF file.
#' 1 - chromosome:location.
#' 2 - chromosome:location:reference_allele:alternate_allele.
#' If snpidformat is 1 and the VCF file uses format 2, an error is
#' generated. Default value is 0.
#' @param bdoptions Character array containing any of the following
#' value, "aaf", "maf", "rsq". The presence of any of these
#' values indicates that the specified values should be
#' calculates and stored in the binary dosage file. These values only
#' apply to format 4.
#'
#' @return
#' None
#' @export
#'
#' @examples
#' # Find the gen file names
#' gen3afile <- system.file("extdata", "set3a.imp", package = "BinaryDosage")
#' gen3asample <- system.file("extdata", "set3a.sample", package = "BinaryDosage")
#' # Get temporary output file name
#' bdfiles <- tempfile()
#' # Convert the file
#' gentobd(genfiles = c(gen3afile, gen3asample),
#'         snpcolumns = c(0L, 2L:5L),
#'         bdfiles = bdfiles)
#' # Verify the file was written correctly
#' bdinfo <- getbdinfo(bdfiles = bdfiles)
gentobd <- function(genfiles,
                    snpcolumns = 1L:5L,
                    startcolumn = 6L,
                    impformat = 3L,
                    chromosome = character(),
                    header = c(FALSE, TRUE),
                    gz = FALSE,
                    sep = "\t",
                    bdfiles,
                    format = 4L,
                    subformat = 0L,
                    snpidformat = 0L,
                    bdoptions = character(0)) {
  if (missing(genfiles) == TRUE)
    stop("No gen file specified")

  if (missing(bdfiles) == TRUE)
    stop("No output files specified")

  validation <- validatebdinput(bdfiles = bdfiles,
                                format = format,
                                subformat = subformat,
                                snpidformat = snpidformat,
                                bdoptions = bdoptions)
  format <- validation$format
  subformat <- validation$subformat

  if (snpidformat == -1)
    readsnpformat = 0L
  else
    readsnpformat = as.integer(snpidformat)
  geninfo <- getgeninfo(genfiles = genfiles,
                        snpcolumns = snpcolumns,
                        startcolumn = startcolumn,
                        impformat = impformat,
                        chromosome = chromosome,
                        header = header,
                        gz = gz,
                        index = FALSE,
                        snpidformat = readsnpformat,
                        sep = sep)
  if (snpidformat == -1)
    geninfo$snpidformat = 1
  else
    geninfo$snpidformat = 0

  if (subformat == 0) {
    if (geninfo$additionalinfo$format == 1L)
      subformat <- 1L
    else
      subformat <- 2L
  }
  WriteBinaryDosageHeader(format = format,
                          subformat = subformat,
                          filename = bdfiles,
                          genefileinfo = geninfo,
                          bdoptions = bdoptions)
  headerinfo <- ReadBinaryDosageHeader(filename = bdfiles)
  bdwriteinfo <- AllocateBinaryDosageWriteMemory(headerinfo = headerinfo)
  genapply(geninfo = geninfo,
           func = WriteBinaryDosageData,
           writeinfo = bdwriteinfo)
  WriteBinaryDosageIndices(writeinfo = bdwriteinfo)
  bdinfo <- getbdinfo(bdfiles = bdfiles)
  if (is.na(match("aaf", bdoptions)) == FALSE)
    updateaaf(bdinfo)
  if (is.na(match("maf", bdoptions)) == FALSE)
    updatemaf(bdinfo)
  if (is.na(match("rsq", bdoptions)) == FALSE)
    updatersq(bdinfo)
}

