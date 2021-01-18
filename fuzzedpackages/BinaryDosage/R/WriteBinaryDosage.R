#***************************************************************************#
#                                                                           #
#                       Writing Binary Dosage files                         #
#                                                                           #
#***************************************************************************#


#***************************************************************************#
#             Support functions for subject and SNP data                    #
#***************************************************************************#

# Formats 1, 2, and 3 all have a separate family and map file
# This routine saves the data frames in RDS format
WriteFamilyAndMapFiles <- function(filename, samples, snps) {
  saveRDS(samples, filename[2])
  saveRDS(snps, filename[3])
  return (md5 <- c(digest(samples, "md5"), digest(snps, "md5")))
}

# Find the groups value in the genetic file info. If it doesn't
# exist it returns the number of samples
FindGroups <- function(geneticfileinfo) {
  x <- match("groups", names(geneticfileinfo$additionalinfo))
  if (is.na(x) == TRUE)
    return (nrow(geneticfileinfo$samples))
  return(geneticfileinfo$additional$groups)
}

# Create the subject and family strings to write to the binary
# dosage header for format 4
SIDandFID4 <- function(genefileinfo) {
  sid <- paste(genefileinfo$samples$sid, collapse = '\t')
  if (genefileinfo$usesfid == TRUE)
    fid <- paste(genefileinfo$samples$fid, collapse = '\t')
  else
    fid <- ""
  return (list(sid = sid, fid = fid))
}

# Find a numeric vector in the SNP info data frame or the
# bdoptions. If option is specified, a vector of zeros is
# return. If it is not found in the options, it
# returns the vector if found in the SNP info, otherwise
# it returns a numeric vector of length 0
FindSNPInfoNumeric <- function(tofind, snpinfo, numsnps, bdoptions) {
  if (is.na(match(tofind, bdoptions)) == FALSE)
    return (numeric(numsnps))
  infocol <- match(tofind, names(snpinfo))
  if (is.na(infocol) == FALSE)
    return (as.vector(snpinfo[[infocol]]))
  return (numeric(0))
}

# Find the SNP information needed for format 4
FindBDSNPInfo <- function(genefileinfo, bdoptions) {
  snpcolnames <- colnames(genefileinfo$snps)
  if (genefileinfo$snpidformat == 0)
    snpid <- genefileinfo$snps$snpid
  else
    snpid <- ""
  chr <- genefileinfo$snps$chromosome
  if (genefileinfo$onechr == TRUE)
    chr <- chr[1]
  loc <- genefileinfo$snps$location
  ref <- genefileinfo$snps$reference
  alt <- genefileinfo$snps$alternate

  aaf <- FindSNPInfoNumeric("aaf", genefileinfo$snpinfo, nrow(genefileinfo$snps), bdoptions)
  maf <- FindSNPInfoNumeric("maf", genefileinfo$snpinfo, nrow(genefileinfo$snps), bdoptions)
  # Average call cannot be calculated, therefore bdoptions are meaningless
  avgcall <- FindSNPInfoNumeric("avgcall", genefileinfo$snpinfo, nrow(genefileinfo$snps), "")
  rsq <- FindSNPInfoNumeric("rsq", genefileinfo$snpinfo, nrow(genefileinfo$snps), bdoptions)

  snpid <- paste0(snpid, collapse = '\t')
  chr <- paste0(chr, collapse = '\t')
  ref <- paste0(ref, collapse = '\t')
  alt <- paste0(alt, collapse = '\t')
  return (list(snpid = snpid,
               chromosome = chr,
               location = loc,
               reference = ref,
               alternate = alt,
               aaf = aaf,
               maf = maf,
               avgcall = avgcall,
               rsq = rsq))
}

#***************************************************************************#
#                                                                           #
#                     Writing the Binary Dosage header                      #
#                                                                           #
#***************************************************************************#


# Writes the header for the various formats of the formats
# of the binary dosage file. These vary for all the different
# formats.
WriteBinaryDosageHeader <- function(format, subformat, filename, genefileinfo, bdoptions) {
  writeHeaderFunc <- list(f1 <- c(WriteBinaryDosageHeader1, WriteBinaryDosageHeader1),
                          f2 <- c(WriteBinaryDosageHeader1, WriteBinaryDosageHeader1),
                          f3 <- c(WriteBinaryDosageHeader31, WriteBinaryDosageHeader31, WriteBinaryDosageHeader33, WriteBinaryDosageHeader34),
                          f4 <- c(WriteBinaryDosageHeader41, WriteBinaryDosageHeader41, WriteBinaryDosageHeader43, WriteBinaryDosageHeader44))
  WriteBinaryDosageBaseHeader(filename[1], format - 1, subformat - 1)
  return (writeHeaderFunc[[format]][[subformat]](filename, genefileinfo, bdoptions))
}

WriteBinaryDosageHeader1 <- function(filename, genefileinfo, bdoptions) {
  md5 <- WriteFamilyAndMapFiles(filename, genefileinfo$samples, genefileinfo$snps)
  return (0)
}

WriteBinaryDosageHeader31 <- function(filename, genefileinfo, bdoptions) {
  md5 <- WriteFamilyAndMapFiles(filename, genefileinfo$samples, genefileinfo$snps)
  WriteBinaryDosageHeader3A(filename[1], nrow(genefileinfo$samples))
  return (0)
}

WriteBinaryDosageHeader33 <- function(filename, genefileinfo, bdoptions) {
  md5 <- WriteFamilyAndMapFiles(filename, genefileinfo$samples, genefileinfo$snps)
  WriteBinaryDosageHeader3B(filename[1], md5[1], md5[2], 0)
  return (0)
}

WriteBinaryDosageHeader34 <- function(filename, genefileinfo, bdoptions) {
  md5 <- WriteFamilyAndMapFiles(filename, genefileinfo$samples, genefileinfo$snps)
  WriteBinaryDosageHeader3B(filename[1], md5[1], md5[2], nrow(genefileinfo$snps))
  return (0)
}

WriteBinaryDosageHeader4 <- function(filename, genefileinfo, bdoptions,
                                     headerEntries, offsets, numindices) {
  subInfo <- SIDandFID4(genefileinfo)
  snpInfo <- FindBDSNPInfo(genefileinfo, bdoptions)

  WriteBinaryDosageHeader4A(filename[1],
                            headerEntries,
                            nrow(genefileinfo$samples),
                            nrow(genefileinfo$snps),
                            FindGroups(genefileinfo),
                            subInfo$sid[1],
                            subInfo$fid[1],
                            snpInfo$snpid[1],
                            snpInfo$chromosome[1],
                            snpInfo$location,
                            snpInfo$reference[1],
                            snpInfo$alternate[1],
                            snpInfo$aaf,
                            snpInfo$maf,
                            snpInfo$avgcall,
                            snpInfo$rsq,
                            offsets,
                            numindices)
  return (0)
}

WriteBinaryDosageHeader41 <- function(filename, genefileinfo, bdoptions) {
  headerEntries <- 8L
  offsets <- c(seq(8L, 36L, 4L), 36L)
  return (WriteBinaryDosageHeader4(filename, genefileinfo, bdoptions,
                                   headerEntries, offsets, 0))
}

WriteBinaryDosageHeader43 <- function(filename, genefileinfo, bdoptions) {
  headerEntries <- 4
  offsets <- c(rep(-1L, 5), seq(8L, 20L, 4L))
  return (WriteBinaryDosageHeader4(filename, genefileinfo, bdoptions,
                                   headerEntries, offsets, 0))
}

WriteBinaryDosageHeader44 <- function(filename, genefileinfo, bdoptions) {
  headerEntries <- 4
  offsets <- c(rep(-1L, 5), seq(8L, 20L, 4L))
  return (WriteBinaryDosageHeader4(filename, genefileinfo, bdoptions,
                                   headerEntries, offsets, nrow(genefileinfo$snps)))
}

#***************************************************************************#
#                                                                           #
#                     Writing the Binary Dosage data                        #
#                                                                           #
#***************************************************************************#

#***************************************************************************#
#                        Allocate memory                                    #
#***************************************************************************#

# Allocates memory needed to write binary dosage files
# This is sufficient for all formats
AllocateBinaryDosageWriteMemory <- function(headerinfo) {
  filename <- headerinfo$filename
  format <- headerinfo$additionalinfo$format
  subformat <- headerinfo$additionalinfo$subformat
  headersize <- headerinfo$additionalinfo$headersize
  snpnumber <- integer(1)
  datasize <- integer(nrow(headerinfo$snps))
  us <- integer(2*nrow(headerinfo$samples))
  return(list(filename = filename,
              format = format,
              subformat = subformat,
              headersize = headersize,
              snpnumber = snpnumber,
              datasize = datasize,
              us = us))
}

#***************************************************************************#
#                        Write the data                                     #
#***************************************************************************#

# Write binary dosage data at the end of the file
# Header has already been written
# writeinfo was already created using AllocateBinaryDosageWriteMemory (see above)
WriteBinaryDosageData <- function(dosage, p0, p1, p2, writeinfo) {
  writeFunc <- list(f1 <- c(WriteBinaryDosageData1, WriteBinaryDosageData2),
                    f2 <- c(WriteBinaryDosageData3, WriteBinaryDosageData4),
                    f3 <- c(WriteBinaryDosageData3, WriteBinaryDosageData5, WriteBinaryDosageData3, WriteBinaryDosageData6),
                    f4 <- c(WriteBinaryDosageData3, WriteBinaryDosageData5, WriteBinaryDosageData3, WriteBinaryDosageData6))
  return (writeFunc[[writeinfo$format]][[writeinfo$subformat]](writeinfo, dosage, p0, p1, p2))
}

WriteBinaryDosageData1 <- function(writeinfo, dosage, p0, p1, p2) {
  return (WriteBinaryDosageDataC(writeinfo$filename, dosage, writeinfo$us, 1))
}

WriteBinaryDosageData2 <- function(writeinfo, dosage, p0, p1, p2) {
  return (WriteBinaryP1P2Data(writeinfo$filename, p1, p2, writeinfo$us,  2))
}

WriteBinaryDosageData3 <- function(writeinfo, dosage, p0, p1, p2) {
  return (WriteBinaryDosageDataC(writeinfo$filename, dosage, writeinfo$us, 3))
}

WriteBinaryDosageData4 <- function(writeinfo, dosage, p0, p1, p2) {
  return (WriteBinaryP1P2Data(writeinfo$filename, p1, p2, writeinfo$us, 3))
}

WriteBinaryDosageData5 <- function(writeinfo, dosage, p0, p1, p2) {
  snpnumber <- -1L
  return (WriteBinaryCompressed(writeinfo$filename,
                                dosage, p0, p1, p2,
                                snpnumber,
                                writeinfo$datasize,
                                writeinfo$us))
}

WriteBinaryDosageData6 <- function(writeinfo, dosage, p0, p1, p2) {
  return (WriteBinaryCompressed(writeinfo$filename,
                                dosage, p0, p1, p2,
                                writeinfo$snpnumber,
                                writeinfo$datasize,
                                writeinfo$us))
}

#***************************************************************************#
#                                                                           #
#                     Writing the Binary Dosage indices                     #
#                                                                           #
#***************************************************************************#

# Write binary dosage indices to  the file
# Header has already been written
# funcData was already created using AllocateBinaryDosageWriteMemory (see above)
WriteBinaryDosageIndices <- function(writeinfo) {
  writeFunc <- list(f1 <- c(WriteBinaryDosageIndices1, WriteBinaryDosageIndices1),
                    f2 <- c(WriteBinaryDosageIndices1, WriteBinaryDosageIndices1),
                    f3 <- c(WriteBinaryDosageIndices1, WriteBinaryDosageIndices1, WriteBinaryDosageIndices1, WriteBinaryDosageIndices2),
                    f4 <- c(WriteBinaryDosageIndices1, WriteBinaryDosageIndices1, WriteBinaryDosageIndices1, WriteBinaryDosageIndices2))
  return (writeFunc[[writeinfo$format]][[writeinfo$subformat]](writeinfo))
}

WriteBinaryDosageIndices1 <- function(writeinfo) {
  return (0)
}

WriteBinaryDosageIndices2 <- function(writeinfo) {
  return(WriteBinaryDosageIndicesC(writeinfo$filename, writeinfo$headersize, writeinfo$datasize))
}

#***************************************************************************#
#                                                                           #
#                    Update snp info, aaf, maf, and rsq                     #
#                                                                           #
#***************************************************************************#
#' Calculate alternate allele frequency
#'
#' Routine to calculate the alternate allele frequency given the dosages.
#' Missing values for dosage ignored. This function is used internally and
#' is exported for use in examples.
#'
#' @param dosage Dosage values
#' @param p0 Pr(g=0) - unused
#' @param p1 Pr(g=1) - unused
#' @param p2 Pr(g=2) - unused
#'
#' @return
#' Alternate allele frequency
#' @export
#'
#' @examples
#' # Get information about binary dosage file
#' bdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
#' bdinfo <- getbdinfo(bdfiles = bdfile)
#' snp1 <- getsnp(bdinfo = bdinfo, 1)
#' aaf <- getaaf(snp1$dosage)
getaaf <- function(dosage, p0, p1, p2) {
  return (mean(dosage, na.rm = TRUE) / 2)
}

#' Calculate minor allele frequency
#'
#' Routine to calculate the minor allele frequency given the dosages.
#' Missing values for dosage ignored. This function is used internally and
#' is exported for use in examples. Note: The minor allele in one data set
#' may be different from another data set. This can make comparing minor
#' allele frequencies between data sets nonsensical.
#'
#' @param dosage Dosage values
#' @param p0 Pr(g=0) - unused
#' @param p1 Pr(g=1) - unused
#' @param p2 Pr(g=2) - unused
#'
#' @return
#' Minor allele frequency
#' @export
#'
#' @examples
#' # Get information about binary dosage file
#' bdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
#' bdinfo <- getbdinfo(bdfiles = bdfile)
#' snp1 <- getsnp(bdinfo = bdinfo, 1)
#' maf <- getmaf(snp1$dosage)
getmaf <- function(dosage, p0, p1, p2) {
  aaf <- mean(dosage, na.rm = TRUE) / 2
  maf <- ifelse(aaf > 0.5, 1. - aaf, aaf)
  return (maf)
}

#' Calculate imputation r squared
#'
#' Routine to calculate the imputation r squared given the dosages
#' and Pr(g=2).
#' This is an estimate for the imputation r squared returned from
#' minimac and impute2. The r squared values are calculated slightly
#' differently between the programs. This estimate is based on the
#' method used by minimac. It does well for minor allele frequencies
#' above 5%. This function is used internally and is exported for
#' use in examples.
#'
#' @param dosage Dosage values
#' @param p0 Pr(g=0) - unused
#' @param p1 Pr(g=1) - unused
#' @param p2 Pr(g=2)
#'
#' @return
#' Imputation r squared
#' @export
#'
#' @examples
#' # Get information about binary dosage file
#' bdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
#' bdinfo <- getbdinfo(bdfiles = bdfile)
#' snp1 <- getsnp(bdinfo = bdinfo, 1, dosageonly = FALSE)
#' rsq <- BinaryDosage:::getrsq(snp1$dosage, p2 = snp1$p2)
getrsq <- function(dosage, p0, p1, p2) {
  q <- numeric(2 * length(dosage))
  d <- dosage * dosage - 4 * p2
  d <- ifelse(d < 0 & d > -0.01, 0., d)
  d <- sqrt(d)
  q[1:length(dosage)] <- 0.5 * (dosage - d)
  q[(length(dosage) + 1):(2*length(dosage))] <- 0.5 * (dosage + d)
  q <- q[is.na(q) == FALSE]
  mu <- mean(q)
  sigma <- mean(q*q) - mu * mu
  rsq <- sigma / (mu * (1. - mu))
  return (rsq)
}

updateaaf <- function (bdinfo) {
  if (bdinfo$additionalinfo$subformat < 3)
    headerinfo <- ReadBinaryDosageHeader4A(bdinfo$filename)
  else
    headerinfo <- ReadBinaryDosageHeader4B(bdinfo$filename)

  aaf <- unlist(bdapply(bdinfo, getaaf))

  updatesnpinfo(bdinfo$filename, headerinfo$snps$aafoffset, aaf)
}

updatemaf <- function (bdinfo) {
  if (bdinfo$additionalinfo$subformat < 3)
    headerinfo <- ReadBinaryDosageHeader4A(bdinfo$filename)
  else
    headerinfo <- ReadBinaryDosageHeader4B(bdinfo$filename)

  maf <- unlist(bdapply(bdinfo, getmaf))

  updatesnpinfo(bdinfo$filename, headerinfo$snps$mafoffset, maf)
}

updatersq <- function (bdinfo) {
  if (bdinfo$additionalinfo$subformat < 3)
    headerinfo <- ReadBinaryDosageHeader4A(bdinfo$filename)
  else
    headerinfo <- ReadBinaryDosageHeader4B(bdinfo$filename)

  rsq <- unlist(bdapply(bdinfo, getrsq))

  updatesnpinfo(bdinfo$filename, headerinfo$snps$rsqoffset, rsq)
}

