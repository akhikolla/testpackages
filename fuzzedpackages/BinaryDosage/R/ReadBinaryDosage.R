#***************************************************************************#
#                                                                           #
#                       Reading Binary Dosage files                         #
#                                                                           #
#***************************************************************************#


#***************************************************************************#
#                                                                           #
#                     Reading the Binary Dosage header                      #
#                                                                           #
#***************************************************************************#

# Reads the header for the various formats of the formats
# of the binary dosage file. These vary for all the different
# formats.
ReadBinaryDosageHeader <- function(filename) {
  ReadHeaderFunc <- list(f1 <- c(ReadBinaryDosageHeader11, ReadBinaryDosageHeader12),
                         f2 <- c(ReadBinaryDosageHeader21, ReadBinaryDosageHeader22),
                         f3 <- c(ReadBinaryDosageHeader31, ReadBinaryDosageHeader32, ReadBinaryDosageHeader33, ReadBinaryDosageHeader34),
                         f4 <- c(ReadBinaryDosageHeader41, ReadBinaryDosageHeader42, ReadBinaryDosageHeader43, ReadBinaryDosageHeader44))
  bdformat <- ReadBinaryDosageBaseHeader(filename[1])
  if (is.na(match("error", names(bdformat))) == FALSE)
    stop(bdformat$error)
  if (bdformat$format == 4) {
    if (length(filename) != 1)
      stop("Binary dosage file format 4 does not use family and map files")
  } else {
    if (length(filename) != 3)
      stop("Binary dosage file format 1, 2, and 3 require family and map files")
  }
  return (ReadHeaderFunc[[bdformat$format]][[bdformat$subformat]](filename))
}

#***************************************************************************#
#             Routines to read the subject and SNP info                     #
#***************************************************************************#

# Function to read subject and map files for formats 1, 2, and 3
# @parameter filename - vector of files names for binary dosage files
# first value is the binary dosage file, the second value is the fam file
# and the third is the map file
# @parameter format - format of the binary dosage file
# @parameter subformat - subformat of the binary dosage file
# @parameter headersize - size of the binary dosage header in bytes
# @return - The subject and family data needed to create a bdinfo list
ReadFamAndMapFiles <- function(filename, format, subformat, headersize) {
  fqfilename <- normalizePath(filename[1], winslash = '/')
  samples <- readRDS(filename[2])
  if (all(samples$fid == "") == TRUE)
    usesfid <- FALSE
  else
    usesfid <- TRUE
  snps <- readRDS(filename[3])
  chr1 <- snps$chromosome[1]
  onechr <- all(snps$chromosome == chr1)
  chrlocid <- paste(snps$chromosome, snps$location, sep = ":")
  if(all(snps$snpid == chrlocid) == TRUE)
    snpidformat <- 1
  else {
    chrlocrefaltid <- paste(snps$chromosome, snps$location,
                            snps$reference, snps$alternate, sep = ":")
    if (all(snps$snpid == chrlocrefaltid) == TRUE)
      snpidformat <- 2
    else
      snpidformat <- 0
  }

  additionalinfo <- list(format = format,
                         subformat = subformat,
                         headersize = headersize,
                         numGroups = 1,
                         groups = nrow(samples))
  class(additionalinfo) <- "bdose-info"

  return (list(filename = fqfilename,
               usesfid = usesfid,
               samples = samples,
               onechr = onechr,
               snpidformat = snpidformat,
               snps = snps,
               snpinfo = list(),
               additionalinfo = additionalinfo))
}

# Function to convert the binary dosage header subject and SNP info
# into the bdinfo format
# @parameter filename - vector of file names for the binary dosage
# file
# @parameter header - header information read from the binary dosage
# file
# @paramter format - format of the binary dosage file - should be 4
# @paramter subformat - subformat of the binary dosage file
# @return - The subject and family data needed to create a bdinfo list
Convert4HeaderToBDInfo <- function(filename, header, format, subformat) {
  sid <- unlist(strsplit(header$samples$sidstring, '\t'))
  if (header$samples$fidsize == 0) {
    usesfid <- FALSE
    fid = rep("", header$numsub)
  } else {
    usesfid <- TRUE
    fid <- unlist(strsplit(header$samples$fidstring, '\t'))
  }
  samples <- data.frame(fid, sid, stringsAsFactors = FALSE)

  if (header$snps$chrsize == 0) {
    onechr <- FALSE
    chromosome <- rep("", header$numSNPs)
  } else {
    chromosome <- unlist(strsplit(header$snps$chrstring, '\t'))
    if (length(chromosome) == 1) {
      onechr <- TRUE
      chromosome <- rep(header$snps$chrstring, header$numSNPs)
    } else {
      onechr <- FALSE
    }
  }
  if (length(header$snps$location) == 0)
    location <- rep(0L, header$numSNPs)
  else
    location <- header$snps$location
  if (header$snps$refsize == 0)
    reference <- rep("", header$numSNPs)
  else
    reference <- unlist(strsplit(header$snps$refstring, '\t'))
  if (header$snps$altsize == 0)
    alternate <- rep("", header$numSNPs)
  else
    alternate <- unlist(strsplit(header$snps$altstring, '\t'))
  if (header$snps$snpsize == 0) {
    if (header$snps$refsize == 0) {
      snpid <- paste(chromosome, location, sep = ':')
      snpidformat <- 1
    } else {
      snpid <- paste(chromosome, location, reference, alternate, sep = ':')
      snpidformat <- 2
    }
  } else {
    snpid <- unlist(strsplit(header$snps$snpstring, '\t'))
    snpidformat <- 0
  }
  snps <- data.frame(chromosome, location, snpid, reference, alternate, stringsAsFactors = FALSE)

  snpinfocol <- match(c("aaf", "maf", "avgcall", "rsq"), names(header$snps))
  snpinfocol <- snpinfocol[sapply(header$snps, function(x) length(x) != 0)[snpinfocol]]
  snpinfo <- lapply(header$snps[snpinfocol], matrix, nrow = header$numSNPs, ncol = header$numgroups)

  additionalinfo <- list(format = format,
                         subformat = subformat,
                         headersize = header$dosageoffset,
                         numgroups = header$numgroups,
                         groups = header$groups)
  class(additionalinfo) <- "bdose-info"

  return (list(filename = normalizePath(filename[1], winslash = "/"),
               usesfid = usesfid,
               samples = samples,
               onechr = onechr,
               snpidformat = snpidformat,
               snps = snps,
               snpinfo = snpinfo,
               additionalinfo = additionalinfo))
}

ReadBinaryDosageHeader11 <- function(filename) {
  return (ReadFamAndMapFiles(filename, 1, 1, 8))
}

ReadBinaryDosageHeader12 <- function(filename) {
  return (ReadFamAndMapFiles(filename, 1, 2, 8))
}

ReadBinaryDosageHeader21 <- function(filename) {
  return (ReadFamAndMapFiles(filename, 2, 1, 8))
}

ReadBinaryDosageHeader22 <- function(filename) {
  return (ReadFamAndMapFiles(filename, 2, 2, 8))
}

ReadBinaryDosageHeader31 <- function(filename) {
  bdInfo <- ReadFamAndMapFiles(filename, 3, 1, 12)
  additionalInfo = ReadBinaryDosageHeader3A(filename[1])
  if (additionalInfo$numsub != nrow(bdInfo$samples))
    stop("Subject file does not line up with binary dosage file")
  return (bdInfo)
}

ReadBinaryDosageHeader32 <- function(filename) {
  bdInfo <- ReadFamAndMapFiles(filename, 3, 2, 12)
  additionalInfo = ReadBinaryDosageHeader3A(filename[1])
  if (additionalInfo$numsub != nrow(bdInfo$samples))
    stop("Subject file does not line up with binary dosage file")
  return (bdInfo)
}

ReadBinaryDosageHeader33 <- function(filename) {
  bdInfo <- ReadFamAndMapFiles(filename, 3, 3, 72)
  additionalInfo = ReadBinaryDosageHeader3B(filename[1])
  if (digest(bdInfo$samples) != additionalInfo$md5[1])
    stop("Subject file does not line up with binary dosage file")
  if (digest(bdInfo$snps) != additionalInfo$md5[2])
    stop("Map file does not line up with binary dosage file")
  return (bdInfo)
}

ReadBinaryDosageHeader34 <- function(filename) {
  bdInfo <- ReadFamAndMapFiles(filename, 3, 4, 72)
  bdInfo$additionalinfo$headersize <- bdInfo$additionalinfo$headersize + 4 * nrow(bdInfo$snps)
  additionalInfo = ReadBinaryDosageHeader3B(filename[1])
  if (digest(bdInfo$samples) != additionalInfo$md5[1])
    stop("Subject file does not line up with binary dosage file")
  if (digest(bdInfo$snps) != additionalInfo$md5[2])
    stop("Map file does not line up with binary dosage file")
  return (bdInfo)
}

ReadBinaryDosageHeader41 <- function(filename) {
  header <- ReadBinaryDosageHeader4A(filename[1])
  return (Convert4HeaderToBDInfo(filename, header, 4, 1))
}

ReadBinaryDosageHeader42 <- function(filename) {
  header <- ReadBinaryDosageHeader4A(filename[1])
  return (Convert4HeaderToBDInfo(filename, header, 4, 2))
}

ReadBinaryDosageHeader43 <- function(filename) {
  header <- ReadBinaryDosageHeader4B(filename[1])
  bdInfo <- Convert4HeaderToBDInfo(filename, header, 4, 3)
  return (bdInfo)
}

ReadBinaryDosageHeader44 <- function(filename) {
  header <- ReadBinaryDosageHeader4B(filename[1])
  bdInfo <- Convert4HeaderToBDInfo(filename, header, 4, 4)
  return (bdInfo)
}

#***************************************************************************#
#                                                                           #
#               Getting the indices for a binary dosage file                #
#                                                                           #
#***************************************************************************#

# Gets the file locations for snps in a binary dosage file
ReadBinaryDosageIndices <- function(bdInfo) {
  ReadIndicesFunc <- list(f1 <- c(ReadIndices1, ReadIndices2),
                          f2 <- c(ReadIndices1, ReadIndices2),
                          f3 <- c(ReadIndices1, ReadIndices3, ReadIndices1, ReadIndices4),
                          f4 <- c(ReadIndices1, ReadIndices3, ReadIndices1, ReadIndices4))
  return (ReadIndicesFunc[[bdInfo$additionalinfo$format]][[bdInfo$additionalinfo$subformat]](bdInfo))
}

# Function to set up the datasize and indices vectors when the
# data size is the same for each SNP.
# @parameter numsub - number of subjects in data set
# @parameter numSNPs - number of snps
FixedIndices <- function(numsub, numSNPs, firstIndex, snpsize) {
  datasize <- snpsize * numsub
  indices <- seq(firstIndex, firstIndex + datasize * (numSNPs - 1), datasize)
  return (list(datasize = datasize, indices = indices))
}

# This routine sets up the indices when only the dosages are
# in the binary dosage file. This is simple because the size
# is fixed 2 bytes per subject per SNP
ReadIndices1 <- function(bdInfo) {
  return (FixedIndices(nrow(bdInfo$samples),
                       nrow(bdInfo$snps),
                       bdInfo$additionalinfo$headersize, 2))
}

# This routine sets up the indices when reading formats 1 and 2
# when there is both dosage and genetic probability data. This
# is also simple because the size is fixed to 4 bytes per subject
# per SNP
ReadIndices2 <- function(bdInfo) {
  return (FixedIndices(nrow(bdInfo$samples),
                       nrow(bdInfo$snps),
                       bdInfo$additionalinfo$headersize, 4))
}

# This routine sets up the indices when reading formats 3 and 4,
# subformat 2. The inidices are stored at the start of each SNP's
# dosage data. This is a time consuming operation.
ReadIndices3 <- function(bdInfo) {
  return (ReadBDIndices3C(bdInfo$filename, nrow(bdInfo$snps), bdInfo$additionalinfo$headersize))
}

# This routine sets up the indices when reading formats 3 and 4,
# subformat 4. The inidices are stored in the header and are easy
# to read in.
ReadIndices4 <- function(bdInfo) {
  return (ReadBDIndices4C(bdInfo$filename, nrow(bdInfo$snps), bdInfo$additionalinfo$headersize))
}

#***************************************************************************#
#                                                                           #
#                              Read the data                                #
#                                                                           #
#***************************************************************************#

# Reads a SNP form the various formats
# of the binary dosage file.
ReadBinaryDosageData <- function(bdInfo, snp, d, p0, p1, p2, us) {
  ReadHeaderFunc <- list(f1 <- c(ReadBinaryDosageData1, ReadBinaryDosageData2),
                         f2 <- c(ReadBinaryDosageData3, ReadBinaryDosageData4),
                         f3 <- c(ReadBinaryDosageData3, ReadBinaryDosageData5, ReadBinaryDosageData3, ReadBinaryDosageData5),
                         f4 <- c(ReadBinaryDosageData3, ReadBinaryDosageData5, ReadBinaryDosageData3, ReadBinaryDosageData5))
  return (ReadHeaderFunc[[bdInfo$additionalinfo$format]][[bdInfo$additionalinfo$subformat]](bdInfo, snp, d, p0, p1, p2, us))
}

ReadBinaryDosageData1 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  ReadBinaryDosageDataC(filename = bdInfo$filename,
                        headersize = bdInfo$additionalinfo$headersize,
                        numsub = nrow(bdInfo$samples),
                        snp = snp,
                        dosage = d,
                        us = us,
                        base = 1)
}

ReadBinaryDosageData2 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  ReadBinaryDosageDataP1P2(filename = bdInfo$filename,
                           headersize = bdInfo$additionalinfo$headersize,
                           numsub = nrow(bdInfo$samples),
                           snp = snp,
                           dosage = d,
                           p0 = p0,
                           p1 = p1,
                           p2 = p2,
                           us = us,
                           base = 2)
}

ReadBinaryDosageData3 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  ReadBinaryDosageDataC(filename = bdInfo$filename,
                        headersize = bdInfo$additionalinfo$headersize,
                        numsub = nrow(bdInfo$samples),
                        snp = snp,
                        dosage = d,
                        us = us,
                        base = 3)
}

ReadBinaryDosageData4 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  ReadBinaryDosageDataP1P2(filename = bdInfo$filename,
                           headersize = bdInfo$additionalinfo$headersize,
                           numsub = nrow(bdInfo$samples),
                           snp = snp,
                           dosage = d,
                           p0 = p0,
                           p1 = p1,
                           p2 = p2,
                           us = us,
                           base = 3)
}

ReadBinaryDosageData5 <- function(bdInfo, snp, d, p0, p1, p2, us) {
  return (ReadBinaryDosageDataCompressed(filename = bdInfo$filename,
                                         index = bdInfo$indices[snp],
                                         datasize = bdInfo$datasize[snp],
                                         numsub = nrow(bdInfo$samples),
                                         dosage = d,
                                         p0 = p0,
                                         p1 = p1,
                                         p2 = p2,
                                         us = us))
}

