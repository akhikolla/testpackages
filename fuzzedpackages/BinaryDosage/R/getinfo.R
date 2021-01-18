###########################################################
#                Binary Dosage                            #
###########################################################

#' Get information about a binary dosage file
#'
#' Routine to return information about a binary dosage file.
#' This information is used by other routines to
#' allow for quicker extraction of values from the
#' file.
#'
#' @param bdfiles Vector of file names. The first is the
#' binary dosage data containing the dosages and genetic
#' probabilities. The second file name is the family information
#' file. The third file name is the SNP information file.
#' The family and SNP information files are not used if the
#' binary dosage file is in format 4. For this format the
#' family and SNP information are in the file with the dosages
#' and genetic probabilities.
#'
#' @return List with information about the binary dosage file.
#' This includes family and subject IDs along with
#' a list of the SNPs in the file. Other information needed
#' to read the file is also included.
#'
#' @export
#'
#' @examples
#' vcf1abdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
#' bdinfo <- getbdinfo(bdfiles = vcf1abdfile)
getbdinfo <- function(bdfiles) {
  if (missing(bdfiles) == TRUE)
    stop("No binary dosage files specified")
  if (is.character(bdfiles) == FALSE)
    stop("bdfiles must be a character vector")
  if (length(bdfiles) != 1 & length(bdfiles) != 3)
    stop("bdfiles must be a character vector of length 1 or 3")
  if (bdfiles[1] == "")
    stop("No binary dosage file specified")
  if (length(bdfiles) == 3) {
    if (bdfiles[2] == "" & bdfiles[3] == "") {
      bdfiles <- bdfiles[1]
    } else if (bdfiles[2] == "" | bdfiles[3] == "") {
      stop("bdfiles contains empty strings")
    }
  }
  headerinfo <- ReadBinaryDosageHeader(bdfiles)
  indices <- ReadBinaryDosageIndices(headerinfo)
  headerinfo$datasize <- indices$datasize
  headerinfo$indices <- indices$indices
  class(headerinfo) <- c("genetic-info")

  return (headerinfo)
}

###########################################################
#                        VCF                              #
###########################################################
summarizevcfadditionalinfo <- function(x) {
  if (length(unique(x)) != 1)
    return (x)
  if (x[1] == '.')
    return (character(0))
  return (x[1])
}

readminimacinfofile <- function(filename) {
  addinfo <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
  if (ncol(addinfo) != 13)
    stop("Error reading information file - Wrong number of columns")
  if (all(colnames(addinfo) == c("SNP", "REF.0.", "ALT.1.", "ALT_Frq", "MAF",
                                 "AvgCall", "Rsq", "Genotyped", "LooRsq",
                                 "EmpR", "EmpRsq", "Dose0", "Dose1")) == FALSE)
    stop("Error reading information file - Wrong column names")
  return(addinfo)
}

#' Get information about a vcf file
#'
#' Routine to return information about a vcf file.
#' This information is used by other routines to
#' allow for quicker extraction of values from the
#' file.
#'
#' @param vcffiles A vector of file names.
#' The first is the name of the vcf file. The
#' second is name of the file that contains information
#' about the imputation of the SNPs. This file is produced
#' by minimac 3 and 4.
#' @param gz Indicator if VCF file is compressed using gzip.
#' Default value is FALSE.
#' @param index Indicator if file should be indexed. This
#' allows for faster reading of the file. Indexing a gzipped
#' file is not supported.
#' Default value is TRUE.
#' @param snpidformat The format that the SNP ID will be saved as.
#' 0 - same as in the VCF file
#' 1 - chromosome:location
#' 2 - chromosome:location:referenceallele:alternateallele
#' If snpidformat is 1 and the VCF file uses format 2, an error is
#' generated. Default value is 0.
#'
#' @return List containing information about the VCF file
#' to include file name, subject IDs, and information about
#' the SNPs. Indices for faster reading will be included
#' if index is set to TRUE
#'
#' @export
#'
#' @examples
#' # Get file names of th vcf and infromation file
#' vcf1afile <- system.file("extdata", "set1a.vcf", package = "BinaryDosage")
#' vcf1ainfo <- system.file("extdata", "set1a.info", package = "BinaryDosage")
#'
#' # Get the information about the vcf file
#' vcf1ainfo <- getvcfinfo(vcffiles = c(vcf1afile, vcf1ainfo))
getvcfinfo <- function(vcffiles,
                       gz = FALSE,
                       index = TRUE,
                       snpidformat = 0L) {
  if (missing(vcffiles) == TRUE)
    stop("No VCF file specified")
  if (is.character(vcffiles) == FALSE)
    stop("vcfiles must be a character value")
  if (length(vcffiles) != 1 & length(vcffiles) != 2)
    stop("vcffiles must be a character vector of length 1 or 2")
  filename = vcffiles[1]
  if (length(vcffiles) != 1)
    infofile <- vcffiles[2]
  else
    infofile <- ""
  if (filename == "")
    stop("No VCF file specified")

  if (is.logical(gz) == FALSE)
    stop("gz must be a logical value")
  if (length(gz) != 1)
    stop("gz must be a logical vector of length 1")

  if (is.logical(index) == FALSE)
    stop("index must be a logical value")
  if (length(index) != 1)
    stop("index must be a logical vector of length 1")
  if (gz == TRUE && index == TRUE)
    stop("Indexing gzipped files is not supported.")

  if (is.numeric(snpidformat) == FALSE)
    stop("snpidformat must be an integer value")
  if (length(snpidformat) != 1)
    stop("snpidformat must be an interger vector of length 1")
  if (floor(snpidformat) != snpidformat)
    stop("snpidformat must be an integer value")
  snpidformat <- as.integer(snpidformat)
  if (snpidformat < 0 || snpidformat > 2)
    stop("snpidformat must have a value of 0, 1, or 2")

  if (gz == FALSE) {
    con <- file(filename, "r")
  } else {
    con <- gzfile(filename, "r")
  }

  fqfilename <- normalizePath(filename, winslash = '/')

  headerlines <- 1L
  headersize <- -1L
  while (TRUE) {
    currentpos <- seek(con, origin = "current")
    line <- readLines(con, n = 1)
    if (substr(line, 1, 1) != '#') {
      close(con)
      stop("Error processing header")
    }
    if (substr(line, 2, 2) != '#') {
      x <- unlist(strsplit(line, "\t"))
      if (x[1] != "#CHROM") {
        close(con)
        stop("Error processing header")
      }
      x[1] = "CHROM"
      begindata <- seek(con, origin = "current")
      break
    }
    headerlines <- headerlines + 1L
  }
  if (index == TRUE)
    headersize <- begindata
  close(con)

  if (all(x[1:9] == c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")) == FALSE)
    stop("Column names incorrect")

  samples = data.frame(fid = rep("", length(x) - 9L),
                       sid = x[10:length(x)],
                       stringsAsFactors = FALSE)
  usesfid = FALSE

  coltypes = c("character", "integer", rep("character", 7),
               rep("NULL", nrow(samples)))
  snps <- read.table(filename,
                     skip = headerlines,
                     colClasses = coltypes,
                     stringsAsFactors = FALSE)
  colnames(snps) <- c("chromosome", "location", "snpid",
                      "reference", "alternate", "quality",
                      "filter", "info", "format")

  vcfinfo <- as.list(snps[,6:9])
  summaryinfo <- lapply(vcfinfo, summarizevcfadditionalinfo)
  datacolumns <- data.frame(numcolumns = rep(0L, length(summaryinfo$format)),
                            dosage = rep(0L, length(summaryinfo$format)),
                            genotypeprob = rep(0L, length(summaryinfo$format)),
                            genotype = rep(0L, length(summaryinfo$format)),
                            stringsAsFactors = FALSE)
  for (i in 1:length(summaryinfo$format)) {
    formatsplit <- unlist(strsplit(summaryinfo$format[i], split = ':'))
    datacolumns$numcolumns[i] <- length(formatsplit)
    datacolumns$dosage[i] <- match("DS", formatsplit)
    datacolumns$genotypeprob[i] <- match("GP", formatsplit)
    datacolumns$genotype[i] <- match("GT", formatsplit)
  }
  additionalinfo <- list(gzipped = gz,
                         headerlines = headerlines,
                         headersize = headersize,
                         quality = summaryinfo$quality,
                         filter = summaryinfo$filter,
                         info = summaryinfo$info,
                         format = summaryinfo$format,
                         datacolumns = datacolumns)
  class(additionalinfo) <- "vcf-info"
  rm(vcfinfo)
  rm(datacolumns)

  if (infofile == "") {
    snpinfo <- list()
  } else {
    minimacinfo <- readminimacinfofile(infofile)
    if (nrow(minimacinfo) == nrow(snps)) {
      if (all(minimacinfo$SNP == snps$snpid) == TRUE &
          all(minimacinfo$REF.0. == snps$reference) == TRUE &
          all(minimacinfo$ALT.1. == snps$alternate) == TRUE) {
        snpinfo <- list(aaf = as.matrix(minimacinfo$ALT_Frq),
                        maf = as.matrix(minimacinfo$MAF),
                        avgcall = as.matrix(minimacinfo$AvgCall),
                        rsq = as.matrix(minimacinfo$Rsq))
      } else {
        stop("Infromation file does not line up with VCF file - different SNPs")
      }
    } else {
      stop("Information file does not line up with VCF file - different number of SNPs")
    }
  }

  snps <- snps[,1:5]
  chr1 <- snps$chromosome[1]
  oneChr <- all(snps$chromosome == chr1)
  chrlocid <- paste(snps$chromosome, snps$location, sep = ":")
  vcfsnpformat1 <- all(snps$snpid == chrlocid)
  chrlocrefaltid <- paste(snps$chromosome, snps$location,
                          snps$reference, snps$alternate, sep = ":")
  vcfsnpformat2 <- all(snps$snpid == chrlocrefaltid)
  if (snpidformat == 0) {
    if (vcfsnpformat1 == TRUE) {
      snps$snpid <- chrlocid
      snpidformat <- 1L
    } else if (vcfsnpformat2 == TRUE) {
      snps$snpid <- chrlocrefaltid
      snpidformat <- 2L
    }
  } else if (snpidformat == 1) {
    if (vcfsnpformat2 == TRUE)
      stop ("snpidformat 1 specified but VCF file uses snpidformat 2")
    if (vcfsnpformat1 == FALSE)
      snps$snpid <- chrlocid
  } else if (snpidformat == 2) {
    if (vcfsnpformat2 == FALSE)
      snps$snpid <- chrlocrefaltid
  }

  if (index == TRUE) {
    datasize <- integer(nrow(snps))
    indices <- numeric(nrow(snps))
    if (gz == FALSE) {
      x <- GetLineLocations(filename)
      indices <- x[(headerlines + 1):(length(x) - 1)]
      for (i in 1:length(datasize))
        datasize[i] <- x[headerlines + i + 1] - x[headerlines + i]
    }
  } else {
    datasize <- integer(0)
    indices <- numeric(0)
  }

  retval = list(filename = fqfilename,
                usesfid = usesfid,
                samples = samples,
                onechr = oneChr,
                snpidformat = snpidformat,
                snps = snps,
                snpinfo = snpinfo,
                datasize = datasize,
                indices = indices,
                additionalinfo = additionalinfo)
  class(retval) <- c("genetic-info")
  return (retval)
}

###########################################################
#                 GEN (Impute2)                           #
###########################################################

#' Get information about a gen, impute2, file
#'
#' Routine to return information about a gen file.
#' This information is used by other routines to
#' allow for quicker extraction of values from the
#' file.
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
#' @param index Indicator if file should be indexed. This
#' allows for faster reading of the file. Indexing a gzipped
#' file is not supported.
#' Default value is TRUE.
#' @param snpidformat Format to change the snpid to.
#' 0 indicates to use the snpid format in the file.
#' 1 indicates to change the snpid into chromosome:location,
#' 2 indicates to change the snpid into chromosome:location:referenceallele:alternateallele,
#' 3 indicates to change the snpid into chromosome:location_referenceallele_alternateallele,
#' Default value is 0.
#' @param sep Separators used in the gen file and sample files,
#' respectively. If only value is provided it is used for both
#' files. Default value is c(`"\t"`, `"\t"`)
#'
#' @return List with information about the gen file.
#' This includes family and subject IDs along with
#' a list of the SNPs in the file. Other information needed
#' to read the file is also included.
#' @export
#'
#' @examples
#' # Get file names of th gen and sample file
#' gen3afile <- system.file("extdata", "set3a.imp", package = "BinaryDosage")
#' gen3ainfo <- system.file("extdata", "set3a.sample", package = "BinaryDosage")
#'
#' # Get the information about the gen file
#' geninfo <- getgeninfo(genfiles = c(gen3afile, gen3ainfo),
#'                       snpcolumns = c(0L, 2L:5L))
getgeninfo <- function(genfiles,
                       snpcolumns = 1L:5L,
                       startcolumn = 6L,
                       impformat = 3L,
                       chromosome = character(),
                       header = c(FALSE, TRUE),
                       gz = FALSE,
                       index = TRUE,
                       snpidformat = 0L,
                       sep = c("\t", "\t")) {
  if (missing(genfiles) == TRUE)
    stop("No gen file specified")
  if (is.character(genfiles) == FALSE)
    stop("genfiles must be a character value")
  if (length(genfiles) != 1 & length(genfiles) != 2)
    stop("genfiles must be a character vector of length 1 or 2")
  genfile <- genfiles[1]
  if (length(genfiles) == 1) {
    samplefile <- character()
  } else {
    samplefile <- genfiles[2]
    if (samplefile == "")
      samplefile <- character()
  }
  if (genfile == "")
    stop("No gen file specified")

  if (is.integer(snpcolumns) == FALSE)
    stop("snpcolumns must be an integer vector")
  if (length(snpcolumns) != 5)
    stop("snpcolumns must be an integer vector of length 5")
  if (min(snpcolumns[2:5]) < 1)
    stop("snpcolumns values other than chromosome must be positive integers")
  if (snpcolumns[1] < -1)
    stop("snpcolumns chromosome value must be -1, or a non-negative integer")

  if (is.integer(startcolumn) == FALSE)
    stop("startcolumn must be an integer value")
  if (length(startcolumn) != 1)
    stop("startcolumn must be an integer vector of length 1")
  if (startcolumn < 1)
    stop("startcolumn must be a positive integer")
  if (startcolumn <= max(snpcolumns))
    stop("startcolumn value must be larger than any value in snpcolumns")

  if (is.integer(impformat) == FALSE)
    stop("impformat must be an integer value")
  if (length(impformat) != 1)
    stop("impformat must be an integer vector of length 1")
  if (impformat < 1 | impformat > 3)
    stop("impformat must have a value of 1, 2, or 3")

  if (is.character(chromosome) == FALSE)
    stop("chromosome must be a character variable")
  if (length(chromosome) > 1)
    stop("chromosome must be a character vector of length 0 or 1")
  if (length(chromosome) == 1) {
    if (chromosome == "")
      chromosome = character()
  }
  if (length(chromosome) == 0) {
    if (snpcolumns[1] == -1)
      stop("No chromosome column or value provided")
  } else {
    if (snpcolumns[1] > -1)
      stop("Both chromosome column and chromosome value provided")
  }

  if (is.logical(header) == FALSE)
    stop("header must be a logical value")
  if (length(header) !=  1 & length(header) != 2)
    stop("header must be a logical vector of length 1 or 2")
  if (length(header) == 1)
    header = c(header, TRUE)
  if (header[1] == FALSE) {
    if (length(samplefile) == 0)
      stop("File has no header and no sample file is provided")
  } else {
    if (length(samplefile) > 0)
      stop("header = TRUE and a sample file is provided")
  }

  if (is.logical(gz) == FALSE)
    stop("gz must be a logical value")
  if (length(gz) != 1)
    stop("gz must be a logical vector of length 1")

  if (is.logical(index) == FALSE)
    stop("index must be a logical value")
  if (length(index) != 1)
    stop("index must be a logical vector of length 1")
  if (gz == TRUE & index == TRUE)
    stop("Indexing gzipped files is not supported.")

  if (is.numeric(snpidformat) == FALSE)
    stop("snpidformat must be an integer value")
  if (length(snpidformat) != 1)
    stop("snpidformat must be an interger vector of length 1")
  if (floor(snpidformat) != snpidformat)
    stop("snpidformat must be an integer value")
  snpidformat <- as.integer(snpidformat)
  if (snpidformat < 0 || snpidformat > 3)
    stop("snpidformat must have a value of 0, 1, 2, or 3")

  if (is.character(sep) == FALSE)
    stop("sep must be a character value")
  if (length(sep) != 1 & length(sep) != 2)
    stop("sep must be a character vector of length 1 or 2")
  if (length(sep) == 1)
    sep = c(sep, sep)
  if (sep[1] == "" | sep[2] == "")
    stop("sep values cannot be empty strings")


  if (header[1] == TRUE) {
    if (gz == TRUE)
      filecon <- gzfile(genfile, "r")
    else
      filecon <- file(genfile, "r")
    headerline <- readLines(filecon, 1)
    close(filecon)
    headervalues <- unlist(strsplit(headerline, sep[1]))
    if (length(headervalues) < startcolumn)
      stop("Number of values in header less than startcolumn")
    headervalues <- headervalues[startcolumn:length(headervalues)]
    if (length(headervalues) != 2 * floor(length(headervalues) / 2))
      stop("Odd number of values for family and subject ID")
    fid <- headervalues[seq(1,length(headervalues) - 1, 2)]
    iid <- headervalues[seq(2,length(headervalues), 2)]
    samples <- data.frame(fid = fid,
                          sid = iid,
                          stringsAsFactors = FALSE)
  } else {
    samples <- read.table(samplefile,
                          header = header[2],
                          sep = sep[2],
                          stringsAsFactors = FALSE)
    if (ncol(samples) == 1)
      samples[,2] <- samples[,1]
    samples <- samples[,1:2]
    samples[,1] <- as.character(samples[,1])
    samples[,2] <- as.character(samples[,2])
    colnames(samples) <- c("fid", "sid")
    if (samples[1,1] == "0" & samples[1,2] == "0")
      samples <- samples[2:nrow(samples),]
  }
  if (all(samples$fid == samples$sid)) {
    usesfid <- FALSE
    samples$fid <- ""
  } else {
    usesfid <- TRUE
  }

  coltypes = rep("NULL", impformat * nrow(samples) + (startcolumn - 1))
  if (snpcolumns[1] > 0)
    coltypes[snpcolumns[1]] <- "character"
  coltypes[snpcolumns[c(2, 4, 5)]] <- "character"
  coltypes[snpcolumns[3]] <- "integer"
  headersize <- 0
  if (header[1] == TRUE)
    headersize <- 1
  snps <- read.table(genfile,
                     skip = headersize,
                     colClasses = coltypes,
                     sep = sep[1],
                     stringsAsFactors = FALSE)
  if (snpcolumns[1] == -1) {
    snps$chromosome <- chromosome
    snpcolumns[1] <- startcolumn
  } else if (snpcolumns[1] == 0) {
    chromosomeid <- unlist(strsplit(snps[,1], ':'))
    snpidentries <- length(chromosomeid) / nrow(snps)
    snps$chromosome <- chromosomeid[seq(1, length(chromosomeid) + 1 - snpidentries, snpidentries)]
    snpcolumns[1] <- startcolumn
  }
  x <- snpcolumns[2]
  snpcolumns[2] <- snpcolumns[3]
  snpcolumns[3] <- x
  snps <- snps[,order(order(snpcolumns))]
  colnames(snps) <- c("chromosome", "location", "snpid", "reference", "alternate")

  chr1 <- snps$chromosome[1]
  onechr <- FALSE
  if (all(snps$chromosome == chr1))
    onechr <- TRUE

  if (snpidformat == 0) {
    chrlocid <- unlist(paste(snps$chromosome, snps$location, sep = ':'))
    if (all(chrlocid == snps$snpid)) {
      snpidformat <- 1
    } else {
      chrlocrefaltid <- unlist(paste(snps$chromosome, snps$location,
                                     snps$reference, snps$alternate, sep = ':'))
      if (all(chrlocrefaltid == snps$snpid)) {
        snpidformat <- 2
      } else {
        chrlocrefaltid <- unlist(paste(snps$chromosome, snps$location, sep = ':'))
        chrlocrefaltid <- unlist(paste(chrlocrefaltid, snps$reference,
                                       snps$alternate, sep = '_'))
        if (all(chrlocrefaltid == snps$snpid))
          snpidformat <- 3
      }
    }
  } else if (snpidformat == 1) {
    chrlocrefaltid <- unlist(paste(snps$chromosome, snps$location,
                                   snps$reference, snps$alternate, sep = ':'))
    if (all(chrlocrefaltid == snps$snpid)) {
      stop("snpidformat 1 specified but GEN file uses snpidformat 2")
    }
    snps$snpid <- unlist(paste(snps$chromosome, snps$location, sep = ':'))
  } else if (snpidformat == 2) {
    snps$snpid <- unlist(paste(snps$chromosome, snps$location,
                               snps$reference, snps$alternate, sep = ':'))
  } else {
    chrlocrefaltid <- unlist(paste(snps$chromosome, snps$location, sep = ':'))
    snps$snpid <- unlist(paste(chrlocrefaltid, snps$reference,
                               snps$alternate, sep = '_'))
  }
  snpinfo <- list()
  if (index == FALSE) {
    datasize <- integer(0)
    indices <- numeric(0)
  } else {
    datasize <- integer(nrow(snps))
    indices <- numeric(nrow(snps))
    if (gz == FALSE) {
      x <- GetLineLocations(genfile)
      if (header[1] == TRUE)
        headerlines <- 1
      else
        headerlines <- 0
      indices <- x[(headerlines + 1):(length(x) - 1)]
      for (i in 1:length(datasize))
        datasize[i] <- x[headerlines + i + 1] - x[headerlines + i]
    }
  }

  additionalinfo <- list(gzipped = gz,
                         headersize = headersize,
                         format = impformat,
                         startcolumn = startcolumn,
                         sep = sep[1])
  class(additionalinfo) <- "gen-info"

  retval <- list(filename = normalizePath(genfile, winslash = '/'),
                 usesfid = usesfid,
                 samples = samples,
                 onechr = onechr,
                 snpidformat = snpidformat,
                 snps = snps,
                 snpinfo = snpinfo,
                 datasize = datasize,
                 indices = indices,
                 additionalinfo = additionalinfo)
  class(retval) <- "genetic-info"

  return(retval)
}
