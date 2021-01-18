
#   Mega2R: Mega2 for R.
#
#   Copyright 2017-2019, University of Pittsburgh. All Rights Reserved.
#
#   Contributors to Mega2R: Robert V. Baron and Daniel E. Weeks.
#
#   This file is part of the Mega2R program, which is free software; you
#   can redistribute it and/or modify it under the terms of the GNU
#   General Public License as published by the Free Software Foundation,
#   either version 2 of the License, or (at your option) any later
#   version.
#
#   Mega2R is distributed in the hope that it will be useful, but WITHOUT
#   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#   for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
#   For further information contact:
#       Daniel E. Weeks
#       e-mail: weeks@pitt.edu
#
# ===========================================================================

#library(GenABEL)


#' generate gwaa.data-class object from a \bold{Mega2R} database
#'
#' @description
#'  Call the \emph{Mega2R} functions to: create a .tped file, a .tfam file and a .phe file.
#'  Then call the GenABEL functions to process these files: the .tped and the .tfam
#'  file are processed by \code{convert.snp.tped} to produce a tped.raw file.  The latter
#'  is combined with a .phe (phenotype) file by \code{load.gwaa.data} to create a gwaa.data-class
#'  object in memory.  All these files are deleted when the exits.
#'
#' @param markers data frame of markers to be processed
#'
#' @param mapno specify which map index to use for physical distances
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return gwaa.data-class object generated from the Mega2R database
#'
#' @export
#'
#' @examples
#'\dontrun{
#' require("GenABEL")
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db)
#' seqsimgwaa = Mega2GenABEL(markers=ENV$markers[1:10,])
#'
#' str(seqsimgwaa)
#' head(summary(seqsimgwaa))
#'}
#'
Mega2GenABEL = function (markers = NULL, mapno = 0, envir = ENV) {
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)

    genABEL.convert.snp.tped = get0("convert.snp.tped", inherits = TRUE)
    genABEL.load.gwaa.data   = get0("load.gwaa.data",   inherits = TRUE)
    if ( is.null(genABEL.convert.snp.tped) ||
         is.null(genABEL.load.gwaa.data)) {
      warning("genABEL has been archived and is not available\n")
      return (NULL)
    }
## print(system.time ({
    if (is.null(markers)) markers = envir$markers

    prefix = file.path(tempdir(), "Mega2GenABEL")

    on.exit(Mega2GenABELClean())

    mkGenABELtped(prefix, markers, mapno = mapno, envir)

    mkGenABELtfam(prefix, envir)

    genABEL.convert.snp.tped(tpedfile=paste0(prefix,".tped"),
                            tfamfile=paste0(prefix,".tfam"),
                            outfile=paste0(prefix, "tped.raw"),
                            strand="u")

    file = paste0(prefix, ".phe")
    unlink(file)
    out = mkGenABELphenotype(envir)
    write.table(out, file=file, sep="\t", quote=FALSE,
                row.names=FALSE, col.names=TRUE)


#x  return (load.gwaa.data(phenofile=paste0(prefix,".phe"),
    ans = genABEL.load.gwaa.data(phenofile=paste0(prefix,".phe"),
                                 genofile=paste0(prefix, "tped.raw"),
                                 force = TRUE)
}

#' delete temporary PLINK tped files processed by GenABEL
#'
#' @description
#'  Delete the PLINK .tped files:  a .tped file, a .tfam file and a .phe file and
#'  the GenABEL tped.raw file.
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#'
#' Mega2GenABELClean()
#'}
Mega2GenABELClean = function () {

    prefix = file.path(tempdir(), "Mega2GenABEL")

    unlink(paste0(prefix, ".tped"))
    unlink(paste0(prefix, ".tfam"))
    unlink(paste0(prefix, "tped.raw"))
    unlink(paste0(prefix, ".phe"))
}

#' generate gwaa.data-class object
#'
#' @description
#'  create a gwaa.data-class object from the data frames in a Mega2 environment.  This
#'  function is a front end that eventually calls a C++ Rcpp function that reads
#'  the genotype data in Mega2 compressed format and converts it to the GenABEL
#'  compressed format.  The results of \code{Mega2ENVGenABEL} are/should be the same
#'  as \code{Mega2GenABEL}, but the calculation is much faster, typically a factor of 10 to 20.
#'
#' @param markers data frame of markers to be processed
#'
#' @param force pass value to gwaa conversion function
#'
#' @param makemap pass value to gwaa conversion function
#'
#' @param sort pass value to gwaa conversion function
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return gwaa.data-class object created from Mega2R database
#'
#' @export
#'
#' @examples
#'\dontrun{
#' require("GenABEL")
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db)
#' gwaa = Mega2ENVGenABEL(markers=ENV$markers[1:10,])
#'
#' str(gwaa)
#' head(summary(gwaa))
#'}
#'
Mega2ENVGenABEL = function (markers = NULL, force = TRUE, makemap = FALSE,
                         sort = TRUE, envir = ENV) {
#browser()
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)

    if (is.null(markers)) markers = envir$markers

## If there is any chance that load.gwaa.data.mega2 exists in GenABEL, ... but that
## may never happen 3/23/18.  Thus the code may always live in mega2genabelexternal.R.
## 
##  fn = get0("load.gwaa.data.mega2", inherits=TRUE, ifnotfound = NULL)
    fn = NULL
    if (is.null(fn)) {
        V2.gwaa.data.mega2(markers, force = force,
            makemap = makemap, sort = sort, envir = envir)
    } else {
        fn(markers, force = force, makemap = makemap, sort = sort, envir = envir)
    }
}

#' generate a PLINK TPED file for GenABEL
#'
#' @description
#'  Generate a PLINK TPED file from the specified Mega2 SQLite database.  The file is named "prefix".tped
#'  If the markers argument is.null(), the entire envir$markers set is include; otherwise the markers argument MUST
#'  be a subset of the envir$markers data.frame -- same columns, but fewer rows.
#'
#' @param prefix prefix for .tped file name
#'
#' @param markers markers selected to be in output file
#'
#' @param mapno specify which map index to use for genetic distances
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return None
#'
#' @importFrom utils write.table
#"
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db)
#' mkGenABELtped("foo", NULL, 0, ENV)
#'}
#
mkGenABELtped = function(prefix, markers=NULL, mapno = 0, envir) {
    file = paste0(prefix, ".tped")

    unlink(file)

    if (is.null(markers)) markers = envir$markers

    allele_table = envir$allele_table[envir$allele_table$locus_link %in% markers$locus_link,]
    map_table = envir$map_table[envir$map_table$marker %in% markers$locus_link,]
    M = nrow(markers)
    C = 1000

    block = data.frame(matrix(0, nrow = C, ncol = nrow(envir$fam) + 4), stringsAsFactors=FALSE)
    blockcol = ncol(block)

    ppl = paste0(envir$fam$PedPre, "_", envir$fam$PerPre)
#   ppl = rbind(ppl, ppl)
    names(block) = c("#CHROM", "ID", "GEN", "POS", ppl)

    j = 0
    while (TRUE) {
        if (M <= 0) break
        R = ((j*C+1):(j*C + C))
        L = length(R)
        if (M < L) {
            L = M
            R = ((j*C+1):(j*C+L))
        }
        M = M - L
        BR = 1:L

## print(system.time ({
        GPos = map_table[map_table$map==mapno, c("position")][R]
        GPosPos = sprintf("%.2f", GPos)

        block[BR , 1] = markers[R , 4]
        block[BR , 2] = markers[R , 3]
        block[BR , 3] = GPosPos
        block[BR , 4] = markers[R , 5]

        cr = getgenotypes(markers[R, ], sepstr = " ", envir = envir ) # 7.17%
        a1 = t(cr)                                      # 0.86%

#       di = dim(a1)
#       line = paste(substr(a1, start=1, stop=1), substr(a1, start=2, stop=2))
#       dim(line) = di
#       block[BR, 5:blockcol] = line

        block[BR, 5:blockcol] = a1

##      }))

## print(system.time ({
        write.table(block[BR, ], file=file, sep="\t", quote=FALSE,     # 48.49%
                    append=TRUE, row.names=FALSE, col.names=FALSE)
## }))
        j = j + 1
##      if (envir$verbose) message(".", appendLF = FALSE)
    }
}

#' generate required fam family for PLINK TPED (.tfam) file
#'
#' @description
#'  Generate the six column .tfam file used with the .tped file.  Note: Only the person id
#'  column appears to be used by GenABEL.
#'
#' @param prefix prefix for generated file name
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return None
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db)
#' mkGenABELtfam("foo", ENV)
#'}
#
mkGenABELtfam = function (prefix, envir) {
    file = paste0(prefix, ".tfam")

    fam = envir$fam
    fam[ , "PerPre"] = paste(fam[ , "PedPre"], fam[ , "PerPre"], sep="_")

    write.table(fam[, -(1:2)], file=file, sep="\t", quote=FALSE,
                row.names=FALSE, col.names=FALSE)
}

#' generate required PLINK (.phe) file
#'
#' @description
#'  Generate the .phe (PLINK phenotype) file needed by GenAbel.  The person
#'  must match that specified in the .tfam file
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return out phenotype data frame
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db)
#' mkGenABELphenotype(ENV)
#'}
#
mkGenABELphenotype = function (envir) {

    fam = envir$fam
    fam$id = paste(fam[ , "PedPre"], fam[ , "PerPre"], sep="_")

# linkage.h:    TYPE_UNSET, QUANT, AFFECTION, BINARY, NUMBERED, XLINKED, YLINKED
#                        0      1          2       3         4        5        6

    out = data.frame(fam$id, stringsAsFactors=FALSE)
    names(out) = "id"
# This (NA / 1 male / 0 female) is different from linkage: NA / 1 male / 2 female
    out$sex = c(NA, 1, 0)[fam$Sex + 1]
    hdr = c("id", "sex")

    phenotype_table = envir$phenotype_table

    raw = unlist(envir$phenotype_table[,4])
    raw = matrix(raw, ncol=8, byrow=T)
    nrows     = nrow(raw)
    nrowpheno = nrow(out)

    for (i in 1:envir$PhenoCnt) {
        hdr = c(hdr, envir$locus_table[i, 2]) # 2 == LocusName

# phenotype_table contains a blob which is a list of entries.  An entry is either an 8 byte
#  double for quant, or two 4 byte ints for affect
        if (envir$locus_table[i, 3] == 2) {              # 3 == Type === AFFECTION
            col = vector("integer", nrowpheno)
            for (j in 1:nrowpheno) {
                col[j] = readBin(raw[envir$PhenoCnt*(j-1)+i, 1:4], integer(), n=1, size=4)
            }
#           col[col==0] = NA
            out$col = col
            names(out) = hdr
        } else if (envir$locus_table[i, 3] == 1) {       # 3 == Type === QUANT
            col = vector("numeric", nrowpheno)
            for (j in 1:nrowpheno) {
                col[j] = readBin(raw[envir$PhenoCnt*(j-1)+i, 1:8], numeric(), n=1, size=8)
            }
            col[col==-99] = NA
            out$col = col
            names(out) = hdr
        }
    }

    out
}

#' generate GenABEL coding vector
#'
#' @description
#'  Each element of the coding vector is the index of the corresponding genotype in the GenABEL:::alleleID.codes() array.
#'  The genotypes are ordered so the allele with the greater frequency appears first.
#'
#' @param markers data frame of markers to be processed
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return None
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' mkGenABELcoding(envir)
#'}
mkGenABELcoding = function(markers = NULL, envir) {
    if (is.null(markers)) markers = envir$markers

    allele_table = envir$allele_table[envir$allele_table$locus_link %in% markers$locus_link,]
    mm = merge(x=allele_table[allele_table$indexX == 1,],
               y=allele_table[allele_table$indexX == 2,],
               by="locus_link")
    nn=ifelse(mm$Frequency.x > mm$Frequency.y,             # before cleaning
              paste0(mm$AlleleName.x, mm$AlleleName.y),
              paste0(mm$AlleleName.y, mm$AlleleName.x))
    envir$xGTy = mm$Frequency.x > mm$Frequency.y
##  nn=ifelse(Freq.x > (1-Freq.x),
##            paste0(mm$AlleleName.x, mm$AlleleName.y),
##            paste0(mm$AlleleName.y, mm$AlleleName.x))
##  envir$xGTy = Freq.x > (1-Freq.x)

    if (any(mm$Frequency.x == .5) && FALSE) {
        w = which(mm$Frequency.x == .5)
        xx = getgenotypesraw(markers[w,], envir)
        yy = xx[1,] == 65537
        nn[w] = ifelse(yy, paste0(mm$AlleleName.x, mm$AlleleName.y)[w], paste0(mm$AlleleName.y, mm$AlleleName.x)[w])
##      if (envir$MARKER_SCHEME == 1) {
##          ms = envir$markerscheme_table[envir$markerscheme_table$key %in% markers$locus_link,]
##          print(mm[w,])
##          print(ms[w,])
##          print(envir$markers[w,])
##          ms1 = ms$allele1[w]
##          nn[w] = ifelse(ms1 == 1, paste0(mm$AlleleName.y[w], mm$AlleleName.x[w]),
##                                   paste0(mm$AlleleName.x[w], mm$AlleleName.y[w]))
##      } else if (envir$MARKER_SCHEME == 2) {
##          print(mm[w,])
##      }
        print(mm[w,])
        print(nn[w])
    }

    nn[mm$Frequency.x == 0 & mm$Frequency.y == 0] = '12'
##  nn[Freq.x == 2] = '12'

    fx = mm$Frequency.x == 1
    nn[fx] = paste0(mm[fx, "AlleleName.x"], mm[fx, "AlleleName.x"])

    fy = mm$Frequency.x == 0
    nn[fy] = paste0(mm[fy, "AlleleName.y"], mm[fy, "AlleleName.y"])

#this is what genabel does for 0/0
    nn[nn=="00"] = "12"
    cc = alleleID.codes()
    ar = match(nn, cc, nomatch=NA)
    if (any(is.na(ar))) {
        warning("Problems mapping genotypes to GenABEL encoding ", mm[is.na(ar)], "/n")
        ar[is.na(ar)] = which(cc == '12')
    }
    as.raw(ar)
}

#' generate GenABEL compressed genotype matrix
#'
#' @description
#'  The matrix is (# of samples / 4 ) x (# of markers).  (# of samples is rounded to a multiple of 4.
#'  Each byte stores data for 4 samples; a byte has 4 - 2 bit encodings.
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @importFrom stats aggregate
#'
#' @return None
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' mkGenABELgenotype(envir)
#'}
mkGenABELgenotype = function(markers = NULL, envir) {
# browser("convert")
    if (is.null(markers)) markers = envir$markers

## print (system.time ({
    raw_mtx = getgenotypesgenabel(markers, envir = envir)
## }))
    return (raw_mtx)
}
