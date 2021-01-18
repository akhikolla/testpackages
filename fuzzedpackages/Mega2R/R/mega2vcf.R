
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


#' generate a VCF file set for a collection of markers
#'
#' @description
#'  Generate a VCF file from the specified Mega2 SQLite database.  The file is named \emph{"prefix".vcf}
#'  If the markers argument is.null(), the entire \bold{envir$markers} set is used, otherwise markers argument MUST
#'  be rows of the markers (\bold{envir$markers}) data frame.
#'  In addition,
#'  several other files are generated to hold additional database information: \emph{"prefix".fam}, \emph{"prefix".freq}, \emph{"prefix".map},
#'  \emph{"prefix".phe}, and \emph{"prefix".pen}, which contain the pedigree, allele frequency, marker genetic and
#'  physical map position, member phenotype and phenotype penetrance data.
#'
#' @param prefix prefix of output files including the VCF file (see Description section above). This prefix can include a path.
#'
#' @param markers markers selected to be in the VCF output file
#'
#' @param mapno specify which map index to use for genetic distances.  The function \code{showMapNames()}
#' will print out the internal map numbers corresponding to all the maps in the Mega2 database.
#'
#' @param alleleOrder how to order alleles in VCF file.
#' 'default' is Mega2order, 'minor' is minor allele freq first, 'major' is major allele freq
#'  first, and 'name' is ascending ascii character order of allele name.
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return None
#'
#' @importFrom utils write.table
#' @export
#'
#' @note This code in this package illustrates how to extract the various kinds of data in the
#'  Mega2 data frames to use for further processing.  Some of the data internal representations
#'  are a bit quirky but the code "explains" it all.
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db)
#' vcfdir = file.path(where_mega2rtutorial_data(), "vcfr")
#' if (!dir.exists(vcfdir)) dir.create(vcfdir)
#' vcffile = file.path(where_mega2rtutorial_data(), "vcfr", "vcf.01")
#' Mega2VCF(vcffile, ENV$markers[ENV$markers$chromosome == 1, ][1:10,], envir = ENV)
#' list.files(vcfdir)
#'
Mega2VCF = function(prefix, markers=NULL, mapno = 0, alleleOrder = 'default', envir = ENV) {
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)

    if (missing(prefix))
        stop("Mega2VCF can not proceed without a filename prefix argument", call. = FALSE)

    file = paste0(prefix, ".vcf")
    unlink(file)

    if (is.null(markers)) markers = envir$markers

    match49r = FALSE
    
    mkVCFhdr(prefix, markers, envir)

    z = c("./.", "./0", "./1", "./2", "0/.", "0/0", "0/1", "0/2",
          "1/.", "1/0", "1/1", "1/2", "2/.", "2/0", "2/1", "2/2")
    zs = sort(z)
    zz = matrix(z, nrow = 4, byrow = TRUE)

    allele_table = envir$allele_table[envir$allele_table$locus_link %in% markers$locus_link,]
    map_table = envir$map_table[envir$map_table$marker %in% markers$locus_link,]
    M = nrow(markers)
    C = 1000

    block = data.frame(matrix(0, nrow = C, ncol = nrow(envir$fam) + 9), stringsAsFactors=TRUE)
    blockcol = ncol(block)

    QUAL   = rep(".",    times=C)
    FILTER = rep("PASS", times=C)
    FORMAT = rep("GT",   times=C)
    block[ , 6] = QUAL
    block[ , 7] = FILTER
    block[ , 9] = FORMAT

    names(block) = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                         paste0(envir$fam$PedPre, "_", envir$fam$PerPre))
    cat(names(block), sep="\t", append=TRUE, file=file)
    cat("\n", append=TRUE, file=file)

    allele_table[allele_table$AlleleName == "dummy", 2] = '.'
    atf = allele_table$AlleleName == ""
    allele_table[atf, 2] = as.character(allele_table[atf, 4])

    n1 = allele_table[allele_table$indexX==1,]
    n2 = allele_table[allele_table$indexX==2,]

    if (envir$MARKER_SCHEME == 2) {
        if (alleleOrder == 'default')
            xn = function (x) x[order(x$indexX),]
        else if (alleleOrder == 'minor')
            xn = function(x) { a=order(x$Frequency); x[c(a[1], sort(a[2:length(a)])),] }
        else if (alleleOrder == 'major')
            xn = function(x) { a=order(x$Frequency, decreasing=TRUE); x[c(a[1], sort(a[2:length(a)])),] }
        else if (alleleOrder == 'name')
            xn = function (x) x[order(x$AlleleName),]
        else
            stop("alleleOrder must be 'default', 'minor' or 'major'\n", call.=FALSE)
        nx = lapply(split(allele_table, allele_table$locus_link), xn)
    }

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

        block[BR , 1] = markers[R , 4]
        block[BR , 2] = markers[R , 5]
        block[BR , 3] = markers[R , 3]
## print(system.time ({

      if (envir$MARKER_SCHEME == 1) {
            REF = n1$AlleleName[R]
            RF  = n1$Frequency[R]
            ALT = n2$AlleleName[R]
            AF  = n2$Frequency[R]
            XXF = REF
            XF  = RF

            if (alleleOrder == 'major' && any(RF < AF)) {
                whichFlip = which(RF < AF)
                doFlip = TRUE
            } else if (alleleOrder == 'minor' && any(AF < RF)) {
                whichFlip = which(AF < RF)
                doFlip = TRUE
            } else if (match49r && any (ALT < REF) ) {
                whichFlip = which(ALT < REF)
                whichFlip = whichFlip[ALT[whichFlip] != "0"]
                whichFlip = whichFlip[AF[whichFlip] != 0]
                doFlip = TRUE
            } else
                doFlip = FALSE
            if (doFlip) {
                REF[whichFlip] = ALT[whichFlip]
                ALT[whichFlip] = XXF[whichFlip]

                RF[whichFlip] = AF[whichFlip]
                AF[whichFlip] = XF[whichFlip]
            }
#           AF = sprintf("%.6f", AF)
            AF = sprintf("%g", round(AF, 6))
        } else if (envir$MARKER_SCHEME == 2) {
            REF = character(L)
            ALT = character(L)
            RF  = numeric(L)
            AF  = character(L)

            k = 0
            for (v in R) {
                k = k + 1
                nxx = nx[[v]]
                nrows = nrow(nxx)
                REF[k] = nxx$AlleleName[1]
                RF[k]  = nxx$Frequency[1]
                ALT[k] = paste0(nxx$AlleleName[2:nrows], collapse=",")
#               AF[k]  = paste0(sprintf("%.6f", nxx$Frequency[2:nrows]), collapse=",")
                AF[k]  = paste0(sprintf("%g", round(nxx$Frequency[2:nrows], 6)), collapse=",")
            }
        }
        GPos = map_table[map_table$map==mapno, c("position", "pos_female", "pos_male")][R, ]
        GPosPos = sprintf("%g", round(GPos$position, 2))
        GPosFem = rep(".", L)
        GPos[is.na(GPos$pos_female),2] = -99.99
        GPosFem[GPos$pos_female != -99.99] = sprintf("%g", round(GPos$pos_female[GPos$pos_female != -99.99], 2))
        GPosMal = rep(".", L)
        GPos[is.na(GPos$pos_male),3] = -99.99
        GPosMal[GPos$pos_male   != -99.99] = sprintf("%g", round(GPos$pos_male[GPos$pos_male != -99.99], 2))

        INFO=paste0("CM=", GPosPos, ",", GPosFem, ",", GPosMal,
#                   ";RF=", sprintf("%.6f", RF),
                    ";RF=", sprintf("%g", round(RF, 6)),
                    ";AF=", AF)
#                   ";")
        block[BR , 4] = REF[BR]
        block[BR , 5] = ALT[BR]
#       block[BR , 6] = QUAL[BR]
#       block[BR , 7] = FILTER[BR]
        block[BR , 8] = INFO
#       block[BR , 9] = FORMAT[BR]

        cr = getgenotypesraw(markers[R, ], envir = envir)   # 7.17%
        a1 = t(cr)                                          # 0.86%
        a2 = a1
        dm = dim(a1)
        a1 = bitwShiftR(a1, 16)
        attr(a1, "dim") = dm
        x3 = a1
        a2 = bitwAnd(a2, 65535)
        attr(a2, "dim") = dm
        if (envir$MARKER_SCHEME == 1) {       # 41.32%
            if (doFlip) {
                a1[whichFlip, ] = match(a1[whichFlip, ], c(2, 1), nomatch=0)
                a2[whichFlip, ] = match(a2[whichFlip, ], c(2, 1), nomatch=0)
            }
        }

#       if ((envir$MARKER_SCHEME > 1) && (a1 > 2 || a2 > 2))        # 41.32%
        if (envir$MARKER_SCHEME > 1) {       # 41.32%
                ##  user  system elapsed
                ##  4.086   0.133   4.251
                ##  user  system elapsed
                ##  1.275   0.053   1.351
            k = 0
            for (v in R) {
                k = k + 1
                nxx = nx[[as.character(markers[v,]$locus_link)]]$indexX
                a1[k,] = match(a1[k,], nxx, nomatch=0)
                a2[k,] = match(a2[k,], nxx, nomatch=0)
            }
            a3 = as.character(a1 - 1)
            a3[a1 == 0] = "."
            a4 = as.character(a2-1)
            a4[a2 == 0] = "."
            a5 = paste0(a3, "/", a4)
            block[BR, 10:blockcol] = a5
        } else
            ##  user  system elapsed
            ##  .774   0.092   0.879
            ##  user  system elapsed
            ##  1.274   0.056   1.351
            block[BR, 10:blockcol] = zz[cbind(as.vector(a1)+1, as.vector(a2)+1)]

##      }))

## print(system.time ({
        write.table(block[BR, ], file=file, sep="\t", quote=FALSE,     # 48.49%
                    append=TRUE, row.names=FALSE, col.names=FALSE)
## }))
        j = j + 1
##      if (envir$verbose) message(".", appendLF = FALSE)
    }
}

#' generate required VCF header
#'
#' @description
#'  Generate the initial boiler plate VCF, then generate ##INFO entries for each entry tag.
#'  Finally, generate the ##contig entries for each chromosome.
#'
#' @param prefix prefix for .vcf file
#'
#' @param markers data.frame of markers being processed
#'
#' @param envir "environment" containing SQLite database and other globals
#'
#' @return None
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' mkVCFhdr(prefix, NULL, envir)
#'}
mkVCFhdr = function (prefix, markers, envir) {
    file = paste0(prefix, ".vcf")

    if (is.null(markers)) markers = envir$markers

    mkVCFfam(prefix,           envir)
    mkVCFfreq(prefix, markers, envir = envir)
    mkVCFmap(prefix,  markers, envir)
    mkVCFpen(prefix,           envir)
    mkVCFphe(prefix,           envir)

    cat('##fileformat=VCFv4.1\n', file=file, append=TRUE)
##  cat('##filedate=19970829\n', file=file, append=TRUE)
    cat('##filedate=', format(Sys.time(), "%Y%m%d"), '\n', sep = "", file=file, append=TRUE)
    cat('##source=MEGA2\n', file=file, append=TRUE)
    cat('##INFO=<ID=CM,Number=3,Type=Float,Description="Genetic Distance in centimorgans (avg, male, female)">\n', file=file, append=TRUE)
    cat('##INFO=<ID=RF,Number=1,Type=Float,Description="Allele Frequency of reference allele">\n', file=file, append=TRUE)
    cat('##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency of alternate allele(s)">\n', file=file, append=TRUE)
    cat('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n', file=file, append=TRUE)
    cat('##FILTER=<ID=PASS,Description="Passed variant FILTERs">\n', file=file, append=TRUE)

    maxes = sapply(split(markers, markers$chromosome), function(x) max(x$position)+1)
    for (i in 1:length(maxes)) {
        cat('##contig=<ID=', names(maxes)[i], ',length=', maxes[i], ',assembly=B37>\n',
            file=file, append=TRUE, sep="")
    }
}

#' generate required VCF family (.fam) file
#'
#' @description
#'  Generate the initial boiler plate VCF, then generate ##INFO entries for each entry tag.
#'  Finally, generate the ##contig entries for each chromosome.
#'
#' @param prefix prefix for .fam file (family pedigree)
#'
#' @param envir "environment" containing SQLite database and other globals
#'
#' @return None
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' mkVCFfam(prefix, envir)
#'}
mkVCFfam = function (prefix, envir) {
    file = paste0(prefix, ".fam")

#   -9 vs 0 for case/control

# mega2 mis
    envir$fam[envir$fam[ , 8] ==0, 8] = -9
    write.table(envir$fam[, -(1:2)], file=file, sep="\t", quote=FALSE,
                    row.names=FALSE, col.names=FALSE)
}

#' generate required VCF frequency (.freq) file
#'
#' @description
#'  Generate the initial boiler plate VCF, then generate ##INFO entries for each entry tag.
#'  Finally, generate the ##contig entries for each chromosome.
#'
#' @param prefix prefix for .freq file (frequency)
#'
#' @param markers data.frame of markers being processed
#'
#' @param recode use 1/2 instead of two given alleles (eg. A/C)
#'
#' @param envir "environment" containing SQLite database and other globals
#'
#' @return None
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' mkVCFfreq(prefix, NULL, FALSE, envir)
#'}
mkVCFfreq = function (prefix, markers, recode = FALSE, envir) {
    file = paste0(prefix, ".freq")

#    unlink(file)

    if (is.null(markers)) markers = envir$markers

    if (recode)
        col = "indexX"
    else
        col = "AlleleName"
    cat('Name\tAllele\tFrequency\n', file=file)
    allele_pheno = merge(envir$locus_table, envir$allele_table[1:(2*envir$PhenoCnt),], by="locus_link")
    alleles = allele_pheno[, c("LocusName", col, "Frequency")]
    alleles[alleles[,col] == "", col] = allele_pheno[alleles[,col] == "", "indexX"]
#std
    alleles$Frequency = sprintf("%.6f", alleles$Frequency)
    write.table(alleles,
                file=file, sep="\t", append=TRUE, quote=FALSE,
                row.names=FALSE, col.names=FALSE)

    allele_table = envir$allele_table[envir$allele_table$locus_link %in% markers$locus_link,
                                    c("locus_link", "AlleleName", "indexX", "Frequency")]
    alleles = merge(markers[, c("locus_link", "MarkerName")],
                    allele_table[, c("locus_link", "AlleleName", "indexX", "Frequency")],
                    by="locus_link")
    alleles[alleles[ , col] == "", col] = alleles[alleles[ , col] == "", "indexX"]
##  alleles = alleles[alleles$Frequency != 0, ]
#std
    alleles$Freq4 = sprintf("%.6f", alleles$Frequency)
    write.table(alleles[ , c(-1, -4, -5)],
                file=file, sep="\t", append=TRUE, quote=FALSE,
                row.names=FALSE, col.names=FALSE)
}

#' generate required Mega2 map (.map) file
#'
#' @description
#'  Generate the initial boiler plate VCF, then generate ##INFO entries for each entry tag.
#'  Finally, generate the ##contig entries for each chromosome.
#'
#' @param prefix prefix for .map file
#'
#' @param markers data.frame of markers being processed
#'
#' @param envir "environment" containing SQLite database and other globals
#'
#' @return None
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' mkVCFmap(prefix, NULL, envir)
#'}
mkVCFmap = function (prefix, markers, envir) {
    file = paste0(prefix, ".map")

    unlink(file)

    if (is.null(markers)) markers = envir$markers

    map_table = envir$map_table[envir$map_table$marker %in% markers$locus_link,]
    mapnames_table = envir$mapnames_table
    TBL = markers[, c("chromosome", "MarkerName")]
    hdr = paste("Chromosome", "Name", sep="\t")

    for (m in mapnames_table$map) {
        if ( (mapnames_table[m+1, "male_sex_map"] == 0) &
             (mapnames_table[m+1, "female_sex_map"] == 0) ) {
            POS  = map_table[map_table$map == m, "position"]

            if (mapnames_table[m+1, "sex_averaged_map"] == 0) {
                TBL  = cbind(TBL, POS)
                hdr = paste0(hdr, "\t", mapnames_table[mapnames_table$map == m, "name"], '.p')
            } else {
#std
                POSS = sprintf("%.6f", POS)
                TBL  = cbind(TBL, POSS)
                hdr = paste0(hdr, "\t", mapnames_table[mapnames_table$map == m, "name"], '.k.a')
            }
          }

        if ( mapnames_table[m+1, "female_sex_map"] != 0) {
            POSF  = map_table[map_table$map == m, "pos_female"]
#std
            POSF[is.na(POSF$pos_female)] = -99.99
            POSFS = sprintf("%.6f", POSF)
            TBL  = cbind(TBL, POSFS)
            hdr = paste0(hdr, "\t", mapnames_table[mapnames_table$map == m, "name"], '.k.f')
        }

        if ( mapnames_table[m+1, "male_sex_map"] != 0) {
            POSM  = map_table[map_table$map == m, "pos_male"]
#std
            POSM[is.na(POSM$pos_male)] = -99.99
            POSMS = sprintf("%.6f", POSM)
            TBL  = cbind(TBL, POSMS)
            hdr = paste0(hdr, "\t", mapnames_table[mapnames_table$map == m, "name"], '.k.m')
        }
    }

    cat(hdr, "\n", file=file, sep="")

    write.table(TBL, file=file, sep="\t", quote=FALSE, append=TRUE,
                row.names=FALSE, col.names=FALSE)
}

#' generate required Mega2 penetrance (.pen) file
#'
#' @description
#'  Generate the initial boiler plate VCF, then generate ##INFO entries for each entry tag.
#'  Finally, generate the ##contig entries for each chromosome.
#'
#' @param prefix prefix for .pen file (penetrance)
#'
#' @param envir "environment" containing SQLite database and other globals
#'
#' @return None
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' mkVCFpen(prefix, envir)
#'}
mkVCFpen = function (prefix, envir) {
    file = paste0(prefix, ".pen")

    unlink(file)

    cat('Name\tClass\tPen.11\tPen.12\tPen.22\tType\n', file=file, append=TRUE)

    all = merge(merge(envir$locus_table[1:envir$PhenoCnt,], envir$traitaff_table, by="locus_link"),
                envir$affectclass_table, by="locus_link")
    ord = order(all$locus_link, all$class_link)
    for (i in ord) {
#       malepen   = all[i, ]$MalePen[[1]]
#       mpen  = readBin(malepen, numeric(), n = length(malepen)/8, size = 8, endian = .Platform$endian)
#       mpens = sprintf("%.4f", mpen)
#       cat(all[i, ]$LocusName, all[i, ]$class_link+1, mpens, file=file, append=TRUE, sep="\t")
#       cat("\t\tmale\n", file=file, append=TRUE)

#       femalepen = all[i, ]$FemalePen[[1]]
#       fpen  = readBin(femalepen, numeric(), n = length(femalepen)/8, size = 8, endian = .Platform$endian)
#       fpens = sprintf("%.4f", fpen)
#       cat(all[i, ]$LocusName, all[i, ]$class_link+1, fpens, file=file, append=TRUE, sep="\t")
#       cat("\tfemale\n", file=file, append=TRUE)

# Each penetrance (for female/male/autosome) is 2 or 3 entries.  (This is for a bialleleic system.)
#  An entry is 8 bytes for a double.
        autopen   = all[i, ]$AutoPen[[1]]
        apen  = readBin(autopen, numeric(), n = length(autopen)/8, size = 8, endian = .Platform$endian)
        apens = sprintf("%.4f", apen)
        cat(all[i, ]$LocusName, all[i, ]$class_link+1, apens, file=file, append=TRUE, sep="\t")
        cat("\tautosomal\n", file=file, append=TRUE)
    }
}

#' generate required PLINK (.phe) file
#'
#' @description
#'  Generate the initial boiler plate VCF, then generate ##INFO entries for each entry tag.
#'  Finally, generate the ##contig entries for each chromosome.
#'
#' @param prefix prefix for .phe file (phenotype)
#'
#' @param envir "environment" containing SQLite database and other globals
#'
#' @return None
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' mkVCFphe(prefix, envir)
#'}
mkVCFphe = function (prefix, envir) {
    file = paste0(prefix, ".phe")

    unlink(file)

# linkage.h:    TYPE_UNSET, QUANT, AFFECTION, BINARY, NUMBERED, XLINKED, YLINKED
#                        0      1          2       3         4        5        6

    out = mkphenotype(envir)

    out$SAMPLEID = paste(out[,1], out[,2], sep="_")

    write.table(out, file=file, sep="\t", quote=FALSE, append=FALSE,
                row.names=FALSE, col.names=TRUE)
}
