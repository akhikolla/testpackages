#   Mega2R: Mega2 for R.
#
#   Copyright 2018-2019, University of Pittsburgh. All Rights Reserved.
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

## For installation do:
## source("http://bioconductor.org/biocLite.R")
## biocLite("gdsfmt")


#library(Mega2R)
#library(gdsfmt)

append_genotype_a=T
append_genotype_b=F
append_genotype_c=F


#' transcode mega2 to gdsfmt/SNP_ARRAY
#'
#' @description
#'   Reads the data frames in "envir" and builds a GDSFMT COREARRAY file from them.
#'
#' @param filename gdsfmt file to create
#'
#' @param markers data frame of markers to be processed
#'
#' @param snp.order TRUE indicates that the "genotype" data matrix has SNP as the first index
#'   which changes more quickly than subsequent indices.  FALSE indicates that SAMPLE is the
#'   the first index.
#'
#' @param SeqArray TRUE uses SeqArray labels for the gdsfmt vector elements.  FALSE it uses labels
#'  shown in SNPRelate
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return writes the "filename" file containing the CoreArray data.  Then returns an internal
#'   pointer, class .gds, to the file data.
#'
#' @export
#'
#' @importFrom gdsfmt createfn.gds openfn.gds closefn.gds cleanup.gds
#' @importFrom gdsfmt add.gdsn addfolder.gdsn compression.gdsn
#' @importFrom gdsfmt index.gdsn read.gdsn write.gdsn append.gdsn
#' @importFrom gdsfmt put.attr.gdsn get.attr.gdsn 
#'
#' @seealso \code{\link{gdsfmt}}
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db)
#' gdsfmtfile = file.path(where_mega2rtutorial_data(), "test.gds")
#' append_genotype_a = TRUE
#' append_genotype_b = append_genotype_c = FALSE
#' gn = Mega2gdsfmt(gdsfmtfile, envir=ENV)
#' gn
#'
Mega2gdsfmt = function(filename = "test.gds", markers = NULL, snp.order = FALSE,  SeqArray = FALSE, envir = ENV)
{
    unlink(filename)
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)

    if (SeqArray)
        Mega2gdsfmtSeq(filename, markers, snp.order, envir)
    else
        Mega2gdsfmtSNP(filename, markers, snp.order, envir)      
}

Mega2gdsfmtSeq = function(filename, markers, snp.order, envir) {

    COMPRESSION = "LZMA_RA"

# print(system.time ({

    g = createfn.gds(filename)
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)
    if (is.null(markers)) markers = envir$markers

    put.attr.gdsn(g$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(g$root, "FileVersion", "v1.0")
    envir$SEQ_ARRAY = TRUE

    desf = addfolder.gdsn(g, "description")

    gsam = add.gdsn(g, "sample.id", envir$fam$PerPre, compress=COMPRESSION)
    nsam = nrow(envir$fam)

    gsnp = add.gdsn(g, "variant.id", as.integer(markers$locus_link - envir$PhenoCnt + 1), compress=COMPRESSION)
    nsnp = nrow(markers)

    gpos = add.gdsn(g, "position", as.integer(markers$position), compress=COMPRESSION)

    gchr = add.gdsn(g, "chromosome", markers$chromosome, compress=COMPRESSION)
    put.attr.gdsn(gchr, "autosome.start", 1)
    put.attr.gdsn(gchr, "autosome.end", 22)
    put.attr.gdsn(gchr, "X",   23)
    put.attr.gdsn(gchr, "XY",  25)
    put.attr.gdsn(gchr, "Y",   24)
    put.attr.gdsn(gchr, "M",   26)
    put.attr.gdsn(gchr, "MT",  26)

#   alleles = mkAlleles(markers, envir = envir)
    alleles = getAlleles(markers, envir = envir)
    gall = add.gdsn(g, "allele", alleles, compress=COMPRESSION)

####
    genf = addfolder.gdsn(g, "genotype")

    if (snp.order) {
        vdim = c(2, nsnp, nsam)
        snpFirst(genf, markers, "data", vdim, envir)
    } else {
        vdim = c(2, nsam, nsnp)
        sampleFirst(genf, markers, "data", vdim, envir)
    }

    geif = add.gdsn(genf, "extra.index", valdim=c(3,0), compress=COMPRESSION, storage="int32")
    gext = add.gdsn(genf, "extra", valdim=c(0), compress=COMPRESSION, storage="int16")

####
    genf = addfolder.gdsn(g, "phase")
    if (snp.order)
        vdim = c(nsnp, nsam)
    else
        vdim = c(nsam, nsnp)
    geno = add.gdsn(genf, "data", rep(0, nsam * nsnp), valdim = vdim,
                    storage="bit1", compress=COMPRESSION)
    peif = add.gdsn(genf, "extra.index", valdim=c(3,0), compress=COMPRESSION, storage="int32")
    pext = add.gdsn(genf, "extra", valdim=c(0), compress=COMPRESSION, storage="bit1")
####
    
    annf  = addfolder.gdsn(g, "annotation")
    ann_s = add.gdsn(annf, "id", markers$MarkerName, compress=COMPRESSION)
    ann_q = add.gdsn(annf, "qual", rep(100.0, nsnp), compress=COMPRESSION, storage="float32")
    filter = factor(c(rep("PASS",nsnp)), levels=c("NA", "PASS"), exclude="")
    ann_ff = add.gdsn(annf, "filter", filter, compress=COMPRESSION, storage="int32")
    inf   = addfolder.gdsn(annf, "info")
    form  = addfolder.gdsn(annf, "format")
        
    samf  = addfolder.gdsn(g, "sample.annotation")
    sam_f = add.gdsn(samf, "family", envir$fam$PedPre, compress=COMPRESSION)

    closefn.gds(g)

    if (append_genotype_a || snp.order) {
        g = openfn.gds(filename, readonly=F)
        ggen = index.gdsn(g, "genotype/data")
        compression.gdsn(ggen, compress=COMPRESSION)
        closefn.gds(g)

        cleanup.gds(filename)
    }
    
# }))

    openfn.gds(filename)
}

Mega2gdsfmtSNP = function(filename, markers, snp.order, envir) {

    COMPRESSION = "LZMA_RA"

# print(system.time ({

    g = createfn.gds(filename)
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)
    if (is.null(markers)) markers = envir$markers

    put.attr.gdsn(g$root, "FileFormat", "SNP_ARRAY")
    put.attr.gdsn(g$root, "FileVersion", "v1.0")
    envir$SEQ_ARRAY = FALSE

    desf = addfolder.gdsn(g, "description")

    gsam = add.gdsn(g, "sample.id", envir$fam$PerPre, compress=COMPRESSION)
    nsam = nrow(envir$fam)

    gsnp = add.gdsn(g, "snp.id", as.integer(markers$locus_link - envir$PhenoCnt + 1), compress=COMPRESSION)
    nsnp = nrow(markers)

    grs = add.gdsn(g, "snp.rs.id", markers$MarkerName, compress=COMPRESSION)

    gpos = add.gdsn(g, "snp.position", as.integer(markers$position), compress=COMPRESSION)

    gchr = add.gdsn(g, "snp.chromosome", markers$chromosome, compress=COMPRESSION)
    put.attr.gdsn(gchr, "autosome.start", 1)
    put.attr.gdsn(gchr, "autosome.end", 22)
    put.attr.gdsn(gchr, "X",   23)
    put.attr.gdsn(gchr, "XY",  25)
    put.attr.gdsn(gchr, "Y",   24)
    put.attr.gdsn(gchr, "M",   26)
    put.attr.gdsn(gchr, "MT",  26)

#   alleles = mkAlleles(markers, envir = envir)
    alleles = getAlleles(markers, envir = envir)
    gall = add.gdsn(g, "snp.allele", alleles, compress=COMPRESSION)

    genf = g

    if (snp.order) {
        vdim = c(nsnp, nsam)
        snpFirst(genf, markers, "genotype", vdim, envir)
    } else {
        vdim = c(nsam, nsnp)
        sampleFirst(genf, markers, "genotype", vdim, envir)
    }
    
    df = data.frame(sample.id = envir$fam$PerPre,
                    family.id = envir$fam$PedPre,
                    father.id = envir$fam$Father,
                    mother.id = envir$fam$Mother,
                    sex.id    = c("male", "female")[envir$fam$Sex],
                    cc.id     = envir$fam$trait,
    #   pop.id = c("CEU", "HCB", "JPT", "YRI")[cc.d + 1]
                    stringsAsFactors = FALSE)

    gann = add.gdsn(g, "sample.annot", val=df, compress=COMPRESSION)

    closefn.gds(g)

    if (append_genotype_a || snp.order) {
        g = openfn.gds(filename, readonly=F)
        ggen = index.gdsn(g, "genotype")
        compression.gdsn(ggen, compress=COMPRESSION)
        closefn.gds(g)
    
        cleanup.gds(filename)
    }

# }))
    
    openfn.gds(filename)
}


## You seem to fill columns at a time, but the column can be samples (of snps) [snpOrder]
##  OR snps (of samples) [sampleFirst]

## Bit2
##  0 = B/B
##  1 = A/B
##  2 = A/A
##  3 = 0/0

snpBlock = function(markers, r, envir) {

    cr = getgenotypesraw(markers[envir$R, ], envir = envir )

    bit2 = envir$mt2[match(cr, envir$mt1)]
    dim(bit2) = dim(cr)

    st = envir$R[1]
    ln = length(envir$R)

    write.gdsn(envir$geno, t(bit2), start=c(st, 1), count=c(ln, -1))
}

seqBlock = function(markers, r, envir) {

    cr = getgenotypesraw(markers[envir$R, ], envir = envir )
    cr = t(cr)

    st = envir$R[1]
    ln = length(envir$R)

    a1a2 = match(cr, envir$mt1)

    x1 = envir$mt3[a1a2]
    x2 = envir$mt4[a1a2]
    x12   = array(c(x1, x2), dim = c(dim(cr), 2))
    bit2  = aperm(x12, perm=c(3, 1, 2))

    write.gdsn(envir$geno, bit2, start=c(1, st, 1), count=c(-1, ln, -1))
}

snpFirst = function(gds, markers, genotype, vdim, envir) {

    geno = add.gdsn(gds, genotype,  valdim = vdim, storage="bit2")
    envir$geno = geno
    envir$mt1 = c(0x10001, 0x10002, 0x20001, 0x20002, 0)
    envir$mt2 = c(      2,       1,       1,       0, 3)

    put.attr.gdsn(geno, "snp.order")

    M = nrow(markers)
    C = 1000

    j = 0
    while (TRUE) {
        if (M <= 0) break
        R = ((j*C+1):(j*C + C))
        L = length(R)
        if (M < L) {
            L = M
            R = ((j*C+1):(j*C+L))
        }
        envir$R  = R
        M = M - L

        if (envir$SEQ_ARRAY)
            tryFn(seqBlock, markers, NULL, envir)
        else
            tryFn(snpBlock, markers, NULL, envir)

        j = j + 1
    }
}

## ALL
# 12.857   0.767  13.711 
# 12.896   0.766  13.863 
## APPEND_GENOTYPE_A
#  5.640   0.520   6.203 
#  5.601   0.527   6.179 
## APPEND_GENOTYPE_B
#  6.609   0.537   7.161 
#  6.624   0.542   7.181 
## APPEND_GENOTYPE_C  
#  6.884   0.507   7.420 
#  6.899   0.508   7.427 
## X0
#  3.213   0.389   3.608 
#  3.469   0.418   3.892 

seqsampleBlock = function(markers, r, envir) {

    cr = getgenotypesraw(markers[envir$R, ], envir = envir )
    st = envir$R[1]
    ln = length(envir$R)

    a1a2 = match(cr, envir$mt1)

    x1 = envir$mt3[a1a2]
    x2 = envir$mt4[a1a2]
    x12   = array(c(x1, x2), dim = c(dim(cr), 2))
    bit2  = aperm(x12, perm=c(3, 1, 2))

    if (append_genotype_a) {
        write.gdsn(envir$geno, bit2, start=c(1,  1, st), count=c(-1, -1, ln))
    }
    if (append_genotype_b) {
        append.gdsn(envir$genx, bit2[,,1:ln])    ## 10.330  0.690  11.084 + geno
    }
    if (append_genotype_c) {
        for (i in 1:ln) {
            bt2 = bit2[,,i]
            append.gdsn(envir$geny, bt2)         ## 12.162  0.667  12.906 + geno
        }
    }
}

snpsampleBlock = function(markers, r, envir) {

    cr = getgenotypesraw(markers[envir$R, ], envir = envir )
    st = envir$R[1]
    ln = length(envir$R)

    a1a2 = match(cr, envir$mt1)

    bit2 = envir$mt2[a1a2]
    dim(bit2) = dim(cr)

    if (append_genotype_a) {
        write.gdsn(envir$geno, bit2, start=c(1, st), count=c(-1, ln))
    }
    if (append_genotype_b) {
        append.gdsn(envir$genx, bit2[,1:ln])    ## 10.330  0.690  11.084 + geno
    }
    if (append_genotype_c) {
        for (i in 1:ln) {
            bt2 = bit2[,i]
             append.gdsn(envir$geny, bt2)         ## 12.162  0.667  12.906 + geno
        }
    }
}

sampleFirst = function(gds, markers, genotype, vdim, envir) {

    if (append_genotype_a) {
        geno = add.gdsn(gds, genotype, valdim = vdim, storage="bit2")
        envir$geno = geno
        put.attr.gdsn(geno, "sample.order")
    }

    if (append_genotype_b) {
        vdim[length(vdim)] = 0
        genx = add.gdsn(gds, genotype, valdim = vdim, storage="bit2", compress="LZMA_RA")
        envir$genx = genx
    }

    if (append_genotype_c) {
        vdim[length(vdim)] = 0
        geny = add.gdsn(gds, genotype, valdim = vdim, storage="bit2", compress="LZMA_RA")
        envir$geny = geny
    }
    #               65537    65538   131073   131074
    #                  ok                         ok
    envir$mt1 = c(0x10001, 0x10002, 0x20001, 0x20002, 0)
    envir$mt2 = c(      2,       1,       1,       0, 3)
    envir$mt3 = c(      0,       0,       1,       1, 3)
    envir$mt4 = c(      0,       1,       0,       1, 3)

    M = nrow(markers)
    C = 1000

    j = 0
    while (TRUE) {
        if (M <= 0) break
        R = ((j*C+1):(j*C + C))
        L = length(R)
        if (M < L) {
            L = M
            R = ((j*C+1):(j*C+L))
        }
        envir$R  = R
        M = M - L

        if (envir$SEQ_ARRAY)
            tryFn(seqsampleBlock, markers, NULL, envir)
        else
            tryFn(snpsampleBlock, markers, NULL, envir)

        j = j + 1
    }
}

#' regenerate fam data frame with unique values in member column
#'
#' @description
#'   Reads the fam data frame in "envir" and returns a new one with unique entries in the member
#'   column
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return a data frame with columns the same as the "fam" data frame but with the member column
#'   containing unique entries
#'
#' @export
#'
#' @seealso \code{\link{mkfam}}
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db)
#' setfam(uniqueFamMember(envir = ENV))
#'
uniqueFamMember = function(envir = ENV) {

    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)

    fam = envir$fam

    mp = fam$Father != 0
    fam$Father = ifelse(mp, paste0(fam$PedPre, "__" , fam$Father), fam$Father)

    mp = fam$Mother != 0
    fam$Mother = ifelse(mp, paste0(fam$PedPre, "__" , fam$Mother), fam$Mother)

    fam$PerPre = paste0(fam$PedPre , "__" , fam$PerPre)

    fam
}

#' generate allele pairs in with MAJ(or) allele first
#'
#' @description
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
#' mkAlleles(envir)
#'}
mkAlleles = function(markers = NULL, separator = "/", envir = ENV) {

    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)
    if (is.null(markers)) markers = envir$markers

    allele_table = envir$allele_table[envir$allele_table$locus_link %in% markers$locus_link,]
    mm = merge(x=allele_table[allele_table$indexX == 1,],
               y=allele_table[allele_table$indexX == 2,],
               by="locus_link")

    nn=ifelse(mm$Frequency.x > mm$Frequency.y,             # before cleaning
              paste0(mm$AlleleName.x, separator, mm$AlleleName.y),
              paste0(mm$AlleleName.y, separator, mm$AlleleName.x))
    envir$xGTy = mm$Frequency.x > mm$Frequency.y

    if (any(mm$Frequency.x == .5)) {
        w = which(mm$Frequency.x == .5)
##      if (envir$MARKER_SCHEME == 1) {
##          ms = envir$markerscheme_table[envir$markerscheme_table$key %in% markers$locus_link,]
##          print(mm[w,])
##          print(ms[w,])
##          print(envir$markers[w,])
##          ms1 = ms$allele1[w]
##          nn[w] = ifelse(ms1 == 1, paste0(mm$AlleleName.y[w], separator, mm$AlleleName.x[w]),
##                                   paste0(mm$AlleleName.x[w], separator, mm$AlleleName.y[w]))
##      } else if (envir$MARKER_SCHEME == 2) {
##          print(mm[w,])
##      }
        print(mm[w,])
        print(nn[w])
    }

    nn[mm$Frequency.x == 0 & mm$Frequency.y == 0] = paste0('0', separator, '0')

    fx = mm$Frequency.x == 1
    nn[fx] = paste0(mm[fx, "AlleleName.x"], separator, mm[fx, "AlleleName.x"])

    fy = mm$Frequency.x == 0
    nn[fy] = paste0(mm[fy, "AlleleName.y"], separator, mm[fy, "AlleleName.y"])

    nn[nn=="00"] = paste0('0', separator, '0')

    nn
}

#' pull allele pairs
#'
#' @description
#'  get alleles in order
#'
#' @param markers data frame of markers to be processed
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return 
#'   string vector of alleles with separator
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' getAlleles(envir)
#'}
getAlleles = function(markers = NULL, separator = "/", envir = ENV) {

    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)
    if (is.null(markers)) markers = envir$markers

    allele_table = envir$allele_table[envir$allele_table$locus_link %in% markers$locus_link,]

    all = vapply(split(allele_table, allele_table$locus_link), FUN.VALUE="A",
                 function(x) { paste0(x$AlleleName, collapse=separator) } )
    names(all) = NULL
    all
#   paste(allele_table[allele_table$indexX == 1,]$AlleleName,
#         allele_table[allele_table$indexX == 2,]$AlleleName,
#         sep = separator)
}
