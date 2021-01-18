
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

#library(pedgene)

#' load Mega2 SQLite database and perform initialization for pedgene usage
#'
#' @description
#'  This populates the \bold{R} data frames from the specified \bold{Mega2} SQLite database.
#'
#' @param db specifies the path of a \bold{Mega2} SQLite database containing study data.
#'
#' @param verbose TRUE indicates that diagnostic printouts should be enabled.
#'  This value is saved in the returned environment.
#'
#' @param traitname Name of the affection status trait to use to set the case/control status; default value = "default".
#'
#' @param ... fed to \emph{dbmega2_import()}; should be bpPosMap= to select from the maps of
#'  base pairs, if the default is not desired.
#'
#' @return "environment" containing data frames from an SQLite database and some computed values.
#'
#' @importFrom utils read.table write.table
#' @export
#'
#' @note
#'  \emph{init_pedgene} calculates schaidPed and pedPer that are used later in the \emph{Dopedgene} calculation.
#'  In addition, it initializes a matrix to aid
#'   in translating a genotype allele matrix to a genotype count matrix.
#'
#'  It also initializes the dataframe \emph{envir$pedgene_results} to zero rows.
#'
#' @seealso \code{\link{DOpedgene}}, \code{\link{Mega2pedgene}}, \code{\link{mkfam}}
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = init_pedgene(db, traitname = "default")
#' ls(ENV)
#'
init_pedgene = function (db = NULL, verbose = FALSE, traitname = "default", ...) {

    if (is.null(db))
        stop("You must specify a database argument!\n", call. = FALSE)

    envir = dbmega2_import(db, verbose = verbose, ...)

    fam = mkfam(envir = envir, traitname = traitname)
#   fam = fam[fam$trait != 0, ]  # b,c vs a
#   vs
    fam$trantrait = NA
    fam$trantrait[fam$trait == 1] = 0
    fam$trantrait[fam$trait == 2] = 1
    fam$trait = NULL

    setfam(fam, envir = envir)  # also updates unified_genotype_table

    envir$schaidPed = envir$fam[ , c(-1, -2)]
    colnames(envir$schaidPed) = c("famid", "person", "father", "mother", "sex", "trait")
    envir$pedPer = envir$schaidPed[ , 1:2]
#   envir$mt = matrix(c(11, 12, 21, 22, 0,    0, 1, 1, 2, 0), nrow = 5, ncol = 2)
    envir$mt1 = c(0x10001, 0x10002, 0x20001, 0x20002, 0)
    envir$mt2 = c(      0,       1,       1,       2, 0)
    envir$pedgene_results <- data.frame(chr = character(0), gene = character(0),
                                        nvariants = numeric(0),
                                        start = numeric(0), end = numeric(0),
                                        sKernel_BT = numeric(0), pKernel_BT = numeric(0),
                                        sBurden_BT = numeric(0), pBurden_BT = numeric(0),
#                                       call_BT    = character(0),

                                        sKernel_MB = numeric(0), pKernel_MB = numeric(0),
                                        sBurden_MB = numeric(0), pBurden_MB = numeric(0),
#                                       call_MB    = character(0),

                                        sKernel_UW = numeric(0), pKernel_UW = numeric(0),
                                        sBurden_UW = numeric(0), pBurden_UW = numeric(0),
#                                       call_UW    = character(0),

                                        stringsAsFactors = FALSE)
    return (envir)
}


#' Execute the pedgene function on a transcript ranges
#'
#' @description
#' Execute the pedgene function on the first \emph{gs} default gene transcript ranges (gs = 1:100).
#'  Update the \emph{envir$pedgene_results} data frame with the results.
#"
#' @param gs a subrange of the default transcript ranges over which to calculate the \emph{Dopedgene} function.
#'
#' @param genes a list of genes over which to calculate the \emph{DOpedgene} function.
#'  The value, "*", means use all the transcripts in the selected Bioconductor database.
#'  If genes is NULL, the gs range of the internal \emph{refRanges} will be used.
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return None
#'  the data frame with the results is stored in the environment and named \emph{pedgene_results},
#'  viz. envir$pedgene_results
#'
#' @export
#'
#' @seealso \code{\link{init_pedgene}}
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = init_pedgene(db)
#' ENV$verbose = TRUE
#' Mega2pedgene(gs = 50:60)
#'
Mega2pedgene = function (gs = 1:100, genes = NULL, envir = ENV) {
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)

    if (is.null(genes))
        applyFnToRanges(DOpedgene, envir$refRanges[gs, ], envir$refIndices, envir = envir)
    else
        applyFnToGenes(DOpedgene, genes, envir = envir)
}

#' pedgene call back function
#'
#' @description
#'  First, ignore call backs that have less than two markers.  Second, convert the genotypesraw()
#'  patterns of 0x10001, 0x10002 (or 0x20001), 0x20002, 0 from the genotype matrix
#'  to the numbers 0, 1, 2, 0 for each marker. (Reverse, the order iff allele "1" has the
#'  minor allele frequency.)  Next, prepend the pedigree and person columns of the family data
#'  to this modified genotype matrix.  Finally, invoke \code{pedgene} with the family data and
#'  genotype matrix for several different weights.  Save the kernel and burden, value and p-value for each
#'  measurement in \emph{envir$pedgene_results}.
#'
#' @param markers_arg a data.frame with the following 5 observations:
#' \describe{
#' \item{locus_link}{is the ordinal ranking of this marker among all loci}
#' \item{locus_link_fill}{is the position of corresponding genotype data in the
#' \emph{unified_genotype_table}}
#' \item{MarkerName}{is the text name of the marker}
#' \item{chromosome}{is the integer chromosome number}
#' \item{position}{is the integer base pair position of marker}
#'  }
#'
#' @param range_arg one row of a ranges_arg.  The latter is a data frame of at least three
#'  integer columns.  The columns indicate a range:
#'  a chromosome number, a start base pair value, and an end base pair value.
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return None
#' @importFrom pedgene pedgene
#' @export
#'
#' @note
#'  This function appends output to the data frame, \emph{envir$pedgene_results}.  It will
#'  print out the lines as they are generated if \emph{envir$verbose} is TRUE. The data frame
#'  \emph{envir$pedgene_results} is initialized by \emph{init_pedgene}, and is appended to
#'  each time \emph{DOpedgene} is run.
#'
#' @seealso \code{\link{init_pedgene}}
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = init_pedgene(db)
#' ENV$verbose = TRUE
#' applyFnToRanges(DOpedgene, ENV$refRanges[50:60,], ENV$refIndices)
#'
#' \donttest{
#' # donttestcheck: try this below if there is time
#' applyFnToGenes(DOpedgene, genes_arg = c("CEP104"))
#' }
#'
DOpedgene = function(markers_arg, range_arg, envir = ENV) {

    if (is.null(range_arg))
        stop("DOpedgene: range is not defined.", calls. = FALSE)

    geno_arg = getgenotypesraw(markers_arg, envir);

    markerNames = markers_arg$MarkerName
    gene  = as.character(range_arg[,envir$refCol[4]])

    di = dim(geno_arg)
    geno = matrix(0, nrow = (di[1]), ncol = di[2])
    for (k in 1:(di[2])) {
        vec = envir$mt2[match(as.integer(geno_arg[ , k]), envir$mt1)]
        g0 = sum(vec == 0)
        g1 = sum(vec == 1)
        g2 = sum(vec == 2)
        if (envir$verbose)
            cat(gene, markerNames[k], g0, g1, g2, "\n")
        if (g0 < g2) {
           geno[ , k] = 2 - vec
        } else {
           geno[ , k] =     vec
        }
    }

    geno = matrix(geno, nrow = di[1])
    maf = colMeans(geno)
    pos = markerNames[maf > 0]
    if (length(pos) >= 2) {       # at least 2 non-polymorphic variants #
        geno <- geno[ , maf > 0]     # remove nonpolymorphic variants #
        nsnp    <- ncol(geno)
        weight <- rep(1, ncol(geno))

        pedgeno <- cbind(envir$pedPer, geno)

        BT <- pedgene(envir$schaidPed, pedgeno, male.dose= 2, checkpeds= FALSE, weights= NULL, weights.mb= FALSE, method= "kounen")
        sKernel_BT <- BT$pgdf$stat.kernel
        pKernel_BT <- BT$pgdf$pval.kernel
        sBurden_BT <- BT$pgdf$stat.burden
        pBurden_BT <- BT$pgdf$pval.burden
#       call_BT    <- BT$call

        MB <- pedgene(envir$schaidPed, pedgeno, male.dose= 2, checkpeds= FALSE, weights= NULL, weights.mb= TRUE, method= "kounen")
        sKernel_MB <- MB$pgdf$stat.kernel
        pKernel_MB <- MB$pgdf$pval.kernel
        sBurden_MB <- MB$pgdf$stat.burden
        pBurden_MB <- MB$pgdf$pval.burden
#       call_MB    <- MB$call

        UW <- pedgene(envir$schaidPed, pedgeno, male.dose= 2, checkpeds= FALSE, weights= weight, weights.mb= TRUE, method= "kounen", acc.davies=1e-9)
        sKernel_UW <- UW$pgdf$stat.kernel
        pKernel_UW <- UW$pgdf$pval.kernel
        sBurden_UW <- UW$pgdf$stat.burden
        pBurden_UW <- UW$pgdf$pval.burden
#       call_UW    <- UW$call

        ## read out the results ##
        chr   <- as.character(range_arg[,envir$refCol[1]])
        start <- range_arg[,envir$refCol[2]]
        end   <- range_arg[,envir$refCol[3]]

        result = list(chr, gene, nsnp, start, end,
                   sKernel_BT, pKernel_BT, sBurden_BT, pBurden_BT,
                   sKernel_MB, pKernel_MB, sBurden_MB, pBurden_MB,
                   sKernel_UW, pKernel_UW, sBurden_UW, pBurden_UW)
        lastp1 = nrow(envir$pedgene_results) + 1
        envir$pedgene_results[lastp1,] = result
        if (envir$verbose) {
            print(envir$pedgene_results[lastp1, ])
        }

    } else {
        if (envir$verbose)
            message("Only one markers in range.  Ignored!\n")
    }
}
