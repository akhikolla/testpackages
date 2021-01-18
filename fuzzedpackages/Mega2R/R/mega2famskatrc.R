
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

#library(famSKATRC)
#library(kinship2)

#' load Mega2 SQLite database and perform initialization for famSKATRC usage
#'
#' @description

#'  This populates the \bold{R} data frames with the specified \bold{Mega2} SQLite database.  It initializes
#'  the fam(ily) table and makes sure the person entries are unique.  Finally, it generates a kinship matrix
#'  from the family data.  It also stores a weighting for the common and rare variant that may be used
#'  later if NULL is specified as a weight in \emph{Mega2famSKATRC}.  The common weighting is the function
#'  \emph{dbeta}(maf, 1, 25).  The rare weighting is the function \emph{dbeta}(maf, 0.5, 0.5).
#'
#' @param db specifies the path of a \bold{Mega2} SQLite database containing study data.
#'
#' @param verbose TRUE indicates that diagnostic printouts should be enabled.
#'  This value is saved in the returned environment.
#'
#' @param ALPHA TRUE indicates that two runs of famSKAT_RC  should be enabled.
#'  One with ALPHA numeric ID's and one with numeric IDs ... this is temporary.
#'  The default is FALSE.
#'
#' @param ... fed to \emph{dbmega2_import()}; should be bpPosMap= to select from the maps of
#'  base pairs, if the default is not desired.
#'
#' @return "environment" containing data frames from an SQLite database and some computed values.
#'
#' @importFrom kinship2 kinship
#' @importFrom stats    dbeta
#' @export
#'
#' @note
#'  \emph{init_famSKATRC} creates a new data frame, \emph{envir$phe}, containing phenotype observations.
#'  In addition, it initializes a matrix to aid
#'  in translating a genotype allele matrix to a genotype count matrix.
#'
#'  It also initializes the data frame \emph{envir$famSKATRC_results} to zero rows.
#'
#' @seealso \code{\link{Mega2famSKATRC}}
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = init_famSKATRC(db, verbose = FALSE)
#' ls(ENV)
#'
init_famSKATRC = function (db = NULL, verbose = FALSE, ALPHA = FALSE, ...) {
# id, fullkins, acc
    if (is.null(db))
        stop("You must specify a database argument!\n", call. = FALSE)

    envir = dbmega2_import(db, verbose = verbose, ...)

    fam = mkfam(envir = envir)
    setfam(fam, envir = envir)  # also updates unified_genotype_table
    fam = uniqueFamMember(envir = envir)
    setfam(fam, envir = envir)  # also updates unified_genotype_table

    envir$phe  = mkphenotype(envir)
    envir$phe[envir$phe == 0] = NA

    envir$famSKATRC_results = data.frame(chr = character(0), gene = character(0),
                                     nvariants = numeric(0), start = integer(0), end = integer(0),
                                     stringsAsFactors = FALSE)

    sample = nrow(fam)
    envir$KIN = kinship(1:sample, match(fam$Father, fam$PerPre, nomatch=0),
                                  match(fam$Mother, fam$PerPre, nomatch=0))
    envir$ALPHA = ALPHA
    if (ALPHA)
        envir$KINL = kinship(fam$PerPre, fam$Father, fam$Mother)

    envir$wuweights_r = function(maf) {ifelse(maf > 0, dbeta(maf, 1, 25), 0) }
    envir$wuweights_c = function(maf) {ifelse(maf > 0, dbeta(maf, 0.5, 0.5), 0) }
    return (envir)
}

#' execute the CRAN famSKAT_RC function on a subset of the gene transcripts
#'
#' @description
#' If the \emph{gene} argument is NULL, execute the famSKAT_RC function on the first \emph{gs}
#'  gene transcripts (default is gs = 1:100).
#'  Update the \emph{envir$famSKATRC_results} data frame with the results.
#'  Otherwise, \emph{gene} is a string vector of genes to process.  The special value '*' stands
#'  for all the known genes.
#'
#' @param gs a subrange of the default transcripts (refRanges) over which to calculate
#'  the \emph{DOfamSKATRC} function.
#'
#' @param genes a list of genes over which to calculate the \emph{DOfamSKATRC} function.
#'  The value, "*", means use all the transcripts in the selected Bioconductor database.
#'  If genes is NULL, the gs range of the internal \emph{refRanges} will be used.
#'
#' @param envir 'environment' containing SQLite database and other globals.
#'
#' @param ... extra arguments that are acceptable to famSKAT_RC.  These are listed with the
#'  \code{\link{DOfamSKATRC}} function.
#'
#' @return
#'  The data frame with the results is stored in the environment and named \emph{famSKATRC_results},
#'  viz. envir$famSKATRC_results
#'
#' @export
#'
#' @note
#'  A helper function
#'  \code{SKAT3arg} is defined for the 3 argument callback function which in turn calls
#'  \code{DOfamSKATRC} with the appropriate arguments (including those specific to the
#'  \code{Mega2famSKATRC} function).
#'
#' @seealso \code{\link{init_famSKATRC}}, \code{\link{DOfamSKATRC}}
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = init_famSKATRC(db, verbose = FALSE)
#' ENV$verbose = FALSE
#' ENV$famSKATRC_results = ENV$famSKATRC_results[0, ]
#' Mega2famSKATRC(gs=50:60, envir=ENV, pheno=3)
#'
#' \donttest{
#' # donttestcheck: try this below if there is time
#'  Mega2famSKATRC(genes=c("CEP104"), envir=ENV, pheno=3 )
#' }
#'
#' ENV$famSKATRC_results
#'
#Mega2famSKATRC = function (gs = 1:100, genes=NULL, envir = ENV,
#                           pheno=3, id=NULL, covariates=NULL,
#                           sqrtweights_c=NULL, sqrtweights_r=NULL,
#                           binomialimpute=TRUE, acc=1e-6, maf=0.05, phi=c(0,0.2,0.5,0.9) )
Mega2famSKATRC = function (gs = 1:100, genes=NULL, envir = ENV, ...) {
# PHENO genotypes covariates sqrtweights_c sqrtweights_r binomialimpute maf phi acc
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)

    famSKAT3arg = function(markers_arg, range_arg, envir) {
        DOfamSKATRC(markers_arg, range_arg, envir, ...)
    }

    if (is.null(genes))
        applyFnToRanges(famSKAT3arg, envir$refRanges[gs, ], envir$refIndices, envir = envir)
    else
        applyFnToGenes(famSKAT3arg, genes, envir = envir)
}

#' DofamSKATRC call back function
#'
#' @description
#'  Convert the genotypesraw() allele patterns of 0x10001, 0x10002 (or 0x20001), 0x20002, 0
#"  from the genotype matrix
#'  to the numbers 0, 1, 2, 9 for each marker. (Reverse, the order iff allele "1" has the
#'  minor allele frequency.)  Ignore markers that have no variants.
#'  Finally, invoke \code{famSKAT_RC} with the converted genotype matrix.
#'  Save information about the range and the p.value calculated by \code{famSKAT_RC}
#'  in \emph{envir$famSKATRC_results}.
#'  If you want to change the argument values to this function they should be changed instead
#'  when calling the \code{Mega2famSKATRC} function.
#'
#' @param markers_arg a data.frame with the following 5 observations:
#' \describe{
#' \item{locus_link}{is the ordinal ranking of this marker among all loci}
#' \item{locus_link_fill}{is the position of corresponding genotype data in the
#' \emph{unified_genotype_table}}
#' \item{MarkerName}{is the text name of the marker}
#' \item{chromosome}{is the integer chromosome number}
#' \item{position}{is the integer base pair position of marker}
#' }
#'
#' @param range_arg one row of a ranges_arg.  The latter is a data frame of at least three
#'  integer columns.  The columns indicate a range:
#'  a chromosome number, a start base pair value, and an end base pair value.
#'
#' @param envir 'environment' containing SQLite database and other globals especially the
#'  phenotype_table, \code{phe}.
#'
#' @param pheno is an index into the phenotypes_table to select the phenotype.  Missing phenotypes
#'  are represented by NA.
#'
#' @param covariates a matrix of covariates for the phenotype.
#'
#' @param id a vector of individuals to be included in the test, a subset of the family members.
#'  If NULL is given, all members will be used.
#'
#' @param binomialimpute if TRUE, impute missing genotypes using a binomial distribution.
#'
#' @param sqrtweights_c weight function for common variants, if NULL use weight set in init_famSKAT
#'
#' @param sqrtweights_r weight function for rare variants, if NULL use weight set in init_famSKAT.
#'
#' @param maf threshold used to separate rare from common variants.
#'
#' @param phi a vector of ratios ratios; each indicates the contribution of rare variants.
#'
#' @param acc accuracy used in Davies approximation.
#'
#' @importFrom famSKATRC famSKAT_RC
#' @return None
#'
#' @export
#'
#' @seealso \code{\link{init_famSKATRC}}, \code{\link{Mega2famSKATRC}}
#'
#' @note
#'  This function accumulates output in the data frame, \emph{envir$famSKATRC_results}.  It will
#'  print out the lines as they are generated if \emph{envir$verbose} is TRUE.  It does not write
#'  the data frame to a file.  You must save the data frame.
#'  You also must initialize the data frame when necessary.
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = init_famSKATRC(db, verbose = TRUE)
#' ENV$famSKATRC_results = ENV$famSKATRC_results[0, ]
#' Mega2famSKATRC(gs=1:1, envir=ENV, pheno=3)
#' # this sets one of the many arguments for DOfamSKATRC
#' # but basically prepares the ENV for the direct use of DOfamSKATRC (below).
#'
#' \donttest{
#' # donttestcheck: try this below instead if there is time
#'  Mega2famSKATRC(genes=c("CEP104"), envir=ENV, pheno=3 )
#' }
#'
#' # DOfamSKATRC is called within Mega2famSKATRC. init_famSKATRC and Mega2famSKATRC need to be
#' # called to set up the environment for famSKAT_RC to run.  BUT, you should ignore DOfamSKATRC
#' # and use Mega2famSKATRC instead.
#' #
#' applyFnToRanges(DOfamSKATRC, ENV$refRanges[50:60, ], ENV$refIndices, envir=ENV)
#' # this will use all the default argument values for DOfamSKATRC
#'
#
DOfamSKATRC = function(markers_arg, range_arg, envir,
                       pheno=3, id=NULL, covariates=NULL,
                       sqrtweights_c=NULL, sqrtweights_r=NULL,
                       binomialimpute=TRUE, acc=1e-6,
                       maf=.05, phi=c(0, .2, .5, .9) ) {

    if (is.null(range_arg))
        stop("DOfamSKATRC: range is not defined.", calls. = FALSE)
    if (pheno > ncol(envir$phe)) stop("pheno index exceeds phenotype array")

    geno = computeDosage(markers_arg, range_arg, envir)
    nsnps = ncol(geno)

    if (nsnps == 0) return (NULL)
    if (nsnps == 1) {
       message("Only 1 marker found in range ",
               range_arg[envir$refCol[4]], ": chr", range_arg[envir$refCol[1]], " ",
               range_arg[envir$refCol[2]], " ", range_arg[envir$refCol[3]])
       return (NULL)
    }

#   The code below is 2 to 3 times faster than that above but is all in Rcpp.  The above
#   code makes it clearer that markers that have no variant are ignored and that the first
#   allele is the one with the most occurrences, and that there must be at least two markers.

#   tim2 = system.time ({
#   genocol = getgenotypesdos(markers_arg, envir)
#   genodos = genocol$geno[ ,1:genocol$ncol]
##  genodos = matrix(genocol$geno, nrow=nrow(genocol$geno), ncol=genocol$ncol)
#   })
#   if (kk != genocol$ncol) {
#       message("kk vs ncol")
#    }
#   if (sum(geno-genodos) != 0) {
#       message("geno vs genodos")
#   }

    chr   = as.character(range_arg[,envir$refCol[1]])
    gene  = as.character(range_arg[,envir$refCol[4]])
    start = range_arg[,envir$refCol[2]]
    end   = range_arg[,envir$refCol[3]]

    if (is.null(id)) {
        id = 1:length(envir$fam$PerPre)
        if (envir$ALPHA)
            idL = envir$fam$PerPre
    }
    if (is.null(sqrtweights_c)) sqrtweights_c = envir$wuweights_c
    if (is.null(sqrtweights_r)) sqrtweights_r = envir$wuweights_r
    argn = list(PHENO=envir$phe[,pheno], genotypes=geno, id=id, fullkins=envir$KIN,
               covariates=covariates, sqrtweights_c=sqrtweights_c, sqrtweights_r=sqrtweights_r,
               binomialimpute=binomialimpute, acc=acc, maf=maf, phi=phi)

    tim = system.time ({
        skat = do.call(famSKAT_RC, argn)
    })

    if (envir$ALPHA) {
##
        argn = list(PHENO=envir$phe[,pheno], genotypes=geno, id=idL, fullkins=envir$KINL,
                   covariates=covariates, sqrtweights_c=sqrtweights_c, sqrtweights_r=sqrtweights_r,
                   binomialimpute=binomialimpute, acc=acc, maf=maf, phi=phi)

        timL = system.time ({
            skatL = do.call(famSKAT_RC, argn)
        })

        result = c(list(), chr=chr, gene=gene, nvariants=nsnps, start=start, end=end,
                   tim[1:3], famSKATRC=skat, timL[1:3], skatL=skatL)
##
    } else
        result = c(list(), chr=chr, gene=gene, nvariants=nsnps, start=start, end=end,
                   tim[1:3], famSKATRC=skat)

    lastp1 = nrow(envir$famSKATRC_results) + 1
    if (lastp1 == 1) {
        envir$famSKATRC_results = do.call(data.frame, c(result, stringsAsFactors = FALSE))
      } else {
        envir$famSKATRC_results[lastp1,] = result
    }

    if (envir$verbose) {
        print(envir$famSKATRC_results[lastp1, ])
    }
}
