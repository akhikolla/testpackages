
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


#' generate the .ped, .fam and .map files of PLINK PED representation of a gwaa.data-class object
#'
#' @description
#'  Use provided gwaa.data-class object and create a PLINK .ped file, PLINK .map file
#'  and a PLINK .phe (phenotypes) file.  By default, \bold{srdta} (a sample from GenABEL)
#'  is used for the gwaa.data-class object.  The files are generated with a prefix of
#'  \emph{srdta} unless a pfx argument is provided.
#'  NOTE: These PLINK files may be used by the Mega2 executable to produce a database.
#'
#' @param gwaa_ name of gwaa.data-class object used as input
#'
#' @param pfx prefix for PLINK .ped/.map/.phe file names
#'
#' @param default name of phenotype used for the 6th column of .ped file
#'
#' @return None
#'
#' @keywords internal
#' @importFrom utils data
#'
#' @examples
#'\dontrun{
#' dmpPed()
#' or
#' dmpPed(mygwaa, "name", "cc")
#'}
dmpPed = function(gwaa_ = srdta, pfx, default = "bt") {

    if (missing(pfx))
        stop("dmpPed can not proceed without a filename prefix argument", call. = FALSE)

    dfphe = data.frame(gwaa_@phdata)
    dfphe$sex = dfphe$id  # don't want sex but need IID and FID; so duplicate id
    names(dfphe)[1:2] = c("FID", "IID")
    write.table(dfphe, file=paste0(pfx, ".phe"), sep="\t", quote=FALSE,
                row.names=FALSE, col.names=TRUE)


    dfmap = data.frame(Chromosome=gwaa_@gtdata@chromosome, Name=gwaa_@gtdata@snpnames,
                       Map.k.a=0, BP.p=gwaa_@gtdata@map)
    write.table(dfmap, file=paste0(pfx, ".map"), sep="\t", quote=FALSE,
                row.names=FALSE, col.names=FALSE)

    dfped = data.frame(pid=gwaa_@gtdata@idnames, person=gwaa_@gtdata@idnames )
    dfped$father = 0
    dfped$mother = 0
    dfped$sex    = gwaa_@phdata$sex
    dfped$sex[dfped$sex == 0] = 2  # here 0 means female
    dfped$default    = gwaa_@phdata[ , default]  # bt is only affection trait


    predfped = sub("/", " ", as.character(gwaa_@gtdata))
    predfped[is.na(predfped)] = "0 0"
    dfped = data.frame(dfped, predfped, stringsAsFactors=FALSE)
    write.table(dfped, file=paste0(pfx, ".ped"), sep="\t", quote=FALSE,
                row.names=FALSE, col.names=FALSE)
}

#' compare two gwaa.data-class objects
#'
#' @description
#'  Verify by fields, all the fields in two gwaa.data-class objects.
#'  Show more detailed marker information iff the coding values are different.  (When comparing
#'  two gwaa.data-class objects, one native and one created via \bold{Mega2R} sometimes
#'  when an allele frequency is .5 for both alleles, the allele order 1/2 vs 2/1 can not be
#'  currently be determined.)
#'
#' @param mega_ name of first gwaa.data-class object
#'
#' @param gwaa_ name of second gwaa.data-class object
#'
#' @param full if TRUE convert genotypes to text as.character(gwaa_@gtdata)\cr and as.character(mega_@gtdata).
#'  Then standardize the order for heterozygous alleles and finally compare.
#'  This step is optional because it can be rather slow.
#'
#' @param envir 'environment' containing SQLite database and other globals
#'
#' @return None
#'
#' @export
#' @importFrom utils data
#'
#' @examples
#'\dontrun{
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' require("GenABEL")
#' ENV = read.Mega2DB(db)
#'
#' y = Mega2ENVGenABEL()
#' Mega2GenABELtst(y, y, full = FALSE)
#'}
#'
#' \dontrun{
#' # donttestcheck: if you have more time, try ...
#' x = Mega2GenABEL()
#' Mega2GenABELtst(x, y, full = FALSE)
#' }
#'
Mega2GenABELtst = function (mega_ = mega, gwaa_ = srdta, full = TRUE, envir = ENV) {
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)

    if (is.null(mega_) || is.null(gwaa_)) {
        warning("One or both GenABEL arguments are NULL.  Aborting.\n")
        return (NULL)
    }
    ANS = TRUE
    phens = names(gwaa_@phdata)
    for (phen in phens[2:length(phens)]) {
        cat("all(mega_@phdata$", phen, " == gwaa_@phdata$", phen, ") ", sep="")
        ans = all(! is.na(mega_@phdata[ , phen]) && ! is.na(gwaa_@phdata[ , phen]) &&
                  mega_@phdata[ , phen] == gwaa_@phdata[ , phen])
        print(ans)
        ANS = ANS && ans
    }

     cat("all(mega_@gtdata@nids == gwaa_@gtdata@nids)")
    ans = all(mega_@gtdata@nids == gwaa_@gtdata@nids)
    print(ans)
    ANS = ANS && ans
    
     cat("all(mega_@gtdata@nsnps == gwaa_@gtdata@nsnps)")
    ans = all(mega_@gtdata@nsnps == gwaa_@gtdata@nsnps)
    print(ans)
    ANS = ANS && ans

     cat("all(mega_@gtdata@nbytes == gwaa_@gtdata@nbytes)")
    ans = all(mega_@gtdata@nbytes == gwaa_@gtdata@nbytes)
    print(ans)
    ANS = ANS && ans

#
     cat("all(mega_@gtdata@idnames == gwaa_@gtdata@idnames)")
    ans = all(mega_@gtdata@idnames == gwaa_@gtdata@idnames)
    print(ans)
    ANS = ANS && ans

     cat("all(mega_@gtdata@snpnames == gwaa_@gtdata@snpnames)")
    ans = all(mega_@gtdata@snpnames == gwaa_@gtdata@snpnames)
    print(ans)
    ANS = ANS && ans

     cat("all(mega_@gtdata@chromosome == gwaa_@gtdata@chromosome)")
    ans = all(mega_@gtdata@chromosome == gwaa_@gtdata@chromosome)
    print(ans)
    ANS = ANS && ans
     cat("all(mega_@gtdata@map == gwaa_@gtdata@map)")
    ans = all(mega_@gtdata@map == gwaa_@gtdata@map)
    print(ans)
    ANS = ANS && ans
     cat("all(mega_@gtdata@male == gwaa_@gtdata@male)")
    ans = all(mega_@gtdata@male == gwaa_@gtdata@male)
    print(ans)
    ANS = ANS && ans
#
     cat("all(mega_@gtdata@coding == gwaa_@gtdata@coding)")
    ansc = all(mega_@gtdata@coding == gwaa_@gtdata@coding)
    print(ansc)
    ANS = ANS && ansc
     cat("all(mega_@gtdata@strand == gwaa_@gtdata@strand)")
    ans = all(mega_@gtdata@strand == gwaa_@gtdata@strand)
    print(ans)
    ANS = ANS && ans

     cat("all(mega_@gtdata@gtps == gwaa_@gtdata@gtps)")
    ans = all(mega_@gtdata@gtps == gwaa_@gtdata@gtps)
    print(ans)
    ANS = ANS && ans
    
    if (full) {
        ms = as.character(gwaa_@gtdata)
        mm = as.character(mega_@gtdata)

        ms[ms == "T/G"] = "G/T"
        ms[ms == "T/C"] = "C/T"
        ms[ms == "T/A"] = "A/T"
        ms[ms == "G/C"] = "C/G"
        ms[ms == "G/A"] = "A/G"
        ms[ms == "C/A"] = "A/C"
        ms[is.na(ms)]   = "0/0"

        mm[mm == "T/G"] = "G/T"
        mm[mm == "T/C"] = "C/T"
        mm[mm == "T/A"] = "A/T"
        mm[mm == "G/C"] = "C/G"
        mm[mm == "G/A"] = "A/G"
        mm[mm == "C/A"] = "A/C"
        mm[is.na(mm)]   = "0/0"

         cat("all(mega_@gtdata == gwaa_@gtdata)")
        ans = all(mm == ms)
        print(ans)
        ANS = ANS && ans
    }

    print(rep(ANS, 10))
    print(rep(ANS, 10))

    if (! ansc) {
        allele_table = envir$allele_table[envir$allele_table$locus_link %in% envir$markers$locus_link,]
        mm = merge(x=allele_table[allele_table$indexX == 1,],
               y=allele_table[allele_table$indexX == 2,],
               by="locus_link")
        cd = which(mega_@gtdata@coding != gwaa_@gtdata@coding)
        print("markers that differ")
        print("markers that differ")
        print(envir$markers[cd,])

        print("allele values for markers that differ")
        print("allele values for markers that differ")
        print(mm[cd,])
    }

}
