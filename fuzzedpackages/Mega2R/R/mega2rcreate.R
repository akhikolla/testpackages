
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
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-101, USA.
#
#   For further information contact:
#       Daniel E. Weeks
#       e-mail: weeks@pitt.edu
#
# ===========================================================================

#' Mega2R package
#'
#' @description This package reads a Mega2 SQLite3 database into data frames and
#'	makes the contained genotypes/phenotypes/linkage data available for analysis.
#'
#' @author Robert V. Baron and Daniel E. Weeks
#' @docType package
#' @name Mega2R
NULL

#library(DBI)
#library(RSQLite)

#' Mega2R version
#'
#' @description This string indicates the current release of Mega2R
#'
#' @author Robert V Baron
#' @docType data
#' @name Mega2RVersion
Mega2RVersion = '6.0.0'

#
# how to load refRanges & refIndices to sysdata
#
# setwd("Mega2R")
# devtools::use_data(refRanges, refIndices, internal = TRUE, overwrite = TRUE)
# setwd("..")
# tools::resaveRdaFiles("Mega2R/R")

#' Mega2R SQLite3 tables
#'
#' @description This character vector indicates the names of the Mega2 SQLite3 database tables
#'  to load.  (Not all of the existing tables are loaded.)
#'
#' @author Robert V Baron
#' @docType data
#' @name Mega2R-TBLS
TBLS = c("int_table",
#        "double_table",
         "charstar_table",
#        "stuff_table",
#        "batch_parameters",
#        "file_table",

         "pedigree_table",
         "person_table",
         "pedigree_brkloop_table",
         "person_brkloop_table",

#        "canonicalallele_table",
         "locus_table",
         "allele_table",
         "marker_table",
#        "markerscheme_table",
         "map_table",
         "mapnames_table",

         "traitaff_table",
         "affectclass_table",
#        "traitquant_table",

         "phenotype_table",
         "genotype_table"

## Derived tables
##        unified_genotype_table
##        markers
  )

#' Mega2R SQLite3 table filter
#'
#' @description This list contains named values.  The name corresponds to an SQLite database table.
#'  The value is a character string of column names from the "named" table that should be fetched.  A table
#'  is in this list, if not all the database table columns are needed.  The columns for each table are separated
#'  by commas.
#'
#' @author Robert V Baron
#' @docType data
#' @name Mega2R-TBLSFilter
#' @note For the data base tables not in this list, all columns are stored in the corresponding data frame.
TBLSFilter = list(
          locus_table    = "pId, LocusName, Type, AlleleCnt, locus_link",

          marker_table   = "pId, MarkerName, pos_avg, pos_female, pos_male, chromosome, locus_link",
          pedigree_table = "pId, Num, EntryCnt, Name, PedPre, OriginalID, origped, pedigree_link",
          person_table   = "pId, UniqueID, OrigID, FamName, PerPre, ID, Father, Mother, Sex, pedigree_link, person_link",
          pedigree_brkloop_table = "pId, Num, EntryCnt, Name, PedPre, OriginalID, origped, pedigree_link",
          person_brkloop_table   = "pId, UniqueID, OrigID, FamName, PerPre, ID, Father, Mother, Sex, pedigree_link, person_link",

          traitaff_table = "pId, ClassCnt, PenCnt, locus_link"
  )

#' make the derived "markers" data frame and store it in the environment.
#'
#' @description Create the markers data frame as a subset of the markers_table data frame.  It contains 5
#'	observations: \describe{
#'	  \item{locus_link:}{locus offset of this marker}
#'	  \item{locus_link_fill:}{locus offset plus an accumulating fudge factor that jumps
#'          with each new chromosome because the count of markers per chromosome is force to be
#'          a multiple of 4.  (This value corresponds to the offset of the marker in the
#'          \code{unified_genotype_table}.)}
#'	  \item{MarkerName:}{name of this marker}
#'	  \item{chromosome:}{chromosome number of this marker}
#'	  \item{position:}{base pair position of this marker (selected by bpPosMap[below])}
#'       }
#'
#' @param bpPosMap An integer that indicates the map (index) to use to fetch the
#'	chromosome/position fields from the map_table data frame to merge with the marker_table.
#'
#' @param envir an environment that contains all the data frames created from the SQLite database.
#'
#' @return None
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' mk_markers_with_skip(1)
#'}
mk_markers_with_skip = function(bpPosMap = 1, envir) {
    if (envir$MARKER_SCHEME == 1) {
        markersPerChr = sapply(split(envir$marker_table$chromosome, envir$marker_table$chromosome), length)
        idx = as.integer(names(markersPerChr))
        m = rep(0, max(idx))
        m[idx] = markersPerChr
        markersPerChr = m
        extra_markers = cumsum(4 * floor((markersPerChr + 3) / 4) - markersPerChr)
        extra_markers = c(0, extra_markers)
        names(extra_markers) = NULL
    } else if (envir$MARKER_SCHEME == 2) {
        extra_markers = vector("integer", length(unique(envir$marker_table$chromosome))+1)
    }
    envir$marker_table$locus_link_fill = envir$marker_table$locus_link + extra_markers[envir$marker_table$chromosome]
    if (any(is.na(envir$marker_table$locus_link_fill))) {
        stop("Internal Error in mk_markers_with_skip.  Bad fill values.", call. = FALSE)
    }
    mkMarkers(bpPosMap = bpPosMap, envir = envir)
}

#' create "markers" data frame
#'
#' @description Create the markers data frame.  It contains 5
#'	observations: \describe{
#'	  \item{locus_link:}{locus offset of this marker}
#'	  \item{locus_link_fill:}{locus offset plus an accumulating fudge factor that jumps
#'          with each new chromosome because the count of markers per chromosome is force to be
#'          a multiple of 4.  (This value corresponds to the offset of the marker in the
#'          \code{unified_genotype_table}.)}
#'	  \item{MarkerName:}{name of the marker}
#'	  \item{chromosome:}{chromosome number of the marker}
#'	  \item{position:}{base pair position of the marker (selected by bpPosMap[below])}
#'       }
#'
#' @details Select a map (index) from the map_table to merge with the select marker_table
#'  data frame to make the marker data frame.  See showMapNames() for the string name to index mapping.
#'
#' @param bpPosMap An integer that indicates the map (index) to use to merge the
#'   chromosome/position fields from the map_table data frame to the marker_table data frame.
#'   See showMapNames() for the string name to index mapping.
#'
#' @param envir an environment that contains all the data frames created from the SQLite database.
#'
#' @return None
#' @export
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db, verbose = FALSE)
#'
#' mkMarkers(1)
#'
#' ENV$markers
#'
mkMarkers = function(bpPosMap = 1, envir = ENV) {
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)

    map_table = envir$map_table[ envir$map_table$map == bpPosMap, c( "marker", "position")]
    if (nrow(map_table) == 0) {
        message("No entry for map == ", bpPosMap, " in map_table.  Using map == 0 instead.")
        map_table = envir$map_table[ envir$map_table$map == 0, c( "marker", "position")]
    }
    envir$markers = merge(envir$marker_table[ , c("locus_link", "locus_link_fill", "MarkerName", "chromosome")],
                        map_table,
                        by.x = "locus_link", by.y = "marker")
}

#' concatenate separate genotype vectors for each chromosome to one extended vector containing all the chromosomes
#'  and store it in the environment.
#'
#' @description The genotype_table data frame contains for each person a separate record for each chromosome.
#'  Concatenate all the separate records for each person in chromosome order and make a data frame of these vectors.
#'  Note:  Genotype vectors are "raw" (byte) vectors.
#'
#' @param envir an environment that contains all the data frames created from the SQLite database.
#'
#' @return None
#'
#' @keywords internal
#'
mk_unified_genotype_table = function(envir) {
    samples = split(envir$genotype_table, envir$genotype_table$person_link)
    samplesize = length(samples)
    person_link = unique(envir$genotype_table$person_link)

    df = data.frame(row.names = 1:samplesize,
                    person_link = person_link,
                    data = vector("raw", samplesize))

    for (i in 1:samplesize) {

        chrOrder = order(samples[[i]][ , 3], decreasing = FALSE)
        v = samples[[i]][chrOrder, 5]
        if (envir$DBcompress && ! is.na(v)) {
            xx = lapply(v, function(a) memDecompress(a, "gzip"))
            uv = unlist(xx)
        } else {
            uv = unlist(v)
        }
        df$data[i] = list(uv)
    }
    envir$unified_genotype_table = df
    envir$genotype_table = NULL
    rm(list = "genotype_table", envir = envir)
}

#' read Mega2 SQLite database into R
#'
#' @description Read the fields of SQLite data base tables that are required for Mega2R into
#'  data frames.
#'  These data frames are stored in an 'environment' which is returned.
#'  This function also adds some state data, extra data frames, and computed data frames
#'  to the 'environment'.
#'
#' @usage
#' dbmega2_import(dbname,
#'                bpPosMap = NULL,
#'                verbose = FALSE)
#'
#' @param dbname file path to SQLite database.
#'
#' @param bpPosMap index that specifies which map in the map_table should be used for
#'  marker chromosome/position.  If it is NULL, the internal variable
#'  \emph{base_pair_position_index} is used instead.
#'  \code{showMapNames()} shows the association between map name and map number.
#'
#' @param verbose print out statistics on the name/size of each table read and show column headers.
#'  Also, save the verbose value for use by other Mega2R functions.
#'
#' @return envir an environment that contains all the data frames made from the SQLite database.
#'
#' @importFrom RSQLite dbConnect dbExistsTable dbReadTable dbListFields SQLITE_RO
#' @importFrom DBI dbGetQuery dbDisconnect
#' @export
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = dbmega2_import(db, verbose = TRUE)
#'
#' ENV = dbmega2_import(db)
#
dbmega2_import = function(dbname,
                          bpPosMap = NULL,
                          verbose = FALSE) {
    con = tryCatch(dbConnect(RSQLite::SQLite(), dbname = dbname, flags = SQLITE_RO),
                   error = function(xx) { stop("DB open failed: ", dbname, call. = FALSE) })

    envir = resetMega2ENV()

    gc(verbose = FALSE)

    envir$verbose = verbose

    for (tbl in TBLS) {
        if (dbExistsTable(con, tbl)) {
            filter = TBLSFilter[tbl][[1]]
            if (is.null(filter))
                assign(tbl, dbReadTable(con, tbl), pos = envir)
            else {
#               assign(tbl, dbReadTable(con, tbl, select.cols = filter), pos = envir)
#               assign(tbl, dbReadTable(con, tbl), pos = envir)
#               print(paste0("select ", filter, " from ", tbl));
                assign(tbl, dbGetQuery(con, paste0("select ", filter, " from ", tbl)), pos = envir)
            }
            if (envir$verbose) {
                cat(tbl, dim(get(tbl, pos=envir)), sep = "\t", end = "\n")
                if (is.null(filter))
                    cat(tbl, dbListFields(con, tbl), sep = "\t", end = "\n")
                else
                  cat(tbl, strsplit(filter, ", ")[[1]], sep = "\t", end = "\n")
                cat(end = "\n")
            }
        }
    }

    dbDisconnect(con)

    envir$PhenoCnt       = envir$int_table[envir$int_table$key == 'PhenoCnt', 3]
    envir$LocusCnt       = envir$int_table[envir$int_table$key == 'LocusCnt', 3]
    envir$MARKER_SCHEME  = envir$int_table[envir$int_table$key == 'MARKER_SCHEME', 3]
    envir$DBcompress     = envir$int_table[envir$int_table$key == 'dbCompression', 3]
    envir$DBMega2Version = envir$charstar_table[envir$charstar_table$key == 'DBMega2Version', 3]
    
    lt3 = function(aa, bb) {
                            ( (aa[1] < bb[1]) || ( (aa[1] == bb[1]) &&
                                            ( (aa[2] < bb[2]) || (aa[2] == bb[2] && (aa[3] < bb[3])) ) ) )
            }
    lt2 = function(aa, bb) {
                            ( (aa[1] < bb[1]) || ( (aa[1] == bb[1]) && (aa[2] < bb[2]) ) ) 
            }
    cc = envir$DBMega2Version != 'X.Y.Z'
    if (cc) {
        aa = as.numeric(strsplit(envir$DBMega2Version, split=".", fixed=T)[[1]])
        bb = as.numeric(strsplit(Mega2RVersion,        split=".", fixed=T)[[1]])
        lt = lt2(aa, bb)
        if (lt) {
           message("NOTE: Mega2R cannot read the Mega2 database because its version (", envir$DBMega2Version, ") is too old.")
           stop("Please recreate the database using the current version of Mega2.", call.=FALSE)
        }
        lt = lt2(bb, aa)
        if (lt) {
           message("NOTE: Mega2R cannot read the Mega2 database because Mega2R version (", Mega2RVersion, ") is too old.")
           stop("Please get the latest version of Mega2R from CRAN.", call.=FALSE)
        }
    }

    if ( length(envir$DBcompress) == 0) envir$DBcompress = 0
    
    if (envir$MARKER_SCHEME > 2) {
        stop("Only compressions levels of 1 or 2 are allowed. (",
             envir$MARKER_SCHEME, ")", call. = FALSE)
    }

    if (is.null(bpPosMap)) {
        bpPosMap = envir$int_table[envir$int_table$key == 'base_pair_position_index', 3]
        if (bpPosMap < 0) {
            message("NOTE: BAD base pair map; trying map 0 ");
            bpPosMap = 0
        }
    }
    mk_unified_genotype_table(envir)
    mk_markers_with_skip(bpPosMap, envir)
    if (envir$MARKER_SCHEME == 2) {
        message("Partitioninging allele_table by locus_link\n");
        envir$locus_allele_table = split(envir$allele_table, envir$allele_table$locus_link)
    } else {
        envir$locus_allele_table = NULL
    }

    gc(verbose = FALSE)
    return(envir)
}

#' generate a phenotype data frame
#'
#' @description
#'  Convert data in phenotype_table to a data frame of columns that are phenotypes.
#'  The columns may be affection status or quantitative values
#'
#' @param envir "environment" containing SQLite database and other globals
#'
#' @return is a data frame with FID column, then IID column, and then
#'  an additional column for each phenotype.
#'
#' @export
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db)
#' out = mkphenotype()
#'
#' out
#
mkphenotype = function (envir) {
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)

# linkage.h:    TYPE_UNSET, QUANT, AFFECTION, BINARY, NUMBERED, XLINKED, YLINKED
#                        0      1          2       3         4        5        6

    out = envir$fam[3:4]
    hdr = c("FID", "IID")

    phenotype_table = envir$phenotype_table

    raw = unlist(envir$phenotype_table[,4])
    raw = matrix(raw, ncol=8, byrow=T)
    nrows     = nrow(raw)
    nrowpheno = nrow(out)

    type = integer(envir$PhenoCnt)
    for (i in 1:envir$PhenoCnt) {
        hdr = c(hdr, envir$locus_table[i, 2]) # 2 == LocusName

# phenotype_table contains a blob which is a list of entries.  An entry is either an 8 byte
#  double for quant, or two 4 byte ints for affect
        if (envir$locus_table[i, 3] == 2) {              # 3 == Type === AFFECTION
            type[i] = 2
            col = vector("integer", nrowpheno)
            for (j in 1:nrowpheno) {
                col[j] = readBin(raw[envir$PhenoCnt*(j-1)+i, 1:4], integer(), n=1, size=4)
            }
#           col[col==0] = NA
            out$col = col
            names(out) = hdr
        } else if (envir$locus_table[i, 3] == 1) {       # 3 == Type === QUANT
            type[i] = 1
            col = vector("numeric", nrowpheno)
            for (j in 1:nrowpheno) {
                col[j] = readBin(raw[envir$PhenoCnt*(j-1)+i, 1:8], numeric(), n=1, size=8)
            }
            col[col==-99] = NA
            out$col = col
            names(out) = hdr
        }
    }
    attr(out, "type") = type
    out
}

#' show Mega2R environment, viz. data frames and related info.
#'
#' Mega2R uses an environment to store the data frames when it reads SQLite database tables.
#'  This function shows the data frames and their sizes; it also
#'  shows the count of samples and markers in the database.
#'  Note: It is not necessary to provide an argument, if the environment is named \emph{ENV}.
#'
#' @param envir an environment that contains all the data frames created from the SQLite database.
#'
#' @return None
#' @export
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db)
#'
#' showMega2ENV()
#'
showMega2ENV = function(envir = ENV) {
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)
    cat("locus count:  ",       envir$LocusCnt)
    cat("; phenotype count: ", envir$PhenoCnt)
    if (envir$MARKER_SCHEME == 1)
        cat("; compression: 2 bits")
    else if (envir$MARKER_SCHEME == 2)
        cat("; compression: 2 bytes")
    else
        cat("; compression: illegal")
    cat("\n")
    cat("marker count: ",      envir$LocusCnt - envir$PhenoCnt)
    cat("; sample count: ",    nrow(envir$fam))
    cat("\n")
    cat("\n")

    cat("genetic and physical maps:\n")
    map = envir$mapnames_table[ , c(6, 2)]
    names(map) = c("map name", "map number")
    print(map)
    cat("\n")

    cat("Phenotypes:\n")
    print(showPhenoNames(envir))
    cat("\n")
    cat("\n")

    cat("basic tables:\n")
    tbls2 = c("fam", "locus_allele_table", "markers", "unified_genotype_table")
    ls = grep("_table", ls(envir), value = TRUE)
    tbls =  sort( ls[! ls %in% tbls2 ] )
    x = sapply(tbls, function(tbl) {dim(get(tbl, envir))})
    x = t(x)
    x = matrix(x, ncol = 2)
    show = data.frame(x)
    names(show) = c("rows", "cols")
    row.names(show) = tbls
    print(show)
    cat("\n")

    cat("derived tables:\n")
    if (exists("fam", envir))
        tbls2 = c("fam", "markers", "unified_genotype_table")
    else
        tbls2 = c("markers", "unified_genotype_table")
    x = sapply(tbls2, function(tbl) {dim(get(tbl, envir))})
    x = t(x)
    x = matrix(x, ncol = 2)
    show = data.frame(x)
    names(show) = c("rows", "cols")
    row.names(show) = tbls2
    print(show)
    cat("\n")
}

#' return an initialized environment
#'
#' Mega2 uses an environment to store the data frames when it reads an SQLite database.
#'  The environment is also used to store metadata.
#'  The function first runs the garbage collector ("gc"), then allocates an empty environment
#'  and finally loads some default data into it.
#'
#' @return an environment that contains a few initial tables read from the Mega2R package.
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' resetMega2ENV()
#'}
resetMega2ENV = function () {

    gc(verbose = FALSE)

    envir = new.env(parent = emptyenv())

    envir$Mega2R     = environment(resetMega2ENV)
    
    envir$refRanges  = refRanges
    envir$refIndices = refIndices

    envir$txdb       = "TxDb.Hsapiens.UCSC.hg19.knownGene"
    envir$entrezGene = "org.Hs.eg.db"

    envir$chr2int = data.frame(chr = c(1:26, 23:26, 1:26, 23:26))
    envir$chr2int = cbind(envir$chr2int, stringsAsFactors = FALSE,
                          string = as.character(envir$chr2int$chr))
    envir$chr2int[27:30,2] = c("X", "Y", "XY", "M")
    envir$chr2int[31:60,2] = paste0("chr", envir$chr2int[1:30,2])

    envir$positionVsName = FALSE

    envir$dosageRaw = c(0x10001, 0x10002, 0x20001, 0x20002, 0)
    envir$dosage    = c(      0,       1,       1,       2, 0)

    return(envir)
}

#' show the association between mapno and mapname
#'
#' Mega2R allows several different physical and genetic maps to be stored and used to select
#'  positions.  This function shows the association between map number and map name.
#'
#' @param envir an environment that contains all the data frames created from the SQLite database.
#'
#' @return None
#' @export
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db)
#'
#' showMapNames()
#'
showMapNames = function (envir = ENV) {
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)
    envir$mapnames_table[ , c(6, 2)]
}

#' show the association between index no and phenotype
#'
#' Mega2R stores several phenotypes, both affective and quantitative. This function displays the
#'  mapping between phenotype (name), index, and the phenotype type (affection or quantitative).
#'
#' @param envir an environment that contains all the data frames created from the SQLite database.
#'
#' @return None
#' @export
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db)
#'
#' showPhenoNames()
#'
showPhenoNames = function (envir = ENV) {
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)

# linkage.h:    TYPE_UNSET, QUANT, AFFECTION, BINARY, NUMBERED, XLINKED, YLINKED
#                        0      1          2       3         4        5        6
    hdr = NULL
    type = NULL
    df = data.frame(Index = 1:envir$PhenoCnt)
    for (i in 1:envir$PhenoCnt) {
        hdr = c(hdr, envir$locus_table[i, 2]) # 2 == LocusName

        if (envir$locus_table[i, 3] == 2) {              # 3 == Type === AFFECTION
            type = c(type, "affection")
        } else if (envir$locus_table[i, 3] == 1) {       # 3 == Type === QUANT
            type = c(type, "quantitative")
        }
    }
    df$Name = hdr
    df$Type = type
    df
}


#' return the genotypes for all markers of a given person
#'
#' @description for selected person, expand genotype bit vector to an integer vector with
#'  one integer per marker.  Then convert the integer vector to genotypes.
#'
#' @param perid person/sample identifier
#'
#' @param envir an environment that contains all the data frames created from the SQLite database.
#'
#' @return None
#' @export
#'
#' @keywords internal
#'
#' @examples
#'\dontrun{
#' # genotypes for all markers for n'th person in genotype table
#' getgenotype_person(1)
#'
#' # genotypes for all markers for range of persons in genotype table
#' getgenotype_person(m:n)
#'}
getgenotype_person = function(perid = 1, envir = ENV) {

  a1 = envir$allele_table[envir$allele_table$indexX==1, 2][-1]
  a2 = envir$allele_table[envir$allele_table$indexX==2, 2][-1]

  rv = envir$genotype_table[perid, 5][[1]]
  rv4 = getgenotypes_forperson(rv)

#  0 1|1
#  1 0|0
#  2 1|2
#  3 2|2

  return(
    ifelse(rv4==0, paste0(a1, a1),
           ifelse(rv4==1, paste0("00"),
                  ifelse(rv4==2, paste0(a1, a2), paste0(a2, a2))
                  )
           )
        )
}

################################################################

#' fetch genotype matrix for specified markers (assemble by rows)
#'
#' @description
#'  This function calls the C++ function that does all the heavy lifting.  It passes the
#'  locus_index and the locus_offset in the \emph{unified_genotype_table} from the
#'  \emph{markers_arg} argument.  It also gathers other data.frames that are in the "global"
#'  \bold{ENV} environment. One frame contains a bit vector of compressed genotype information,
#'  another contains the alleles for each marker, and finally there are some bookkeeping related
#'  data.  Note this function is for Testing only and is not exported.
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
#' @param sepstr separator string for alleles (default is none)
#'
#' @param envir an environment that contains all the data frames created from the SQLite database.
#'
#' @return a matrix of genotypes represented as a nucleotide pair.  There is one column for each
#'  marker in \emph{markers_arg} argument.  There is one row for each person in the family
#'  (\emph{fam}) table.
#'
#' @keywords internal
#' @useDynLib Mega2R
#'
#' @details
#'  The \emph{unified_genotype_table} contains one raw vector for each person.  In the vector
#'  there are two bits for each genotype.  This function creates an output matrix by selecting
#'  from each row all the needed markers and then repeating for each person.
#'
#' @examples
#'\dontrun{
#' # genotypes for all persons in markers data.frame argument
#' getgenotypes_R(ENV$markers)
#'
#' # genotypes for all persons in chromosome n
#' getgenotypes_R(ENV$markers[ENV$markers$chromosome == n,])
#'}
getgenotypes_R = function(markers_arg, sepstr = "", envir = ENV) {

  return(
    if (envir$MARKER_SCHEME == 1) {
        getgenotypes_Ri(markers_arg$locus_link, markers_arg$locus_link_fill,
                        envir$unified_genotype_table, envir$allele_table,
                        sepstr, envir$PhenoCnt)
    } else {  # must be == 2

    }
   )
}


#' fetch genotype character matrix for specified markers
#'
#' @description
#'  This function calls a C++ function that does all the heavy lifting.  It passes the arguments
#'  necessary for the C++ function: some from the caller's arguments and some from data frames
#'  that are in the "global" environment, \bold{envir}.  From the markers_arg argument, it fetches
#'  the locus_index and the index in the \emph{unified_genotype_table}. It also passes the allele
#'  nucleotide separator argument.
#'  From the "global" environment, \bold{envir}, it gets a bit vector of compressed genotype information,
#'  the alleles for each marker, and some bookkeeping related data.
#'  Note: This function also contains a dispatch/switch on the type of compression in the genotype
#'  vector.  A different C++ function is called when there is compression versus when there is no
#'  compression.
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
#' @param sepstr separator string inserted between the alleles (default is none).  When present, this
#'  is typically a space, a tab or "/".
#'
#' @param envir an environment that contains all the data frames created from the SQLite database.
#'
#' @return a matrix of genotypes represented as two allele pairs.  The matrix has one column for each
#'  marker in \emph{markers_arg} argument.  There is one row for each person in the family
#'  (\emph{fam}) table.
#'
#' @export
#' @useDynLib Mega2R
#'
#' @details
#'  The \emph{unified_genotype_table} contains one raw vector for each person.  In the vector
#'  there are two bits for each genotype.  This function creates an output matrix by fixing
#'  the marker and collecting genotype information for each person and then repeating for
#'  all the needed markers.  (Currently, this appears slightly faster than a scan
#'  which is fixes the person and iterates over markers.)
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db)
#'
#' getgenotypes(ENV$markers)
#'
getgenotypes = function(markers_arg, sepstr = "", envir = ENV) {
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)

  return(
    if (envir$MARKER_SCHEME == 1) {
        getgenotypes_1(markers_arg$locus_link, markers_arg$locus_link_fill,
                       envir$unified_genotype_table, envir$allele_table,
                       sepstr, envir$PhenoCnt)
    } else {  # must be == 2
        getgenotypes_2(markers_arg$locus_link,
                       envir$unified_genotype_table, envir$locus_allele_table,
                       sepstr, envir$PhenoCnt)
    }
   )
}

getgenotypes_C = getgenotypes

#' fetch genotype integer matrix for specified markers
#'
#' @description
#'  This function calls a C++ function that does all the heavy lifting.  It passes the arguments
#'  necessary for the C++ function: some from the caller's arguments and some from data frames
#'  that are in the "global" environment, \bold{envir}.  From its markers_arg argument, it gets
#'  the locus_index and the index in the \emph{unified_genotype_table}.
#'  From the "global" environment, \bold{envir}, it gets a bit vector of compressed genotype information,
#'  and some bookkeeping related data.
#'  Note: This function also contains a dispatch/switch on the type of compression in the genotype
#'  vector.  A different C++ function is called when there is compression versus when there is no
#'  compression.
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
#' @param envir an environment that contains all the data frames created from the SQLite database.
#'
#' @return a matrix of genotypes represented as integers.  Each 32 bit integer represents contains
#'  two allele values: the high 16 bits contains the index of allele1 and the low 16 bits contains
#'  the index of allele2.  In the matrix, there is one column for each
#'  marker in the \emph{markers_arg} argument.  There is one row for each person in the family
#'  (\emph{fam}) table.
#'
#' @export
#' @useDynLib Mega2R
#'
#' @details
#'  The \emph{unified_genotype_table} contains one raw vector for each person.  In the vector,
#'  there are two bits for each genotype.  This function creates an output matrix by fixing
#'  the marker and collecting genotype information for each person and then repeating for
#'  all the needed markers.
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db)
#'
#' # two ints in upper/lower half integer representing allele # for all persons in chromosome 1
#' getgenotypesraw(ENV$markers[ENV$markers$chromosome == 1,])
#'
getgenotypesraw = function(markers_arg, envir = ENV) {
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)

  return(
    if (envir$MARKER_SCHEME == 1) {
        getgenotypesraw_1(markers_arg$locus_link, markers_arg$locus_link_fill,
                          envir$unified_genotype_table, envir$allele_table,
                          envir$PhenoCnt)
    } else {  # must be == 2
        getgenotypesraw_2(markers_arg$locus_link,
                          envir$unified_genotype_table, envir$locus_allele_table,
                          envir$PhenoCnt)
    }
   )
}

#' process the genotype matrix for specified markers and return the corresponding GenABEL genotype matrix
#'
#' @description
#'  This function calls a C++ function that does all the heavy lifting.  It passes the arguments
#'  necessary for the C++ function: some from the caller's arguments and some from data frames
#'  that are in the "global" environment, \bold{envir}.  From its markers_arg argument, it gets
#'  the locus_index and the index in the \emph{unified_genotype_table}.
#'  From the "global" environment, \bold{envir}, it gets a bit vector of compressed genotype information,
#'  allele information, and some bookkeeping related data.
#'  Note: This function also contains a dispatch/switch on the type of compression in the genotype
#'  vector.  A different C++ function is called when there is compression versus when there is no
#'  compression.
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
#' @param envir an environment that contains all the data frames created from the SQLite database.
#'
#' @return the GenABEL gwaa.data-class object component that contains the genotype data.
#'
#' @export
#' @useDynLib Mega2R
#'
#' @details
#'  This function reads the genotype data in Mega2 compressed format and converts it to the GenABEL
#'  compressed format.
#'  The \emph{unified_genotype_table} contains one raw vector for each person.  In the vector,
#'  there are two bits for each genotype;  each byte has the data for 4 markers.  In GenABEL,
#'  there is one raw vector per marker, and each byte has the data for 4 persons.  The C++
#'  function does the conversion as well as adjusts the bits' contents.  For example, in GenABEL
#'  the genotype represented by bits == 0, is what Mega2 represents with 2.
#'  Doing the conversion in C++ is 10 - 20 times faster than converting the Mega2 data to
#'  PLINK .tped files and then having GenABEL read in and process/convert those files.
#' @note This function is called from \code{Mega2ENVGenABEL}; it is not intended to be called
#'  by the programmer.
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db)
#'
#' aa = getgenotypesgenabel(ENV$markers[ENV$markers$chromosome == 1,])
#'
#' aa
#'
getgenotypesgenabel = function(markers_arg, envir = ENV) {
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)

  return(
    if (envir$MARKER_SCHEME == 1) {
        getgenotypesgenabel_1(markers_arg$locus_link, markers_arg$locus_link_fill,
                              envir$unified_genotype_table, envir$allele_table,
                              envir$PhenoCnt)
    } else {  # must be == 2
        getgenotypesgenabel_2(markers_arg$locus_link,
                              envir$unified_genotype_table, envir$locus_allele_table,
                              envir$PhenoCnt)
    }
   )
}

#' fetch dosage integer matrix for specified markers
#'
#' @description
#'  This function calls a C++ function that does all the heavy lifting.  It passes the arguments
#'  necessary for the C++ function: some from the caller's arguments and some from data frames
#'  that are in the "global" environment, \bold{envir}.  From its markers_arg argument, it gets
#'  the locus_index and the index in the \emph{unified_genotype_table}.
#'  From the "global" environment, \bold{envir}, it gets a bit vector of compressed genotype information,
#'  and some bookkeeping related data.
#'  Note: This function also contains a dispatch/switch on the type of compression in the genotype
#'  vector.  A different C++ function is called when there is compression versus when there is no
#'  compression.
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
#' @param envir an environment that contains all the data frames created from the SQLite database.
#'
#' @return a list of 3 values, named "ncol", "zero", "geno".
#'  \describe{
#'  \item{geno}{is a matrix of dosages as integers.  The value 0 is given to the Major allele 
#'  value, 1 is given to the heterozygote value, and 2 is given to the Minor allele.
#'  In the matrix, there is usually one column for each  marker in the \emph{markers_arg} argument.
#'  But if there would be only the one allele 0 or 2 in the column, the column is ignorednot present.
#'  There is one row for each person in the family (\emph{fam}) table.}
#'  \item{ncol}{Is the count of the actual number of columns in the geno matrix.}
#'  \item{zero}{Is a vector with one entry per marker.  The value will be 0 if the marker is
#'  not in the geno matrix.  Otherwise the value is the column number in the geno matrix where
#'  the marker data appears.}}
#'
#' @export
#' @useDynLib Mega2R
#'
#' @details
#'  The \emph{unified_genotype_table} contains one raw vector for each person.  In the vector,
#'  there are two bits for each genotype.  This function creates an output matrix by fixing
#'  the marker and collecting genotype information for each person and then repeating for
#'  all the specified markers.
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = read.Mega2DB(db)
#'
#' getgenotypesdos(ENV$markers[ENV$markers$chromosome == 1,])
#'
getgenotypesdos = function(markers_arg, envir = ENV) {
    if (missing(envir)) envir = get("ENV", parent.frame(), inherits = TRUE)

  return(
    if (envir$MARKER_SCHEME == 1) {
        getgenotypesdos_1(markers_arg$locus_link, markers_arg$locus_link_fill,
                          envir$unified_genotype_table, envir$allele_table,
                          envir$PhenoCnt)
    } else {  # must be == 2
        getgenotypesdos_2(markers_arg$locus_link,
                          envir$unified_genotype_table, envir$locus_allele_table,
                          envir$PhenoCnt)
    }
  )
}

#' computeDosage function
#'
#' @description
#'  Convert the genotypesraw() allele patterns of 0x10001, 0x10002 (or 0x20001), 0x20002, 0
#"  from the genotype matrix
#'  to the numbers 0, 1, 2, 9 for each marker. (Reverse, the order iff allele "1" has the
#'  minor allele frequency.)  
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
#' @return a matrix of samples X markers for all the markers that have nonzero changes.
#'
#' @export
#'
#' @seealso \code{\link{DOfamSKATRC}}
#'
#' @examples
#' db = system.file("exdata", "seqsimm.db", package="Mega2R")
#' ENV = init_famSKATRC(db, verbose = TRUE)
#' dimDosage = function(m, r, e) {print(dim(computeDosage(m, r, e)))}
#' applyFnToRanges(dimDosage, ENV$refRanges[50:60, ], ENV$refIndices, envir=ENV)
#' # This will use return dosage matrices for the markers in the ranges 50 - 60,
#' # but is basically ignores the results.
#'
#
computeDosage = function (markers_arg, range_arg, envir) {

#   tim0 = system.time ({

    geno_arg = getgenotypesraw(markers_arg, envir);

    markerNames = markers_arg$MarkerName
    gene  = as.character(range_arg[,envir$refCol[4]])
    
    di = dim(geno_arg)
    geno = matrix(0, nrow = (di[1]), ncol = di[2])
    kk = 0
    for (k in 1:(di[2])) {
        vec = envir$dosage[match(as.integer(geno_arg[ , k]), envir$dosageRaw)]
        g0 = sum(vec == 0)
        g1 = sum(vec == 1)
        g2 = sum(vec == 2)
        if (envir$verbose)
            cat(gene, markerNames[k], g0, g1, g2, "\n")

        if (g0 == di[1] || g1 == di[1] || g2 == di[1])
            next

        if (g0 < g2) {
           kk = kk + 1
           geno[ , kk] = 2 - vec
        } else {
           kk = kk + 1
           geno[ , kk] =     vec
        }
    }

    if (kk < di[2])
        geno = geno[,1:kk]

#   })

    geno
}


Rcpp::sourceCpp("src/getgenotypes.cpp")
