# Gene usage analysis

#### Calculation functions ####

#' Tabulates V(D)J allele, gene or family usage.
#' 
#' Determines the count and relative abundance of V(D)J alleles, genes or families within
#' groups.
#'
#' @param    data    data.frame with AIRR-format or Change-O style columns.
#' @param    gene    column containing allele assignments. Only the first allele in the
#'                   column will be considered when \code{mode} is "gene", "family" or 
#'                   "allele". The value will be used as it is with \code{mode="asis"}. 
#' @param    groups  columns containing grouping variables. If \code{NULL} do not group.
#' @param    copy    name of the \code{data} column containing copy numbers for each 
#'                   sequence. If this value is specified, then total copy abundance
#'                   is determined by the sum of copy numbers within each gene.
#'                   This argument is ignored if \code{clone} is specified.
#' @param    clone   name of the \code{data} column containing clone identifiers for each 
#'                   sequence. If this value is specified, then one gene will be considered 
#'                   for each clone. Note, this is accomplished by using the most 
#'                   common gene within each \code{clone} identifier. As such,
#'                   ambiguous alleles within a clone will not be accurately represented.
#' @param    mode    one of \code{c("gene", "family", "allele", "asis")} defining
#'                   the degree of specificity regarding allele calls. Determines whether 
#'                   to return counts for genes (calling \code{getGene}), 
#'                   families (calling \code{getFamily}), alleles (calling 
#'                   \code{getAllele}) or using the value as it is in the column
#'                   \code{gene}, without any processing.
#' @param    fill    logical of \code{c(TRUE, FALSE)} specifying when if groups (when specified)
#'                   lacking a particular gene should be counted as 0 if TRUE or not (omitted) 
#' 
#' @return   A data.frame summarizing family, gene or allele counts and frequencies 
#'           with columns:
#'           \itemize{
#'             \item \code{gene}:         name of the family, gene or allele
#'             \item \code{seq_count}:    total number of sequences for the gene.
#'             \item \code{seq_freq}:     frequency of the gene as a fraction of the total
#'                                        number of sequences within each grouping.
#'             \item \code{copy_count}:   sum of the copy counts in the \code{copy} column.
#'                                        for each gene. Only present if the \code{copy} 
#'                                        argument is specified.
#'             \item \code{copy_freq}:    frequency of the gene as a fraction of the total
#'                                        copy number within each group. Only present if 
#'                                        the \code{copy} argument is specified.
#'             \item \code{clone_count}:  total number of clones for the gene.
#'             \item \code{clone_freq}:   frequency of the gene as a fraction of the total
#'                                        number of clones within each grouping.
#'           }
#'           Additional columns defined by the \code{groups} argument will also be present.
#'
#' @examples
#' # Without copy numbers
#' genes <- countGenes(ExampleDb, gene="v_call", groups="sample_id", mode="family")
#' genes <- countGenes(ExampleDb, gene="v_call", groups="sample_id", mode="gene")
#' genes <- countGenes(ExampleDb, gene="v_call", groups="sample_id", mode="allele")
#'
#' # With copy numbers and multiple groups
#' genes <- countGenes(ExampleDb, gene="v_call", groups=c("sample_id", "c_call"), 
#'                     copy="duplicate_count", mode="family")
#' 
#' # Count by clone
#' genes <- countGenes(ExampleDb, gene="v_call", groups=c("sample_id", "c_call"), 
#'                     clone="clone_id", mode="family")
#'
#' # Count absent genes 
#' genes <- countGenes(ExampleDb, gene="v_call", groups="sample_id", 
#'                     mode="allele", fill=TRUE)
#'
#'@export
countGenes <- function(data, gene, groups=NULL, copy=NULL, clone=NULL, fill=FALSE,
                       mode=c("gene", "allele", "family", "asis")) {
    ## DEBUG
    # data=ExampleDb; gene="c_call"; groups=NULL; mode="gene"; clone="clone_id"
    # data=subset(db, clond_id == 3138)
    # Hack for visibility of dplyr variables
    . <- NULL
    
    # Check input
    mode <- match.arg(mode)
    check <- checkColumns(data, c(gene, groups, copy))
    if (check != TRUE) { stop(check) }

    # Extract gene, allele or family assignments
    if (mode != "asis") {
        gene_func <- switch(mode,
                            allele=getAllele,
                            gene=getGene,
                            family=getFamily)
        data[[gene]] <- gene_func(data[[gene]], first=TRUE)
    }
    
    # Tabulate abundance
    if (is.null(copy) & is.null(clone)) {
        # Tabulate sequence abundance
        gene_tab <- data %>% 
            group_by(!!!rlang::syms(c(groups, gene))) %>%
            dplyr::summarize(seq_count=n()) %>%
            mutate(., seq_freq=!!rlang::sym("seq_count")/sum(!!rlang::sym("seq_count"), na.rm=TRUE)) %>%
            arrange(desc(!!rlang::sym("seq_count")))
    } else if (!is.null(clone) & is.null(copy)) {
        # Find count of genes within each clone and keep first with maximum count
        gene_tab <- data %>%
            group_by(!!!rlang::syms(c(groups, clone, gene))) %>%
            dplyr::mutate(clone_gene_count=n()) %>%
            ungroup() %>%
            group_by(!!!rlang::syms(c(groups, clone))) %>%
            slice(which.max(!!rlang::sym("clone_gene_count"))) %>%
            ungroup() %>%
            group_by(!!!rlang::syms(c(groups, gene))) %>%
            dplyr::summarize(clone_count=n()) %>%
            mutate(clone_freq=!!rlang::sym("clone_count")/sum(!!rlang::sym("clone_count"), na.rm=TRUE)) %>%
            arrange(!!rlang::sym("clone_count"))
    } else {
        if (!is.null(clone) & !is.null(copy)) {
            warning("Specifying both 'copy' and 'clone' columns is not meaningful. ",
                    "The 'clone' argument will be ignored.")
        }
        # Tabulate copy abundance
        gene_tab <- data %>% 
            group_by(!!!rlang::syms(c(groups, gene))) %>%
            summarize(seq_count=length(!!rlang::sym(gene)),
                       copy_count=sum(!!rlang::sym(copy), na.rm=TRUE)) %>%
            mutate(seq_freq=!!rlang::sym("seq_count")/sum(!!rlang::sym("seq_count"), na.rm=TRUE),
                   copy_freq=!!rlang::sym("copy_count")/sum(!!rlang::sym("copy_count"), na.rm=TRUE)) %>%
            arrange(desc(!!rlang::sym("copy_count")))
    }

    # If a gene is present in one GROUP but not another, will fill the COUNT and FREQ with 0s
    if (fill) {
        gene_tab <- gene_tab %>%
            ungroup() %>%
            tidyr::complete(!!!rlang::syms(as.list(c(groups, gene)))) %>%
            replace(is.na(.), 0)
    }

    # Rename gene column
    gene_tab <- rename(gene_tab, "gene"=gene)
    
    return(gene_tab)
}


#### Annotation functions ####

#' Get Ig segment allele, gene and family names
#' 
#' \code{getSegment} performs generic matching of delimited segment calls with a custom 
#' regular expression. \link{getAllele}, \link{getGene} and \link{getFamily} extract 
#' the allele, gene and family names, respectively, from a character vector of 
#' immunoglobulin (Ig) or TCR segment allele calls in IMGT format.
#'
#' @param     segment_call    character vector containing segment calls delimited by commas.
#' @param     segment_regex   string defining the segment match regular expression.
#' @param     first           if \code{TRUE} return only the first call in 
#'                            \code{segment_call}; if \code{FALSE} return all calls 
#'                            delimited by commas.
#' @param     collapse        if \code{TRUE} check for duplicates and return only unique 
#'                            segment assignments; if \code{FALSE} return all assignments 
#'                            (faster). Has no effect if \code{first=TRUE}.
#' @param     strip_d         if \code{TRUE} remove the "D" from the end of gene annotations 
#'                            (denoting a duplicate gene in the locus); 
#'                            if \code{FALSE} do not alter gene names.
#' @param     omit_nl         if \code{TRUE} remove non-localized (NL) genes from the result.
#'                            Only applies at the gene or allele level.
#' @param     sep             character defining both the input and output segment call 
#'                            delimiter.
#'
#' @return    A character vector containing allele, gene or family names.
#' 
#' @references
#'   \url{http://imgt.org}
#'
#' @seealso  \link{countGenes}
#'
#' @examples
#' # Light chain examples
#' kappa_call <- c("Homsap IGKV1D-39*01 F,Homsap IGKV1-39*02 F,Homsap IGKV1-39*01",
#'                 "Homsap IGKJ5*01 F")
#'
#' getAllele(kappa_call)
#' getAllele(kappa_call, first=FALSE)
#' getAllele(kappa_call, first=FALSE, strip_d=FALSE)
#' 
#' getGene(kappa_call)
#' getGene(kappa_call, first=FALSE)
#' getGene(kappa_call, first=FALSE, strip_d=FALSE)
#' 
#' getFamily(kappa_call)
#' getFamily(kappa_call, first=FALSE)
#' getFamily(kappa_call, first=FALSE, collapse=FALSE)
#' getFamily(kappa_call, first=FALSE, strip_d=FALSE)
#' 
#' # Heavy chain examples
#' heavy_call <- c("Homsap IGHV1-69*01 F,Homsap IGHV1-69D*01 F", 
#'                 "Homsap IGHD1-1*01 F", 
#'                 "Homsap IGHJ1*01 F")
#' 
#' getAllele(heavy_call, first=FALSE)
#' getAllele(heavy_call, first=FALSE, strip_d=FALSE)
#'
#' getGene(heavy_call, first=FALSE)
#' getGene(heavy_call, first=FALSE, strip_d=FALSE)
#'
#' # Filtering non-localized genes
#' nl_call <- c("IGHV3-NL1*01,IGHV3-30-3*01,IGHV3-30*01", 
#'              "Homosap IGHV3-30*01 F,Homsap IGHV3-NL1*01 F",
#'              "IGHV1-NL1*01")
#'              
#' getAllele(nl_call, first=FALSE, omit_nl=TRUE)
#' getGene(nl_call, first=FALSE, omit_nl=TRUE)
#' getFamily(nl_call, first=FALSE, omit_nl=TRUE)
#'
#' @export
getSegment <- function(segment_call, segment_regex, first=TRUE, collapse=TRUE, 
                       strip_d=TRUE, omit_nl=FALSE, sep=",") {
    # Define boundaries of individual segment calls
    edge_regex <- paste0("[^", sep, "]*")
    
    # Extract calls
    r <- gsub(paste0(edge_regex, "(", segment_regex, ")", edge_regex), "\\1", 
              segment_call, perl=T)
    
    # Remove NL genes
    if (omit_nl) {
        nl_regex <- paste0('(IG[HLK]|TR[ABGD])[VDJ][0-9]+-NL[0-9]([-/\\w]*[-\\*][\\.\\w]+)*(', 
                           sep, "|$)")
        r <- gsub(nl_regex, "", r, perl=TRUE)
    }
    
    # Strip D from gene names if required
    if (strip_d) {
        strip_regex <- paste0("(?<=[A-Z0-9])D(?=\\*|-|", sep, "|$)")
        r <- gsub(strip_regex, "", r, perl=TRUE)
    }
    
    # Collapse to unique set if required
    if (first) {
        r <- gsub(paste0(sep, ".*$"), "", r)
    } else if (collapse) {
        r <- sapply(strsplit(r, sep), function(x) paste(unique(x), collapse=sep))
    }
    
    return(r)
}


#' @rdname getSegment
#' @export
getAllele <- function(segment_call, first=TRUE, collapse=TRUE, 
                      strip_d=TRUE, omit_nl=FALSE, sep=",") {    
    allele_regex <- '((IG[HLK]|TR[ABGD])[VDJ][A-Z0-9\\(\\)]+[-/\\w]*[-\\*]*[\\.\\w]+)'
    r <- getSegment(segment_call, allele_regex, first=first, collapse=collapse, 
                    strip_d=strip_d, omit_nl=omit_nl, sep=sep)
    
    return(r)
}


#' @rdname getSegment
#' @export
getGene <- function(segment_call, first=TRUE, collapse=TRUE, 
                    strip_d=TRUE, omit_nl=FALSE, sep=",") {
    gene_regex <- '((IG[HLK]|TR[ABGD])[VDJ][A-Z0-9\\(\\)]+[-/\\w]*)'
    r <- getSegment(segment_call, gene_regex, first=first, collapse=collapse, 
                    strip_d=strip_d, omit_nl=omit_nl, sep=sep)
    
    return(r)
}


#' @rdname getSegment
#' @export
getFamily <- function(segment_call, first=TRUE, collapse=TRUE, 
                      strip_d=TRUE, omit_nl=FALSE, sep=",") {
    family_regex <- '((IG[HLK]|TR[ABGD])[VDJ][A-Z0-9\\(\\)]+)'
    r <- getSegment(segment_call, family_regex, first=first, collapse=collapse, 
                    strip_d=strip_d, omit_nl=omit_nl, sep=sep)
    
    return(r)
}


#### Utility functions ####

# Get all VJ(L) combinations from one or more chains of the same type
#
# Input: 
# V annotation, J annotation, and optionally junction length of
# one of more chains of the same type
# - v: V annotation
# - j: J annotation
# - l: junction length (optional)
# - sep_chain: character separting multiple chains
# - sep_anno: character separating multiple/ambiguous annotations within each chain
# - first: to be passed to getGene()
# 
# Output:
# A vector containing all unique VJ(L) combinations represented
#
# Assumption: 
# 1) number of chains match across v, j, l
# 2) if length value is supplied, each chain has a single length value
#
# Example input
# v <- "Homsap IGLV2-20*09 F,Homsap IGLV3-30*09 F;Homsap IGKV1-27*01 F;Homsap IGKV1-25*01 F,Homsap IGKV1-25*38 F,Homsap IGKV1-32*02 F"
# j <- "Homsap IGLJ3*02 F;Homsap IGKJ5*02 F;Homsap IGKJ3*03 F,Homsap IGKJ9*03 F"
# l <- "36;39;60"
# getAllVJL(v, j, l, ";", ",", FALSE)
# - 3 light chains
# - number of V annotations per chain: 2/1/3
# - number of J annotations per chain: 1/1/2
# - Expected supremum of the number of VJL combinations: 2*1 + 1*1 + 3*2 = 9
# - Note that this is the supremum because the ACTUAL number (7) CAN be lower due to presence
#   of different alleles from the SAME gene (since computation is performed at the gene level)
getAllVJL <- function(v, j, l, sep_chain, sep_anno, first) {
    # is l NULL?
    l_NULL <- is.null(l)
    
    # are there multiple chains?
    # assumes that number of chains match across v, j, l
    # (pre-checked in groupGenes)
    multi_chain <- stringi::stri_detect_fixed(str=v, pattern=sep_chain)
    
    # are there multiple annotations per chain?
    multi_anno_v <- stringi::stri_detect_fixed(str=v, pattern=sep_anno)
    multi_anno_j <- stringi::stri_detect_fixed(str=j, pattern=sep_anno)
    
    # separate chains
    # gets a vector of strings
    # each vector entry corresponds to a chain
    if (multi_chain) {
        v <- stringi::stri_split_fixed(str=v, pattern=sep_chain)[[1]]
        j <- stringi::stri_split_fixed(str=j, pattern=sep_chain)[[1]]
        
        if (!l_NULL) { l <- stringi::stri_split_fixed(str=l, pattern=sep_chain)[[1]] }
    }
    
    # separate annotations
    # gets a list
    # each list entry has a vector of one or more strings
    if (multi_anno_v) {
        v <- stringi::stri_split_fixed(str=v, pattern=sep_anno)
        # take care of 'first' here because getGene will not split by ","
        if (first) { v <- lapply(v, function(x){x[1]}) }
    }
    if (multi_anno_j) {
        j <- stringi::stri_split_fixed(str=j, pattern=sep_anno)
        if (first) { j <- lapply(j, function(x){x[1]}) }
    }
    
    # get gene
    # gets a list
    # each list entry corresponds to a chain, and is a vector of one or more strings
    v <- sapply(v, getGene, collapse=TRUE, simplify=FALSE, USE.NAMES=FALSE)
    j <- sapply(j, getGene, collapse=TRUE, simplify=FALSE, USE.NAMES=FALSE)
    
    # if there is multiple chains and/or multiple annotations
    if ( multi_chain | (multi_anno_v & first) | (multi_anno_j & first) ) {
        # gets a list
        # each list entry is a vector of one or more strings
        # do if/else outside sapply so it only gets evaluated once
        if (!l_NULL) {
            exp <- sapply(1:length(v), function(i){
                #eg_df <- expand.grid(v[[i]], j[[i]], l[i])
                #eg_vec <- apply(eg_df, 1, stringi::stri_paste, collapse="@")
                n_v <- length(v[[i]])
                n_j <- length(j[[i]])
                eg_vec = stringi::stri_paste(rep.int(v[[i]], times=n_j),
                                             rep(j[[i]], each=n_v),
                                             rep.int(l[i], times=n_v*n_j),
                                             sep="@")
                return(eg_vec)
            }, simplify=FALSE, USE.NAMES=FALSE)
        } else {
            exp <- sapply(1:length(v), function(i){
                #eg_df = expand.grid(v[[i]], j[[i]])
                #eg_vec = apply(eg_df, 1, stringi::stri_paste, collapse="@")
                n_v <- length(v[[i]])
                n_j <- length(j[[i]])
                eg_vec <- stringi::stri_paste(rep.int(v[[i]], times=n_j),
                                              rep(j[[i]], each=n_v),
                                              sep="@")
                return(eg_vec)
            }, simplify=FALSE, USE.NAMES=FALSE)
        }
        
        # concat and convert to vector; keep distinct values
        exp <- unique(unlist(exp, use.names=FALSE))
    } else {
        
        n_v <- length(v[[1]])
        n_j <- length(j[[1]])
        
        if (!l_NULL) {
            exp <- stringi::stri_paste(rep.int(v[[1]], times=n_j),
                                       rep(j[[1]], each=n_v),
                                       rep.int(l[1], times=n_v*n_j),
                                       sep="@")
        } else {
            exp <- stringi::stri_paste(rep.int(v[[1]], times=n_j),
                                       rep(j[[1]], each=n_v),
                                       sep="@")
        }
        
    }
    
    return(exp)
}


#' Group sequences by gene assignment
#'
#' \code{groupGenes} will group rows by shared V and J gene assignments, 
#' and optionally also by junction lengths. IGH:IGK/IGL, TRB:TRA, and TRD:TRG 
#' paired single-cell BCR/TCR sequencing and unpaired bulk sequencing 
#' (IGH, TRB, TRD chain only) are supported. In the case of ambiguous (multiple) 
#' gene assignments, the grouping may be specified to be a union across all 
#' ambiguous V and J gene pairs, analogous to single-linkage clustering 
#' (i.e., allowing for chaining).
#'
#' @param    data          data.frame containing sequence data.
#' @param    v_call        name of the column containing the heavy/long chain 
#'                         V-segment allele calls.
#' @param    j_call        name of the column containing the heavy/long chain 
#'                         J-segment allele calls.
#' @param    junc_len      name of column containing the junction length.
#'                         If \code{NULL} then 1-stage partitioning is perform
#'                         considering only the V and J genes is performed. 
#'                         See Details for further clarification.
#' @param    cell_id       name of the column containing cell identifiers or barcodes. 
#'                         If specified, grouping will be performed in single-cell mode
#'                         with the behavior governed by the \code{locus} and 
#'                         \code{only_heavy} arguments. If set to \code{NULL} then the 
#'                         bulk sequencing data is assumed.
#' @param    locus         name of the column containing locus information. 
#'                         Only applicable to single-cell data.
#'                         Ignored if \code{cell_id=NULL}.
#' @param    only_heavy    use only the IGH (BCR or TRB/TRD (TCR) sequences 
#'                         for grouping. Only applicable to single-cell data.
#'                         Ignored if \code{cell_id=NULL}.
#' @param    first         if \code{TRUE} only the first call of the gene assignments 
#'                         is used. if \code{FALSE} the union of ambiguous gene 
#'                         assignments is used to group all sequences with any 
#'                         overlapping gene calls.
#'
#' @return   Returns a modified data.frame with disjoint union indices 
#'           in a new \code{vj_group} column. 
#'           
#'           If \code{junc_len} is supplied, the grouping this \code{vj_group} 
#'           will have been based on V, J, and junction length simultaneously. However, 
#'           the output column name will remain \code{vj_group}.
#'           
#'           The output \code{v_call}, \code{j_call}, \code{cell_id}, and \code{locus}
#'           columns will be converted to type \code{character} if they were of type 
#'           \code{factor} in the input \code{data}.
#'
#' @details
#' To invoke single-cell mode the \code{cell_id} argument must be specified and the \code{locus} 
#' column must be correct. Otherwise, \code{groupGenes} will be run with bulk sequencing assumptions, 
#' using all input sequences regardless of the values in the \code{locus} column.
#' 
#' Values in the \code{locus} column must be one of \code{c("IGH", "IGI", "IGK", "IGL")} for BCR 
#' or \code{c("TRA", "TRB", "TRD", "TRG")} for TCR sequences. Otherwise, the function returns an 
#' error message and stops.
#' 
#' Under single-cell mode with paired chained sequences, there is a choice of whether 
#' grouping should be done by (a) using IGH (BCR) or TRB/TRD (TCR) sequences only or
#' (b) using IGH plus IGK/IGL (BCR) or TRB/TRD plus TRA/TRG (TCR). 
#' This is governed by the \code{only_heavy} argument.
#' 
#' Specifying \code{junc_len} will force \code{groupGenes} to perform a 1-stage partitioning of the 
#' sequences/cells based on V gene, J gene, and junction length simultaneously. 
#' If \code{junc_len=NULL} (no column specified), then \code{groupGenes} performs only the first 
#' stage of a 2-stage partitioning in which sequences/cells are partitioned in the first stage 
#' based on V gene and J gene, and then in the second stage further splits the groups based on 
#' junction length (the second stage must be performed independently, as this only returns the
#' first stage results).
#' 
#' In the input \code{data}, the \code{v_call}, \code{j_call}, \code{cell_id}, and \code{locus} 
#' columns, if present, must be of type \code{character} (as opposed to \code{factor}). 
#' 
#' It is assumed that ambiguous gene assignments are separated by commas.
#' 
#' All rows containing \code{NA} values in any of the \code{v_call}, \code{j_call}, and \code{junc_len} 
#' (if \code{junc_len != NULL}) columns will be removed. A warning will be issued when a row 
#' containing an \code{NA} is removed.
#' 
#' @section Expectations for single-cell data:
#' 
#' Single-cell paired chain data assumptions:
#'   \itemize{
#'      \item every row represents a sequence (chain).
#'      \item heavy/long and light/short chains of the same cell are linked by \code{cell_id}.
#'      \item the value in \code{locus} column indicates whether the chain is the heavy/long or light/short chain.
#'      \item each cell possibly contains multiple heavy/long and/or light/short chains.
#'      \item every chain has its own V(D)J annotation, in which ambiguous V(D)J 
#'            annotations, if any, are separated by a comma.
#'   }
#'   
#' Single-cell example:
#'   \itemize{
#'      \item A cell has 1 heavy chain and 2 light chains.
#'      \item There should be 3 rows corresponding to this cell.
#'      \item One of the light chains may have an ambiguous V annotation which looks like \code{"Homsap IGKV1-39*01 F,Homsap IGKV1D-39*01 F"}.
#'   }
#' 
#' @examples
#' # Group by genes
#' db <- groupGenes(ExampleDb)
#' head(db$vj_group)
#'  
#' @export
groupGenes <- function(data, v_call="v_call", j_call="j_call", junc_len=NULL,
                       cell_id=NULL, locus="locus", only_heavy=TRUE,
                       first=FALSE) {
    # Check base input
    check <- checkColumns(data, c(v_call, j_call, junc_len))
    if (check != TRUE) { stop(check) }
    
    # Check single-cell input
    if (!is.null(cell_id)) {
        check <- checkColumns(data, c(cell_id, locus))
        if (check != TRUE) { stop(check) }
    }
    
    # if necessary, cast select columns to character (factor not allowed later on)
    if (!is(data[[v_call]], "character")) { data[[v_call]] <- as.character(data[[v_call]]) }
    if (!is(data[[j_call]], "character")) { data[[j_call]] <- as.character(data[[j_call]]) }
    
    # e.g.: "Homsap IGHV3-7*01 F,Homsap IGHV3-6*01 F;Homsap IGHV1-4*01 F"
    separator_within_seq <- ","
    separator_between_seq <- ";"
    
    # single-cell mode?
    if (!is.null(cell_id) & !is.null(locus)) {
        single_cell <- TRUE
        
        if (!is(data[[cell_id]], "character")) { data[[cell_id]] <- as.character(data[[cell_id]]) }
        if (!is(data[[locus]], "character")) { data[[locus]] <- as.character(data[[locus]]) }
        
        # check locus column
        valid_loci <- c("IGH", "IGI", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG")
        check <- !all(unique(data[[locus]]) %in% valid_loci)
        if (check) {
            stop("The locus column contains invalid loci annotations.")
        }
    } else {
        single_cell <- FALSE
    }
    
    # only set if `single_cell` & `only_heavy`
    v_call_light <- NULL
    j_call_light <- NULL
    junc_len_light <- NULL
    
    # single-cell mode
    if (single_cell) {
        
        # regardless of using heavy only, or using both heavy and light
        # for each cell
        # - index wrt data of heavy chain
        # - index wrt data of light chain(s)
        cell_id_uniq <- unique(data[[cell_id]])
        cell_seq_idx <- sapply(cell_id_uniq, function(x){
            # heavy chain
            idx_h <- which( data[[cell_id]]==x & data[[locus]] %in% c("IGH", "TRB", "TRD")) 
            # light chain
            idx_l <- which( data[[cell_id]]==x & data[[locus]] %in% c("IGK", "IGL", "TRA", "TRG") )
            
            return(list(heavy=idx_h, light=idx_l))
        }, USE.NAMES=FALSE, simplify=FALSE)
        
        # make a copy
        data_orig <- data; rm(data)
        
        if (only_heavy) {
            
            # use heavy chains only
            
            # Straightforward subsetting like below won't work in cases 
            #     where multiple HCs are present for a cell 
            # subset to heavy only
            # data <- data_orig[data_orig[[locus]]=="IGH", ]
            
            # flatten data
            cols <- c(cell_id, v_call, j_call, junc_len)
            data <- data.frame(matrix(NA, nrow=length(cell_seq_idx), ncol=length(cols)))
            colnames(data) <- cols
            
            for (i_cell in 1:length(cell_seq_idx)) {
                i_cell_h <- cell_seq_idx[[i_cell]][["heavy"]]
                
                data[[cell_id]][i_cell] <- cell_id_uniq[i_cell]
                
                # heavy chain V, J, junc_len
                data[[v_call]][i_cell] <- paste0(data_orig[[v_call]][i_cell_h], 
                                               collapse=separator_between_seq)
                data[[j_call]][i_cell] <- paste0(data_orig[[j_call]][i_cell_h], 
                                               collapse=separator_between_seq)
                if (!is.null(junc_len)) {
                    data[[junc_len]][i_cell] <- paste0(data_orig[[junc_len]][i_cell_h], 
                                                     collapse=separator_between_seq)
                }
            }
            
            
        } else {
            
            # use heavy AND light chains for grouping
            v_call_light <- "v_call_light"
            j_call_light <- "j_call_light"
            # ifelse won't return NULL
            if (is.null(junc_len)) {
                junc_len_light <- NULL
            } else {
                junc_len_light <- "len_light"
            }
            
            # flatten data
            cols <- c(cell_id, v_call, j_call, junc_len, v_call_light, j_call_light, junc_len_light)
            data <- data.frame(matrix(NA, nrow=length(cell_seq_idx), ncol=length(cols)))
            colnames(data) <- cols
            
            for (i_cell in 1:length(cell_seq_idx)) {
                i_cell_h <- cell_seq_idx[[i_cell]][["heavy"]]
                i_cell_l <- cell_seq_idx[[i_cell]][["light"]]
                
                data[[cell_id]][i_cell] <- cell_id_uniq[i_cell]
                
                # heavy chain V, J, junc_len
                data[[v_call]][i_cell] <- paste0(data_orig[[v_call]][i_cell_h], 
                                               collapse=separator_between_seq)
                data[[j_call]][i_cell] <- paste0(data_orig[[j_call]][i_cell_h], 
                                               collapse=separator_between_seq)
                if (!is.null(junc_len)) {
                    data[[junc_len]][i_cell] <- paste0(data_orig[[junc_len]][i_cell_h], 
                                                     collapse=separator_between_seq)
                }
                
                # light chain V, J, junc_len
                data[[v_call_light]][i_cell] <- paste0(data_orig[[v_call]][i_cell_l], 
                                                     collapse=separator_between_seq)
                data[[j_call_light]][i_cell] <- paste0(data_orig[[j_call]][i_cell_l], 
                                                     collapse=separator_between_seq)
                if (!is.null(junc_len_light)) {
                    data[[junc_len_light]][i_cell] <- paste0(data_orig[[junc_len]][i_cell_l], 
                                                           collapse=separator_between_seq)
                }
                
            }
            
            # It cannot be the case that there are 2 V annotations but only 1 J annotation for 
            # the two light chains. Both J annotations must be spelled out for each light chain,
            # separated by \code{separator_between_seq}, even if the annotated alleles are the same.
            # 
            # This one-to-one annotation-to-chain correspondence for both V and J is explicitly
            # checked and an error raised if the requirement is not met. 
            
            # one-to-one annotation-to-chain correspondence for both V and J (light)
            # for each cell/row, number of bewteen_seq separators in light V annotation and in light J annotation must match
            n_separator_btw_seq_v_light <- stringi::stri_count_fixed(str=data[[v_call_light]], pattern=separator_between_seq)
            n_separator_btw_seq_j_light <- stringi::stri_count_fixed(str=data[[j_call_light]], pattern=separator_between_seq)
            if (any( n_separator_btw_seq_v_light != n_separator_btw_seq_j_light )) {
                stop("Requirement not met: one-to-one annotation-to-chain correspondence for both V and J (light)")
            }
            
        }
    } 
    
    # one-to-one annotation-to-chain correspondence for both V and J (heavy)
    # for each cell/row, number of bewteen_seq separators in heavy V annotation and in heavy J annotation must match
    # (in theory, there should be 1 heavy chain per cell; but 10x can return cell with >1 heavy chains and 
    #  you never know if the user will supply this cell as input)
    n_separator_btw_seq_v_heavy <- stringi::stri_count_fixed(str=data[[v_call]], pattern=separator_between_seq)
    n_separator_btw_seq_j_heavy <- stringi::stri_count_fixed(str=data[[j_call]], pattern=separator_between_seq)
    if (any( n_separator_btw_seq_v_heavy != n_separator_btw_seq_j_heavy )) {
        stop("Requirement not met: one-to-one annotation-to-chain correspondence for both V and J (heavy)")
    }
    
    
    # NULL will disappear when doing c()
    # c(NULL,NULL) gives NULL still
    cols_for_grouping_heavy <- c(v_call, j_call, junc_len)
    cols_for_grouping_light <- c(v_call_light, j_call_light, junc_len_light)
    
    # cols cannot be factor
    if (any( sapply(cols_for_grouping_heavy, function(x){class(data[[x]]) == "factor"}) )) {
        stop("one or more of { ", v_call, ", ", j_call,  
             ifelse(is.null(junc_len), " ", ", "), junc_len, 
             "} is factor. Must be character.\nIf using read.table(), make sure to set stringsAsFactors=FALSE.\n")
    }
    if (single_cell & !only_heavy) {
        if (any( sapply(cols_for_grouping_light, function(x) {class(data[[x]]) == "factor"}) )) {
            stop("one or more of { ", v_call_light, ", ", j_call_light,  
                 ifelse(is.null(junc_len_light), " ", ", "), junc_len_light, 
                 "} is factor. Must be character.\nIf using read.table(), make sure to set stringsAsFactors=FALSE.\n")
        }  
    }
    
    # Check NA(s) in columns
    bool_na <- rowSums( is.na( data[, c(cols_for_grouping_heavy, cols_for_grouping_light)] ) ) >0
    if (any(bool_na)) {
        entityName <- ifelse(single_cell, " cell(s)", " sequence(s)")
        msg <- paste0("NA(s) found in one or more of { ", 
                      v_call, ", ", j_call, 
                      ifelse(is.null(junc_len), "", ", "), junc_len,
                      ifelse(is.null(v_call_light), "", ", "), v_call_light,
                      ifelse(is.null(j_call_light), "", ", "), j_call_light,
                      ifelse(is.null(junc_len_light), "", ", "), junc_len_light,
                      " } columns. ", sum(bool_na), entityName, " removed.\n")
        warning(msg)
        data <- data[!bool_na, ]
        if (single_cell) {
            # maintain one-to-one relationship between 
            # rows of data, cell_id_uniq, and cell_seq_idx
            cell_id_uniq <- cell_id_uniq[!bool_na]
            cell_seq_idx <- cell_seq_idx[!bool_na]
        }
    }
    
    ### expand
    
    # speed-up strategy
    # compute expanded VJL combos for unique rows
    # then distribute back to all rows
    
    # unique combinations of VJL
    # heavy chain seqs only
    if ( (!single_cell) | (single_cell & only_heavy) ) {
        combo_unique <- unique(data[, cols_for_grouping_heavy])
        
        # unique components
        v_unique <- unique(combo_unique[[v_call]])
        j_unique <- unique(combo_unique[[j_call]])
        
        # map each row in full data to unique combo
        m_v <- match(data[[v_call]], v_unique)
        m_j <- match(data[[j_call]], j_unique)
        
        if (is.null(junc_len)) {
            combo_unique_full_idx <- sapply(1:nrow(combo_unique), function(i) {
                idx_v <- which (v_unique == combo_unique[[v_call]][i])
                idx_j <- which (j_unique == combo_unique[[j_call]][i])
                idx <- which(m_v==idx_v & m_j==idx_j)
                return(idx)
            }, simplify=FALSE, USE.NAMES=FALSE) 
        } else {
            l_unique <- unique(combo_unique[[junc_len]])
            m_l <- match(data[[junc_len]], l_unique)
            combo_unique_full_idx <- sapply(1:nrow(combo_unique), function(i) {
                idx_v <- which(v_unique == combo_unique[[v_call]][i])
                idx_j <- which(j_unique == combo_unique[[j_call]][i])
                idx_l <- which(l_unique == combo_unique[[junc_len]][i])
                idx <- which(m_v==idx_v & m_j==idx_j & m_l==idx_l)
                return(idx)
            }, simplify=FALSE, USE.NAMES=FALSE)
        }
        
        # expand combo_unique
        if (is.null(junc_len)) {
            exp_lst <- sapply(1:nrow(combo_unique), function(i){
                getAllVJL(v=combo_unique[[v_call]][i], j=combo_unique[[j_call]][i], 
                          l=NULL, first=first,
                          sep_anno=separator_within_seq, sep_chain=separator_between_seq)
            }, simplify=F, USE.NAMES=FALSE)
        } else {
            exp_lst <- sapply(1:nrow(combo_unique), function(i){
                getAllVJL(v=combo_unique[[v_call]][i], j=combo_unique[[j_call]][i], 
                          l=combo_unique[[junc_len]][i], first=first,
                          sep_anno=separator_within_seq, sep_chain=separator_between_seq)
            }, simplify=F, USE.NAMES=FALSE)
        }
        
    } else {
        # single_cell & !only_heavy
        
        # important: do not do this separately for heavy and light
        # must keep the pairing structure
        combo_unique <- unique(data[, c(cols_for_grouping_heavy, cols_for_grouping_light)])
        
        # unique components
        v_unique_h <- unique(combo_unique[[v_call]])
        j_unique_h <- unique(combo_unique[[j_call]])
        v_unique_l <- unique(combo_unique[[v_call_light]])
        j_unique_l <- unique(combo_unique[[j_call_light]])
        
        # map each row in full data to unique combo
        m_v_h <- match(data[[v_call]], v_unique_h)
        m_j_h <- match(data[[j_call]], j_unique_h)
        m_v_l <- match(data[[v_call_light]], v_unique_l)
        m_j_l <- match(data[[j_call_light]], j_unique_l)
        
        if (!is.null(junc_len)) {
            l_unique_h <- unique(combo_unique[[junc_len]])
            m_l_h <- match(data[[junc_len]], l_unique_h)
        }
        if (!is.null(junc_len_light)) {
            l_unique_l <- unique(combo_unique[[junc_len_light]])
            m_l_l <- match(data[[junc_len_light]], l_unique_l)
        }
        
        # expand combo_unique
        if (is.null(junc_len) & is.null(junc_len_light)) {
            
            # map
            combo_unique_full_idx <- sapply(1:nrow(combo_unique), function(i){
                idx_v_h <- which( v_unique_h == combo_unique[[v_call]][i] )
                idx_j_h <- which( j_unique_h == combo_unique[[j_call]][i] )
                idx_v_l <- which( v_unique_l == combo_unique[[v_call_light]][i] )
                idx_j_l <- which( j_unique_l == combo_unique[[j_call_light]][i] )
                idx <- which(m_v_h==idx_v_h & m_j_h==idx_j_h & m_v_l==idx_v_l & m_j_l==idx_j_l)
                return(idx)
            }, simplify=FALSE, USE.NAMES=FALSE) 
            
            # heavy
            exp_h <- sapply(1:nrow(combo_unique), function(i){
                getAllVJL(v=combo_unique[[v_call]][i], j=combo_unique[[j_call]][i], 
                          l=NULL, first=first,
                          sep_anno=separator_within_seq, sep_chain=separator_between_seq)
            }, simplify=FALSE, USE.NAMES=FALSE)
            # light
            exp_l <- sapply(1:nrow(combo_unique), function(i){
                getAllVJL(v=combo_unique[[v_call_light]][i], j=combo_unique[[j_call_light]][i], 
                          l=NULL, first=first,
                          sep_anno=separator_within_seq, sep_chain=separator_between_seq)
            }, simplify=FALSE, USE.NAMES=FALSE)
            
            
        } else if (!is.null(junc_len) & !is.null(junc_len_light)) {
            
            # map
            combo_unique_full_idx <- sapply(1:nrow(combo_unique), function(i){
                idx_v_h <- which( v_unique_h == combo_unique[[v_call]][i] )
                idx_j_h <- which( j_unique_h == combo_unique[[j_call]][i] )
                idx_l_h <- which( l_unique_h == combo_unique[[junc_len]][i] )
                idx_v_l <- which( v_unique_l == combo_unique[[v_call_light]][i] )
                idx_j_l <- which( j_unique_l == combo_unique[[j_call_light]][i] )
                idx_l_l <- which( l_unique_l == combo_unique[[junc_len_light]][i] )
                idx <- which(m_v_h==idx_v_h & m_j_h==idx_j_h & m_l_h==idx_l_h & 
                                m_v_l==idx_v_l & m_j_l==idx_j_l & m_l_l==idx_l_l)
                return(idx)
            }, simplify=FALSE, USE.NAMES=FALSE) 
            
            # heavy
            exp_h <- sapply(1:nrow(combo_unique), function(i){
                getAllVJL(v=combo_unique[[v_call]][i], j=combo_unique[[j_call]][i], 
                          l=combo_unique[[junc_len]][i], first=first,
                          sep_anno=separator_within_seq, sep_chain=separator_between_seq)
            }, simplify=FALSE, USE.NAMES=FALSE)
            # light
            exp_l <- sapply(1:nrow(combo_unique), function(i){
                getAllVJL(v=combo_unique[[v_call_light]][i], j=combo_unique[[j_call_light]][i], 
                          l=combo_unique[[junc_len_light]][i], first=first,
                          sep_anno=separator_within_seq, sep_chain=separator_between_seq)
            }, simplify=FALSE, USE.NAMES=FALSE)
            
        } else if (is.null(junc_len) & !is.null(junc_len_light)) {
            
            # map
            combo_unique_full_idx <- sapply(1:nrow(combo_unique), function(i){
                idx_v_h <- which( v_unique_h == combo_unique[[v_call]][i] )
                idx_j_h <- which( j_unique_h == combo_unique[[j_call]][i] )
                idx_v_l <- which( v_unique_l == combo_unique[[v_call_light]][i] )
                idx_j_l <- which( j_unique_l == combo_unique[[j_call_light]][i] )
                idx_l_l <- which( l_unique_l == combo_unique[[junc_len_light]][i] )
                idx <- which(m_v_h==idx_v_h & m_j_h==idx_j_h & 
                                m_v_l==idx_v_l & m_j_l==idx_j_l & m_l_l==idx_l_l)
                return(idx)
            }, simplify=FALSE, USE.NAMES=FALSE) 
            
            # heavy
            exp_h <- sapply(1:nrow(combo_unique), function(i){
                getAllVJL(v=combo_unique[[v_call]][i], j=combo_unique[[j_call]][i], 
                          l=NULL, first=first,
                          sep_anno=separator_within_seq, sep_chain=separator_between_seq)
            }, simplify=FALSE, USE.NAMES=FALSE)
            # light
            exp_l <- sapply(1:nrow(combo_unique), function(i){
                getAllVJL(v=combo_unique[[v_call_light]][i], j=combo_unique[[j_call_light]][i], 
                          l=combo_unique[[junc_len_light]][i], first=first,
                          sep_anno=separator_within_seq, sep_chain=separator_between_seq)
            }, simplify=FALSE, USE.NAMES=FALSE)
            
        } else if (!is.null(junc_len) & is.null(junc_len_light)) {
            
            # map
            combo_unique_full_idx <- sapply(1:nrow(combo_unique), function(i){
                idx_v_h <- which( v_unique_h == combo_unique[[v_call]][i] )
                idx_j_h <- which( j_unique_h == combo_unique[[j_call]][i] )
                idx_l_h <- which( l_unique_h == combo_unique[[junc_len]][i] )
                idx_v_l <- which( v_unique_l == combo_unique[[v_call_light]][i] )
                idx_j_l <- which( j_unique_l == combo_unique[[j_call_light]][i] )
                idx <- which(m_v_h==idx_v_h & m_j_h==idx_j_h & m_l_h==idx_l_h & 
                                m_v_l==idx_v_l & m_j_l==idx_j_l)
                return(idx)
            }, simplify=FALSE, USE.NAMES=FALSE) 
            
            # heavy
            exp_h <- sapply(1:nrow(combo_unique), function(i){
                getAllVJL(v=combo_unique[[v_call]][i], j=combo_unique[[j_call]][i], 
                          l=combo_unique[[junc_len]][i], first=first,
                          sep_anno=separator_within_seq, sep_chain=separator_between_seq)
            }, simplify=FALSE, USE.NAMES=FALSE)
            # light
            exp_l <- sapply(1:nrow(combo_unique), function(i){
                getAllVJL(v=combo_unique[[v_call_light]][i], j=combo_unique[[j_call_light]][i], 
                          l=NULL, first=first,
                          sep_anno=separator_within_seq, sep_chain=separator_between_seq)
            }, simplify=FALSE, USE.NAMES=FALSE)
            
        }
        
        # pair heavy & light
        stopifnot( length(exp_h) == length(exp_l) )
        exp_lst <- sapply(1:length(exp_h), function(i){
            n_h <- length(exp_h[[i]])
            n_l <- length(exp_l[[i]])
            stringi::stri_paste(rep.int(exp_h[[i]], times=n_l),
                                rep(exp_l[[i]], each=n_h),
                                sep="@")
        }, simplify=FALSE, USE.NAMES=FALSE)
        
    }
    
    # one-to-one correspondence btw exp_lst and combo_unique_full_idx
    # exp_lst: VJL combinations
    # combo_unique_full_idx: rows in data carrying each exp_lst
    # exp_lst may not be all unique because gene-level info is kept instead of allele-level
    # make exp_lst unique
    
    exp_lst_uniq <- unique(exp_lst)
    exp_lst_uniq_full_idx <- sapply(exp_lst_uniq, function(x){
        # wrt exp_lst, therefore also wrt combo_unique_full_idx
        idx_lst <- which(unlist(lapply(exp_lst, function(y){ 
            length(y)==length(x) && all(y==x) 
        })))
        # merge
        
        unlist(combo_unique_full_idx[idx_lst], use.names=FALSE)
    }, simplify=FALSE, USE.NAMES=FALSE)
    
    stopifnot( length(unique(unlist(exp_lst_uniq_full_idx, use.names=FALSE))) == nrow(data) )
    
    # tip: unlist with use.names=F makes it much faster (>100x)
    # https://www.r-bloggers.com/speed-trick-unlist-use-namesfalse-is-heaps-faster/
    exp_uniq <- sort(unique(unlist(exp_lst_uniq, use.names=FALSE)))
    n_cells_or_seqs <- nrow(data)
    
    # notes on implementation
    
    # regular/dense matrix is more straightforward to implement but very costly memory-wise
    # sparse matrix is less straightforward to implement but way more memory efficient
    
    # sparse matrix is very slow to modify to on-the-fly (using a loop like for dense matrix)
    # way faster to construct in one go
    
    # (DO NOT DELETE)
    # for illustrating the concept 
    # this is the way to go if using regular matrix (memory-intensive)
    # same concept implemented using sparse matrix
    
    # mtx_cell_VJL <- matrix(0, nrow=nrow(data), ncol=length(exp_uniq))
    # colnames(mtx_cell_VJL) <- exp_uniq
    # 
    # mtx_adj <- matrix(0, nrow=length(exp_uniq), ncol=length(exp_uniq))
    # rownames(mtx_adj) <- exp_uniq
    # colnames(mtx_adj) <- exp_uniq
    # 
    # outdated:
    # for (i_cell in 1:length(exp_lst)) { 
    #     #if (i_cell %% 1000 == 0) { cat(i_cell, "\n") }
    #     cur_uniq <- unique(exp_lst[[i_cell]])
    #     mtx_cell_VJL[i_cell, cur_uniq] <- 1
    #     mtx_adj[cur_uniq, cur_uniq] <- 1
    # }
    
    # actual implementation using sparse matrix from Matrix package
    
    ### matrix indicating relationship between cell and VJ(L) combinations
    # row: cell
    # col: unique heavy VJ(L) (and light VJ(L))
    
    # row indices
    m1_i <- lapply(1:length(exp_lst_uniq), function(i){
        rep(exp_lst_uniq_full_idx[[i]], each=length(exp_lst_uniq[[i]]))
    })
    m1_i_v <- unlist(m1_i, use.names=FALSE)
    
    # column indices
    m1_j <- lapply(1:length(exp_lst_uniq), function(i){
        # wrt exp_uniq
        idx <- match(exp_lst_uniq[[i]], exp_uniq)
        #stopifnot( all.equal( exp_uniq[idx], exp_lst_uniq[[i]] ) )
        
        rep.int(idx, length(exp_lst_uniq_full_idx[[i]]))
    })
    m1_j_v <- unlist(m1_j, use.names=FALSE)
    
    stopifnot( length(m1_i_v) == length(m1_j_v) )
    
    # no particular need for this to be not of class "nsparseMatrix"
    # so no need to specify x=rep(1, length(m1_i))
    # not specifying makes it even more space-efficient
    mtx_cell_VJL <- Matrix::sparseMatrix(i=m1_i_v, j=m1_j_v, 
                                         dims=c(n_cells_or_seqs, length(exp_uniq)), 
                                         symmetric=F, triangular=F, index1=T, 
                                         dimnames=list(NULL, exp_uniq))
    
    ### adjacency matrix
    # row and col: unique heavy VJ(L) (and light VJ(L))
    
    # row indices
    m2_i <- lapply(1:length(exp_lst_uniq), function(i){
        # wrt exp_uniq
        idx <- match(exp_lst_uniq[[i]], exp_uniq)
        #stopifnot( all.equal( exp_uniq[idx], exp_lst_uniq[[i]] ) )
        
        rep(idx, each=length(exp_lst_uniq[[i]]))
    })
    m2_i_v <- unlist(m2_i, use.names=FALSE)
    
    # col indices
    m2_j <- lapply(1:length(exp_lst_uniq), function(i){
        # wrt exp_uniq
        idx <- match(exp_lst_uniq[[i]], exp_uniq)
        #stopifnot( all.equal( exp_uniq[idx], exp_lst_uniq[[i]] ) )
        
        rep.int(idx, length(exp_lst_uniq[[i]]))
    })
    m2_j_v <- unlist(m2_j, use.names=FALSE)
    
    stopifnot( length(m2_i_v) == length(m2_j_v) )
    
    # important: x must be specified for mtx_adj in order to make it not of class "nsparseMatrix"
    # this is because igraph accepts sparse matrix from Matrix but not the "pattern" matrices variant
    mtx_adj <- Matrix::sparseMatrix(i=m2_i_v, j=m2_j_v, x=rep(1,length(m2_i_v)), 
                                    dims=c(length(exp_uniq), length(exp_uniq)), 
                                    symmetric=F, triangular=F, index1=T, 
                                    dimnames=list(exp_uniq, exp_uniq))
    
    rm(m1_i, m1_j, m2_i, m2_j, m1_i_v, m1_j_v, m2_i_v, m2_j_v, exp_lst)
    
    ### identify connected components based on adjcencey matrix
    # this is the grouping
    # source: https://stackoverflow.com/questions/35772846/obtaining-connected-components-in-r
    
    g <- igraph::graph_from_adjacency_matrix(adjmatrix=mtx_adj, mode="undirected", diag=FALSE)
    #plot(g, vertex.size=10, vertex.label.cex=1, vertex.color="skyblue", vertex.label.color="black", vertex.frame.color="transparent", edge.arrow.mode=0)
    
    connected <- igraph::components(g)
    VJL_groups <- igraph::groups(connected)
    names(VJL_groups) <- paste0("G", 1:length(VJL_groups))
    
    ### identify cells associated with each connected component (grouping)
    
    # each entry corresponds to a group/partition
    # each element within an entry is a cell
    
    cellIdx_byGroup_lst <- lapply(VJL_groups, function(x){ 
        if (length(x)>1) {
            # matrix
            # important to specify rowSums from Matrix package
            # base::rowSums will NOT work
            cell_idx <- which(Matrix::rowSums(mtx_cell_VJL[, x])>0)
        } else {
            # vector
            cell_idx <- which(mtx_cell_VJL[, x]>0)
        }
        return(cell_idx)
    })
    
    # sanity check: there should be perfect/disjoint partitioning 
    # (each cell has exactly one group assignment)
    stopifnot( n_cells_or_seqs == length(unique(unlist(cellIdx_byGroup_lst, use.names=FALSE))) )
    
    # assign
    data$vj_group <- NA
    for (i in 1:length(cellIdx_byGroup_lst)) {
        data[["vj_group"]][cellIdx_byGroup_lst[[i]]] <- names(VJL_groups)[i]
    }
    stopifnot(!any(is.na(data[["vj_group"]])))
    
    if (!single_cell) {
        return(data)
    } else {
        data_orig$vj_group <- NA
        
        # map back to data_orig
        for (i_cell in 1:nrow(data)) {
            # wrt data_orig
            i_orig_h <- cell_seq_idx[[i_cell]][["heavy"]]
            i_orig_l <- cell_seq_idx[[i_cell]][["light"]]
            # sanity check
            stopifnot( all( data_orig[[cell_id]][c(i_orig_h, i_orig_l)] == cell_id_uniq[i_cell] ) )
            # grouping
            data_orig$vj_group[c(i_orig_h, i_orig_l)] <- data$vj_group[i_cell]
        }
        
        # remove rows with $vj_group values of NA
        # these had already been removed by the NA check for `data`
        data_orig <- data_orig[!is.na(data_orig$vj_group), ]
        
        return(data_orig)
    }
    
}


#' Sort V(D)J genes
#'
#' \code{sortGenes} sorts a vector of V(D)J gene names by either lexicographic ordering 
#' or locus position. 
#' 
#' @param    genes    vector of strings respresenting V(D)J gene names.
#' @param    method   string defining the method to use for sorting genes. One of:
#'                    \itemize{
#'                      \item \code{"name"}:      sort in lexicographic order. Order is by 
#'                                                family first, then gene, and then allele. 
#'                      \item \code{"position"}:  sort by position in the locus, as
#'                                                determined by the final two numbers 
#'                                                in the gene name. Non-localized genes 
#'                                                are assigned to the highest positions.
#'                    }
#'                    
#' @return   A sorted character vector of gene names.
#' 
#' @seealso  See \code{getAllele}, \code{getGene} and \code{getFamily} for parsing
#'           gene names.
#' 
#' @examples
#' # Create a list of allele names
#' genes <- c("IGHV1-69D*01","IGHV1-69*01","IGHV4-38-2*01","IGHV1-69-2*01",
#'            "IGHV2-5*01","IGHV1-NL1*01", "IGHV1-2*01,IGHV1-2*05", 
#'            "IGHV1-2", "IGHV1-2*02", "IGHV1-69*02")
#' 
#' # Sort genes by name
#' sortGenes(genes)
#' 
#' # Sort genes by position in the locus
#' sortGenes(genes, method="pos")
#' 
#' @export
sortGenes <- function(genes, method=c("name", "position")) { 
    ## DEBUG
    # method="name"
    
    # Check arguments
    method <- match.arg(method)

    # Build sorting table
    sort_tab <- tibble(CALL=sort(getAllele(genes, first=FALSE, strip_d=FALSE))) %>%
        # Determine the gene and family
        mutate(FAMILY=getFamily(!!rlang::sym("CALL"), first=TRUE, strip_d=FALSE),
               GENE=getGene(!!rlang::sym("CALL"), first=TRUE, strip_d=FALSE),
               ALLELE=getAllele(!!rlang::sym("CALL"), first=TRUE, strip_d=FALSE)) %>%
        # Identify first gene number, second gene number and allele number
        mutate(G1=gsub("[^-]+-([^-\\*D]+).*", "\\1", !!rlang::sym("GENE")),
               G1=as.numeric(gsub("[^0-9]+", "99", !!rlang::sym("G1"))),
               G2=gsub("[^-]+-[^-]+-?", "", !!rlang::sym("GENE")),
               G2=as.numeric(gsub("[^0-9]+", "99", !!rlang::sym("G2"))),
               A1=as.numeric(sub("[^\\*]+\\*|[^\\*]+$", "", !!rlang::sym("ALLELE")))
        )

    # Convert missing values to 0
    sort_tab[is.na(sort_tab)] <- 0
    
    # Sort
    if (method == "name") {  
        sorted_genes <- arrange(sort_tab, !!!rlang::syms(c("FAMILY", "G1", "G2", "A1")))[["CALL"]]
    } else if (method == "position") {
        sorted_genes <- arrange(sort_tab, 
                                desc(!!rlang::sym("G1")), 
                                desc(!!rlang::sym("G2")), 
                                !!rlang::sym("FAMILY"), 
                                !!rlang::sym("A1"))[["CALL"]]
    }
    
    return(sorted_genes)
}
