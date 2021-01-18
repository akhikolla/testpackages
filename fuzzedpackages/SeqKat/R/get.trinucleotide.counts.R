library(doParallel);
globalVariables("chr");

#' Get Trinucleotide Counts
#'
#' Aggregates the total counts of each possible trinucleotide.
#'
#'
#' @param ref.dir Path to a directory containing the reference genome.
#' @param ref.name Name of the reference genome being used (i.e. hg19, GRCh38, etc)
#' @param output.dir Path to a directory where output will be created.
#'
#' @examples
#' \dontrun{
#' get.trinucleotide.counts(ref.dir, ref.name, output.dir)
#' }
#'
#' @author Fan Fan
#' @author Fouad Yousif

get.trinucleotide.counts <- function(ref.dir, ref.name, output.dir) {

    # for multi-core parallel
    doMC::registerDoMC(cores = 24);

    chrs <- paste0('chr', c(as.character(1:22), 'X', 'Y'));

    bases.raw <- c('A', 'C', 'G', 'N', 'T');
    tri.types.raw <- c(outer( c(outer(bases.raw, bases.raw, function(x, y) paste0(x,y))), bases.raw, function(x, y) paste0(x,y)));
    tri.types.raw <- sort(tri.types.raw);
    tri.types.filter <- !grepl('N', tri.types.raw)
    tri.types <- tri.types.raw[tri.types.filter];

    # multi-core parallel. In tri.counts, row=tri.types, col=chromosomes

    tri.counts <- foreach(chr = chrs, .combine = cbind) %dopar% {

        c.counts <- cget_trinucleotide_counts(key = tri.types.raw, chr = file.path(ref.dir, paste0(chr, '.fa')));
        c.counts[tri.types.filter]
        
    };
        
    rownames(tri.counts) <- tri.types;

    write.table(x = tri.counts, file = file.path(output.dir, paste0('trinucleotide_counts_on_', ref.name, '_by_chr.txt')), sep = '\t', quote = FALSE);

    # collapsed counts
    base.pairs <- c('T', 'G', 'C', 'A');
    names(base.pairs) <- base.pairs;

    tri.count.genomewise <- data.frame(trinucleotide = rownames(tri.counts), count = rowSums(tri.counts), stringsAsFactors = FALSE); 

    tri.count.genomewise$trinucleotide <- sapply(
        X   = tri.count.genomewise$trinucleotide, 
        FUN = function(x) {
            ifelse(
                test = substr(x, 2, 2) %in% c('C', 'T'), 
                yes  = x, 
                no   = paste0(base.pairs[rev(unlist(strsplit(x,'')))], collapse = ''))
            }
        );

    tri.count.genomewise.collapsed <- tapply( X = tri.count.genomewise$count, INDEX = tri.count.genomewise$trinucleotide, sum);
    tri.count.genomewise.collapsed <- tri.count.genomewise.collapsed[order(names(tri.count.genomewise.collapsed))];

    write.table(
        x = data.frame(trinucleotide = names(tri.count.genomewise.collapsed), count = tri.count.genomewise.collapsed),
        file = file.path(output.dir, paste0('trinucleotide_', ref.name, '_whole_genome.txt')), 
        sep = '\t', quote = FALSE, row.names = FALSE
        );

    }

# function to collapse count matrix, i.e., ACG is equivalent to CAT.
tricount.collapse <- function(counts) {
    base.pairs <- c('A'='T', 'C'='G', 'G'='C', 'T'='A');
    counts.df <- data.frame(triname = names(counts), count = counts, stringsAsFactors = FALSE);

    counts.df$triname <- sapply(
        X = counts.df$triname,
        FUN = function(x) {
            ifelse(substr(x, 2, 2) %in% c('C', 'T'), x, paste0(base.pairs[rev(unlist(strsplit(x,'')))], collapse = ''))
            }
        );
    tapply(
        X = counts.df$count,
        INDEX = counts.df$triname,
        FUN = sum
        );
    }
