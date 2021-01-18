library(Rcpp);

#' Get Nucleotide Chunk Counts
#'
#' Obtain counts for all possible trinucleotides within a specified genomic region
#'
#' 
#' @param key List of specify trinucleotides to count
#' @param chr Chromosome
#' @param upstream Length upstream to read
#' @param downstream Length downstream to read
#' @param start Starting position
#' @param end Ending position
#'
#'
#' @author Fouad Yousif

get.nucleotide.chunk.counts  <-  function (key, chr, upstream = 1, downstream = 1, start = 1, end = -1){
	cget_nucleotide_chunk_counts(key, chr, upstream, downstream, start, end);
	};
