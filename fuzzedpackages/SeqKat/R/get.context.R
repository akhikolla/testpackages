library(Rcpp);

#' Get Context
#'
#' Gets the 5' and 3' neighboring bases to the mutated base
#'
#'
#' @param file Reference files directory
#' @param start The position of the mutation gene

#' @return The trinucleotide context.
#'
#' @examples
#' \dontrun{
#' get.context(file.path(referencie.genome.dir, 'chr1.fa'), c(158297133, 161176181))
#' }
#'
#' @author Fouad Yousif
#' @author Fan Fan


get.context  <-  function (file,start){
	cpp_get_context(file, start, length(start));
	};
