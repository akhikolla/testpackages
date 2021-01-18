#' Get Pair
#'
#' Generates the reverse compliment of a nucleotide sequence
#'
#' Reverses and compliments the bases of the input string. Bases must be (A, C, G, T, or N).
#' 
#' @param x asdf
#'
#' @examples
#' \dontrun{
#' get.pair("GATTACA")
#' }
#'
#' @author Fouad Yousif

get.pair <- function(x) {
	dict <- c('A' = 'T', 'C' = 'G', 'G' = 'C', 'T' = 'A', 'N' = 'N');
	neu <- strsplit(x, '');
	sapply(neu, function(z) paste(rev(dict[z]), collapse=''))
	};
