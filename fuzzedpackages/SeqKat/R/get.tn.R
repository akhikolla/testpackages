#' Get Trinucleotides
#'
#' Count the frequencies of 32 trinucleotide in a region respectively
#'
#'
#' @param chr Chromosome
#' @param start.bp Starting position
#' @param end.bp Ending position
#' @param ref.dir Path to a directory containing the reference genome.
#'
#' @examples
#' \dontrun{
#' get.tn(chr,start.bp,end.bp,ref.dir)
#' }
#'
#' @author Fan Fan

### FUNCTION: GET.TN  ##############################################################
# count the frequencies of 32 trinucleotide in a region respectively 

get.tn<- function(chr,start.bp,end.bp,ref.dir){
	
	bases.raw <- c('A','C','G','T','N');
	tri.types.raw <- c(outer( c(outer(bases.raw, bases.raw, function(x, y) paste0(x,y))), bases.raw, function(x, y) paste0(x,y)));
	tri.types.raw <- sort(tri.types.raw);
	tri.types.filter <- !grepl('N', tri.types.raw);
	tri.types <- tri.types.raw[tri.types.filter];

	tri.counts <- get.nucleotide.chunk.counts(
		key = tri.types.raw, 
		chr = file.path(ref.dir, paste0( paste('chr',chr,sep = ''), '.fa')),
		start = start.bp,
		end = end.bp
		);
	tri.counts <- tri.counts[tri.types.filter];
	names(tri.counts) <- tri.types;


# collapse counts
	base.pairs <- c('T' = 'A', 'G' = 'C', 'C' = 'G', 'A' = 'T');

	tri.count <- data.frame(
		trinucleotide = names(tri.counts), 
		count = tri.counts, 
		stringsAsFactors = FALSE
		);


	tri.count$trinucleotide <- sapply(
		X   = tri.count$trinucleotide, 
		FUN = function(x) {
			ifelse(
				test = substr(x, 2, 2) %in% c('C', 'T'), 
				yes  = x, 
				no   = paste(base.pairs[rev(unlist(strsplit(x,'')))], collapse = ''))
			}
		);

	tri.count.collapsed <- tapply( 
		X = tri.count$count, 
		INDEX = tri.count$trinucleotide, 
		sum
		);

	tri.count.collapsed <- tri.count.collapsed[order(names(tri.count.collapsed))];
	return(tri.count.collapsed);

	}
