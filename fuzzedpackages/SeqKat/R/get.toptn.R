#' Get Top Trinucleotides
#'
#' Generate a tri-nucleotide summary for each sliding window
#'
#' 
#' @param somatic.subset Data frame of somatic variants subset for a specific chromosome
#' @param chr Chromosome
#' @param start.bp Starting position
#' @param end.bp Ending position
#' @param ref.dir Path to a directory containing the reference genome.
#'
#' @examples
#' \dontrun{
#' get.toptn(somatic.subset, chr, start.bp, end.bp, ref.dir)
#' }
#'
#' @author Fan Fan
#' @author Fouad Yousif

get.toptn <- function(somatic.subset, chr, start.bp, end.bp, ref.dir){
	dict <- c('A' = 'T', 'C' = 'G', 'G' = 'C', 'T' = 'A', 'N' = 'N');
	
	index <- which(somatic.subset$POS %in% c(start.bp, end.bp));
	index <- index[1]:index[2]
	mut.index <- somatic.subset$POS[index];
	mut.tn <- toupper(get.context(
			file.path(ref.dir,paste(paste('chr',chr,sep=''),'fa',sep='.')),
			mut.index
			));
	mut.snv <- toupper(somatic.subset$ALT[index]);
	mut.tn.adj <- ifelse(substr(mut.tn,2,2) %in% c('A','G'), get.pair(mut.tn),mut.tn);
	mut.snv.adj <- ifelse(substr(mut.tn,2,2) %in% c('A','G'),dict[mut.snv],mut.snv);
	changes <- paste(mut.tn.adj, mut.snv.adj, sep ='>');
	tmp <- table(changes);
	return(paste(paste(dimnames(tmp)[[1]],as.vector(tmp),sep = ':'), collapse = ', '));
	}
