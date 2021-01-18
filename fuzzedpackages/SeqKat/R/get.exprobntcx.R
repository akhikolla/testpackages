library(doParallel);
globalVariables("chr");

#' get.exprobntcx
#'
#' Gets the expected probability for each trinucleotide and total number of tcx
#'
#'
#' @param somatic Data frame of somatic variants
#' @param ref.dir Path to a directory containing the reference genome.
#' @param trinucleotide.count.file A tab seprarated file containing a count of all trinucleotides present in the reference genome. This can be generated with the get.trinucleotide.counts() function in this package.
#'
#' @examples
#' \dontrun{
#' get.exprobntcx(somatic, ref.dir, trinucleotide.count.file)
#' }
#'
#' @author Fan Fan
#' @author Fouad Yousif

get.exprobntcx <- function(somatic, ref.dir, trinucleotide.count.file){
	ref.tn.total <- vector();
	tn.tcx <- vector();
	foreach( chr = 1:24) %do% {
		chromosome.subset <- subset(
			somatic,
			somatic$CHR == paste('chr',chr,sep='')
			);
		tn.index <- chromosome.subset$POS; 
		ref.tn <- get.context(
			file.path(ref.dir,paste(paste('chr',chr,sep=''),'fa',sep='.')),
			tn.index
			);
		ref.tn <- toupper(ref.tn);
		ref.tn.rev <-  ifelse(substr(ref.tn,2,2)%in%c("A","G"),get.pair(ref.tn),ref.tn);
		ref.tn.total <- c(ref.tn.total,ref.tn.rev);
		tcx.chr <- sum(substr(ref.tn.rev,1,2) == 'TC');
		# get TCX in all chromosomes
		tn.tcx <- c(tn.tcx,tcx.chr);
		}
	ref.mut <- as.data.frame(table(ref.tn.total));
	ref.genome <- read.table(
		file=trinucleotide.count.file,
		header = TRUE,
		);
	ref.mut[,1] <- as.character(ref.mut[,1]);
	ref.genome[,1] <- as.character(ref.genome[,1]);
	# assign 0s to the trinucleotide without mutations
	if(	dim(ref.mut)[1] < 32) {
		ref.mut <- 	data.frame(
			trinucleotide = c(ref.mut[,1],ref.genome[,1][is.na(pmatch(ref.genome[,1],ref.mut[,1]))]),
			count = c(ref.mut[,2],rep(0,sum(is.na(pmatch(ref.genome[,1],ref.mut[,1])))))
			)
		}
	ref.mut <- ref.mut[order(ref.mut[,1]),];
	exp.prob <- ref.mut[,2]/ref.genome[,2];
	names(exp.prob) <- ref.genome[,1];
	# result include 32 probabilities and count of tcx
	return(c(exp.prob,tn.tcx));
	}
