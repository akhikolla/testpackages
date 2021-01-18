#' Final Score
#'
#' Assigns hypermutation score (hm.score) and kataegic score (k.score)
#'
#' 
#' @param test.table Data frame of kataegis test scores
#' @param cutoff The minimum hypermutation score used to classify the windows in the sliding binomial test as significant windows. The score is calculated per window as follows: -log10(binomial test p-value). Recommended value: 5
#' @param somatic Data frame of somatic variants
#' @param output.name Name of the generated output directory.
#'
#' @examples
#' \dontrun{
#' final.score(test.table, cutoff, somatic, output.name)
#' }
#'
#' @author Fan Fan
#' @author Fouad Yousif

final.score <- function(test.table, cutoff, somatic, output.name){
	if (dim(test.table)[1] > 0){
		temp <- sapply( 
			as.character(test.table$base.change[as.character(test.table$base.change) != 'base change']),
			function(x) strsplit(x, ', ')
			);
		temp2 <- sapply(temp,
			function(x) {
				y <-sapply(x,function(z) strsplit(z, ':')); 
				sum(as.numeric(sapply(y,'[[',2))[ifelse(substr(sapply(y,'[[',1),1,2)=='TC',TRUE,FALSE)])
				}
			);	
		test.table$num.tcx <- rep(0, dim(test.table)[1]);
		test.table$num.tcx[as.character(test.table$base.change) != 'base change'] <- as.vector(temp2);
		test.table[,c(1:3,5:8,10:14)] <- apply(test.table[,c(1:3,5:8,10:14)],2,function(x) as.numeric(as.character(x)));
		test.table$exp.tcx <- ifelse(test.table$exp.tcx == 0, 1, test.table$exp.tcx);
		test.table$hm.score <- -log10(as.numeric(test.table$p.value))*(test.table$num.snvs/test.table$exp.snvs);
		#make sure that both num.tcx and exp.tcx are numeric
		if(is.numeric(test.table$num.tcx) & is.numeric(test.table$exp.tcx)){
			test.table$k.score <- test.table$hm.score*(test.table$num.tcx/test.table$exp.tcx);
			}
		else{
			#message(paste("pvalue is: ",test.table$p.value,sep=""))
			#message(paste("num snvs is: ",test.table$num.snvs,sep=""))
			#message(paste("exp snvs is: ",test.table$exp.snvs,sep=""))
			test.table$k.score <- NA;
			}
		test.table$significance <- as.numeric(test.table$hm.score >= cutoff) + 1;
		test.table$p.value[test.table$p.value == 1000] <- 1;
		test.table$window.medians <- sapply(
			1:dim(test.table)[1],
			function(x) median(test.table[x,2], test.table[x,3])
			);
		test.table$cutoff <- rep(cutoff,length(test.table$hm.score));
		if (typeof(test.table$num.tcx) == "list") {
			test.table$num.tcx <- unlist(test.table$num.tcx);
			}
		write.table(
			test.table,
			file = paste0(output.name,'/finaltable_chr',unique(test.table$chromosome),'_cutoff',cutoff,'.txt'),
			row.names = FALSE,
			sep = '\t'
			);
		return(test.table);
		} else { test.table <- NULL; return(test.table); } 
	}
