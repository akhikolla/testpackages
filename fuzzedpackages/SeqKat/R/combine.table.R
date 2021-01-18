#' Combine Table
#'
#' Merges overlapped windows to identify genomic boundaries of kataegic events. This function also assigns hypermuation and kataegic score for combined windows
#'
#'
#' @param test.table Data frame of kataegis test scores
#' @param somatic Data frame of somatic variants
#' @param mutdistance The maximum intermutational distance allowed for SNVs to be grouped in the same kataegic event. Recommended value: 3.2
#' @param segnum Minimum mutation count. The minimum number of mutations required within a cluster to be identified as kataegic. Recommended value: 4
#' @param output.name Name of the generated output directory.
#'
#' @examples
#' \dontrun{
#' combine.table(test.table, somatic, mutdistance,segnum,output.name)
#' }
#'
#' @author Fouad Yousif
#' @author Fan Fan

combine.table <- function(test.table, somatic, mutdistance,segnum,output.name){
	sig.table <- test.table[test.table$significance > 1,];
	# need a ifelse here, if no rows in sigtable, then OK, no potential kataegis events
	# if yes, perform following tests
	if (!is.null(sig.table)){
		if (dim(sig.table)[1] >0){
			# combine the table here.
			start.raw <- sig.table[,11];
			end.raw <- sig.table[,12];
			score.raw <- sig.table[,16];
			score.k.raw <- sig.table[,17];
			start.mut <- start.raw[1];
			end.mut <- end.raw[1];
			score.mut <- score.raw[1];
			score.k.mut <- score.k.raw[1];

			if (length(start.raw) == 1){
				start.mut <- start.mut;
				end.mut <- end.mut;
				score.mut <- score.mut;
				score.k.mut <- score.k.mut;
				} 
			else {
				# combine the windows with any overlapped elements
				for ( i in 2:length(start.raw)){
					if (start.raw[i] < tail(end.mut,1)){
						end.mut[length(end.mut)] <- end.raw[i]
						if(score.raw[i]>score.raw[i-1]){
						score.mut[length(score.mut)]<-score.raw[i]
						score.k.mut[length(score.k.mut)]<-score.k.raw[i]
						}else{
						score.mut[length(score.mut)]<-score.raw[i-1]
						score.k.mut[length(score.k.mut)]<-score.k.raw[i-1]
						}
						#score.mut[length(score.mut)]<-max(score.raw[2:length(start.raw)])
						#score.k.mut[length(score.k.mut)]<-max(score.k.raw[2:length(start.raw)])
						} 
					else {
						start.mut <- c(start.mut,start.raw[i]);
						end.mut <- c(end.mut,end.raw[i]);
						score.mut <-c(score.mut,score.raw[i]);
						score.k.mut <- c(score.k.mut,score.k.raw[i]);
						}
					}
				}
			somatic.subset <- somatic[somatic$CHR == paste0('chr', unique(sig.table$chromosome)),];
			raw.table <- vector();
			score <- 0;
			score.k <- 0;
			for( x in 1:length(start.mut)) {
				index1 <- which(somatic.subset$POS == start.mut[x]);
				index2 <- which(somatic.subset$POS == end.mut[x]);
				distance <- log10(diff(somatic.subset$POS[index1:index2]));
				check <- distance < mutdistance;
				status <- c(rle(check)$values,as.logical(1-tail(rle(check)$values,1)));
				length.seg <- c(rle(check)$lengths,0);
				start.position <- cumsum(c(1,rle(check)$lengths))[status&length.seg>segnum];
				num.mutations <- length.seg[status&length.seg>segnum]+1;
				start.pos <- somatic.subset$POS[index1:index2][start.position];
				end.pos <- somatic.subset$POS[index1:index2][start.position+length.seg[status&length.seg>segnum]];

				score <- score.mut[x];
				score.k <- score.k.mut[x];

				raw.table <- rbind(
					raw.table,
					cbind(
						rep(unique(as.character(sig.table$sample.name)),length(start.pos)),
						rep(unique(sig.table$chromosome),length(start.pos)),
						start.pos,
						end.pos,
						num.mutations,
						rep(unique(score),length(start.pos)),
						rep(unique(score.k),length(start.pos))
						)
					);
				}
			if(dim(raw.table)[1] > 0){
				final.result <- as.data.frame(raw.table);
				names(final.result) <- c('sample','chr','start','end','variants','score.hm','score.kat');
				write.table(
					final.result,
					file = paste0(output.name,'/',unique(test.table$sample.name),'_chr',unique(test.table$chromosome),'_cutoff',unique(test.table$cutoff),'_mutdist',mutdistance,'_segnum',segnum,'.txt'),
					sep="\t",
					row.names=FALSE,
					quote=FALSE
					);
				} else { final.result <- NULL;}
			} else { final.result <- NULL;}
		} else { final.result <- NULL; }
	}
