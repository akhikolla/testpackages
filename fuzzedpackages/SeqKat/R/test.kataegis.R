#' Test Kataegis
#'
#' Performs exact binomial test to test the deviation of the 32 tri-nucleotides counts from expected
#'
#' 
#' @param chromosome.num Chromosome
#' @param somatic Data frame of somatic variants
#' @param units Base window size
#' @param exprobntcx Expected probability for each trinucleotide and total number of tcx
#' @param output.name Name of the generated output directory.
#' @param ref.dir Path to a directory containing the reference genome.
#' @param chromosome.length.file A tab separated file containing the lengths of all chromosomes in the reference genome.
#'
#' @examples
#' \dontrun{
#' test.kataegis(
#'	chromosome.num,
#'	somatic,units,
#'	exprobntcx,
#'	output.name,
#'	ref.dir,
#'	chromosome.length.file
#'	);
#' }
#'
#' @author Fouad Yousif

test.kataegis <- function(
	chromosome.num,
	somatic,
	units,
	exprobntcx,
	output.name,
	ref.dir,
	chromosome.length.file
	){
	# suppose all are female patients here: use the length without chrY
	chromosome.subset <- somatic[somatic$CHR == paste0('chr',chromosome.num),];
	chr.length <- read.table(file=chromosome.length.file, header = TRUE, stringsAsFactors = FALSE);
	average.length <- floor(chr.length[25,2]/dim(somatic)[1]);
	# unit is the expected mutation number the user wanna include in each window, defult = 2
	window.size <- average.length*(units - 1);
	# make sure the window size is a even number
	window.size <- ifelse(window.size %% 2 ==0, window.size, window.size+1);
	chromosome.subset.length=dim(chromosome.subset)[1];
	#check if there are any variants actually present in that chromosome
	if (chromosome.subset.length>0){
	#check if the number of variants is larger than the sliding window.
	#if(window.size+1<chromosome.subset.length){
	message(paste("Testing Chromosome ",chromosome.num,sep=''));
	#update the chromosome subset table to include variant order, REF/ALT, and intermutation distance
	chromosome.subset<-data.frame(
		1:(chromosome.subset.length),
		chromosome.subset,
		paste(chromosome.subset[1:chromosome.subset.length,c('REF')],chromosome.subset[1:chromosome.subset.length,c('ALT')],sep='/'),
		c(0,log10(apply(as.matrix(chromosome.subset$POS),2,diff)))
		);
	names(chromosome.subset)<-c('VAR','CHR','POS','REF','ALT','BASE','DIST'); # DIST: distance
	#measure the distance between the first and last position on the chromosome, this will probably need to be changed to the actual chromosome length in bases (either use supplied or we add it to the beginning of the code)
	chromosome.bp=max(chromosome.subset$POS)-min(chromosome.subset$POS);
	#number of full size windows
	segment.rawnum <- chr.length[as.character(chromosome.num),2]/window.size;


	if (segment.rawnum < 1){
		cutpoint <- c(1, chr.length[as.character(chromosome.num),2]);
		loops <- 1;
		grouping <- split(
			chromosome.subset,
			cut(
				chromosome.subset$POS,
				breaks <- (c(0,chr.length[as.character(chromosome.num),2]))
			   )
			);
		} 
	else {
		segment.num <- ifelse(
			segment.rawnum - floor(segment.rawnum) < 0.5,
			floor(segment.rawnum),
			floor(segment.rawnum) + 0.5
			);
		
		# flag the position of the start points of all windows
		cutpoint <- c(seq(0,segment.num,0.5)*window.size + 1, chr.length[as.character(chromosome.num),2]);
		# total number of windows
		loops <- length(cutpoint) - 2;
	
		# group the mutations into each half window
		grouping <- split(
			chromosome.subset,
			cut(
				chromosome.subset$POS,
				breaks = (c(seq(0, segment.num,0.5)*window.size,chr.length[as.character(chromosome.num),2]))
				)
			);
		}
	# indicator indicates how many mutations in each interval
	indicator <- as.vector(sapply(grouping,dim)[1,]);
	# loops = length(indicator) - 1
	# check total mutation number of the current window and the next window, if >1, will apply binomial test
	test.indicator <- c();
 	# message(length(indicator))
	for (i in 1:loops){
		if (length(indicator)>1){
			if (sum(indicator[i:(i+1)]) < 2){
				test.indicator[i] = 0;
				}
			else {
				test.indicator[i] = 1;
				}
			}
		else {
			test.indicator[i] = 0;
			}
		}

	# if want to check the mutations stored in each interval, use grouping[[i]]

	pvalue <- matrix(
		nrow = loops,
		ncol = 14
		); 
	# prepare a new table to save the sliding window test results. add one more col for trinucleotide.
	# apply biomial test for the selected windows
	for(i in 1:loops){
		if (test.indicator[i] == 1){
			mutation.inrange <- do.call(rbind,grouping[i:(i + 1)]);
			start.bp <- cutpoint[i]; #pull the starting position for the sliding window
			end.bp <- cutpoint[i + 2] - 1; 
			window.bp <- abs(end.bp-start.bp); #calculate the #basepairs that span the sliding window
			tn.table <- get.tn(chr = chromosome.num, start.bp = start.bp, end.bp = end.bp, ref.dir = ref.dir); 
			total.tn <-sum(tn.table);
			exp.prob <- exprobntcx[1:32];
			num.snvs <- dim(mutation.inrange)[1];
			exp.snvs <- tn.table %*% exp.prob;
			exp.tcx <- exprobntcx[32+chromosome.num]*(end.bp - start.bp)/chr.length[as.character(chromosome.num),2];
			expected.freq <- exp.snvs/total.tn;
			observed.freq <- num.snvs/total.tn;
			begin.snv <- mutation.inrange$POS[1];
			end.snv <- tail(mutation.inrange$POS,1);
			tn_summary <- get.toptn(chromosome.subset, chromosome.num, begin.snv, end.snv, ref.dir);
			raw.pvalue <- binom.test(dim(mutation.inrange)[1],total.tn,p=expected.freq,alternative="greater")$p.value; 
			pvalue[i,1:13]=c(chromosome.num,start.bp,end.bp,tn_summary,observed.freq,expected.freq,raw.pvalue,num.snvs,output.name,exp.snvs,begin.snv,end.snv,exp.tcx);
			}
		else{
			start.bp <- cutpoint[i];
			end.bp  <- cutpoint[i + 2] - 1;
			window.bp  <- abs(end.bp - start.bp);
			tn.table <- -100;
			total.tn <- -100;
			expected.freq <- -100;
			observed.freq <- -100;
			raw.pvalue <- -100;
			base.change <- -100; 
			pvalue[i,1:13]=c(chromosome.num,start.bp,end.bp,'base change',observed.freq,expected.freq,raw.pvalue,0,output.name,1000,-100,-100,1000 );
			}
		}


	adjusted.pvalue <- vector(length = loops);
	adjusted.pvalue[which(test.indicator == 1)]	<-p.adjust(as.numeric(pvalue[which(test.indicator == 1),7]),method="fdr"); #adjust the p-value from the sliding window test for multiple hypothesis testing using method FDR
	adjusted.pvalue[which(test.indicator == 0)] <- 1000;
	adjusted.pvalue.binary <- as.numeric(adjusted.pvalue<0.05)+1; #convert the p-value into a binary score (<0.05 =1 and >0.05 =0)
	# do the following to those with tests  
	# this binary value could be ignored, will replace by another significant indicator in another function
	pvalue[,14] <- adjusted.pvalue.binary;
	pvalue[,7] <- adjusted.pvalue;
	chromosome.score = na.omit(pvalue);
	#message(paste(length(which(pvalue[,14]>1))," potential kataegis events identified ",sum(test.indicator)," tests done",sep=''));
	#       }else{
	#           chromosome.score = NA;
	#			message(paste("Skipping Chromosome ",chromosome.num," Your window size is large. Not sufficient somatic events.",sep=''));
	#       }
	}
	else{
		chromosome.score = NA;
		message(paste("Skipping Chromosome ",chromosome.num," Not present in your bed file",sep=''));
		}
	combined <-data.frame(na.omit(rbind(chromosome.score))) ;
	names(combined) <- c("chromosome","window.start","window.end","base.change","observed.frequency","expected.frequency","p.value","num.snvs","sample.name","exp.snvs","1st.snvPOS", "last.snvPOS","exp.tcx","significance");

	# combined<-data.frame(rbind(na.omit(chromosome.score))) #combined all the chromosomes data into one table
	# names(combined)<-c("chromosome","window.start","window.end","base.change","observed.frequency","expected.frequency","p.value","num.snvs","sample.name","exp.snvs","1st.snvPOS", "last.snvPOS","exp.tcx","significance");
	# write.table(combined,paste(output.name,"/kataegisPvalues_units", units, "_weight", weight, "_chr",chromosome.num,".txt",sep=''),sep="\t",row.names=FALSE,quote=FALSE)
	return(combined)
	};
