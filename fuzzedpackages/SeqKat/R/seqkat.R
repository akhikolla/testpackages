library(Rcpp);
library(foreach);

#' SeqKat
#'
#' Kataegis detection from SNV BED files
#'
#' The default paramters in SeqKat have been optimized using Alexanrov's "Signatures of mutational processes in human cancer" dataset. SeqKat accepts a BED file and outputs the results in TXT format. A file per chromosome is generated if a kataegic event is detected, otherwise no file is generated. SeqKat reports two scores per kataegic event, a hypermutation score and an APOBEC mediated kataegic score.
#' 
#' @param sigcutoff The minimum hypermutation score used to classify the windows in the sliding binomial test as significant windows. The score is calculated per window as follows: -log10(binomial test p-value). Recommended value: 5
#' @param mutdistance The maximum intermutational distance allowed for SNVs to be grouped in the same kataegic event. Recommended value: 3.2
#' @param segnum Minimum mutation count. The minimum number of mutations required within a cluster to be identified as kataegic. Recommended value: 4
#' @param ref.dir Path to a directory containing the reference genome. Each chromosome should have its own .fa file and chromosomes X and Y are named as chr23 and chr24. The fasta files should contain no header
#' @param bed.file Path to the SNV BED file. The BED file should contain the following information: Chromosome, Position, Reference allele, Alternate allele
#' @param output.dir Path to a directory where output will be created.
#' @param chromosome The chromosome to be analysed. This can be (1, 2, ..., 23, 24) or "all" to run sequentially on all chromosomes.
#' @param chromosome.length.file A tab separated file containing the lengths of all chromosomes in the reference genome.
#' @param trinucleotide.count.file A tab seprarated file containing a count of all trinucleotides present in the reference genome. This can be generated with the get.trinucleotide.counts() function in this package.
#'
#' @examples
#' example.bed.file <- paste0(
#'	path.package("SeqKat"),
#'	"/extdata/test/PD4120a-chr4-1-2000000_test_snvs.bed"
#'	);
#' example.ref.dir <- paste0(
#'	path.package("SeqKat"),
#'	"/extdata/test/ref/"
#'	);
#' example.chromosome.length.file <- paste0(
#'	path.package("SeqKat"),
#'	"/extdata/test/length_hg19_chr_test.txt"
#'	);
#' seqkat(
#'		5,
#'		3.2,
#'		2,
#'		bed.file = example.bed.file,
#'		output.dir = ".",
#'		chromosome = "4",
#'		ref.dir = example.ref.dir,
#'		chromosome.length.file = example.chromosome.length.file
#'		);
#'
#' @author Fouad Yousif
#' @author Fan Fan
#' @author Christopher Lalansingh

seqkat <- function(
	sigcutoff = 5,
	mutdistance = 3.2,
	segnum = 4,
	ref.dir  = NULL,
	bed.file = "./",
	output.dir = "./",
	chromosome = "all",
	chromosome.length.file = NULL,
	trinucleotide.count.file = NULL
	) {

	# Validate sigcutoff, mutdistance, segnum
	stopifnot(is.numeric(sigcutoff), is.numeric(mutdistance), is.numeric(segnum));

	# Validate the bed.file
	if (!file.exists(bed.file)) {
		stop("The bed.file provided does not exist.");
		}

	# Validate the chromosome
	if ((chromosome != "all") & (!grepl('^(\\d*)$', chromosome))) {
		stop("The provided value for chromosome should be of the form (1, 2, ..., 23, 24), without a prepended 'chr'.");
		}
	else if (grepl('^(\\d*)$', chromosome)) {
		stopifnot(as.numeric(chromosome) >= 1, as.numeric(chromosome) <= 24, as.numeric(chromosome)%%1 == 0);
		}

	# Validate the ref.dir
	if (is.null(ref.dir)) {
		stop("Please supply a path to the reference genome with the ref.dir argument.")
		}
	if (!all(file.exists(sapply(1:24, function(x) { paste0(ref.dir, "/chr", x, ".fa") })))) {
		stop("The ref.dir directory must contain separate fasta files for each chromosome of the form (chr1.fa, chr2.fa, ..., chr23.fa, chr24.fa).")
		}

	# Validate the chromosome.length.file
	if (is.null(chromosome.length.file)) {
		warning("No chromosome.length.file provided, using hg19 lengths by default.");
		chromosome.length.file <- paste0(path.package("SeqKat"),"/extdata/length_hg19_chr.txt")
		}
	else {
		chr.length <- read.table(file=chromosome.length.file, header = TRUE, stringsAsFactors = FALSE);
		if (is.null(chr.length$num)) {
			stop("The supplied chromosome.length.file is missing the 'num' column.");
			}
		if (is.null(chr.length$length)) {
			stop("The supplied chromosome.length.file is missing the 'length' column.");
			}
		if (!all.equal(c(as.character(1:24), "sum.f", "sum.m"), chr.length$num)) {
			stop("The supplied chromosome.length.file is missing required rows.");
			}
		}

	# Validate the trinucleotide.count.file
	if (is.null(trinucleotide.count.file)) {
		warning("No trinucleotide.count.file provided, using hg19 counts by default. This file can be generated using the get.trinucleotide.counts() function in this package.");
		trinucleotide.count.file <- paste0(path.package("SeqKat"),"/extdata/tn_count.txt")
		}

	if (!file.exists(output.dir)) {
		dir.create(output.dir);
		}

	sample.name <- paste(strsplit(basename(bed.file),'_')[[1]][strsplit(basename(bed.file),'_')[[1]]!='snvs.bed'],collapse='_');
	somatic.directory <- paste0(
		output.dir,
		"/",
		sample.name
		);
	dir.create(somatic.directory);
	somatic.file <- bed.file;
	
	print(somatic.file);
	somatic<-read.delim(
		somatic.file,
		header = TRUE,
		stringsAsFactors = FALSE
		);

	names(somatic) <- c("CHR","POS","REF","ALT");

	#Rename chrX and chrY
	somatic$CHR <- gsub("chrX","chr23",somatic$CHR);
	somatic$CHR <- gsub("chrY","chr24",somatic$CHR);

	output.name <- somatic.directory;
	exprobntcx <- get.exprobntcx(somatic, ref.dir, trinucleotide.count.file);
	if (chromosome == "all") {
		chrs <- sort(na.omit(as.numeric(unlist(strsplit(unique(somatic$CHR),'[^0-9]+')))));
		}
	else {
		chrs <- c(as.numeric(chromosome));
		}
	combined <- sapply(
		chrs,
		function(x){
			combine.table(
				 final.score(
					 test.table =  test.kataegis(
						 chromosome.num = x, 
						 somatic, 
						 units = 2,
						 exprobntcx = exprobntcx,
						 output.name = sample.name,
						 ref.dir = ref.dir,
						 chromosome.length.file = chromosome.length.file
						 ), 
					 cutoff = sigcutoff, 
					 somatic,
					 output.name = output.name
					 ),
				somatic, 
				mutdistance, 
				segnum,
				output.name = output.name
				);
			}
		);
	}
