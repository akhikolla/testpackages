get_sample_names <- function(file, connection_type=base::file, description_columns=9) {
	cnx <- connection_type(file, open="r")
	on.exit(close(cnx))
	headings <- NULL
	while (is.null(headings)) {
		l <- readLines(cnx, n=1)
		if (grepl(x=l, pattern="^#CHROM")) headings <- strsplit(l, split="\t")[[1]]
	}
	headings[(description_columns+1):length(headings)]
}

compressed_vcf_sample_names <- function(vcf_file_name) {
	get_sample_names(paste0("zcat ", vcf_file_name), connection_type=pipe)
}

just_counts <- function(parts, description_columns=9) {
	y <- sapply(parts, "[", -(1:description_columns))
	structure(grepl(x=y, pattern="^[^0.][/|].") + grepl(x=y, pattern="^.[/|][^0.]"), dim=c(if (length(parts) > 0) length(parts[[1]])-description_columns else 0, length(parts)))
}

var_info <- function(parts, description_columns=9) structure(dimnames=list(NULL, c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","FORMAT","INFO")), structure(dim=c(length(parts),description_columns), t(sapply(parts, "[", seq(length.out=description_columns)))))

get_block_parts <- function(vcf_file_name, chr, from, to) {
	cmd <- paste("tabix ", vcf_file_name, " ", chr, ":", from, "-", to, sep="")
	z <- pipe(cmd)
	lines <- grep(value=TRUE, pattern="^#", invert=TRUE, x=readLines(z))
	close(z)
	strsplit(lines, split="\t")
}

vcf2matrix <- function(vcf_file_name, chr, from, to, samples=compressed_vcf_sample_names(vcf_file_name), include_variant_info=FALSE, description_columns=9, warn_if_AF_greater_than=0.1) {
	file_samples <- compressed_vcf_sample_names(vcf_file_name) 


	sample_inds <- match(samples, file_samples)
	parts <- get_block_parts(vcf_file_name, chr, from, to)

	info <- var_info(parts, description_columns)

	if (length(parts) == 0)
		stop("No data")

	if (any(grepl(",",info[,"ALT"])))
		stop("Region contains variants with multiple alternate alleles: please use 'bcftools norm' to split such variants onto separate rows")

	counts <- just_counts(parts, description_columns)[sample_inds,,drop=FALSE]
	if (any(apply(counts > 0, 2, mean) > warn_if_AF_greater_than))
		warning(paste0("Allele count matrix contains variants with allele frequency greater than ", warn_if_AF_greater_than))

	if (include_variant_info)
		list( info=info, G=counts)
	else 
		counts
}
