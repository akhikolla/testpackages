rbm.haplos.freqs <- function(haplos, freqs, size, replicates) {
  bed <- .Call('rbm_haplos_freqs', PACKAGE = "Ravages", haplos, freqs, size, replicates)

  nb_inds <- sum(size);

  ids <- sprintf("A%0*d", log10(nb_inds) + 1, 1:nb_inds)

  if( ncol(freqs) != length(size) )
    stop("Dimensions mismatch")

  if(is.null(colnames(freqs)))
    lev <- (1:ncol(freqs)) - 1
  else
    lev <- colnames(freqs)

  pheno <- factor( unlist(mapply(rep, lev, each = size, SIMPLIFY = FALSE)) , levels = lev )

  ped <- data.frame(famid = ids,  id = ids, father = 0, mother = 0, sex = 0,
            pheno = pheno, stringsAsFactors = FALSE)

  snps <- data.frame(chr = NA, id = NA, dist = NA, pos = NA, A1 = NA, A2 = NA, 
               genomic.region = factor( rep(sprintf("R%0*d", log10(replicates) + 1, 1:replicates), each = ncol(haplos)) ),
               stringsAsFactors = FALSE)
 
  if( is.null(colnames(haplos)) )
    snps$id <- paste( snps$genomic.region, sprintf("SNP%0*d", log10(ncol(haplos)) + 1, 1:ncol(haplos)), sep = "_")
  else
    snps$id <- paste( snps$genomic.region, colnames(haplos), sep = "_")

  x <- new("bed.matrix", bed = bed, snps = snps, ped = ped, p = NULL, mu = NULL,
           sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats(x, verbose = FALSE)
  x
}
