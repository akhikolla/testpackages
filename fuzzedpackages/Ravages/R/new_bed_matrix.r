# mafs = une matrice de mafs comme renvoyée par group.mafs
# ie autant de colonnes que de SNP, autant de ligne que de groupes
# size = un vecteur donnant la taille de chaque groupe

new.bed.matrix <- function(nb_inds, nb_snps) {
  bed <- .Call('oz_new_bed_matrix', PACKAGE = "Ravages", nb_snps, nb_inds)

  ids <- sprintf("A%0*d", log10(nb_inds) + 1, 1:nb_inds)
  ped <- data.frame(famid = ids,  id = ids, father = 0, mother = 0, sex = 0,
            pheno = 0, stringsAsFactors = FALSE)

  ids <- sprintf("m%0*d", log10(nb_snps) + 1, 1:nb_snps)
  snps <- data.frame(chr = NA, id = ids, dist = NA, pos = NA,
               A1 = NA, A2 = NA, stringsAsFactors = FALSE)

  x <- new("bed.matrix", bed = bed, snps = snps, ped = ped, p = NULL, mu = NULL,
           sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )
  if(getOption("gaston.auto.set.stats", TRUE))
    x <- set.stats(x, verbose = FALSE)
  x
}

