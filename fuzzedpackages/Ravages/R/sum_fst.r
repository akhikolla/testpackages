Sum.Fst <- function(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region, 
                   which.snps = rep(TRUE, ncol(x)), target = 10, B.max = 1e6) {
  r <- .Call('oz_sum_fst', PACKAGE = "Ravages", x@bed, which.snps, as.factor(genomic.region), as.factor(group), target, B.max) 
  res <- data.frame( genomic.region = levels(genomic.region), r )
  return(res)
}

Sum.Fst.higher.perms <-  function(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region,
                   which.snps = rep(TRUE, ncol(x)), n.keep = 500, B = 1e6) {
  r <- .Call('oz_sum_fst_max_perm', PACKAGE = "Ravages", x@bed, which.snps, as.factor(genomic.region), as.factor(group), n.keep, B)
  names(r) <- levels(genomic.region)
  r[ sapply(r, length) > 0 ] 
}

# pour comparer tps d'exec
Sum.Fst1 <- function(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region, 
                   which.snps = rep(TRUE, ncol(x)), target = 10, B.max = 1e6) {
  r <- .Call('oz_sum_fst1', PACKAGE = "Ravages", x@bed, which.snps, as.factor(genomic.region), as.factor(group), target, B.max) 
  res <- data.frame( genomic.region = levels(genomic.region), r )
  return(res)
}

Sum.Fst1.higher.perms <-  function(x, group = x@ped$pheno, genomic.region = x@snps$genomic.region,
                   which.snps = rep(TRUE, ncol(x)), n.keep = 500, B = 1e6) {
  r <- .Call('oz_sum_fst1_max_perm', PACKAGE = "Ravages", x@bed, which.snps, as.factor(genomic.region), as.factor(group), n.keep, B)
  names(r) <- levels(genomic.region)
  r[ sapply(r, length) > 0 ] 
}

