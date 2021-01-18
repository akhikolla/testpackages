readVCF <- function(filename) {
  vcf <- vcfR::read.vcfR(filename, verbose = FALSE)
  foo <- vcfR::extract.haps(vcf, unphased_as_NA = TRUE, verbose = FALSE)

  if (any(is.na(foo))) {
    stop("VCF data must be phased and without missing values")
  }

  swapallele <- function(x) {
    x[x == vcfR::getREF(vcf)] <- 0
    x[x == vcfR::getALT(vcf)] <- 1

    return(x)
  }

  foo <- apply(foo, 2, swapallele)
  foo <- matrix(as.numeric(foo), nrow = nrow(foo))

  retval <- list()

  rows <- ncol(foo) / 2
  hap <- list()
  hap[[1]] <- matrix(0, nrow = rows, ncol = nrow(foo))
  hap[[2]] <- hap[[1]]

  for (i in 1:rows) {
    hap[[1]][i, ] <- foo[, i * 2 - 1]
    hap[[2]][i, ] <- foo[, i * 2]
  }

  retval$genotypes <- hap2geno(hap)

  retval$map <- data.frame(SNP = vcfR::getID(vcf), chr = vcfR::getCHROM(vcf), pos = vcfR::getPOS(vcf))

  return(retval)
}
