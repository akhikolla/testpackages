#' @title Estimation of recombination rate and maternal LD
#' @name hsrecombi
#' @description Wrapper function for estimating recombination rate and maternal
#'   linkage disequilibrium between intra-chromosomal SNP pairs by calling EM
#'   algorithm
#' @details Paternal recombination rate and maternal linkage disequilibrium (LD)
#'   are estimated for pairs of biallelic markers (such as single nucleotide
#'   polymorphisms; SNPs) from progeny genotypes and sire haplotypes. At least
#'   one sire has to be double heterozygous at the investigated pairs of SNPs.
#'   All progeny are merged in two genomic families: (1) coupling phase family
#'   if sires are double heterozygous 0-0/1-1 and (2) repulsion phase family if
#'   sires are double heterozygous 0-1/1-0.
#'   So far it is recommended processing the chromosomes separately.
#'   If maternal half-sib families are used, the roles of sire/dam are swapped.
#'   Multiple families can be considered.
#' @param hap list (LEN 2) of lists
#' \describe{
#'   \item{famID}{list (LEN number of sires) of vectors (LEN n.progeny) of
#'    progeny indices relating to lines in genotype matrix}
#'   \item{sireHap}{list (LEN number of sires) of matrices (DIM 2 x p)
#'    of sire haplotypes (0, 1) on investigated chromosome}
#' }
#' @param genotype.chr matrix (DIM n x p) of all progeny genotypes (0, 1, 2) on
#'   a chromosome with p SNPs
#' @param snp.chr vector(LEN p) of SNP indices as in physical map
#' @param only.adj logical; if \code{TRUE}, recombination rate is calculated
#'   only between neighbouring markers
#' @param prec scalar; precision of estimation
#' @return list (LEN p - 1) of data.frames; for each SNP, parameters are
#'   estimated with all following SNPs; two solutions (prefix sln1 and sln2) are
#'   obtained for two runs of the EM algorithm
#' \describe{
#'   \item{\code{SNP1}}{index 1. SNP}
#'   \item{\code{SNP2}}{index 2. SNP}
#'   \item{\code{D}}{maternal LD}
#'   \item{\code{fAA}}{frequency of maternal haplotype 1-1}
#'   \item{\code{fAB}}{frequency of maternal haplotype 1-0}
#'   \item{\code{fBA}}{frequency of maternal haplotype 0-1}
#'   \item{\code{fBB}}{frequency of maternal haplotype 0-0}
#'   \item{\code{p1}}{Maternal allele frequency (allele 1)}
#'   \item{\code{p2}}{Maternal allele frequency (allele 0)}
#'   \item{\code{nfam1}}{size of genomic family 1}
#'   \item{\code{nfam2}}{size of genomic family 2}
#'   \item{\code{error}}{0 if computations were without error; 1 if EM algorithm
#'     did not converge}
#'   \item{\code{iteration}}{number of EM iterations}
#'   \item{\code{theta}}{paternal recombination rate}
#'   \item{\code{r2}}{\eqn{r^2} of maternal LD}
#'   \item{\code{logL}}{value of log likelihood function}
#'   \item{\code{unimodal}}{1 if likelihood is unimodal; 0 if likelihood is
#'     bimodal}
#'   \item{\code{critical}}{0 if parameter estimates are unique; 1 if parameter
#'     estimates at both solutions are valid, then decision process follows in
#'     post-processing function "editraw"}
#' }
#'   Afterwards, solutions are compared and processed with function
#'   \code{editraw}, yielding the final estimates for each valid pair of SNPs.
#' @examples
#'   ### test data
#'   data(targetregion)
#'   ### make list for paternal half-sib families
#'   hap <- makehaplist(daughterSire, hapSire)
#'   ### parameter estimates on a chromosome
#'   res <- hsrecombi(hap, genotype.chr, map.chr$SNP)
#'   ### post-processing to achieve final and valid set of estimates
#'   final <- editraw(res, map.chr)
#' @references
#'   Hampel, A., Teuscher, F., Gomez-Raya, L., Doschoris, M. & Wittenburg, D.
#'    (2018) Estimation of recombination rate and maternal linkage
#'    disequilibrium in half-sibs. Frontiers in Genetics 9:186.
#'    \url{https://doi.org/10.3389/fgene.2018.00186}
#'
#'   Gomez-Raya, L. (2012) Maximum likelihood estimation of linkage
#'    disequilibrium in half-sib families. Genetics 191:195-213.
#' @importFrom Rcpp evalCpp
#' @useDynLib hsrecombi
#' @export
hsrecombi <- function(hap, genotype.chr, snp.chr, only.adj = FALSE, prec = 1e-6){
  if((length(snp.chr) != ncol(genotype.chr)) | (length(snp.chr) != ncol(hap$sireHap[[1]]))) stop('ERROR: inconsistency in number of SNPs')
  ls <- list()
  for(j in 1:(length(snp.chr) - 1)){
    out <- c()
    if(j %% 100 == 1) message(paste('Processing SNP', j))
    q <- ifelse(!only.adj, length(snp.chr), j + 1)
    for (i in (j + 1):q){

      ## set-up genomic families for snp pairs
      GenFam1 <- GenFam2 <- c()
      for (l in 1:length(hap$sireHap)){

        if((hap$sireHap[[l]][1, i] != hap$sireHap[[l]][2, i]) && (hap$sireHap[[l]][1, j] != hap$sireHap[[l]][2, j])){
          if((hap$sireHap[[l]][1, j] == 0 & hap$sireHap[[l]][1, i] == 0) | (hap$sireHap[[l]][1, j] == 1 & hap$sireHap[[l]][1, i] == 1)){
            GenFam1 <- rbind(GenFam1, genotype.chr[hap$famID[[l]], c(j, i)])

          } else GenFam2 <- rbind(GenFam2, genotype.chr[hap$famID[[l]], c(j, i)])
        }
      }

      if(!is.null(GenFam1) | !is.null(GenFam2)){
        ## start values 1. run
        theta.start <- prec
        start <- startvalue(GenFam1, GenFam2, 0)
        ## 1.run
        sln1 <- LDHScpp(GenFam1, GenFam2, start$fAA.start, start$fAB.start, start$fBA.start, theta.start, F, prec)

        ## start values 2. run
        Ds <- sln1$D + (1 - 2 * start$p1) * (1 - 2 * start$p2) / 4
        theta.start <- (1 - 4 * Ds) / 2
        Dd <- (1 - 2 * sln1$theta) / 4 - (1 - 2 * start$p1) * (1 - 2 * start$p2) / 4
        fAA.start <- Dd + start$p1 * start$p2
        fAB.start <- -Dd + start$p1 * (1 - start$p2)
        fBA.start <- -Dd + (1 - start$p1) * start$p2

        ## 2. run if second mode exists
        if((!(theta.start > 1 | theta.start < 0)) & ((Dd >= start$L1) & (Dd <= start$L2))){
          sln2 <- LDHScpp(GenFam1, GenFam2, fAA.start, fAB.start, fBA.start, theta.start, F, prec)
          unimodal <- 0
        } else {
          sln2 <- sln1
          unimodal <- 1
        }

        out <- rbind(out, c(SNP1 = snp.chr[j], SNP2 = snp.chr[i], sln1 = unlist(sln1), sln2 = unlist(sln2),
                            unimodal = unimodal, critical = start$critical))
      }
    }
    ls[[j]] <- data.frame(out)
  }
  return(ls)
}



#' @title Estimation of genetic position
#' @name geneticPosition
#' @description Estimation of genetic positions (in centi Morgan)
#' @details Smoothing of recombination rates (theta) <= 0.05 via quadratic
#'   optimization provides an approximation of genetic distances (in Morgan)
#'   between SNPs. The cumulative sum * 100 yields the genetic positions in cM.
#'
#'   The minimization problem \code{(theta - D d)^2} is solved s.t. d > 0 where
#'   d is the vector of genetic distances between adjacent markers but theta is
#'   not restricted to adjacent markers. The incidence matrix D contains 1's for
#'   those intervals contributing to the total distance relevant for each theta.
#'
#'   Estmates of theta = 1e-6 are neglected as these values coincide with start
#'   values and indicate that (because of a very flat likelihood surface) no
#'   meaningful estimate of recombination rate has been obtained.
#' @param final table of results produced by \code{editraw} with pairwise
#'   estimates of recombination rate between p SNPs within chromosome; minimum
#'   required data frame with columns \code{SNP1}, \code{SNP2} and \code{theta}
#' @param exclude optional vector (LEN q) of SNPs to be excluded (e.g.,
#'   candidates of misplaced SNPs)
#' @param threshold optional value; recombination rates <= threshold are
#'   considered for smoothing
#' @return vector (LEN p) of genetic positions of SNPs (in cM)
#' @examples
#'   ### test data
#'   data(targetregion)
#'   ### make list for paternal half-sib families
#'   hap <- makehaplist(daughterSire, hapSire)
#'   ### parameter estimates on a chromosome
#'   res <- hsrecombi(hap, genotype.chr, map.chr$SNP)
#'   ### post-processing to achieve final and valid set of estimates
#'   final <- editraw(res, map.chr)
#'   ### approximation of genetic positions
#'   pos <- geneticPosition(final)
#' @importFrom quadprog solve.QP
#' @export
geneticPosition <- function(final, exclude = NULL, threshold = 0.05){
  part <- final[(final$theta <= threshold) & (final$theta >= 1e-5) &
                  !(final$SNP1 %in% exclude) & !(final$SNP2 %in% exclude), ]
  p <- max(part$SNP2) - min(part$SNP1) + 1

  # quadratic optimization min(theta - D d)² s.t. d > 0
  D <- matrix(0, ncol = p - 1, nrow = nrow(part))
  for (i in 1:nrow(part)) {
    D[i, (part$SNP1[i] - min(part$SNP1) + 1):(part$SNP2[i] - min(part$SNP1))] <- 1
  }

  # components for quadprog, make D p.d.
  dvec <- crossprod(D, part$theta)
  Amat <- diag(1, p - 1)
  Dmat <- crossprod(D) + diag(1e-6, p - 1)

  sln <- solve.QP(Dmat, dvec, Amat)
  gen <- c(0, cumsum(sln$solution)) * 100

  if(min(part$SNP1) > 1) gen <- c(rep(NA, min(part$SNP1) - 1), gen)
  names(gen) <- 1:max(part$SNP2)
  id <- unique(c(part$SNP1, part$SNP2))
  gen[!(names(gen) %in% id)] <- NA

  return(gen)
}

