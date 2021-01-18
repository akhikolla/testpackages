#' @title Make list of sire haplotypes
#' @name makehaplist
#' @description List of sire haplotypes is set up in the format required for
#'   hsrecombi. Haplotypes (obtained by external software) are provided.
#' @param daughterSire vector (LEN n) of sire ID for each progeny
#' @param hapSire matrix (DIM 2N x p + 1) of sire haplotype at p SNPs; 2 lines
#'  per sire, 1. columns contains sire ID
#' @param nmin scalar, minimum number of progeny required, default 1
#' @return hap list (LEN 2) of lists. For each sire:
#' \describe{
#'   \item{\code{famID}}{list (LEN N) of vectors (LEN n.progeny) of progeny
#'    indices relating to lines in genotype matrix}
#'   \item{\code{sireHap}}{list (LEN N) of matrices (DIM 2 x p) of sire
#'    haplotypes (0, 1) on investigated chromosome}
#' }
#' @examples
#'   data(targetregion)
#'   hap <- makehaplist(daughterSire, hapSire)
#' @export
makehaplist <- function(daughterSire, hapSire, nmin = 1){

  sireID <- unique(hapSire[, 1])
  ListFam <- ListHap <- list()

  for (i in sireID){
    # index of progeny
    hsID <- which(daughterSire == i)

    # filter only sufficiently large half-sib families
    if(length(hsID) >= nmin){
      ListFam[[as.character(i)]] <- hsID
      ListHap[[as.character(i)]] <- hapSire[hapSire[, 1] == i, ][, -1]
    }
  }
  return(list(famID = ListFam, sireHap = ListHap))
}


#' @title Make list of imputed sire haplotypes
#' @name makehap
#' @description List of sire haplotypes is set up in the format required for
#'   hsrecombi. Sire haplotypes are imputed from progeny genotypes using R
#'   package \code{hsphase}.
#' @param sireID vector (LEN N) of IDs of all sires
#' @param daughterSire vector (LEN n) of sire ID for each progeny
#' @param genotype.chr matrix (DIM n x p) of progeny genotypes on a single
#'   chromosome with p SNPs
#' @param nmin scalar, minimum required number of progeny for proper imputation,
#'   default 30
#' @return list (LEN 2) of lists. For each sire:
#' \describe{
#'   \item{\code{famID}}{list (LEN N) of vectors (LEN n.progeny) of progeny
#'    indices relating to lines in genotype matrix}
#'   \item{\code{sireHap}}{list (LEN N) of matrices (DIM 2 x p) of sire
#'    haplotypes (0, 1) on investigated chromosome}
#' }
#' @examples
#'  data(targetregion)
#'  hap <- makehap(unique(daughterSire), daughterSire, genotype.chr)
#' @references
#'  Ferdosi, M., Kinghorn, B., van der Werf, J., Lee, S. & Gondro, C. (2014)
#'   hsphase: an R package for pedigree reconstruction, detection of
#'   recombination events, phasing and imputation of half-sib family groups
#'   BMC Bioinformatics 15:172.
#'   \url{https://CRAN.R-project.org/package=hsphase}
#' @import hsphase
#' @export
makehap <- function(sireID, daughterSire, genotype.chr, nmin = 30){

  ListFam <- ListHap <- list()

  for (i in sireID){
    # index of progeny
    hsID <- which(daughterSire == i)

    # filter only sufficiently large half-sib families
    if(length(hsID) >= nmin){
      genotype.chr.fam <- genotype.chr[hsID, ]
      ListFam[[as.character(i)]] <- hsID

      b <- NULL
      try(b <- hsphase::bmh(genotype.chr.fam), silent = T)

      if(is.null(b)){
        ListHap[[as.character(i)]] <- matrix(9, nrow = 2, ncol = ncol(genotype.chr.fam))
      } else{
        ListHap[[as.character(i)]] <- hsphase::ssp(b, genotype.chr.fam)
      }

    }
  }
  return(list(famID = ListFam, sireHap = ListHap))
}


#' @title Make list of imputed haplotypes and recombination rate
#' @name makehappm
#' @description List of sire haplotypes is set up in the format required for
#'   hsrecombi. Sire haplotypes are imputed from progeny genotypes using R
#'   package \code{hsphase}. Furthermore, recombination rate estimates between
#'   adjacent SNPs from hsphase are reported.
#' @param sireID vector (LEN N) of IDs of all sires
#' @param daughterSire vector (LEN n) of sire ID for each progeny
#' @param genotype.chr matrix (DIM n x p) of progeny genotypes on a single
#'   chromosome with p SNPs
#' @param nmin scalar, minimum number of progeny required, default 1
#' @return hap list (LEN 2) of lists. For each sire:
#' \describe{
#'   \item{\code{famID}}{list (LEN N) of vectors (LEN n.progeny) of progeny
#'    indices relating to lines in genotype matrix}
#'   \item{\code{sireHap}}{list (LEN N) of matrices (DIM 2 x p) of sire
#'    haplotypes (0, 1) on investigated chromosome}
#'  \item{probRec}{vector (LEN p - 1) of proportion of recombinant progeny over
#'    all families between adjacent SNPs}
#'  \item{numberRec}{list (LEN N) of vectors (LEN n.progeny) of number of
#'    recombination events per animal}
#' }
#' @examples
#'   data(targetregion)
#'   hap <- makehappm(unique(daughterSire), daughterSire, genotype.chr)
#' @references
#'  Ferdosi, M., Kinghorn, B., van der Werf, J., Lee, S. & Gondro, C. (2014)
#'   hsphase: an R package for pedigree reconstruction, detection of
#'   recombination events, phasing and imputation of half-sib family groups
#'   BMC Bioinformatics 15:172.
#'   \url{https://CRAN.R-project.org/package=hsphase}
#' @import hsphase
#' @importFrom rlist list.rbind
#' @export
makehappm <- function(sireID, daughterSire, genotype.chr, nmin = 30){

  ListFam <- ListHap <- ListPm <- ListRec <- list()

  for (i in sireID){
    # index of progeny
    hsID <- which(daughterSire == i)

    # filter only sufficiently large half-sib families
    if(length(hsID) >= nmin){
      genotype.chr.fam <- genotype.chr[hsID, ]
      ListFam[[as.character(i)]] <- hsID

      b <- NULL
      try(b <- hsphase::bmh(genotype.chr.fam), silent = T)

      if(is.null(b)){
        ListHap[[as.character(i)]] <- matrix(9, nrow = 2, ncol = ncol(genotype.chr.fam))
        ListPm[[as.character(i)]] <- matrix(NA, nrow = 1, ncol = ncol(genotype.chr.fam) - 1)
        ListRec[[as.character(i)]] <- rep(NA, length(hsID))
      } else{
        ListHap[[as.character(i)]] <- hsphase::ssp(b, genotype.chr.fam)
        ListPm[[as.character(i)]] <- hsphase::pm(b, method = 'relative')
        ListRec[[as.character(i)]] <- rowSums(ListPm[[as.character(i)]], na.rm = TRUE)
      }
    }
    big <- rlist::list.rbind(ListPm)

  }
  return(list(famID = ListFam, sireHap = ListHap, probRec = colMeans(big, na.rm = TRUE), numberRec = ListRec))
}



#' @title Start value for maternal allele and haplotype frequencies
#' @name startvalue
#' @description Determine default start values for Expectation Maximisation (EM)
#'  algorithm that is used to estimate paternal recombination rate and maternal
#'  haplotype frequencies
#' @param Fam1 matrix (DIM n.progeny x 2) of progeny genotypes of genomic family
#'   with coupling phase sires (1) at SNP pair
#' @param Fam2 matrix (DIM n.progeny x 2) of progeny genotypes of genomic family
#'   with repulsion phase sires (2) at SNP pair
#' @param Dd maternal LD, default 0
#' @param prec minimum accepted start value for fAA, fAB, fBA; default
#'  \code{1e-6}
#' @return list (LEN 8)
#' \describe{
#'  \item{\code{fAA.start}}{frequency of maternal haplotype 1-1}
#'  \item{\code{fAB.start}}{frequency of maternal haplotype 1-0}
#'  \item{\code{fBA.start}}{frequency of maternal haplotype 0-1}
#'  \item{\code{p1}}{estimate of maternal allele frequency (allele 1) when sire
#'   is heterozygous at \code{SNP1}}
#'  \item{\code{p2}}{estimate of maternal allele frequency (allele 1) when sire
#'   is heterozygous at \code{SNP2}}
#'  \item{\code{L1}}{lower bound of maternal LD}
#'  \item{\code{L2}}{upper bound for maternal LD}
#'  \item{\code{critical}}{0 if parameter estimates are unique; 1 if parameter
#'   estimates at both solutions are valid}
#' }
#' @examples
#'  n1 <- 100
#'  n2 <- 20
#'  G1 <- matrix(ncol = 2, nrow = n1, sample(c(0:2), replace = TRUE,
#'   size = 2 * n1))
#'  G2 <- matrix(ncol = 2, nrow = n2, sample(c(0:2), replace = TRUE,
#'   size = 2 * n2))
#'  startvalue(G1, G2)
#' @export
startvalue <- function(Fam1, Fam2, Dd = 0, prec = 1e-6){
  allFam <- rbind(Fam1, Fam2)

  ## allele frequency first locus
  n2 <- sum(allFam[, 1] == 2)
  n0 <- sum(allFam[, 1] == 0)
  if(n2 + n0 > 0) p1 <- n2 / (n2 + n0) else p1 <- 0.5

  ## allele frequency second locus
  n2 <- sum(allFam[, 2] == 2)
  n0 <- sum(allFam[, 2] == 0)
  if(n2 + n0 > 0) p2 <- n2 / (n2 + n0) else p2 <- 0.5

  ## limits for Ddam
  L1 <- max(-p1 * p2, -(1 - p1) * (1 - p2))
  L2 <- min(p1 * (1 - p2), (1 - p1) * p2)

  crit = ifelse((p1 > 0.52) | (p2 > 0.52) | (p1 < 0.48) | (p2 < 0.48), 0, 1)
  fAA <- max(Dd + p1 * p2 - prec, prec)
  fAB <- max(-Dd + p1 * (1 - p2) - prec, prec)
  fBA <- max(-Dd + (1 - p1) * p2 - prec, prec)

  return(list(fAA.start =  fAA, fAB.start = fAB, fBA.start = fBA,
              p1 = p1, p2 = p2, L1 = L1, L2 = L2, critical = crit))
}


#' @title Editing results of hsrecombi
#' @name editraw
#' @description Process raw results from \code{hsrecombi}, decide which out of
#'   two sets of estimates is more likely and prepare list of final results
#' @param Roh list of raw results from \code{hsrecombi}
#' @param map1 data.frame containing information on physical map, at least:
#' \describe{
#'  \item{\code{SNP}}{SNP ID}
#'  \item{\code{locus_Mb}}{physical position in Mbp of SNP on chromosomes}
#'  \item{\code{Chr}}{chromosome of SNP}
#' }
#' @return final table of results
#' \describe{
#'   \item{\code{SNP1}}{index 1. SNP}
#'   \item{\code{SNP2}}{index 2. SNP}
#'   \item{\code{D}}{maternal LD}
#'   \item{\code{fAA}}{frequency of maternal haplotype 1-1}
#'   \item{\code{fAB}}{frequency of maternal haplotype 1-0}
#'   \item{\code{fBA}}{frequency of maternal haplotype 0-1}
#'   \item{\code{fBB}}{frequency of maternal haplotype 0-0}
#'   \item{\code{p1}}{Maternal allele frequency (allele 1) SNP1}
#'   \item{\code{p2}}{Maternal allele frequency (allele 1) SNP2}
#'   \item{\code{nfam1}}{size of genomic family 1}
#'   \item{\code{nfam2}}{size of genomic family 2}
#'   \item{\code{error}}{0 if computations were without error; 1 if EM algorithm
#'     did not converge}
#'   \item{\code{iteration}}{number of EM iterations}
#'   \item{\code{theta}}{paternal recombination rate}
#'   \item{\code{r2}}{\eqn{r^2} of maternal LD}
#'   \item{\code{logL}}{value of log likelihood function}
#'   \item{\code{unimodal}}{1 if likelihood is unimodal; 0 if likelihood is
#'    bimodal}
#'   \item{\code{critical}}{0 if parameter estimates were unique; 1 if parameter
#'     estimates were obtained via decision process}
#'   \item{\code{locus_Mb}}{physical distance between SNPs in Mbp}
#' }
#' @examples
#'   ### test data
#'   data(targetregion)
#'   ### make list for paternal half-sib families
#'   hap <- makehaplist(daughterSire, hapSire)
#'   ### parameter estimates on a chromosome
#'   res <- hsrecombi(hap, genotype.chr, map.chr$SNP)
#'   ### post-processing to achieve final and valid set of estimates
#'   final <- editraw(res, map.chr)
#' @importFrom dplyr select starts_with
#' @importFrom data.table rbindlist
#' @export
editraw <- function(Roh, map1){

  jointset <- rbindlist(Roh)

  ## 1: final estimate for non-critical SNPs based on loglik
  idx1 <- which(jointset$critical == 0)

  ## idx1.2 -> 2. set is more likely
  idx1.2 <- (jointset$unimodal[idx1] == 0) & (jointset$sln1.logL[idx1] < jointset$sln2.logL[idx1])
  idx1.1 <- !idx1.2

  part1 <- dplyr::select(jointset, dplyr::starts_with('sln1'))
  part2 <- dplyr::select(jointset, dplyr::starts_with('sln2'))
  text <- unlist(lapply(colnames(part1), function(x){strsplit(x, '.', fixed = T)[[1]][2]}))
  colnames(part1) <- colnames(part2) <- text

  binder <- rbind(cbind(SNP1 = jointset$SNP1[idx1[idx1.1]], SNP2 = jointset$SNP2[idx1[idx1.1]],
                        unimodal = jointset$unimodal[idx1[idx1.1]], critical = jointset$critical[idx1[idx1.1]],
                        part1[idx1[idx1.1], ]),
                  cbind(SNP1 = jointset$SNP1[idx1[idx1.2]], SNP2 = jointset$SNP2[idx1[idx1.2]],
                        unimodal = jointset$unimodal[idx1[idx1.2]], critical = jointset$critical[idx1[idx1.2]],
                        part2[idx1[idx1.2], ]))


  ## 2: final estimates for crital SNPs based on smoothing spline for non-critical SNPs
  idx2 <- which(jointset$critical != 0)
  crit <- unique(jointset$SNP1[idx2])
  jointset$sp <- NA

  if(length(crit) > 0){
    ## relative physical position
    map1$d <- map1$locus_Mb / max(map1$locus_Mb)

    for(j in 1:length(crit)){
      ## verified information of non-critical SNP pairs
      vec1 <- which(binder$SNP1 == crit[j])
      vec2 <- which(binder$SNP2 == crit[j])
      snp1 <- binder$SNP2[vec1]
      snp2 <- binder$SNP1[vec2]
      d1 <- map1$d[match(snp1, map1$SNP)] - map1$d[match(crit[j], map1$SNP)]
      d2 <- map1$d[match(snp2, map1$SNP)] - map1$d[match(crit[j], map1$SNP)]
      sp <- stats::smooth.spline(x = c(d1, d2), y = binder$theta[c(vec1, vec2)], df = 5)

      snp2.id <- which(jointset$SNP1[idx2] == crit[j])
      for(i in snp2.id){
        d <- map1$d[match(jointset$SNP2[idx2[i]], map1$SNP)] - map1$d[match(crit[j], map1$SNP)]
        est <- stats::predict(sp, x = d)$y
        jointset$sp[idx2[i]] <- ifelse((est - jointset$sln1.theta[idx2[i]]) ^ 2 <= (est - jointset$sln2.theta[idx2[i]]) ^ 2, T, F)
      }
    }
  }

  ## idx2.1 -> 1. set is closer to the smoothing curve
  idx2.1 <- jointset$sp[idx2]
  idx2.2 <- !idx2.1

  binder2 <- rbind(cbind(SNP1 = jointset$SNP1[idx2[idx2.1]], SNP2 = jointset$SNP2[idx2[idx2.1]],
                        unimodal = jointset$unimodal[idx2[idx2.1]], critical = jointset$critical[idx2[idx2.1]],
                        part1[idx2[idx2.1], ]),
                  cbind(SNP1 = jointset$SNP1[idx2[idx2.2]], SNP2 = jointset$SNP2[idx2[idx2.2]],
                        unimodal = jointset$unimodal[idx2[idx2.2]], critical = jointset$critical[idx2[idx2.2]],
                        part2[idx2[idx2.2], ]))
  binder <- rbind(binder, binder2)

  ## sort by SNP1/SNP2
  binder <- binder[with(binder, order(SNP1, SNP2)), ]
  binder$dist_Mb <- map1$locus_Mb[match(binder$SNP2, map1$SNP)] - map1$locus_Mb[match(binder$SNP1, map1$SNP)]

  return(binder)
}


#' @title Candidates for misplacement
#' @name checkCandidates
#' @description Search for SNPs with unusually large estimates of recombination
#'   rate
#' @details Markers with unusually large estimates of recombination rate to
#'   close SNPs are candidates for misplacements in the underlying assembly. The
#'   mean of recombination rate estimates with \code{win} subsequent or
#'   preceeding markers is calculated and those SNPs with mean value exceeding
#'   the \code{quant} quantile are denoted as candidates which have to be
#'   manually curated!
#'   This can be done, for instance, by visual inspection of a correlation plot
#'   containing estimates of recombination rate in a selected region.
#' @param final table of results produced by \code{editraw} with pairwise
#'   estimates of recombination rate between p SNPs within chromosome; minimum
#'   required data frame with columns \code{SNP1}, \code{SNP2} and \code{theta}
#' @param win optional value for window size; default value 30
#' @param quant optional value; default value 0.99, see details
#' @return vector of SNP indices for further verification
#' @examples
#'   ### test data
#'   data(targetregion)
#'   ### make list for paternal half-sib families
#'   hap <- makehaplist(daughterSire, hapSire)
#'   ### parameter estimates on a chromosome
#'   res <- hsrecombi(hap, genotype.chr, map.chr$SNP)
#'   ### post-processing to achieve final and valid set of estimates
#'   final <- editraw(res, map.chr)
#'   ### check for candidates of misplacement
#'   snp <- checkCandidates(final)
#' @importFrom stats quantile
#' @export
checkCandidates <- function(final, win = 30, quant = 0.99){
  p <- max(final$SNP2)
  meanval <- c()
  for (i in min(final$SNP1):(p - win)){
    meanval[i] <- mean(final$theta[(final$SNP1 == i) & (final$SNP2 <= i + win)])
  }

  for (i in (p - win + 1):p) {
    meanval[i] <- mean(final$theta[(final$SNP1 >= i - win) & (final$SNP2 == i)])
  }

  q <- quantile(meanval, prob = quant, na.rm = T)
  cand <- which(meanval >= q)

  return(cand)
}
