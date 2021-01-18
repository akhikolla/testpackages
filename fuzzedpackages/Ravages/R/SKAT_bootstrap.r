SKAT.bootstrap <- function(x, NullObject, genomic.region = x@snps$genomic.region, 
                 weights = (1-x@snps$maf)**24, maf.threshold = 0.5, 
                 perm.target = 100, perm.max = 5e4, debug = FALSE,
                 estimation.pvalue = "kurtosis") {

  if(!is.factor(genomic.region)) stop("'genomic.region' should be a factor")
  genomic.region <- droplevels(genomic.region)

  which.snps <- (x@snps$maf <= maf.threshold) & (x@snps$maf > 0)

  group <- NullObject$group
  Pi <- NullObject$Pi.data
  X <- NullObject$X 

  # matrice U = Id - WX (X'WX)^{-1} X'W pour le bootstrap
  W <- W.mat(Pi[,-1,drop=FALSE])
  XX <- block.diag( rep(list(X), ncol(Pi) - 1) )
  WX <- W %*% XX
  XWX <- t(XX) %*% WX
  U <- -(WX %*% solve(XWX, t(XX)))
  diag(U) <- diag(U)+1

  B <- .Call('skat_bootstrap', PACKAGE = "Ravages", x@bed, which.snps, genomic.region, group, x@p, Pi, weights, U, perm.target, perm.max);

  names(B)[5] <- "p.perm"

  # pour la procédure d'approx par un chi2
  M1 <- B$M1;
  M2 <- B$M2;
  M3 <- B$M3;
  M4 <- B$M4;
  # retrouver moments centrés
  S2 <- M2 - M1**2 # variance
  m3 <- M3 - 3*S2*M1 - M1**3 # 3e moment
  m4 <- M4 - 4*m3*M1 - 6*S2*M1**2 - M1**4 # 4e moment
  B$sigma  <- sqrt(S2)
  B$skewness <- m3/S2**1.5
  B$kurtosis <- m4/S2**2
  
  B$p.chi2 <- as.vector(mapply(p.valeur.moments.liu, Q = B$stat, mu = B$M1, sigma = B$sigma, skewness = B$skewness, kurtosis = B$kurtosis - 3, estimation.pvalue = estimation.pvalue))

  names(B)[6] <- "mean"
  B$p.value <- ifelse(B$nb.perm < perm.max, B$p.perm, B$p.chi2) 
  #If p.chi2 is NA, return NA in p.value
  B$p.value <- ifelse(is.na(B$p.chi), NA, B$p.value)
 
  B <- as.data.frame(B, row.names = levels(genomic.region))
  if(debug) 
    B[,!(names(B) %in% c("M2", "M3", "M4"))]
  else
    B[ , c("stat", "p.perm", "p.chi2", "p.value") ]
}

