#Function for the simulations
#Sample causal variants, compute haplotypes burden
#and thresholds corresponding to the desired prevalence
#p.causal=P(causal variant) ; p.protect=P(protective variant | causal variant)

rbm.haplos.thresholds <- function(haplos, weights = -0.4*log10(colMeans(haplos)), maf.threshold = 0.01, 
             nb.causal, p.protect = 0, h2, prev, normal.approx = TRUE, size, replicates, rep.by.causal, verbose = TRUE) {

  if( (replicates %% rep.by.causal) != 0 ) 
    stop("replicates should be a multiple of rep.by.causal")

  if(length(h2) != length(prev) | length(h2) != length(size) | length(prev) != length(size)) stop("h2 and prev should have same size as size")
  if(length(maf.threshold)==1){
    maf.threshold <- rep(maf.threshold, length(h2))
  }else{
    if(length(maf.threshold) != length(h2)) stop("maf.threshold should be of size 1 or of same size as h2 and prev")
  }
  if(length(nb.causal)==1){
    nb.causal <- rep(nb.causal, length(h2))
  }else{
    if(length(nb.causal) != length(h2)) stop("nb.causal should be of size 1 or of same size as h2 and prev")
  }
  if(length(p.protect)==1){
    p.protect <- rep(p.protect, length(h2))
  }else{
    if(length(p.protect) != length(h2)) stop("p.protect should be of size 1 or of same size as h2 and prev")
  }
  
  if(verbose) cat(length(h2), " groups of individuals will be simulated \n")
    
  
  x <- new.bed.matrix(nb_inds=sum(size), nb_snps=ncol(haplos)*replicates)
 
  for(i in 1:length(seq(1,replicates, by=rep.by.causal))){
    burdens <- mapply(get.causal.burden, weights=list(weights), haplos = list(haplos), maf.threshold, nb.causal, p.protect, h2, SIMPLIFY=FALSE)
  
    ##Computation of thresholds to respect 'prev' 
    ##normal.approx = TRUE : ok if small h2
    if(normal.approx) {
      s <- qnorm(prev, lower.tail = FALSE)
    } else {  
      BRD <- sample(burdens, 1e6, TRUE) + sample(burdens, 1e6, TRUE) + rnorm(1e6, sd = sqrt(1-h2))
      s <- quantile(BRD, 1 - prev)
    }

    .Call('rbm_haplos_thresholds_filling', PACKAGE = 'Ravages', x@bed,  haplos, burdens, sd = sqrt(1-h2), thr1 = s, thr2 = rep(Inf, length(size)), size, repNumber = i-1, reps=rep.by.causal)
  }

  x@ped$pheno <- factor(rep.int( 1:length(size) - 1, size))
  x@snps$genomic.region <- factor( rep( sprintf("R%0*d", log10(replicates) + 1, 1:replicates), each = ncol(haplos)) )
  
 
  if( is.null(colnames(haplos)) )
    x@snps$id <- paste( x@snps$genomic.region, x@snps$id, sep="_")
  else
    x@snps$id <- paste( x@snps$genomic.region, colnames(haplos), sep = "_")

  x <- set.stats(x, verbose = FALSE)
  x
}



get.causal.burden <- function(weights, haplos, maf.threshold, nb.causal, p.protect, h2){
  weights[colMeans(haplos) == 0 | colMeans(haplos) > maf.threshold ] <- 0
  w <- which(weights > 0) # parmi les SNPs potentiellement causaux
  
  if(nb.causal > length(w)) 
    stop("There are not enough positively weighted variants to pick ", nb.causal, " ones")

  nb.pro <- round(nb.causal * p.protect)
  nb.del <- nb.causal - nb.pro
  
  w.causal <- sample(w, nb.causal)
  directions <- numeric(ncol(haplos))
  directions[w.causal] <- c(rep(-1,nb.pro), rep(1,nb.del))

  ##Computation of burden (genetic component of liability)
  ##Based on the vector 'weights'
  we <- directions * weights
  burdens <- haplos %*% we
  #Standardisation of burdens
  burdens <- burdens - mean(burdens)

  ##Respect h2: adjust burdens' variance of the two haplotypes 
  tau <- as.numeric(2*var(burdens))
  burdens <- burdens / sqrt(tau) * sqrt(h2)
  #To generate the liability, a gaussian variable of standard deviation sqrt(1-h2) should be added to the burdens
  burdens
}
