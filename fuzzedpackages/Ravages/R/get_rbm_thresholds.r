#Function for the simulations
#Sample causal variants, compute haplotypes burden
#and thresholds corresponding to the desired prevalence
#p.causal=P(causal variant) ; p.protect=P(protective variant | causal variant)
simus.parameters <- function(haplos, weights = -0.4*log10(colMeans(haplos)), maf.threshold = 0.01, nb.causal, p.protect, h2, prev, normal.approx) {
  weights[colMeans(haplos) == 0 | colMeans(haplos) > maf.threshold ] <- 0
  w <- which(weights > 0) # parmi les SNPs potentiellement causaux
  if(nb.causal > length(w)) 
    stop("There are not enough positively weighted variants to pick", nb.causal, "ones")

  nb.pro <- round(nb.causal * p.protect)
  nb.del <- nb.causal - nb.pro

  
  w.causal <- sample(w, nb.causal)
  directions <- numeric(length(weights))
  directions[w.causal] <- c(rep(-1,nb.pro), rep(1,nb.del))


  ## calcul des fardeaux (composante gÃ©nÃ©tique de la liabilitÃ©) 
  ## Ã  partir des poids du vecteur 'weights'
  we <- directions * weights
  burdens <- haplos %*% we
  # on centre les fardeaux
  burdens <- burdens - mean(burdens)

  # respecter h2 :
  # on ajuste la variance des fardeaux des deux haplotypes
  tau <- as.numeric(2*var(burdens))
  burdens <- burdens / sqrt(tau) * sqrt(h2)
  # pour gÃ©nÃ©rer la liabilitÃ©, il faudra ajouter aux burdens gÃ©notypiques
  # une v.a. normale d'e.t. sqrt(1-h2)


  ## calcul des seuils pour respecter les prÃ©valences donnÃ©es dans 'prev'
  ## si h2 est petit on peut utiliser normal.approx = TRUE sans problÃ¨me
  if(normal.approx) {
    s <- qnorm(prev, lower.tail = FALSE)
  } else {  
    BRD <- sample(burdens, 1e6, TRUE) + sample(burdens, 1e6, TRUE) + rnorm(1e6, sd = sqrt(1-h2))
    s <- quantile(BRD, 1 - prev)
  }

  list(burdens = burdens, sd = sqrt(1-h2), thr1 = s)
}


get.rbm.thresholds <- function(haplos, weights = -0.4*log10(colMeans(haplos)), maf.threshold = 0.01, nb.causal, p.protect, h2, prev, normal.approx = TRUE){
  RR <- mapply(simus.parameters, list(haplos), list(weights), maf.threshold, nb.causal, p.protect, h2, prev, normal.approx, SIMPLIFY=FALSE)
  list(burdens = lapply(RR, function(z) z$burdens), sd = unlist(lapply(RR, function(z) z$sd)), thr1 = unlist(lapply(RR, function(z) z$thr1)) )
}