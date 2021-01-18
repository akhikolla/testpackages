permute_select_surv <- function(xdata, Y.names, P.names, T.names, importance = 'IIS', method = "R2", Bag, args.rpart, args.parallel = list(numWorkers = 1), nperm = 50)
{
  survbag <- Bagg_Surv(xdata, Y.names, P.names, T.names, method, args.rpart, args.parallel, Bag) ;
  IMPH1 <- rep(0, length(T.names)) ;
  names(IMPH1) <- T.names ;
  IMPH1[names(survbag[[importance]])] <- survbag[[importance]] ;
  PERMH0 <- matrix(0, nrow = length(T.names), ncol = nperm ) ;
  for(i in 1:nperm)
  {
    xdataH0 <- xdata ;
    xdataH0[, c(Y.names, P.names)] = xdata[, c(Y.names, P.names)][sample(nrow(xdata)), ]
    survbagH0 <- Bagg_Surv(xdataH0, Y.names, P.names, T.names, method, args.rpart, args.parallel, Bag) ;
    tempimp <- rep(0, length(T.names)) ;
    names(tempimp) <- T.names ;
    impH0 <- survbagH0[[importance]] ; 
    tempimp[names(impH0)] <- impH0 ;
    PERMH0[, i] = tempimp ;
  }
  
  pvalperm1 <- sapply(1:length(T.names), function(u) mean(PERMH0[u,] > IMPH1[u])) ;
  names(pvalperm1) <- T.names ;
  pvalperm1 <- sort(pvalperm1) ;
  
  
  ## normal law adequation under the null hypothesis
  meansH0 <- apply(PERMH0, 1, mean) ;
  sddH0 <- apply(PERMH0, 1, sd) ;
  sddH0[sddH0 <= mean(sddH0)]  <- mean(sddH0) ;
  pvalperm2 <- sort(sapply(1:length(T.names), function(uu) 1-pnorm(IMPH1[uu], meansH0[uu], sddH0[uu])))
  
  pvalKS <- sapply(1:length(T.names), function(u) ks.test(PERMH0[u,], rnorm(nperm, meansH0[u], sddH0[u]))$p.value) ;
  
  
  return(list(pvalperm1 = pvalperm1, pvalperm2 = pvalperm2, pvalKS = pvalKS, IMPH1 = IMPH1, PERMH0 = PERMH0)) ;
}
