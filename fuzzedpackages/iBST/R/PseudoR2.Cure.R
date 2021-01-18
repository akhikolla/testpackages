PseudoR2.Cure <-
function(ygene, ydelai, yetat, strate, ordered = FALSE) {
  n <- length(ydelai) ;
  
  if(is.vector(strate)) strate <- matrix(strate, ncol = 1)
  
  surv.object <- Surv(ydelai,yetat)
  n_times_events <- length(unique( ydelai[ yetat == max(yetat) ] ))
  
  alpha <- if(is.null(strate)) 0 else coxph(surv.object ~ strate)$coeff 
  
  if(!ordered) { # on ordonne tout
    o <- order(ydelai) ;
    strate <- strate[o,,drop=FALSE] 
    ydelai <- ydelai[o]
    yetat  <- yetat[o]
    ygene  <- ygene[o]
  }
  # cureE = exp( beta Z ) 
  cureE <- if(is.null(strate)) rep(1,n) else exp(strate %*% alpha )
  
  .Call("iBST_COX", ydelai, cureE, ygene, yetat, PACKAGE = "iBST")/n_times_events
  # COX(ydelai, cureE, ygene, yetat)/n_times_events
  
  
}
