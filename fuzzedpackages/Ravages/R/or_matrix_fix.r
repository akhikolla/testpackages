OR.matrix <- function (n.variants, n.groups, OR.del, OR.pro = 1/OR.del, p.causal, prob.pro) {
    if (length(OR.del) != length(OR.pro))
        stop("Dimensions mismatch")
    
    if(!is.matrix(OR.del)){
      if(length(OR.del)==n.variants){
		OR.del <- matrix(OR.del, nrow=1)
		OR.pro <- matrix(OR.pro, nrow=1)
	}
	else if(length(OR.del)==n.groups){
		OR.del <- matrix( rep_len(OR.del, n.variants*length(OR.del)), nrow=n.groups)
		OR.pro <- matrix( rep_len(OR.pro, n.variants*length(OR.pro)), nrow=n.groups)
	}
	else
		stop("OR Dimensions mismatch")
    }

   nb.causal <- n.variants*p.causal
   #Sample of causal variants
   v.causal <- lapply(1:nrow(OR.del), function(z) if(nb.causal>0) sample(1:n.variants, nb.causal))
   v.protect <- lapply(v.causal, function(z) if(prob.pro>0) sample(z, nb.causal*prob.pro))

   #Matrix of 1 to be change with GRR values of causal variants	
   OR.tot <- matrix(rep(1, n.variants*nrow(OR.del)), nrow=nrow(OR.del))
   
   OR.tot <- t(sapply(1:nrow(OR.tot), function(z){ OR.tot[z,v.causal[[z]]] <- OR.del[z, v.causal[[z]]] ; return(OR.tot[z,])}))
   OR.tot <- t(sapply(1:nrow(OR.tot), function(z){ OR.tot[z,v.protect[[z]]] <- OR.pro[z, v.protect[[z]]] ; return(OR.tot[z,])})) 
   return(OR.tot)
}


