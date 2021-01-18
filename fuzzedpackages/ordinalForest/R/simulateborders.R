simulateborders <-
  function(J, nsets=300, npermtrial=500, permperdefault=FALSE) {
    
    # Number of possible permutations of the ordering of the widths of the
    # J intervals corresponding to the J classes of the ordinal target
    # variable:
    nperm <- factorial(J)
    
    # If the number of possible permutations of the ordering
    # of the widths of the intervals 'nperm' is smaller than the number
    # of divisions, than consider all possible permutations.
    # Otherwise heterogeneous as possible rankings are considered:
    if((nperm < nsets) & !permperdefault) {
      # Matrix with all possible permutations:
      permmat <- combinat::permn(1:J)
      permmat <- matrix(nrow=length(permmat), ncol=J, data=unlist(permmat), byrow=TRUE)
      permmat <- permmat[sample(1:nrow(permmat)),]
      
      # Repeat individual permutation tupels until they have the number 'nsets':
      repind <- rep(1:nperm, each=floor(nsets/nperm))
      repind <- c(repind, sample(1:nperm, size=nsets - length(repind)))
      permind <- permmat[repind,]
    } 
    else {
      # Initialize matrix of rankings:
      permind <- matrix(nrow=nsets, ncol=J)
      
      # First ranking is random:
      permind[1,] <- sample(1:J)
      
      # Simulate other rankings:
      for(i in 2:nsets) {
        # Simulate 'npermtrial' (e.g. 500) permutations to try
        # for the i-th ranking:
        permindtemp <- t(replicate(npermtrial, sample(1:J)))
        # Use that of the 'npermtrial' rankings that is most
        # dissimilar to its precessor:
        dists <- apply(permindtemp, 1, function(x) mean((x-permind[i-1,])^2))
        permind[i,] <- permindtemp[which.max(dists),]
      }
    }
    
    # Simulate unordered borders:
    bordersbraw <- t(rbind(0, replicate(nsets, sort(runif(J-1))), 1))
    
    # Calculate and order differences between borders:
    diffmat <- t(apply(bordersbraw, 1, function(x) sort(diff(x))))
    
    # Order differences according to entries of the permutation matrix:
    permdiff <- t(sapply(1:nsets, function(x) diffmat[x,][permind[x,]]))
    
    # Sum up differences in order to obtained the borders:
    permdiffcumsum <- t(apply(permdiff, 1, cumsum))
    
    # Make border matrix:
    bordersb <- matrix(nrow=nsets, ncol=J+1)
    bordersb[,2:ncol(permdiffcumsum)] <- permdiffcumsum[,-ncol(permdiffcumsum)]
    bordersb[,1] <- 0
    bordersb[,J+1] <- 1
    
    # Return borders:
    bordersb
    
  }
