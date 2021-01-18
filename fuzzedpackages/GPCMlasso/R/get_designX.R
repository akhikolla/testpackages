get_designX <- function(X, DSF, m, I, q, n){
  if(!DSF){
    designX <- matrix(0, ncol = m * I, nrow = n * sum(q))
    
    pos_u <- 1
    for (u in 1:n) {
      for(uu in 1:I){
        designX[pos_u:(pos_u+q[uu]-1), ((uu - 1) * m + 1):(uu * m)] <-  
          matrix(rep(X[u,], q[uu]), byrow = TRUE, ncol = m, nrow = q[uu])
        pos_u <- pos_u+q[uu]
      }
    }
  }else{
    designX <- matrix(0, ncol = m * sum(q), nrow = n * sum(q))
    
    pos_u <- 1
    for (u in 1:n) {
      pos_uu <- 1
      for(uu in 1:I){
        for(uuu in 1:q[uu]){
          designX[pos_u, pos_uu:(pos_uu+m-1)] <-  
            X[u,]
          pos_u <- pos_u+1
          pos_uu <- pos_uu+m
        }
      }
    }
  }
  return(designX)
}


get_acoefs_old <- function(RSM, DSF, m, I, q, n_sigma){

  if(!DSF){
    pen1 <- diag(m*I)
  }else{
    pen1 <- matrix(0,nrow=m*sum(q),ncol=m*sum(choose(q,2)))
    pos1 <-1

    pos_pos <- 1
    for(u in 1:I){
      n_comb <- choose(q[u],2)
      if(n_comb>0){
      combis <- combn(q[u],2)-1
      for(uuu in 1:m){
        for(uu in 1:n_comb){
          pen1[combis[1,uu]*m+pos_pos,pos1] <- 1
          pen1[combis[2,uu]*m+pos_pos,pos1] <- -1
          pos1 <- pos1+1
        }
        pos_pos <- pos_pos+1
      }
      pos_pos <- pos_pos+(q[u]-1)*m
    }
    }
    pen1 <- cbind(diag(m*sum(q)),pen1)
  }
    if(RSM){
      acoefs <- rbind(matrix(0,nrow=q[1]+I-1,ncol=ncol(pen1)),pen1,
                      matrix(0, ncol = ncol(pen1), nrow = n_sigma))
      
    }else{
      acoefs <- rbind(matrix(0,nrow=sum(q),ncol=ncol(pen1)),pen1,
                      matrix(0, ncol = ncol(pen1), nrow = n_sigma))
    }
  }



get_acoefs <- function(RSM, DSF, m, I, q, n_sigma){

  if(!DSF){
    pen1 <- diag(m*I)
  }else{
    pen1 <- matrix(0,nrow=m*sum(q),ncol=m*sum(q-1))
    pos1 <-1
    
    pos_pos <- 1
    for(u in 1:I){
      n_comb <- q[u] - 1
      if(n_comb>0){
        combis <- rep(1:q[u], each = 2)
        combis <- matrix(combis[-c(1,length(combis))]-1, nrow = 2)
        for(uuu in 1:m){
          for(uu in 1:n_comb){
            pen1[combis[1,uu]*m+pos_pos,pos1] <- 1
            pen1[combis[2,uu]*m+pos_pos,pos1] <- -1
            pos1 <- pos1+1
          }
          pos_pos <- pos_pos+1
        }
        pos_pos <- pos_pos+(q[u]-1)*m
      }
    }
    pen1 <- cbind(diag(m*sum(q)),pen1)
  }
  if(RSM){
    acoefs <- rbind(matrix(0,nrow=q[1]+I-1,ncol=ncol(pen1)),pen1,
                    matrix(0, ncol = ncol(pen1), nrow = n_sigma))
    
  }else{
    acoefs <- rbind(matrix(0,nrow=sum(q),ncol=ncol(pen1)),pen1,
                    matrix(0, ncol = ncol(pen1), nrow = n_sigma))
  }
}


