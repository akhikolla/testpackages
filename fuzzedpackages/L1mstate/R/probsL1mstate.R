##############################################################
###Compute subject-specific or overall transition probabilities
###based on 'probtrans' function of 'mstate' package

probs.l1mstate = function(object, longdt, tmat, predt, direction=c("forward","fixedhorizon")){
  
  transit = tmat2(tmat)
  numtrans = nrow(transit)
  
  #process to get cumHaz for all unique time points from all trans
  x1 = object$Haz
  unt = sort(unique(x1$time))
  K = unique(x1$trans)
  k=min(K)
  xk = x1[x1$trans==k,]
  xk_new = matrix(rep(0,3*length(unt)), nrow = length(unt))
  xk_new[,1] = unt
  xk_new[,2] = ifelse(xk_new[,1] %in% xk$time, xk$Haz, 0)
  xk_new[,2] = cumsum(xk_new[,2])
  xk_new[,3] = rep(k, length(unt))
  x1_new = xk_new
  for(k in K){
    if(k!=min(K)){
      xk = x1[x1$trans==k,]
      xk_new = matrix(rep(0,3*length(unt)), nrow = length(unt))
      xk_new[,1] = unt
      xk_new[,2] = ifelse(xk_new[,1] %in% xk$time, xk$Haz, 0)
      xk_new[,2] = cumsum(xk_new[,2])
      xk_new[,3] = rep(k, length(unt))
      x1_new = rbind(x1_new, xk_new)
    }
  }
  x1_new = data.frame(x1_new)
  names(x1_new)=c("time", "cumHaz", "trans")
  stackhaz = x1_new
  
  
  for (i in 1:numtrans)
    stackhaz$dhaz[stackhaz$trans==i] = diff(c(0,stackhaz$cumHaz[stackhaz$trans==i]))
  if (direction=="forward"){
    stackhaz = stackhaz[stackhaz$time > predt,]
  }else{
    stackhaz = stackhaz[stackhaz$time <= predt,]
  } 
  
  untimes = sort(unique(stackhaz$time))
  TT = length(untimes)
  S = nrow(tmat)
  
  if (direction=="forward") {
    res = array(0,c(TT+1,S+1,S)) # S+1 for time, probs (S)
    # first line (in case of forward) contains values for t=predt
    res[1,1,] = predt
    for (j in 1:S) res[1, 1+j,] = rep(c(0,1,0), c(j-1,1,S-j))
  }else{
    # situation for backward is different from forward,
    # depends on whether predt is an event time
    if (predt %in% untimes) {
      res = array(0,c(TT+1,S+1,S))
      res[TT+1,1,] = predt
      for (j in 1:S) res[TT+1, 1+j,] = rep(c(0,1,0), c(j-1,1,S-j))
    }else{
      res = array(0,c(TT+2,S+1,S))
      res[TT+1,1,] = max(untimes)
      for (j in 1:S) res[TT+1, 1+j,] = rep(c(0,1,0), c(j-1,1,S-j))
      
      res[TT+2,1,] = predt
      for (j in 1:S) res[TT+2, 1+j,] = rep(c(0,1,0), c(j-1,1,S-j))
    }
  }
  
  P = diag(S)
  
  for (i in 1:TT)
  {
    idx = ifelse(direction=="forward",i,TT+1-i)
    tt = untimes[idx]
    Haztt = stackhaz[stackhaz$time==tt,]
    lHaztt = nrow(Haztt)
    # build S x S matrix IplusdA
    IplusdA = diag(S)
    for (j in 1:lHaztt){
      from = transit$from[transit$transno==Haztt$trans[j]]
      to = transit$to[transit$transno==Haztt$trans[j]]
      IplusdA[from, to] = Haztt$dhaz[j]
      IplusdA[from, from] = IplusdA[from, from] - Haztt$dhaz[j]
    }
    
    if (any(diag(IplusdA)<0))
      warning("Warning! Negative diagonal elements of (I+dA); the estimate may not be meaningful. \n")
    
    if (direction=="forward") 
    {  
      P = P %*% IplusdA
    }else{
      P = IplusdA %*% P
    }
    
    if (direction=="forward") {
      res[idx+1,1,] = tt
      res[idx+1,2:(S+1),] = t(P)
    }else{
      res[idx,1,] = ifelse(i==TT,0,untimes[TT-i])
      res[idx,2:(S+1),] = t(P)
    }
  }
  
  ### res[,,s] contains prediction from state s, convert to list of dataframes
  res2 = vector("list", S)
  
  for (s in 1:S) {
    tmp = as.data.frame(res[,,s])
    if (min(dim(tmp))==1) tmp = res[,,s]
    names(tmp) = c("time",paste("pstate",1:S,sep=""))
    res2[[s]] = tmp
  }
  res2$trans = x1$trans
  res2$tmat = tmat
  return(res2)
}

tmat2 = function(tmat)
{
  dm = dim(tmat)
  mx = max(tmat,na.rm=TRUE)
  res = matrix(NA,mx,3)
  res[,1] = 1:mx
  transvec = as.vector(tmat)
  for (i in 1:mx) {
    idx = which(transvec==i)
    res[i,2:3] = c(idx %% dm[2],idx %/% dm[2] + 1)
  }
  res = data.frame(res)
  names(res) = c("transno","from","to")
  if (!is.null(dimnames(tmat))) {
    states = dimnames(tmat)[[1]]
    res$fromname = states[res$from]
    res$toname = states[res$to]
  }
  return(res)
}
