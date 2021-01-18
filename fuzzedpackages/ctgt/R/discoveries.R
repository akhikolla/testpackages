## discoveries for standardized data

discoveries = function(y, X, xs, hyps, maxit=0, alpha=0.05){
  rh = actgt(y,X,xs,hyps,maxit=maxit,alpha=alpha)[1]
  if(rh == "not reject" ){
    dis=0
  } else if(rh=="unsure"){
    dis = "unsure"
  } else{
    sqrW = sqrt(mean(y)*(1-mean(y)) ) # sqrt of covariance of y
    
    WIHZ = sqrW *(sweep(X,2,colMeans(X))) ## W^{1/2}*(I-H)Z
    IHZ = WIHZ/sqrW  ##(I-H)Z
    
    itest = sapply(hyps, function(i) sum(colSums(y*IHZ[,i,drop=FALSE])^2))
    candidates = names(sort(itest))
    if(length(candidates)==1){
      dis=1
    } else{
      for(i in (length(candidates)-1):1){
        rc = actgt(y=y, X=X, xs=xs,hyps = candidates[1:i],maxit=maxit,alpha=alpha)
        if(rc[1]=="not reject" || rc[1]=="unsure"){
          dis = length(candidates)-i
          break
        }
        if(i==1 && rc[1]=="reject"){
          dis = length(candidates)-i
        }
      }
    }

  }
  return(dis)
}