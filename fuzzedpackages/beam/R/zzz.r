#######################
## Internal R functions
#######################

.upperTriIdxs <- function(p){
  z <- sequence(p)
  cbind(
    row = unlist(lapply(1:(p-1), function(x) 1:x), use.names = FALSE),
    col = rep(z[-1], times = tail(z, -1)-1))
}

.Idx2RowCol <- function(idx){
  col = ceiling((1 + sqrt(1+8*idx))/2)
  row = idx - (col^2 - 3*col +4)/2 + 1
  return(data.frame(row=row, col=col))
}

# Shrinkage estimation of variances
.shrinkvars <- function(s2jk, nk, method=c("eb", "mean", "median", "none")){
  
  if(method!="none"){
    
    # Functions needed
    logmlyk <- function(p, n, sj){
      a <- p[1]
      b <- p[2]
      astar <- a+0.5*n
      sum(-0.5*n*log(2*pi) + a*log(a-1) + a*log(b) - lgamma(a) + lgamma(astar) - astar*log((a-1)*b+0.5*n*sj))
    }
    logmlyka <- function(a, b, n, sj){
      astar <- a+0.5*n
      sum(-0.5*n*log(2*pi) + a*log(a-1) + a*log(b) - lgamma(a) + lgamma(astar) - astar*log((a-1)*b+0.5*n*sj))
    }
    atoalpha <- function(a, n){
      (a - 1)/(a + n - 1)
    }
    alphatoa <- function(alpha, n){
      ((n-1)*alpha + 1)/(1-alpha)
    }
    
    # Optimal shrinkage
    if(method=="eb"){
      resOpt <- optim(c(2,2), logmlyk, n=nk, sj=s2jk, control=list(fnscale=-1))
      a <- resOpt$par[1]
      b <- resOpt$par[2]
    }
    
    if(method=="mean" || method=="median"){
      
      if(method == "mean"){
        b <- mean(s2jk)
      }else{
        b <- median(s2jk)
      }
      
      lb <- alphatoa(1-0.999, nk)
      ub <- alphatoa(0.999, nk)
      resOpt <- optim(1, lower=lb , upper=ub, logmlyka, b = b, n = nk, sj = s2jk, method="Brent", control=list(fnscale=-1))
      a <- resOpt$par
      
    }
    
    # Shrunken variances
    alphaOpt <- atoalpha(a, nk)
    out <- alphaOpt*b + (1-alphaOpt)*s2jk
    
    print(resOpt)
    cat("convergence = ", resOpt$convergence, "\n")
    cat("alphaOpt = ", alphaOpt, "\n")
    cat("a = ", a, "\n")
    cat("b = ", b, "\n")
    
    
  }else{
    
    out <- s2jk
  }
  
  return(out)
}

# .logCPO <- function(mydelta, myp, myn, myX, myeigs, myXXT, myS, myD, mylogdetD){
#   part1 <- -0.5*myp*log(pi) + .lpvarGamma((mydelta+myn)*0.5, p=myp) - .lpvarGamma((mydelta+myn-1)*0.5, p=myp) 
#   part2 <- -0.5*sum(log((mydelta-myp-1)+myeigs))
#   mycpo <- part1 + part2
#   if(!is.null(mylogdetD)){
#     mycpo <- mycpo - 0.5*mylogdetD
#     TinvDelta <- chol2inv(chol((mydelta-myp-1)*myD + myS))
#   }else{
#     if(is.null(myXXT)){
#       TinvDelta <- chol2inv(chol((mydelta-myp-1)*diag(myp) + myS))
#     }else{
#       TinvDelta <- chol2inv(chol(diag(myn)+(1/(mydelta-myp-1))*myXXT))
#       TinvDelta <- crossprod(myX, TinvDelta)/((mydelta-myp-1)^2)
#       TinvDelta <- TinvDelta%*%myX
#       TinvDelta <- (1/(mydelta-myp-1))*diag(myp) - TinvDelta
#     }
#   }
#   mycpo <- mycpo + 0.5*(mydelta+myn-1)*log(1-diag((myX%*%TinvDelta)%*%t(myX)))
#   return(sum(mycpo))
# }

