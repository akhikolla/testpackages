# functions 
lp <- function(p,i,r,theta,delta,gamma,k,scaling=FALSE){
  
  if(k%%2==0){ # k even, odd number of categories 
    m <- k/2
    if(!scaling){
      predictor <- theta[p]+sign(m-r+0.5)*gamma[p]-delta[i,r]
    } else{
      predictor <- theta[p]+(m-r+0.5)*gamma[p]-delta[i,r]
    }
  } else{ # k odd, even number of categories 
    m <- floor(k/2)+1
    if(!scaling){
      predictor <- theta[p]+sign(m-r)*gamma[p]-delta[i,r]
    } else{
      predictor <- theta[p]+(m-r)*gamma[p]-delta[i,r]
    }
  }
  return(predictor)
}

prob_pir <- function(p,i,r,theta,delta,gamma,k,...){
  
  lps <- sapply(1:k,function(l) lp(p,i,l,theta,delta,gamma,k,...))
  
  counter <- exp(sum(lps[1:r]))
  denomi  <- 1+sum(exp(cumsum(lps)))
  prob    <- counter/denomi
  
  return(prob)
}

prob_pi <- function(p,i,theta,delta,gamma,k,...){
  
  probs <- numeric(k+1)
  probs[c(2:(k+1))] <- sapply(1:k,function(l) prob_pir(p,i,l,theta,delta,gamma,k,...))
  probs[1] <- 1-sum(probs[c(2:(k+1))])
  
  return(probs)
}

y_pi <- function(p,i,theta,delta,gamma,k,...){
  
  prob <- prob_pi(p,i,theta,delta,gamma,k,...)
  y    <- sample(0:k,1,prob=prob)
  
  return(y)
}

y_i <- function(P,i,theta,delta,gamma,k,...){
  
  y_i <- sapply(1:P,function(l) y_pi(l,i,theta,delta,gamma,k,...) )
  
  return(y_i)
}

sim_Y <- function(P,I,theta,delta,gamma,k,...){
  
  Y <- sapply(1:I,function(l) y_i(P,l,theta,delta,gamma,k,...))
  
  return(Y)
}
