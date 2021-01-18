`p_to_e` <- function(p){c(sum(p), p[-length(p)]/rev(cumsum(rev(p))[-1]))}

`e_to_p` <- function(e){e[1]* cumprod(c(1,1-e[-1]))* c(e[-1],1)}

`B` <- function(H, disallowed=NULL, give=FALSE, ...){  # modelled on hyperdirichlet::calculate_B()
  n <- size(H)-1  # "-1" because this is the dimension of the integrand, not the number of p's.

  ## No special dispensation for the Dirichlet!

  if(is.null(disallowed)){
    f <- function(e){ dhyper2_e(e,H, include.Jacobian=TRUE) }  # integrate over the whole simplex
  } else {
    f <- function(e){as.numeric(!disallowed(e_to_p(c(1,e)))) * dhyper2_e(e, H, include.Jacobian=TRUE)}
  }

  out <- adaptIntegrate(f,lowerLimit=rep(0,n),upperLimit=rep(1,n), ...)
  if(give){
    return(out)
  } else {
    return(out$integral)
  }
}

`Jacobian` <-
function(e){
  n <- length(e)
  e1 <- e[1]
  e <- e[-c(1,length(e))]
  prod(cumprod(1-e))*e1^n
}

`dhyper2` <- function(P,H,...){
  P <- rbind(P)
  loglik(H,P,log=FALSE)/B(H,...)
}

`dhyper2_e` <-  # analogous to dhyperdirichlet_e()
function(e, H, include.Jacobian=TRUE){
  
  e <- c(1,e)
  p <- e_to_p(e)

  out <- loglik(H,indep(p),log=FALSE)
  if(include.Jacobian){out <- Jacobian(e)*out}
  return(out)
}

`mgf` <- function(H,powers,...){ B(H + dirichlet(powers = powers), ...)/B(H, ...) }

`mean_hyper2` <- function(H, normalize=TRUE, ...){
  f <- function(i){
    jj <- rep(0,size(H))
    jj[i] <- 1
    return(mgf(H,jj, ...))
  }

  out <- sapply(seq_len(size(H)),f)
  if(!identical(pnames(H),NA)){names(out) <- pnames(H)}
  if(normalize){
    out <- out/sum(out)
  }
  return(out)
}

"probability" <- function(H , disallowed=NULL, ...){
  B(H , disallowed=disallowed , give=FALSE, ...) / B(H, ...)
}
