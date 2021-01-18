globalVariables(c(".", ".SD"))
  
aggr <- function(x, clusters){
  #seems like lapply(.SD, sum) on data.table drops 
  #columns with identical names, thus ensure that no identical names
  colnames(x) <- NULL 
  temp <- data.table(x)
  temp <- as.matrix(temp[, j=lapply(.SD, sum), by=clusters][, -1])

} 

cglm <- function(method, formula, data, id, link, ...){
  
  dots <- list(...)
   
  #preparation
  data <- data[order(data[, id]), ]
  id <- data[, id]
  X <- model.matrix(object=formula, data=data)[, -1, drop=FALSE]
  y <- data[, as.character(formula)[2]] 
  ni <- as.vector(tapply(X=id, INDEX=id, FUN=length))
  n <- length(unique(id))
  nbeta <- ncol(X)
  
  if(method=="cml"){  
  
    nidcumsum <- c(0, cumsum(ni))  
    args <- list(par=rep(0, ncol(X)), fn=llfun, y=y, design=X, 
      nidcumsum=nidcumsum, hessian=TRUE)
    args[names(dots)] <- dots
    if(is.na(match("control", names(args))))
      args$control <- list(maxit=1000)
    fit <- do.call("optim", args=args)
    convergence <- fit$convergence==0
    coefficients <- fit$par
    var <- try(solve(fit$hessian))
    if(inherits(x=var, what="try-error")){
      var <- matrix(nrow=nbeta, ncol=nbeta)
      convergence <- FALSE
    }    
    coeff.names <- attr(terms(formula), "term.labels")
    names(coefficients) <- coeff.names
    colnames(var) <- coeff.names
    rownames(var) <- coeff.names  
    
  }
  
  if(method=="ts"){
  
    X.mean <- aggr(x=X, clusters=id)/ni   
    X.mean <- X.mean[rep(1:n, ni), ]
    X.cent <- X-X.mean
    
    resfun <- function(betastar){
      pred <- as.vector(X%*%matrix(betastar))
      if(link=="identity")
        res <- y-pred
      if(link=="log")
        res <- y*exp(-pred)
      return(res)
    }
    Ubetastarfun <- function(betastar){
      res <- resfun(betastar)  
      U <- aggr(x=X.cent*res, clusters=id)
      return(U)
    } 
    Ubetastarfunsums <- function(betastar){
      return(colSums(Ubetastarfun(betastar)))
    }
    
    #---ESTIMATES---
     
    #estimate of beta
  
    if(link=="identity"){
      betastar <- as.vector(solve(t(X.cent)%*%X.cent)%*%(t(X.cent)%*%matrix(y))) 
      convergence <- TRUE
    }
    if(link=="log"){
      args <- list(x=rep(0, nbeta), fn=Ubetastarfunsums)
      args[names(dots)] <- dots
      if(is.na(match("control", names(args))))
        args$control <- list(maxit=100)
      fit <- do.call("nleqslv", args=args)
      betastar <- fit$x
      convergence <- fit$termcd==1
    }
    
      
    #estimate of phi
    
    res <- resfun(betastar)
    pred <- as.vector(X%*%matrix(betastar))
    if(link=="identity"){
      b <- rep(as.vector(aggr(x=res, clusters=id))/ni, ni)
      m <- b+pred
      App <- 1
      Bi <- ni>1
    }   
    if(link=="log"){
      b <- log(rep(as.vector(aggr(x=res, clusters=id))/ni, ni))
      m <- exp(b+pred)
      App <- m
      all.zeros <- tapply(abs(y), id, FUN=sum)==0
      Bi <- ni>1 & !all.zeros
    }
    eps <- y-m 
    phiij <- eps^2/App
    phii <- tapply(phiij, id, FUN=mean)    
    phii <- phii*ni/(ni-1)
    phii[!Bi] <- 0
    phi <- mean(phii[Bi])
    
    #estimate of beta
    
    beta <- betastar/phi
    names(beta) <- colnames(X)
    coefficients <- beta
    
    #---VARIANCE---

    #meat
    Ubetastar <- Ubetastarfun(betastar)
    Uphi <- Bi*(phii-phi)
    Ubeta <- matrix(rep(betastar/phi-beta, each=n), nrow=n, ncol=nbeta)
    U <- cbind(Ubetastar, Uphi, Ubeta) 
    J <- var(U)
    
    #bread
    if(link=="identity")
      dUbetastar.dbetastar <- -t(X.cent)%*%X/n
    if(link=="log")
      dUbetastar.dbetastar <- -t(X.cent)%*%(X*res)/n  
    dUbetastar.dphi <- 0
    dUbetastar.dbeta <- matrix(0, nrow=nbeta, ncol=nbeta)
    dUbetastar <- cbind(dUbetastar.dbetastar, dUbetastar.dphi, dUbetastar.dbeta)
    if(link=="identity"){
      dm.dbetastar <- X.cent
      dUphi.dbetastar <- -2*eps*dm.dbetastar
    }
    if(link=="log"){
      tmp1 <- X*res
      tmp1 <- aggr(x=tmp1, clusters=id)/ni 
      tmp1 <- tmp1[rep(1:n, ni), ]
      res <- rep(as.vector(aggr(x=res, clusters=id))/ni, ni)
      res <- res[rep(1:n, ni)]
      tmp2 <- X*res
      dm.dbetastar <- exp(pred)*(tmp2-tmp1) 
      dUphi.dbetastar <- -dm.dbetastar*eps/m*(1+eps/m)
    } 
    dUphi.dbetastar <- aggr(x=dUphi.dbetastar, clusters=id)/(ni-1)
    dUphi.dbetastar[!Bi, ] <- 0 
    dUphi.dbetastar <- colMeans(dUphi.dbetastar)
    dUphi.dphi <- -mean(Bi)
    dUphi.dbeta <- rep(0, nbeta)
    dUphi <- c(dUphi.dbetastar, dUphi.dphi, dUphi.dbeta)
    dUbeta.dbetastar <- diag(1/phi, nbeta)
    dUbeta.dphi <- -betastar/phi^2
    dUbeta.dbeta <- -diag(nbeta)
    dUbeta <- cbind(dUbeta.dbetastar, dUbeta.dphi, dUbeta.dbeta)
    I <- rbind(dUbetastar, dUphi, dUbeta)
 
    #sandwich formula
    var <- solve(I)%*%J%*%t(solve(I))/n
    var <- var[(nbeta+2):(2*nbeta+1), (nbeta+2):(2*nbeta+1), drop=FALSE] 
    rownames(var) <- colnames(X)
    colnames(var) <- colnames(X)
      
  }  
  
  out <- list(call=match.call(), coefficients=coefficients, var=var,
    convergence=convergence) 
  class(out) <- "cglm"
  return(out) 
  
}

summary.cglm <- function(object, ...) {
  s.err <- sqrt(diag(as.matrix(object$var)))
  zvalue <- object$coefficients/s.err
  pvalue <- 2*pnorm(-abs(zvalue))
  coef.table <- as.matrix(cbind(object$coefficients, s.err, zvalue, 
        pvalue))
  dimnames(coef.table) <- list(names(object$coefficients), c("Estimate", 
      "Std. Error", "z value", "Pr(>|z|)"))
  ans <- list(call=object$call, coefficients=coef.table, 
    convergence=object$convergence)
  class(ans) <- "summary.cglm"
  return(ans)
}

print.summary.cglm <- function(x, digits=max(3L, getOption("digits")-3L), 
  signif.stars=getOption("show.signif.stars"), ...) {
  cat("\nCall:  ", "\n", paste(deparse(x$call), sep="\n", collapse="\n"), 
      "\n", sep="")
  if(!x$convergence)
    cat("\nWarning: no solution to the estimating equations was found,
      consider trying different values of control parameters", "\n")  
  cat("\nCoefficients:", "\n")
  printCoefmat(x$coefficients, digits=digits, signif.stars=signif.stars, 
    na.print = "NA", ...)
  cat("\n") 
}
