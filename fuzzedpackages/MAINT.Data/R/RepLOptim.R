RepLOptim <- function(start, parsd, fr, gr=NULL, inphess=NULL, ..., method="nlminb",
  lower=NULL, upper=NULL, rethess=FALSE, parmstder=FALSE, control=list()) 

{
  if (!is.null(control$EnfCnstrs)) EnfCnstrs <- control$EnfCnstrs 
  else EnfCnstrs <- FALSE
# EnfCnstrs ("Enforce Constraints") -- non-documented boolean argument. 
# Needs to be explicitly set to TRUE if non-linear constraints are to be
# strictly enforced (instead of using BigM penalties) by local optimizers.   

  maxrepet_default <- 50
  maxnoimprov_default <- 250
  maxreplic_default <- 1000
  maxiter_default <- 1500
  maxSANNiter_default <- 3000
  maxeval_default <- 2000
  srMachinetol <- sqrt(.Machine$double.eps)
#  RLOtol_default <- srMachinetol
  RLOtol_default <- 1e-3
  reltol_default <- srMachinetol
  
  if (!is.null(control$maxrepet)) maxrepet <- control$maxrepet else maxrepet <- maxrepet_default  
  if (!is.null(control$maxnoimprov)) maxnoimprov <- control$maxnoimprov else maxnoimprov <- maxnoimprov_default
  if (!is.null(control$maxreplic)) maxreplic <- control$maxreplic else maxreplic <- maxreplic_default
  if (!is.null(control$maxiter)) maxiter <- control$maxiter 
  else if (class(method)[1]!="character") maxiter <- maxiter_default
    else if (method!="SANN") maxiter <- maxiter_default
      else maxiter <-  maxSANNiter_default
  if (!is.null(control$maxeval)) maxeval <- control$maxeval else maxeval <- maxeval_default
  if (!is.null(control$objbnd)) objbnd <- control$objbnd else objbnd <- Inf
  if (!is.null(control$RLOtol)) RLOtol <- control$RLOtol else RLOtol <- RLOtol_default
  if (!is.null(control$rel.tol)) reltol <- control$reltol else reltol <- reltol_default
  maxiter1 <- control$maxiter1
  maxeval1 <- control$maxeval1
  allrep <- control$allrep
  if (!is.null(control$start)) start <- control$start 
  if (!is.null(control$sdfactor)) 
    if (!is.null(control$parsd)) parsd <- control$sdfactor * control$parsd 
    else parsd <- control$sdfactor * parsd 
  else if (!is.null(control$parsd)) parsd <- control$parsd 
  if (!is.null(control$inphess)) inphess <- control$inphess 
  if (!is.null(control$lower)) lower <- control$lower 
  if (!is.null(control$upper)) upper <- control$upper 
  if (!is.null(control$rethess)) lower <- control$rethess 
  if (!is.null(control$parmstder)) parmstder <- control$parmstder 
  if (!is.null(control$EnfCnstrs)) EnfCnstrs <- control$EnfCnstrs

  if (is.na(as.integer(maxrepet)) || length(maxrepet) > 1 || maxrepet<1) 
    stop("Wrong value (=",maxrepet,") for the control$maxrepet parameter\n")   

  npar <- length(start)
  values <- NULL
  if (is.null(lower)) lower <- rep(-Inf,npar)
  if (is.null(upper)) upper <- rep(Inf,npar)
  if (is.finite(objbnd))  { if (is.null(allrep)) allrep <- 10*maxreplic }
  else allrep <- maxreplic

  if (length(lower)!=npar) { stop("Incorrect length of the lower limits vector\n") }
  if (length(upper)!=npar) { stop("Incorrect length of the upper limits vector\n") }

  if (length(lower)!=npar) { 
    cat("npar =",npar," -- start =",start,"\nlower =",lower,"\n")
    stop("Incorrect length of the lower limits vector\n") 
  }
  if (length(upper)!=npar) { 
    cat("npar =",npar," -- start =",start,"\nupper =",upper,"\n")
    stop("Incorrect length of the upper limits vector\n") 
  }

  bestres <- NULL	
  bestval <- Inf
  initpar <- bestpar <- start
  cnt <- repcnt <- noimpcnt <- 0 
  if (!is.null(maxiter1) || !is.null(maxeval1)) nsteps <- 2
  else nsteps <- 1

  for (steps in 1:nsteps) for (i in 1:maxreplic)
  {
    if (steps==2) { maxiter <- maxiter1 ; maxeval <- maxeval1 ; cnt <- allrep-1 ; repcnt <- noimpcnt <- 0 }
    if (cnt > allrep || repcnt >= maxrepet || noimpcnt >= maxnoimprov) break
    value <- Inf
    while (value >= objbnd && cnt < allrep && repcnt < maxrepet && noimpcnt < maxnoimprov)
    {
      if (is.function(method)) {
#        tmpres <- method(initpar,gr=gr,lbound=lower, ubound=upper, control=control, hessian=rethess, ...)
        tmpres <- method(initpar,gr=gr,lbound=lower, ubound=upper, hessian=rethess, ...)
      } else {
        if (method == "nlminb") {
          tmpres <- nlminb(start=initpar,fr,gradient=gr,hessian=inphess,lower=lower,upper=upper,
            control=list(iter.max=maxiter,eval.max=maxeval,rel.tol=reltol),...)
        }
        else if (method == "nlm") 
          tmpres <- nlm(fr,p=initpar,lbound=lower,ubound=upper,iterlim=maxiter,...)
        else if (method == "L-BFGS-B") {
          tmpres <- optim(initpar,fr,gr=gr,method=method,lower=lower,upper=upper,
            control=list(maxit=maxiter),hessian=rethess,...)
        } 
        else if (method == "Nelder-Mead" || method == "BFGS" || method == "CG" || method == "SANN")
          tmpres <- optim(initpar,fr,gr=gr,method=method,control=list(maxit=maxiter),
            lbound=lower,ubound=upper,hessian=rethess,...)
        else stop(paste("The argument 'method' is neither a function object\n",
          "or a string describing any of the available local optimizers\n")) 
      }
      if (EnfCnstrs)
      {
        if (is.function(method) || method != "nlm") par <- tmpres$par
        else par <- tmpres$estimate
        value <- fr(par,nopenalty=TRUE,...)
      }
      else  {  
        if (!is.function(method))  {
          if (method == "nlminb") value <- tmpres$objective
          else if (method == "nlm") value <- tmpres$minimum
               else value <- tmpres$value
        }
        else value <- tmpres$value
      }
      if (is.null(value) || is.na(value) || value==.Machine$double.xmax) value <- objbnd
      cnt <- cnt+1
      values <- c(values,value)
      if (is.na(value)) { noimpcnt <- noimpcnt + 1 ; repcnt <- 0 }
      else
      {
        if (is.finite(bestval))
          if (abs((value-bestval)/bestval) < RLOtol) repcnt <- repcnt + 1 
          else repcnt <- 0 
        if (value < bestval)
        {
          bestval <- value
          if (is.function(method) || method != "nlm") bestpar <- tmpres$par
          else bestpar <- tmpres$estimate
          bestres <- tmpres
          noimpcnt <- 0 
        }
        else noimpcnt <- noimpcnt + 1
      }
      if (cnt < allrep)
      { 
        hlfrng <- sqrt(3)*parsd  #  generate new parameters from an uniform distribution
        initpar <- runif(npar,min=pmax(lower,bestpar-hlfrng),max=pmin(upper,bestpar+hlfrng))
      #  u <- runif(n=npar)     # generate npar uniform random numbers
      #  initpar <- qnorm(u,mean=bestpar,sd=parsd) #  generate new parameters from a normal distribution
      #  lbndind <- which(initpar<lower)   #  identify indices of parameters that fell below their lower bounds
      #  ubndind <- which(initpar>upper)   #  identify indices of parameters that fell above their upper bounds
      #  initpar[lbndind] <- lower[lbndind] + u[lbndind] * (bestpar[lbndind]-lower[lbndind]) # and correct them
      #  initpar[ubndind] <- upper[ubndind] - u[ubndind] * (upper[ubndind]-bestpar[ubndind])
      }
    } 
  }
 
  if (!is.function(method) && method == "nlminb")
  {
    iterations <- bestres$iterations
    counts <- bestres$evaluations
    hess <- NULL
    egval <- NULL
    parstd <- NULL
  } 
  else
  {
    iterations <- NULL
    counts <- bestres$counts
    if (rethess==TRUE)
    {
      hess <- bestres$hessian
      egval <- eigen(hess,symmetric=TRUE,only.values=TRUE)$values
      if (parmstder==TRUE)
        if (egval[npar] < RLOtol) parstd <- "Not computed because the hessian is not positive definite"
        else parstd <- sqrt(diag(solve(hess)))
      else parstd <- NULL
    }
    else {
      hess <- NULL
      egval <- NULL
      parstd <- NULL
    }
  } 

  if (!is.null(bestres))
    return( list(par=bestpar,val=bestval,iterations=iterations,vallist=values,counts=counts,
      convergence=bestres$convergence,message=bestres$message,hessian=hess,hessegval=egval,stderrors=parstd) )
  else
    return( list(par=NULL,val=Inf,iterations=NULL,vallist=NULL,counts=NULL,convergence=NULL,
      message="RepLOptim was unable to find any valid solution",hessian=NULL,hessegval=NULL,stderrors=NULL) )

}
