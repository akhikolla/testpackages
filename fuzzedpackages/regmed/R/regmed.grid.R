regmed.grid <-
function(x, mediator, y, lambda.vec, frac.lasso = 0.8, max.outer=5000, max.inner=100,
                        x.std=TRUE, med.std=TRUE,  step.multiplier = 0.5, wt.delta = .5, print.iter = FALSE){
 
  zed<-match.call()

  
  ## ----------input parameters
  ## x:  vector of exposure variable
  ## mediator: matrix of possible mediators
  ## y: vector of outcome
  ## lambda.vec: vector of penalty parmeters (preferable in decreasing order for
  ##             using parameters estimates from one penalty output as warm start
  ##             initial values for next penalty parameter)
  ## frac.lasso: fraction of penalty assigned to lasso penalty (in contrast to
  ##             group penalty). Simulations suggest frac.lasso = 0.8 offers
  ##             adequate control of false-positive results without much loss in power
  ## max.outer: maximum number of iterations in outer loop of optimization steps
  ## max.inner: maximum number of iterations in inner loop of optimization steps
  ## x.std: if TRUE, standardize x before analyses (center and scale by standard deviation)
  ## med.std: if TRUE, standardize mediator before analyses
  ## print.iter: if TRUE, print when each iteration of optimization is conducted (verbose output)

  tol <- 1e-4 


  ### check bounded parameters  
  if(any(lambda.vec < 0)) stop("invalid lambda, lambda must be >= 0")
  if((frac.lasso < 0) | (frac.lasso > 1)) stop("invalid frac.lasso, must be [0,1]")
  if(length(frac.lasso)!=1) stop("invalid frac.lasso, must be scalar")
  if(length(max.outer)!=1) stop("invalid max.outer, must be scalar")
  if(max.outer <= 0) stop("invalid max.outer, must be > 0")
  if(length(max.inner)!=1) stop("invalid max.inner, must be scalar")
  if(max.inner <= 0) stop("invalid max.inner, must be > 0")



  ### check x,y and mediator ###

  checked.dat<-regmed_dat_check(x=x,y=y,mediator=mediator)

  ### scale and center x,y and mediator as appropriate, initialize variables for rcpp_regmed ###

  inits<-regmed_init(dat.obj=checked.dat,x.std=x.std, med.std=med.std)
 
   
  ### loop through lambda values ###

  fit.lst <- vector("list", length=length(lambda.vec))

  for(i in 1:length(lambda.vec)){
  
    if(print.iter) cat("lambda i = ", i, "\n")
  
    save <- rcpp_regmed(alpha=inits$Alpha, beta=inits$Beta, delta=inits$Delta,
                        vary = inits$vary, varx = inits$varx, SampCov = inits$SampCov,
                        inits$MedCov, inits$sampleSize, fracLasso = frac.lasso,
                        lambda = lambda.vec[i], wt_delta = wt.delta,
                        max_iter= max.outer, max_iter_inner=max.inner,
                        tol=1e-6,vary_step_size = inits$vary.step.size,
                        step_multiplier = step.multiplier,
                        verbose=print.iter)


    ### zero out alphas/betas if they are sufficiently close to 0, based on tol parameter ###    

    save$alpha <- ifelse(abs(save$alpha) < tol, 0.0, save$alpha)
    save$beta  <- ifelse(abs(save$beta)  < tol, 0.0, save$beta)

    rownames(save$alpha)<-checked.dat$mediator.names
    rownames(save$beta)<-checked.dat$mediator.names

    ### check on invalid alpha/beta ###

    if(any(is.nan(save$alpha)) | any(is.nan(save$beta))) stop("invalid alpha or beta estimates")
    if(any(is.na(save$alpha)) | any(is.na(save$beta))) stop("invalid alpha or beta estimates")
   
    fit.lst[[i]] <- save

    ### update initialization parameters: Alpha, Delta, Beta, varx and vary ###

    inits <- regmed_grid_update(regmed.fit.obj=save,regmed.inits=inits)

  }
   
  gridData <- regmed_gridData(fit.lst=fit.lst,lambda.vec=lambda.vec)
  
  out<-list(fit.list=fit.lst, grid.data=gridData, sample.size=inits$sampleSize, MedCov=inits$MedCov, call=zed)

    out$frac.lasso <- frac.lasso
    
  class(out)<-"regmed.grid"

  return(out)

}
