regmed.fit <-
function(x, mediator, y, lambda, frac.lasso, 
         x.std =TRUE, med.std=TRUE, max.outer=5000, max.inner=100,
         step.multiplier = 0.5, wt.delta = .5,
         print.iter=FALSE){
  
  zed<-match.call()
 
  ## ----------input parameters
  ## x:  vector of exposure variable
  ## mediator: matrix of possible mediators
  ## y: vector of outcome
  ## labmda: a single penalty parmeter
  ## frac.lasso: fraction of penalty assigned to lasso penalty (in contrast to
  ##             group penalty). Simulations suggest frac.lasso = 0.8 offers
  ##             adequate control of false-positive results without much loss in power
  ## max.outer: maximum number of iterations in outer loop of optimization steps
  ## max.inner: maximum number of iterations in inner loop of optimization steps
  ## x.std: if TRUE, standardize x before analyses (center and scale by standard deviation)
  ## med.std: if TRUE, standardize mediator before analyses
  ## print.iter: if TRUE, print when each iteration of optimization is conducted (verbose output)

  ### check bounded parameters  
  if(lambda < 0 | lambda > 1) stop("invalid lambda, must be [0,1]")
  if(frac.lasso < 0 | frac.lasso > 1) stop("invalid frac.lasso, must be [0,1]")
  if(length(frac.lasso)!=1) stop("invalid frac.lasso, must be scalar")
  if(length(max.outer)!=1) stop("invalid max.outer, must be scalar")
  if(max.outer <= 0) stop("invalid max.outer, must be > 0")
  if(length(max.inner)!=1) stop("invalid max.inner, must be scalar")
  if(max.inner <= 0) stop("invalid max.inner, must be > 0")


  ### check x,y and mediator ###

  checked.dat<-regmed_dat_check(x=x,y=y,mediator=mediator)

  ### scale and center x,y and mediator as appropriate, initialize variables for rcpp_regmed ###

  inits<-regmed_init(dat.obj=checked.dat,x.std=x.std, med.std=med.std)
 
  ### fit regmed ###


  save<- rcpp_regmed(alpha=inits$Alpha, beta=inits$Beta, delta=inits$Delta,
                     vary = inits$vary, varx = inits$varx, SampCov = inits$SampCov,
                     inits$MedCov, inits$sampleSize, fracLasso = frac.lasso,
                     lambda = lambda, wt_delta = wt.delta,
                     max_iter=max.outer, max_iter_inner=max.inner,
                     tol=1e-6, vary_step_size = inits$vary.step.size,
                     step_multiplier = step.multiplier,
                     verbose=print.iter)

  ### G: this is done in regmed.grid, do you want to do this here as well?
  ### zero out alphas/betas if they are sufficiently close to 0, based on tol parameter ###    

  #savealpha <- ifelse(abs(save$alpha) < tol, 0.0, save$alpha)
  #save$beta  <- ifelse(abs(save$beta)  < tol, 0.0, save$beta)

  rownames(save$alpha)<-checked.dat$mediator.names
  rownames(save$beta)<-checked.dat$mediator.names

  ### check on invalid alpha/beta ###

  if(any(is.nan(save$alpha)) | any(is.nan(save$beta))) stop("invalid alpha or beta estimates")
  if(any(is.na(save$alpha)) | any(is.na(save$beta))) stop("invalid alpha or beta estimates")

  
  save$MedCov <- inits$MedCov

  save$call <- zed
  save$frac.lasso <- frac.lasso
  class(save)<-"regmed"

  return(save)

}
