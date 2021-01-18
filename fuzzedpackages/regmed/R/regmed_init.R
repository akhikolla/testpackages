regmed_init <-
function(dat.obj,x.std=TRUE, med.std=TRUE){
	
  ## Note that y should always be centered, so always do this, even if 
  ## input y was already standardized
  
  y <- scale(dat.obj$y,center=TRUE,scale=FALSE)
  x <- scale(dat.obj$x,center=TRUE,scale=x.std)
  mediator <- scale(dat.obj$mediator,center=TRUE,scale=med.std) 

  n.med <- ncol(dat.obj$mediator)

  dat <- cbind(x, mediator, y)
  SampCov <- var(dat)
 
  ## regress each med on x to get residuals, using seemingly unrelated regression

  res<-apply(mediator,2,function(a,b) residuals(lm(a~0+b)),b=x)

  ## estimate penalized var matrix of residuals

  smat <- var(res)
  save.glasso <- glasso(smat, rho=.02, penalize.diagonal = FALSE)

  MedCov <- save.glasso$w
  dimnames(MedCov)<-list(colnames(mediator),colnames(mediator))
 
  sampleSize <- nrow(dat)

  Alpha <-  Beta <- rep(0, n.med)
  Delta <- 0.0
  varx <- var(x)
  vary <- var(y)

  vary.step.size <- 0.10 * vary
  if(vary.step.size < 0.01) vary.step.size <- .01
    
  return(list(Alpha=Alpha,Beta=Beta,Delta=Delta,varx=varx,vary=vary,SampCov=SampCov,
              MedCov=MedCov,vary.step.size=vary.step.size, sampleSize=sampleSize))

}
