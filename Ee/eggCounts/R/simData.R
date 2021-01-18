simData1s <- function(n=10, mean=500, kappa=0.5, phi=1, f=50, rounding=TRUE, seed = NULL){
  # error checking
  if (n%%1!=0) {stop("sample size is not an integer")}
  if (kappa<0) {stop("overdispersion parameter is not positive")}
  if ((phi<0)||(phi>1)) {stop("prevalence is not between 0 and 1")}
  if (sum(f%%1!=0)>0) {stop("correction factor is not integer(s)")}
  if (!is.logical(rounding)) {stop("rounding argument is not a logical")}
  if (length(f)>1 && length(f)!=n) {stop("the length of correction factors and sample size do not match")}
  # random seed
  if (!is.null(seed)) {set.seed(seed)}
   # sample true egg counts
   trueMean <- rgamma(n,shape=kappa, rate=kappa/mean)
   infected <- sample(c(FALSE,TRUE),size=n, replace=TRUE,prob=c(1-phi,phi))
   trueMean[!infected] <- 0
   trueEPG <- rpois(n,trueMean)
   # take a subsample and count eggs y, such that epg = f*y
   if (rounding){
   fec <- rpois(n,lambda=round(trueEPG/f))
   } else {
     fec <- rpois(n,lambda=trueEPG/f)}
   # fec <- rbinom(n,trueEPG,1/f)
   data <- data.frame(obs=fec*f, master=fec, true=trueEPG)
   return(data)
}

simData2s <- function(n=10, preMean=500, delta=0.1, kappa=0.5, deltaShape = NULL,
                      phiPre=1, phiPost=phiPre, f=50, paired=TRUE, rounding=TRUE, seed = NULL){
  # error checking
  if (n%%1!=0) {stop("sample size is not an integer")}
  if (kappa<0) {stop("overdispersion parameter is not positive")}
  if ((phiPre<0)||(phiPre>1)) {stop("pre-treatment prevalence is not between 0 and 1")}
  if ((phiPost<0)||(phiPost>1)) {stop("post-treatment prevalence is not between 0 and 1")}
  if (sum(f%%1!=0)>0) {stop("correction factor is not integer(s)")}
  if (!is.logical(paired)) {stop("paired is not a logical")}
  if (!is.logical(rounding)) {stop("rounding argument is not a logical")}
  if (length(f)>1 && length(f)!=n) {stop("the length of correction factors and sample size do not match")}
  # random seed
  if (!is.null(seed)) {set.seed(seed)}
  if(paired){
	# sample true egg counts before treatment
	truePreMean <- rgamma(n,shape=kappa, rate=kappa/preMean)
	infected <- sample(c(FALSE,TRUE),size=n, replace=TRUE,prob=c(1-phiPre,phiPre))
	truePreMean[!infected] <- 0
	truePreEPG <- rpois(n,truePreMean)
	# take a subsample and count eggs y, such that epg = f*y
	if (rounding){
	fec <- rpois(n,lambda=round(truePreEPG/f))
	} else {fec <- rpois(n,lambda=truePreEPG/f)}
	# fec <- rbinom(n,truePreEPG,1/f)	
	# now add sample egg counts after treatment
	infectedA <- infected
	infectedA[infected] <- sample(c(FALSE,TRUE),size=sum(infected), replace=TRUE,prob=c(1-phiPost/phiPre,phiPost/phiPre))

	if (is.null(deltaShape)){
	  deltas <- delta
	} else {
	  deltas <- rgamma(n, shape = deltaShape, rate = deltaShape/delta)
	}
	truePostMean <- truePreMean*deltas
	truePostMean[!infectedA] <- 0
	truePostEPG <- rpois(n,truePostMean)
	# take a subsample and count eggs y, such that epg = f*y
	if (rounding){
	postFEC <- rpois(n,lambda=round(truePostEPG/f))
	} else {postFEC <- rpois(n,lambda=truePostEPG/f)}
	#postFEC <- rbinom(n,truePostEPG,1/f)	
	data <- data.frame(obsPre =fec*f, masterPre=fec, truePre=truePreEPG, obsPost=postFEC*f, masterPost=postFEC, truePost=truePostEPG)

  } else {
    
    if (is.null(deltaShape)){
      deltas <- delta
    } else {
      deltas <- rgamma(n, shape = deltaShape, rate = deltaShape/delta)
    }
    
	dataPre <- simData1s(n=n, mean=preMean, kappa=kappa, phi=phiPre,rounding=rounding)
	dataPost <- simData1s(n=n, mean=preMean*deltas, kappa=kappa, phi=phiPost,rounding=rounding)

	data <- data.frame(dataPre,dataPost)
    colnames(data) <- paste(colnames(data),rep(c("pre","post"),each=3),sep=".")
	  
  }
   return(data)
}

