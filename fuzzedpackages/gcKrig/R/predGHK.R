
#######################################################################################################################
#######################################################################################################################
#### Inner function only used to be called by another functions inside package
#######################################################################################################################
#######################################################################################################################

#### param: Last one (or two) elements are: Range and Nugget (if applicable)
#### The first (1 + number of covariates) parameters are: regression parameters
#### Before Range parameter: overdispersion parameter (if applicable)
#### y: counts.
#### x: covariates
#### kappa: smoothness parameter in matern marginal of correlations





ghkgcmr_R <- function(R, pdf, cdf, nrep){
  .Call('ghkgcmr', PACKAGE = 'gcKrig', R, pdf, cdf, nrep)
}

ghkgcmr2_R <- function(R, pdf, cdf, nrep){
  .Call('ghkgcmr2', PACKAGE = 'gcKrig', R, pdf, cdf, nrep)
}




likGHK <- function(
  pars, y, x = NULL, locs, marginal, corr, effort, longlat = FALSE, distscale = 1,
  nrep = 1000, reorder = FALSE, seed = 12345
  ){

  if(longlat == FALSE){
   matD <- as.matrix(dist(locs, method = "euclidean", diag = TRUE, upper = TRUE))*distscale
 }else if(longlat == TRUE){
   if (requireNamespace("sp", quietly = TRUE)) {
   matD <- sp::spDists(x = as.matrix(locs), longlat = TRUE)*distscale
   }else{
     stop("Please install {sp} first!")
   }
 }
  npars <- length(pars)
  nparcorr <- corr$npar.cor

  if(marginal$nod == 1 & pars[ncol(x)+1] < 0){
    pars[ncol(x)+1] <- 0
  }
  if(corr$nug == 0 & pars[npars] < 0){
    pars[npars] <- 0
  }else if(corr$nug == 1){
    pars[npars] <- ifelse(pars[npars] > 1, 1, pars[npars])
    pars[npars] <- ifelse(pars[npars] < 0, 0, pars[npars])
    pars[npars-1] <- ifelse(pars[npars-1] < 0, 0, pars[npars-1])
  }

  R <- corr$corr(pars[(npars-nparcorr+1):npars], D = matD)
  pdf = marginal$pdf(y = y, x = x, pars = pars, effort = effort)
  cdf = marginal$cdf(y = y, x = x, pars = pars, effort = effort)
  set.seed(seed)
  if(reorder == FALSE){
    loglik <- ghkgcmr_R(R = R, pdf = pdf, cdf = cdf, nrep = nrep)
  } else{
    loglik <- ghkgcmr2_R(R = R, pdf = pdf, cdf = cdf, nrep = nrep)
  }
  if(is.nan(loglik) | loglik == -Inf) loglik = -1e6
  if (loglik == Inf) loglik = 1e6
  return(-loglik)
}


profilelikGHK <- function(
    theta, par.index, fitted, fixvalue = NULL, single = TRUE, alpha, nrep = 1000, seed = 12345
){
  if(single == TRUE){
    optpar <- append(fixvalue, theta, after = par.index-1)
    llik <- likGHK(pars = optpar, y = fitted$args$y, x = fitted$args$x, locs = fitted$args$locs,
                   marginal = fitted$args$marginal, corr = fitted$args$corr, effort = fitted$args$effort,
                   longlat = fitted$args$longlat, distscale = fitted$args$distscale,
                   nrep = nrep, reorder = FALSE, seed = seed)
    llikout <- (-llik - (fitted$log.lik-qchisq(1-alpha,1)/2) )^2
  }else{
    optpar <- append(theta, fixvalue, after = par.index-1)
    llikout <- likGHK(pars = optpar, y = fitted$args$y, x = fitted$args$x, locs = fitted$args$locs,
                      marginal = fitted$args$marginal, corr = fitted$args$corr, effort = fitted$args$effort,
                      longlat = fitted$args$longlat, distscale = fitted$args$distscale,
                      nrep = nrep, reorder = FALSE, seed = seed)
  }
  return(llikout)
}


mleProfilelikGHK <- function(
    fitted, par.index, fixvalue = NULL, single = TRUE, start, alpha, nrep = 1000, seed = 12345
){
  if(single == TRUE){
    est <- optim(start, profilelikGHK, par.index = par.index, fitted = fitted,fixvalue = fixvalue,
                 single = TRUE, alpha = alpha, nrep = nrep, seed = seed, method = 'L-BFGS-B',
                 lower = fitted$optlb[par.index], upper = fitted$optub[par.index])$par

  }else{
    est <- optim(start, profilelikGHK, par.index = par.index, fitted = fitted, fixvalue = fixvalue,
                 single = FALSE, alpha = alpha, nrep = nrep, seed = seed, method = 'L-BFGS-B',
                 lower = fitted$optlb[-par.index], upper = fitted$optub[-par.index])$par
  }
  return(est)
}


#######################################################################################################################
#######################################################################################################################
#### Inner function only used to be called by another functions inside package
#######################################################################################################################
#######################################################################################################################


mleGHK <- function(
  y, x = NULL, locs, marginal, corr, effort = 1, longlat = FALSE, distscale = 1, corrpar0 = NULL,
  ghkoptions = list(nrep = c(100,1000), reorder = FALSE, seed = 12345)
  ){
  if(longlat == FALSE){
    D <- as.matrix(dist(locs, method = "euclidean", diag = TRUE, upper = TRUE))*distscale
  }else if(longlat == TRUE){
    if (requireNamespace("sp", quietly = TRUE)) {
      D <- sp::spDists(x = as.matrix(locs), longlat = TRUE)*distscale
     }else{
      stop("Please install {sp} first!")
    }
  }
  N <- length(y)
  if(length(effort) == 1 & effort[1] == 1) effort <- rep(1,N)
  if(!is.null(x)){
  x <- as.data.frame(x)
  x <- cbind(rep(1,N), x)
  }else{
    x <- cbind(rep(1,N), x)
  }
  marg0 <- marginal$start(y = y, x = x, effort = effort)
  n.nugget0 <- corr$nug
  n.reg0 <- ncol(x)
  n.od0 <- marginal$nod
  if(is.null(corrpar0)){
    corpar0 <- corr$start(D)
  }else{
    corpar0 <- corrpar0
    names(corpar0) <- names(corr$start(D))
  }
   est <- c(marg0, corpar0)
   optlb <- c(rep(-Inf,n.reg0), rep(0, n.od0), 0, rep(0, n.nugget0))
   optub <- c(rep(Inf,n.reg0), rep(Inf, n.od0), Inf, rep(1, n.nugget0))

  for(j in 1:length(ghkoptions$nrep)){
    fit <- optim(par = est, fn = likGHK, y = y, x = x, locs = locs, marginal = marginal, corr = corr,
                 effort = effort, longlat = longlat, distscale = distscale, nrep = ghkoptions$nrep[j],
                 reorder = ghkoptions$reorder, seed = ghkoptions$seed, method = "L-BFGS-B",
                 lower = optlb, upper = optub)

    est <- fit$par
  }
   kmarg <- length(marg0)
   k <- length(est)

   if (is.null(fit$convergence) || fit$convergence != 0)
    warnings("Maximum likelihood estimation failed. Algorithm does not converge or MLEs do not exist.")

    result.list <- list(MLE = est,
                        x = x,
                        nug = n.nugget0,
                        nreg = n.reg0,
                        log.lik = -fit$value,
                        AIC = 2*k+2*fit$value,
                        AICc = 2*k+2*fit$value + 2*k*(k+1)/(N-k-1),
                        BIC = 2*fit$value + k*log(N),
                        kmarg = kmarg,
                        par.df = k,
                        N = N,
                        D = D,
                        optlb = optlb,
                        optub = optub,
                        args = mget(names(formals()),sys.frame(sys.nframe())))
  return(result.list)
}


likGHKXX <- function(pars, y, XX, locs, marginal, corr, effort, longlat, distscale, nrep, seed)
{
  likGHK(pars = pars, y = y, x = cbind(rep(1, length(y)), XX), locs = locs, marginal = marginal, corr = corr,
         effort = effort, longlat = longlat, distscale = distscale, nrep = nrep, reorder = FALSE, seed = seed)
}


#######################################################################################################################
#######################################################################################################################
#### OUTPUT FUNCTION: Prediction using GHK simulator: serial version
#######################################################################################################################
#######################################################################################################################

#### Output: n(number of prediction locations)*2 matrix
#### First Column: predicting value
#### Second Column: Estimated MSPE

predGHK <- function(
   obs.y, obs.x = NULL, obs.locs, pred.x = NULL, pred.locs, longlat = FALSE, distscale = 1,
    marginal, corr, obs.effort = 1, pred.effort = 1, estpar = NULL,
    corrpar0 = NULL, pred.interval = NULL,
     ghkoptions = list(nrep = c(100, 1000), reorder = FALSE, seed = 12345)
){

  x <- as.data.frame(cbind(rep(1,length(obs.y)), obs.x))
  colnames(x)[1] <- c("Intercept")

  if(!is.matrix(pred.locs) & !is.data.frame(pred.locs))
    stop("Input 'pred.locs' must be a data frame or matrix!")

  if(!is.matrix(obs.locs) & !is.data.frame(obs.locs))
    stop("Input 'obs.locs' must be a data frame or matrix!")

  if(length(obs.effort) == 1) obs.effort <- rep(obs.effort, nrow(obs.locs))
  if(!length(obs.effort) == nrow(obs.locs))
    stop("Sampling Effort must be equal to the number of sampling locations!")

  if(length(pred.effort) == 1) pred.effort <- rep(pred.effort, nrow(pred.locs))

  if(!length(pred.effort) == nrow(pred.locs))
    stop("Prediction Effort must be equal to the number of prediction locations!")

  #### First calculate the MLE and output the log-likelihood as denominator

  if(is.null(estpar)){

  MLE.est <- mleGHK(y = obs.y, x = x[,-1], locs = obs.locs, marginal = marginal, corr = corr, effort = obs.effort,
                    longlat = longlat, distscale = distscale, corrpar0 = corrpar0, ghkoptions = ghkoptions)

  loglik <- MLE.est$log.lik; estpar <- MLE.est$MLE

  }else{
  loglik <- -likGHK(pars = estpar, y = obs.y, x = x, locs = obs.locs, marginal = marginal, corr = corr,
           effort = obs.effort, longlat = longlat, distscale = distscale,
           nrep = ghkoptions$nrep[length(ghkoptions$nrep)],
           reorder = ghkoptions$reorder, seed = ghkoptions$seed)
  }

  if(is.null(pred.x)) {
    pred.x <- matrix(1, nrow = nrow(pred.locs), ncol = 1)
  }else {
    pred.x <- cbind(rep(1, nrow(pred.locs)) , pred.x)
  }
  pred.x <- as.data.frame(pred.x)
  names(pred.x) <- names(x)

  if(nrow(pred.x)!= nrow(pred.locs))
    stop("Number of prediction locations did not match rows of covariates")

  if (requireNamespace("FNN", quietly = TRUE)) {
   if(nrow(pred.locs) == 1){
    indexloc <- which.min(FNN::get.knnx(pred.locs, obs.locs, 1)$nn.dist)
    m0 <- n0 <- round(unique(pred.effort)*obs.y[indexloc]/obs.effort[indexloc])+1
  }else{
    m0 <- n0 <- round(pred.effort*apply(pred.locs, 1, function(x) obs.y[which.min(FNN::get.knnx(t(as.matrix(x)),
             obs.locs, 1)$nn.dist)]/obs.effort[which.min(FNN::get.knnx(t(as.matrix(x)), obs.locs, 1)$nn.dist)]))+1
   }
  }else{
    stop("Please install {FNN} first!")
  }
  NPL <- length(m0)
  #### m0 and n0: initial prediction values for probability search. Scalor or Vector
  max.count <- ceiling(5*max(obs.y/obs.effort)*max(pred.effort))

  #### A for loop which cannot be avoided

  if(is.null(pred.interval)){
    ans <- matrix(NA, nrow = NPL, ncol = 2)
  } else if(pred.interval >=0 & pred.interval<=1 ) {
    ans <- matrix(NA, nrow = NPL, ncol = 6)
  }else stop("Input pred.interval must be a number between 0 and 1!")
  nnrep <- length(ghkoptions$nrep)

  for(j in 1:NPL){
    tmpfun <- function(xtmp) {
      exp(-likGHK(pars = estpar, y = c(obs.y, xtmp), x = rbind(as.matrix(x), pred.x[j,]),
                  locs = rbind(obs.locs, pred.locs[j,]), marginal = marginal, corr = corr,
                  effort = c(obs.effort, pred.effort[j]), longlat = longlat, distscale = distscale,
                  nrep = ghkoptions$nrep[nnrep], reorder = ghkoptions$reorder, seed = ghkoptions$seed) - loglik)
    }

    p.m0 <- p.n0 <- tmpfun(m0[j]); mu.m0 <- mu.n0 <- p.m0*m0[j];  mu2.m0 <- mu2.n0 <- p.m0*m0[j]^2

    MM1 <- matrix(0, nrow = 2, ncol = 4); MM2 <- matrix(0, nrow = 2, ncol = 4)

    MM1[1,] <- c(p.m0, mu.m0, mu2.m0, m0[j]);  MM2[1,] <- c(p.n0, mu.n0, mu2.n0, n0[j])

    p.m0 <- tmpfun(m0[j]-1); mu.m0 <- p.m0*(m0[j]-1)
    mu2.m0 <- p.m0*(m0[j]-1)^2;  MM1[2,] <- c(p.m0, mu.m0, mu2.m0, m0[j]-1)

    while( (p.m0 > sqrt(.Machine$double.eps) | MM1[nrow(MM1), 2] > MM1[nrow(MM1)-1, 2]) & m0[j] > 1)
    {
      p.m0 <- tmpfun(m0[j]-2); mu.m0 <- p.m0*(m0[j]-2); mu2.m0 <- p.m0*(m0[j]-2)^2
      MM1 <- rbind(MM1, c(p.m0, mu.m0, mu2.m0, m0[j]-2)); m0[j] <- m0[j]-1
    }
    #### Search from n0 to the right
    p.n0 <- tmpfun(n0[j]+1); mu.n0 <- p.n0*(n0[j]+1); mu2.n0 <- p.n0*(n0[j]+1)^2
    MM2[2, ] <- c(p.n0, mu.n0, mu2.n0, n0[j]+1)

    while( (p.n0 > sqrt(.Machine$double.eps) | MM2[nrow(MM2), 2] > MM2[nrow(MM2)-1, 2]) & n0[j] < max.count)
    {
      p.n0 <- tmpfun(n0[j]+2); mu.n0 <- p.n0*(n0[j]+2); mu2.n0 <- p.n0*(n0[j]+2)^2
      MM2 <- rbind(MM2, c(p.n0, mu.n0, mu2.n0, n0[j]+2)); n0[j] <- n0[j]+1
    }

    MM2 <- MM2[-1, ];  MM.all <- rbind(MM1, MM2); weight <- 1/sum(MM.all[,1])

    if(!is.null(pred.interval)){
      #### Equal Tail (alpha/2) Prediction Interval
      pd <- cbind(MM.all[,4], MM.all[,1]*weight )
      pd1 <- pd[order(pd[,1]), ];pd1 <- cbind(pd1, cumsum(pd1[,2]))
      id1 <- suppressWarnings(ifelse(max(which(pd1[,3] <= (1-pred.interval)/2 ))==-Inf, 0,
                                     max(which(pd1[,3] <= (1-pred.interval)/2 ))))

      id2 <- min(which(pd1[,3] >= 1-(1-pred.interval)/2 ))
      L1 <- id1; U1 <- id2-1

      #### Shortest Length Prediction Interval
      pd2 <- pd[order(pd[,2], decreasing = TRUE),]
      pd2 <- cbind(pd2, cumsum(pd2[,2]))
      id3 <- which(pd2[,3] >= pred.interval)[1]
      L2 <- min(pd2[1:id3,1]); U2 <- max(pd2[1:id3,1])
      ans[j, ] <- c(sum(MM.all[,2])*weight, sum(MM.all[,3]*weight)-(sum(MM.all[,2])*weight)^2,
                    L1, U1, L2, U2)
    }else{
      ans[j, ] <- c(sum(MM.all[,2])*weight, sum(MM.all[,3]*weight)-(sum(MM.all[,2])*weight)^2)
    }
  }
  if(!is.null(pred.interval)){
    anslist <- (list(obs.locs = obs.locs,
                     obs.y = obs.y,
                     pred.locs = pred.locs,
                     predValue = ans[,1],
                     predCount = round(ans[,1]),
                     predVar = ans[,2],
                     ConfidenceLevel = pred.interval,
                     predInterval.EqualTail = ans[,3:4],
                     predInterval.Shortest = ans[,5:6]))
  }else{
    anslist <- (list(obs.locs = obs.locs,
                     obs.y = obs.y,
                     pred.locs = pred.locs,
                     predValue = ans[,1],
                     predCount = round(ans[,1]),
                     predVar = ans[,2]))
  }
  return(anslist)
}


#######################################################################################################################
#######################################################################################################################
#### OUTPUT FUNCTION: Prediction using GHK simulator: parallel version via snowfall.
#######################################################################################################################
#######################################################################################################################

predGHK.sf <- function(
    obs.y, obs.x = NULL, obs.locs, pred.x = NULL, pred.locs, longlat = FALSE, distscale = 1,
      marginal, corr, obs.effort = 1, pred.effort = 1, estpar = NULL,
      corrpar0 = NULL, pred.interval = NULL,
       n.cores = 2, cluster.type="SOCK", ghkoptions = list(nrep = c(100,1000), reorder = FALSE, seed = 12345)
){
  x <- as.data.frame(cbind(rep(1,length(obs.y)), obs.x))
  colnames(x)[1] <- c("Intercept")

  if(is.null(ghkoptions[["nrep"]])) ghkoptions$nrep = c(100,1000)
  if(is.null(ghkoptions[["reorder"]])) ghkoptions$reorder = FALSE
  if(is.null(ghkoptions[["seed"]])) ghkoptions$seed = 12345

  if(!is.matrix(pred.locs) & !is.data.frame(pred.locs))
    stop("Input 'pred.locs' must be a data frame or matrix!")

  if(!is.matrix(obs.locs) & !is.data.frame(obs.locs))
    stop("Input 'obs.locs' must be a data frame or matrix!")

  if(length(obs.effort) == 1) obs.effort <- rep(obs.effort, nrow(obs.locs))
  if(!length(obs.effort) == nrow(obs.locs))
    stop("Sampling Effort must be equal to the number of sampling locations!")

  if(length(pred.effort) == 1) pred.effort <- rep(pred.effort, nrow(pred.locs))

  if(!length(pred.effort) == nrow(pred.locs))
    stop("Prediction Effort must be equal to the number of prediction locations!")

  #### First calculate the MLE and output the log-likelihood as denominator
  if(is.null(estpar)){

  MLE.est <- mleGHK(y = obs.y, x = x[,-1], locs = obs.locs, marginal = marginal, corr = corr,
                    effort = obs.effort, longlat = longlat, distscale = distscale,
                    corrpar0 = corrpar0, ghkoptions = ghkoptions)

  loglik <- MLE.est$log.lik; estpar <- MLE.est$MLE

  }else{
  loglik <- -likGHK(pars = estpar, y = obs.y, x = x, locs = obs.locs, marginal = marginal, corr = corr,
                    effort = obs.effort, longlat = longlat, distscale = distscale,
                    nrep = ghkoptions$nrep[length(ghkoptions$nrep)],
                    reorder = ghkoptions$reorder, seed = ghkoptions$seed)
  }

  if(is.null(pred.x)) {
    pred.x <- matrix(1, nrow = nrow(pred.locs), ncol = 1)
  }else {
    pred.x <- cbind(rep(1, nrow(pred.x)) , pred.x)
  }
  pred.x <- as.data.frame(pred.x)
  colnames(pred.x)[1] <- c("Intercept")
  names(pred.x) <- names(x)

  if(nrow(pred.x)!= nrow(pred.locs))
    stop("Number of prediction locations did not match the number of covariates")

  if (requireNamespace("FNN", quietly = TRUE)) {
    if(nrow(pred.locs) == 1){
      indexloc <- which.min(FNN::get.knnx(pred.locs, obs.locs, 1)$nn.dist)
     m0 <- n0 <- round(unique(pred.effort)*obs.y[indexloc]/obs.effort[indexloc])+1
  }else{
    m0 <- n0 <- round(pred.effort*apply(pred.locs, 1, function(x) obs.y[which.min(FNN::get.knnx(t(as.matrix(x)),
                obs.locs, 1)$nn.dist)]/obs.effort[which.min(FNN::get.knnx(t(as.matrix(x)), obs.locs, 1)$nn.dist)]))+1
   }
  }else {
    stop("Please install {FNN} first!")
  }
  NPL <- length(m0)

  max.count <- ceiling(5*max(obs.y/obs.effort)*max(pred.effort))
  nnrep <- length(ghkoptions$nrep)

  #### Begin to parallel
  if (requireNamespace("snowfall", quietly = TRUE)) {

    snowfall::sfInit(parallel =TRUE, cpus = n.cores, type = cluster.type)
    suppressMessages(snowfall::sfExportAll(except = NULL, debug = FALSE))
    suppressMessages(snowfall::sfLibrary("gcKrig", character.only= TRUE))
    #snowfall::sfClusterSetupRNG()
    #  sfClusterEvalQ( ls() )


  par.pred.inner <- function(j){

    tmpfun <- function(xtmp) {
      exp(-likGHK(pars = estpar, y = c(obs.y, xtmp), x = rbind(as.matrix(x), pred.x[j,]),
                  locs = rbind(obs.locs, pred.locs[j,]), marginal = marginal, corr = corr,
                  effort = c(obs.effort, pred.effort[j]), longlat = longlat, distscale = distscale,
                  nrep = ghkoptions$nrep[nnrep], reorder = ghkoptions$reorder, seed = ghkoptions$seed) - loglik)
    }
    p.m0 <- p.n0 <- tmpfun(m0[j]); mu.m0 <- mu.n0 <- p.m0*m0[j];  mu2.m0 <- mu2.n0 <- p.m0*m0[j]^2

    MM1 <- matrix(0, nrow = 2, ncol = 4); MM2 <- matrix(0, nrow = 2, ncol = 4)

    MM1[1,] <- c(p.m0, mu.m0, mu2.m0, m0[j]);  MM2[1,] <- c(p.n0, mu.n0, mu2.n0, n0[j])

    p.m0 <- tmpfun(m0[j]-1); mu.m0 <- p.m0*(m0[j]-1);
    mu2.m0 <- p.m0*(m0[j]-1)^2;  MM1[2,] <- c(p.m0, mu.m0, mu2.m0, m0[j]-1)

    while( (p.m0 > sqrt(.Machine$double.eps) | MM1[nrow(MM1), 2] > MM1[nrow(MM1)-1, 2]) & m0[j] > 1)
    {
      p.m0 <- tmpfun(m0[j]-2); mu.m0 <- p.m0*(m0[j]-2); mu2.m0 <- p.m0*(m0[j]-2)^2
      MM1 <- rbind(MM1, c(p.m0, mu.m0, mu2.m0, m0[j]-2));   m0[j] <- m0[j]-1
    }
    #### Search from n0 to the right
    p.n0 <- tmpfun(n0[j]+1); mu.n0 <- p.n0*(n0[j]+1);  mu2.n0 <- p.n0*(n0[j]+1)^2
    MM2[2, ] <- c(p.n0, mu.n0, mu2.n0, n0[j]+1)

    while( (p.n0 > sqrt(.Machine$double.eps) | MM2[nrow(MM2), 2] > MM2[nrow(MM2)-1, 2]) & n0[j] < max.count)
    {
      p.n0 <- tmpfun(n0[j]+2); mu.n0 <- p.n0*(n0[j]+2);   mu2.n0 <- p.n0*(n0[j]+2)^2
      MM2 <- rbind(MM2, c(p.n0, mu.n0, mu2.n0, n0[j]+2));  n0[j] <- n0[j]+1
    }

    MM2 <- MM2[-1, ];  MM.all <- rbind(MM1, MM2); weight <- 1/sum(MM.all[,1])


    if(!is.null(pred.interval)){
      #### Equal Tail (alpha/2) Prediction Interval
      pd <- cbind(MM.all[,4], MM.all[,1]*weight )
      pd1 <- pd[order(pd[,1]), ];pd1 <- cbind(pd1, cumsum(pd1[,2]))
      id1 <- suppressWarnings(ifelse(max(which(pd1[,3] <= (1-pred.interval)/2 ))==-Inf, 0,
                                     max(which(pd1[,3] <= (1-pred.interval)/2 ))))

      id2 <- min(which(pd1[,3] >= 1-(1-pred.interval)/2 )); L1 <- id1; U1 <- id2-1

      #### Shortest Length Prediction Interval
      pd2 <- pd[order(pd[,2], decreasing = TRUE),]
      pd2 <- cbind(pd2, cumsum(pd2[,2]))
      id3 <- which(pd2[,3] >= pred.interval)[1]
      L2 <- min(pd2[1:id3,1]); U2 <- max(pd2[1:id3,1])
      ans <- c(sum(MM.all[,2])*weight, sum(MM.all[,3]*weight)-(sum(MM.all[,2])*weight)^2,
                    L1, U1, L2, U2)
    }else{
      ans <- c(sum(MM.all[,2])*weight, sum(MM.all[,3]*weight)-(sum(MM.all[,2])*weight)^2)
    }
    return(ans)
  }
  out = snowfall::sfSapply(1:NPL, par.pred.inner)
  snowfall::sfStop()
  }else{
      stop("Please install {snowfall} first before using this function!")
  }
  ans <- t(out)

  if(!is.null(pred.interval)){
    anslist <- (list(obs.locs = obs.locs,
                     obs.y = obs.y,
                     pred.locs = pred.locs,
                     predValue = ans[,1],
                     predCount = round(ans[,1]),
                     predVar = ans[,2],
                     ConfidenceLevel = pred.interval,
                     predInterval.EqualTail = ans[,3:4],
                     predInterval.Shortest = ans[,5:6]))
  }else{
    anslist <- (list(obs.locs = obs.locs,
                     obs.y = obs.y,
                     pred.locs = pred.locs,
                     predValue = ans[,1],
                     predCount = round(ans[,1]),
                     predVar = ans[,2]))
  }
  return(anslist)
}
