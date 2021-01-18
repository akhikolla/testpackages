
#' High-powered detection of genetic effects on DNA methylation
#'
#' Perform high-powered detection of genetic effects on DNA
#' methylation using integrated methylation QTL (methylation quantitative-trait locus)
#' mapping and allele-specific analysis.
#'
#'
#' @param geno 	a data list containing the genotype data
#' @param data a data list containing the methylation data
#' @param K a known kinship matrix. This matrix should be a positive semi-definite
#' matrix with dimensions equal to the sample size
#' @param Covariates a matrix containing the covariates subject to adjustment (Default = NULL)
#' @param numCore a positive integer specifying the number of cores for parallel computing (default = 1)
#' @param fit.maxiter a positive integer specifying the maximum number of iterations
#' when fitting the generalized linear mixed model (default = 500)
#' @param fit.tol a positive number specifying tolerance, the difference threshold
#' for parameter estimates below which iterations should be stopped (default = 1e-5)
#' @param verbose include verbose output
#'
#'
#' @return A \code{data.frame} containing the following named elements:
#' \itemize{
#'   \item \code{loc:} ordinal number of SNP-CpG pair being analyzed
#'   \item \code{numIDV:} number of observations of SNP-CpG pair being analyzed
#'   \item \code{beta:} the fixed effect parameter estimate for the predictor of interest
#'   \item \code{se_beta:} the standard deviation of fixed effect
#'   \item \code{pvalue:} P value for the fixed effect, based on the Wald test
#'   \item \code{h2:} heritability of the transformed rate
#'   \item \code{sigma2:} total variance component
#'   \item \code{converged:} a logical indicator for convergence
#'   \item \code{time:} time to converge
#' }
#'
#' @references Fan, Y., Vilgalys, T.P., Sun, S., Peng, Q., Tung, J. and Zhou, X.,
#' 2019. High-powered detection of genetic effects on DNA methylation using integrated
#' methylation QTL mapping and allele-specific analysis. bioRxiv, p.615039.
#'
#' @examples
#'
#' # This example demonstrates IMAGE:
#' \donttest{
#'   data(ExampleData)
#'   geno <- ExampleData$geno
#'   K <- ExampleData$K
#'   data <- ExampleData$data
#'   res=image(geno,data,K)
#' }
#'
#' \donttest{
#'   # We've saved the results of the example above to show an example of
#'   # the outputs IMAGE produces:
#'   data(example_results)
#' }
#'
#' # Toy example for testing purposes only:
#'
#' geno <- list()
#' geno$hap1 = matrix(sample(c(0,1),25, replace = TRUE, prob = c(0.6,0.4)),
#'                     nrow = 5, ncol = 5)
#' geno$hap2 = matrix(sample(c(0,1),25, replace = TRUE, prob = c(0.6,0.4)),
#'                     nrow = 5, ncol = 5)
#'
#' data <- list()
#' data$r = matrix(sample(0:50,25, replace = TRUE), nrow = 5, ncol = 5)
#' data$y = matrix(sample(0:50,25, replace = TRUE), nrow = 5, ncol = 5)
#' data$r1 = matrix(sample(0:50,25, replace = TRUE), nrow = 5, ncol = 5)
#' data$r2 = matrix(sample(0:50,25, replace = TRUE), nrow = 5, ncol = 5)
#' data$y1 = matrix(sample(0:50,25, replace = TRUE), nrow = 5, ncol = 5)
#' data$y2 = matrix(sample(0:50,25, replace = TRUE), nrow = 5, ncol = 5)
#'
#' K = matrix(runif(25,-0.1,0.1), nrow = 5, ncol = 5)
#'
#' res=image(geno,data,K)
#'
#' \dontshow{
#'   closeAllConnections()
#' }
#'
#' @import parallel
#' @import foreach
#' @import doParallel
#' @import Matrix
#' @importFrom Rcpp evalCpp
#' @importFrom Rcpp sourceCpp
#' @useDynLib IMAGE
#' @importFrom stats binomial
#' @importFrom stats glm
#' @importFrom stats lm
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats na.omit
#' @importFrom stats na.pass
#' @importFrom stats pchisq
#' @importFrom stats poisson
#' @importFrom stats var
#' @export image
#'

image <- function(geno,data,K,Covariates=NULL,numCore=1,fit.maxiter=500,fit.tol=1e-5,verbose=TRUE) {

  if(numCore > 1){
    if(numCore>detectCores()){warning("PQLseq:: the number of cores you're setting is larger than detected cores!");numCore = detectCores()-1}
  }

  m=ncol(geno$hap1)
  n=nrow(geno$hap1)
  registerDoParallel(cores=numCore)

  # cleanup connections
  #on.exit(expr = closeAllConnections(), add = FALSE)

  if(is.null(Covariates)){
    numCov <- 0
  }else{
    numCov     <- dim(Covariates)[2]
    Covariates <- as.matrix(Covariates)
  }
  iVar <- NULL # made explicit to apease R CMD CHECK
  resBMM <- foreach(iVar=1:m, .combine=rbind, .errorhandling = 'remove')%dopar%{

    res <- data.frame()


    if(iVar%%100==0&&verbose)
    {
    cat(iVar, '\n')
    }
    geno1 = geno$hap1[, iVar]
    geno2 = geno$hap2[, iVar]
    heter  = which(geno1 + geno2 == 1)
    if (length(heter) != 0) {
      homo = (1:n)[-heter]
      ratio <- c(data$y[homo, iVar])/c(data$r[homo, iVar])
      ratio[is.na(ratio)] <- 1.1
      idx <- which(ratio > 1.0)
      if (length(idx) > 0) {
        homo = homo[-idx]
      }
      ratio1 <- c(data$y1[heter, iVar])/c(data$r1[heter, iVar])
      ratio1[is.na(ratio1)] <- 1.1
      ratio2 <- c(data$y2[heter, iVar])/c(data$r2[heter, iVar])
      ratio2[is.na(ratio2)] <- 1.1
      idx <- union( which(ratio1 > 1.0), which(ratio2 > 1.0) )
      if (length(idx) > 0) {
        heter = heter[-idx]
      }
      n1 <- length(homo)
      n2 <- length(heter)
      gheter = c(rbind(geno1[heter], geno2[heter]))
      genotypes = c(geno1[homo], gheter)
      if(numCov>0)
      {
        if(length(heter)>1)
        {
          covariates=rbind(Covariates[homo,],t(kronecker(t(Covariates[heter,]),matrix(1, nrow=1, ncol=2))))
        }else{
          covariates=rbind(Covariates[homo,],t(kronecker((Covariates[heter,]),matrix(1, nrow=1, ncol=2))))
        }
      }
      CountData <- c(data$y[homo, iVar], c(rbind(data$y1[heter, iVar], data$y2[heter, iVar])))
      LibSize <- c(data$r[homo, iVar], c(rbind(data$r1[heter, iVar], data$r2[heter, iVar])))
    } else {
      homo <- 1:n
      ratio <- c(data$y[homo, iVar])/c(data$r[homo, iVar])
      ratio[is.na(ratio)] <- 1.1
      idx <-which(ratio > 1.0)
      if (length(idx) > 0) {
        homo = homo[-idx]
      }
      n1 <- length(homo)
      n2 <- 0
      CountData <- c(data$y[homo, iVar])
      LibSize <- c(data$r[homo, iVar])
      pheno1 <- NULL # made explicit to apease R CMD CHECK
      genotypes <- c(pheno1[homo])
      if(numCov>0)
      {
        covariates=Covariates[homo,]
      }
    }

    numIDV <- length(homo)+length(heter)
    ratio <- CountData/LibSize

    modified_data=data_modify(genotypes,ratio,LibSize,CountData)
    LibSize=modified_data$r
    CountData=modified_data$y
    ratio <- CountData/LibSize

    if(numCov>0){
      idxcov=numeric()
      for(icov in 1:numCov)
      {
        if(length(unique(covariates[,icov])) == 1){
          idxcov=c(idxcov,icov)

        }
      }
      covariates<- covariates[,-idxcov]

      numCov     <- dim(covariates)[2]
      covariates=as.matrix(covariates)

      model0<-lm(formula=genotypes~covariates)
      genotypes=genotypes-cbind(rep(1,length(genotypes)),covariates)%*%coef(summary(model0))[,1]

    }


    t1 <- system.time( model0 <- glm(formula = ratio~genotypes,
                                     family = binomial(link = "logit"),
                                     weights = LibSize) )
    coef <- coef(summary(model0))

    ###############################################################


    ###############################################################

    K11 <- K[homo, homo]
    K22 <- K[heter, heter]
    K12 <- K[homo, heter]
    K21 <- K[heter, homo]
    if(length(homo)==1)
    {
      RelatednessMatrix <- cbind(K11, kronecker(t(K12), matrix(1, nrow=1, ncol=2)))
    }else{
      RelatednessMatrix <- cbind(K11, kronecker(K12, matrix(1, nrow=1, ncol=2)))
    }
    if(length(heter)==1)
    {
      RelatednessMatrix <- rbind(RelatednessMatrix, cbind(kronecker(t(K21), matrix(1, nrow=2, ncol=1)),
                                                          kronecker(K22, matrix(1, nrow=2, ncol=2))))
    }else{
      RelatednessMatrix <- rbind(RelatednessMatrix, cbind(kronecker(K21, matrix(1, nrow=2, ncol=1)),
                                                          kronecker(K22, matrix(1, nrow=2, ncol=2))))
    }
    eig               <- eigen(RelatednessMatrix)
    eigval            <- eig$value
    eigvector         <- eig$vectors
    if(any(eigval<1e-10)){
      #   warning("PQLseq::the relatedness matrix is singular, it has been modified!")
      RelatednessMatrix <- as.matrix(nearPD(RelatednessMatrix,corr=T)$mat)
    }
    rm(eig)
    rm(eigval)
    rm(eigvector)
    RelatednessMatrix <- list( RelatednessMatrix, diag(c(rep(0.5, n1), rep(1, 2*n2))))
    model0$numTotal <- LibSize
    model0$numSucc  <- CountData
    tmpRelatednessMatrix=RelatednessMatrix
    names(tmpRelatednessMatrix) <- paste("kins", 1:length(tmpRelatednessMatrix), sep="")



    #############################################################################################
    GRM2 <- as.matrix(bdiag(diag(1, nrow=n1),
                            kronecker(diag(nrow=n2), matrix(1, nrow=2, ncol=2)) ))
    RelatednessMatrix[[3]]=GRM2
    beta <- tau1 <- tau2 <- se_beta <- pvalue <- converged <- h2 <- sigma2 <- overdisp <- NA

    model0$numTotal <- LibSize
    model0$numSucc  <- CountData
    tmpRelatednessMatrix=RelatednessMatrix
    names(tmpRelatednessMatrix) <- paste("kins", 1:length(tmpRelatednessMatrix), sep="")
    t1 <- system.time(model1 <- try( PQLseq.fit(model0, tmpRelatednessMatrix,maxiter = fit.maxiter, tol = fit.tol) ))
    if(class(model1) != "try-error"){
      # numIDV <- length(idx)
      beta        <- model1$coefficients[ length(model1$coefficients) ]# the last one
      se_beta     <- sqrt( diag(model1$cov)[ length(model1$coefficients) ] )
      pvalue      <- pchisq( (beta/se_beta)^2, 1, lower.tail = F)
      sigma2      <- model1$theta[2]+model1$theta[3]
      h2          <- model1$theta[2]/(sigma2)
      tau1        <- model1$theta[2]
      tau2        <- model1$theta[3]
      rho<-model1$theta[4]/(model1$theta[3]+model1$theta[4])
      converged   <- model1$converged
    }else{converged <- FALSE}

    if(!converged)
    {
      tmpRelatednessMatrix=tmpRelatednessMatrix[-3]
      names(tmpRelatednessMatrix) <- paste("kins", 1:length(tmpRelatednessMatrix), sep="")
      t1 <- system.time(model1 <- try( PQLseq2.fit(model0, tmpRelatednessMatrix,maxiter = fit.maxiter, tol = fit.tol) ))
      if(class(model1) != "try-error"){
        # numIDV <- length(idx)
        beta        <- model1$coefficients[ length(model1$coefficients) ]# the last one
        se_beta     <- sqrt( diag(model1$cov)[ length(model1$coefficients) ] )
        pvalue      <- pchisq( (beta/se_beta)^2, 1, lower.tail = F)
        sigma2      <- model1$theta[2]+model1$theta[3]
        h2          <- model1$theta[2]/(sigma2)
        tau1        <- model1$theta[2]
        tau2        <- model1$theta[3]
        #  rho<-model1$theta[4]/(model1$theta[3]+model1$theta[4])
        rho=model1$tau.rho
        converged   <- model1$converged
      }else{converged <- FALSE}
    }

    res <- rbind(res,data.frame(Loc=iVar,numIDV = numIDV, beta = beta, se_beta = se_beta,
                                pvalue = pvalue, h2 = h2, sigma2 = sigma2,
                                converged = converged,time=t1[3]))





    return(res)

  }
}


data_modify <- function(genotypes,ratio,LibSize,CountData) {
  data=list()
  if(mean(ratio[which(genotypes==1)])==0)
  {
    idx1=which(genotypes==1)
    idx2=which(LibSize[idx1]==max(LibSize[idx1]))
    idx2=idx2[1]
    CountData[idx1[idx2]]=CountData[idx1[idx2]]+1
    ratio <- CountData/LibSize
  }else if((mean(ratio[which(genotypes==1)])==1)){
    idx1=which(genotypes==1)
    idx2=which(LibSize[idx1]==min(LibSize[idx1]))
    idx2=idx2[1]
    CountData[idx1[idx2]]=CountData[idx1[idx2]]-1
    ratio <- CountData/LibSize
  }else if((mean(ratio[which(genotypes==0)])==0)){
    idx1=which(genotypes==0)
    idx2=which(LibSize[idx1]==max(LibSize[idx1]))
    idx2=idx2[1]
    CountData[idx1[idx2]]=CountData[idx1[idx2]]+1
    ratio <- CountData/LibSize
  }else if((mean(ratio[which(genotypes==0)])==1)){
    idx1=which(genotypes==0)
    idx2=which(LibSize[idx1]==max(LibSize[idx1]))
    idx2=idx2[1]
    CountData[idx1[idx2]]=CountData[idx1[idx2]]-1
    ratio <- CountData/LibSize
  }
  data$y=CountData
  data$r=LibSize
  return(data)
}



##########################################################
#           	   PQLseq FIT FUNCTION					 #
##########################################################
PQLseq.fit <- function(model0, RelatednessMatrix, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, verbose = FALSE) {

  names(RelatednessMatrix) <- paste("kins", 1:length(RelatednessMatrix), sep="")
  # if((method.optim == "AI")&(!sum(model0$fitted.values<1e-5))) {
  if(method.optim == "AI") {
    fixtau.old <- rep(0, length(RelatednessMatrix)+1)
    # to use average information method to fit alternative model
    model1 <- PQLseq.AI(model0, RelatednessMatrix, maxiter = maxiter, tol = tol, verbose = verbose)
    fixtau.new <- 1*(model1$theta < 1.01 * tol)

    while(any(fixtau.new != fixtau.old)) {
      fixtau.old <- fixtau.new
      # to use average information method to fit alternative model
      model1 <- PQLseq.AI(model0, RelatednessMatrix, fixtau = fixtau.old, maxiter = maxiter, tol = tol, verbose = verbose)
      fixtau.new <- 1*(model1$theta < 1.01 * tol)
    }
    return(model1)
  }
}

##########################################################
#       PQLseq FIT AVERAGE INFORMATION FUNCTION			 #
##########################################################
PQLseq.AI <- function(model0, RelatednessMatrix, tau = rep(0, length(RelatednessMatrix)+1), fixtau = rep(0, length(RelatednessMatrix)+1), maxiter = 500, tol = 1e-5, verbose = FALSE) {

  if(model0$family$family %in% c("binomial")){
    y <- model0$numSucc
  }else{
    y <- model0$y
  }
  numIDV <- length(y)
  offset <- model0$offset
  if(is.null(offset)) {offset <- rep(0, numIDV)}

  family <- model0$family
  eta <- model0$linear.predictors
  mu <- model0$fitted.values
  mu.eta <- family$mu.eta(eta)
  D <- mu.eta/sqrt(model0$family$variance(mu))

  if(family$family %in% c("binomial")){
    mu.eta <- model0$numTotal*mu.eta
    D <- mu.eta/sqrt(model0$numTotal*model0$family$variance(mu))
    mu <- model0$numTotal*mu
  }

  Y <- eta - offset + (y - mu)/mu.eta
  X <- model.matrix(model0)
  alpha <- model0$coef

  if(family$family %in% c("poisson", "binomial")) {
    tau[1] <- 1
    fixtau[1] <- 1
  }
  numK <- length(RelatednessMatrix)
  idxtau <- which(fixtau == 0)
  numK2 <- sum(fixtau == 0)
  if(numK2 > 0) {
    tau[fixtau == 0] <- rep(min(0.9,var(Y)/(numK+1)), numK2)

    H <- tau[1]*diag(1/D^2)
    for(ik in 1:numK) {H <- H + tau[ik+1]*RelatednessMatrix[[ik]]}

    Hinv <- chol2inv(chol(H))
    HinvX <- crossprod(Hinv, X)
    XHinvX <- crossprod(X, HinvX)

    P <- try(Hinv - tcrossprod(tcrossprod(HinvX, chol2inv(chol( XHinvX ))), HinvX))

    if(class(P) == "try-error"){
      stop("Error in P matrix calculation!")
    }

    PY <- crossprod(P, Y)
    tau0 <- tau
    for(ik in 1:numK2) {
      if(ik == 1 && fixtau[1] == 0) tau[1] <- max(0, tau0[1] + tau0[1]^2 * (sum((PY/D)^2) - sum(diag(P)/D^2))/numIDV)
      else {
        PAPY <- crossprod(P, crossprod(RelatednessMatrix[[idxtau[ik]-1]], PY))
        tau[idxtau[ik]] <- max(0, tau0[idxtau[ik]] + tau0[idxtau[ik]]^2 * (crossprod(Y, PAPY) - sum(P*RelatednessMatrix[[idxtau[ik]-1]]))/numIDV)
      }
    }
  }

  for (iter in seq_len(maxiter)) {
    alpha0 <- alpha
    tau0 <- tau
    model1 <- AI(Y, X, length(RelatednessMatrix), RelatednessMatrix, D^2, tau, fixtau, tol)

    tau <- as.numeric(model1$tau)
    cov <- as.matrix(model1$cov)
    alpha <- as.numeric(model1$alpha)
    eta <- as.numeric(model1$eta) + offset


    mu <- family$linkinv(eta)
    mu.eta <- family$mu.eta(eta)
    D <- mu.eta/sqrt(family$variance(mu))

    if(family$family %in% c("binomial")){
      mu.eta <- model0$numTotal*mu.eta
      D <- mu.eta/sqrt(model0$numTotal*family$variance(mu))
      mu <- model0$numTotal*mu
    }

    Y <- eta - offset + (y - mu)/mu.eta

    if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) {break}
    if(max(tau) > tol^(-2)|any(is.infinite(D))|any(is.infinite(mu))|any(is.infinite(eta)) ) {

      iter <- maxiter
      break
    }
  }

  converged <- ifelse(iter < maxiter, TRUE, FALSE)

  if(converged==FALSE)
  {
    tau = rep(0, length(RelatednessMatrix)+1)
    tau[1] <- 1
    if(model0$family$family %in% c("binomial")){
      y <- model0$numSucc
    }else{
      y <- model0$y
    }
    numIDV <- length(y)
    offset <- model0$offset
    if(is.null(offset)) {offset <- rep(0, numIDV)}

    family <- model0$family
    eta <- model0$linear.predictors
    mu <- model0$fitted.values
    mu.eta <- family$mu.eta(eta)
    D <- mu.eta/sqrt(model0$family$variance(mu))

    if(family$family %in% c("binomial")){
      mu.eta <- model0$numTotal*mu.eta
      D <- mu.eta/sqrt(model0$numTotal*model0$family$variance(mu))
      mu <- model0$numTotal*mu
    }

    Y <- eta - offset + (y - mu)/mu.eta
    alpha <- model0$coef



    if(numK2 > 0) {

      H <- tau[1]*diag(1/D^2)
      for(ik in 1:numK) {H <- H + tau[ik+1]*RelatednessMatrix[[ik]]}

      Hinv <- chol2inv(chol(H))
      HinvX <- crossprod(Hinv, X)
      XHinvX <- crossprod(X, HinvX)

      P <- try(Hinv - tcrossprod(tcrossprod(HinvX, chol2inv(chol( XHinvX ))), HinvX))

      if(class(P) == "try-error"){
        stop("Error in P matrix calculation!")
      }
      PY <- crossprod(P, Y)
      tau0 <- tau
    }

    for (iter in seq_len(maxiter)) {
      alpha0 <- alpha
      tau0 <- tau
      model1 <- AI(Y, X, length(RelatednessMatrix), RelatednessMatrix, D^2, tau, fixtau, tol)

      tau <- as.numeric(model1$tau)
      cov <- as.matrix(model1$cov)
      alpha <- as.numeric(model1$alpha)
      eta <- as.numeric(model1$eta) + offset


      mu <- family$linkinv(eta)
      mu.eta <- family$mu.eta(eta)
      D <- mu.eta/sqrt(family$variance(mu))

      if(family$family %in% c("binomial")){
        mu.eta <- model0$numTotal*mu.eta
        D <- mu.eta/sqrt(model0$numTotal*family$variance(mu))
        mu <- model0$numTotal*mu
      }

      Y <- eta - offset + (y - mu)/mu.eta

      if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tau - tau0)/(abs(tau) + abs(tau0) + tol)) < tol) {break}
      if(max(tau) > tol^(-2)|any(is.infinite(D))|any(is.infinite(mu))|any(is.infinite(eta)) ) {

        iter <- maxiter
        break
      }
    }

  }

  converged <- ifelse(iter < maxiter, TRUE, FALSE)

  res <- y - mu
  P <- model1$P
  return(list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, P = P, residuals = res, cov = cov, converged = converged))
}# end function


#########################################
#             CODE END                  #
#########################################



##########################################################
#           	   PQLseq2 FIT FUNCTION					 #
##########################################################
PQLseq2.fit <- function(model0, RelatednessMatrix, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, verbose = FALSE) {

  names(RelatednessMatrix) <- paste("kins", 1:length(RelatednessMatrix), sep="")
  # if((method.optim == "AI")&(!sum(model0$fitted.values<1e-5))) {
  if(method.optim == "AI") {
    fixtau.old <- rep(0, length(RelatednessMatrix)+1)
    # to use average information method to fit alternative model
    model1 <- PQLseq.AI2(model0, RelatednessMatrix, maxiter = maxiter, tol = tol, verbose = verbose)
    fixtau.new <- 1*(model1$theta < 1.01 * tol)




    ##########################




    while(any(fixtau.new != fixtau.old)) {
      fixtau.old <- fixtau.new
      # to use average information method to fit alternative model
      model1 <- PQLseq.AI2(model0, RelatednessMatrix, fixtau = fixtau.old, maxiter = maxiter, tol = tol, verbose = verbose)






      fixtau.new <- 1*(model1$theta < 1.01 * tol)
    }
    return(model1)
  }
}

##########################################################
#       PQLseq2 FIT AVERAGE INFORMATION FUNCTION			 #
##########################################################
PQLseq2.AI <- function(model0, RelatednessMatrix, tau = rep(0, length(RelatednessMatrix)+1), fixtau = rep(0, length(RelatednessMatrix)+1), maxiter = 500, tol = 1e-5, verbose = FALSE) {

  n=nrow(RelatednessMatrix[[1]])

  if(model0$family$family %in% c("binomial")){
    y <- model0$numSucc
  }else{
    y <- model0$y
  }
  numIDV <- length(y)
  offset <- model0$offset
  if(is.null(offset)) {offset <- rep(0, numIDV)}

  family <- model0$family
  eta <- model0$linear.predictors
  mu <- model0$fitted.values
  mu.eta <- family$mu.eta(eta)
  D <- mu.eta/sqrt(model0$family$variance(mu))

  if(family$family %in% c("binomial")){
    mu.eta <- model0$numTotal*mu.eta
    D <- mu.eta/sqrt(model0$numTotal*model0$family$variance(mu))
    mu <- model0$numTotal*mu
  }

  Y <- eta - offset + (y - mu)/mu.eta
  X <- model.matrix(model0)
  alpha <- model0$coef

  if(family$family %in% c("poisson", "binomial")) {
    tau[1] <- 1
    fixtau[1] <- 1
  }
  numK <- length(RelatednessMatrix)
  idxtau <- which(fixtau == 0)
  numK2 <- sum(fixtau == 0)
  if(numK2 > 0) {
    tau[fixtau == 0] <- rep(min(0.9,var(Y)/(numK+1)), numK2)

    H <- tau[1]*diag(1/D^2)
    for(ik in 1:numK) {H <- H + tau[ik+1]*RelatednessMatrix[[ik]]}

    Hinv <- chol2inv(chol(H))
    HinvX <- crossprod(Hinv, X)
    XHinvX <- crossprod(X, HinvX)

    P <- try(Hinv - tcrossprod(tcrossprod(HinvX, chol2inv(chol( XHinvX ))), HinvX))

    if(class(P) == "try-error"){
      stop("Error in P matrix calculation!")
    }

    PY <- crossprod(P, Y)
    tau0 <- tau
    for(ik in 1:numK2) {
      if(ik == 1 && fixtau[1] == 0) tau[1] <- max(0, tau0[1] + tau0[1]^2 * (sum((PY/D)^2) - sum(diag(P)/D^2))/numIDV)
      else {
        PAPY <- crossprod(P, crossprod(RelatednessMatrix[[idxtau[ik]-1]], PY))
        tau[idxtau[ik]] <- max(0, tau0[idxtau[ik]] + tau0[idxtau[ik]]^2 * (crossprod(Y, PAPY) - sum(P*RelatednessMatrix[[idxtau[ik]-1]]))/numIDV)
      }
    }
  }

  d=0.1
  rho=seq(0,1,d)
  GRM1 <- as.matrix(bdiag(diag(0.5, nrow=n)))
  GRM2 <- as.matrix(bdiag(diag(0.5, nrow=n)))
  tmpGRM<-list(RelatednessMatrix[[1]],GRM1,GRM2)
  tau.rho=0

  for (iter in seq_len(maxiter)) {
    alpha0 <- alpha
    tau0 <- tau
    rho0<-tau.rho
    model1 <- AI(Y, X, length(RelatednessMatrix), RelatednessMatrix, D^2, tau, fixtau, tol)

    tau <- as.numeric(model1$tau)
    cov <- as.matrix(model1$cov)
    alpha <- as.numeric(model1$alpha)
    eta <- as.numeric(model1$eta) + offset


    mu <- family$linkinv(eta)
    mu.eta <- family$mu.eta(eta)
    D <- mu.eta/sqrt(family$variance(mu))

    if(family$family %in% c("binomial")){
      mu.eta <- model0$numTotal*mu.eta
      D <- mu.eta/sqrt(model0$numTotal*family$variance(mu))
      mu <- model0$numTotal*mu
    }

    Y <- eta - offset + (y - mu)/mu.eta


    numK=length(tmpGRM)
    X <- model.matrix(model0)

    likelihood=numeric(length=length(rho))
    H <- tau[1]*diag(1/D^2)
    for(ik in 1:(numK-1)) {H <- H + tau[ik+1]*tmpGRM[[ik]]}
    for(rhoi in 1:length(rho))
    {
      H1<-H+rho[rhoi]*tau[3]*tmpGRM[[3]]
      #  Px=model1$P
      Hinv <- chol2inv(chol(H1))
      HinvX <- crossprod(Hinv, X)
      XHinvX <- crossprod(X, HinvX)

      P <- try(Hinv - tcrossprod(tcrossprod(HinvX, chol2inv(chol( XHinvX ))), HinvX))
      likelihood[rhoi]=-0.5*(log(det(H1))+log(det(diag(D^2)))+t(Y)%*%P%*%Y)
    }

    tau.rho=rho[which(likelihood==max(likelihood))][1]
    #  if(d>=0.025)
    #  {
    #  d=d*0.5
    #  }
    #  rho=seq(tau.rho-5*d,tau.rho+5*d,d)
    #  rho=rho[which(rho<=1&rho>=-1)]


    tautmp=c(tau,tau.rho)
    tau0tmp=c(tau0,rho0)

    RelatednessMatrix<-list(RelatednessMatrix[[1]],GRM1+tau.rho*GRM2)
    names(RelatednessMatrix) <- paste("kins", 1:length(RelatednessMatrix), sep="")
    if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tautmp - tau0tmp)/(abs(tautmp) + abs(tau0tmp) + tol)) < tol) {break}
    if(max(tau) > tol^(-2)|any(is.infinite(D))|any(is.infinite(mu))|any(is.infinite(eta)) ) {

      iter <- maxiter
      break
    }
  }

  converged <- ifelse(iter < maxiter, TRUE, FALSE)

  if(converged==FALSE)
  {
    tau = rep(0, length(RelatednessMatrix)+1)
    tau[1] <- 1
    if(model0$family$family %in% c("binomial")){
      y <- model0$numSucc
    }else{
      y <- model0$y
    }
    numIDV <- length(y)
    offset <- model0$offset
    if(is.null(offset)) {offset <- rep(0, numIDV)}

    family <- model0$family
    eta <- model0$linear.predictors
    mu <- model0$fitted.values
    mu.eta <- family$mu.eta(eta)
    D <- mu.eta/sqrt(model0$family$variance(mu))

    if(family$family %in% c("binomial")){
      mu.eta <- model0$numTotal*mu.eta
      D <- mu.eta/sqrt(model0$numTotal*model0$family$variance(mu))
      mu <- model0$numTotal*mu
    }

    Y <- eta - offset + (y - mu)/mu.eta
    X <- model.matrix(model0)
    alpha <- model0$coef

    if(family$family %in% c("poisson", "binomial")) {
      tau[1] <- 1
      fixtau[1] <- 1
    }
    numK <- length(RelatednessMatrix)
    idxtau <- which(fixtau == 0)
    numK2 <- sum(fixtau == 0)
    if(numK2 > 0) {


      H <- tau[1]*diag(1/D^2)
      for(ik in 1:numK) {H <- H + tau[ik+1]*RelatednessMatrix[[ik]]}

      Hinv <- chol2inv(chol(H))
      HinvX <- crossprod(Hinv, X)
      XHinvX <- crossprod(X, HinvX)

      P <- try(Hinv - tcrossprod(tcrossprod(HinvX, chol2inv(chol( XHinvX ))), HinvX))

      if(class(P) == "try-error"){
        stop("Error in P matrix calculation!")
      }

      PY <- crossprod(P, Y)
      tau0 <- tau
    }

    d=0.1
    rho=seq(0,0,d)
    GRM1 <- as.matrix(bdiag(diag(0.5, nrow=n) ))
    GRM2 <- as.matrix(bdiag(diag(0.5, nrow=n) ))
    tmpGRM<-list(RelatednessMatrix[[1]],GRM1,GRM2)
    tau.rho=0

    for (iter in seq_len(maxiter)) {
      alpha0 <- alpha
      tau0 <- tau
      rho0<-tau.rho
      model1 <- AI(Y, X, length(RelatednessMatrix), RelatednessMatrix, D^2, tau, fixtau, tol)

      tau <- as.numeric(model1$tau)
      cov <- as.matrix(model1$cov)
      alpha <- as.numeric(model1$alpha)
      eta <- as.numeric(model1$eta) + offset


      mu <- family$linkinv(eta)
      mu.eta <- family$mu.eta(eta)
      D <- mu.eta/sqrt(family$variance(mu))

      if(family$family %in% c("binomial")){
        mu.eta <- model0$numTotal*mu.eta
        D <- mu.eta/sqrt(model0$numTotal*family$variance(mu))
        mu <- model0$numTotal*mu
      }

      Y <- eta - offset + (y - mu)/mu.eta


      numK=length(tmpGRM)
      X <- model.matrix(model0)

      likelihood=numeric(length=length(rho))
      H <- tau[1]*diag(1/D^2)
      for(ik in 1:(numK-1)) {H <- H + tau[ik+1]*tmpGRM[[ik]]}
      for(rhoi in 1:length(rho))
      {
        H1<-H+rho[rhoi]*tau[3]*tmpGRM[[3]]
        #  Px=model1$P
        Hinv <- chol2inv(chol(H1))
        HinvX <- crossprod(Hinv, X)
        XHinvX <- crossprod(X, HinvX)

        P <- try(Hinv - tcrossprod(tcrossprod(HinvX, chol2inv(chol( XHinvX ))), HinvX))
        likelihood[rhoi]=-0.5*(log(det(H1))+log(det(diag(D^2)))+t(Y)%*%P%*%Y)
      }

      tau.rho=rho[which(likelihood==max(likelihood))][1]



      tautmp=c(tau,tau.rho)
      tau0tmp=c(tau0,rho0)

      RelatednessMatrix<-list(RelatednessMatrix[[1]],GRM1+tau.rho*GRM2)
      names(RelatednessMatrix) <- paste("kins", 1:length(RelatednessMatrix), sep="")
      if(2*max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol), abs(tautmp - tau0tmp)/(abs(tautmp) + abs(tau0tmp) + tol)) < tol) {break}
      if(max(tau) > tol^(-2)|any(is.infinite(D))|any(is.infinite(mu))|any(is.infinite(eta)) ) {

        iter <- maxiter
        break
      }
    }
  }
  converged <- ifelse(iter < maxiter, TRUE, FALSE)

  res <- y - mu
  P <- model1$P
  return(list(theta = tau, coefficients = alpha, linear.predictors = eta, fitted.values = mu, Y = Y, P = P, residuals = res, cov = cov, converged = converged))
}# end function



#########################################
#             CODE END                  #
#########################################

## quiets concerns of R CMD check re: the .'s that appear in pipelines
utils::globalVariables(c("PQLseq.AI2", "iVar"))
