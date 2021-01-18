# Joint Modeling Main Function with NMRE (nonparametric Multiplicative random effects)

jmodelMult <- function (fitLME, fitCOX, data, model = 1, rho = 0, timeVarY = NULL, timeVarT = NULL, 
                        control = list(), ...) 
{
  cat("Running jmodelMult(), may take some time to finish.\n")
  call <- match.call()

  CheckInputs(fitLME, fitCOX, rho)

  cntrlLst <- GenerateControlList(control, 1)
  
  ID1 <- as.vector(unclass(fitLME$groups[[1]]))
  uniqueID <- !duplicated(ID1)
  tempID <- which(uniqueID)
  tempID <- c(tempID, length(ID1) + 1)
  ni <- diff(tempID) 
  ID <- rep(1:sum(uniqueID), times = ni)           
  nLong <- length(ni)
  
  if (ncol(fitCOX$y) != 3)
    stop("\n must fit time-dependent Cox model in coxph().")
  start <- as.vector(fitCOX$y[, 1])
  stop <- as.vector(fitCOX$y[, 2])
  event <- as.vector(fitCOX$y[, 3])
  Time <- stop[cumsum(ni)]
  d<- event[cumsum(ni)] 
  nSurv <- length(Time) 
  if (sum(d) < 5)
    warning("\n more than 5 events are required.")
  if (nLong != nSurv)
    stop("\n sample sizes in the longitudinal and event processes differ.")
  
  Z <- as.matrix(fitCOX$x)
  ncz <- ncol(Z)
  
  phi.names <- colnames(Z)
  formSurv <- formula(fitCOX)
  TermsSurv <- fitCOX$terms
  mfSurv <- model.frame(TermsSurv, data)[cumsum(ni), ]
  if (!is.null(timeVarT)) {
    if (!all(timeVarT %in% names(mfSurv)))
      stop("\n'timeVarT' does not correspond columns in the fixed-effect design matrix of 'fitCOX'.")
    mfSurv[timeVarT] <- Time
  }
  if (ncz > 0) {
    Ztime <- as.matrix(model.matrix(formSurv, mfSurv))
    if(attr(TermsSurv, 'intercept')) Ztime <- as.matrix(Ztime[, - 1])
    # design matrix of the covariates in survival part, one row for each subject, excluding intercept #
  } else Ztime <- matrix(, ncol = 0, nrow = nSurv)

  
  TermsLongX <- fitLME$terms
  mydata <- fitLME$data[all.vars(TermsLongX)] 
  formLongX <- formula(fitLME) 
  mfLongX <- model.frame(TermsLongX, data = mydata) 
  B <- as.matrix(model.matrix(formLongX, mfLongX)) # may include a column of intercept #
  Y <- as.vector(model.response(mfLongX, "numeric"))
  # give the column in mfLongX which is considered as response, may be transformed, of length N #
  alpha.name <- rownames(attr(TermsLongX, "factors"))[attr(TermsLongX, "response")]
  
  data.id <- mydata[uniqueID, ] # pick the first row of each subject in mydata, nrow=n #
  
  if (!is.null(timeVarY)) {
    if (!all(timeVarY %in% names(mydata)))
      stop("\n'timeVarY' does not correspond to columns in the fixed-effect design matrix of 'fitLME'.")
    data.id[timeVarY] <- Time
  }
  
  mfLongX.id <- model.frame(TermsLongX, data = data.id)
  Btime <- as.matrix(model.matrix(formLongX, mfLongX.id)) # same structure with B, but with only n rows #
  
  U <- sort(unique(Time[d == 1])) # ordered uncensored observed event time #
  tempU <- lapply(Time, function(t) U[t >= U]) 
  times <- unlist(tempU) # vector of length M #
  nk <- sapply(tempU, length)  # length of each element in times, vector of length n #
  M <- sum(nk)
  Index <- rep(1:nLong, nk) # repeat 1:n by nk, length M #
  Index0 <- match(Time, U)
  # vector of length n, for the uncensored subjects, return the value as nk, otherwise return NA # 
  Index1 <- unlist(lapply(nk[nk != 0], seq, from = 1)) # vector of length M #
  Index2 <- colSums(d * outer(Time, U, "==")) # vector of length nu #
  
  data.id2 <- data.id[Index, ]
  if (!is.null(timeVarY)) {
    data.id2[timeVarY] <- times
  }
  mfLongX2 <- model.frame(TermsLongX, data = data.id2)
  Btime2 <- as.matrix(model.matrix(formLongX, mfLongX2))

  mfSurv2 <- mfSurv[Index, ]
  if (!is.null(timeVarT)) {
    mfSurv2[timeVarT] <- times
  }
  if (ncz > 0) {
    Ztime2 <- as.matrix(model.matrix(formSurv, mfSurv2))
    if(attr(TermsSurv, 'intercept')) Ztime2 <- as.matrix(Ztime2[, - 1]) # excluding intercept #
  } else Ztime2 <- matrix(, ncol = 0, nrow = M)

  
  n <- nLong
  N <- length(Y)
  nu <- length(U)
  ncb <- ncol(B)
  
  GHQ <- gauss.quad(cntrlLst$nknot, kind = "hermite")
  b <- GHQ$nodes
  wGQ <- GHQ$weights
  
  Y.st <- split(Y, ID)
  B.st <- lapply(split(B, ID), function(x) matrix(x, ncol = ncb))
  Ztime22 <- if (ncz > 1) t(apply(Ztime2, 1, function(x) tcrossprod(x))) else Ztime2 ^ 2
  Btime22 <- if (ncb > 1) t(apply(Btime2, 1, function(x) tcrossprod(x))) else Btime2 ^ 2
  B2 <- if (ncb > 1) t(apply(B, 1, function(x) tcrossprod(x))) else B ^ 2
    
  tempResp <- strsplit(toString(formLongX), ", ")[[1]][c(2, 1)]
  tempResp <- paste(tempResp, collapse = "")
  tempForm <- strsplit(toString(splitFormula(formLongX)[[1]]), ", ")[[1]]
  tempForm <- tempForm[-1]
  tempForm[length(tempForm) + 1] <- "data = fitLME$data"
  tempForm <- paste(tempForm, collapse = ",")
  tempForm <- paste("lm(", tempResp, tempForm, ")", sep = "")
  fitLM <- eval(parse(text = tempForm))
  gamma <- as.vector(fitLM$coefficients)
  gamma.names <- names(fitLM$coefficients)
  gamma.names <- gsub("bs\\(.*\\)", "bs", gamma.names)
  
  surv.init <-  InitValMultGeneric(gamma = gamma,  B.st = B.st, n = n, Y.st = Y.st, ni = ni, model = model, ID = ID, Index = Index, start = start, stop = stop, B = B, Btime = Btime, 
                                   Btime2 = Btime2, event = event, Z = Z, ncz = ncz, Ztime2 = Ztime2, Index2 = Index2, Index1 = Index1, rho = rho, nk = nk, d = d, Ztime22 = Ztime22, 
                                   Ztime = Ztime, tol.P = cntrlLst$tol.P, iter = cntrlLst$max.iter)  
  phi <- surv.init$phi
  alpha <- surv.init$alpha
  lamb <- surv.init$lamb
  Ysigma <- surv.init$Ysigma
  Bsigma <- surv.init$Bsigma
  
  theta.old <- list(gamma = gamma, phi = phi,  alpha = alpha, Ysigma = Ysigma, Bsigma = Bsigma, 
                   lamb = lamb, lgLik = 0)
  err.P <- err.L <- step <- 1
  
  while (step <= cntrlLst$max.iter) {
    
    if (err.P < cntrlLst$tol.P | err.L < cntrlLst$tol.L) break
    
    theta.new <-  EMiterMultGeneric(theta.old, B.st, n, Y.st, b, model, Btime, Btime2, Index, Index0, Ztime, Ztime2, nknot = cntrlLst$nknot, nk, Index1, rho, d, wGQ, ID, ncb, B, Y, N, 
                                    ncz, Ztime22, Index2, B2, Btime22)
    
    new.P <- c(theta.new$gamma, theta.new$phi, theta.new$alpha, theta.new$Ysigma, theta.new$Bsigma)
    old.P <- c(theta.old$gamma, theta.old$phi, theta.old$alpha, theta.old$Ysigma, theta.old$Bsigma)
    err.P <- max(abs(new.P - old.P) / (abs(old.P) + .Machine$double.eps * 2))
    # add .Machine$double.eps *2 to avoid zero value of the estimated parameters #
    
    new.L <- theta.new$lgLik
    old.L <- theta.old$lgLik
    err.L <- abs(new.L - old.L) / (abs(old.L) + .Machine$double.eps * 2)

    step <- step + 1
    theta.old <- theta.new
  }
  converge <- as.numeric(err.P < cntrlLst$tol.P | err.L < cntrlLst$tol.L)
  

  if (cntrlLst$SE.method == 'PFDS') {
    time.SE <- system.time(Vcov <- PFDSMult(model, theta.new, delta = cntrlLst$delta, ncz = ncz, ncb = ncb, B.st = B.st, n =n, Y.st = Y.st, b = b, Btime = Btime, Btime2 = Btime2, 
                                            Index = Index, Ztime = Ztime, Ztime2 = Ztime2, Index0 = Index0, nk = nk, Index1 = Index1, rho = rho, d = d, wGQ = wGQ, Index2 = Index2, 
                                            alpha.name = alpha.name, phi.names = phi.names,N = N, Y = Y, B = B, ID = ID, nknot = cntrlLst$nknot, iter = cntrlLst$max.iter, 
                                            tol = min(cntrlLst$tol.P, cntrlLst$delta) / 100))[3]
    if (any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
      warning("NA's present in StdErr estimation due to numerical error!\n")
  } else if (cntrlLst$SE.method == 'PRES') {
    if (CheckDeltaMult(theta.new, cntrlLst$delta)) {
      time.SE <- system.time(Vcov <- PRESMult(model, theta.new, delta = cntrlLst$delta, ncz = ncz, ncb = ncb, B.st = B.st, n = n, Y.st = Y.st, b =b, Btime = Btime, Btime2 = Btime2, 
                                              Index = Index, Ztime = Ztime, Ztime2 = Ztime2, Index0 = Index0 , nk = nk, Index1 = Index1, rho = rho, d = d, wGQ = wGQ, Index2 =Index2, 
                                              alpha.name =alpha.name, phi.names = phi.names,N = N, Y = Y, B = B, ID = ID, nknot = cntrlLst$nknot, iter = cntrlLst$max.iter, 
                                              tol = min(cntrlLst$tol.P, cntrlLst$delta) / 100  ))[3]
      if (any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
        warning("NA's present in StdErr estimation due to numerical error!\n")
    } else {
      Vcov <- time.SE <- NA
      warning("\n 'delta' is too large, use smaller 'delta'!")
    }
  } else if (cntrlLst$SE.method == 'PLFD') {
    time.SE <- system.time(Vcov <- PLFDMult( model = model, theta.new, delta = cntrlLst$delta, B.st = B.st, n = n, Y.st = Y.st, b = b, Btime = Btime, Btime2 = Btime2, Index = Index, 
                                             Index0 = Index0, Ztime = Ztime, Ztime2 = Ztime2, nk = nk, Index1 = Index1, rho = rho, d = d, wGQ = wGQ, ncz = ncz, ncb = ncb, 
                                             Index2 = Index2, alpha.name = alpha.name, phi.names = phi.names, nknot = cntrlLst$nknot, iter = cntrlLst$max.iter, 
                                             tol = min(cntrlLst$tol.P, cntrlLst$delta) / 100))[3]
    if (any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
      warning("NA's present in StdErr estimation due to numerical error!\n")
  } else {
    Vcov <- time.SE <- NA
    warning("\n Standard error estimation method should be either 'PFDS', 'PRES' or 'PLFD'.")
  }
  
  theta.new$lamb <- data.frame("time" = U, "bashaz" = theta.new$lamb)
  names(theta.new$gamma) <- gamma.names
  names(theta.new$phi) <- phi.names
  names(theta.new$alpha) <- if (model == 1) alpha.name else "alpha"
  names(theta.new$Ysigma) <- "sigma.e"
  names(theta.new$Bsigma) <- "sigma.b"
  
  theta.new$est.bi <- matrix(theta.new$est.bi, ncol = 1)
  colnames(theta.new$est.bi) <- "bi"
  rownames(theta.new$est.bi) <- (fitLME$groups[[1]])[uniqueID]
  
  result <- list()
  result$coefficients <- theta.new
  result$logLik <- theta.new$lgLik
  result$call <- call
  result$numIter <- step
  result$Vcov <- Vcov
  result$est.bi <- theta.new$est.bi
  result$coefficients$est.bi <- NULL
  result$convergence <- if (converge == 1) "success" else "failure"
  result$control <- cntrlLst
  result$time.SE <- time.SE
  result$N <- N
  result$n <- n
  result$d <- d
  result$rho <- rho
  result$dataMat <- list(Y = Y, B = B, ID = ID, IDName = fitLME$groups[[1]])
  class(result) <-  unlist(strsplit(deparse(sys.call()), split = '\\('))[1]

  return(result)
}
