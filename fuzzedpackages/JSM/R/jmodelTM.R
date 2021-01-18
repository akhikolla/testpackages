# Joint Modeling Main Function with LME (linear mixed effects)

jmodelTM <- function (fitLME, fitCOX, data, model = 1, rho = 0, timeVarY = NULL, timeVarT = NULL, 
                      control = list(), ...)
{
  cat("Running jmodelTM(), may take some time to finish.\n")
  call <- match.call()

  CheckInputs(fitLME, fitCOX, rho)
  
  ID1 <- as.vector(unclass(fitLME$groups[[1]]))
  uniqueID <- !duplicated(ID1)
  tempID <- which(uniqueID)
  tempID <- c(tempID, length(ID1) + 1)
  ni <- diff(tempID) 
  ID <- rep(1:sum(uniqueID), times = ni)
  bBLUP1 <- data.matrix(ranef(fitLME))
  if (ncol(bBLUP1) == 1) {
    bBLUP <- matrix(bBLUP1[ID1[uniqueID], ], ncol = 1)
  } else {
    bBLUP <- bBLUP1[ID1[uniqueID], ]
    dimnames(bBLUP) <- NULL
  }
  nLong <- nrow(bBLUP)
  if (ncol(fitCOX$y) != 3)
    stop("\n must fit time-dependent Cox model in coxph().")
  start <- as.vector(fitCOX$y[, 1])
  stop <- as.vector(fitCOX$y[, 2])
  event <- as.vector(fitCOX$y[, 3])
  Time <- stop[cumsum(ni)]
  d <- event[cumsum(ni)] 
  nSurv <- length(Time) 
  if (sum(d) < 5)
    warning("\n more than 5 events are required.")
  if (nLong != nSurv)
    stop("\n sample sizes in the longitudinal and event processes differ.")
  
  W <- as.matrix(fitCOX$x)
  ncw <- ncol(W)

  varNames <- list()
  varNames$phi.names <- colnames(W)
  formSurv <- formula(fitCOX)
  TermsSurv <- fitCOX$terms
  mfSurv <- model.frame(TermsSurv, data)[cumsum(ni), ]
  if (!is.null(timeVarT)) {
    if (!all(timeVarT %in% all.vars(TermsSurv)))
      stop("\n'timeVarT' does not correspond columns in the fixed-effect design matrix of 'fitCOX'.")
    mfSurv[timeVarT] <- Time
  }
  if (ncw > 0) {
    Wtime <- as.matrix(model.matrix(formSurv, mfSurv))
    if(attr(TermsSurv, 'intercept')) Wtime <- as.matrix(Wtime[, - 1])
    # design matrix in survival part, one row for each subject, excluding intercept #
  } else Wtime <- matrix(, ncol = 0, nrow = nSurv)

  
  TermsLongX <- fitLME$terms
  mydata <- fitLME$data[all.vars(TermsLongX)] 
  formLongX <- formula(fitLME) 
  mfLongX <- model.frame(TermsLongX, data = mydata) 
  X <- as.matrix(model.matrix(formLongX, mfLongX))
  varNames$alpha.name <- rownames(attr(TermsLongX, "factors"))[attr(TermsLongX, "response")]
  
  formLongZ <- formula(fitLME$modelStruct$reStruct[[1]]) 
  mfLongZ <- model.frame(terms(formLongZ), data = mydata)
  TermsLongZ <- attr(mfLongZ, "terms") 
  Z <- as.matrix(model.matrix(formLongZ, mfLongZ)) 
  Y <- as.vector(model.response(mfLongX, "numeric"))
  # give the column in mfLongX which is considered as response, may be transformed #
  
  data.id <- mydata[uniqueID, ] # pick the first row of each subject in mydata, nrow=n #
  
  if (!is.null(timeVarY)) {
    if (!all(timeVarY %in% names(mydata)))
      stop("\n'timeVarY' does not correspond to columns in the fixed-effect design matrix of 'fitLME'.")
    data.id[timeVarY] <- Time
  }
  
  mfLongX.id <- model.frame(TermsLongX, data = data.id)
  Xtime <- as.matrix(model.matrix(formLongX, mfLongX.id)) # same structure with X, but with only n rows #
  
  mfLongZ.id <- model.frame(TermsLongZ, data = data.id)
  Ztime <- as.matrix(model.matrix(formLongZ, mfLongZ.id)) # same structure with Z, but with only n rows #
  
  U <- sort(unique(Time[d == 1])) # ordered uncensored observed event time #
  tempU <- lapply(Time, function(t) U[t >= U]) 
  times <- unlist(tempU) # vector of length M #
  nk <- sapply(tempU, length)  # length of each element in times, vector of length n #
  M <- sum(nk)
  
  Indcs <- list(); 
  
  Indcs$Index <- rep(1:nLong, nk) # repeat 1:n by nk, length M #
  Indcs$Index0 <- match(Time, U)
  Indcs$Index1 <- unlist(lapply(nk[nk != 0], seq, from = 1)) # vector of length M #
  Indcs$Index2 <- colSums(d * outer(Time, U, "==")) # vector of length nu #
  
  data.id2 <- data.id[Indcs$Index, ]
  if (!is.null(timeVarY)) {
    data.id2[timeVarY] <- times
  }
  mfLongX2 <- model.frame(TermsLongX, data = data.id2)
  Xtime2 <- as.matrix(model.matrix(formLongX, mfLongX2))
  mfLongZ2 <- model.frame(TermsLongZ, data = data.id2)
  Ztime2 <- as.matrix(model.matrix(formLongZ, mfLongZ2))
  
  mfSurv2 <- mfSurv[Indcs$Index, ]
  if (!is.null(timeVarT)) {
    mfSurv2[timeVarT] <- times
  }
  if (ncw > 0) {
    Wtime2 <- as.matrix(model.matrix(formSurv, mfSurv2))
    if(attr(TermsSurv, 'intercept')) Wtime2 <- as.matrix(Wtime2[, - 1]) # excluding intercept #
  } else Wtime2 <- matrix(, ncol = 0, nrow = M)

  
  n <- nLong
  N <- length(Y)
  nu <- length(U)
  ncz <- ncol(Z)
  ncx <- ncol(X)
  ncz2 <- ncz ^ 2
  p <- ncz * (ncz + 1) / 2
  
  cntrlLst <- GenerateControlList(control, ncz)  
  
  GHQ <- gauss.quad(cntrlLst$nknot, kind = "hermite")
  b <- as.matrix(expand.grid(rep(list(GHQ$nodes), ncz)))
  wGQ <- as.matrix(expand.grid(rep(list(GHQ$weights), ncz)))
  wGQ <- apply(wGQ, 1, prod)
  GQ <- nrow(b)
  
  Z.st <- lapply(split(Z, ID), function(x) matrix(x, ncol = ncz))
  Y.st <- split(Y, ID)
  X.st <- lapply(split(X, ID), function(x) matrix(x, ncol = ncx))
  Ztime2.st <- vector('list', n)
  for (i in (1:n)[nk != 0]) { Ztime2.st[[i]] <- matrix(Ztime2[Indcs$Index == i, ], ncol = ncz) }
  Wtime22 <- if(ncw > 1) t(apply(Wtime2, 1, function(x) tcrossprod(x))) else Wtime2 ^ 2
  Xtime22 <- if(ncx > 1) t(apply(Xtime2, 1, function(x) tcrossprod(x))) else Xtime2 ^ 2
  X2 <- if(ncx > 1) t(apply(X, 1, function(x) tcrossprod(x))) else X ^ 2
  X2.sum <- matrix(colSums(X2), nrow = ncx)  
      
  Bsigma <- lapply(lapply(fitLME$modelStruct$reStruct, as.matrix), 
                   function(x) x * fitLME$sigma ^ 2)[[1]]

  # the estimated variance-covariance matrix for the random effects  #
  beta <- as.vector(fixef(fitLME))
  varNames$beta.names <- names(fixef(fitLME))
  Ysigma <- fitLME$sigma
  
  surv.init <- InitValTMGeneric(beta, model = model, n = n, X = X, Z = Z, bBLUP = bBLUP, ID = ID, Xtime = Xtime, Ztime = Ztime, Xtime2 = Xtime2, Ztime2 = Ztime2, Indcs = Indcs, 
                                start = start, event = event, stop = stop, W = W, ncw = ncw, Wtime2 = Wtime2, rho= rho, nk = nk, Wtime22 = Wtime22, d = d,  Wtime = Wtime, cvals = cntrlLst)
  phi <- surv.init$phi
  alpha <- surv.init$alpha
  lamb <- surv.init$lamb
   
  theta.old <- list(beta = beta, phi = phi, alpha = alpha, Ysigma = Ysigma, Bsigma = Bsigma,  
                    lamb = lamb, lgLik = 0)
  err.P <- err.L <- step <- 1
  
  while (step <= cntrlLst$max.iter) {
     
    if (err.P < cntrlLst$tol.P | err.L < cntrlLst$tol.L) break
     
    theta.new <- EMiterTMGeneric(theta.old, n = n, Z.st = Z.st, Ztime = Ztime, Ztime2.st = Ztime2.st, nk = nk, Indcs = Indcs, Wtime2 = Wtime2, Xtime2 = Xtime2, GQ = GQ, 
                                 rho = rho, wGQ = wGQ, d = d, Y.st = Y.st, X.st = X.st, ncz = ncz, ncz2 = ncz2, b = b, model =  model, Wtime = Wtime, Xtime = Xtime, X = X, 
                                 Y = Y, ID = ID, N = N, ncw = ncw, Wtime22 = Wtime22, ncx = ncx, Xtime22 = Xtime22, Z = Z, X2.sum = X2.sum)
    new.P <- c(theta.new$beta, theta.new$phi, theta.new$alpha, theta.new$Ysigma, theta.new$Bsigma)
    old.P <- c(theta.old$beta, theta.old$phi, theta.old$alpha, theta.old$Ysigma, theta.old$Bsigma)
    err.P <- max(abs(new.P - old.P) / (abs(old.P) + .Machine$double.eps * 2))

    new.L <- theta.new$lgLik
    old.L <- theta.old$lgLik
    err.L <- abs(new.L - old.L) / (abs(old.L) + .Machine$double.eps * 2)
    step <- step + 1
    theta.old <- theta.new
  }
  converge <- as.numeric(err.P < cntrlLst$tol.P | err.L < cntrlLst$tol.L)
     
  if (cntrlLst$SE.method == 'PFDS') {
    if (CheckDeltaFD(theta.new, ncz, cntrlLst$delta)) {
      time.SE <- system.time(Vcov <- PFDS(model, theta.new, ncx = ncx, ncz = ncz, ncw = ncw, p = p, cvals = cntrlLst, varNames = varNames, Indcs = Indcs, n = n, Z.st = Z.st, 
                                          Y.st = Y.st, X.st = X.st, Ztime = Ztime, nk = nk, Wtime = Wtime, Wtime2 = Wtime2, Xtime = Xtime, Xtime2 = Xtime2, GQ = GQ, rho = rho, 
                                          d = d, wGQ = wGQ, ncz2 = ncz2, b = b, Ztime2.st = Ztime2.st, X = X, Y = Y, ID = ID, N = N, Z = Z))[3]
      if (any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
        warning("NA's present in StdErr estimation due to numerical error!\n")
    } else {
      Vcov <- time.SE <- NA
      warning("\n 'delta' is too large, use smaller 'delta'!")
    }
  } else if (cntrlLst$SE.method == 'PRES') {
    if (CheckDeltaRE(theta.new, ncz, cntrlLst$delta)) {
      time.SE <- system.time(Vcov <- PRES(model, theta.new, ncz = ncz, ncx = ncx, ncw = ncw, n = n, Z.st = Z.st, Y.st = Y.st, X.st = X.st, b = b, Ztime = Ztime, Ztime2.st = Ztime2.st, 
                                          nk = nk, Wtime = Wtime, Xtime = Xtime, Wtime2 = Wtime2, Xtime2 = Xtime2, rho = rho, Indcs = Indcs, wGQ =wGQ, GQ = GQ, d = d, p = p, ncz2 = ncz2, 
                                          X = X, Y = Y, Z = Z, ID = ID, N = N, cvals = cntrlLst, varNames = varNames))[3]
      if (any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
        warning("NA's present in StdErr estimation due to numerical error!\n")
    } else {
      Vcov <- time.SE <- NA
      warning("\n 'delta' is too large, use smaller 'delta'!")
    }
  } else if (cntrlLst$SE.method == 'PLFD') {
    if (CheckDeltaFD(theta.new, ncz, cntrlLst$delta)) {
      time.SE <- system.time(Vcov <- PLFD(model, theta.new, n= n, ncx = ncx, ncz = ncz, ncw = ncw, p = p, cvals = cntrlLst, varNames = varNames, Z.st = Z.st, Y.st =Y.st, X.st = X.st, 
                                          b = b, Ztime = Ztime, nk = nk, Indcs = Indcs, Wtime = Wtime, Xtime = Xtime, Wtime2 = Wtime2, Xtime2 = Xtime2, GQ = GQ, rho = rho, d = d, 
                                          wGQ = wGQ, Ztime2.st = Ztime2.st ))[3]
      if (any(is.na(suppressWarnings(sqrt(diag(Vcov))))))
        warning("NA's present in StdErr estimation due to numerical error!\n")
    } else {
      Vcov <- time.SE <- NA
      warning("\n 'delta' is too large, use smaller 'delta'!")
    }
  } else {
    Vcov <- time.SE <- NA
    warning("\n Standard error estimation method should be either 'PFDS', 'PRES' or 'PLFD'.")
  }

  theta.new$lamb <- data.frame("time" = U, "bashaz" = theta.new$lamb)
  names(theta.new$beta) <- varNames$beta.names
  names(theta.new$phi) <- varNames$phi.names
  names(theta.new$alpha) <- if(model == 1) varNames$alpha.name else "alpha"
  names(theta.new$Ysigma) <- "sigma.e"
  if (ncz > 1) dimnames(theta.new$Bsigma) <- dimnames(Bsigma) 
  else names(theta.new$Bsigma) <- "sigma.b"
  colnames(theta.new$est.bi) <- colnames(Bsigma)
  rownames(theta.new$est.bi) <- (fitLME$groups[[1]])[uniqueID]
  
  result <- list()
  result$coefficients <- theta.new
  result$logLik <- theta.new$lgLik
  result$call <- call
  result$numIter <- step
  result$Vcov <- Vcov
  result$est.bi <- theta.new$est.bi
  result$coefficients$est.bi <- NULL
  result$convergence <- if(converge == 1) "success" else "failure"
  result$control <- cntrlLst
  result$time.SE <- time.SE
  result$N <- N
  result$n <- n
  result$d <- d
  result$rho <- rho
  result$dataMat <- list(Y = Y, X = X, Z = Z, ID = ID, IDName = fitLME$groups[[1]])
  class(result) <-  unlist(strsplit(deparse(sys.call()), split = '\\('))[1]

  return(result)
}
