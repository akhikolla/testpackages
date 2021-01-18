
transfo <- function(X, type = "YJ", robust = TRUE, lambdarange = NULL,
                    prestandardize = TRUE, prescaleBC = F, scalefac = 1,
                    quant = 0.99, nbsteps = 2, checkPars = list()) {
  #
  # Function for fitting the Box-Cox (BC) or Yeo-Johnson (YJ) 
  # transformation to the columns of X, as described in:
  # 
  # J. Raymaekers and P.J. Rousseeuw, "Transforming variables to 
  #   central normality," https://arxiv.org/abs/2005.07946 .
  #
  # Arguments: 
  #   X              : a data matrix of dimensions n x d. Its columns
  #                    are the variables to be transformed.
  #   type           : one of the following:
  #                    "BC" : Box-Cox power transformation. Only works
  #                           for strictly positive variables. If this type
  #                           is given but a variable is not strictly
  #                           positive, the function stops with a
  #                           message about that variable.
  #                    "YJ" : Yeo-Johnson power transformation. The data
  #                           may have positive as well as negative values.
  #                    "bestObj" : for strictly positive variables both BC
  #                                and YJ are run, and the solution with
  #                                lowest objective is kept. On the other
  #                                variables YJ is run.
  #   robust         : if TRUE the proposed RewML (Reweighted Maximum 
  #                    Likelihood) method is used, which first
  #                    computes a robust initial estimate of the
  #                    transformation parameter lambda. If FALSE
  #                    the classical ML method is used.
  #   lambdarange    : range of lambda values that will be optimized over. 
  #                    If NULL, the range goes from -4 to 6.    
  #   prestandardize : whether to standardize the variables _before_ the
  #                    power transformation.
  #                    For BC the variable is divided by its median.
  #                    For YJ it subtracts a location estimate and 
  #                    divides by a scale estimate.
  #   prescaleBC     : for BC only. This standardizes the logarithm of
  #                    original variable by subtracting its median and
  #                    dividing by its mad, after which exp() turns
  #                    the result into a positive variable again.
  #   scalefac       : when YJ is fit and prestandardize = T, the 
  #                    standardized data is multiplied by scalefac. 
  #                    When BC is fit and prescaleBC = T the same happens 
  #                    to the standardized log of the original variable.
  #   quant          : quantile for determining the weights in the 
  #                    reweighting step (ignored when robust=F).
  #   nbsteps        : number of reweighting steps (ignored when 
  #                    robust=F).
  # 
  # Results:
  #   lambdahats  : the estimated lambda parameter for each column of x
  #   Xt          : each column is the transformed version of the
  #                 corresponding column of x
  #   muhat       : the estimated location of each column of Xt
  #   sigmahat    : the estimated scale of each column of Xt
  #   Zt          : Xt poststandardized by the centers in muhat and the
  #                 scales in sigmahat. Is always provided.
  #   weights     : the final weights from the reweighting.
  #   ttypes      : the type of transform used in each column.
  
  
  if (is.vector(X)) {X <- matrix(X, ncol = 1)}
  X <- as.matrix(X)
  
  # parameters for checkDataSet
  if (!"coreOnly" %in% names(checkPars)) {
    checkPars$coreOnly <- FALSE
  }
  if (!"silent" %in% names(checkPars)) {
    checkPars$silent <- FALSE
  }
  if (!"numDiscrete" %in% names(checkPars)) {
    checkPars$numDiscrete <- 5
  }
  if (!"precScale" %in% names(checkPars)) {
    checkPars$precScale <- 1e-12
  }
  
  # Start with checkDataSet to avoid transforming too discrete variables
  if (!checkPars$coreOnly) {
    # Check the data set and set aside columns and rows that do
    # not satisfy the conditions:
    out <- checkDataSet(X,
                        fracNA = 1,
                        numDiscrete = checkPars$numDiscrete,
                        precScale = checkPars$precScale,
                        silent = checkPars$silent)
    
    X <- out$remX 
  } 
  
  # The code of the main function starts here.
  #
  n <- nrow(X)
  d <- ncol(X)
  n ; d
  Xt         <- X
  lambdahats <- rep(1, d)
  ttypes     <- rep(NA, d)
  objective  <- rep(NA, d)
  weights    <- matrix(0, n, d)
  muhat      <- rep(0, d)
  sigmahat   <- rep(1, d)
  if (is.null(lambdarange)) { lambdarange <- c(-4, 6) }
  for (j in 1:d) { # loop over the columns of X
    x         <- X[, j]
    goodInds  <- which(!is.na(x))
    x         <- x[goodInds]
    ttype     <- getTtype(x, type, j)   
    ttypes[j] <- ttype
    if (robust == TRUE) { 
      order.x <- order(x)
      xsort   <- x[order.x] # the estimation code uses the sorted data
      # throughout, to avoid sorting many times.
      if (ttype == "BC" || ttype == "bestObj") {
        est.out.BC <- RewML_BC(xsort, lambdarange = lambdarange,
                               prestandardize = prestandardize, 
                               init = "BCr", quant = quant,
                               nbsteps = nbsteps,
                               prescaleBC = prescaleBC,
                               scalefac = scalefac)
        lambdahats[j]                 <- est.out.BC$lambdahat.rew
        objective[j]                  <- est.out.BC$critval.rew
        Xt[goodInds[order.x], j]      <- est.out.BC$yt.rew
        weights[goodInds[order.x], j] <- est.out.BC$weights
        muhat[j]                      <- est.out.BC$muhat
        sigmahat[j]                   <- est.out.BC$sigmahat
      }
      if (ttype == "YJ" || ttype == "bestObj") {
        est.out.YJ <- RewML_YJ(xsort, lambdarange = lambdarange,
                               prestandardize = prestandardize,
                               init = "YJr", quant = quant,
                               nbsteps = nbsteps,
                               scalefac = scalefac)
        lambdahats[j]                 <- est.out.YJ$lambdahat.rew
        objective[j]                  <- est.out.YJ$critval.rew
        Xt[goodInds[order.x], j]      <- est.out.YJ$yt.rew
        weights[goodInds[order.x], j] <- est.out.YJ$weights
        muhat[j]                      <- est.out.YJ$muhat
        sigmahat[j]                   <- est.out.YJ$sigmahat
      }
      if (ttype == "bestObj") {
        if (est.out.BC$critval.rew < est.out.YJ$critval.rew) {
          # RewML minimizes its objective
          lambdahats[j]                 <- est.out.BC$lambdahat.rew
          objective[j]                  <- est.out.BC$critval.rew
          Xt[goodInds[order.x], j]      <- est.out.BC$yt.rew
          weights[goodInds[order.x], j] <- est.out.BC$weights
          ttypes[j]                     <- "BC"
          muhat[j]                      <- est.out.BC$muhat
          sigmahat[j]                   <- est.out.BC$sigmahat
          ttypes[j] <- "BC"
        } else {
          ttypes[j] <- "YJ"
        }
      }
    } else {# robust = F so classical ML
      if (ttype == "BC" || ttype == "bestObj") {
        est.out.BC           <- estML(x, type = "BC", lambdarange = lambdarange,
                                      prestandardize = prestandardize,
                                      prescaleBC = prescaleBC,
                                      scalefac = scalefac)
        lambdahats[j]        <- est.out.BC$lambda
        objective[j]         <- est.out.BC$objective
        Xt[goodInds, j]      <- est.out.BC$xt
        weights[goodInds, j] <- est.out.BC$weights
      }
      if (ttype == "YJ" || ttype == "bestObj") {
        est.out.YJ           <- estML(x, type = "YJ", lambdarange = lambdarange,
                                      prestandardize = prestandardize,
                                      scalefac = scalefac)
        lambdahats[j]        <- est.out.YJ$lambda
        objective[j]         <- est.out.YJ$objective
        Xt[goodInds, j]      <- est.out.YJ$xt
        weights[goodInds, j] <- est.out.YJ$weights        
      }
      if (ttype == "bestObj") {
        if (est.out.BC$objective > est.out.YJ$objective) {
          # ML maximizes its objective
          lambdahats[j]        <- est.out.BC$lambda
          objective[j]         <- est.out.BC$objective
          Xt[goodInds, j]      <- est.out.BC$xt
          weights[goodInds, j] <- est.out.BC$weights
          ttypes[j] <- "BC"
        } else {
          ttypes[j] <- "YJ"
        }
      }
      muhat[j]    <- mean(Xt[,j], na.rm = TRUE)
      sigmahat[j] <- sd(Xt[,j], na.rm = TRUE)
    } # ends classical ML
  } # ends loop over columns of X
  Zt <- scale(Xt, center = muhat, scale = sigmahat)
  return(c(list(lambdahats = lambdahats, objective = objective,
              Xt = Xt, Zt = Zt, weights = weights, ttypes = ttypes, 
              muhat = muhat, sigmahat = sigmahat), out))
}


getTtype <- function(x, type, j=NULL) {
  # Gets the transformation type (BC, YJ, or bestObj) and checks 
  # whether the domain of the variable is compatible with that type.
  # args:
  #   x    : vector with the values of the variable
  #   type : one of "BC", "YJ", "bestObj"
  #   j    : number of the variable
  # 
  if (!type %in% c("BC", "YJ", "bestObj")) {
    stop(" type should be 'BC' or 'YJ' or 'bestObj'.")
  }
  minx <- min(x)
  if (type == "BC") {
    if (minx <= 0) {
      stop(paste("The minimum of variable ", j,
                 " is not strictly positive.\n",
                 "Please choose a type other than BC.",sep = ""))
    } else {
      ttype <- "BC"
    }
  } else if (type == "YJ") {
    ttype <- "YJ"
  } else if (type == "bestObj") {
    if (minx <= 0) {
      ttype <- "YJ"
    } else {
      ttype <- "bestObj"
    }
  }
  return(ttype)
}


RewML_BC <- function(x,
                     lambdarange = NULL,
                     prestandardize = TRUE,
                     prescaleBC = F,
                     scalefac = 1,
                     init = "BCr",
                     quant = 0.99,
                     nbsteps = 2) {
  # The reweighted ML estimator for the Box-Cox transformation parameter.
  #
  # args: 
  #   x:              vector of _sorted_ observations
  #   lambdarange:        grid of lambda values. If NULL, a grid between 
  #                   -2 and 4 is chosen
  #   prestandardize: whether to start by robustly standardizing the data
  #   init:           initial estimator. should be "BCr" or "BC"
  #   quant:          quantile for determining the weights in the 
  #                   reweighting step
  #   nbsteps:        number of reweighting steps
  #
  if (!init %in% c("BC", "BCr")) {
    stop("init should be either 'BC' or 'BCr'")
  }
  x <- na.omit(x)
  
  if (is.unsorted(x)) {
    x.order <- order(x, decreasing = FALSE)
    x       <- x[x.order]
  } else {
    x.order <- NULL
  }
  
  if (prestandardize) {
    x <- x / median(x) # so median is 1
    if (prescaleBC){ 
      if(is.null(scalefac)) { scalefac = 1 }
      x <- exp( scalefac*log(x)/mad(log(x)) ) # median stays 1
      if(min(x) < 1e-6) { x <- x + 1e-6 }
    }
  }
  
  # Range of lambda over which to optimize:
  if (is.null(lambdarange)) { lambdarange <- c(-4, 6) }
  if (init == "BCr") {
    tempfunc <- function(lambdatemp) {
      getCritval(x, lambdatemp, tfunc = BCr)
    }
    opt.out <- optimize(tempfunc, range(lambdarange))
    lambdahat.raw <- opt.out$minimum
    BCr.out.raw   <- BCr(x, lambdahat.raw)
    yt.raw        <- BCr.out.raw$yt
    zt.raw        <- BCr.out.raw$zt
  } else { 
    tempfunc <- function(lambdatemp) {
      getCritval(x, lambdatemp, tfunc = BC)
    }
    opt.out <- optimize(tempfunc, range(lambdarange))
    lambdahat.raw <- opt.out$minimum
    BC.out.raw    <- BC(x, lambdahat.raw)
    yt.raw        <- BC.out.raw$yt
    zt.raw        <- BC.out.raw$zt
  }
  rew.out <- reweightBCr(x = x, zt.raw = zt.raw,
                         lambdahat.raw = lambdahat.raw, 
                         lambdarange = lambdarange, quant = quant,
                         nbsteps = nbsteps,
                         prestandardize = prestandardize)
  BC.out.rew    <- rew.out$BC.out.rew
  yt.rew        <- BC.out.rew$yt
  zt.rew        <- BC.out.rew$zt
  critval.rew   <- rew.out$critval.rew
  lambdahat.rew <- rew.out$lambdahat.rew
  weights       <- rew.out$wgts
  muhat         <- mean(yt.rew[weights == 1])
  sigmahat      <- sd(yt.rew[weights == 1])
  
  if (!is.null(x.order)) {
    yt.raw[x.order] <- yt.raw
    zt.raw[x.order] <- zt.raw
    yt.rew[x.order] <- yt.rew
    zt.rew[x.order] <- zt.rew
    weights[x.order] <- weights
  }
  return(list(lambdahat.raw = lambdahat.raw,
              yt.raw = yt.raw,
              zt.raw = zt.raw,
              weights = weights,
              lambdahat.rew = lambdahat.rew,
              yt.rew = yt.rew,
              zt.rew = zt.rew,
              critval.rew = critval.rew, 
              muhat = muhat,
              sigmahat = sigmahat))
}


RewML_YJ <- function(x, 
                     lambdarange = NULL, 
                     quant = 0.99, 
                     init = "YJr",
                     prestandardize = TRUE, 
                     scalefac = 1,
                     nbsteps = 2) {
  # The reweighted ML estimator for the Yeo-Johnson transformation parameter.
  #
  # args: 
  #   x             : vector of _sorted_ observations
  #   lambdarange       : grid of lambda values. If NULL, a grid between 
  #                   -2 and 4 is chosen
  #   prestandardize: whether to start by robustly standardizing the data
  #   init          : initial estimator. should be "YJr" or "YJ"
  #   quant         : quantile for determining the weights in the
  #                   reweighting step
  #   nbsteps       : number of reweighting steps
  #
  if (!init %in% c("YJ", "YJr")) {
    stop("init should be either 'YJ' or 'YJr'")
  }
  
  x <- na.omit(x)
  
  if(is.unsorted(x)) {
    x.order <- order(x, decreasing = FALSE)
    x       <- x[x.order]
  } else {
    x.order <- NULL
  }
  
  if (prestandardize) {
    if(is.null(scalefac)) { scalefac = 1 }
    x <- scalefac*(x - median(x))/mad(x)
  }
  
  # Range of lambda over which to optimize:
  if (is.null(lambdarange)) { lambdarange <- c(-4,6) }
  if (init == "YJr") {
    tempfunc <- function(lambdatemp) {
      getCritval(x, lambdatemp, tfunc = YJr)
    }
    opt.out <- optimize(tempfunc, range(lambdarange))
    lambdahat.raw <- opt.out$minimum
    YJr.out.raw   <- YJr(x, lambdahat.raw)
    yt.raw        <- YJr.out.raw$yt
    zt.raw        <- YJr.out.raw$zt
  } else { 
    tempfunc <- function(lambdatemp) {
      getCritval(x, lambdatemp, tfunc = YJ)
    }
    opt.out <- optimize(tempfunc, range(lambdarange))
    lambdahat.raw <- opt.out$minimum
    YJr.out.raw   <- YJ(x, lambdahat.raw)
    yt.raw        <- YJr.out.raw$yt
    zt.raw        <- YJr.out.raw$zt
  }
  # Reweighting:
  rew.out <- reweightYJr(x, zt.raw, lambdahat.raw, lambdarange, 
                         quant = quant, nbsteps = nbsteps)
  weights       <- rew.out$wgts
  YJ.out.rew    <- rew.out$YJ.out.rew
  yt.rew        <- YJ.out.rew$yt
  zt.rew        <- YJ.out.rew$zt
  critval.rew   <- rew.out$critval.rew
  lambdahat.rew <- rew.out$lambdahat.rew
  muhat         <- mean(yt.rew[weights == 1])
  sigmahat      <- sd(yt.rew[weights == 1])
  
  if (!is.null(x.order)) {
    yt.raw[x.order] <- yt.raw
    zt.raw[x.order] <- zt.raw
    yt.rew[x.order] <- yt.rew
    zt.rew[x.order] <- zt.rew
    weights[x.order] <- weights
  }
  
  return(list(lambdarange = lambdarange,
              lambdahat.raw = lambdahat.raw,
              yt.raw = yt.raw,
              zt.raw = zt.raw,
              weights = weights,
              critval.rew = critval.rew,
              lambdahat.rew = lambdahat.rew,
              yt.rew = yt.rew,
              zt.rew = zt.rew, 
              muhat = muhat,
              sigmahat = sigmahat))
}


reweightYJr <- function(x, 
                        zt.raw, 
                        lambdahat.raw,
                        lambdarange,
                        quant = 0.99, 
                        nbsteps = 2) {
  # Function for reweighted maximum likelihood, based on an initial estimate.
  # args: 
  #   x             : vector of sorted original observations
  #   zt.raw        : vector of sorted poststandardized transformed data 
  #                   (with initial lambdahat)
  #   lambdahat.raw : initial estimate for transformation parameter lambda
  #   lambdarange   : range of lambda values.
  #   quant         : quantile for determining the weights in the 
  #                   reweighting step
  #   nbsteps       : number of reweighting steps
  #
  if (is.null(lambdarange)) { lambdarange = c(-4,6) }
  # initial reweighting:
  wgts    <- abs(zt.raw) <= sqrt(qchisq(quant,1))
  rewinds <- which(wgts == 1)
  x.rew   <- x[rewinds]
  x.all   <- x
  # YJ transform with weighted ML:
  for (k in 1:nbsteps) {
    lambdahat.rew <- estML(x = x.rew, lambdarange = lambdarange, 
                           type = "YJ", prestandardize = F)$lambda
    YJ.out.rew    <- YJ(x.all, lambdahat.rew)
    wgts    <- abs(YJ.out.rew$zt) <= sqrt(qchisq(quant,1))
    rewinds <- which(wgts == 1)
    x.rew   <- x[rewinds]
  }
  critval.rew <- getCritval(x.all, lambdahat.rew, tfunc = YJ,
                            quant=quant)
  return(list(wgts = wgts,
              critval.rew  = critval.rew,
              lambdahat.rew = lambdahat.rew,
              YJ.out.rew = YJ.out.rew))
}


getChangepointYJr <- function(y, lambda, fac = 1.5, eps = 1e-5) {
  # Function to caculate the "changepoint" for the rectified
  # Yeo-Johnson transformation
  # 
  quarts = localFivenum(y)[2:4]
  if (lambda < 1){
    chg = YJ(quarts[3], lambda, stdToo = FALSE)$yt * 1.5
  } else if (lambda > 1){
    chg = YJ(quarts[1], lambda, stdToo = FALSE)$yt * 1.5
  }
  if (lambda < 0) {
    chg <- min(chg, abs(1 / lambda) - eps)
  } else if (lambda > 2) {
    chg <- max(chg, (1 / (2 - lambda)) + eps)
  }
  chg <- iYJ(chg, lambda, stdToo = FALSE)$yt
  chg <- min(max(chg, y[1]), y[length(y)])
  return(chg)
}


getChangepointBCr <- function(y, lambda, fac = 1.5, eps = 1e-5) {
  # Function to caculate the "changepoint" for the rectified
  # Box Cox transformation
  # 
  quarts = localFivenum(y)[2:4]
  if (lambda < 1){
    chg = BC(quarts[3], lambda, stdToo = FALSE)$yt * 1.5
  } else if (lambda > 1){
    chg = BC(quarts[1], lambda, stdToo = FALSE)$yt * 1.5
  }
  #correct for image of BC
  if (lambda < 0) {
    chg <- min(chg, abs(1 / lambda) - eps)
  } else  if (lambda > 0) {
    chg <- max(chg, -(1 / lambda) + eps)
  }
  chg <- iBC(chg, lambda, stdToo = FALSE)$yt
  chg <- min(max(chg, y[1]), y[length(y)])
  return(chg)
}


reweightBCr <- function(x, zt.raw, lambdahat.raw, lambdarange, quant = 0.99,
                        nbsteps = 2, prestandardize = TRUE) {
  # function for reweighted maximum likelihood, based on an initial estimate.
  # args: 
  #   x              : vector of sorted original observations
  #   zt.raw         : vector of sorted poststandardized transformed data 
  #                    (from the initial lambdahat)
  #   lambdahat.raw  : initial estimate for transformation parameter lambda
  #   lambdarange    : range of lambda values.
  #   prestandardize : whether to start by robustly standardizing the data
  #   quant          : quantile for determining the weights in the 
  #                    reweighting step
  #   nbsteps        : number of reweighting steps
  #
  if (is.null(lambdarange)) { lambdarange = c(-4,6) }
  # Initial weights:
  weights <- abs(zt.raw) <= sqrt(qchisq(quant,1))
  x.all   <- x
  x.rew   <- x[which(weights == 1)]
  # 
  # BC transform with weighted ML
  #
  for (k in 1:nbsteps) {
    lambdahat.rew <- estML(x = x.rew, lambdarange = lambdarange, 
                           type = "BC", prestandardize = F)$lambda
    BC.out.rew <- BC(x.all, lambdahat.rew)
    wgts       <- abs(BC.out.rew$zt) <= sqrt(qchisq(quant,1))
    rewinds    <- which(wgts == 1)
    x.rew      <- x[rewinds]
  }
  critval.rew <- getCritval(x.all, lambdahat.rew, tfunc = BC,quant = quant)
  return(list(BC.out.rew = BC.out.rew,
              lambdahat.rew = lambdahat.rew, 
              critval.rew = critval.rew,
              wgts = wgts))
}



estML <- function(x, lambdarange = NULL, type = "BC",
                  prestandardize = TRUE,
                  prescaleBC = F, scalefac = 1) {
  # Computes the ML estimator of the parameter lambda in the BC 
  # and YJ transformation, for a vector x (corresponding to a
  # column of the input data X).
  # Assumes x has no NA's, as they were already taken out in the main. 
  #
  n <- length(x)
  if (is.null(lambdarange)) { lambdarange <- c(-4, 6) }
  if (type == "YJ") { # Yeo-Johnson
    if (prestandardize) {
      if(is.null(scalefac)) { scalefac = 1 }
      x <- scalefac*scale(x)
    }
    tempfunc <- function(lambdatemp) {
      xyj    <- YJ(x, lambdatemp, stdToo = FALSE)$yt
      mu     <- mean(xyj)
      sigma2 <- mean((xyj - mu)^2)
      result <- -0.5 * n * log(2 * pi) -
        0.5 * n * log(sigma2) -
        0.5/sigma2 * sum( (xyj - mu)^2) +
        (lambdatemp - 1) * sum(sign(x) * log(1 + abs(x)))
      return(result)
    }
    opt.out <- optimize(tempfunc, range(lambdarange), maximum = TRUE)
    objective = tempfunc(opt.out$maximum) 
    xt <- YJ(x, opt.out$maximum, stdToo = FALSE)$yt
    zt <- scale(xt)
  } else {# Box-Cox
    if (prestandardize) {
      x <- x/median(x)
      if (prescaleBC){ 
        if(is.null(scalefac)) { scalefac = 1 }
        x <- exp( scalefac*log(x)/sd(log(x)) ) # median stays 1
        if(min(x) < 1e-6) { x <- x + 1e-6 }
      }
    }
    tempfunc <- function(lambdatemp) {
      xyj    <- BC(x, lambdatemp, stdToo = FALSE)$yt
      mu     <- mean(xyj)
      sigma2 <- mean((xyj - mu)^2)
      result <- -n/2 * log(sigma2) + (lambdatemp - 1) * sum(log(x))
      return(result)
    }
    opt.out <- optimize(tempfunc, range(lambdarange), maximum = TRUE)
    objective = tempfunc(opt.out$maximum) 
    xt <- BC(x, opt.out$maximum, stdToo = FALSE)$yt
    zt <- scale(xt)
  }
  return(list(lambda = opt.out$maximum,
              lambdarange = lambdarange, objective = objective,
              xt = xt, zt = zt, weights = rep(1, n)))
}


getCritval <- function(x, lambda, tfunc = YJr, chg2=NULL,  quant=0.99) {
  # Calculates the value of the robust fitting criterion.
  # Arguments:
  # x      : the univariate data
  # lambda : the value of lambda
  # tfunc  : either BC, BCr, YJ, YJr
  # chg2   : 2 changepoints: one for lambda >= 1 and one if lambda < 1
  #          so chg[1] < chg[2]
  #          If NULL, uses Q1 for lambda > 1 and Q3 for lambda > 1
  #
  if (is.null(chg2)) { chg = NULL } else {
    chg <- if (lambda >= 1) {chg2[1]} else {chg2[2]}
  }
  robnormality(tfunc(x,lambda,chg,stdToo = FALSE)$yt, b = 0.5)
}


BC <- function(y, lambda, chg = NULL, stdToo = TRUE) {
  # The Box-Cox transformation.
  #
  if (lambda == 0) {
    yt <- log(y)
  } else {
    yt <- (y^lambda - 1) / lambda
  }
  if (stdToo) {
    if (length(y) > 1) {
      locScale <- cellWise::estLocScale(matrix(yt, ncol = 1), 
                                        type = "hubhub")
      zt <- (yt - locScale$loc) / locScale$scale
    } else { 
      zt <- yt
    }
  } else {
    zt <- NULL
  }
  return(list(yt = yt, zt =zt))
}


iBC <- function(y, lambda, chg = NULL, stdToo = TRUE) {
  # Inverse of the Box-Cox transformation.
  #
  if (lambda == 0) {
    yt <- exp(y)
  } else {
    yt <- (1 + y * lambda)^(1/lambda)
  }
  
  if (stdToo) {
    zt = (yt - median(yt)) / mad(yt)
  } else {
    zt = NULL
  }
  return(list(yt = yt, zt = zt))
}


YJ <- function(y, lambda, chg = NULL, stdToo = TRUE) {
  # The Yeo-Johnson transformation.
  #
  indlow  <- which(y < 0)
  indhigh <- which(y >= 0)
  if (lambda != 0) {
    y[indhigh] = ((1 + y[indhigh])^(lambda) - 1) / lambda
  } else { 
    y[indhigh] = log(1 + y[indhigh])
  }
  if (lambda != 2) {
    y[indlow] = -((1 - y[indlow])^(2 - lambda) - 1) / (2 - lambda)
  } else {
    y[indlow] = -log(1 - y[indlow])
  }
  if (stdToo) {
    if (length(y) > 1) {
      locScale <- cellWise::estLocScale(matrix(y, ncol = 1), 
                                        type = "hubhub")
      zt <- (y - locScale$loc) / locScale$scale
    } else {
      zt <- y
    }
  } else {
    zt <- NULL
  }
  return(list(yt = y, zt = zt))
}


iYJ <- function(y, lambda, stdToo = TRUE) {
  # Inverse of the Yeo-Johnson transformation.
  #
  indlow  <- which(y < 0)
  indhigh <- which(y >= 0)
  if (lambda != 0) {
    y[indhigh] = (1 + lambda * y[indhigh])^(1 / lambda) - 1
  } else { 
    y[indhigh] = exp(y[indhigh]) - 1
  }
  if (lambda != 2) {
    y[indlow] = -((1 + (lambda - 2) * y[indlow])^(1/(2-lambda))) + 1
  } else {
    y[indlow] = 1 - exp(-y[indlow])
  }
  if (stdToo) {
    zt = (y - median(y)) / mad(y)
  } else {
    zt = NULL
  }
  return(list(yt = y, zt = zt))
}


YJprimex <- function(x, lambda) {
  # partial derivative of YJ with respect to the argument x.
  (1 + abs(x))^(sign(x) * (lambda - 1))
}


estMTL <- function(x, alpha = 0.95,
                   lambdagrid = NULL, type = "YJ",
                   prestandardize = TRUE) {
  # The Maximum Trimmed Likelihood estimator.
  # 
  if (is.null(lambdagrid)) {
    lambdagrid <- ((0:200) - 100)/50
  }
  result  <- rep(0, length(lambdagrid))
  n       <- length(x)
  h       <- ceiling(alpha * n)
  if (prestandardize) {
    if (type == "YJ") {
      x <- (x - median(x))/mad(x)
    }
    if (type == "BC") {
      x <- x / median(x)
    }
  }
  x_order <- order(x) # ranks of the data
  if (type == "YJ") { # Yeo-Johnson
    for (i in 1:length(lambdagrid)) { # loop over lambdagrid
      lambda <- lambdagrid[i]
      xyj <- YJ(x, lambda)$yt
      nbSubsets <- n - h + 1
      sds   <- rep(1, nbSubsets)
      means <- rep(0, nbSubsets)
      for (j in 1:nbSubsets) { # loop over subsets
        datatemp <- xyj[x_order[j:(j + h - 1)]]
        # subset of ordered xyj
        means[j] <- mean(datatemp)
        sds[j]   <- sd(datatemp)
      } # ends loop over subsets
      bestSubset <- which.min(sds) 
      # this gives the LTS=MCD subset
      idx   <- x_order[bestSubset:(bestSubset + h - 1)] 
      # ranks of LTS subset
      mu    <- means[bestSubset] # LTS=MCD location
      sigma <- sds[bestSubset]   # LTS=MCD scale
      result[i] <- -0.5 * h * log(2 * pi) -
        0.5 * h * log(sigma^2) -
        0.5  / sigma^2 * sum( (xyj[idx] - mu)^2) + 
        (lambda - 1) * sum(sign(x[idx])*log(1 + abs(x[idx])))
    } # ends loop over lambdagrid
    xt <- YJ(x, lambda = lambdagrid[which.max(result)])$yt
    zt <- (xt - median(xt)) / mad(xt)
  } else { # Box-Cox
    for (i in 1:length(lambdagrid)) {
      lambda  <- lambdagrid[i]
      xyj     <- BC(x, lambda)$yt
      nbSubsets <- n - h + 1
      sds   <- rep(1, nbSubsets)
      means <- rep(0, nbSubsets)
      for (j in 1:nbSubsets) { # loop over subsets
        datatemp <- xyj[x_order[j:(j + h - 1)]]
        # subset of ordered xyj
        means[j] <- mean(datatemp)
        sds[j]   <- sd(datatemp)
      } # ends loop over subsets
      bestSubset <- which.min(sds) 
      # this gives the LTS=MCD subset
      idx   <- x_order[bestSubset:(bestSubset + h - 1)] 
      # ranks of LTS subset
      mu    <- means[bestSubset] # LTS=MCD location
      sigma <- sds[bestSubset]   # LTS=MCD scale
      result[i] <- -n / 2 * log(sigma^2) + (lambda - 1) * sum(log(x))
    }
    xt <- BC(x, lambdagrid[which.max(result)])$yt
    zt <- (xt - median(xt)) / mad(xt)
  }
  return(list(lambda = lambdagrid[which.max(result)],
              lambdagrid = lambdagrid, LL = result,
              xt = xt, zt = zt))
}


localFivenum = function(x, na.rm = TRUE) {
  # Our version of stats::fivenum where now the quartiles
  # are order statistics, with symmetric ranks.
  # Assumes x is _sorted_ already.
  xna <- is.na(x)
  if (any(xna)) {
    if (na.rm) 
      x <- x[!xna]
    else return(rep.int(NA, 5))
  }
  n <- length(x)
  if (n == 0) 
    rep.int(NA, 5)
  else {
    med <- if ((n %% 2) == 0) {
      (x[(n/2)] + x[(n/2) + 1]) / 2
    } else {
      x[(n - 1)/2 + 1]
    }
    result <- c(x[1], x[ceiling(n / 4)], med,
                x[n - ceiling(n / 4) + 1], x[n])
    result
  }
}


BCr = function(y, lambda, chg = NULL, stdToo = TRUE) {
  # Rectified Box Cox transformation.
  # Like the Classical Box Cox transformation but with a linear 
  # tail on the shrinking side, and a continuous derivative.
  #
  # Arguments:
  # y      : univariate data
  # lambda : transformation parameter for BC
  # chg    : change point: where the linear part begins
  #          when lambda < 1, or ends when lambda > 1.
  #          If NULL, uses the getChangepoint function
  # stdToo : also output poststandardized transformed data.
  #
  if(min(y) <= 0) stop("Data values should be strictly positive")
  yt = rep(NA,times = length(y))
  chgt <- NULL
  if(lambda == 1){ # is already linear
    yt = y-1
    chgt = 0
  }
  if (lambda > 1){
    if(is.null(chg)) { chg = getChangepointBCr(y, lambda = lambda) }
    indl = which(y < chg)
    if (length(indl) > 0) {
      yt[-indl] = (y[-indl]^lambda-1)/lambda
      chgt      = (chg^lambda-1)/lambda
      yt[indl]  = chgt + (y[indl]-chg)*chg^(lambda-1)
    } else {
      yt = (y^lambda-1)/lambda
    }
  } else if(lambda < 1){
    if(is.null(chg)) { chg = getChangepointBCr(y,  lambda = lambda) }
    indu = which(y > chg)    
    if (lambda == 0) {
      if (length(indu) > 0) {
        yt[-indu] = log(y[-indu])
        chgt      = log(chg)
        yt[indu]  = chgt + (y[indu]-chg)/chg
      } else {
        yt = log(y)
      }
    } else {
      if (length(indu) > 0) {
        yt[-indu] = (y[-indu]^lambda-1)/lambda
        chgt      = (chg^lambda-1)/lambda
        yt[indu]  = chgt + (y[indu]-chg)*chg^(lambda-1)
      } else {
        yt = (y^lambda-1)/lambda
      }
    }
  }
  if (stdToo) {
    if (length(y) > 1) {
      locScale <- cellWise::estLocScale(matrix(yt, ncol = 1), 
                                        type = "hubhub")
      zt <- (yt - locScale$loc)/locScale$scale
    } else {
      zt <- yt
    }
  } else {
    zt = NULL
  }
  return(list(yt=yt, chg=chg, chgt=chgt, zt=zt))
}


YJr = function(y, lambda, chg = NULL, prec = 1e-10, stdToo = TRUE) {
  # Rectified Yeo-Johnson transformation.
  # Like the classical YJ transformation but with a linear tail 
  # on the shrinking side, and a continuous derivative.
  #
  # Arguments:
  # y      : univariate data. Assumed shifted to median zero.
  # lambda : transformation parameter for YJ.
  # chg    : change point: where the linear part begins
  #          when lambda < 1, or ends when lambda > 1.
  #          If NULL, uses Q3 for lambda > 1 and
  #          Q1 for lambda > 1.
  # stdToo : also output poststandardized transformed data.
  #
  quarts = localFivenum(y)[2:4]
  yt = rep(NA,times = length(y))
  if (lambda == 1){ # is already linear
    yt = y
    chgt = 0
  }
  indneg = which(y < 0)
  indpos = which(y >= 0)
  #
  if(lambda > 1){
    if(is.null(chg)) { chg = getChangepointYJr(y, lambda = lambda) }
    indl = which(y < chg) # to be linearized
    indb = which(y < 0 & y >= chg) # in between
    yt[indpos] = ((1 + y[indpos])^(lambda) - 1)/lambda
    #
    if(lambda == 2){
      yt[indb] = -log(1-y[indb])
      chgt     = -log(1-chg)
      yt[indl] = chgt + (y[indl]-chg) * YJprimex(chg,lambda)
    }
    else {
      yt[indb] = -((1-y[indb])^(2-lambda)-1)/(2-lambda)
      chgt     = -((1-chg)^(2-lambda)-1)/(2-lambda)
      yt[indl] = chgt + (y[indl]-chg) * YJprimex(chg,lambda)
    }
  }  
  if(lambda < 1){
    if(is.null(chg)) { chg = getChangepointYJr(y, lambda = lambda) }
    # if(chg <= 0) stop("chg should be positive")
    indu = which(y > chg) # to be linearized
    indb = which(y >= 0 & y <= chg) # in between
    yt[indneg] = -((1-y[indneg])^(2-lambda)-1)/(2-lambda)
    #
    if(lambda == 0){
      yt[indb] = log(1 + y[indb])
      chgt     = log(1 + chg)
      yt[indu] = chgt + (y[indu]-chg)*YJprimex(chg,lambda)
    }
    else {
      yt[indb] = ((1+y[indb])^(lambda)-1)/lambda
      chgt     = ((1+chg)^(lambda)-1)/lambda
      yt[indu] = chgt + (y[indu]-chg)*YJprimex(chg,lambda)
    }
  } 
  if (stdToo) {
    if (length(y) > 1) {
      locScale <- cellWise::estLocScale(matrix(yt, ncol = 1), 
                                        type = "hubhub")
      zt <- (yt - locScale$loc) / locScale$scale
    } else { 
      zt <- yt
    } 
  } else {
    zt <- NULL
  }
  return(list(yt=yt, chg=chg, chgt=chgt, zt=zt))
}


rhoBiweight <- function(x, c = 0.5) {
  xbw <- c * (1 - (1 - (x/c)^2)^3)
  xbw[which(abs(x) > c)] <- c
  return(xbw)
}


robnormality <- function(y, b = 0.5) {
  # Assess normality of a dataset after robustly standardizing.
  # Outliers only have a limited effect on this criterion.
  # y is assumed to be _sorted_.
  #
  # args: 
  #   y : vector of observations
  #   b : tuning parameter of biweight rho
  # 
  sy       <- y[!is.na(y)]
  n        <- length(sy)
  locScale <- cellWise::estLocScale(matrix(sy, ncol = 1), 
                                    type = "hubhub")
  sy       <- (sy - locScale$loc) / locScale$scale
  # Theoretical quantiles and differences:
  hprobs     <- ((1:n) - 1/3)/(n + 1/3)
  theoQs     <- qnorm(hprobs)
  diffs      <- sy - theoQs
  # apply Biweight rho
  rhodiffs <- rhoBiweight(diffs, b)
  crit     <- mean(abs(rhodiffs))
  return(crit)
}