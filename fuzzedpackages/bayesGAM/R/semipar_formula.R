
#' \code{bayesGAM} fits a variety of regression models using Hamiltonian Monte Carlo
#'
#' Based on \code{\link[stats:glm]{glm}}.  \code{bayesGAM} is used to fit a variety of statistical models, including linear models, generalized lienar models, mixed effect models with random intercept, and semiparametric regression models.
#'
#' @export
#' @param formula a \code{\link[stats:formula]{formula}} object describing the model to be fitted.
#' @param random (optional) specify a random intercept in the form '~var'
#' @param family distribution and link function for the model
#' @param data (optional) data frame containing the variables in the model.
#' @param offset Same as \code{\link[stats:glm]{glm}}
#' @param beta (optional) list of priors for the fixed effects parameters.  Sensible priors are selected as a default.
#' @param eps (optional) list of priors for the error term in linear regression.  Sensible priors are selecteda as a default.
#' @param lambda (optional) list of priors for random effects variance parameters.  Sensible priors are selected as a default.
#' @param a (optional) list of priors for the off diagonal of the LDLT decomposed covariance matrix for multivariate response models.  Vague normal priors are used as a default.
#' @param spcontrol a list of control parameters for fitting the model in STAN.  See 'details'
#' @param store_plot_data a logical indicator for storing the plot data frame after simulation. Defaults to \code{FALSE}
#' @param method default currently set to 'bayesGAMfit'.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @details Similar to \code{glm}, models are typically specified by formula.
#' The formula typically takes the form \code{response ~ terms}, where the response is numeric
#' and terms specify the linear predictor for the response.  The terms may be numeric variables or factors.
#' @details The link function for the Generalized Linear Model is specified with a \code{\link[stats:family]{family}} object.
#' Currently, this package supports gaussian, binomial, and poisson families with all available link functions.
#' @details The list \code{spcontrol} currently supports additional parameters to facilitate fitting models.
#' \code{qr} is a logical indicator specifying whether the design matrix should be transformed via QR decomposition prior to HMC sampling.
#' QR decomposition often improves the efficiency with which HMC samples, as the MCMC chain navigates an orthogonal
#' space more easily than highly correlated parameters.  
#' \code{mvindep} is a logical indicator for multivariate response models with random intercepts. This indicates whether 
#' the multivariate responses should be considered independent. Defaults to \code{FALSE}
#' @references Hastie, T. J. (1992) Generalized additive models. Chapter 7 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie, Wadsworth & Brooks/Cole.
#' @references Dobson, A. J. (1990) \emph{An Introduction to Generalized Linear Models}. London: Chapman and Hall.
#' @return An object of class `bayesGAMfit`.  Includes slots:
#' @return
#' @return \code{results}:  `stanfit` object returned by `rstan::sampling`
#' @return \code{model}:  `glmModel` object 
#' @return \code{offset}:  offset vector from the input parameter
#' @return \code{spcontrol}:  list of control parameters from input
#' @examples
#' ## Dobson (1990) Page 93: Randomized Controlled Trial :
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' fpois<- bayesGAM(counts ~ outcome + treatment, family = poisson(),
#'                  spcontrol = list(qr = TRUE))
#' summary(fpois)
bayesGAM <- function (formula, random=NULL,
                         family = gaussian, data, offset, 
                         beta = list(),
                         eps = list(),
                         lambda = list(),
                         a = list(),
                         spcontrol = list(qr=TRUE, mvindep=FALSE, ...),
                         store_plot_data = FALSE, 
                         method = "bayesGAMfit", ...)
{

  # extract np terms
  getnp <- parse_formula(formula)
  
  Senv <- new.env(parent = environment(formula)) 
  assign("L", L, envir = Senv)
  
  if (missing(data)) {
    data <- environment(formula)
  }
  
  tt <- terms(getnp$sub_form, 
              data=data)
  
  has_intercept <- attr(terms(getnp$sub_form, 
                              data=data), "intercept") == 1
  
  # formula for everything but nl terms.  parse like glm
  call <- match.call()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  # from glm()
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  # family numbers
  # 1 = gaussian
  # 2 = binomial
  # 3 = poisson
  #

  famnum <- -1L
  linknum <- -1L
  if (family$family == "gaussian") {
    famnum <- 1L
  } else if (family$family == "binomial") {
    famnum <- 2L
  } else if (family$family == "poisson") {
    famnum <- 3L
  }

  # link numbers
  # 1 = identity
  # 2 = log
  # 3 = inverse
  # 4 = logit
  # 5 = probit
  # 6 = cauchit
  # 7 = cloglog
  # 8 = sqrt

  if (family$link == "identity") {
    linknum <- 1L
  } else if (family$link == "log") {
    linknum <- 2L
  } else if (family$link == "inverse") {
    linknum <- 3L
  } else if (family$link == "logit") {
    linknum <- 4L
  } else if (family$link == "probit") {
    linknum <- 5L
  } else if (family$link == "cauchit") {
    linknum <- 6L
  } else if (family$link == "log") {
    linknum <- 2L
  } else if (family$link == "cloglog") {
    linknum <- 7L
  } else if (family$link == "sqrt") {
    linknum <- 8L
  }


  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  # remove np() terms
  mf$formula <- getnp$sub_form
  
  mf2 <- mf
  mf2[[1L]] <- quote(function(...) { 
    stats::model.frame(na.action=na.pass, ...)
  })

  mf2 <- eval(mf2, Senv)
  
  not.na.obs <- complete.cases(mf2)
  
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, Senv)
  
  mt <- attr(mf, "terms")

  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  
  if (famnum == 2L) {
    if (is.factor(Y))
      Y <- as.integer(Y) -1
  }
  
  # add offset code
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y))
      stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }
  
  X <- if (!is.empty.model(mt)) {
    model.matrix(mt, mf)
  } else {
    matrix(, NROW(Y), 0L)
  }

  # check for random intercept
  e1 <- environment(mf$formula)
  assign("L", L, envir = e1)
  

  # create list for Z
  Z1 <- NULL
  Zlst <- list()
  random_intercept <- FALSE
  if (!is.null(random)) {
    # update to remove global intercept
    random <- update(random, ~ . +0)
    Z1 <- eval(model.matrix(random, data=data[not.na.obs, ]), envir = e1)
    Zlst <- c(Zlst, list(Z1))
    random_intercept <- TRUE
  }


  # get np terms
  npt <- getnp$npterms

  Xnp <- NULL
  Znp <- NULL
  allZ <- list()
  allknots <- list()
  allbasis <- character()
  alldegree <- integer()
  npargs <- list()

  if(length(npt) > 0) {

    allnp <- lapply(npt, function(xx) {
      ncall <- parse(text=xx)
      with(e1$data[not.na.obs, ], eval(ncall))
    })

    allX <- lapply(allnp, function(xx) {
      res <- xx$X
      if (xx$basis == "tps") {
        colnames(res) <- xx$Xnms
      } else if (xx$basis == "trunc.poly") {
        colnames(res) <- paste0(xx$Xnms, 1:xx$degree)
      }
      res
    })
    
    npargs <- lapply(allnp, function(xx) xx$Xnms)

    allZ <- lapply(allnp, function(xx) xx$Z)

    Xnp <- do.call(cbind, allX)

    allknots <- lapply(allnp, function(xx) xx$knots)
    allbasis <- sapply(allnp, function(xx) xx$basis)
    alldegree <- sapply(allnp, function(xx) xx$degree)
    names(allknots) <- npt
  }
  
  X <- cbind(X, Xnp)
  Zlst <- c(Zlst, allZ)

  if (length(Zlst) == 0) {
    q <- 0L
    zvars <- c(0L, 0L)
    Z <- matrix(, nrow=0, ncol=0)
  } else {
    q <- length(Zlst)
    zvars <- sapply(Zlst, ncol)
    Z <- do.call(cbind, Zlst)
  }
  
  mixed <- length(Z) > 0

  # name for dependent variable
  ynm <- getResponse(formula, data=data)

  # number of multresponse Y
  r <- ncol(as.matrix(Y))
  
  # remove NA
  if (mixed) {
    whichcomplete <- complete.cases(X, Z, Y)
    if (r > 1) {
      Y <- Y[whichcomplete, ]
    } else {
      Y <- Y[whichcomplete]
    }
    X <- X[whichcomplete, ,drop=FALSE]
    Z <- Z[whichcomplete, ]
    
  } else {
    whichcomplete <- complete.cases(X, Y)
    if (r > 1) {
      Y <- Y[whichcomplete, ]
    } else {
      Y <- Y[whichcomplete]
    }
    X <- X[whichcomplete, ,drop=FALSE]
  }
  
  # offset
  nobs <- NROW(Y)
  if (is.null(offset))
    offset <- rep.int(0, nobs)
  
  if (random_intercept) {
    Zint <- Z[, 1:zvars[1]]
  } else {
    Zint <- matrix(, nrow=0, ncol=0)
  }
  
  if (ncol(Z) > 0 & random_intercept) {
    Znp <- Z[, -c(1:zvars[1])]
  } else if (ncol(Z) > 0 & !random_intercept) {
    Znp <- Z
  } else {
    Znp <- matrix(, nrow=0, ncol=0)
  }
  
  spMCMCmodel <- new("glmModel",
                     famnum = famnum,
                     famname = family$family,
                     linknum = linknum,
                     linkname = family$link,
                     has_intercept = has_intercept,
                     p = ncol(X),
                     q = q,
                     r = r,
                     X = X,
                     Z = Z,
                     Zint = Zint, 
                     Znp = Znp, 
                     y = Y,
                     names_y = ynm,
                     zvars = zvars,
                     mixed = mixed,
                     knots = allknots,
                     basis = allbasis,
                     npdegree = as.integer(alldegree),
                     npargs = npargs,
                     npterms = npt,
                     sub_form = getnp$sub_form, 
                     random_intercept = random_intercept, 
                     call = call, 
                     offset = offset)

  # set priors if provided
  if (!all(length(beta)==0, length(eps)==0, length(lambda)==0, length(a)==0)) {

    spMCMCmodel <- setPrior(spMCMCmodel,
                            beta = beta,
                            eps = eps,
                            lambda = lambda,
                            a = a,
                            random_intercept = random_intercept)
  }

  # offset
  nobs <- NROW(Y)
  if (is.null(offset))
    offset <- rep.int(0, nobs)

  # call bayesGAM fit
  if (method == "bayesGAMfit") {

    stanResults <- bayesGAMfit(spMCMCmodel,
                                  offset = offset,
                                  spcontrol = spcontrol,
                                  ...)
  } else if (method == "predict") {
    return(spMCMCmodel)
  } else {
    stop("invalid method selected")
  }

  spResults <- new("bayesGAMfit",
                    results = stanResults,
                    model = spMCMCmodel,
                    offset = offset,
                    spcontrol = spcontrol)

  if (store_plot_data) {
    spResults@mcmcres <- as.matrix(spResults@results)
    spResults@pdata <- createPlotData(spResults, type="smooth")  
  }

  return(spResults)
}

