#' fitClere function
#' 
#' This function runs the CLERE Model. It returns an object of class
#' \code{\linkS4class{Clere}}. For more details please refer to
#' \code{\link{clere}}.
#' 
#' 
#' @param y [numeric]: The vector of observed responses - size \code{n}.
#' @param x [matrix]: The matrix of predictors - size \code{n} rows and
#' \code{p} columns.
#' @param g [integer]: Either the number or the maximum of groups for fitting
#' CLERE. Maximum number of groups is considered when model selection is
#' required.
#' @param nItMC [numeric]: Number of Gibbs iterations to generate the
#' partitions. After the \code{nBurn} iterations, this number is automatically
#' set to \code{1}.
#' @param nItEM [numeric]: Number of SEM iterations.
#' @param nBurn [numeric]: Number of SEM iterations discarded before
#' calculating the MLE which is averaged over SEM draws.
#' @param dp [numeric]: Number of iterations between sampled partitions when
#' calculating the likelihood at the end of the run.
#' @param nsamp [numeric]: Number of sampled partitions for calculating the
#' likelihood at the end of the run.
#' @param maxit [numeric]: An EM algorithm is used inside the SEM to maximize
#' the complete log-likelihood p(y, Z|theta). \code{maxit} stands as the
#' maximum number of EM iterations for the internal EM.
#' @param tol [numeric]: Maximum increased in complete log-likelihood for the
#' internal EM (stopping criterion).
#' @param nstart [integer]: Number of random starting points to be used for
#' fitting the model.
#' @param parallel [logical]: Should the estimation from \code{nstart} random
#' starting points run in parallel?
#' @param seed [integer]: An integer given as a seed for random number
#' generation. If set to \code{NULL}, then a random seed is generated between
#' \code{1} and \code{1000}.
#' @param plotit [logical]: Should a summary plot (base plot) be drawn after
#' the run?
#' @param sparse [logical]: Should a \code{0} class be imposed to the model?
#' @param analysis [character]: Which analysis is to be performed. Values are
#' \code{"fit"}, \code{"bic"}, \code{"aic"} and \code{"icl"}.
#' @param algorithm [character]: The algorithm to be chosen to fit the model.
#' Either the SEM-Gibbs algorithm or the MCEM algorithm. The most efficient
#' algorithm being the SEM-Gibbs approach. MCEM is not available for binary
#' response.
#' @param theta0 [vector(numeric)]: An initial guess of the model parameters.
#' When considering g components, the length of \code{theta0} must be
#' \code{2*g+3} and \code{theta0} should be filled as intercept, the b_k's (g
#' real numbers), the pi_k's (g real numbers summing to 1), sigma^2 and gamma^2
#' (two positive numbers).
#' @param Z0 [vector(integer)]: A vector of integers representing an initial
#' partition for the variables. For 10 variables and 3 groups \code{Z0} can be
#' defined as \ \code{Z0 = c(rep(0, 2), rep(1, 3), rep(2, 5))}.
#' 
#' @return Object of class \code{\linkS4class{Clere}.}
#' 
#' @export
#' 
#' @seealso Overview : \code{\link{clere-package}} \cr 
#' Classes : \code{\linkS4class{Clere}}, \code{\linkS4class{Pacs}} \cr 
#' Methods : \code{\link{plot}}, \code{\link{clusters}}, \code{\link{predict}}, \code{\link{summary}} \cr 
#' Functions : \code{\link{fitClere}}, \code{\link{fitPacs}} 
#' Datasets : \code{\link{numExpRealData}}, \code{\link{numExpSimData}}, \code{\link{algoComp}}
#' 
#' @examples
#' 
#' library(clere)
#' plotit    <- FALSE
#' sparse    <- FALSE
#' nItEM     <- 100
#' nBurn     <- nItEM / 2
#' nsamp     <- 100
#' analysis  <- "fit"
#' algorithm <- "SEM"
#' nItMC     <- 1
#' dp        <- 2
#' maxit     <- 200
#' tol       <- 1e-3
#' 
#' n         <- 50
#' p         <- 50
#' intercept <- 0
#' sigma     <- 10
#' gamma     <- 10
#' rho       <- 0.5
#' 
#' g         <- 5 
#' probs     <- c(0.36, 0.28, 0.20, 0.12, 0.04)
#' Eff       <- p * probs
#' a         <- 5
#' B         <- a**(0:(g-1))-1
#' Z         <- matrix(0, nrow = p, ncol = g)
#' imax      <- 0
#' imin      <- 1
#' 
#' for (k in 1:g) {
#'     imin <- imax+1
#'     imax <- imax+Eff[k]
#'     Z[imin:imax, k] <- 1
#' }
#' Z <- Z[sample(1:p, p), ]
#' if (g>1) {
#'     Beta <- rnorm(p, mean = c(Z%*%B), sd = gamma)
#' } else {
#'     Beta <- rnorm(p, mean = B, sd = gamma)
#' }
#' 
#' theta0 <- NULL # c(intercept, B, probs, sigma^2, gamma^2)
#' Z0     <- NULL # apply(Z, 1, which.max)-1
#' 
#' gmax <- 7
#' 
#' ## Prediction
#' eps  <- rnorm(n, mean = 0, sd = sigma)
#' X    <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' Y    <- as.numeric(intercept+X%*%Beta+eps)
#' tt   <- system.time(mod <- fitClere(y = Y, x = X, g = gmax, 
#'                         analysis = analysis,algorithm = algorithm,
#'                         plotit = plotit, 
#'                         sparse = FALSE,nItEM = nItEM, 
#'                         nBurn = nBurn, nItMC = nItMC, 
#'                         nsamp = nsamp, theta0 = theta0, Z0 = Z0) )
#' plot(mod)
#' Yv <- predict(object = mod, newx = X)
#' 
fitClere <- function(
  y,
  x,
  g = 1,
  nItMC = 50,
  nItEM = 1000,
  nBurn = 200,
  dp = 5,
  nsamp = 200,
  maxit = 500,
  tol = 1e-3,
  nstart = 2,
  parallel = FALSE,
  seed = NULL,
  plotit = FALSE,
  sparse = FALSE,
  analysis = "fit",
  algorithm = "SEM",
  theta0 = NULL,
  Z0 = NULL
) {
  toyfit_fun <- function(ns) {
    if (is.null(seed)) {
      iseed <- sample(1:1000, 1)
    } else {
      iseed <- as.integer(seed * ns)
    }
    ## Added here since seed is no longer really passed as an argument
    set.seed(iseed)
    clereObj <- methods::new("Clere",
      y = y, x = x, g = g, nItMC = nItMC, nItEM = nItEM, nBurn = nBurn,
      dp = dp, nsamp = nsamp, sparse = sparse, analysis = analysis,
      algorithm = algorithm, initialized = FALSE, maxit = maxit, tol = tol, seed = iseed
    )
    ## Dimensions
    clereObj@n <- nrow(clereObj@x)
    clereObj@p <- ncol(clereObj@x)
    if (analysis == "fit" & !is.null(theta0)) {
      K <- length(theta0)
      if (K == 2 * g + 3) {
        clereObj@intercept <- theta0[1]
        clereObj@b <- theta0[2:(g + 1)]
        clereObj@pi <- theta0[(g + 2):(2 * g + 1)]
        clereObj@sigma2 <- theta0[2 * g + 2]
        clereObj@gamma2 <- theta0[2 * g + 3]
        if (abs(sum(clereObj@pi) - 1) >= 1e-10) {
          stop("[Clere:fitClere] Invalid proportions for initial vector of parameters.\n\tPlease check your initial guess of parameters", call. = FALSE)
        }
        if (clereObj@sigma2 <= 0 | clereObj@gamma2 <= 0) {
          stop("[Clere:fitClere] Non negative variances are not allowed!\n\tPlease check your initial guess of parameters", call. = FALSE)
        }
      } else {
        stop("[Clere:fitClere] The required number of groups is not consistent with your initial guess of parameters.\n\tPlease check your initial guess of parameters", call. = FALSE)
      }

      ## Check the initial partition
      if (is.null(Z0)) {
        stop("[Clere:fitClere] An initial partition should also be given with the initial set of parameters", call. = FALSE)
      } else {
        if (length(Z0) != clereObj@p) {
          stop(paste0("[Clere:fitClere] The length of the initial partition should be ", clereObj@p), call. = FALSE)
        }
        tz <- table(Z0)
        nz <- as.numeric(names(tz))
        if (min(nz) > 0 | max(nz) >= g) {
          stop(paste0("[Clere:fitClere] The partition only assigns numbers between 0 and ", g - 1, ", not between ", min(nz), " and ", max(nz), " as specifed by the user"), call. = FALSE)
        }
      }
      clereObj@Z0 <- Z0
      clereObj@initialized <- TRUE
    }
    .Call("clere", clereObj, PACKAGE = "clere")
    ## Last refinement
    ## Clean the output
    clereObj@Bw[!is.finite(clereObj@Bw)] <- NA
    clereObj@Zw[!is.finite(clereObj@Zw)] <- NA

    clereObj@g <- length(clereObj@b)
    ## Rename matrix P
    colnames(clereObj@P) <- paste("Group", 1:clereObj@g)
    rownames(clereObj@P) <- colnames(x)

    ## Rename output columns
    K <- seq_along(clereObj@b)
    colnames(clereObj@theta) <- c("intercept", paste0("b", K), paste0("pi", K), "sigma2", "gamma2", "CLL")
    if (plotit) {
      graphics::plot(clereObj)
    }

    clereObj
  }

  analysis <- match.arg(arg = analysis, choices = c("fit", "aic", "bic", "icl"))
  ## replication now
  if (parallel) {
    Models <- parallel::mclapply(seq_len(nstart), function(startingPoint) {
      mod <- toyfit_fun(startingPoint)
    }, mc.cores = min(parallel::detectCores(), nstart))
  } else {
    Models <- sapply(seq_len(nstart), function(startingPoint) {
      mod <- toyfit_fun(startingPoint)
    })
  }

  switch(EXPR = analysis,
    "fit" = {
      Models[[which.max(sapply(Models, function(m) {
        summary(m)@likelihood
      }))]]
    },
    "aic" = {
      Models[[which.min(sapply(Models, function(m) {
        summary(m)@AIC
      }))]]
    },
    "bic" = {
      Models[[which.min(sapply(Models, function(m) {
        summary(m)@BIC
      }))]]
    },
    "icl" = {
      Models[[which.min(sapply(Models, function(m) {
        summary(m)@ICL
      }))]]
    },
    stop(paste0('[Clere:fitClere] Please select one of the available option for "analysis": "fit", "aic", "bic" or "icl"'), call. = FALSE)
  )
}
