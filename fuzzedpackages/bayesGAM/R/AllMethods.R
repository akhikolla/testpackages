# -----------------GLM family methods---------------------------- #

# ----------------- number of parameters for a distribution ---------------------------- #

setGeneric("numParams", function(object, ...) {
  standardGeneric("numParams")
})

setMethod("numParams", "distribution",
          function(object) {
  object@num_params
})


# -----------------normal distribution class---------------------------- #


setMethod("initialize", "normalDistribution",
          function(.Object, name="Normal",
                   distnum = 1L,
                   num_params = 2L,
                   param_names=c("mu", "sigma"),
                   param_values=c(0, 1e6), ...) {
            
            if (missing(name))
              name <- "Normal"
            if (missing(distnum))
              distnum <- 1L
            if (missing(param_names))
              param_names <- c("mu", "sigma")
            if (missing(param_values))
              param_values <- c(0, 1e6)
            if (missing(num_params))
              num_params <- length(param_values)

            if (!all.equal(param_names, c("mu", "sigma")))
              stop("invalid param_names for Normal.  Should be mu sigma")
            if (param_values[2] <= 0)
              stop("sigma parameter must be positive for Normal distribution")
            if (name != "Normal")
              stop("must be named Normal")
            .Object@param_names <- param_names
            .Object@param_values <- param_values
            .Object@name <- name
            .Object@distnum <- distnum
            .Object@num_params <- num_params
            return(.Object)
          })



# -----------------Student t distribution class---------------------------- #

setMethod("initialize", "tDistribution",
          function(.Object, name="tdist",
                   distnum=2L,
                   num_params=3L,
                   param_names=c("nu", "mu", "sigma"),
                   param_values=c(1, 0, 25), ...) {
            
            if (missing(name))
              name <- "tdist"
            if (missing(distnum))
              distnum <- 2L
            if (missing(param_names))
              param_names <- c("nu", "mu", "sigma")
            if (missing(param_values))
              param_values <- c(1, 0, 25)
            if (missing(num_params))
              num_params <- length(param_values)

            if (!all.equal(param_names, c("nu", "mu", "sigma")))
              stop("invalid param_names for student-t  Should be nu, mu, sigma")
            if (param_values[1] <= 0)
              stop("nu parameter must be positive for student-t distribution")
            if (param_values[3] <= 0)
              stop("sigma parameter must be positive for student-t distribution")
            if (name != "tdist")
              stop("must be named tdist")

            .Object@param_names <- param_names
            .Object@param_values <- param_values
            .Object@name <- name
            .Object@distnum <- distnum
            .Object@num_params <- num_params
            return(.Object)
          })

# -----------------GLM family class---------------------------- #

setMethod("initialize", "glmFamily",
          function(.Object,
                   famnum,
                   famname,
                   linknum,
                   linkname,
                   mixed, ...) {
            
            if (missing(famnum))
              # default to Gaussian
              famnum <- 1L
            if (missing(famname)) {
              if (famnum == 1L) {
                famname <- "gaussian"
              } else if (famnum == 2L) {
                famname <- "binomial"
              } else if (famnum == 3L) {
                famname <- "poisson"
              } else {
                stop("family not supported")
              }
            }

            if (missing(linkname)) {
              if (linknum == 1L) {
                linkname <- "identity"
              } else if (linknum == 2L) {
                linkname <- "log"
              } else if (linknum == 3L) {
                linkname <- "inverse"
              } else if (linknum == 4L) {
                linkname <- "logit"
              } else if (linknum == 5L) {
                linkname <- "probit"
              } else if (linknum == 6L) {
                linkname <- "cauchit"
              } else if (linknum == 7L) {
                linkname <- "cloglog"
              } else if (linknum == 8L) {
                linkname <- "sqrt"
              } else {
                stop("link not supported")
              }
            }

            if (missing(mixed)) {
              mixed <- FALSE
            }

            .Object@famnum <- famnum
            .Object@famname <- famname
            .Object@linknum <- linknum
            .Object@linkname <- linkname
            .Object@mixed <- mixed
            return(.Object)
          })


# -----------------GLM model class---------------------------- #

setMethod("initialize", "glmModel",
          function(.Object, X, Z, Zint, Znp, y, p, r, q, n, has_intercept, zvars, names_beta, names_u, names_y,
                   prior_selections, knots, npargs, npterms, random_intercept, basis, npdegree, sub_form, 
                   call, offset, ...) {
            .Object <- callNextMethod(.Object, ...)
            

            if (missing(y)) {
              stop("y must be provided to fit model")
            } else {
              # n <- length(y)
              y <- as.matrix(y)
              n <- nrow(y)

              multresponse <- ncol(y) > 1
            }

            if (missing(X)) {
              stop("design matrix X needed to fit model")
            } else {
              if (nrow(X) != n) stop("number of rows in X must match length of y")
            }
            if (missing(p)) {
              p <- ncol(X)
            } else if (p != ncol(X)) {
              stop("p must match the number of columns in X")
            }

            if (missing(Z)) {
              Z <- Zint <- Znp <- matrix(0, nrow=0, ncol=0)
              #Zlst <- list(Z)
              max_col <- 0L
            } 

            nk <- max_col <- ncol(Z)
            if (missing(zvars)) {
              zvars <- ncol(Z)
            }
            if (missing(random_intercept)) {
              if (length(Z) == 0) {
                random_intercept <- FALSE
                Zint <- matrix(0, nrow=0, ncol=0)
              } else {
                stop("must specify if random intercept is present")
              }
            }
            
            # special case for random intercept
            if (missing(Zint) & random_intercept) {
              Zint <- Z[, 1:zvars[1]]
            } 
            
            if (missing(Znp) & random_intercept) {
              Znp <- Z[, -c(1:zvars[1])]
            } else if (missing(Znp) & !random_intercept) {
              Znp <- Z
            }
            
            if (missing(names_beta)) {
                if (length(colnames(X)) == p) {
                  names_beta <- colnames(X)
                } else {
                  names_beta <- paste0("beta", 1:p)
                }
            } else if (length(names_beta) != p) {
              stop("length of names_beta must match number of columns in X")
            }
            if (missing(names_u)) {
              if (nk == 0) {
                names_u <- character(0)
              } else {
                names_u <- paste0("u", 1:nk)
              }
            } else {
              if (length(names_u) != nk) stop("length of names_u must match number of columns in Z")
            }
            if(.Object@mixed) {
              if (nk == 0) stop("Z matrix not supplied for mixed model")
            }
            if (missing(names_y)) {
              names_y <- "depvar"
            }

            # assign prior if missing. Otherwise, verify that prior is valid for model
            if (missing(q)) {
              q <- 0L
            }
            if (missing(r)) {
              if (.Object@famnum == 1L) {
                # r <- 1L
                r <- ncol(as.matrix(y))
              } else {
                r <- 0L
              }
            }

            if (missing(basis)) {
              basis <- character()
            }
            if (missing(npdegree)) {
              npdegree <- integer(0L)
            }

            if (missing(prior_selections)) {

              prior_selections <- new("prior",
                                      famnum = .Object@famnum,
                                      mixed = .Object@mixed,
                                      p = p,
                                      q = q,
                                      r = r,
                                      random_intercept = random_intercept)
            }

            # check matching family and mixed
            if (.Object@famnum != prior_selections@famnum) {
              stop("family selected for prior does not match family for model")
            }
            if (.Object@mixed != prior_selections@mixed) {
              stop("mixed model form (including Z) must be consistent between priors and model selection")
            }

            if (missing(knots)) {
              knots <- list()
            }

            if (missing(npterms)) {
              npterms <- character()
            }

            if (missing(npargs)) {
              npargs <- list()
            }
            
            .Object@X <- X
            .Object@Z <- Z
            .Object@Zint <- Zint
            .Object@Znp <- Znp
            .Object@max_col <- max_col
            .Object@y <- y
            .Object@p <- p
            .Object@r <- r
            .Object@q <- q
            .Object@n <- n
            .Object@has_intercept <- has_intercept
            .Object@zvars <- zvars
            .Object@names_beta <- names_beta
            .Object@names_u <- names_u
            .Object@names_y <- names_y
            .Object@prior <- prior_selections
            .Object@knots <- knots
            .Object@npargs <- npargs
            .Object@npterms <- npterms
            .Object@sub_form <- sub_form
            .Object@random_intercept <- random_intercept
            .Object@basis <- basis
            .Object@npdegree <- npdegree
            .Object@multresponse <- multresponse
            .Object@call <- call
            .Object@offset <- offset
            return(.Object)
          })


# ----------------- set priors ---------------------------- #

# initialize prior object
# famnum: family number
# mixed:  logical whether to include Z
# beta:  list of priors for beta
# eps:  list of priors for eps
# lambda:  list of priors for lambda
# a: list of priors for a
# p:  number of beta params
# r:  number of eps params
# q:  number of lambda params
setMethod("initialize", "prior",
          function(.Object, famnum, mixed, beta, eps, lambda, a, p, r, q, random_intercept, ...) {

            # for a off-diagonal multivariate response
            # random_intercept <- .Object@random_intercept

            if (missing(famnum)) {
              famnum <- 1L
            }
            if (missing(mixed)) {
              mixed <- FALSE
            }

            # use eps for Normal only
            if (missing(r)) {
              if (famnum == 1L) {
                r <- 1L
              } else {
                r <- 0L
              }
            }

            # mixed
            if (mixed) {
              if (missing(q)) {
                stop("must supply number of parameters q")
              }
            } else {
              q <- 0L
            }

            # always supply p
            if (missing(p)) {
              stop("must supply number of fixed effect parameters p")
            }

            # assign default distributions
            if (famnum == 1L & missing(eps)) {
              eps <- replicate(r, new("tDistribution", param_values = c(3, 0, 1)), simplify=FALSE)
            }
            if (famnum %in% c(2L, 3L)) {
              eps <- replicate(r, new("distribution"), simplify=FALSE)
            }

            if (r == 0) {
              names(eps) <- paste0("eps", pmax(1L, r))
            } else {
              names(eps) <- paste0("eps", 1:r)
            }

            if (missing(beta)) {
              # check for multiple response
              if (r <= 1L) {
                beta <- replicate(p, new("normalDistribution", param_values = c(0, 1e6)), simplify=FALSE)
                names(beta) <- paste0("beta", 1:p)
              } else {
                beta <- replicate(p*r, new("normalDistribution", param_values = c(0, 1e6)), simplify=FALSE)
                names(beta) <- paste0("beta", rep(1:r, each=p), "_", rep(1:p))
              }
            } else {
              if (length(beta) != r*p) {
                stop("invalid number of beta parameters")
              }
            }

            if (missing(lambda) | q == 0) {
              if (mixed) {
                if (r <= 1L) {
                  lambda <- replicate(q, new("tDistribution", param_values = c(3, 0, 1)), simplify=FALSE)
                  names(lambda) <- paste0("lambda", 1:pmax(1L, q))
                } else {
                  lambda <- replicate(q*r, new("tDistribution", param_values = c(3, 0, 1)), simplify=FALSE)
                  names(lambda) <- paste0("lambda", rep(1:r, each=q), "_", rep(1:q))
                }
              } else {
                lambda <- replicate(1, new("distribution"), simplify=FALSE)
                names(lambda) <- paste0("lambda", 1:pmax(1L, q))
              }
            } else {
              if (length(lambda) != r*q) {
                stop("invalid number of lambda parameters")
              }
            }

            if (missing(a)) {
              if (random_intercept & (r>1)) {
                a <- replicate(r*(r-1)/2, new("normalDistribution", param_values = c(0, 1e6)), simplify=FALSE)
                names(a) <- paste0("a", 1:(r*(r-1)/2))
              } else {
                a <- replicate(1, new("distribution"), simplify=FALSE)
              }
            }

            .Object@famnum <- famnum
            .Object@mixed <- mixed
            .Object@beta <- beta
            .Object@eps <- eps
            .Object@lambda <- lambda
            .Object@a <- a
            return(.Object)
          })

# methods to get priors
setGeneric("getPrior", function(object, ...) {
  standardGeneric("getPrior")
})

setMethod("getPrior", signature(object = "prior"),
          function(object, params="beta", type="param_values", transpose=TRUE, FUN=sapply, ...) {

            if (!type %in% slotNames("distribution")) {
              stop("invalid slot name for distribution object")
            }
            if (!params %in% c("beta", "eps", "lambda", "a")) {
              stop("invalid parameter selection.  must be beta, eps, lambda, or a")
            }

            if (params == "beta") {
              retval <- FUN(object@beta, slot, type, ...)
            } else if (params == "eps") {
              retval <- FUN(object@eps, slot, type, ...)
            } else if (params == "lambda") {
              retval <- FUN(object@lambda, slot, type, ...)
            } else if (params == "a") {
              retval <- FUN(object@a, slot, type, ...)
            } else {
              stop("invalid prior selections")
            }

            if (transpose) {
              retval <- t(retval)
            }

            # if (type == "param_values") {
            #   retval <- do.call(rbind, retval)
            # }

            return(retval)
          })

# ----------------- parse list/distribution for priors---------------------------- #


setGeneric("parsePrior", function(object, ...) {
  standardGeneric("parsePrior")
})

# parsing prior parameter when given a list
setMethod("parsePrior", signature(object="list"),
          function(object, num_dist, ...) {
            num_dist_from_list <- length(object)
            if (num_dist_from_list == 1) {
              if (is(object[[1]], "distribution")) {
                res <- replicate(num_dist, object[[1]], simplify=FALSE)
              } else {
                stop("must provide a list of distribution objects")
              }
            } else if (num_dist_from_list == num_dist) {
              if (!all(sapply(object, is, "distribution"))) {
                stop("all elements in list must be distributions")
              } else {
                res <- object
              }
            } else {
              stop("number of elements in list must equal one or correct number of parameters")
            }
            return(res)
          })

# parsing prior parameter when given a distribution
setMethod("parsePrior", signature(object="distribution"),
          function(object, num_dist, ...) {
            replicate(num_dist, object, simplify=FALSE)
          })


# ---------------- set prior for glm model ---------------------------- #

setGeneric("setPrior", function(object, ...) {
  standardGeneric("setPrior")
})

# beta, eps, lambda are lists of distributions
# added a prior for offdiagonal multivariate modeling
setMethod("setPrior", signature(object="glmModel"),
          function(object, famnum, mixed, p, q, r, random_intercept, beta, eps, lambda, a, ...) {
            if (missing(famnum)) {
              famnum <- object@famnum
            }
            if (missing(mixed)) {
              mixed <- object@mixed
            }
            if (missing(p)) {
              p <- object@p
            }
            if (missing(q)) {
              q <- object@q
            }
            if (missing(r)) {
              r <- object@r
            }

            # current priors
            prior <- object@prior

            if (missing(beta) | length(beta) == 0) {
              beta <- prior@beta
            } else {
              beta <- parsePrior(beta, p*r)
              names(beta) <- paste0("beta", rep(1:r, each=p), "_", rep(1:p))
            }

            if (missing(eps) | length(eps) == 0) {
              eps <- prior@eps
            } else {
              eps <- parsePrior(eps, r)
            }

            if (missing(lambda) | length(lambda) == 0) {
              lambda <- prior@lambda
            } else {
              lambda <- parsePrior(lambda, q*r)
              names(lambda) <- paste0("lambda", rep(1:r, each=q), "_", rep(1:q))
            }

            if (missing(a) | length(a) == 0) {
              a <- prior@a
            } else {
              a <- parsePrior(a, r*(r-1)/2)
              names(a) <- paste0("a", 1:(r*(r-1)/2))
            }

            # check family of glm model
            prior_selections <- new("prior",
                                    famnum = famnum,
                                    mixed = mixed,
                                    p = p,
                                    q = q,
                                    r = r,
                                    beta = beta,
                                    eps = eps,
                                    lambda = lambda,
                                    a = a,
                                    random_intercept = random_intercept)

            object@prior <- prior_selections
            return(object)
})



# ---------------- show prior from model ---------------------------- #

setGeneric("getDistName", function(object, ...) {
  standardGeneric("getDistName")
})

setMethod("getDistName", signature(object="integer"),
          function(object, ...) {

  res <- sapply(object, function(xx) {
    if (xx == 1L) {
      nm <- "Normal"
    } else if (xx == 2L) {
      nm <- "t"
    } else {
      stop("unknown distribution")
    }
    nm
  })
  res

})

#' Display the priors used in \code{bayesGAM}
#' 
#' Prints a list of priors for \eqn{\beta}, \eqn{\lambda}, \eqn{\epsilon}, and \eqn{a}, where applicable. 
#' 
#' @param object Object of type \code{bayesGAMfit} generated from \code{bayesGAM}
#' @param params character vector of the names of parameters to return
#' \itemize{
#' \item \eqn{\beta} beta
#' \item \eqn{\epsilon} eps
#' \item \eqn{\lambda} lambda
#' \item \eqn{a] a}
#' }
#' @return none
#' @export
#' @name showPrior
#' @examples
#' require(stats); require(graphics)
#' f <- bayesGAM(weight ~ np(height), data = women, 
#'               family = gaussian, iter = 500, chains = 1)
#' showPrior(f)
NULL

#' @rdname showPrior
#' @param ... Additional arguments for \code{showPrior}
#' @export
setGeneric("showPrior", function(object, ...) {
  standardGeneric("showPrior")
})

#' @rdname showPrior
#' @export
setMethod("showPrior", signature(object="bayesGAMfit"),
  function(object, params = "all") {

    prior_obj <- object@model@prior
    famnum <- prior_obj@famnum
    mixed <- prior_obj@mixed
    random_intercept <- object@model@random_intercept

    # for multresponse with multiple y
    multresponse <- object@model@multresponse
    r <- object@model@r

    beta_distnum <- getPrior(prior_obj, "beta", "distnum", transpose=F)
    beta_distname <- getDistName(beta_distnum)
    beta_param <- getPrior(prior_obj, "beta", "param_values", transpose=F, FUN=lapply)
    beta_num_params <- getPrior(prior_obj, "beta", "num_params", transpose=F)

    beta_nms <- names(prior_obj@beta)
    beta_param_str <- mapply(pastelim, x=beta_param, nparam=beta_num_params, MoreArgs=list(collapse=', '))
    beta_ret <- paste0(beta_nms, " ~ ", beta_distname, "(", beta_param_str, ")", "\n")

    if (famnum == 1L) {
      eps_distnum <- getPrior(prior_obj, "eps", "distnum", transpose=F)
      eps_distname <- getDistName(eps_distnum)
      eps_param <- getPrior(prior_obj, "eps", "param_values", transpose=F, FUN=lapply)
      eps_num_params <- getPrior(prior_obj, "eps", "num_params", transpose=F)

      eps_nms <- names(prior_obj@eps)
      eps_param_str <- mapply(pastelim, x=eps_param, nparam=eps_num_params, MoreArgs=list(collapse=', '))
      eps_ret <- paste0(eps_nms, " ~ ", eps_distname, "(", eps_param_str, ")", "\n")

    }

    if (mixed) {
      lambda_distnum <- getPrior(prior_obj, "lambda", "distnum", transpose=F)
      lambda_distname <- getDistName(lambda_distnum)
      lambda_param <- getPrior(prior_obj, "lambda", "param_values", transpose=F, FUN=lapply)
      lambda_num_params <- getPrior(prior_obj, "lambda", "num_params", transpose=F)

      lambda_nms <- names(prior_obj@lambda)
      lambda_param_str <- mapply(pastelim, x=lambda_param, nparam=lambda_num_params, MoreArgs=list(collapse=', '))
      lambda_ret <- paste0(lambda_nms, " ~ ", lambda_distname, "(", lambda_param_str, ")", "\n")

    }

    # off diagonal multresponse
    if (multresponse & (r > 1) & random_intercept) {
      a_distnum <- getPrior(prior_obj, "a", "distnum", transpose=F)
      a_distname <- getDistName(a_distnum)
      a_param <- getPrior(prior_obj, "a", "param_values", transpose=F, FUN=lapply)
      a_num_params <- getPrior(prior_obj, "a", "num_params", transpose=F)

      a_nms <- names(prior_obj@a)
      a_param_str <- mapply(pastelim, x=a_param, nparam=a_num_params, MoreArgs=list(collapse=', '))
      a_ret <- paste0(a_nms, " ~ ", a_distname, "(", a_param_str, ")", "\n")
    }

    # strings to return
    if (params %in% c("beta", "all")) {
      cat(beta_ret)
    }
    if (famnum == 1L & params %in% c("eps", "epsilon", "all")) {
      cat(eps_ret)
    }
    if(mixed & params %in% c("lambda", "all")) {
      cat(lambda_ret)
    }
    if (multresponse & (r > 1) & random_intercept & params %in% c("a", "all")) {
      cat(a_ret)
    }
})


# ---------------- store fit ---------------------------- #


setMethod("initialize", "bayesGAMfit",
          function(.Object, results, model, offset, spcontrol, ...) {

            if (missing(offset)) {
              offset <- numeric(0)
            }
            if (missing(spcontrol)) {
              spcontrol <- list(qr=TRUE, mvindep=FALSE)
            }

            .Object@results <- results
            .Object@model <- model
            .Object@offset <- offset
            .Object@spcontrol <- spcontrol
            .Object@mcmcres <- matrix(, nrow=0L, ncol=0L)
            .Object@pdata <- data.frame()
            return(.Object)
          })


# ----------------- get distribution details ---------------------------- #

setGeneric("getDistribution", function(object, ...) {
  standardGeneric("getDistribution")
})

# extract distribution from distribution object
setMethod("getDistribution", signature(object="distribution"),
          function(object, max_length=3) {
            vals <- object@param_values
            vals <- c(vals, rep(-Inf, max_length - length(vals)))
            retval <- list(param_nums = vals)
            retval$distnum <- object@distnum
            return(retval)
          })

# extract distribution matrix from list of distributions
setMethod("getDistribution", signature(object="list"),
          function(object, max_length=3) {
            res1 <- sapply(object, function(xx) {
              getDistribution(xx)$param_nums
            })

            retval <- list(param_nums = t(res1))

            res2 <- sapply(object, function(xx) {
              getDistribution(xx)$distnum
            }, simplify = "integer")
            retval$distnum <- as.integer(res2)
            return(retval)
          })




# ----------------- Fit MCMC model---------------------------- #

setGeneric("bayesGAMfit", function(object, ...) {
  standardGeneric("bayesGAMfit")
})

setMethod("bayesGAMfit", signature(object="glmModel"),
          function(object, linknum, offset, spcontrol, ...) {

            use_Z <- length(object@Z) > 0
            
            if (is.null(spcontrol$mvindep) | length(spcontrol$mvindep) == 0) {
              mvindep <- 0
            } else {
              mvindep <- spcontrol$mvindep*1
            }
            
            if (is.null(spcontrol$qr) | length(spcontrol$qr) == 0) {
              qr <- 1
            } else {
              qr <- spcontrol$qr*1
            }

            # convert to numeric if single response
            if (!object@multresponse) {
              y <- as.numeric(object@y)
            } else {
              y <- object@y
            }

            beta_distnum <- getPrior(object@prior, "beta", "distnum", transpose=F, simplify=TRUE)
            beta_param <- getPrior(object@prior, "beta", "param_values", transpose=F, FUN=lapply)
            beta_param <- t(sapply(beta_param, '[', 1:max(sapply(beta_param, length))))
            beta_param[is.na(beta_param)] <- 0

            eps_distnum <- getPrior(object@prior, "eps", "distnum", transpose=F, simplify=TRUE)
            eps_param <- getPrior(object@prior, "eps", "param_values", transpose=F, FUN=lapply)
            eps_param <- t(sapply(eps_param, '[', 1:max(sapply(eps_param, length))))
            eps_param[is.na(eps_param)] <- 0

            lambda_distnum <- getPrior(object@prior, "lambda", "distnum", transpose=F, simplify=TRUE)
            lambda_param <- getPrior(object@prior, "lambda", "param_values", transpose=F,FUN=lapply)
            lambda_param <- t(sapply(lambda_param, '[', 1:max(sapply(lambda_param, length))))
            lambda_param[is.na(lambda_param)] <- 0

            a_distnum <- getPrior(object@prior, "a", "distnum", transpose=F, simplify=TRUE)
            a_param <- getPrior(object@prior, "a", "param_values", transpose=F,FUN=lapply)
            a_param <- t(sapply(a_param, '[', 1:max(sapply(a_param, length))))
            a_param[is.na(a_param)] <- 0
            
            # Zint dimensions
            if (object@random_intercept) {
              Zint <- object@Zint
              a_num_offdiagonal <- object@r*(object@r-1)/2
              anum <- c(a_distnum, 99999)
            } else {
              Zint <- matrix(0, nrow=0, ncol=0)
              a_num_offdiagonal <- 0L
              anum <- integer(0L)
              a_param <- matrix(0, nrow=0, ncol=0)
            }
            
            # Znp dim
            if (ncol(object@Znp) > 0) {
              Znp <- object@Znp
            } else {
              Znp <- matrix(0, nrow=0, ncol=0)
            }
            
            # Z dim
            if (ncol(object@Z) > 0) {
              lambdanum <- c(lambda_distnum, 99999)
              lambda_max_params <- ncol(lambda_param)
              zvars <- c(object@zvars, -1e6)
              Z <- object@Z
            } else {
              lambdanum <- integer(0L)
              lambda_max_params <- 0L
              lambda_param <- matrix(0, nrow=0, ncol=0)
              zvars <- integer(0L)
              Z <- matrix(0, nrow=0, ncol=0)
            }
            
            dat <- list(N = nrow(object@X),
                        p = ncol(object@X),
                        y = y,
                        X = object@X,
                        r = object@r,
                        q = object@q,
                        nk = ncol(object@Z),
                        ny = ncol(y),
                        zvars = zvars,
                        Z = Z,
                        Zint = Zint, 
                        Znp = Znp, 
                        nnp = ncol(object@Znp), 
                        nrandint = ncol(object@Zint), 
                        max_col = object@max_col,
                        qr = qr + 0,
                        qrsplit = 0,
                        mvindep = mvindep + 0,
                        randint = object@random_intercept + 0, 
                        randeff = (ncol(object@Znp) > 0) + 0, 
                        famnum = object@famnum,
                        offset = offset,
                        linknum = object@linknum,
                        lambdanum = lambdanum,
                        lambda_max_params = lambda_max_params,
                        lambda_param = lambda_param,
                        epsnum = c(eps_distnum, 99999),
                        eps_max_params = ncol(eps_param),
                        eps_param = eps_param,
                        betanum = c(beta_distnum, 99999),
                        beta_max_params = ncol(beta_param),
                        beta_param = beta_param,
                        a_num_offdiagonal = a_num_offdiagonal, # number of off-diagonal for multresponse
                        anum = anum,
                        a_max_params = ncol(a_param),
                        a_param = a_param)
            
            if (!object@multresponse & object@famnum == 1) {
              res <- rstan::sampling(stanmodels$glmm_continuous, data = dat, ...)
            } else if (!object@multresponse & object@famnum %in% c(2, 3)) {
              res <- rstan::sampling(stanmodels$glmm_discrete, data = dat, ...)
            } else if (object@multresponse & object@famnum == 1) {
              res <- rstan::sampling(stanmodels$multresponse_continuous, data = dat, ...)
            } else if (object@multresponse & object@famnum %in% c(2, 3)) {
              res <- rstan::sampling(stanmodels$multresponse_discrete, data = dat, ...)
            }

            else {
              stop("model not supported")
            }

            return(res)
          })

# ----------------- get design matrix  ---------------------------- #

#' @rdname getDesign
#' @param ... Additional arguments for \code{getDesign}
setGeneric("getDesign", function(object, ...) {
  standardGeneric("getDesign")
})

#' Design matrices from a \code{bayesGAMfit} object
#'
#' Contains the design matrices produced for model fitting. The fixed effects design matrix \code{X} or random effects design matrix \code{Z} can be specified.
#'
#' @param object Object of type \code{bayesGAMfit} generated from \code{bayesGAM}
#' @param type Character for fixed effect design matrix \code{X} or random effects design matrix \code{Z}
#' @return Contents of \code{stanfit} results
#' @name getDesign
#' @export
#' @examples
#' require(stats); require(graphics)
#' f <- bayesGAM(weight ~ np(height), data = women, family = gaussian, 
#'               iter = 500, chains = 1)
#' getDesign(f, "Z")
NULL

#' @rdname getDesign
#' @export
setMethod("getDesign", signature(object="bayesGAMfit"),
          function(object, type="X") {
            if (type == "X") {
              retval <- object@model@X
            } else if (type == "Z") {
              retval <- object@model@Z
            } else {
             retval <- matrix(, nrow=0L, ncol=0L)
            }
            return(retval)
          })

#' @rdname getDesign
setMethod("getDesign", signature(object="glmModel"),
          function(object, type="X") {
            if (type == "X") {
              retval <- object@X
            } else if (type == "Z") {
              retval <- object@Z
            } else {
              retval <- matrix(, nrow=0L, ncol=0L)
            }
            return(retval)
          })

# ------------------ get model slots ---------------------------- # 

#' @param ... Additional arguments for \code{getModelSlots}
#' @rdname getModelSlots
setGeneric("getModelSlots", function(object, ...) {
  standardGeneric("getModelSlots")
})

#' Return one or slots from the \code{Stan} model in \code{bayesGAM}
#'
#' Contains the objects and parameters passed to Stan in object of type \code{glmModel}, contained in object type \code{bayesGAMfit}
#'
#' @param object Object of type \code{bayesGAMfit}
#' @param name Character name of slot in \code{glmModel}
#' \itemize{
#' \item X Fixed effects design matrix
#' \item Z Random effects design matrix
#' \item Zlst list of individual random effects design matrices that, combined, form \code{Z}
#' \item Zarray array of individual random effects design matrices.  Used for multiple response models
#' \item max_col maximum number of columns of an individual \code{Z} matrix.  Padding for STAN
#' \item y numeric response matrix
#' \item p number of beta parameters
#' \item r number of eps parameters
#' \item q number of lambda parameters
#' \item n number of records in the dataset
#' \item has_intercept logical of whether the model includes an intercept term
#' \item zvars number of random effects variables
#' \item names_beta parameter names for beta
#' \item names_u parameter names for the random effects
#' \item names_y response names
#' \item prior `prior` object with priors used in the model
#' \item knots list of knots used in non-parametric functions
#' \item basis character indicating basis function.  \code{tps} for thin-plate splines and \code{trunc.poly} for truncated polynomial
#' \item npargs arguments passed to non-parametric functions in the model
#' \item npterms variables used in non-parametric functions
#' \item sub_form formula with the \code{np} terms removed
#' \item random_intercept logical indicator of whether a random effects intercept is used
#' \item multresponse logical indicator of whether the model is multiple response
#' }
#' @return Contents of slot in \code{glmModel}
#' @name getModelSlots
#' @export
#' @examples
#' require(stats); require(graphics)
#' f <- bayesGAM(weight ~ np(height), data = women, family = gaussian, 
#'               iter = 500, chains = 1)
#' getModelSlots(f, "X")
NULL

#' @rdname getModelSlots
#' @export
setMethod("getModelSlots", signature(object="bayesGAMfit"), 
          function(object, name="X") {
            if (name=="X") {
              object@model@X           
            } else if (name == "Z") {
              object@model@Z
            } else if (name == "max_col") {
              object@model@max_col
            } else if (name == "y") {
              object@model@y
            } else if (name == "p") {
              object@model@p
            } else if (name == "r") {
              object@model@r
            } else if (name == "q") {
              object@model@q
            } else if (name == "n") {
              object@model@n
            } else if (name == "has_intercept") {
              object@model@has_intercept
            } else if (name == "zvars") {
              object@model@zvars
            } else if (name == "names_beta") {
              object@model@names_beta
            } else if (name == "names_u") {
              object@model@names_u
            } else if (name == "names_y") {
              object@model@names_y
            } else if (name == "prior") {
              object@model@prior
            } else if (name == "knots") {
              object@model@knots
            } else if (name == "basis") {
              object@model@basis
            } else if (name == "npargs") {
              object@model@npargs
            } else if (name == "npterms") {
              object@model@npterms
            } else if (name == "random_intercept") {
              object@model@random_intercept
            } else if (name == "multresponse") {
              object@model@multresponse
            } else if (name == "sub_form") {
              object@model@sub_form
            }
          })

# ------------------ get model slots ---------------------------- # 

#' @rdname getStanResults
setGeneric("getStanResults", function(object) {
  standardGeneric("getStanResults")
})

#' Returns the \code{stanfit} object generated by \pkg{rstan}
#'
#' Contains the full content of the \code{stanfit} object
#'
#' @param object Object of type \code{bayesGAMfit} returned from \code{bayesGAM}
#' @return Contents of \code{stanfit} results
#' @name getStanResults
#' @export
#' @examples
#' require(stats); require(graphics)
#' f <- bayesGAM(weight ~ np(height), data = women, family = gaussian, 
#'               iter = 500, chains = 1)
#' sres <- getStanResults(f)
#' plot(sres) # rstan method
NULL

#' @rdname getStanResults
#' @export
setMethod("getStanResults", signature(object="bayesGAMfit"), 
          function(object) {
            object@results
          })

# ---------------- rstan extract  ---------------------------- #

setMethod("extract", signature("bayesGAMfit"),
          function(object, ...) {
            rstan::extract(as(object, "stanfit"), ...)
          })

# ---------------- fitted values ---------------------------- #

#' Extract Model Coefficients
#'
#' Method for \code{bayesGAMfit} objects.  Extracts the specified quantile of the posterior. 
#' The user may specify all or some of the parameters \eqn{\beta}, \eqn{\epsilon}, \eqn{\lambda}, \eqn{u}, \eqn{sigma}, \eqn{a}.
#'
#' @param object an object of class \code{bayesGAMfit}, usually a result of a call to \code{bayesGAM}.
#' @param params character vector of the names of parameters to return
#' \itemize{
#' \item \eqn{\beta} beta
#' \item \eqn{\epsilon} eps
#' \item \eqn{\lambda} lambda
#' \item \eqn{a] a}
#' }
#' @param FUN function from which to estimate coefficients. Default is \code{median}
#' @return Numeric vector of parameter point estimates based on the given \code{prob}, with a default of the median estimate.
#' @export
#' @name coefficients
#' @examples
#' require(stats); require(graphics)
#' f <- bayesGAM(weight ~ np(height), data = women, family = gaussian, 
#'               iter = 500, chains = 1)
#' coef(f, params=c("beta", "eps"))
NULL

#' @rdname coefficients
setMethod("coefficients", signature("bayesGAMfit"),
          function(object, params=c("beta", "eps", "lambda", "u", "sigma", "a"), FUN=median) {

            if (length(object@mcmcres) > 0) {
              sims <- object@mcmcres  
            } else {
              sims <- as.matrix(object@results)
            }
            
            vals <- apply(sims, 2, FUN)

            pars <- paste0("^", params)

            whichpars <- unique (grep(paste(pars,collapse="|"),
                                                 names(vals), value=TRUE))
            return(vals[whichpars])
          })

#' @rdname coefficients
#' @export
setMethod("coef", signature("bayesGAMfit"),
          function(object, params=c("beta", "eps", "lambda", "u", "sigma", "a"), FUN=median) {
            coefficients(object, params=params, FUN=FUN)
          })

#' Extract fitted values from a model fit by \code{bayesGAM}
#'
#' Method for \code{bayesGAMfit} objects.  Extracts the fitted values based on a specified quantile for the posterior distribution. The median is the default. 
#'
#' @param object an object of class \code{bayesGAMfit}, usually a result of a call to \code{bayesGAM}.
#' @param ... additional arguments to pass to \code{coefficients}
#' @return Numeric vector of fitted values
#' @export
#' @name fitted
#' @examples
#' require(stats); require(graphics)
#' f <- bayesGAM(weight ~ np(height), data = women, family = gaussian, 
#'               iter = 500, chains = 1)
#' plot(fitted(f), women$weight, type='o', xlab="fitted", ylab="actual")
NULL

#' @rdname fitted
setMethod("fitted", signature("bayesGAMfit"),
          function(object, ...) {

            X <- object@model@X
            Z <- object@model@Z
            multresponse <- object@model@multresponse
            
            # get inverse link function
            linkinfo <- make.link(object@model@linkname)
            invlink <- linkinfo$linkinv
            
            coef_beta <- coefficients(object, params="beta", ...)
            
            if (multresponse) {
              
              betacols <- strsplit(names(coef_beta),  
                                    split="\\[|\\,|\\]")
              nyvalsbeta <- sapply(betacols, function(xx) {
                as.numeric(xx[2])
              })
              
              if (object@model@mixed) {
                coef_u <- coefficients(object, params="u", ...)
                ucols <- strsplit(names(coef_u), 
                                  split="\\[|\\,|\\]")
                nyvalsu <- sapply(ucols, function(xx) {
                  as.numeric(xx[2])
                })
                
                coef_betau <- c(coef_beta, coef_u)
                nyvalsbetau <- c(nyvalsbeta, nyvalsu)
                Cgrid <- cbind(X, Z)
              } else {
                coef_betau <- coef_beta
                nyvalsbetau <- nyvalsbeta
                Cgrid <- X
              }
              
              fitvals <- lapply(sort(unique(nyvalsbetau)), function(i) {
                coef_temp <- coef_betau[nyvalsbetau == i]
                fv <- Cgrid %*% coef_temp
                invlink(fv)
              })
              
              fitvals <- do.call(cbind, fitvals)
              
              
            } else {
              # pull coefficients
              fitvals <- X%*%coef_beta
              
              if (object@model@mixed) {
                coef_u <- coefficients(object, params="u", ...)
                fitvals <- fitvals + Z%*%coef_u
              }
              
              invlink(fitvals)
            }
           

          })


# ---------------- bayesGAM plotting methods ---------------------------- #

setAs("bayesGAMfit", "stanfit",
      function(from, to) {
        new("stanfit", from@results)
      })

#' Additional plotting for MCMC visualization and diagnostics. 
#'
#' Marginal response smooth plot functions for parametric and nonparametric associations. 
#'
#' @seealso \code{\link{mcmc_plots}}
#'
#' @param x an object of class \code{hmclearn}, usually a result of a call to \code{mh} or \code{hmc}
#' @param y unused
#' @param applylink logical to indicate whether the inverse link function should be applied to the plots
#' @param ... optional additional arguments to pass to the \code{ggplot2}
#' @references H. Wickham. \emph{ggplot2: Elegant Graphics for Data Analysis}. Springer-Verlag New York, 2016.
#' @export
#' @name plot
#' @return A list of \emph{univariate} and \emph{bivariate} plots generated by plot functions based on \code{ggplot2}
#' @examples
#' f <- bayesGAM(weight ~ np(height), data = women, 
#'               family = gaussian, iter=500, chains = 1)
#' plot(f)
#' 
NULL

#' @rdname plot
setMethod("plot", signature(x="bayesGAMfit", y="missing"),
          function(x, y, applylink=TRUE, ...) {
            smooth(x, applylink=applylink, ...)
          })


#' @rdname plot
setMethod("plot", signature(x="predictPlotObject", y="missing"),
          function(x, ...) {
            object <- x
            ypred_plot(object, ...)
          })

#' @rdname plot
setMethod("plot", signature(x="posteriorPredictObject", y="missing"),
          function(x, ...) {
            object <- x
            ppc_dens_overlay(object, ...)
          })

setGeneric("createPlotData", function(object, ...) {
  standardGeneric("createPlotData")
})

setMethod("createPlotData", signature("bayesGAMfit"),
          function(object, type="standard", applylink=TRUE, nvals=100, ...) {
            
            model <- object@model
            X <- getDesign(model, "X")
            Z <- getDesign(model, "Z")
            y <- model@y
            p <- model@p
            q <- model@q
            has_intercept <- model@has_intercept
            zvars <- model@zvars
            random_intercept <- model@random_intercept
            npargs <- model@npargs
            npterms <- model@npterms
            npdegree <- model@npdegree
            npbasis <- model@basis
            npknots <- model@knots
            betavals <- coefficients(object, params="beta")
            betanms <- object@model@names_beta
            uvals <- coefficients(object, params="u")
            unms <- object@model@names_u
            ynm <- object@model@names_y
            linkname <- object@model@linkname
            names_y <- get_multy_names(object@model@names_y)
            results <- object@results
            mcmcres <- object@mcmcres
            multresponse <- object@model@multresponse
            
            if (type == "standard") {
              
              if (q > 0) {
                Z <- getDesign(model, "Z")
                
                allnms <- c(model@names_y,
                            model@names_beta,
                            model@names_u)
                
                pdata <- data.frame(y = y,
                                    X = X,
                                    Z = Z)
              } else {
                allnms <- c(model@names_y,
                            model@names_beta)
                pdata <- data.frame(y = y,
                                    X = X)
              }
              
              colnames(pdata) <- allnms
            } else if (type == "smooth") {
              smoothobj <- new("smoothPlotObject",
                               betanms = betanms,
                               unms = unms,
                               betavals = betavals,
                               uvals = uvals,
                               npargs = npargs,
                               npbasis = npbasis,
                               xvars = betanms,
                               xvars_static = "",
                               npdegree = npdegree,
                               has_intercept= has_intercept,
                               random_intercept = random_intercept,
                               Xorig = X,
                               knots = npknots,
                               zvars = zvars,
                               linkname = linkname,
                               names_y = names_y,
                               results = results,
                               multresponse =  multresponse, 
                               mcmcres = mcmcres)
              
              pdata <-  createSmoothPdata(smoothobj, applylink=applylink, ...)
            }
    
            return(pdata)
          })


# fitted v actual plot
setGeneric("fitvact", function(object, ...) {
  standardGeneric("fitvact")
})

setMethod("fitvact", signature("bayesGAMfit"),
          function(object, ...) {
            fitvals <- fitted(object, ...)
            actvals <- object@model@y

            # TODO:  replace with ggplot2
            plot(fitvals, actvals, xlab="predicted", ylab="actual")
            abline(a=0,b=1)
          })


# smoothed response plot by independent variables
setGeneric("smooth", function(object, applylink=TRUE, ...) {
  standardGeneric("smooth")
})

setMethod("smooth", signature("bayesGAMfit"),
          function(object, applylink=TRUE, ...) {

            pdata <- object@pdata
            if (length(pdata) == 0) {
              pdata <- createPlotData(object, type="smooth", applylink=applylink, ...)
            }
            
            if (is.null(pdata)) {
              stop("plot not available for this model")
            }
            
            colnames(pdata) <- make.names(colnames(pdata))

            # split dataframe by bivariate
            pdata$isbivariate <- factor(ifelse(pdata$type == 'bivariate',
                                               'bivariate', 'univariate'),
                                        levels=c("univariate", "bivariate"), ordered=T)

            multresponse <- object@model@multresponse
            if (multresponse) {
              ynm_all <- trimws(get_multy_names(object@model@names_y))
              colnames(pdata)[1:length(ynm_all)] <- ynm_all
              ylower_all <- colnames(pdata)[grepl("ylower.*", colnames(pdata))]
              yupper_all <- colnames(pdata)[grepl("yupper.*", colnames(pdata))]
            } else {
              ynm_all <- make.names(object@model@names_y)
              ylower_all <- "ylower"
              yupper_all <- "yupper"
            }

            res <- lapply(split(pdata, pdata$isbivariate), function(xx) {
              if (nrow(xx) > 0) {

                if (xx$isbivariate[1] == 'univariate') {

                  # include multresponse option
                  p1all <- lapply(seq_along(ynm_all), function(yy) {
                    p1 <- ggplot2::ggplot(data=xx,
                                          ggplot2::aes_string(x="xplotvals",
                                                              y=ynm_all[yy]))
                    p1 <- p1 + ggplot2::xlab("")

                    p1 <- p1 + ggplot2::geom_ribbon(data=xx,
                                                    ggplot2::aes_string(ymin=ylower_all[[yy]],
                                                                        ymax=yupper_all[[yy]]),
                                                    alpha=0.2)
                    p1 <- p1 + ggplot2::geom_line(data=xx,
                                                  colour="blue",
                                                  inherit.aes=TRUE)


                    p1 <- p1 + ggplot2::facet_wrap( ~ grouping, scales="free_x")
                    p1 <- p1 + ggplot2::theme_bw()
                  })

                  p1all

                  # bivariate plot
                } else {
                  p2all <- lapply(seq_along(ynm_all), function(yy) {
                    
                    # get x and y labels
                    gnames <- unique(xx$grouping)
                    gkeep <- gnames[grepl("~", gnames)]
                    xyparam <- trimws(unlist(strsplit(gkeep, "~")))
                    labels.list <- as.list(xyparam)
                    
                    # for facets
                    xx$grouping[xx$grouping == gkeep] <- ynm_all[yy]
                    
                    
                    # carveout portion of plot not in data
                    # lapply to apply logit to each facet
                    Xorig <- object@model@X[, rev(xyparam)]
                    xx <- lapply(split(xx, xx$grouping), function(zz) {
                      ch <- geometry::convhulln(Xorig)
                      ptsinhull <- geometry::inhulln(ch, as.matrix(zz[, c("xplotvals", "yplotvals")]))
                      zz[!ptsinhull, which(colnames(xx) == ynm_all[yy])] <- NA
                      zz
                    })
                    xx <- as.data.frame(do.call(rbind, xx))
                    
           
                    
                    p2 <- ggplot2::ggplot(data=xx,
                                          ggplot2::aes_string(x="xplotvals",
                                                              y="yplotvals",
                                                              z=ynm_all[yy]))
                    p2 <- p2 + ggplot2::geom_raster(data=xx,
                                                    ggplot2::aes_string(x="xplotvals",
                                                                        y="yplotvals",
                                                                        fill = ynm_all[yy])) +
                                                    ggplot2::geom_contour(colour="white", na.rm=TRUE)

                    p2 <- p2 + ggplot2::xlab(labels.list[[2]]) + ggplot2::ylab(labels.list[[1]])
                    p2 <- p2 + ggplot2::scale_fill_distiller(palette="Spectral", na.value="transparent")
                    p2 <- p2 + ggplot2::facet_wrap( ~ grouping, scales="free_x")
                    p2 <- p2 + ggplot2::guides(fill=ggplot2::guide_legend(title=NULL))
                    p2 <- p2 + ggplot2::theme_bw()
                  })

                  p2all
                }
              }
            })
            
            res
    })

# initialize smooth plot obj
setGeneric("createSmoothPlotObject", function(object, ...) {
  standardGeneric("createSmoothPlotObject")
})

setMethod("initialize", "smoothPlotObject",
          function(.Object, betanms, unms, betavals, uvals, npargs, npbasis, npdegree,
                   xvars, xvars_static, zvals,
                   xvars_npargs, xvars_np, xvars_basis, has_intercept, random_intercept,
                   Xorig, knots, zvars, linkname, names_y, results, multresponse,
                   mcmcres, pdata, ...) {


            if (has_intercept) {
              xvars <- xvars[-1]
            }

            # only univariate np terms
            check_bivariate <- sapply(npargs, length) == 2
            npargs_v <- npargs[!check_bivariate]
            npargs_v <- unlist(npargs_v)
            xvars_bivariate <- npargs[check_bivariate]

            which_truncpoly <- which(npbasis == "trunc.poly")

            npargs_v_truncpoly <- character()
            npargs_v_exclude <- character()
            if (length(npargs_v) > 0 & length(which_truncpoly) > 0) {

              npargs_v_truncpoly <- paste0(npargs_v, 1)
              deg <- npdegree[which_truncpoly]

              npargs_v_exclude <- mapply(function(nms, maxdegree) {
                paste0(nms, 2:maxdegree)
              }, nms = as.list(npargs_v), maxdegre=as.list(deg), SIMPLIFY=FALSE)
              npargs_v_exclude <- unlist(npargs_v_exclude)
            }

            xvars_static <- xvars[!xvars %in% npargs_v &
                                    !xvars %in% npargs_v_truncpoly &
                                    !xvars %in% npargs_v_exclude &
                                    !xvars %in% unlist(xvars_bivariate)]

            xvars_np <- npargs_v
            
            # check trunc poly
            ind_trunc_poly <- npbasis == "trunc.poly"
            ind_trunc_poly_univ <- ind_trunc_poly[!check_bivariate]
            xvars_np[ind_trunc_poly_univ] <- 
              paste0(xvars_np[ind_trunc_poly_univ], 1)
            
            if (length(npargs_v) == 0) {
              xvars_np <- character()
            }
            
            .Object@betanms <- betanms
            .Object@unms <- unms
            .Object@betavals <- betavals
            .Object@uvals <- uvals
            .Object@xvars <- xvars
            .Object@xvars_static <- xvars_static
            .Object@xvars_npargs <- npargs
            .Object@xvars_np <- xvars_np
            .Object@xvars_bivariate <- xvars_bivariate
            .Object@xvars_basis <- npbasis
            .Object@npdegree <- as.integer(npdegree)
            .Object@has_intercept <- has_intercept
            .Object@random_intercept <- random_intercept
            .Object@Xorig <- Xorig
            .Object@knots <- knots
            .Object@zvars <- zvars
            .Object@linkname <- linkname
            .Object@names_y <- names_y
            .Object@results <- results
            .Object@multresponse <- multresponse
            .Object@mcmcres <- mcmcres
            return(.Object)
          })



setGeneric("createSmoothPdata", function(smoothobj, ...) {
  standardGeneric("createSmoothPdata")
})


setMethod("createSmoothPdata", signature("smoothPlotObject"),
          function(smoothobj, nvals=100, applylink=TRUE, rg=c(-1.96, 1.96),
                   probs = c(0.025, 0.975),
                   interval = "simulation", ...) {
            
            betanms <- smoothobj@betanms
            unms <- smoothobj@unms
            betavals <- smoothobj@betavals
            uvals <- smoothobj@uvals
            xvars <- smoothobj@xvars
            xvars_static <- smoothobj@xvars_static
            xvars_npargs <- smoothobj@xvars_npargs
            xvars_np <- smoothobj@xvars_np
            xvars_basis <- smoothobj@xvars_basis
            xvars_bivariate <- smoothobj@xvars_bivariate
            npdegree <- smoothobj@npdegree
            has_intercept <- smoothobj@has_intercept
            random_intercept <- smoothobj@random_intercept
            Xorig <- smoothobj@Xorig
            knots <- smoothobj@knots
            zvars <- smoothobj@zvars
            linkname <- smoothobj@linkname
            names_y <- smoothobj@names_y
            results <- smoothobj@results
            multresponse <- smoothobj@multresponse
            
            if (length(smoothobj@mcmcres) > 0) {
              mcmcres <- smoothobj@mcmcres  
            } else {
              mcmcres <- as.matrix(smoothobj@results)
            }
            
            Xmeans <- create_xmeans(xvars_static, Xorig, nvals, xvars_np, xvars_npargs,
                                    xvars_bivariate, xvars_basis,
                                    npdegree, has_intercept, betanms)

            allvars <- as.list(c(xvars_static, xvars_np))
            if (length(xvars_bivariate) > 0) {
              allvars <- c(allvars, xvars_bivariate)
            }
            
            res <- lapply(allvars, create_single_smooth_data,
                          Xorig=Xorig,
                          Xmeans=Xmeans,
                          xvars_static=xvars_static,
                          xvars_np=xvars_np,
                          xvars_npargs=xvars_npargs,
                          xvars_basis=xvars_basis,
                          xvars_bivariate=xvars_bivariate,
                          npdegree=npdegree,
                          random_intercept=random_intercept,
                          nvals=nvals,
                          knots=knots,
                          zvars=zvars,
                          betanms=betanms,
                          unms=unms,
                          betavals=betavals,
                          uvals=uvals,
                          linkname=linkname,
                          applylink=applylink,
                          rg=rg,
                          probs=probs,
                          names_y = trimws(names_y), 
                          simresults = results,
                          interval = interval,
                          multresponse = multresponse, 
                          mcmcres = mcmcres)

            pdata <- do.call(rbind, res)
            return(pdata)
          })


# ---------------- multivariate correlation matrix plot  ---------------------------- #
setGeneric("mvcorrplot", function(object, ...) {
  standardGeneric("mvcorrplot")
})

#' Multivariate response correlation plot for \code{bayesGAMfit} objects
#'
#' Creates a correlation plot of the multivariate responses based on \code{corrplot} 
#'
#' @param object model object of class \code{bayesGAMfit}
#' @param ... Additional parameters passed to \code{corrplot.mixed}
#' @references Taiyun Wei and Viliam Simko (2017). R package \emph{corrplot}: Visualization of a Correlation Matrix (Version 0.84).
#' @return \code{corrplot} object
#' @name mvcorrplot
#' @export
#' @examples
#' 
#' require(MASS)
#' sig <- matrix(c(1, 0.5, 0.5, 1), ncol=2)
#' set.seed(123)
#' Y <- mvrnorm(50, mu=c(-2, 2), Sigma=sig)
#' dat <- data.frame(id = rep(1:5, each=10),
#'                   y1 = Y[, 1], 
#'                   y2 = Y[, 2])
#' 
#' f <- bayesGAM(cbind(y1, y2) ~ 1, random = ~factor(id), 
#'               data=dat, 
#'               a = normal(c(0, 5)), 
#'               chains = 1, iter = 500)
#' mvcorrplot(f)
#' 
NULL

#' @rdname mvcorrplot
#' @export
setMethod("mvcorrplot", "bayesGAMfit",
          function(object, ...) {
            if (!object@model@multresponse | !object@model@random_intercept) {
              stop("mvcorrplot is valid for multivariate response models with random intercepts only")
            }
            # get correlations from simulation
            allcoef <- coefficients(object)
            nms <- names(allcoef)
            sigucoef <- allcoef[grepl("^sigma_u_correlation", nms)]
            
            # sort coefficients column major
            nmsplit <- strsplit(names(sigucoef), split="\\[|\\,|\\]")
            nmdf <- data.frame(do.call(rbind, nmsplit)[, 2:3])
            nmdf$corr <- sigucoef
            colnames(nmdf) <- c("var1", "var2", "corr")
            nmdf <- nmdf[order(nmdf$var2, nmdf$var1), ]
            
            # create matrix
            sigu_correlation <- matrix(nmdf$corr, ncol=sqrt(nrow(nmdf)))

            # plot matrix
            corrplot::corrplot.mixed(sigu_correlation, ...)
            
          })


# get MCMC samples for marginal plotting
#' @rdname getSamples
#' @export
setGeneric("getSamples", function(object, ...) {
  standardGeneric("getSamples")
})

#' Extract the MCMC samples from an object of type \code{bayesGAMfit}
#'
#' Returns an array of the posterior simulation from \code{Stan}. Optionally, may return a subsample from the full MCMC simulation.
#'
#' @param object model object of class \code{bayesGAMfit}
#' @param nsamp Optional number of samples to return
#' @param results Matrix of HMC posterior samples
#' @param seednum Optional integer for seed number when selecting a random sample
#' @param ... Additional parameters passed to \code{corrplot.mixed}
#' @return array of the posterior simulation, or subsample of the array
#' @name getSamples
#' @section 
#' 
#' @export
#' @examples
#' require(stats); require(graphics)
#' f <- bayesGAM(weight ~ np(height), data = women, family = gaussian, 
#'               iter = 500, chains = 1)
#' allres <- getSamples(f)             
NULL

#' @rdname getSamples
#' @export
setMethod("getSamples", signature("bayesGAMfit"),
          function(object, nsamp=NULL, seednum=NULL, ...) {
            
            famnum <- object@model@famnum

            results <- object@results

            if (is.numeric(seednum)) {
              set.seed(seednum)
            }

            # get samples from rstan object
            if (length(object@mcmcres) > 0) {
              mcmcres <- object@mcmcres  
            } else {
              mcmcres <- as.matrix(results)
            }
            
            nsim <- nrow(mcmcres)

            if (is.null(nsamp)) {
              index <- 1:nrow(mcmcres)
            } else {
              index <- sample(1:nsim, nsamp, replace=TRUE)              
            }

            # condition on famnum
            xnms <- colnames(mcmcres)[grepl("^beta", colnames(mcmcres), ignore.case = TRUE)]
            has_random_effects <- grepl("^u", colnames(mcmcres), ignore.case = TRUE)
            
            if (any(has_random_effects)) {
              znms <- colnames(mcmcres)[has_random_effects]
              allnms <- c(xnms, znms)
            } else {
              allnms <- xnms
            }
            
            # epsilon for Gaussian distribution
            if (famnum == 1) {
              has_eps <- grepl("^eps", colnames(mcmcres), ignore.case=TRUE)
              epsnms <- colnames(mcmcres)[has_eps]
              allnms <- c(allnms, epsnms)
            }

            mcmcres[index, allnms]

          })


#' @rdname getSamples
setMethod("getSamples", signature("stanfit"),
          function(object, nsamp=1000, seednum=NULL, results=NULL, ...) {
            
            # results <- object
            
            if (is.numeric(seednum)) {
              set.seed(seednum)
            }
            
            # get samples from rstan object
            if (is.null(results) | length(results) == 0) {
              mcmcres <- as.matrix(object)  
            } else {
              mcmcres <- results
            }
            
            nsim <- nrow(mcmcres)
            
            index <- sample(1:nsim, nsamp, replace=T)
            
            # condition on famnum
            xnms <- colnames(mcmcres)[grepl("^beta", colnames(mcmcres), ignore.case = T)]
            
            has_random_effects <- grepl("^u", colnames(mcmcres), ignore.case = T)
            
            if (any(has_random_effects)) {
              znms <- colnames(mcmcres)[has_random_effects]
              allnms <- c(xnms, znms)
            } else {
              allnms <- xnms
            }
            
            mcmcres[index, allnms]
            
          })


#' @rdname getSamples
#' @export
setMethod("getSamples", signature("glmModel"),
 function(object, nsamp=NULL, seednum=NULL, results=NULL, ...) {
  
  famnum <- object@famnum
  
  # results <- object@results
  
  if (is.numeric(seednum)) {
    set.seed(seednum)
  }
  
  # get samples from rstan object
  mcmcres <- results
  nsim <- nrow(mcmcres)
  
  if (is.null(nsamp)) {
    index <- 1:nrow(mcmcres)
  } else {
    index <- sample(1:nsim, nsamp, replace=TRUE)              
  }
  
  # condition on famnum
  xnms <- colnames(mcmcres)[grepl("^beta", colnames(mcmcres), ignore.case = TRUE)]
  has_random_effects <- grepl("^u", colnames(mcmcres), ignore.case = TRUE)
  
  if (any(has_random_effects)) {
    znms <- colnames(mcmcres)[has_random_effects]
    allnms <- c(xnms, znms)
  } else {
    allnms <- xnms
  }
  
  # epsilon for Gaussian distribution
  if (famnum == 1) {
    has_eps <- grepl("^eps", colnames(mcmcres), ignore.case=TRUE)
    epsnms <- colnames(mcmcres)[has_eps]
    allnms <- c(allnms, epsnms)
  }
  
  mcmcres[index, allnms]
})


# ---------------- show bayesGAMfit object  ---------------------------- #

#' @rdname bayesGAMfit
setMethod("show", "bayesGAMfit",
          function(object) {
            stanobj <- object@results
            multresponse <- object@model@multresponse
            xbetanms <- object@model@names_beta

            # rename beta param
            if (multresponse) {
              ynm_all <- trimws(get_multy_names(object@model@names_y))
              beta_param_nms <- names(stanobj)[grepl("^beta", names(stanobj))]
              xnums <- get_multresponse_xnums(beta_param_nms)
              xnum_uniq <- sort(unique(xnums))

              beta_param_nms <- lapply(seq_along(xbetanms), function(yy) {
                beta_nms_temp <- beta_param_nms[xnums == yy]
                beta_nms_temp <- gsub("beta", paste("beta", xbetanms[yy], sep="_"),
                                      beta_nms_temp)
                beta_nms_temp
              })

              beta_param_nms <- unlist(beta_param_nms)
              beta_param_nms <- gsub(pattern="\\,*.\\]", "]", beta_param_nms)

              names(stanobj)[grepl("^beta", names(stanobj))] <- beta_param_nms

            } else {
              names(stanobj)[grepl("^beta", names(stanobj))] <-
                object@model@names_beta
            }

            default_pars <- stanobj@sim$pars_oi
            pars <- default_pars[default_pars %in% c("beta", "eps", "lambda", object@model@names_beta, "a")]
            print(x=stanobj, pars=pars)
          })

#' Summarizing Model Fits from \code{bayesGAM}
#'
#' summary method for class \code{bayesGAMfit}
#'
#' @param object an object of class \code{hmclearn}, usually a result of a call to \code{mh} or \code{hmc}
#' @returns Returns a matrix with posterior quantiles and the posterior scale reduction factor statistic for each parameter.
#' @references Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.1.
#' @name summary
#' @export
#' @examples
#' f <- bayesGAM(weight ~ np(height), data = women, 
#'               family = gaussian, iter=500, chains = 1)
#'
#' summary(f)
NULL

#' @rdname summary
#' @export
setMethod("summary", "bayesGAMfit",
          function(object) {
            show(object)
          })


# ---------------- bayesplot functionality  ---------------------------- #

#' Plotting for MCMC visualization and diagnostics provided by \code{bayesplot} package
#'
#' Plots of Rhat statistics, ratios of effective sample size to total sample
#' size, and autocorrelation of MCMC draws.
#'
#' @name mcmc_plots
#'
#' @param object an object of class \code{bayesGAMfit}
#' @param regex_pars character vector of regular expressions of variable names to plot
#' @param ... optional additional arguments to pass to the \code{bayesplot} functions
#'
#' @section Plot Descriptions from the \code{bayesplot} package documentation:
#' \itemize{
#'   \item{`mcmc_hist(object, ...)`}{
#'    Default plot called by `plot` function.  Histograms of posterior draws with all chains merged.
#'   }
#'   \item{`mcmc_dens(object, ...)`}{
#'    Kernel density plots of posterior draws with all chains merged.
#'   }
#'   \item{`mcmc_hist_by_chain(object, ...)`}{
#'    Histograms of posterior draws with chains separated via faceting.
#'   }
#'   \item{`mcmc_dens_overlay(object, ...)`}{
#'    Kernel density plots of posterior draws with chains separated but
#'    overlaid on a single plot.
#'   }
#'   \item{`mcmc_violin(object, ...)`}{
#'    The density estimate of each chain is plotted as a violin with
#'    horizontal lines at notable quantiles.
#'   }
#'   \item{`mcmc_dens_chains(object, ...)`}{
#'    Ridgeline kernel density plots of posterior draws with chains separated
#'    but overlaid on a single plot. In `mcmc_dens_overlay()` parameters
#'    appear in separate facets; in `mcmc_dens_chains()` they appear in the
#'    same panel and can overlap vertically.
#'   }
#'   \item{`mcmc_intervals(object, ...)`}{
#'    Plots of uncertainty intervals computed from posterior draws with all
#'    chains merged.
#'   }
#'   \item{`mcmc_areas(object, ...)`}{
#'    Density plots computed from posterior draws with all chains merged,
#'    with uncertainty intervals shown as shaded areas under the curves.
#'   }
#'   \item{`mcmc_scatter(object, ...)`}{
#'    Bivariate scatterplot of posterior draws. If using a very large number of
#'    posterior draws then `mcmc_hex()` may be preferable to avoid
#'    overplotting.
#'   }
#'   \item{`mcmc_hex(object, ...)`}{
#'    Hexagonal heatmap of 2-D bin counts. This plot is useful in cases where
#'    the posterior sample size is large enough that `mcmc_scatter()` suffers
#'    from overplotting.
#'   }
#'   \item{`mcmc_pairs(object, ...)`}{
#'    A square plot matrix with univariate marginal distributions along the
#'    diagonal (as histograms or kernel density plots) and bivariate
#'    distributions off the diagonal (as scatterplots or hex heatmaps).
#'
#'    For the off-diagonal plots, the default is to split the chains so that
#'    (roughly) half are displayed above the diagonal and half are below (all
#'    chains are always merged together for the plots along the diagonal). Other
#'    possibilities are available by setting the `condition` argument.
#'   }
#' \item{`mcmc_rhat(object, ...)`, `mcmc_rhat_hist(object, ...)`}{
#'   Rhat values as either points or a histogram. Values are colored using
#'   different shades (lighter is better). The chosen thresholds are somewhat
#'   arbitrary, but can be useful guidelines in practice.
#'   * _light_: below 1.05 (good)
#'   * _mid_: between 1.05 and 1.1 (ok)
#'   * _dark_: above 1.1 (too high)
#'  }
#'  \item{`mcmc_neff(object, ...)`, `mcmc_neff_hist(object, ...)`}{
#'   Ratios of effective sample size to total sample size as either points or a
#'   histogram. Values are colored using different shades (lighter is better).
#'   The chosen thresholds are somewhat arbitrary, but can be useful guidelines
#'   in practice.
#'   * _light_: between 0.5 and 1 (high)
#'   * _mid_: between 0.1 and 0.5 (good)
#'   * _dark_: below 0.1 (low)
#'  }
#'  \item{`mcmc_acf(object, ...)`, `mcmc_acf_bar(object, ...)`}{
#'   Grid of autocorrelation plots by chain and parameter. The `lags` argument
#'   gives the maximum number of lags at which to calculate the autocorrelation
#'   function. `mcmc_acf()` is a line plot whereas `mcmc_acf_bar()` is a
#'   barplot.
#'  }
#' }
#' @return These functions call various plotting functions from the \code{bayesplot} package, which returns a list including \code{ggplot2} objects.
#' @references Gabry, Jonah and Mahr, Tristan (2019).  \emph{bayesplot:  Plotting for Bayesian Models}.  \url{https://mc-stan.org/bayesplot/}
#' @references Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A (2019).  \emph{Visualization in Bayesian Workflow}.  Journal of the Royal Statistical Society: Series A. Vol 182.  Issue 2.  p.389-402.
#' @references Gelman, A. and Rubin, D. (1992) \emph{Inference from Iterative Simulation Using Multiple Sequences}.  Statistical Science 7(4) 457-472.
#' @references Gelman, A., et. al. (2013) \emph{Bayesian Data Analysis}.  Chapman and Hall/CRC.
NULL

#' @rdname mcmc_plots
#' @export
#' @examples
#' f <- bayesGAM(weight ~ np(height), data = women, 
#'               family = gaussian, iter=1000, chains = 1)
#' mcmc_trace(f)
setGeneric("mcmc_intervals", function(object, ...) {
  standardGeneric("mcmc_intervals")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_intervals", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {
            stanobj <- set_varnms(object)
            bayesplot::mcmc_intervals(stanobj, regex_pars=regex_pars, ...)
          })

#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_areas", function(object, ...) {
  standardGeneric("mcmc_areas")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_areas", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {
            stanobj <- set_varnms(object)
            bayesplot::mcmc_areas(stanobj, regex_pars=regex_pars, ...)
          })

#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_hist", function(object, ...) {
  standardGeneric("mcmc_hist")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_hist", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {
            stanobj <- set_varnms(object)
            bayesplot::mcmc_hist(stanobj, regex_pars=regex_pars, ...)
          })

#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_hist_by_chain", function(object, ...) {
  standardGeneric("mcmc_hist_by_chain")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_hist_by_chain", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {
            stanobj <- set_varnms(object)
            bayesplot::mcmc_hist_by_chain(stanobj, regex_pars=regex_pars, ...)
          })

#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_dens", function(object, ...) {
  standardGeneric("mcmc_dens")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_dens", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {
            stanobj <- set_varnms(object)
            bayesplot::mcmc_dens(stanobj, regex_pars=regex_pars, ...)
          })

#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_scatter", function(object, ...) {
  standardGeneric("mcmc_scatter")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_scatter", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {
            stanobj <- set_varnms(object)
            bayesplot::mcmc_scatter(stanobj, regex_pars=regex_pars, ...)
          })

#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_hex", function(object, ...) {
  standardGeneric("mcmc_hex")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_hex", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {
            stanobj <- set_varnms(object)
            bayesplot::mcmc_hex(stanobj, regex_pars=regex_pars, ...)
          })

#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_pairs", function(object, ...) {
  standardGeneric("mcmc_pairs")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_pairs", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {
            stanobj <- set_varnms(object)
            bayesplot::mcmc_pairs(stanobj, regex_pars=regex_pars, ...)
          })

#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_acf", function(object, ...) {
  standardGeneric("mcmc_acf")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_acf", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {
            stanobj <- set_varnms(object)
            bayesplot::mcmc_acf(stanobj, regex_pars=regex_pars, ...)
          })

#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_acf_bar", function(object, ...) {
  standardGeneric("mcmc_acf_bar")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_acf_bar", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {
            stanobj <- set_varnms(object)
            bayesplot::mcmc_acf_bar(stanobj, regex_pars=regex_pars, ...)
          })

#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_trace", function(object, ...) {
  standardGeneric("mcmc_trace")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_trace", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {
            stanobj <- set_varnms(object)
            
            # bayesplot::mcmc_trace(object@results, regex_pars=regex_pars, ...)
            bayesplot::mcmc_trace(stanobj, regex_pars=regex_pars, ...)
          })

#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_rhat", function(object, ...) {
  standardGeneric("mcmc_rhat")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_rhat", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {
            stanobj <- set_varnms(object)

            rhat_vals <- bayesplot::rhat(stanobj)
            nms <- names(rhat_vals)
            rhat_vals <- rhat_vals[which(grepl(paste(regex_pars, collapse="|"), nms))]

            bayesplot::mcmc_rhat(rhat_vals, ...)
          })

#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_rhat_hist", function(object, ...) {
  standardGeneric("mcmc_rhat_hist")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_rhat_hist", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {

            stanobj <- set_varnms(object)
            rhat_vals <- bayesplot::rhat(stanobj)
            nms <- names(rhat_vals)
            rhat_vals <- rhat_vals[which(grepl(paste(regex_pars, collapse="|"), nms))]

            bayesplot::mcmc_rhat_hist(rhat_vals, ...)
         })

#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_rhat_data", function(object, ...) {
  standardGeneric("mcmc_rhat_data")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_rhat_data", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {

            stanobj <- set_varnms(object)
            rhat_vals <- bayesplot::rhat(stanobj)
            nms <- names(rhat_vals)
            rhat_vals <- rhat_vals[which(grepl(paste(regex_pars, collapse="|"), nms))]

            bayesplot::mcmc_rhat_data(rhat_vals, ...)
          })

#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_neff", function(object, ...) {
  standardGeneric("mcmc_neff")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_neff", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {

            stanobj <- set_varnms(object)
            neffr <- bayesplot::neff_ratio(stanobj)
            nms <- names(neffr)
            neffr <- neffr[which(grepl(paste(regex_pars, collapse="|"), nms))]

            bayesplot::mcmc_neff(neffr, ...)
         })


#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_neff_hist", function(object, ...) {
  standardGeneric("mcmc_neff_hist")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_neff_hist", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {

            stanobj <- set_varnms(object)
            neffr <- bayesplot::neff_ratio(stanobj)
            nms <- names(neffr)
            neffr <- neffr[which(grepl(paste(regex_pars, collapse="|"), nms))]

            bayesplot::mcmc_neff_hist(neffr, ...)
          })

#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_neff_data", function(object, ...) {
  standardGeneric("mcmc_neff_data")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_neff_data", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {

            stanobj <- set_varnms(object)
            neffr <- bayesplot::neff_ratio(stanobj)
            nms <- names(neffr)
            neffr <- neffr[which(grepl(paste(regex_pars, collapse="|"), nms))]

            bayesplot::mcmc_neff_data(neffr, ...)
          })

#' @export
#' @rdname mcmc_plots
setGeneric("mcmc_violin", function(object, ...) {
  standardGeneric("mcmc_violin")
})

#' @export
#' @rdname mcmc_plots
setMethod("mcmc_violin", "bayesGAMfit",
          function(object,
                   regex_pars = c("^beta", "^lambda", "^eps", "^a", "^sigma_u_correlation"),
                   ...) {
            
            stanobj <- set_varnms(object)
            bayesplot::mcmc_violin(stanobj, regex_pars=regex_pars, ...)
          })


#' Posterior predictive samples from models fit by \code{bayesGAM}
#' 
#' Draw from the posterior predictive distribution
#' 
#' @param object Object of type \code{bayesGAMfit} generated from \code{bayesGAM}.
#' @param draws An integer indicating the number of draws to return. The default and maximum number of draws is the size of the posterior sample.
#' @param results Matrix of HMC posterior samples
#' @return a list of \emph{D} by \emph{N} matrices, where \emph{D} is the number of draws from the posterior predictive distribution and \emph{N} is the number of data points being predicted per draw.
#' @export
#' @references Goodrich B, Gabry J, Ali I & Brilleman S. (2020). rstanarm: Bayesian applied regression modeling via Stan. R package version 2.19.3 https://mc-stan.org/rstanarm.
#' @references Jonah Gabry, Ben Goodrich and Martin Lysy (2020). rstantools: Tools for Developing R Packages Interfacing with 'Stan'. https://mc-stan.org/rstantools/, https://discourse.mc-stan.org/.
#' @name posterior_predict
NULL

#' @rdname posterior_predict
#' @param ... Additional arguments for \code{postrior_predict}
#' @examples
#' f <- bayesGAM(weight ~ np(height), data = women, 
#'               family = gaussian, iter=1000, chains = 1)
#' res <- posterior_predict(f, draws=100)
#' @export
setGeneric("posterior_predict", function(object, ...) {
  standardGeneric("posterior_predict")
})

#' @rdname posterior_predict
#' @export
setMethod("posterior_predict", signature(object="bayesGAMfit"),
          function(object, draws=NULL, ...) {
            
      famnum <- object@model@famnum
      linknum <- object@model@linknum
      offset <- object@offset
      mcmcres <- object@mcmcres
      q <- object@model@q
      r <- object@model@r
      X <- getDesign(object, "X")
      Z <- getDesign(object, "Z")
      if (q > 0) {
        W <- cbind(X, Z)  
      } else {
        W <- X
      }
      
      multresponse <- object@model@multresponse
            
      # retrieve samples from posterior
      thetasamp <- getSamples(object)
      nsamp <- nrow(thetasamp)
      if (is.null(draws) | draws > nsamp) {
        draws <- nsamp
      }
      useallsamples <- draws == nsamp
      
      # get sample indices
      ivals <- 1:nsamp
      if (useallsamples) {
        index <- ivals
      } else {
        index <- sample(ivals, size=draws)
      }

      if (famnum == 1) {
        betaunms <- grepl("^beta|^u", colnames(thetasamp))
        thetasamp_betau <- thetasamp[, betaunms]
        thetasamp_eps <- thetasamp[, !betaunms]
        yrep <- pp_sim(index, r, W, famnum, linknum, offset,  
                           thetasamp_betau, thetasamp_eps, multresponse)
      } else if (famnum %in% c(2, 3)) {
        thetasamp_betau <- thetasamp
        thetasamp_eps <- NULL
        yrep <- pp_sim(index, r, W, famnum, linknum, offset,  
                       thetasamp_betau, thetasamp_eps, multresponse)
      }
      
      mdl <- object@model
      
      obj <- new("posteriorPredictObject", 
                 model=mdl, 
                 pp = yrep, 
                 thetasamp = thetasamp)
      
    return(obj)
})

#' @rdname posterior_predict
#' @export
setMethod("posterior_predict", signature(object="glmModel"),
     function(object, draws=NULL, results=NULL, ...) {
  
  famnum <- object@famnum
  linknum <- object@linknum
  offset <- object@offset
  q <- object@q
  r <- object@r
  X <- object@X
  Z <- object@Z
  
  if (q > 0) {
    W <- cbind(X, Z)  
  } else {
    W <- X
  }

  multresponse <- object@multresponse
  
  # retrieve samples from posterior
  thetasamp <- getSamples(object, results=results)
  nsamp <- nrow(thetasamp)
  if (is.null(draws) | draws > nsamp) {
    draws <- nsamp
  }
  useallsamples <- draws == nsamp
  
  # get sample indices
  ivals <- 1:nsamp
  if (useallsamples) {
    index <- ivals
  } else {
    index <- sample(ivals, size=draws)
  }
  
  if (famnum == 1) {
    betaunms <- grepl("^beta|^u", colnames(thetasamp))
    thetasamp_betau <- thetasamp[, betaunms]
    thetasamp_eps <- thetasamp[, !betaunms]
    yrep <- pp_sim(index, r, W, famnum, linknum, offset,  
                   thetasamp_betau, thetasamp_eps, multresponse)
  } else if (famnum %in% c(2, 3)) {
    thetasamp_betau <- thetasamp
    thetasamp_eps <- NULL
    yrep <- pp_sim(index, r, W, famnum, linknum, offset,  
                   thetasamp_betau, thetasamp_eps, multresponse)
  }
  
  # for predict method
  names(yrep) <- gsub("yrep", "ypred", names(yrep))
  
  return(yrep)
})


#' Plotting for MCMC visualization and diagnostics provided by \code{bayesplot} package
#'
#' Plots of Rhat statistics, ratios of effective sample size to total sample
#' size, and autocorrelation of MCMC draws.
#'
#' @name ppc_plots
#'
#' @param object an object of class \code{bayesGAMfit}
#' @param draws An integer indicating the number of draws to return. The default and maximum number of draws is the size of the posterior sample.
#' @param ... optional additional arguments to pass to the \code{bayesplot} functions
#'
#' @section Plot Descriptions from the \code{bayesplot} package documentation:
#' \itemize{
#'   \item{`ppc_hist(object, draws=NULL, ...)`}{
#'    A separate histogram estimate is displayed for y and each dataset (row) in yrep. For these plots yrep should therefore contain only a small number of rows. 
#'   }
#'   \item{`ppc_boxplot(object, draws=NULL, ...)`}{
#'    A separate box and whiskers plot is displayed for y and each dataset (row) in yrep. For these plots yrep should therefore contain only a small number of rows. 
#'   }
#'   \item{`ppc_freqpoly(object, draws=NULL, ...)`}{
#'    A separate shaded frequency polygon is displayed for y and each dataset (row) in yrep. For these plots yrep should therefore contain only a small number of rows. 
#'   }
#'   \item{`ppc_dens(object, draws=NULL, ...)`}{
#'    A separate smoothed kernel density estimate is displayed for y and each dataset (row) in yrep. For these plots yrep should therefore contain only a small number of rows. 
#'   }
#'   \item{`ppc_dens_overlay(object, draws=NULL, ...)`}{
#'    Kernel density estimates of each dataset (row) in \code{yrep} are overlaid, with the distribution of \code{y} itself on top (and in a darker shade).
#'   }
#'   \item{`ppc_ecdf_overlay(object, draws=NULL, ...)`}{
#'    Empirical CDF estimates of each dataset (row) in \code{yrep} are overlaid, with the distribution of \code{y} itself on top (and in a darker shade).
#'   }
#' }
#' @return These functions call various plotting functions from the \code{bayesplot} package, which returns a list including \code{ggplot2} objects.
#' @references Gabry, Jonah and Mahr, Tristan (2019).  \emph{bayesplot:  Plotting for Bayesian Models}.  \url{https://mc-stan.org/bayesplot/}
#' @references Gabry, J., Simpson, D., Vehtari, A., Betancourt, M., and Gelman, A (2019).  \emph{Visualization in Bayesian Workflow}.  Journal of the Royal Statistical Society: Series A. Vol 182.  Issue 2.  p.389-402.
#' @references Gelman, A. and Rubin, D. (1992) \emph{Inference from Iterative Simulation Using Multiple Sequences}.  Statistical Science 7(4) 457-472.
#' @references Gelman, A., et. al. (2013) \emph{Bayesian Data Analysis}.  Chapman and Hall/CRC.
#' @references Gabry, J. , Simpson, D. , Vehtari, A. , Betancourt, M. and Gelman, A. (2019), Visualization in Bayesian workflow. J. R. Stat. Soc. A, 182: 389-402. doi:10.1111/rssa.12378.
NULL

#' @export
#' @examples
#' f <- bayesGAM(weight ~ np(height), data = women, 
#'               family = gaussian, iter=500, chains = 1)
#' ppc_dens(f, draws=2)
#' @rdname ppc_plots
setGeneric("ppc_dens", function(object, ...) {
  standardGeneric("ppc_dens")
})

#' @export
#' @rdname ppc_plots
setMethod("ppc_dens", "bayesGAMfit",
          function(object, draws=NULL, ...) {
            pp <- posterior_predict(object, draws=draws)@pp
            
            y <- object@model@y
            ylst <- unlist(apply(y, 2, list), recursive = FALSE)
            
            pfun <- function(y, yrep, nms, ...) {
              bayesplot::ppc_dens(y, yrep, ...) + 
                ggplot2::ggtitle(nms) 
            }
            
            bp <- mapply(pfun, 
                         y = ylst, 
                         yrep = pp, 
                         nms = as.list(names(pp)), 
                         MoreArgs=list(... = ...), 
                         SIMPLIFY=FALSE)
            bayesplot::bayesplot_grid(plots=bp)
            
          })


#' @export
#' @rdname ppc_plots
setMethod("ppc_dens", "posteriorPredictObject",
          function(object, ...) {
            pp <- object@pp
            
            y <- object@model@y
            ylst <- unlist(apply(y, 2, list), recursive = FALSE)
            
            pfun <- function(y, yrep, nms, ...) {
              bayesplot::ppc_dens(y, yrep, ...) + 
                ggplot2::ggtitle(nms) 
            }
            
            bp <- mapply(pfun, 
                         y = ylst, 
                         yrep = pp, 
                         nms = as.list(names(pp)), 
                         MoreArgs=list(... = ...), 
                         SIMPLIFY=FALSE)
            bayesplot::bayesplot_grid(plots=bp)
            
          })

#' @export
#' @rdname ppc_plots
setGeneric("ppc_dens_overlay", function(object, ...) {
  standardGeneric("ppc_dens_overlay")
})

#' @export
#' @rdname ppc_plots
setMethod("ppc_dens_overlay", "bayesGAMfit",
          function(object, draws=NULL, ...) {
            
            if (is.null(draws)) {
              draws <- pmin(100, nrow(object@model@y))
            }
            
            
            pp <- posterior_predict(object, draws=draws)@pp
            
            y <- object@model@y
            ylst <- unlist(apply(y, 2, list), recursive = FALSE)
            
            pfun <- function(y, yrep, nms, ...) {
              bayesplot::ppc_dens_overlay(y, yrep, ...) + 
                ggplot2::ggtitle(nms) 
            }
            
            bp <- mapply(pfun, 
                          y = ylst, 
                          yrep = pp, 
                          nms = as.list(names(pp)), 
                        MoreArgs=list(... = ...), 
                      SIMPLIFY=FALSE)
            bayesplot::bayesplot_grid(plots=bp)
            
          })

#' @export
#' @rdname ppc_plots
setMethod("ppc_dens_overlay", "posteriorPredictObject",
          function(object, ...) {
            pp <- object@pp
            
            y <- object@model@y
            ylst <- unlist(apply(y, 2, list), recursive = FALSE)
            
            pfun <- function(y, yrep, nms, ...) {
              bayesplot::ppc_dens_overlay(y, yrep, ...) + 
                ggplot2::ggtitle(nms) 
            }
            
            bp <- mapply(pfun, 
                         y = ylst, 
                         yrep = pp, 
                         nms = as.list(names(pp)), 
                         MoreArgs=list(... = ...), 
                         SIMPLIFY=FALSE)
            bayesplot::bayesplot_grid(plots=bp)
            
          })


#' @export
#' @rdname ppc_plots
setGeneric("ppc_hist", function(object, ...) {
  standardGeneric("ppc_hist")
})

#' @export
#' @rdname ppc_plots
setMethod("ppc_hist", "bayesGAMfit",
          function(object, draws=NULL, ...) {
            pp <- posterior_predict(object, draws=draws)@pp
            
            y <- object@model@y
            ylst <- unlist(apply(y, 2, list), recursive = FALSE)
            
            pfun <- function(y, yrep, nms, ...) {
              bayesplot::ppc_hist(y, yrep, ...) + 
                ggplot2::ggtitle(nms) 
            }
            
            bp <- mapply(pfun, 
                         y = ylst, 
                         yrep = pp, 
                         nms = as.list(names(pp)), 
                         MoreArgs=list(... = ...), 
                         SIMPLIFY=FALSE)
            bayesplot::bayesplot_grid(plots=bp)
            
          })

#' @export
#' @rdname ppc_plots
setMethod("ppc_hist", "posteriorPredictObject",
          function(object, ...) {
            pp <- object@pp
            
            y <- object@model@y
            ylst <- unlist(apply(y, 2, list), recursive = FALSE)
            
            pfun <- function(y, yrep, nms, ...) {
              bayesplot::ppc_hist(y, yrep, ...) + 
                ggplot2::ggtitle(nms) 
            }
            
            bp <- mapply(pfun, 
                         y = ylst, 
                         yrep = pp, 
                         nms = as.list(names(pp)), 
                         MoreArgs=list(... = ...), 
                         SIMPLIFY=FALSE)
            bayesplot::bayesplot_grid(plots=bp)
            
          })


#' @export
#' @rdname ppc_plots
setGeneric("ppc_boxplot", function(object, ...) {
  standardGeneric("ppc_boxplot")
})

#' @export
#' @rdname ppc_plots
setMethod("ppc_boxplot", "bayesGAMfit",
          function(object, draws=NULL, ...) {
            pp <- posterior_predict(object, draws=draws)@pp
            
            y <- object@model@y
            ylst <- unlist(apply(y, 2, list), recursive = FALSE)
            
            pfun <- function(y, yrep, nms, ...) {
              bayesplot::ppc_boxplot(y, yrep, ...) + 
                ggplot2::ggtitle(nms) 
            }
            
            bp <- mapply(pfun, 
                         y = ylst, 
                         yrep = pp, 
                         nms = as.list(names(pp)), 
                         MoreArgs=list(... = ...), 
                         SIMPLIFY=FALSE)
            bayesplot::bayesplot_grid(plots=bp)
            
          })

#' @export
#' @rdname ppc_plots
setMethod("ppc_boxplot", "posteriorPredictObject",
          function(object, ...) {
            pp <- object@pp
            
            y <- object@model@y
            ylst <- unlist(apply(y, 2, list), recursive = FALSE)
            
            pfun <- function(y, yrep, nms, ...) {
              bayesplot::ppc_boxplot(y, yrep, ...) + 
                ggplot2::ggtitle(nms) 
            }
            
            bp <- mapply(pfun, 
                         y = ylst, 
                         yrep = pp, 
                         nms = as.list(names(pp)), 
                         MoreArgs=list(... = ...), 
                         SIMPLIFY=FALSE)
            bayesplot::bayesplot_grid(plots=bp)
            
          })

#' @export
#' @rdname ppc_plots
setGeneric("ppc_freqpoly", function(object, ...) {
  standardGeneric("ppc_freqpoly")
})

#' @export
#' @rdname ppc_plots
setMethod("ppc_freqpoly", "bayesGAMfit",
          function(object, draws=NULL, ...) {
            pp <- posterior_predict(object, draws=draws)@pp
            
            y <- object@model@y
            ylst <- unlist(apply(y, 2, list), recursive = FALSE)
            
            pfun <- function(y, yrep, nms, ...) {
              bayesplot::ppc_freqpoly(y, yrep, ...) + 
                ggplot2::ggtitle(nms) 
            }
            
            bp <- mapply(pfun, 
                         y = ylst, 
                         yrep = pp, 
                         nms = as.list(names(pp)), 
                         MoreArgs=list(... = ...), 
                         SIMPLIFY=FALSE)
            bayesplot::bayesplot_grid(plots=bp)
            
          })

#' @export
#' @rdname ppc_plots
setMethod("ppc_freqpoly", "posteriorPredictObject",
          function(object, ...) {
            pp <- object@pp
            
            y <- object@model@y
            ylst <- unlist(apply(y, 2, list), recursive = FALSE)
            
            pfun <- function(y, yrep, nms, ...) {
              bayesplot::ppc_freqpoly(y, yrep, ...) + 
                ggplot2::ggtitle(nms) 
            }
            
            bp <- mapply(pfun, 
                         y = ylst, 
                         yrep = pp, 
                         nms = as.list(names(pp)), 
                         MoreArgs=list(... = ...), 
                         SIMPLIFY=FALSE)
            bayesplot::bayesplot_grid(plots=bp)
            
          })



#' @export
#' @rdname ppc_plots
setGeneric("ppc_ecdf_overlay", function(object, ...) {
  standardGeneric("ppc_ecdf_overlay")
})

#' @export
#' @rdname ppc_plots
setMethod("ppc_ecdf_overlay", "bayesGAMfit",
          function(object, draws=NULL, ...) {
            pp <- posterior_predict(object, draws=draws)@pp 
            
            y <- object@model@y
            ylst <- unlist(apply(y, 2, list), recursive = FALSE)
            
            pfun <- function(y, yrep, nms, ...) {
              bayesplot::ppc_ecdf_overlay(y, yrep, ...) + 
                ggplot2::ggtitle(nms) 
            }
            
            bp <- mapply(pfun, 
                         y = ylst, 
                         yrep = pp, 
                         nms = as.list(names(pp)), 
                         MoreArgs=list(... = ...), 
                         SIMPLIFY=FALSE)
            bayesplot::bayesplot_grid(plots=bp)
            
          })

#' @export
#' @rdname ppc_plots
setMethod("ppc_ecdf_overlay", "posteriorPredictObject",
          function(object, ...) {
            pp <- object@pp 
            
            y <- object@model@y
            ylst <- unlist(apply(y, 2, list), recursive = FALSE)
            
            pfun <- function(y, yrep, nms, ...) {
              bayesplot::ppc_ecdf_overlay(y, yrep, ...) + 
                ggplot2::ggtitle(nms) 
            }
            
            bp <- mapply(pfun, 
                         y = ylst, 
                         yrep = pp, 
                         nms = as.list(names(pp)), 
                         MoreArgs=list(... = ...), 
                         SIMPLIFY=FALSE)
            bayesplot::bayesplot_grid(plots=bp)
            
          })

setMethod("initialize", "predictPlotObject",
          function(.Object, model, pp, ...) {
            .Object@model <- model
            .Object@pp <- pp
          
            return(.Object)
          })
            
setMethod("initialize", "posteriorPredictObject",
          function(.Object, model, pp, ...) {
            .Object@model <- model
            .Object@pp <- pp
            
            return(.Object)
          })

#########################################################
# predict method
#

#' Posterior predictive samples from models fit by \code{bayesGAM}, but with new data
#' 
#' Draw from the posterior predictive distribution applied to new data
#' 
#' @param object Object of type \code{bayesGAMfit} generated from \code{bayesGAM}.
#' @param newdata A data frame with new data applied to the \code{bayesGAMfit} object
#' @param draws An integer indicating the number of draws to return. The default and maximum number of draws is the size of the posterior sample.
#' @return a list of \emph{D} by \emph{N} matrices, where \emph{D} is the number of draws from the posterior predictive distribution and \emph{N} is the number of data points being predicted per draw.
#' @export
#' @references Goodrich B, Gabry J, Ali I & Brilleman S. (2020). rstanarm: Bayesian applied regression modeling via Stan. R package version 2.19.3 https://mc-stan.org/rstanarm.
#' @references Jonah Gabry, Ben Goodrich and Martin Lysy (2020). rstantools: Tools for Developing R Packages Interfacing with 'Stan'. https://mc-stan.org/rstantools/, https://discourse.mc-stan.org/.
#' @name predict
NULL

#' @rdname predict
#' @param ... Additional arguments for \code{postrior_predict}
#' @examples
#' set.seed(432)
#' f <- bayesGAM(weight ~ np(height), data = women,
#'               family = gaussian, iter=500, chains = 1)
#' newheights <- with(women, rnorm(10, mean = mean(height)), sd=sd(height))
#' women2 <- data.frame(height=newheights)
#' 
#' pred <- predict(f, women2, draws=100)
#' @export
setMethod("predict", signature(object="bayesGAMfit"),     
          function(object, newdata, draws=NULL, ...) {
  
  if (missing(newdata) | !is.data.frame(newdata)) {
    stop("newdata must be a data.frame")
  }
  
  ynms <- trimws(get_multy_names(object@model@names_y))
  
  nvals <- nrow(newdata)
  r <- object@model@r
  y <- object@model@y
  if (length(object@mcmcres) > 0) {
    mcmcres <- object@mcmcres  
  } else {
    mcmcres <- as.matrix(object@results)
  }
  
  newdata_y <- matrix(0, nrow=nvals, ncol=r)
  newdata_y <- as.data.frame(newdata_y)
  colnames(newdata_y) <- ynms
  
  # check if newdata already has these cols
  ycolstokeep <- ynms %in% colnames(newdata)
  newdata_y <- data.frame(newdata_y[, !ycolstokeep])
  if (ncol(newdata_y) > 0) {
    colnames(newdata_y)[!ycolstokeep] <- 
      ynms[!ycolstokeep]
    newdata <- cbind(newdata, newdata_y)
  }

  # np values from current model
  npargs <- object@model@npargs
  knots <- object@model@knots
  basis <- object@model@basis
  npdegree <- object@model@npdegree
  npterms <- object@model@npterms
  
  # get list of new npargs
  npnew <- mapply(append_knots_to_formula, 
                  nparg=npargs, 
                  kvals=knots,
                  basis=basis,
                  npdegree=npdegree, SIMPLIFY=FALSE)
  
  # change formula
  cl <- getCall(object@model)
  oldform <- cl$formula
  oldformchar <- as.character(oldform)
  
  newformchar <- oldformchar
  
  numnp <- length(npterms)
  iter <- 1
  
  while (iter <= numnp) {
    newformchar[3] <- gsub(npterms[iter], 
                           npnew[iter], 
                           newformchar[3], fixed=TRUE)
    iter <- iter + 1
  }
  
  newform <- as.formula(paste(newformchar[2], 
                              newformchar[1], 
                              newformchar[3]))  
  
  # get new glmModel from newdata
  cl$formula <- newform
  mod <- object@model
  mod@call <- cl
  newmodel <- update(mod, 
                     method="predict", 
                     data = newdata)
  
  # predictions
  res <- posterior_predict(newmodel, draws, results=mcmcres)
  
  # export object
  obj <- new("predictPlotObject", 
             model = newmodel, 
             pp = res)
  
  return(obj)
})


#' Extract the log likelihood from models fit by \code{bayesGAM}
#' 
#' Convenience function for extracting the pointwise log-likelihood matrix
#' or array from a model fit by \code{bayesGAM}. Calls the \code{extract_log_lik} method
#' from the \code{loo} package
#' 
#' @param object Object of type \code{bayesGAMfit} generated from \code{bayesGAM}.
#' @param ... Additional parameters to pass to \code{loo::extract_log_lik}
#' @return A matrix with the extracted log likelihood values post-warmup
#' @export
#' @references Stan Development Team (2017). The Stan C++ Library, Version 2.16.0. https://mc-stan.org/
#' @references Stan Development Team (2017). RStan: the R interface to Stan, Version 2.16.1. https://mc-stan.org/
#' @references Vehtari A, Gabry J, Magnusson M, Yao Y, Gelman A (2019). “loo: Efficient leave-one-out cross-validation and WAIC for Bayesian models.” R package version 2.2.0, <URL: https://mc-stan.org/loo>.
#' @references Vehtari A, Gelman A, Gabry J (2017). “Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC.” _Statistics and Computing_, *27*, 1413-1432. doi:10.1007/s11222-016-9696-4 (URL: https://doi.org/10.1007/s11222-016-9696-4).
#' @name extract_log_lik_bgam
#' @examples 
#' f <- bayesGAM(weight ~ np(height), data = women,
#'               family = gaussian, iter=500, chains = 1)
#' ll <- extract_log_lik_bgam(f)
NULL

#' @export
#' @rdname extract_log_lik_bgam
setGeneric("extract_log_lik_bgam", function(object, ...) {
  standardGeneric("extract_log_lik_bgam")
})

#' @export
#' @rdname extract_log_lik_bgam
setMethod("extract_log_lik_bgam", "bayesGAMfit", 
          function(object, ...) {
            loo::extract_log_lik(object@results, ...)
          })


#' Calls the \code{loo} package to perform efficient approximate leave-one-out cross-validation on models fit with \code{bayesGAM}
#' 
#' Computes PSIS-LOO CV, efficient approximate leave-one-out (LOO) cross-validation for 
#' Bayesian models using Pareto smoothed importance sampling (PSIS). This calls the implementation 
#' from the \code{loo} package of the methods described in Vehtari, Gelman, and Gabry (2017a, 2017b).
#' 
#' @param object Object of type \code{bayesGAMfit} generated from \code{bayesGAM}.
#' @param ... Additional parameters to pass to pass to \code{loo::loo}
#' @return a named list of class \code{c("psis_loo", "loo")}
#' \describe{
#'  \item{`estimates`}{
#'   A matrix with two columns (`Estimate`, `SE`) and three rows
#'   (`elpd_loo`, `p_loo`, `looic`). This
#'   contains point estimates and standard errors of the expected log pointwise
#'   predictive density (`elpd_loo`), the effective number of parameters
#'   (`p_loo`) and the LOO information criterion `looic` (which is
#'   just `-2 * elpd_loo`, i.e., converted to deviance scale).
#'  }
#'
#'  \item{`pointwise`}{
#'   A matrix with five columns (and number of rows equal to the number of
#'   observations) containing the pointwise contributions of the measures
#'   (`elpd_loo`, `mcse_elpd_loo`, `p_loo`, `looic`, `influence_pareto_k`).
#'   in addition to the three measures in \code{estimates}, we also report
#'   pointwise values of the Monte Carlo standard error of \code{elpd_loo}
#'   (\code{mcse_elpd_loo}), and statistics describing the influence of
#'   each observation on the posterior distribution (\code{influence_pareto_k}).
#'   These are the estimates of the shape parameter \eqn{k} of the
#'   generalized Pareto fit to the importance ratios for each leave-one-out
#'   distribution. See the [pareto-k-diagnostic] page for details.
#'  }
#'
#'
#'  \item{`diagnostics`}{
#'  A named list containing two vectors:
#'    * `pareto_k`: Importance sampling reliability diagnostics. By default,
#'      these are equal to the \code{influence_pareto_k} in \code{pointwise}.
#'      Some algorithms can improve importance sampling reliability and
#'      modify these diagnostics. See the [pareto-k-diagnostic] page for details.
#'    * `n_eff`: PSIS effective sample size estimates.
#'  }
#'
#'  \item{`psis_object`}{
#'  This component will be `NULL` unless the `save_psis` argument is set to
#'  `TRUE` when calling `loo()`. In that case `psis_object` will be the object
#'  of class `"psis"` that is created when the `loo()` function calls [psis()]
#'  internally to do the PSIS procedure.
#'  }
#' }
#'
#' @export
#' @references Vehtari, A., Gelman, A., and Gabry, J. (2017a). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing. 27(5), 1413–1432. doi:10.1007/s11222-016-9696-4 (journal version, preprint arXiv:1507.04544).
#' @references Vehtari, A., Gelman, A., and Gabry, J. (2017b). Pareto smoothed importance sampling. preprint arXiv:1507.02646
#' @references Vehtari A, Gabry J, Magnusson M, Yao Y, Gelman A (2019). “loo: Efficient leave-one-out cross-validation and WAIC for Bayesian models.” R package version 2.2.0, <URL: https://mc-stan.org/loo>.
#' @name loo_bgam
#' @examples
#' f <- bayesGAM(weight ~ np(height), data = women,
#'               family = gaussian, iter=500, chains = 1)
#' loo_bgam(f)
NULL

#' @export
#' @rdname loo_bgam
setGeneric("loo_bgam", function(object, ...) {
  standardGeneric("loo_bgam")
})

#' @export
#' @rdname loo_bgam
setMethod("loo_bgam", "bayesGAMfit", 
          function(object, ...) {
            ll <- extract_log_lik_bgam(object)
            loo::loo(ll, ...)
          })

#' @export
#' @rdname loo_bgam
setMethod("loo_bgam", "array", 
          function(object, ...) {
            loo::loo(object, ...)
          })

#' Calls the \code{loo} package to calculate the widely applicable information criterion (WAIC)
#' 
#' Computes WAIC by calling the appropriate function from the \code{loo} package
#' 
#' @param object Object of type \code{bayesGAMfit} generated from \code{bayesGAM}.
#' @param ... Additional parameters to pass to pass to \code{loo::waic}
#' @return a named list of class \code{c("waic", "loo")}
#' \describe{
#'  \item{`estimates`}{
#'  A matrix with two columns (`"Estimate"`, `"SE"`) and three
#'  rows (`"elpd_waic"`, `"p_waic"`, `"waic"`). This contains
#'  point estimates and standard errors of the expected log pointwise predictive
#'  density (`elpd_waic`), the effective number of parameters
#'  (`p_waic`) and the information criterion `waic` (which is just
#'  `-2 * elpd_waic`, i.e., converted to deviance scale).
#'  }
#'  \item{`pointwise`}{
#'  A matrix with three columns (and number of rows equal to the number of
#'  observations) containing the pointwise contributions of each of the above
#'  measures (`elpd_waic`, `p_waic`, `waic`).
#'  }
#' }
#'
#' @export
#' @references Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely application information criterion in singular learning theory. Journal of Machine Learning Research 11, 3571-3594.
#' @references Vehtari, A., Gelman, A., and Gabry, J. (2017a). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing. 27(5), 1413–1432. doi:10.1007/s11222-016-9696-4 (journal version, preprint arXiv:1507.04544).
#' @references Vehtari, A., Gelman, A., and Gabry, J. (2017b). Pareto smoothed importance sampling. preprint arXiv:1507.02646
#' @references Vehtari A, Gabry J, Magnusson M, Yao Y, Gelman A (2019). “loo: Efficient leave-one-out cross-validation and WAIC for Bayesian models.” R package version 2.2.0, <URL: https://mc-stan.org/loo>.
#' @name waic_bgam
#' @examples
#' f <- bayesGAM(weight ~ np(height), data = women,
#'               family = gaussian, iter=500, chains = 1)
#' waic_bgam(f)
NULL

#' @export
#' @rdname waic_bgam
setGeneric("waic_bgam", function(object, ...) {
  standardGeneric("waic_bgam")
})

#' @export
#' @rdname waic_bgam
setMethod("waic_bgam", "bayesGAMfit", 
          function(object, ...) {
            ll <- extract_log_lik_bgam(object)
            loo::waic(ll, ...)
          })

#' @export
#' @rdname waic_bgam
setMethod("waic_bgam", "array", 
          function(object, ...) {
            loo::waic(object, ...)
          })

#' Calls the \code{loo} package to compare models fit by \code{bayesGAMfit}
#' 
#' Compares fitted models based on ELPD, the expected log pointwise predictive 
#' density for a new dataset.  
#' 
#' 
#' @param object Object of type \code{bayesGAMfit} generated from \code{bayesGAM}.
#' @param ... Additional objects of type \code{bayesGAMfit}
#' @return a matrix with class \code{compare.loo} that has its own print method from the \code{loo} package
#' @export
#' @references Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely application information criterion in singular learning theory. Journal of Machine Learning Research 11, 3571-3594.
#' @references Vehtari, A., Gelman, A., and Gabry, J. (2017a). Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. Statistics and Computing. 27(5), 1413–1432. doi:10.1007/s11222-016-9696-4 (journal version, preprint arXiv:1507.04544).
#' @references Vehtari, A., Gelman, A., and Gabry, J. (2017b). Pareto smoothed importance sampling. preprint arXiv:1507.02646
#' @references Vehtari A, Gabry J, Magnusson M, Yao Y, Gelman A (2019). “loo: Efficient leave-one-out cross-validation and WAIC for Bayesian models.” R package version 2.2.0, <URL: https://mc-stan.org/loo>.
#' @references Gabry, J. , Simpson, D. , Vehtari, A. , Betancourt, M. and Gelman, A. (2019), Visualization in Bayesian workflow. J. R. Stat. Soc. A, 182: 389-402. doi:10.1111/rssa.12378
#' @examples
#' f1 <- bayesGAM(weight ~ height, data = women,
#'               family = gaussian, iter=500, chains = 1)
#' f2 <- bayesGAM(weight ~ np(height), data=women, 
#'               family = gaussian, iter=500, chains = 1)
#' loo_compare_bgam(f1, f2)
#' @name loo_compare_bgam
NULL


#' @export
#' @rdname loo_compare_bgam
setGeneric("loo_compare_bgam", function(object, ...) {
  standardGeneric("loo_compare_bgam")
})

#' @export
#' @rdname loo_compare_bgam
setMethod("loo_compare_bgam", "bayesGAMfit", 
          function(object, ...) {
              dots <- list(...)
              mvfit <- c(list(object), dots)
              
              all_loo <- lapply(mvfit, loo_bgam)
              loo::loo_compare(all_loo)
            }
          ) 

