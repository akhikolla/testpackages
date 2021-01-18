
parse_priors <- function(priors, famnum, X, Z) {
  # default priors
  # families:  1 = gaussian, 2 = binomial, 3 = poisson
  if (famnum %in% c(1, 2, 3)) {
    if (is.null(priors$beta))
      priors$beta <- c(0, 1e6)
  } else {
    stop("Invalid family selected")
  }

  # random effects variance
  if (!is.null(Z)) {
    if (is.null(priors$lambda))
      priors$lambda <- c(35, 0, 1)
  }

  # residual variance
  if (famnum == 1) {
    if (is.null(priors$eps))
      priors$eps <- c(35, 0, 1)
  }

  # convert beta priors to matrix
  p <- ncol(X)
  if (p > 1 && !is.matrix(priors$beta)) {
    priors$beta <- matrix(rep(priors$beta, p), ncol=2, byrow=TRUE)
  }

  return(priors)
}

parse_formula <- function(form)
{
  formtochar <- as.character(form)
  rhs <- formtochar[[3]]    
  
  # extract all terms separated by +
  allterms <- trimws(strsplit(rhs, "\\+")[[1]])

  # check pattern match
  checknp <- sapply(allterms, grepl, pattern="np\\((.*)\\)")

  # extract nonparametric portions of formula
  npterms <- allterms[checknp]
  otherterms <- allterms[!checknp]
  otherterms <- paste(otherterms, collapse = " + ")

  # include intercept if needed
  if (otherterms == "") {
    otherterms <- "1"
  }

  # create new formula after extracting nonlinear
  sub_form <- paste(formtochar[[2]], formtochar[[1]],
                    otherterms, collapse=" ")
  sub_form <- as.formula(sub_form)
  retval <- list(npterms = npterms,
                 sub_form = sub_form)

  return(retval)
}

# function to get multiple response names
get_multy_names <- function(yvar) {
  checky <- sapply(yvar, grepl, pattern="cbind\\((.*)\\)")
  yvals <- yvar
  if (checky) {
    temp <- gsub("cbind", "", yvar)
    temp <- gsub("\\(", "", temp)
    temp <- gsub("\\)", "", temp)
    yvals <- unlist(strsplit(temp, ",", fixed=TRUE))
  }

  return(yvals)

}


# function to get multresponse values
get_multresponse_ynums <- function(ynm) {
  mult_pattern <- "\\[.*\\,"
  result <- regmatches(ynm, regexec(mult_pattern, ynm))
  result <- gsub("\\[", "", result)
  result <- as.integer(gsub("\\,", "", result))
  return(result)
}

get_multresponse_xnums <- function(xnm) {
  mult_pattern <- "\\,.*\\]"
  result <- regmatches(xnm, regexec(mult_pattern, xnm))
  result <- gsub("\\]", "", result)
  result <- as.integer(gsub("\\,", "", result))
  return(result)
}

# function to apply given vector length
pastelim <- function(x, nparam, ...) {
  paste(x[1:nparam], ...)
}

# get name of y
getResponse <- function(formula, data=NULL) {
  tt <- terms(formula, data=data)
  vars <- as.character(attr(tt, "variables"))[-1]
  response <- attr(tt, "response")
  vars[response]
}

# get np names.  returns list
getNPnames <- function(txt) {
  txt <- gsub("^np\\(", "", txt)
  txt <- gsub("\\)", "", txt)
  ss <- strsplit(txt,",")
  n <- sapply(ss,length)
  ss
}


# pad list of matrices
pad_z <- function(Zlst, padval=0) {
  if (length(Zlst) == 0) {
    retval <- list(max_col=0L,
         num_z=0L,
         Zarray = array(0))
  } else {
    max_col <- max(sapply(Zlst, ncol))
    num_row <- sapply(Zlst, nrow)[1]

    new_z <- lapply(Zlst, function(xx) {
      res <- xx
      if (ncol(xx) < max_col) {
        tmp <- matrix(padval, nrow=num_row, ncol=max_col - ncol(res) )
        res <- cbind(res, tmp)
      }
      res
    })

    # coerce to array
    num_z <- length(Zlst)
    Zarray <- simplify2array(new_z)
    if (!is.array(Zarray)) stop("Z not converted to array")
    Zarray <- aperm(Zarray, c(3, 1, 2))

    list(max_col=as.integer(max_col),
         num_z = length(Zlst),
         Zarray = Zarray)
  }

}


#' Lag function for autoregressive models
#' 
#' Creates lagged variables for use with \code{bayesGAM}, including the functionality
#' to create lags for each specified subject if desired. The input data must be pre-
#' sorted according by time, and within each subject id if specified. 
#' 
#' @export
#' @param x numeric vector
#' @param k integer vector of lagged variables to create
#' @param id optional identification number for each subject
#' @references Zeileis A (2019). dynlm: Dynamic Linear Regression. R package version 0.3-6
#' @return numeric vector or matrix of the lagged variable(s)
#' @examples
#' x <- rnorm(20)
#' id <- rep(1:4, each=5)
#' L(x, 1:2, id)
#' 
#' # autoregressive
#' ar.ols(lh, demean = FALSE, intercept=TRUE, order=1)
#' f <- bayesGAM(lh ~ L(lh), family=gaussian)
#' coef(f)
L <- function(x, k=1, id=NULL) {
  arg <- deparse(substitute(x))
  
  if (is.null(id)) {
    res <- Lbase(x=x, k=k) 
    if (length(k) > 1) {
      res <- as.matrix(res)
      colnames(res) <- paste0("L", arg, k)
    }
  } else {
    if (is.factor(id)) {
      id <- as.character(id)
    }
    
    fid <- split(x, factor(id))
    
    if (length(k) == 1) {
      res <- lapply(fid, function(zz) {
        rv <- Lbase(zz, k=k)
        rv[1:length(zz)]
      })
      res <- do.call(c, res) 
    } else {
      res <- lapply(fid, function(zz) {
        rv <- Lbase(zz, k=k)
        rv[1:length(zz), ]
      })
      
      res <- as.matrix(do.call(rbind, res))
      colnames(res) <- paste0("L", arg, k)
    }
    
  }
  return(res)
}

Lbase <- function(x, k = 1) {
  x <- as.ts(x)
  
  if (length(k) > 1) {
    res <- lapply(k, function(i) lag(x, k=-i))
    res <- do.call(ts.union, c(x=list(x), res)) 
    res <- as.matrix(res[1:length(x), -1])
  } else {
    res <- ts.union(x, lag(x, k = -k), dframe=TRUE)[, -1]
    res <- as.numeric(res)[1:length(x)]
  }
  return(res)
}

# helper function for posterior predictive sampling
pp_sim <- function(index, r, W, famnum, linknum, offset, mcmcres_betau, mcmcres_eps=NULL, multresponse=FALSE) {
  # identity
  if (linknum == 1) {
    linkinvFUN <- identity
    # log
  } else if (linknum == 2) {
    linkinvFUN <- exp
    # inverse
  } else if (linknum == 3) {
    linkinvFUN <- function(x) { 1/x }
  } else if (linknum == 4) {
    linkinvFUN <- boot::inv.logit
  } else if (linknum == 5) {
    linkinvFUN <- pnorm
  } else if (linknum == 6) { 
    linkinvFUN <- atan
  } else if (linknum == 7) {
    linkinvFUN <- function(x) { 1 - exp(-exp(x))}
  } else if (linknum == 8) {
    linkinvFUN <- function(x) { x^2 }
  } else {
    stop("invalid link number")
  }
 
  if (multresponse) {
    # lists of sample matrices, for each dependent variable
    betaucols <- strsplit(colnames(mcmcres_betau), 
                          split="\\[|\\,|\\]")
    nyvals <- sapply(betaucols, function(xx) {
      as.numeric(xx[2])
    })
    
    mcmcres_betau_lst <- split(as.data.frame(t(mcmcres_betau)), 
                               nyvals) 
    mcmcres_betau_lst <- lapply(mcmcres_betau_lst, function(xx) {
      as.matrix(t(xx))
    })
    
    # gaussian multivariate response has additional param
    if (famnum == 1) {
      epscols <- strsplit(colnames(mcmcres_eps), 
                          split="\\[|\\,|\\]")
      nepsvals <- sapply(epscols, function(xx) {
        as.numeric(xx[2])
      })
      mcmcres_eps_lst <- split(as.data.frame(t(mcmcres_eps)), 
                               nepsvals)
      mcmcres_eps_lst <- lapply(mcmcres_eps_lst, function(xx) {
        as.matrix(t(xx))
      })
      
      allres <- mapply(pp_samples, 
                        mcmcres_betau=mcmcres_betau_lst, 
                        mcmcres_eps=mcmcres_eps_lst, 
                       MoreArgs=list(index=index, 
                                     W=W, 
                                     famnum=famnum, 
                                     offset=offset,
                                     linkinvFUN=linkinvFUN), 
                       SIMPLIFY=FALSE)
      
    } else {
      allres <- lapply(mcmcres_betau_lst, FUN=pp_samples, 
                       index=index, 
                       W=W, 
                       famnum=famnum, 
                       offset=offset, 
                       mcmcres_eps=mcmcres_eps, 
                       linkinvFUN=linkinvFUN)
    }
    
    names(allres) <- paste0("yrep", 1:r)
    
    
  } else {
    res <- pp_samples(index, W, famnum, mcmcres_betau, 
                      mcmcres_eps, offset, linkinvFUN)
    allres <- list(yrep = res)
  }
  
  return(allres)
}

pp_samples <- function(index, W, famnum, mcmcres_betau, mcmcres_eps, 
                       offset, linkinvFUN) {
  
  # gaussian
  if (famnum == 1) {
    res <- lapply(index, function(i) {
      betau <- mcmcres_betau[i, ]
      eps <- mcmcres_eps[i]
      meanvals <- linkinvFUN(W %*% betau + offset)
      sapply(meanvals, FUN=rnorm, n=1, sd=eps)
    })    
    # binomial
  } else if (famnum == 2) {
    res <- lapply(index, function(i) {
      betau <- mcmcres_betau[i, ]
      meanvals <- linkinvFUN(W %*% betau + offset)
      sapply(meanvals, FUN=rbinom, n=1, size=1)
    })
    # poisson
  } else if (famnum == 3) {
    res <- lapply(index, function(i) {
      betau <- mcmcres_betau[i, ]
      meanvals <- linkinvFUN(W %*% betau + offset) 
      sapply(meanvals, FUN=rpois, n=1)
    })
  }
  res <- as.matrix(do.call(rbind, res))
  return(res)
}


# ensure same knots for prediction with new data
append_knots_to_formula <- function(nparg, kvals, basis, npdegree) {
  knstring <- paste(kvals, collapse=',')
  
  # check for bivariate
  if (length(nparg) > 1) {
    nparg <- paste(nparg, collapse=", ")
  }
  
  newargs <- paste0(nparg, ", ", 
                    "knots=c(", knstring, "), ", 
                    "basis=", "\'", basis, "\'")
  
  newcall <- paste0("np(", newargs, ")")
  
  return(newcall)
}

# set variable names for bayesplot. object is type bayesGAMfit, returns stan object
set_varnms <- function(object) {
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
      paste("beta", object@model@names_beta, sep="_")
  }
  
  return(stanobj)
}
