###############################################
#
# Control function
#
###############################################



#' @export 
#' @title Control Function for Fitting a Multi-Type Regularized GLM Using the SMuRF Algorithm.
#' 
#' @description Control function to handle parameters for fitting a multi-type regularized generalized linear model (GLM) using the SMuRF algorithm.
#' The function sets defaults and performs input checks on the provided parameters. 
#' 
#' @param epsilon Numeric tolerance value for stopping criterion. A numeric strictly larger than 0, default is \code{1e-8}.
#' @param maxiter Maximum number of iterations of the SMuRF algorithm. A numeric larger than or equal to 1, default is \code{10 000}.
#' @param step Initial step size, a numeric strictly larger than 0 or \code{NULL}. When \code{NULL} (default), it is equal to \code{0.1} times the sample size.
#' @param tau Parameter for backtracking the step size. A numeric strictly between 0 and 1, default is 0.5.
#' @param reest A logical indicating if the obtained (reduced) model is re-estimated using \code{\link[speedglm]{speedglm}} or \code{\link[stats]{glm}}. Default is \code{TRUE}.
#' @param lambda.vector Values of lambda to consider when selecting the optimal value of lambda. A vector of strictly positive numerics (which is preferably a decreasing sequence as we make use of warm starts) or \code{NULL} (default).
#' When \code{NULL}, it is set to an exponential decreasing sequence of length \code{lambda.length} between \code{lambda.max} and \code{lambda.min}.
#' @param lambda.min Minimum value of lambda to consider when selecting the optimal value of lambda. A strictly positive numeric or \code{NULL} (default).
#' When \code{NULL}, it is equal to \code{0.0001} times \code{lambda.max}. This argument is ignored when \code{lambda.vector} is not \code{NULL}.
#' @param lambda.max Maximum value of lambda to consider when selecting the optimal value of lambda. A strictly positive numeric larger than \code{lambda.min} or \code{NULL} (default).
#' In the latter case, \code{lambda.max} will be determined based on the used penalty types such that it is one of the smallest values of lambda that results in an intercept-only model. This argument is ignored when \code{lambda.vector} is not \code{NULL}.
#' @param lambda.length Number of lambda values to consider when selecting the optimal value of lambda. A strictly positive integer, default is 50. This argument is ignored when \code{lambda.vector} is not \code{NULL}.
#' @param lambda.reest Logical indicating if the re-estimated coefficients are used when selecting lambda, default is \code{FALSE}.
#' This argument is only used if \code{reest} is \code{TRUE}.
#' @param k Number of folds when selecting lambda using cross-validation. A strictly positive integer, default is 5 (i.e. five-fold cross-validation). This number cannot be larger than the number of observations. Note that cross-validation with one fold (\code{k=1}) is the same as in-sample selection of \code{lambda}.
#' @param oos.prop Proportion of the data that is used as the validation sample when selecting \code{lambda} out-of-sample. A numeric strictly between 0 and 1, default is 0.2.
#' This argument is ignored when \code{validation.index} is not \code{NULL}. 
#' @param validation.index Vector containing the row indices of the data matrix corresponding to the observations that are used as the validation sample.
#' This argument is only used when \code{lambda} is selected out-of-sample. Default is \code{NULL} meaning that randomly 100*\code{oos.prop}\% of the data are used as validation sample.
#' @param ncores Number of cores used when performing cross-validation. A strictly positive integer or \code{NULL} (default). 
#' When \code{NULL}, \code{max(nc-1,1)} cores are used where \code{nc} is the number of cores as determined by \code{\link{detectCores}}.
#' @param po.ncores Number of cores used when computing the proximal operators. A strictly positive integer or \code{NULL} (default).
#' When \code{NULL} or \code{ncores > 1}, \code{po.ncores} is set to 1.
#' @param print A logical indicating if intermediate results need to be printed, default is \code{FALSE}.
#'
#' @details More details on the selection of lambda can be found in the package vignette.
#' @return A list with elements named as the arguments.
#'
#' @seealso Fitting procedures: \code{\link{glmsmurf}} and \code{\link{glmsmurf.fit}} (given design matrix). \code{\link[stats]{glm.control}}
#' 
#' @examples 
#' ## See example(plot_lambda) for examples 
glmsmurf.control <- function(epsilon = 1e-8, maxiter = 1e4, step = NULL, tau = 0.5, reest = TRUE,
                             lambda.vector = NULL, lambda.min = NULL, lambda.max = NULL, lambda.length = 50L, lambda.reest = FALSE, 
                             k = 5L, oos.prop = 0.2, validation.index = NULL, ncores = NULL, po.ncores = NULL, print = FALSE) {
  
  
  # step
  if (!is.null(step)) {
    
    if (!is.numeric(step) | length(step) != 1) {
      stop("'step' must be a numeric of length 1.")
    }
    
    if (is.nan(step) | is.infinite(step)) {
      stop("'step' must be a finite numeric of length 1.")
    }
    
    if (step <= 0) {
      stop("'step' must be strictly larger than 0.")
    }
  }

  
  # tau
  if (!is.numeric(tau) | length(tau) != 1) {
    stop("'tau' must be a numeric of length 1.")
  } 
  
  if (is.nan(tau) | is.infinite(tau)) {
    stop("'tau' must be a finite numeric of length 1.")
  }
  
  if (tau <= 0 | tau >= 1) {
    stop("'tau' must be strictly between 0 and 1.")
  }

  
  # epsilon
  if (!is.numeric(epsilon) | length(epsilon) != 1) {
    stop("'epsilon' must be a numeric of length 1.")
  }
  
  if (is.nan(epsilon) | is.infinite(epsilon)) {
    stop("'epsilon' must be a finite numeric of length 1.")
  }
  
  if (epsilon <= 0) {
    stop("'epsilon' must be strictly larger than 0.")
  }
  
  
  # maxiter
  if (!is.numeric(maxiter) | length(maxiter) != 1) {
    stop("'maxiter' must be a numeric of length 1.")
  }
  
  if (is.nan(maxiter) | is.infinite(maxiter)) {
    stop("'maxiter' must be a finite numeric of length 1.")
  }
  
  if (maxiter < 1) {
    stop("'maxiter' must be at least 1.")
  }
  
  
  # reest    
  reest_stop <- "'reest' must be a logical of length 1."
  
  if (is.null(reest) | length(reest) != 1) {
    stop(reest_stop)
    
  } else if (is.na(reest) | (!is.logical(reest) & !(reest %in% c(0, 1)))) {
    stop(reest_stop)
  }
  
  
  # k
  if (!is.numeric(k) | length(k) != 1) {
    stop("'k' must be a strictly positive integer of length 1.")
  }
  
  if (is.nan(k) | is.infinite(k)) {
    stop("'k' must be a finite strictly positive integer of length 1.")
  }
  
  if (k <= 0 | k %% 1 != 0) {
    stop("'k' must be a strictly positive integer.")
  }
  
  
  # oos.prop
  if (!is.numeric(oos.prop) | length(oos.prop) != 1) {
    stop("'oos.prop' must be a numeric of length 1.")
  }
  
  if (is.nan(oos.prop) | is.infinite(oos.prop)) {
    stop("'oos.prop' must be a finite numeric of length 1.")
  }
  
  if (oos.prop <= 0 | oos.prop >= 1) {
    stop("'oos.prop' must be a numeric strictly between 0 and 1.")
  }
  
  
  # validation.index
  if (!is.null(validation.index)) {
   
    if (!is.numeric(validation.index)) {
      stop("'validation.index' should be a numeric vector or NULL.")
    }

    if (!all(.is.wholenumber(validation.index)) | any(validation.index <= 0) | 
        any(is.nan(validation.index)) | any(is.infinite(validation.index))) {
      stop("All elements of 'validation.index' should be strictly positive integers.")
    }
  }
  
  
  # lambda.vector
  if (!is.null(lambda.vector)) {
    
    if (!is.numeric(lambda.vector) | length(lambda.vector) == 0L) {
      stop("'lambda.vector' must be a vector of strictly positive numbers or 'NULL'.")
    }
    
    if (any(is.nan(lambda.vector) | is.infinite(lambda.vector))) {
      stop("All elements of 'lambda.vector' should be strictly positive numbers.")
    }
    
    if (min(lambda.vector) <= 0) {
      stop("All elements of 'lambda.vector' should be strictly positive numbers.")
    }
  }
  
  
  # lambda.min
  if (!is.null(lambda.min)) {
    
    if (!is.numeric(lambda.min) | length(lambda.min) != 1) {
      stop("'lambda.min' must be a strictly positive number of length 1 or 'NULL'.")
    }
    
    if (is.nan(lambda.min) | is.infinite(lambda.min)) {
      stop("'lambda.min' must be a finite numeric of length 1 or 'NULL'.")
    }
    
    if (lambda.min <= 0) {
      stop("'lambda.min' must be a strictly positive number of length 1 or 'NULL'.")
    }
  }

  
  # lambda.max
  if (!is.null(lambda.max)) {
    
    if (!is.numeric(lambda.max) | length(lambda.max) != 1) {
      stop("'lambda.max' must be a strictly positive number of length 1 or 'NULL'.")
    }
    
    if (is.nan(lambda.max) | is.infinite(lambda.max)) {
      stop("'lambda.max' must be a finite numeric of length 1 or 'NULL'.")
    }
    
    if (lambda.max <= 0) {
      stop("'lambda.max' must be a strictly positive number of length 1 or 'NULL'.")
    }
    
    if (!is.null(lambda.min)) {
      if (lambda.max <= lambda.min) {
        stop("'lambda.max' must be strictly larger than 'lambda.min'.")
      }
    }
  }
  
  
  # lambda.length
  if (!is.numeric(lambda.length) | length(lambda.length) != 1) {
    stop("'lambda.length' must be a strictly positive integer of length 1.")
  }
  
  if (is.nan(lambda.length) | is.infinite(lambda.length)) {
    stop("'lambda.length' must be a strictly positive integer of length 1.")
  }
  
  if (lambda.length <= 0 | lambda.length %% 1 != 0) {
    stop("'lambda.length' must be a strictly positive integer of length 1.")
  }
  
  
  # lambda.reest
  if (is.null(lambda.reest) | length(lambda.reest) != 1) {
    stop("'lambda.reest' must be a logical of length 1.")
    
  } else if (is.na(lambda.reest) | (!is.logical(lambda.reest) & !(lambda.reest %in% c(0, 1)))) {
    stop("'lambda.reest' must be a logical of length 1.")
  }

  
  # Always set lambda.reest to FALSE if reest is FALSE
  if (!reest) {
    lambda.reest <- FALSE
  }
  
  
  # ncores
  if (is.null(ncores)) {
    
    # Use total number of cores minus one (or one if total number of cores is 1)
    ncores <- max(detectCores() - 1, 1)
    
  } else {
    
    if (!is.numeric(ncores) | length(ncores) != 1) {
      stop("'ncores' must be a positive integer of length 1 or 'NULL'.")
    }
    
    if (is.nan(ncores) | is.infinite(ncores)) {
      stop("'ncores' must be a positive integer of length 1 or 'NULL'.")
    }
    
    if (ncores <= 0 | ncores %% 1 != 0) {
      stop("'ncores' must be a positive integer or 'NULL'.")
    }
    
  }
  
  
  # po.ncores
  if (is.null(po.ncores)) {
    
    # Use one core to compute proximal operators
    po.ncores <- 1L
    
  } else {
    
    if (!is.numeric(po.ncores) | length(po.ncores) != 1) {
      stop("'po.ncores' must be a positive integer of length 1 or 'NULL'.")
    } 
    
    if (is.nan(po.ncores) | is.infinite(po.ncores)) {
      stop("'po.ncores' must be a positive integer of length 1 or 'NULL'.")
    }
    
    if (po.ncores <= 0 | po.ncores %% 1 != 0) {
      stop("'po.ncores' must be a positive integer or 'NULL'.")
    }

    if (ncores > 1 & po.ncores > 1) {
      po.ncores <- 1L
      warning("'po.ncores' is set to 1 since 'ncores' is already larger than 1.")
    }
    
    if (po.ncores > 1) {
      warning("Setting 'po.ncores' larger than one is not advised unless at least one of the ",
              "proximal operators takes long to compute.")
    }
    
  }
  
  
  # print
  print_stop <- "'print' must be a logical of length 1."
  
  if (is.null(print) | length(print) != 1) {
    stop(print_stop)
    
  } else if (is.na(print) | (!is.logical(print) & !(print %in% c(0, 1)))) {
    stop(print_stop)
  }
  
  
  return(list(step = step, tau = tau, epsilon = epsilon, maxiter = maxiter, k = k, oos.prop = oos.prop, 
              validation.index = validation.index, reest = reest,
              lambda.vector = lambda.vector, lambda.min = lambda.min, lambda.max = lambda.max, lambda.length = lambda.length, 
              lambda.reest = lambda.reest, ncores = ncores, po.ncores = po.ncores, print = print))
}

# Check if integer
#
# x: A number or a vector of numbers
# tol: Numerical tolerance value
.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {
  abs(x - round(x)) < tol
}
