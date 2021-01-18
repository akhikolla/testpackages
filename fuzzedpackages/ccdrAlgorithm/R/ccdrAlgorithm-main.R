#
#  ccdrAlgorithm-main.R
#  ccdrAlgorithm
#
#  Created by Bryon Aragam (local) on 1/22/16.
#  Copyright (c) 2014-2017 Bryon Aragam. All rights reserved.
#

#
# PACKAGE CCDRALGORITHM: Main CCDr methods
#
#   CONTENTS:
#     ccdr.run
#     ccdr_call
#     ccdr_gridR
#     ccdr_singleR
#

###--- These two lines are necessary to import the auto-generated Rcpp methods in RcppExports.R---###
#' @useDynLib ccdrAlgorithm, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Main CCDr Algorithm
#'
#' Estimate a Bayesian network (directed acyclic graph) from observational data using the
#' CCDr algorithm as described in \href{http://jmlr.org/papers/v16/aragam15a.html}{Aragam and Zhou (2015)}.
#'
#' Instead of producing a single estimate, this algorithm computes a solution path of estimates based
#' on the values supplied to \code{lambdas} or \code{lambdas.length}. The CCDr algorithm approximates
#' the solution to a nonconvex optimization problem using coordinate descent. Instead of AIC or BIC,
#' CCDr uses continuous regularization based on concave penalties such as the minimax concave penalty
#' (MCP).
#'
#' This implementation includes two options for the penalty: (1) MCP, and (2) L1 (or Lasso). This option
#' is controlled by the \code{gamma} argument.
#'
#' @param data Data as \code{\link[sparsebnUtils]{sparsebnData}} object. Must be numeric and contain no missing values.
#' @param lambdas Numeric vector containing a grid of lambda values (i.e. regularization
#'                parameters) to use in the solution path. If missing, a default grid of values will be
#'                used based on a decreasing log-scale  (see also \link{generate.lambdas}).
#' @param lambdas.length Integer number of values to include in the solution path. If \code{lambdas}
#'                       has also been specified, this value will be ignored. Note also that the final
#'                       solution path may contain fewer estimates (see
#'                       \code{alpha}).
#' @param whitelist A two-column matrix of edges that are guaranteed to be in each
#'                  estimate (a "white list"). Each row in this matrix corresponds
#'                  to an edge that is to be whitelisted. These edges can be
#'                  specified by node name (as a \code{character} matrix), or by
#'                  index (as a \code{numeric} matrix).
#' @param blacklist A two-column matrix of edges that are guaranteed to be absent
#'                  from each estimate (a "black list"). See argument
#'                  "\code{whitelist}" above for more details.
#' @param gamma Value of concavity parameter. If \code{gamma > 0}, then the MCP will be used
#'              with \code{gamma} as the concavity parameter. If \code{gamma < 0}, then the L1 penalty
#'              will be used and this value is otherwise ignored.
#' @param error.tol Error tolerance for the algorithm, used to test for convergence.
#' @param max.iters Maximum number of iterations for each internal sweep.
#' @param alpha Threshold parameter used to terminate the algorithm whenever the number of edges in the
#'              current DAG estimate is \code{> alpha * ncol(data)}.
#' @param betas Initial guess for the algorithm. Represents the weighted adjacency matrix
#'              of a DAG where the algorithm will begin searching for an optimal structure.
#' @param sigmas Numeric vector of known values of conditional variances for each node in the network. If this is
#'               set by the user, these parameters will not be computed and the input will
#'               be used as the "true" values of the variances in the algorithm. Note that setting
#'               this to be all ones (i.e. \code{sigmas[j] = 1} for all \code{j}) is
#'               equivalent to using the least-squares loss.
#' @param verbose \code{TRUE / FALSE} whether or not to print out progress and summary reports.
#'
#' @return A \code{\link[sparsebnUtils]{sparsebnPath}} object.
#'
#' @examples
#'
#' ### Generate some random data
#' dat <- matrix(rnorm(1000), nrow = 20)
#' dat <- sparsebnUtils::sparsebnData(dat, type = "continuous")
#'
#' # Run with default settings
#' ccdr.run(data = dat, lambdas.length = 20)
#'
#' ### Optional: Adjust settings
#' pp <- ncol(dat$data)
#'
#' # Initialize algorithm with a random initial value
#' init.betas <- matrix(0, nrow = pp, ncol = pp)
#' init.betas[1,2] <- init.betas[1,3] <- init.betas[4,2] <- 1
#'
#' # Run with adjusted settings
#' ccdr.run(data = dat, betas = init.betas, lambdas.length = 20, alpha = 10, verbose = TRUE)
#'
#' @export
ccdr.run <- function(data,
                     lambdas = NULL,
                     lambdas.length = NULL,
                     whitelist = NULL,
                     blacklist = NULL,
                     gamma = 2.0,
                     error.tol = 1e-4,
                     max.iters = NULL,
                     alpha = 10,
                     betas,
                     sigmas = NULL,
                     verbose = FALSE
){
    ### Check data format
    if(!sparsebnUtils::is.sparsebnData(data)) stop(sparsebnUtils::input_not_sparsebnData(data))

    ### Extract the data and ivn
    ### CCDr now works on both observational data and interventional data, and a mixture of both
    data_matrix <- data$data
    ivn_list <- data$ivn

    ### If ivn_list contains character names, convert to indices
    if("character" %in% sparsebnUtils::list_classes(ivn_list)){
        ivn_list <- lapply(ivn_list, function(x){
            idx <- match(x, names(data_matrix))
            if(length(idx) == 0) NULL # return NULL if no match (=> observational)
            else idx
        })
    }

    ### Call the CCDr algorithm
    ccdr_call(data = data_matrix,
              ivn = ivn_list,
              betas = betas,
              sigmas = sigmas,
              lambdas = lambdas,
              lambdas.length = lambdas.length,
              whitelist = whitelist,
              blacklist = blacklist,
              gamma = gamma,
              error.tol = error.tol,
              rlam = NULL,
              max.iters = max.iters,
              alpha = alpha,
              verbose = verbose)
} # END CCDR.RUN

### Maximum number of nodes allowed
MAX_CCS_ARRAY_SIZE <- function() 10000

# ccdr_call
#
#   Handles most of the bookkeeping for CCDr. Sets default values and prepares arguments for
#    passing to ccdr_gridR and ccdr_singleR. Some type-checking as well, although most of
#    this is handled internally by ccdr_gridR and ccdr_singleR.
#
ccdr_call <- function(data,
                      ivn = NULL,
                      betas,
                      sigmas,
                      lambdas,
                      lambdas.length,
                      whitelist = NULL,
                      blacklist = NULL,
                      gamma,
                      error.tol,
                      rlam,
                      max.iters,
                      alpha,
                      verbose = FALSE
){
    node_names <- names(data)
#     ### Allow users to input a data.frame, but kindly warn them about doing this
#     if(is.data.frame(data)){
#         warning(sparsebnUtils::alg_input_data_frame())
#         data <- sparsebnUtils::sparsebnData(data)
#     }
#
#     ### Check data format
#     if(!sparsebnUtils::is.sparsebnData(data)) stop(sparsebnUtils::input_not_sparsebnData(data))
#
#     ### Extract the data (CCDr only works on observational data, so ignore the intervention part)
#     data_matrix <- data$data

    ### Check data format
    if(!sparsebnUtils::check_if_data_matrix(data)) stop("'data' argument must be a data.frame or matrix!")

    # Could use check_if_complete_data here, but we avoid this in order to save a (small) amount of computation
    #  and give a more informative error message
    num_missing_values <- sparsebnUtils::count_nas(data)
    if(num_missing_values > 0) stop(sprintf("%d missing values detected!", num_missing_values))

    ### Get the dimensions of the data matrix
    nn <- as.integer(nrow(data))
    pp <- as.integer(ncol(data))

    if(pp > MAX_CCS_ARRAY_SIZE()){
        stop(max_nodes_warning(pp))
    }

    if(is.null(ivn)) ivn <- vector("list", nn) # to pass testthat for observational data cases
    ### Check ivn
    if(!check_if_ivn_list(ivn)) stop("ivn must be a list of NULLs or vectors!")
    if(!check_ivn_size(ivn, data)) stop(sprintf("Length of ivn is %d, expected to match the number of rows in data: %d.", length(ivn), nn))
    check_ivn_label(ivn, data)
    ### if(!check_ivn_label(ivn, data)) stop("Intervention labels are incorrect.")

    ### use a vector nj to count how many times each node is under intervention
    ### refer to nj as "intervention times vector"
    nj <- rep(0, pp)
    for(j in 1:pp) { ## include 0 here or not?
        nj[j] <- sum(!sapply(lapply(ivn, is.element, j), any)) ## optimize for sorted column?
    }

    ### Set default for sigmas (negative values => ignore initial value and update as usual)
    if(is.null(sigmas)){
        sigmas <- rep(-1., pp)
    }

    ### Use default values for lambda if not specified
    if(is.null(lambdas)){
        if(is.null(lambdas.length)){
            stop("Both lambdas and lambdas.length unspecified: Must specify a value for at least one of these arguments!")
        } else{
            ### Check lambdas.length if specified
            if(!is.numeric(lambdas.length)) stop("lambdas.length must be numeric!")
            if(lambdas.length <= 0) stop("lambdas.length must be positive!")
        }

        if(missing(rlam)){
            ### Even though ccdr_call should never be called on its own, this behaviour is left for testing backwards-compatibility
            stop("rlam must be specified if lambdas is not explicitly specified.")
        } else if(is.null(rlam)){
            ### rlam = NULL is used as a sentinel value to indicate a default value should be used
            rlam <- 1e-2
        } else{
            ### Check rlam if specified
            if(!is.numeric(rlam)) stop("rlam must be numeric!")
            if(rlam < 0) stop("rlam must be >= 0!")
        }

        # If no grid of lambdas is passed, then use the standard log-scale that starts at
        #  max.lam = sqrt(nn) and descends to min.lam = rlam * max.lam
        lambdas <- sparsebnUtils::generate.lambdas(lambda.max = sqrt(nn),
                                                   lambdas.ratio = rlam,
                                                   lambdas.length = as.integer(lambdas.length),
                                                   scale = "log")
    }

    ### Check lambdas
    if(!is.numeric(lambdas)) stop("lambdas must be a numeric vector!")
    if(any(lambdas < 0)) stop("lambdas must contain only nonnegative values!")

#     if(length(lambdas) != lambdas.length){
#         warning("Length of lambdas vector does not match lambdas.length. The specified lambdas vector will be used and lambdas.length will be overwritten.")
#     }

    ### By default, set the initial guess for betas to be all zeroes

    if(missing(betas)){
        betas <- matrix(0, nrow = pp, ncol = pp)
        # betas <- SparseBlockMatrixR(betas) # 2015-03-26: Deprecated and replaced with .init_sbm below
        betas <- .init_sbm(betas, rep(0, pp))

        # If the initial matrix is the zero matrix, indexing does not matter so we don't need to use reIndexC here
        #   Still need to set start = 0, though.
        betas$start <- 0
    } # Type-checking for betas happens in ccdr_singleR

    # This parameter can be set by the user, but in order to prevent the algorithm from taking too long to run
    #  it is a good idea to keep the threshold used by default which is O(sqrt(pp))
    if(is.null(max.iters)){
        max.iters <- sparsebnUtils::default_max_iters(pp)
    }

    ### White/black lists
    # Be careful about handling various NULL cases
    if(!is.null(whitelist)) whitelist <- bwlist_check(whitelist, node_names)
    if(!is.null(blacklist)) blacklist <- bwlist_check(blacklist, node_names)

    if(!is.null(whitelist) && !is.null(blacklist)){
        if(length(intersect(whitelist, blacklist)) > 0){
            badinput <- vapply(intersect(whitelist, blacklist), function(x) sprintf("\t[%s]\n", paste(x, collapse = ",")), FUN.VALUE = "vector")
            badinput <- paste(badinput, collapse = "")
            msg <- sprintf("Duplicate entries found in blacklist and whitelist: \n%s", badinput)
            stop(msg)
        }
    }

    weights <- bwlist_to_weights(blacklist, whitelist, nnode = pp)

    ### Pre-process correlation data
    t1.cor <- proc.time()[3]
    #     cors <- cor(data)
    #     cors <- cors[upper.tri(cors, diag = TRUE)]
    corlist <- sparsebnUtils::cor_vector_ivn(data, ivn)
    cors <- corlist$cors
    indexj <- corlist$indexj
    t2.cor <- proc.time()[3]

    fit <- ccdr_gridR(cors,
                      as.integer(pp),
                      as.integer(nn),
                      as.integer(nj),
                      as.integer(indexj),
                      betas,
                      as.numeric(sigmas),
                      as.numeric(lambdas),
                      as.integer(weights),
                      as.numeric(gamma),
                      as.numeric(error.tol),
                      as.integer(max.iters),
                      as.numeric(alpha),
                      verbose)

    #
    # Output DAGs as edge lists (i.e. edgeList objects).
    #  This is NOT the same as sbm$rows since some of these rows may correspond to edges with zero coefficients.
    #  See docs for SparseBlockMatrixR class for details.
    #
    for(k in seq_along(fit)){
        ### Coerce sbm output to edgeList
        names(fit[[k]])[1] <- "edges" # rename 'sbm' slot to 'edges': After the next line, this slot will no longer be an SBM object
        fit[[k]]$edges <- sparsebnUtils::as.edgeList(fit[[k]]$edges) # Before coercion, li$edges is actually an SBM object
        names(fit[[k]]$edges) <- names(data)

        ### Add node names to output
        fit[[k]] <- append(fit[[k]], list(node_names), after = 1) # insert node names into second slot
        names(fit[[k]])[2] <- "nodes"
    }

    fit <- lapply(fit, sparsebnUtils::sparsebnFit)    # convert everything to sparsebnFit objects
    sparsebnUtils::sparsebnPath(fit)                  # wrap as sparsebnPath object
} # END CCDR_CALL

# ccdr_gridR
#
#   Main subroutine for running the CCDr algorithm on a grid of lambda values.
ccdr_gridR <- function(cors,
                       pp, nn,
                       nj = NULL,
                       indexj = NULL,
                       betas,
                       sigmas,
                       lambdas,
                       weights,
                       gamma,
                       eps,
                       maxIters,
                       alpha,
                       verbose
){

    ### Check alpha
    if(!is.numeric(alpha)) stop("alpha must be numeric!")
    if(alpha < 0) stop("alpha must be >= 0!")

    ### nlam is now set automatically
    nlam <- length(lambdas)

    ### Check indexj
    if(is.null(indexj)) indexj <- rep(0L, pp + 1)
    ### Check nj
    if(is.null(nj)) nj <- as.integer(rep(nn, pp))

    ccdr.out <- list()
    for(i in 1:nlam){

        if(verbose) message("Working on lambda = ", round(lambdas[i], 5), " [", i, "/", nlam, "]")

        t1.ccdr <- proc.time()[3]
        ccdr.out[[i]] <- ccdr_singleR(cors,
                                      pp, nn,
                                      nj,
                                      indexj,
                                      betas,
                                      sigmas,
                                      lambdas[i],
                                      weights,
                                      gamma = gamma,
                                      eps = eps,
                                      maxIters = maxIters,
                                      alpha = alpha,
                                      verbose = verbose
        )
        t2.ccdr <- proc.time()[3]

        betas <- ccdr.out[[i]]$sbm
        betas <- sparsebnUtils::reIndexC(betas) # use C-friendly indexing

        if(verbose){
            test.nedge <- sum(as.matrix(betas) != 0)
            message("  Estimated number of edges: ", ccdr.out[[i]]$nedge)
            # message("  Estimated total variance: ", sum(1 / (betas$sigmas)^2))
        }

        # 7-16-14: Added code below to check edge threshold via alpha parameter
        if(ccdr.out[[i]]$nedge > alpha * pp){
            if(verbose) message("Edge threshold met, terminating algorithm with ", ccdr.out[[i-1]]$nedge, " edges.")
            ccdr.out <- ccdr.out[1:(i-1)] # only return up to i - 1 since the last (ith) model did not finish
            break
        }
    }

    ccdr.out
} # END CCDR_GRIDR

# ccdr_singleR
#
#   Internal subroutine for handling calls to singleCCDr: This is the only place where C++ is directly
#    called. Type-checking is strongly enforced here.
ccdr_singleR <- function(cors,
                         pp, nn,
                         nj = NULL,
                         indexj = NULL,
                         betas,
                         sigmas,
                         lambda,
                         weights,
                         gamma,
                         eps,
                         maxIters,
                         alpha,     # 2-9-15: No longer necessary in ccdr_singleR, but needed since the C++ call asks for it
                         verbose = FALSE
){

    ### Check dimension parameters
    if(!is.integer(pp) || !is.integer(nn)) stop("Both pp and nn must be integers!")
    if(pp <= 0 || nn <= 0) stop("Both pp and nn must be positive!")

    ### These variables, if NULL, need to be initialized before checking anything
    if(is.null(indexj)) indexj <- rep(0L, pp + 1) # initialize indexj
    if(is.null(nj)) nj <- as.integer(rep(nn, pp)) # initialize nj

    ### Check indexj
    if(!is.vector(indexj)) stop("Index vector for cors is not a vector.")
    if(length(indexj) > pp + 1) stop(sprintf("Index vector for cors is too long, expected to be no greater than %d, the number of columns of data.", pp))
    if(!is.integer(indexj)) stop("Index vector for cors has non-integer component(s).")
    if(any(is.na(indexj) | is.null(indexj))) stop("Index vector cannot have missing or NULL values.")
    if(any(indexj < 0 | indexj > pp + 1)) stop(sprintf("Index vector for cors has out-of-range component(s), expected to be between 0 and %d.", pp))

    ### Check nj
    if(!is.vector(nj)) stop("Intervention times vector is not a vector.")
    if(length(nj) != pp) stop(sprintf("Length of intervention times vector is %d, expected to match the number of columns of data = %d", length(nj), pp))
    if(!is.integer(nj)) stop("Intervention times vector has non-integer component(s).")
    if(any(is.na(nj) | is.null(nj))) stop("Intervention times vector cannot have missing or NULL values.")
    if(any(nj < 0 | nj > nn)) stop(sprintf("Intervention times vector has out-of-range component(s), expected to be between 0 and %d.", nn))

    ### Check cors
    ### This check must come after the checks for indexj, nj since these values are used to check cors
    if(!is.numeric(cors)) stop("cors must be a numeric vector!")
    if(length(cors) != length(unique(indexj))*pp*(pp+1)/2) stop(paste0("cors has incorrect length: Expected length = ", length(unique(indexj))*pp*(pp+1)/2, " input length = ", length(cors)))

    ### add a weight a_j to penalty on beta_{ij}
    ### since now with intervention data, beta_{ij} only appears n_j times out of total nn samples
    aj <- nj / nn

    ### Check betas
    if(sparsebnUtils::check_if_matrix(betas)){ # if the input is a matrix, convert to SBM object
        betas <- .init_sbm(betas, rep(0, pp)) # if betas is non-numeric, SparseBlockMatrixR constructor should throw error
        betas <- reIndexC(betas) # use C-friendly indexing
    } else if(!is.SparseBlockMatrixR(betas)){ # otherwise check that it is an object of class SparseBlockMatrixR
        stop("Incompatible data passed for betas parameter: Should be either matrix or list in SparseBlockMatrixR format.")
    }

    ### Check sigmas
    if(!is.numeric(sigmas)) stop("sigmas must be numeric!")
    if(length(sigmas) != pp) stop(sprintf("sigmas must have length = %d!", pp))
    if(any(sigmas < 0)){
        # -1 is a sentinel value for updating sigmas via the CD updates
        if(any(sigmas != -1.)){
            stop("sigmas must be > 0!")
        }
    }

    ### Check lambda
    if(!is.numeric(lambda)) stop("lambda must be numeric!")
    if(lambda < 0) stop("lambda must be >= 0!")

    ### Check weights
    if(length(weights) != pp*pp) stop(sprintf("weights must have length p^2 = %d!", pp*pp))
    if(!is.numeric(weights)) stop("weights must be numeric!")
    if(weights < -1 || weights > 1) stop("weights out of bounds!")

    ### Check gamma
    if(!is.numeric(gamma)) stop("gamma must be numeric!")
    if(gamma < 0 && gamma != -1) stop("gamma must be >= 0 (MCP) or = -1 (Lasso)!")

    ### Check eps
    if(!is.numeric(eps)) stop("eps must be numeric!")
    if(eps <= 0){
        if(eps < 0) stop("eps must be >= 0!")
        if(eps == 0) warning("eps is set to zero: This may cause the algorithm to fail to converge, and maxIters will be used to terminate the algorithm.")
    }

    ### Check maxIters
    if(!is.integer(maxIters)) stop("maxIters must be an integer!")
    if(maxIters <= 0) stop("maxIters must be > 0!")

    ### alpha check is in ccdr_gridR

    # if(verbose) cat("Opening C++ connection...")
    t1.ccdr <- proc.time()[3]
    ccdr.out <- singleCCDr(cors,
                           betas,
                           sigmas,
                           nj,
                           indexj,
                           aj,
                           lambda,
                           weights,
                           c(gamma, eps, maxIters, alpha),
                           verbose = verbose)
    t2.ccdr <- proc.time()[3]
    # if(verbose) cat("C++ connection closed. Total time in C++: ", t2.ccdr-t1.ccdr, "\n")

    #
    # Convert output back to SBM format
    #
    ccdr.out <- list(sbm = SparseBlockMatrixR(list(rows = ccdr.out$rows, vals = ccdr.out$vals, blocks = ccdr.out$blocks, sigmas = ccdr.out$sigmas, start = 0)),
                     lambda = ccdr.out$lambda,
                     nedge = ccdr.out$length,
                     pp = pp,
                     nn = nn,
                     time = t2.ccdr - t1.ccdr)
    ccdr.out$sbm <- sparsebnUtils::reIndexR(ccdr.out$sbm)

    # sparsebnFit(ccdr.out)
    ccdr.out
} # END CCDR_SINGLER
