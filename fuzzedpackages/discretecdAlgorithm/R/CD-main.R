# You can learn more about package authoring with RStudio at: http://r-pkgs.had.co.nz/
# Some useful keyboard shortcuts for package authoring: Build and Reload Package: 'Cmd +
# Shift + B' Check Package: 'Cmd + Shift + E' Test Package: 'Cmd + Shift + T'

# ========================================================
# Note that this is only a trial version.
# There will be no such choice that we can input the initial beta.
# No documentation yet
# ========================================================

#' @useDynLib discretecdAlgorithm, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#========================================================
# The main function CD.run
# ========================================================
## will be exported

#' cd.run
#'
#' Structure learning of discrete Bayesian network
#'
#' Estimate structure of a discrete Bayesian network from observational/interventional data using the CD algorithm described in \href{http://arxiv.org/abs/1403.2310}{Gu et al. (2016)}.
#'
#' Instead of producing a single estimate, this algorithm computes a solution path of estimates based
#' on the values supplied to \code{lambdas} or \code{lambdas.length}.
#' This package do not provide a model selection method in this version, users can choose their own model selection criterion.
#' In later version of this package we will provide an empirical model selection method.
#'
#' This package can handle interventional data by input a list of intervention. See example for more detail.
#'
#' @param indata A sparsebnData object.
#' @param weights Weight matrix. Weight can be the \code{l_2} norm of a consistent estimate of \code{beta_{j.i}}.
#'                See paper \href{http://arxiv.org/abs/1403.2310}{Gu et al. (2016)} chapter 3.3 for more details.
#'                A weight matrix that is set improperly may cause convergence issues and lead to a suboptimal solution.
#' @param lambdas Numeric vector containing a grid of lambda values (i.e. regularization parameters)
#'                to use in the solution path. If missing, a default grid of values will be used based on a decreasing log-scale.
#'                To generate a sequence of lambdas see \code{\link[sparsebnUtils]{generate.lambdas}}.
#'                For discrete network, the paper provided a way to calculate a maximum lambda that penalizes all parameters to zero,
#'                \href{http://arxiv.org/abs/1403.2310}{Gu et al. (2016)} chapter 3.4.
#'                See function \code{\link{max_lambda}} for details.
#' @param lambdas.length Integer number of values to include in the solution path.
#' @param whitelist A two-column matrix of edges that are guaranteed to be in each
#'                  estimate (a "white list"). Each row in this matrix corresponds
#'                  to an edge that is to be whitelisted. These edges can be
#'                  specified by node name (as a \code{character} matrix), or by
#'                  index (as a \code{numeric} matrix).
#' @param blacklist A two-column matrix of edges that are guaranteed to be absent
#'                  from each estimate (a "black list"). See argument
#'                  "\code{whitelist}" above for more details.
#' @param error.tol Error tolerance for the algorithm, used to test for convergence.
#' @param convLb Small positive number used in Hessian approximation.
#' @param weight.scale A positive number to scale weight matrix.
#' @param upperbound A large positive value used to truncate the adaptive weights.
#'                   A -1 value indicates that there is no truncation.
#' @param alpha Threshold parameter used to terminate the algorithm whenever the number of edges in the
#'              current DAG estimate is \code{> alpha * ncol(data)}.
#' @param permute A bool parameter, default value is FALSE. If TRUE, will randomize order of going through blocks.
#' @param adaptive A bool parameter, default value is FALSE. If FALSE, a regular lasso algorithm will be run.
#'                 If TRUE, an adaptive lasso algorithm will be run.
#' @return A \code{\link[sparsebnUtils]{sparsebnPath}} object.
#'         The CD Algorithm will be stopped if the number of edges exceeds 3 times of number of variables.
#' @examples
#'
#' \dontrun{
#'
#' ### Generate some random data
#' dat <- matrix(rbinom(200, size = 3, prob = 0.4), nrow = 20)
#' # for observational data
#' dat_obs <- sparsebnUtils::sparsebnData(dat, type = "discrete")
#' # for interventional data
#' data_size <- nrow(dat)
#' ivn <- lapply(1:data_size, function(x){return(as.integer(x/10))})
#' # if there is no intervention for an observation, use 0.
#' # cd algorithm can handle multiple interventions for a single observation.
#' dat_int <- sparsebnUtils::sparsebnData(dat, ivn = ivn, type = "discrete")
#'
#' # Run with default settings for observational data
#' cd.run(indata = dat_obs)
#' # Run with default settings for interventional data
#' cd.run(indata = dat_int)
#' # Run adaptive algorithm for observational data
#' cd.run(indata = dat_obs, adaptive = TRUE)
#'
#' ### Optional: Adjust settings
#' n_node <- ncol(dat)
#'
#' # Run algorithm with a given weight
#' # Careful with this option.
#' weights <- matrix(1, nrow = n_node, ncol = n_node)
#'
#' # Run with adjusted settings
#' cd.run(indata = dat_obs, weights = weights, lambdas.length = 10)
#'
#' # Note: Normally, users do not need to change default settings.
#' }
#' @export
cd.run <- function(indata,
                   weights=NULL,
                   lambdas=NULL,
                   lambdas.length=30,
                   whitelist = NULL,
                   blacklist = NULL,
                   error.tol=0.0001,
                   convLb=0.01,
                   weight.scale=1.0,
                   upperbound = 100.0,
                   alpha = 3,
                   permute = FALSE,
                   adaptive = FALSE) {

  cd_adaptive_run(indata = indata,
                  eor = NULL,
                  weights = weights,
                  lambda_seq = lambdas,
                  fmlam = 0.01,
                  nlam = lambdas.length,
                  whitelist = whitelist,
                  blacklist = blacklist,
                  eps = error.tol,
                  convLb = convLb,
                  qtol = error.tol,
                  gamma = weight.scale,
                  upperbound = upperbound,
                  threshold = alpha,
                  permute = permute,
                  adaptive = adaptive)

}

cd_adaptive_run <- function(indata,
                            eor,
                            weights,
                            lambda_seq,
                            fmlam,
                            nlam,
                            whitelist,
                            blacklist,
                            eps,
                            convLb,
                            qtol,
                            gamma,
                            upperbound,
                            threshold,
                            permute,
                            adaptive)
{
  if (adaptive == FALSE) {
    return(CD_call(indata, eor, permute, weights, lambda_seq, fmlam, nlam, whitelist, blacklist, eps, convLb, qtol, gamma, upperbound, threshold)$fit)
  }
  else {
    cd_call_out <- CD_call(indata, eor, permute, weights, lambda_seq, fmlam, nlam, whitelist, blacklist, eps, convLb, qtol, gamma, upperbound, threshold)
    adaptive_weights <- cd_call_out$adaptive_weights
    return(CD_call(indata, eor, permute, adaptive_weights, lambda_seq = NULL, fmlam, nlam, whitelist, blacklist, eps, convLb, qtol, gamma, upperbound, threshold)$fit)
  }
}

# Convert input to the right form.
CD_call <- function(indata,
                    eor,
                    permute,
                    weights,
                    lambda_seq,
                    fmlam,
                    nlam,
                    whitelist,
                    blacklist,
                    eps,
                    convLb,
                    qtol,
                    gamma,
                    upperbound,
                    threshold) {

  # Allow users to input a data.frame, but kindly warn them about doing this.
  # if the input is a dataframe, the data set is treated as an observational data set.
  # ivn will be initialized to be a list of length dataSize, and every element is 0.
  if(is.data.frame(indata)){
    warning(sparsebnUtils::alg_input_data_frame())
    dataSize <- nrow(indata)
    ivn <- vector("list", length = dataSize)
    ivn <- lapply(ivn, function(x){
      return(c(0L))
    })
    data <- sparsebnUtils::sparsebnData(indata, ivn = ivn, type = "discrete")
  }
  else {
    data <- indata
  }

  # Check data format
  if(!sparsebnUtils::is.sparsebnData(data)) stop(sparsebnUtils::input_not_sparsebnData(data))

  # Extract the data and the intervention list.
  data_matrix <- dat_transform(data)
  data_ivn <- data$ivn
  if (is.null(data_ivn)) {
    data_ivn <- as.list(rep(0L, nrow(data_matrix)))
  }
  if (length(data_ivn)!= nrow(data_matrix)) stop("length of ivn should be equals to number of observations!")

  node_index <- 0:ncol(data_matrix)
  data_level <- data$levels
  data_names <- names(data$data)

  for(i in 1:length(data_ivn)){
    if(is.character(data_ivn[[i]])) {
      data_ivn[[i]] = match(data_ivn[[i]], data_names)
    }
  }

  # Get the dimensions of the data matrix
  dataSize <- nrow(data_matrix)
  node <- ncol(data_matrix)

  # the input data_matrix should be a matrix
  data_matrix <- as.matrix(data_matrix)

  # get n_levels.
  n_levels <- as.integer(sapply(data$levels, function(x){length(x)}))
  if(sum(n_levels<2)) stop("Some nodes has only one level! There must be at least two levels for each node! Remove nodes with one level!")

  # get observational index (obsIndex_R) from interventional list (ivn)
  obsIndex_R <- get_obsIndex(data_ivn, node)

  # make sure that for each observation, at least on node is not under intervention. If all nodes are under intervention, stop and require user to remove that observation.
  ind <- 1:node
  is_obs_zero <- sapply(obsIndex_R, function(x){(length(x)==1 && x == 0)})
  if(length(ind[is_obs_zero])!=0) {
    stop(sprintf("%d th node has been intervened in all observations, remove this node \n", ind[is_obs_zero]))
  }

  # minus 1 from all elements in obsIndex_R to incorporate with C++.
  obsIndex_R <- lapply(obsIndex_R, function(x) {as.integer(x-1)})

  # check/generate eor and eor_nr
  if(is.null(eor)) {
    eor_nr <- node*(node-1)/2
    eor <- matrix(0, nrow=eor_nr, ncol=2)
    cnt1=1
    for (i in 1:(node-1)) {
      for (j in (i+1):node) {
        eor[cnt1, 1] = i;
        eor[cnt1, 2] = j;
        cnt1 = cnt1+1;
      }
    }
  }

  if (permute) {
    eor <- eor[sample(1:eor_nr), ] # run with random order of blocks
  }

  eor_nr <- as.integer(eor_nr)
  eor <- matrix(as.integer(eor), ncol = 2)

  # check/generate weight matrix
  if(is.null(weights)) {
    weights <- matrix(1, node, node)
  }
  if(ncol(weights)!=nrow(weights)) stop("weights should be a square matrix!")
  if(ncol(weights)!=node) stop("wrong dimension for weights, number of colmn of weight matrix should be node!")

  # type conversion for tunning parameters
  node = as.integer(node)
  dataSize = as.integer(dataSize)
  eps = as.numeric(eps)
  convLb = as.numeric(convLb)
  qtol = as.numeric(qtol)
  gamma = as.numeric(gamma)
  upperbound = as.numeric(upperbound)
  threshold = as.integer(threshold)

  # check/generate lambda sequence
  if(is.null(lambda_seq)) {
    lambda_m <- max_lambda(indata,
                           weights,
                           gamma,
                           upperbound)
    lambda_seq <- sparsebnUtils::generate.lambdas(lambda.max = lambda_m, lambdas.ratio = fmlam, lambdas.length = nlam, scale = "log")
  }
  else {
    nlam = length(lambda_seq)
  }
  lambda_seq <- as.numeric(lambda_seq)
  # fmlam = as.numeric(fmlam)
  nlam = as.integer(nlam)

  # white list should be calculated after lambda sequence is generated
  if(!is.null(whitelist)) {
    if(!is.matrix(whitelist) || ncol(whitelist) != 2){
      stop("whitelist must be a matrix with exactly 2 columns!")
    }
    if(any(is.na(whitelist))){
      stop("whitelist cannot have missing values!")
    }
    if(!(all(is.numeric(whitelist))||all(is.character(whitelist)))) {
      stop("whitelis must be a list of all integers or characters!")
    }
    if(all(is.character(whitelist))){
      whitelist <- get_bwlist(whitelist, data_names)
    }
    for(i in 1:nrow(whitelist)) {
      weights[whitelist[i, 1], whitelist[i, 2]] = 0.0
    }
  }

  if(!is.null(blacklist)) {
    if(!is.matrix(blacklist) || ncol(blacklist) != 2){
      stop("blacklist must be a matrix with exactly 2 columns!")
    }
    if(any(is.na(blacklist))){
      stop("blacklist cannot have missing values!")
    }
    if(!(all(is.numeric(blacklist))||all(is.character(blacklist)))) {
      stop("blacklist must be a list of all integers or characters!")
    }
    if(all(is.character(blacklist))){
      blacklist <- get_bwlist(blacklist, data_names)
    }
    for(i in 1:nrow(blacklist)) {
      weights[blacklist[i, 1], blacklist[i, 2]] = upperbound
    }
  }

  # check white and black list
  # check white and black list do not have overlap
  if(!is.null(whitelist) && !is.null(blacklist)) {
    for(i in 1:nrow(whitelist)) {
      if(sum(apply(blacklist, 1, function(x, w_row) {identical(w_row, x)}, whitelist[i, ]))) stop("blacklist and whitelist cannot have common edges!")
    }
  }

  weights <- matrix(as.numeric(weights), ncol = node)

  # run CD algorithm
  estimate <- CD_path(node,
                      dataSize,
                      data_matrix,
                      n_levels,
                      obsIndex_R,
                      eor_nr,
                      eor,
                      lambda_seq,
                      nlam,
                      eps,
                      convLb,
                      qtol,
                      weights,
                      gamma,
                      upperbound,
                      threshold)

  # extract lambdas
  # lambda <- estimate$lambdas
  # extract adjacency matrix
  estimateG <- estimate$estimateG
  time <- estimate$time
  beta_l2 <- estimate$adaptive_weights

  adaptive_weights <- get_adaptWeights(beta_l2)

  # convert each element in fit to sparsebnFit object
  fit <- get.edgeList(estimateG, dataSize, lambda_seq, time)

  # delete null graphs along the solution path due to upper bound on the number of edges. Now the upper bound is fixed, will let the user to decide the number of maximum number of edges in the future
  if_remove <- sapply(fit, function(x) {x$nedge == 0})
  if(if_remove[1]==TRUE) {if_remove[1] = FALSE}
  fit[if_remove] = NULL

  # add node names to output
  for(k in seq_along(fit)){
    fit[[k]] <- append(fit[[k]], list(data_names), after = 1) # insert node names into second slot
    names(fit[[k]])[2] <- "nodes"
    names(fit[[k]]$edges) <- data_names
  }

  # convert element of fit to sparsebnFit object
  fit <- lapply(fit, sparsebnUtils::sparsebnFit)

  # convert fit to sparsebnPath object
  fit <- sparsebnUtils::sparsebnPath(fit)

  return(list(fit = fit, adaptive_weights = adaptive_weights))
  # return(fit)
}

# a function that directly calls from cpp
# type check, no converting type of an input at this point
CD_path <- function(node,
                    dataSize,
                    data_matrix,
                    n_levels,
                    obsIndex_R,
                    eor_nr,
                    eor,
                    lambda_seq,
                    nlam,
                    eps,
                    convLb,
                    qtol,
                    weights,
                    gamma,
                    upperbound,
                    threshold) {
  # check node parameter
  if(!is.integer(node)) stop("node must be a integer!")
  if(node <= 0) stop("node must be a positive integer!")

  # check dataSize parameter
  if(!is.integer(dataSize)) stop("dataSize must be a integer!")
  if(dataSize <= 0) stop("dataSize must be a positive integer!")

  # check data_matrix
  if(node!=ncol(data_matrix) || dataSize!=nrow(data_matrix)) stop("dimension does not match. node should be the number of columns of data matrix, and dataSize should be numbe of rows of data matrix.")
  if(sum(sapply(data_matrix, function(x){!is.integer(x)}))!=0) stop ("data_matrix has to be a data.frame with integer entries!")

  # check n_levels
  if (!is.integer(n_levels)) stop("n_levels must be a vector of integers!")
  if (length(n_levels)!=node) stop("Length of n_levels does not compatible with the input data set. n_levels must be a vector of length equals to the number of node!")
  if(sum(n_levels<2)) stop("Some node has only one level! There must be at least two levels for each node! Remove nodes with one level!")
  max_levels <- sapply(as.data.frame(data_matrix), function(x){length(unique(x))})
  if (sum(max_levels>n_levels)) stop("The number of levels and the data set is not compatible! Check data set and the input ivn list!")

  # check obsIndex_R
  if (!is.list(obsIndex_R)) stop("obsIndex_R must be a list!")
  if (sum(sapply(obsIndex_R, function(x){!is.integer(x)}))!=0) stop("element of obsIndex_R must be a vector of integers!")
  if (sum(sapply(obsIndex_R, function(x){sum(x<0)}))!=0) stop("element of obsIndex_R must be a positive number!")
  if (length(obsIndex_R)!=node) stop("obsIndex_R must be a list of length equals to the number of nodes!")

  # check eor_nr and eor
  if (!is.integer(eor)) stop("eor must be a vector of integers!")
  if (nrow(eor)==0) stop("eor cannot be empty!")
  if (!is.integer(eor_nr)) stop("eor_nr must be an integer!")
  if (eor_nr != nrow(eor)) stop("eor_nr must be the number of rows of eor!")

  # check lambda_seq
  if(!is.numeric(lambda_seq)) stop("lambda_seq must be a numeric vector!")
  if(sum(lambda_seq<0)>0) stop("lambda_seq must be a non-negative numeric vector!")

  # check nlam
  if(!is.integer(nlam)) stop("nlam must be an integer!")
  if(nlam<=0) stop("nlam must be a positive integer!")
  if(nlam!=length(lambda_seq)) stop("nlam must be the length of lambda_seq!")

  # check eps
  if(!is.numeric(eps)) stop("eps must be a numeric number!")
  if(eps<=0) stop("eps must be a positive number!")

  # check convLb
  if (!is.numeric(convLb)) stop("convLb must be a numeric number!")
  if (convLb<=0) stop("convLb must be a positive number!")

  # check qtol
  if (!is.numeric(qtol)) stop("qtol must be a numeric number!")
  if (qtol<=0) stop("qtol must be a positive number!")

  # check weights
  if (!is.numeric(weights)) stop("weights must be a numeric matrix!")
  if (is.integer(weights)) stop("weights must not be an integer matrix!")
  if (ncol(weights)!=node || nrow(weights)!=node) stop("weigths must be a matrix with both number of rows and columns equal to number of variables!")

  # check gamma
  if (!is.numeric(gamma)) stop("gamma must be a numeric number!")
  if (gamma<=0) stop ("gamma must be a positive number!")

  # check upperbound
  if (!is.numeric(upperbound)) stop("upperbound must be a numeric number!")
  if (upperbound <= 0 && upperbound!=-1) stop("upperbound must be a large positive integer to truncate the adaptive weights. Or it can be -1, which indicates that there is no truncation!")

  # check threshold
  if (!is.integer(threshold)) stop("threshold (alpha) must be an integer number!")
  if (threshold <= 0) stop("threshold (alpha must be a positive integer!")

  CD.out <- CD(node,
               dataSize,
               data_matrix,
               n_levels,
               obsIndex_R,
               eor_nr,
               eor,
               lambda_seq,
               nlam,
               eps,
               convLb,
               qtol,
               weights,
               gamma,
               upperbound,
               threshold)

  return(CD.out)
}

