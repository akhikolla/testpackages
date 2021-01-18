#' Solve PU problem with lasso or group lasso penalty.
#' 
#' Fit a model using PUlasso algorithm over a regularization path. The regularization path is computed at a grid of values for the regularization parameter lambda. 
#' 
#'@importFrom Rcpp evalCpp
#'@importFrom methods as
#'@import Matrix
#'@useDynLib PUlasso
#'@param X Input matrix; each row is an observation. Can be a matrix or a sparse matrix.
#'@param z Response vector representing whether an observation is labeled or unlabeled.
#'@param py1 True prevalence Pr(Y=1)
#'@param initial_coef A vector representing an initial point where we start PUlasso algorithm from.
#'@param group A vector representing grouping of the coefficients. For the least ambiguity, it is recommended if group is provided in the form of vector of consecutive ascending integers.
#'@param penalty penalty to be applied to the model. Default is sqrt(group size) for each of the group.
#'@param lambda A user supplied sequence of lambda values. If unspecified, the function automatically generates its own lambda sequence based on nlambda and lambdaMinRatio.
#'@param nlambda The number of lambda values. 
#'@param lambdaMinRatio Smallest value for lambda, as a fraction of lambda.max which leads to the intercept only model.
#'@param maxit Maximum number of iterations. 
#'@param maxit_inner Maximum number of iterations for a quadratic sub-problem for CD.
#'@param weights observation weights. Default is 1 for each observation.
#'@param eps Convergence threshold for the outer loop. The algorithm iterates until the maximum change in coefficients is less than eps in the outer loop.
#'@param inner_eps Convergence threshold for the inner loop. The algorithm iterates until the maximum change in coefficients is less than eps in the inner loop.
#'@param verbose A logical value. if TRUE, the function prints out the fitting process.
#'@param stepSize A step size for gradient-based optimization. if NULL, a step size is taken to be stepSizeAdj/mean(Li) where Li is a Lipschitz constant for ith sample
#'@param stepSizeAdjustment A step size adjustment. By default, adjustment is 1 for GD and SGD, 1/8 for SVRG and 1/16 for SAG.
#'@param batchSize A batch size. Default is 1.
#'@param updateFrequency An update frequency of full gradient for method =="SVRG"
#'@param samplingProbabilities sampling probabilities for each of samples for stochastic gradient-based optimization. if NULL, each sample is chosen proportionally to Li.
#'@param method Optimization method. Default is Coordinate Descent. CD for Coordinate Descent, GD for Gradient Descent, SGD for Stochastic Gradient Descent, SVRG for Stochastic Variance Reduction Gradient, SAG for Stochastic Averaging Gradient.
#'@param trace An option for saving intermediate quantities. All intermediate standardized-scale parameter estimates(trace=="param"), objective function values at each iteration(trace=="fVal"), or both(trace=="all") are saved in optResult. Since this is computationally very heavy, it should be only used for decently small-sized dataset and small maxit. A default is "none".
#'@return coef A p by length(lambda) matrix of coefficients
#'@return std_coef A p by length(lambda) matrix of coefficients in a standardized scale
#'@return lambda The actual sequence of lambda values used.
#'@return nullDev Null deviance defined to be 2*(logLik_sat -logLik_null)
#'@return deviance Deviance defined to be 2*(logLik_sat -logLik(model))
#'@return optResult A list containing the result of the optimization. fValues, subGradients contain objective function values and subgradient vectors at each lambda value. If trace = TRUE, corresponding intermediate quantities are saved as well.
#'@return iters Number of iterations(EM updates) if method = "CD". Number of steps taken otherwise.
#'@examples
#'data("simulPU")
#'fit<-grpPUlasso(X=simulPU$X,z=simulPU$z,py1=simulPU$truePY1)
#'@export
#'
grpPUlasso <- function(X,
                       z,
                       py1,
                       initial_coef = NULL,
                       group = 1:ncol(X),
                       penalty = NULL,
                       lambda = NULL,
                       nlambda = 100,
                       lambdaMinRatio = ifelse(N < p, 0.05, 0.005),
                       maxit = ifelse(method == "CD", 1000, N * 10),
                       maxit_inner = 100000,
                       weights = NULL,
                       eps = 1e-04,
                       inner_eps = 1e-02,
                       verbose = FALSE,
                       stepSize = NULL,
                       stepSizeAdjustment = NULL,
                       batchSize = 1,
                       updateFrequency = N,
                       samplingProbabilities = NULL,
                       method = c("CD", "GD", "SGD", "SVRG", "SAG"),
                       trace = c("none", "param", "fVal", "all")
)
{
  N = nrow(X); p = ncol(X)
  # X, z input check
  input_check(X, z, group, penalty, stepSize, samplingProbabilities, weights)
  if (is.null(colnames(X))) {
    colnames(X) <- paste("V", 1:ncol(X), sep = "")
  }
  
  # Reorder samples of X,z
  row_ordering = order(z, decreasing = T)
  col_ordering = order(group)
  ordering_res = ordering_data(row_ordering, col_ordering, X, z, group, weights)
  X_lu = ordering_res$X_lu
  z_lu = ordering_res$z_lu
  w_lu = ordering_res$w_lu
  
  group = ordering_res$group
  group0 = ordering_res$group0
  
  remove(X, z, ordering_res)
  
  # Check the type of X
  is.sparse = FALSE
  if (inherits(X_lu, "sparseMatrix")) {
    is.sparse = TRUE
    X_lu = as(X_lu, "CsparseMatrix")
    X_lu = as(X_lu, "dgCMatrix")
  } else if (inherits(X_lu, "dgeMatrix")) {
    X_lu = as.matrix(X_lu)
  }
  if (!(class(X_lu) == "matrix" ||
        class(X_lu) == "dgCMatrix")) {
    stop("X must be a matrix or a sparse matrix")
  }
  if (typeof(X_lu) != "double") {
    X_lu <- X_lu + 0.0
  } # Ensure type of X is double
  
  # Normalize weights
  if (!is.null(w_lu)) {
    weiOption <- TRUE
    w_lu <- w_lu / sum(w_lu) * length(w_lu)
  } else{
    weiOption <- FALSE
    w_lu <- rep(1, N)
  }
  
  # Apply strong set screening if p >N
  usestrongSet = ifelse(N < p, FALSE, TRUE)
  
  # Match arguments
  method = match.arg(method, choices = c("CD", "GD", "SGD", "SVRG", "SAG"))
  trace = match.arg(trace, choices = c("none", "param", "fVal", "all"))
  
  # Fitting and opt_option setup
  fitting_ls = fitting_setup(
    py1 = py1,
    lambda = lambda,
    lambdaMinRatio = lambdaMinRatio,
    nlambda = nlambda,
    initial_coef = initial_coef,
    group = group,
    penalty = penalty,
    p = p
  )
  
  opt_ls = opt_option_setup(
    method = method,
    trace = trace,
    stepSize = stepSize,
    stepSizeAdjustment = stepSizeAdjustment,
    samplingProbabilities = samplingProbabilities
  )
  
  if(weiOption&&
     (method != "CD")) {
    opt_ls$method = "CD"
    message("Currently the weight option is available for method == CD. Method switched to CD")
  }
  skip_fitting = getOption('PUlasso.skip_fitting')
  
  if(!is.sparse) {
    g <- LU_dense_cpp(
      X_ = X_lu,
      z_ = z_lu,
      icoef_ = fitting_ls$icoef,
      gsize_ = fitting_ls$gsize,
      pen_ = fitting_ls$pen,
      lambdaseq_ = fitting_ls$lambdaseq,
      user_lambdaseq_ = fitting_ls$user_lambdaseq,
      pathLength_ = nlambda,
      lambdaMinRatio_ = lambdaMinRatio,
      pi_ = py1,
      max_nUpdates_ = maxit,
      maxit_ = maxit_inner,
      wei_ = w_lu,
      weiOption_ = weiOption,
      tol_ = eps,
      inner_tol_ = inner_eps,
      useStrongSet_ = usestrongSet,
      verbose_ = verbose,
      stepSize_ = opt_ls$stepSize,
      stepSizeAdj_ = opt_ls$stepSizeAdjustment,
      batchSize_ = batchSize,
      updateFreq_ = updateFrequency,
      samplingProbabilities_ = opt_ls$samplingProbabilities,
      useLipschitz_ = opt_ls$use_Lipschitz_for_ss_or_sProb,
      method_ = method,
      trace_ = opt_ls$trace,
      skipFitting_ = skip_fitting
    )
  } else{
    g <- LU_sparse_cpp(
      X_ = X_lu,
      z_ = z_lu,
      icoef_ = fitting_ls$icoef,
      gsize_ = fitting_ls$gsize,
      pen_ = fitting_ls$pen,
      lambdaseq_ = fitting_ls$lambdaseq,
      user_lambdaseq_ = fitting_ls$user_lambdaseq,
      pathLength_ = nlambda,
      lambdaMinRatio_ = lambdaMinRatio,
      pi_ = py1,
      max_nUpdates_ = maxit,
      maxit_ = maxit_inner,
      wei_ = w_lu,
      weiOption_ = weiOption,
      tol_ = eps,
      inner_tol_ = inner_eps,
      useStrongSet_ = usestrongSet,
      verbose_ = verbose,
      stepSize_ = opt_ls$stepSize,
      stepSizeAdj_ = opt_ls$stepSizeAdjustment,
      batchSize_ = batchSize,
      updateFreq_ = updateFrequency,
      samplingProbabilities_ = opt_ls$samplingProbabilities,
      useLipschitz_ = opt_ls$use_Lipschitz_for_ss_or_sProb,
      method_ = method,
      trace_ = opt_ls$trace,
      skipFitting_ = skip_fitting
    )
  }
  
  
  cpp_results = summary_cpp_results(g, method, trace, colnames = colnames(X_lu), group0 =
                                      group0)
  
  optResult = list(
    method = method,
    convergence = g$convFlag,
    fValues = g$fVals,
    subGradients = g$subgrads,
    stepSize = g$stepSize,
    samplingProbabilities = g$samplingProbabilities,
    std_coef_all = cpp_results$std_coef_all,
    fValues_all = cpp_results$fVals_all,
    maxit = maxit
    
  )
  
  # warning
  if(method %in% c("CD","GD")){
    widx<-which(g$convFlag==1)
    if(length(widx)>0){
      for(i in 1:length(widx)){
        warning(paste("convergence failed at ",widx[i],"th lambda, ", cpp_results$iters[widx[i]],"th iterations",sep=""))
      }
    }
  }else{
    if(verbose){
      widx<-which(g$convFlag==0)
      if(length(widx)>0){
        for(i in 1:length(widx)){
          cat('|param.diff| < eps at',widx[i],'th lambda,', cpp_results$iters[widx[i]],'th iterations\n')
        }
      }
    }
  }
  
  result <- structure(
    list(
      coef = cpp_results$coef,
      std_coef = cpp_results$std_coef,
      lambda = g$lambda,
      nullDev = g$nullDev,
      deviance = g$deviance,
      optResult = optResult,
      iters = cpp_results$iters,
      call = match.call()
    ),
    class = "PUfit"
  )
  
  return(result)
}
