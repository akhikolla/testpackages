################################################################################
# File:             ddalpha.train.r
# Created by:       Pavlo Mozharovskyi, Oleksii Pokotylo
# First published:  28.02.2013
# Last revised:     20.02.2019
# 
# Contains the training function of the DDalpha-classifier.
# 
# For a description of the algorithm, see:
#   Lange, T., Mosler, K. and Mozharovskyi, P. (2012). Fast nonparametric 
#     classification based on data depth. Statistical Papers.
#   Mozharovskyi, P., Mosler, K. and Lange, T. (2013). Classifying real-world 
#     data with the DDalpha-procedure. Mimeo.
################################################################################

ddalpha.train <- function(formula, data, subset,
                          depth = "halfspace", 
                          separator = "alpha", 
                          outsider.methods = "LDA", 
                          outsider.settings = NULL, 
                          aggregation.method = "majority",
                          pretransform = NULL,
                          use.convex = FALSE,
                          seed = 0,
                          
                          ...
                          
                          #                           # knn
                          #                           knnrange = NULL, 
                          #                           # alpha
                          #                           num.chunks = 10, 
                          #                           max.degree = 3, 
                          #                          
#                           # halfspace depth
                          #                           num.directions = 1000,                          
                          #                           # Mahalanobis depth
                          #                           mah.estimate = "moment", 
                          #                           mah.parMcd = 0.75, 
                          #                           mah.priors = NULL
                          
){
  
 
  # Raw data processing
  ddalpha <- .ddalpha.create.structure(formula, data, subset, ...) 

  # Check for data consistency #1 (moved to .ddalpha.create.structure)
  #if (!(is.matrix(data) && is.numeric(data)
  #    || is.data.frame(data) && prod(sapply(data[,-ncol(data)], is.numeric)))){
  #  stop("Argument data has unacceptable format. Classifier can not be trained!!!")
  #}
  
  # Check for data consistency #2
  if (ddalpha$numPatterns < 2){
    stop("Number of patterns is < 2. Classifier can not be trained!!!")
  }
  
  # TODO ddalpha$numPatterns should be restricted from above as well
  if (ddalpha$dimension < 2){
    stop("Data dimension is < 2. Classifier can not be trained!!!")
  }
  
  if(separator == "Dknn"){
    dknn = dknn.train(ddalpha$raw, depth = depth, seed = seed, ...)
    dknn$call = ddalpha$call
    dknn$colnames = ddalpha$colnames
    dknn$classif.formula = ddalpha$classif.formula
    return(dknn)
  }
  
  for (i in 1:ddalpha$numPatterns){
    if (ddalpha$patterns[[i]]$cardinality < ddalpha$dimension + 1){
      stop("At least in one patern number of the points < (dimension + 1). Classifier can not be trained!!!")
    }
  }
  if (seed != 0){
    set.seed(seed)
  }
  ddalpha$seed = seed
  
  ## Validating the properties
  depthsThatNeedNoScaling = c("zonoid", "halfspace", "Mahalanobis", "projection", "spatial", "spatialLocal", "simplicial", "simplicialVolume", "ddplot") # note: "spatialLocal" thansforms the data inside, by itself
  supportedDepths = c(depthsThatNeedNoScaling, "potential")
  if (is.null(depth) || toupper(depth) %in% c("", "NULL", "NONE")){
      ddalpha$methodDepth <- NULL
      warning("Argument \"depth\" specified as NULL.")    
  } else
  if (!is.character(depth)
      || length(depth) != 1){
    stop("Argument \"depth\" not specified correctly.")
  } else
  if(!(depth %in% supportedDepths)){
    .check.depth.exists(depth)
    ddalpha$methodDepth <- depth
  }else{
    ddalpha$methodDepth <- depth
  }
  if (!is.character(separator)
      || length(separator) != 1){
    stop("Argument \"separator\" not specified correctly.")
  } else
  if(!(separator %in% c("alpha", "polynomial", "knnlm", "maxD"))){
    fname = paste0(".", separator, "_validate")
    f <- try(match.fun(fname), silent = T)
    if (!is.function(f))
      warning(paste0("No validator function: ", fname, ". Cannot set 'methodSeparatorBinary'."))
    fname = paste0(".", separator, "_learn")
    f <- (match.fun(fname))
    if (!is.function(f))
      stop(paste0("No function: ", fname))
    fname = paste0(".", separator, "_classify")
    f <- (match.fun(fname))
    if (!is.function(f))
      stop(paste0("No function: ", fname))
    ddalpha$methodSeparator <- separator
  }else{
    ddalpha$methodSeparator <- separator
  }
  if (!is.character(aggregation.method)
      || length(aggregation.method) != 1
      || !(aggregation.method %in% c("majority", "sequent", "none"))){
    ddalpha$methodAggregation <- "majority"
    warning("Argument \"aggregation.method\" not specified correctly. \"majority\" is used as a default value")
  }else{
    ddalpha$methodAggregation <- aggregation.method
    
    ddalpha$methodSeparatorBinary = (aggregation.method !="none")
  }
  
  ddalpha$needtransform = 0
  if (!is.null(pretransform))
    if (ddalpha$methodDepth %in% depthsThatNeedNoScaling){
      warning("The used depth method is affine-invariant and pretransform doesn't influence the result. The data won't be transformed.")
    } else if (pretransform == "1Mom" || pretransform == "1MCD"){
      ddalpha$needtransform = 1
      
      if (pretransform == "1Mom")
        mm <- mah.moment(data[,-ncol(data)])
      else   # "1MCD"
        mm <- mah.mcd(data[,-ncol(data)], .mah.parMcd.fromDots(...))              
      
      for (i in 1:ddalpha$numPatterns){
        ddalpha$patterns[[i]]$transformer <- MahMomentTransformer(mm$mu, mm$b)
        ddalpha$patterns[[i]]$points <- ddalpha$patterns[[i]]$transformer(ddalpha$patterns[[i]]$points)
      }
    } else if (pretransform == "NMom" || pretransform == "NMCD"){
      ddalpha$needtransform = 2
      
      for (i in 1:ddalpha$numPatterns){
        if (pretransform == "NMom")
          mm <- mah.moment(ddalpha$patterns[[i]]$points)
        else   # "NMCD"
          mm <- mah.mcd(ddalpha$patterns[[i]]$points, .mah.parMcd.fromDots(...))
        
        ddalpha$patterns[[i]]$transformer <- MahMomentTransformer(mm$mu, mm$b)
      }
    } else {
      warning("Argument pretransform has unacceptable format. The data won't be transformed.")
    }  
  
  # appends ddalpha with the values from the given list (lst)
  ddalpha.append <- function(lst){ if(is.list(lst)) for(k in names(lst))  ddalpha[[k]] <<- lst[[k]] }
  # calls the validation method for the selected separator || depth
  # NO errors or warnings if the function doesn't exist!!
  # ddalpha is READONLY inside the validators
  validate <- function(method_name){
    f <- try(match.fun(paste0(".", method_name, "_validate")), silent = T)
    if (is.function(f)){
      lst  <- f(ddalpha, ...)
      ddalpha.append(lst) 
    }
  }
  
  ## Separator parameters validation
  
  if (!is.null(ddalpha$methodDepth))  
  validate(ddalpha$methodSeparator)
  
  ## Depth parameters validation
  
  if (!is.logical(use.convex) 
      || length(use.convex) != 1 
      || !(use.convex %in% c(TRUE, FALSE))){
    ddalpha$useConvex <- FALSE
    warning("Argument \"use.convex\" not specified correctly. FALSE is used as a default value")
  }else{
    ddalpha$useConvex <- use.convex
  }
  
  if (!is.null(ddalpha$methodDepth))
  validate(ddalpha$methodDepth)
  
  ## The learning procedures
  
  if (!is.null(ddalpha$methodDepth)){
    # Calculate depths
    if(ddalpha$methodDepth == "ddplot"){
      for(i in 1:ddalpha$numPatterns)
        ddalpha$patterns[[i]]$depths <- ddalpha$patterns[[i]]$points
    } else{
      ddalpha <- .ddalpha.learn.depth(ddalpha)
      if(is.null(ddalpha))
        stop("The depth method did not return the 'ddalpha' object.")
    }
    
    # Learn classification machine
    if (ddalpha$methodSeparatorBinary){
      ddalpha <- .ddalpha.learn.binary(ddalpha)
    } else {
      ddalpha <- .ddalpha.learn.multiclass(ddalpha)
    }
    if(is.null(ddalpha))
      stop("The separator method did not return the 'ddalpha' object.")
    
    # if (ddalpha$methodSeparator == "alpha"){
    #   ddalpha <- .ddalpha.learn.alpha(ddalpha)
    # } else
    # if (ddalpha$methodSeparator == "polynomial"){
    #   ddalpha <- .ddalpha.learn.polynomial(ddalpha)
    # } else
    # if (ddalpha$methodSeparator == "knnlm"){
    #   ddalpha <- .ddalpha.learn.knnlm(ddalpha)
    # } else
    #   stop("Define custom classifier")
  } else {
    ddalpha$numClassifiers = 0
  }
  
  # Learn outsider treatments if needed
  if (is.null(ddalpha$methodDepth) || !(ddalpha$methodDepth %in% 
          c("Mahalanobis", "projection", "spatial", "simplicialVolume", "potential"))){#, "simplicial" (may obtain too small values)
    ddalpha <- .ddalpha.learn.outsiders(ddalpha = ddalpha, 
                                        methodsOutsider = outsider.methods, 
                                        settingsOutsider = outsider.settings)
  }
  class(ddalpha) <- "ddalpha"
  
  return (ddalpha)
}

################################################################################
# Validation functions
################################################################################

.mah.parMcd.fromDots <- function(mah.parMcd = 0.75, ...) {
  return(mah.parMcd)
}

.alpha_validate  <- function(ddalpha, num.chunks = 10, max.degree = 3, debug = F,...){
  
  if (ddalpha$methodAggregation == "majority"){
    maxChunks <- ddalpha$patterns[[ddalpha$numPatterns]]$cardinality + ddalpha$patterns[[ddalpha$numPatterns - 1]]$cardinality
  }else{
    maxChunks <- ddalpha$numPoints
  }
  
  if (is.character(num.chunks) && toupper(num.chunks)=="MAX")
    num.chunks <- maxChunks
  else if (!is.numeric(num.chunks) 
           || is.na(num.chunks) 
           || length(num.chunks) != 1 
           || !.is.wholenumber(num.chunks) 
           || !(num.chunks > 0 && num.chunks <= maxChunks)){
    if (!missing(num.chunks))
      warning("Argument \"num.chunks\" not specified correctly. ", maxChunks, " is used instead")
    num.chunks <- maxChunks
  }
  
  if(!is.numeric(max.degree) 
     || is.na(max.degree) 
     || length(max.degree) != 1 
     || !.is.wholenumber(max.degree) 
     || !(max.degree %in% 1:10)){
    max.degree <- 3
    warning("Argument \"max.degree\" not specified correctly. 3 is used as a default value")
  }
  
  return (list(numChunks = num.chunks, maxDegree = max.degree, debug = (debug == T), methodSeparatorBinary = T))  
}

.polynomial_validate  <- .alpha_validate  # the same  

.knnlm_validate  <- function(ddalpha, knnrange = 10*( (ddalpha$numPoints)^(1/ddalpha$numPatterns) ) + 1,...){
  isnull = missing(knnrange) || is.null(knnrange)
  
  if (is.character(knnrange) && toupper(knnrange)=="MAX")
    knnrange = ceiling(ddalpha$numPoints/2)
  else if(is.null(knnrange)
          || !is.numeric(knnrange) 
          || is.na(knnrange) 
          || length(knnrange) != 1 
          || !.is.wholenumber(knnrange) 
          || !(knnrange >=2 && knnrange <= ceiling(ddalpha$numPoints/2))){    
    knnrange <- 10*( (ddalpha$numPoints)^(1/ddalpha$numPatterns) ) + 1   # Default
    
    knnrange <- min(knnrange, ceiling(ddalpha$numPoints/2))
    knnrange <- max(knnrange, 2)
    
    if (!isnull) warning("Argument \"knnrange\" not specified correctly. ", knnrange, " is used instead")
  }
  return (list(knnrange = knnrange, methodSeparatorBinary = F))
}

.maxD_validate  <- function(ddalpha,...){
  return(list(methodSeparatorBinary = F))
}

.ddplot_validate <- function(ddalpha, ...){
  if(ddalpha$dimension!=ddalpha$numPatterns)
    stop("You must pass a DD-plot, with the number of columns equal to the number of classes as data.")
}

.halfspace_validate  <- function(ddalpha, exact, method, num.directions = 1000,...){
  method = .parse_HSD_pars(exact, method)
  if(method == 0)
    if(!is.numeric(num.directions) 
     || is.na(num.directions) 
     || length(num.directions) != 1 
     || !.is.wholenumber(num.directions) 
     || !(num.directions > 1 && num.directions < 10000000) ){
    num.directions <- 1000
    warning("Argument \"num.directions\" not specified correctly. 1000 is used as a default value")
  }
  return (list(dmethod = method, numDirections = num.directions))
}

.projection_validate <- function(ddalpha, method = "random", num.directions = 1000,...){
  if (!(method %in% c("random","linearize")))
    stop("Wrong method")
  if(method == "random")
  if(!is.numeric(num.directions) 
     || is.na(num.directions) 
     || length(num.directions) != 1 
     || !.is.wholenumber(num.directions) 
     || !(num.directions > 1 && num.directions < 10000000) ){
    num.directions <- 1000
    warning("Argument \"num.directions\" not specified correctly. 1000 is used as a default value")
  }
  return (list(dmethod = method, numDirections = num.directions))
}

.simplicial_validate <- function(ddalpha, exact = F, k = 0.05, ...){
  if (exact)
    return(list(d_exact = exact))
  
  if (k <= 0) stop("k must be positive")
  else if (k < 1) k = choose(ddalpha$numPoints, ddalpha$dimension)*k
  
  return(list(d_exact = exact, d_k = k))
}

.simplicialVolume_validate <- function(ddalpha, exact = F, k = 0.05, mah.estimate = "moment", mah.parMcd = 0.75, ...){
  if (toupper(mah.estimate) == "NONE"){
    useCov <- 0
  } else if (toupper(mah.estimate) == "MOMENT"){
    useCov <- 1
  } else if (toupper(mah.estimate) == "MCD"){
    useCov <- 2
  } else {stop("Wrong argument \"mah.estimate\", should be one of \"moment\", \"MCD\", \"none\"")}
  
  if (exact){
    return(list(d_exact = exact, d_k = k, d_useCov = useCov, 
                d_parMcd = mah.parMcd))
  }
  
  if (k <= 0) stop("k must be positive")
  else if (k < 1) k = choose(ddalpha$numPoints, ddalpha$dimension)*k
  
  return(list(d_exact = exact, d_k = k, d_useCov = useCov, 
              d_parMcd = mah.parMcd))
}

.Mahalanobis_validate  <- function(ddalpha, mah.estimate = "moment", mah.priors = NULL, mah.parMcd = 0.75, ...){ 
  if (!is.character(mah.estimate) 
      || length(mah.estimate) != 1 
      || !(mah.estimate %in% c("moment", "MCD"))){
    mah.estimate <- "moment"
    warning("Argument \"mah.estimate\" not specified correctly. \"moment\" is used as a default value")
  }
  
  if (!is.vector(mah.priors, mode = "double") 
      || is.na(min(mah.priors)) 
      || length(mah.priors) != ddalpha$numPatterns 
      || min(mah.priors) <= 0 
      || max(mah.priors) <= 0){
    if (!is.null(mah.priors)){
      warning("Argument \"mah.priors\" not specified correctly. Defaults in the form of class portions are applied")
    }
    mah.priors <- NULL
  }else{
    mah.priors <- mah.priors/sum(mah.priors)
    
  }
  
  ret <- list(mahEstimate = mah.estimate, mahPriors = mah.priors)
  
  if (mah.estimate == "MCD"){
    if (is.null(mah.parMcd) 
        || !is.vector(mah.parMcd, mode = "double") 
        || is.na(min(mah.parMcd)) 
        || length(mah.parMcd) != 1 
        || mah.parMcd < 0.5 
        || mah.parMcd > 1){
      mah.parMcd <- 0.75
      warning("Argument \"mah.parMcd\" not specified correctly. 0.75 is used as a default value")
    }
    ret$mahParMcd = mah.parMcd
  }
  return (ret)
}

.spatial_validate <- function(ddalpha, mah.estimate = "moment", mah.parMcd = 0.75, ...){ 
  if(mah.estimate == "none")
    return(list(mahEstimate = "none"))
  
  return(.Mahalanobis_validate(ddalpha, mah.estimate = mah.estimate, mah.parMcd = mah.parMcd, ...))
}

.spatialLocal_validate  <- function(ddalpha, kernel = "GKernel", kernel.bandwidth = 1, ...)
{  
  # validate paraameters mah.estimate, mah.parMcd
  spatial_val = .spatial_validate(ddalpha, ...)
  
  if (is.null(kernel)
      || suppressWarnings (
        !((kernel  %in% .potentialKernelTypes) || !is.na(as.numeric(kernel))&&(as.numeric(kernel) %in% c(1:length(.potentialKernelTypes))))
      )){
    stop("Argument \"Kernel\" has invaid format.")
  }  
  if (is.null(kernel.bandwidth) 
      || !is.numeric(kernel.bandwidth)){
    stop("Argument \"kernel.bandwidth\" has invaid format.")
  }
  if (length(kernel.bandwidth) == 1){
    if (length(kernel.bandwidth) || is.na(kernel.bandwidth) || kernel.bandwidth == 0){
      stop("Argument \"kernel.bandwidth\" is Zero or NA.")
    } 
    kernel.bandwidth = rep(kernel.bandwidth, ddalpha$numPatterns)
  }
  else {
    if (sum(!is.na(kernel.bandwidth)) != ddalpha$numPatterns || sum(kernel.bandwidth != 0) != ddalpha$numPatterns){
      stop("Argument \"kernel.bandwidth\" has invaid length, Zero or NA elements.")
    }
     
    
    # order bandwidths the same as the classes
    names = sapply(ddalpha$patterns, FUN=function(X) X$name)
    kernel.bandwidth = kernel.bandwidth[order(names)]
  }
  
  spatial_val$kernel = kernel
  spatial_val$kernel.bandwidth = kernel.bandwidth
  
  return(spatial_val)
}

.potential_validate  <- function(ddalpha, kernel = "GKernel", kernel.bandwidth = NULL, ignoreself = FALSE, ...)
{  
  # if kernel.bandwidth is a vector - the values are given in the alphabetical order of the classes nemes
  if (ddalpha$needtransform == 0)
    stop("'pretransform' must be set for depth = 'potential'")
  
  if (is.null(kernel)
      || suppressWarnings (
        !((kernel  %in% .potentialKernelTypes) || !is.na(as.numeric(kernel))&&(as.numeric(kernel) %in% c(1:length(.potentialKernelTypes))))
        )){
    stop("Argument \"Kernel\" has invaid format.")
  }  
  if (is.null(kernel.bandwidth)) {    # use the rule of thumb
    if (ddalpha$needtransform == 2){
      kernel.bandwidth = sapply(ddalpha$patterns, FUN=function(X) nrow(X$points)) ^ (-2/(ddalpha$dimension+4))
    } else {
      kernel.bandwidth = ddalpha$numPoints ^ (-2/(ddalpha$dimension+4))
    }
  }
  else{
  if (#is.null(kernel.bandwidth) ||
        !is.numeric(kernel.bandwidth)
      ||!(is.vector(kernel.bandwidth) || is.list(kernel.bandwidth))){
    stop("Argument \"kernel.bandwidth\" has invaid format.")
  }
  
  if (ddalpha$needtransform == 2){
    if (length(kernel.bandwidth) == 1)
      kernel.bandwidth = rep(kernel.bandwidth, ddalpha$numPatterns)
    if (sum(!is.na(kernel.bandwidth)) != ddalpha$numPatterns || sum(kernel.bandwidth != 0) != ddalpha$numPatterns){
      stop("Argument \"kernel.bandwidth\" has invaid length, Zero or NA elements.")
    } 
    
    # order bandwidths the same as the classes
    names = sapply(ddalpha$patterns, FUN=function(X) X$name)
    kernel.bandwidth = kernel.bandwidth[order(names)]
  } else if (length(kernel.bandwidth) != 1 || is.na(kernel.bandwidth) || kernel.bandwidth == 0){
    stop("Argument \"kernel.bandwidth\" has invaid length, Zero or NA elements.")
  }  
  }
  
  if (is.null(ignoreself) || !is.logical(ignoreself))
    warning ("Argument \"ignoreself\" has invaid format. FALSE used.")
  
  return(list("kernel" = kernel, "kernel.bandwidth" = kernel.bandwidth, "ignoreself" = ignoreself))
}

