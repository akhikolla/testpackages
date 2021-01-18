estLocScale <- function(X, 
                        type = "wrap", precScale = 1e-12,
                        center = TRUE,
                        alpha = 0.5,
                        nLocScale = 25000,
                        silent = FALSE) {
  # Estimates location and scale for every column in X
  #
  # args:
  #   X: data matrix
  #   nLocScale: if < n, then rows are sampled tocalculate loc/scale
  #   type: "wrap" or "mcd", the location/scale estimators used
  #   precScale: precision scale used throughout the algorithm
  #   alpha = value of alpha for the unimcd, h = ceil(alpha * n)
  # Returns: 
  #   loc: the locations of the columns in X
  #   scale: the scales of the columns in X
  
  # Check inputs
  
  # The random seed is retained when leaving the function
  if (exists(".Random.seed",envir = .GlobalEnv,inherits = FALSE)) {
    seed.keep <- get(".Random.seed", envir = .GlobalEnv, 
                     inherits = FALSE)
    on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
  }
  set.seed(0)
  
  if (!is.data.frame(X) & !is.matrix(X) & !is.vector(X)) {
    stop("The input data must be a vector, matrix or a data frame")
  }
  type <- match(type, c("1stepM","hubhub","wrap","mcd", "rawmcd", "wrapmedmad")) - 1
  if (is.na(type)) {
    stop(paste("Invalid \"type\" argument. Should be \"wrap\", \"mcd\" or \"1stepM\""))
  }
  
  if (is.na(alpha)) {
    alpha <- 0.5
  }
  
  if (nLocScale == 0) {
    nLocScale = dim(X)[1]
  }
  
  # Estimate location/scale
  res <- tryCatch( .Call('_cellWise_estLocScale_cpp', as.matrix(X), nLocScale,
                         type,  precScale,
                         center, alpha, PACKAGE = 'cellWise'),
                   "std::range_error" = function(e){
                     conditionMessage( e ) })
  zeroscales <- which(res$scale <= precScale)
  if (!silent) {
    if ( length(zeroscales) > 0) {
      warning(paste(length(zeroscales)," out of ", dim(X)[2], " variables have an estimated scale <= 
\"precScale\" = ", precScale, "."))
    }
  }
  
  loc.out <- drop(res$loc)
  scale.out <- drop(res$scale)
  names(loc.out) <- names(scale.out) <- colnames(X)
  return(list(loc = loc.out, scale = scale.out))
}