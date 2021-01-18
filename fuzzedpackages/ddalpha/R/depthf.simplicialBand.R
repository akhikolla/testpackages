depthf.simplicialBand <- function(objectsf, dataf, modified = TRUE, J = NULL, 
                                  range = NULL, d = 101){
  # Calculate the simplicial band depth
  # Follows the article:
  # Lopez-Pintado, Sun, Lin, Genton (2014). 
  # "Simplicial band depth for multivariate data", 
  # Advances in Data Analysis and Classification, 8(3), 321-338.
  # Args:
  #   objectsf: functoins for which the depth should be computed; a list 
  #             containing lists (functions) of two vectors of equal length, 
  #             named "args" and "vals": arguments sorted in ascending order 
  #             and corresponding them values respectively.
  #   dataf:    data sample of functoins w.r.t. which the depth should be 
  #             computed; structure as for "objectsf".
  #   modified: whether modified simplicial band depth should be computed; 
  #             logical, TRUE by default.
  #   J:        how many functions to consider in each tuple of the 
  #             U-statistics; integer, d+1 by default.
  #   range:    The common range of the domain where the functions of objectsf 
  #             and dataf are observed. Vector of length 2 with the left and 
  #             the right end of the interval. Must contain all arguments given 
  #             in objectsf and dataf.
  #   d:        Grid size to which all the functional data are transformed. For 
  #             depth computation, all functional observations are first 
  #             transformed into vectors of their functional values of 
  #             length d corresponding to equi-spaced points in the domain 
  #             given by the interval range. Functional values in these points 
  #             are reconstructed using linear interpolation, and extrapolation.
  # Returns:
  #   A vector of depths of each of "objectsf" w.r.t. "dataf".
  
  # Check input data for consistency:
  if (length(objectsf) < 1){
    stop("Number of functions for which the depth should be computed is < 1.")
  }
  #if (J < 2){
  #  stop("Impossible to calculate depth with 'J' < 2.")
  #}
  p <- ifelse(is.null(dim(objectsf[[1]]$vals)), 1, dim(objectsf[[1]]$vals)[2])
  J <- p + 1
  if (length(dataf) < J){
    stop("Number of functions w.r.t. which the depth should be computed 
         is < dimension + 1.")
  }
  m <- length(objectsf)
  n <- length(dataf)
  objArgs <- unique(unlist(lapply(objectsf, function(x){return(x$args)})))
  datArgs <- unique(unlist(lapply(dataf, function(x){return(x$args)})))
  if (length(objArgs) != length(datArgs) || 
      sum(objArgs == datArgs) != length(objArgs)){
    stop("Not the same arguments for 'objectsf' and 'dataf'.")
  }
  l <- length(objArgs)
  numObjArgs <- unlist(lapply(objectsf, function(x){return(length(x$args))}))
  if (sum(numObjArgs == length(objArgs)) != m){
    stop("Not the same arguments for all functions in 'objectsf'.")
  }
  numDatArgs <- unlist(lapply(dataf, function(x){return(length(x$args))}))
  if (sum(numDatArgs == length(datArgs)) != n){
    stop("Not the same arguments for all functions in 'dataf'.")
  }
  if (p == 1){
    numObjVals <- unlist(lapply(objectsf, function(x){return(length(x$vals))}))
    numDatVals <- unlist(lapply(dataf, function(x){return(length(x$vals))}))
  }else{
    numObjVals <- unlist(lapply(objectsf, function(x){return(nrow(x$vals))}))
    numDatVals <- unlist(lapply(dataf, function(x){return(nrow(x$vals))}))
  }
  if (sum(numObjVals == numObjArgs) != m){
    stop("Number of arguments and values for (some) functions in 'objectsf' differ.")
  }
  if (sum(numDatVals == numDatArgs) != n){
    stop("Number of arguments and values for (some) functions in 'dataf' differ.")
  }
  # Interpolate the data and prepare it for the C++ transfer
  A <- dataf2rawfd(objectsf, range = range, d = d)
  B <- dataf2rawfd(dataf, range = range, d = d)
  At <- apply(A, 1, function(x) t(x))
  Bt <- apply(B, 1, function(x) t(x))
  fArgs <- approx(dataf[[1]]$args, n = d)$x
  # fArgs <- objectsf[[1]]$args
  #fObjVals <- unlist(lapply(A, function(x){return(t(x$vals))}))
  #fDatVals <- unlist(lapply(B, function(x){return(t(x$vals))}))
  # Call the C++ routine
  #print(as.double(fObjVals))
  #print(as.double(fDatVals))
  #print(as.integer(m))
  #print(as.integer(n))
  #print(as.integer(l))
  #print(as.integer(d))
  ds <- .C("SimplicialBandDepthF", 
           as.double(At),
           as.double(Bt), 
           as.double(fArgs), 
           as.integer(m), 
           as.integer(n), 
           as.integer(d), 
           as.integer(p), 
           as.integer(modified), 
           as.integer(J), 
           depths = double(m))$depths
  return(ds)
}
