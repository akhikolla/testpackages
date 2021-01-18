################################################################################
# File:             dknn.r
# Created by:       Oleksii Pokotylo
# First published:  
# Last revised:     
# 
# Contains the realization of the Depth-based KNN classifier of Paindaveine and Van Bever (2015).
################################################################################

dknn.train <- function(data, kMax = -1, depth = "halfspace", seed = 0){

  dpth = switch (depth,
                 "halfspace" = 1,
                 "Mahalanobis" = 2,
                 "simplicial" = 3, 0)
  if(dpth == 0)
    stop("Wrong depth: ", depth)
  
  n = nrow(data)
  chunkNumber = n #10
  dimension = ncol(data)-1
  labs <- data[,dimension+1]
  labels <- integer(length(labs))
  
  uniquelab = unique(labs)
  cardinalities = unlist(lapply(uniquelab, function(l)sum(labs == l)))
  uniquelab = uniquelab[order(cardinalities)]
  cardinalities = cardinalities[order(cardinalities)]
  for (i in seq_along(uniquelab)){
    labels[labs == uniquelab[i]] <- i
  }
  if(length(uniquelab)<2)
    stop("There is only one class")
  
  if (!is.numeric(kMax) 
      || is.na(kMax) 
      || length(kMax) != 1 
      || !.is.wholenumber(kMax) 
      || !(kMax >= 1 
           && kMax <= (cardinalities[1]+rev(cardinalities)[1]) 
           || kMax == -1)){
    warning("In treatment number ", i, ": Argument \"k\" not specified 
            correctly. Defaults are applied")
    kMax <- - 1
  }
  if(kMax == -1)    kMax <- n/2

  kMax <- min(kMax, n - 1)
  kMax <- max(kMax, 2)
  
  points <- as.vector(t(data[,1:dimension]))
  k <- as.integer(.C("DKnnLearnCv", 
          as.double(points), 
          as.integer(labels),
          as.integer(n), 
          as.integer(dimension), 
          as.integer(kMax), 
          as.integer(dpth),
          k = integer(1),
          as.integer(chunkNumber),
          as.integer(seed))$k)
  
  dknn <- list(data = data, 
               n = n, 
               dimension = dimension,
               labels = labels,
               uniquelab = uniquelab,
               methodSeparator = "Dknn",
               k = k,
               depth = dpth,#depth,
               seed = seed)
  
  return (dknn)
}

dknn.classify.trained <- function(objects, dknn){
  
  # Correct input data
  if(is.data.frame(objects))
    objects = as.matrix(objects)
  if (!is.matrix(objects)){
    objects <- matrix(objects, nrow=1)
  }
  if(ncol(objects)!=dknn$dimension)
    stop("Parameter 'objects' has wrong dimension")
  
  points <- as.vector(t(dknn$data[,1:dknn$dimension]))
  obj <- as.vector(t(objects))
  output <- .C("DKnnClassify", 
          as.double(obj),
          as.integer(nrow(objects)),
          as.double(points), 
          as.integer(dknn$labels),
          as.integer(dknn$n), 
          as.integer(dknn$dimension), 
          as.integer(dknn$k), 
          as.integer(dknn$depth),
          as.integer(dknn$seed),
          output=integer(nrow(objects)))$output

  results = dknn$uniquelab[output]
  return (results)
}

dknn.classify <- function(objects, data, k, depth = "halfspace", seed = 0){
  
  n = nrow(data)
  dimension = ncol(data)-1
  labs <- data[,dimension+1]
  labels <- integer(length(labs))
  
  uniquelab = unique(labs)
  cardinalities = unlist(lapply(uniquelab, function(l)sum(labs == l)))
  uniquelab = uniquelab[order(cardinalities)]
  cardinalities = cardinalities[order(cardinalities)]
  for (i in seq_along(uniquelab)){
    labels[labs == uniquelab[i]] <- i
  }
  if(length(uniquelab)<2)
    stop("There is only one class")
  
  if (!is.numeric(k) 
      || is.na(k) 
      || length(k) != 1 
      || !.is.wholenumber(k) 
      || k < 1
      || k > nrow(data)){
    stop("Argument \"k\" not specified correctly.")
  }
  
  dpth = switch (depth,
                 "halfspace" = 1,
                 "Mahalanobis" = 2,
                 "simplicial" = 3, 0)
  if(dpth == 0)
    stop("Wrong depth: ", depth)
  
  dknn <- list(data = data, 
               n = n, 
               dimension = dimension,
               labels = labels,
               uniquelab = uniquelab,
               methodSeparator = "Dknn",
               k = k,
               depth = dpth,#depth,
               seed = seed)
  
  return (dknn.classify.trained(objects, dknn))
}