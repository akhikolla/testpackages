################################################################################
# File:             ddalpha-internal.r
# Created by:       Pavlo Mozharovskyi, Oleksii Pokotylo
# First published:  28.02.2013
# Last revised:     20.02.2019
# 
# Contains the internal functions of the DDalpha-classifier.
# 
# For a description of the algorithm, see:
#   Lange, T., Mosler, K. and Mozharovskyi, P. (2014). Fast nonparametric 
#     classification based on data depth. Statistical Papers.
#   Mozharovskyi, P., Mosler, K. and Lange, T. (2015). Classifying real-world 
#     data with the DDalpha-procedure. Mimeo.
################################################################################

.ddalpha.create.structure <- function(formula, data, subset, ...){
  
  # if the user calls ddalpha(data, <other named parameters>, ...)
  # for backward compatibility
  if(!missing(formula) && missing(data) && (is.data.frame(formula) || is.matrix(formula))){
    data = formula
    formula = NULL
  }
  
  # formula is present
  if(!missing(formula) && !is.null(formula)){
    
    needed.frame <- 1 # PM(2018-06-22)

    cl <- match.call(call = sys.call(sys.parent(n = needed.frame)))
    mf <- match.call(expand.dots = FALSE, call = sys.call(sys.parent(n = needed.frame)))
    m <- match(c("formula", "data", "subset"#, "weights", "na.action", "offset"
                 ), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    
    clres = colnames(mf)
    cat("Selected columns: ", paste(clres, collapse = ",   "), "\n")
    classif.formula = delete.response(terms(mf))
    
    data = cbind(mf[,-1,drop=F], mf[,1,drop=F])
    
  #  y = model.response(mf)
  #  mm = model.matrix(formula, data = mf)
  #  mm = as.data.frame(mm[,-1,drop=F])
    
  #  colnames(mm) <- paste0("P_",1:ncol(mm))
  #  data = cbind(mm, y)
    
    #colnames(mf) <- paste0("P_",1:ncol(mm))
                           
  } else {  # no formula
    # Check for data consistency
    if (!(is.matrix(data) && is.numeric(data)
          || is.data.frame(data) && prod(sapply(data[,-ncol(data)], is.numeric)))){
      stop("Argument data has unacceptable format. Classifier can not be trained!!!")
    }
    cl = NULL
    classif.formula = NULL
    clres = colnames(data)
    
    if(!is.data.frame(data))
      data = as.data.frame(data)
    if(!missing(subset))
      data = data[subset,]
  }
  
  names(data)[ncol(data)] <- "CLASS"
  
  # Elemantary statistics
  dimension <- ncol(data) - 1
  numOfPoints <- nrow(data)
  classNames <- unique(data[,dimension + 1])
  numOfClasses <- length(classNames)

  if(!is.data.frame(data))
    data = as.data.frame(data)
  names(data)[ncol(data)] <- "CLASS"
  
  # Creating overall structure
  ddalpha <- list(
         call = cl,
         colnames = clres, # colnames after formula is applied
         classif.formula = classif.formula, # for classification
         raw = data, # after formula is applied
         dimension = dimension, 
         numPatterns = numOfClasses, 
         numPoints = numOfPoints, 
         patterns = list(), 
         needtransform = 0,    # 0 - no transform, 1 - transform new points before classification, all classes the same, 2 - transform differently w.r.t. to classes.
         classifiers = list(), 
#         numClassifiers = 0, 
#         methodDepth = "halfspace", 
#         methodSeparator = "alpha", 
#         methodAggregation = "majority", 
#         methodsOutsider = NULL, 
#         numDirections = 1000, 
#         directions = NULL, 
#         projections = NULL, 
         sameDirections = TRUE, 
         useConvex = FALSE, 
         maxDegree = 3, 
         numChunks = numOfPoints, 
#         knnrange = NULL,
#         mahEstimate = "moment", 
#         mahParMcd = 0.75, 
#         mahPriors = NULL, 
#         knnK = 1, 
#         knnX = NULL, 
#         knnY = NULL, 
#         knnD = NULL, 
         treatments = c("LDA", "KNNAFF", "KNN", "DEPTH.MAHALANOBIS", "RANDEQUAL", "RANDPROP", "IGNORE"))

  # Ordering patterns according to their cardinalities
  classCardinalities <- rep(0, numOfClasses)
  for (i in 1:numOfClasses){
    classCardinalities[i] <- nrow(data[data[,dimension + 1] == classNames[i],])
  }
  # Creating pattern templates
  patterns <- as.list("")
  for (i in 1:numOfClasses){
    maxCarIndex <- which.max(classCardinalities)
    # Creating a single template
    ddalpha$patterns[[i]] <- list(
      index = i, 
      points = data[data[,dimension + 1] == classNames[maxCarIndex],1:dimension], 
      name = classNames[maxCarIndex], 
      cardinality = classCardinalities[maxCarIndex], 
      depths = matrix(rep(0, numOfClasses*classCardinalities[maxCarIndex]), nrow = classCardinalities[maxCarIndex], ncol = numOfClasses), 
      votes = 0#, 
      #         center = 0, 
      #         cov = 0, 
      #         sigma = 0, 
      #         centerMcd = 0, 
      #         covMcd = 0, 
      #         sigmaMcd = 0
    )    

    # Adding pattern template to the list of patterns
    class(ddalpha$patterns[[i]])<-"ddalpha.pattern"
    # Deleting processed pattern
    classCardinalities[maxCarIndex] <- -1    
  }

  return (ddalpha)
}

.check.depth.exists <- function(depth) {
  fname = paste0(".", depth, "_validate")
  f <- try(match.fun(fname), silent = T)
  if (!is.function(f))
    warning(paste0("No validator function: ", fname))
  fname = paste0(".", depth, "_learn")
  
  f <- (match.fun(fname))
  if (!is.function(f))
    stop(paste0("No function: ", fname))
  fname = paste0(".", depth, "_depths")
  f <- (match.fun(fname))
  if (!is.function(f))
    stop(paste0("No function: ", fname))
}

.ddalpha.learn.depth <- function(ddalpha){
  
  # try to find a custom depth
  fname = paste0(".", ddalpha$methodDepth, "_learn")
  f <- try(match.fun(fname), silent = T)
  if (is.function(f)){
    ddalpha = f(ddalpha)
    return(ddalpha)
  }
  
  # If it's the random Tukey depth, compute it first
  if (ddalpha$methodDepth == "halfspace"){
    dSpaceStructure <- .halfspace_space(ddalpha)
    ddalpha$directions <- dSpaceStructure$directions
    ddalpha$projections <- dSpaceStructure$projections
    tmpDSpace <- dSpaceStructure$dspace
  }
  if (ddalpha$methodDepth == "projection"){
    dSpaceStructure <- .projection_space(ddalpha)
    tmpDSpace <- dSpaceStructure$dspace
    if (ddalpha$dmethod == "random"){
      ddalpha$directions <- dSpaceStructure$directions
      ddalpha$projections <- dSpaceStructure$projections
    }
  }
  
  classBegin = 1
  if (ddalpha$methodDepth == "potential"){
    tmpDSpace <- .ddalpha.count.depths(ddalpha, NULL)
  }
  
  # Calculating depths in each pattern
  for (i in 1:ddalpha$numPatterns){
    if (   ddalpha$methodDepth == "halfspace" 
        || ddalpha$methodDepth == "projection"
        || ddalpha$methodDepth == "potential"){
      # Random depth is already calculated, just distribute
      ddalpha$patterns[[i]]$depths <- tmpDSpace[classBegin:(classBegin+ddalpha$patterns[[i]]$cardinality-1),]
      classBegin = classBegin+ddalpha$patterns[[i]]$cardinality
    } else
    if (ddalpha$methodDepth == "zonoid"){
      # Calculate depths for the class w.r.t all classes, saying to which of the classes the chunk belongs
      ddalpha$patterns[[i]]$depths <- .zonoid_depths(ddalpha, ddalpha$patterns[[i]]$points, i)
    } else
    if (ddalpha$methodDepth == "Mahalanobis"){
      ddalpha$patterns[[i]]$depths <- .Mahalanobis_depths(ddalpha, ddalpha$patterns[[i]]$points)
    } else
    if (ddalpha$methodDepth == "spatial"){
      # Calculate depths for the class w.r.t all classes, saying to which of the classes the chunk belongs
      ddalpha$patterns[[i]]$depths <- .spatial_depths(ddalpha, ddalpha$patterns[[i]]$points)
    } else
    if (ddalpha$methodDepth == "spatialLocal"){
      # Calculate depths for the class w.r.t all classes, saying to which of the classes the chunk belongs
      ddalpha$patterns[[i]]$depths <- .spatialLocal_depths(ddalpha, ddalpha$patterns[[i]]$points)
    } else
    if (ddalpha$methodDepth == "simplicial"){
      # Calculate depths for the class w.r.t all classes, saying to which of the classes the chunk belongs
      ddalpha$patterns[[i]]$depths <- .simplicial_depths(ddalpha, ddalpha$patterns[[i]]$points)
    } else
    if (ddalpha$methodDepth == "simplicialVolume"){
      # Calculate depths for the class w.r.t all classes, saying to which of the classes the chunk belongs
      ddalpha$patterns[[i]]$depths <- .simplicialVolume_depths(ddalpha, ddalpha$patterns[[i]]$points)
    } else
      stop("Unknown depth ", ddalpha$methodDepth)
  }

  return (ddalpha)
}

.getFunction <- function(fname) {
  f <- try(match.fun(fname), silent = T)
  if (!is.function(f))
    stop("Wrong or absent function: ", fname)
  f
}

.ddalpha.learn.binary <- function(ddalpha){
  
  fname = paste0(".", ddalpha$methodSeparator, "_learn")
  learn <- try(match.fun(fname), silent = T)
  if (!is.function(learn))
    stop("Wrong or absent function: ", fname)
  
  # Separating (calculating extensions and normals)
  counter <- 1
  # Determining multi-class behaviour
  if (ddalpha$methodAggregation == "majority"){
    for (i in 1:(ddalpha$numPatterns - 1)){
      for (j in (i + 1):ddalpha$numPatterns){
        # Creating a classifier
        classifier <- learn(ddalpha, i, j,
                            ddalpha$patterns[[i]]$depths, 
                            ddalpha$patterns[[j]]$depths)
        
        classifier$index          = counter
        classifier$index1         = i
        classifier$index2         = j
        
        if(class(classifier)=="list")
          class(classifier) <- paste0("ddalpha.", ddalpha$methodSeparator)
        else
          class(classifier) <- c(class(classifier), paste0("ddalpha.", ddalpha$methodSeparator))
        
        # Adding the classifier to the list of classifiers
        ddalpha$classifiers[[counter]] <- classifier
        
        counter <- counter + 1
      }
    }
    ddalpha$numClassifiers <- counter - 1
  }
  if (ddalpha$methodAggregation == "sequent"){
    for (i in 1:ddalpha$numPatterns){
      anotherClass <- NULL
      for (j in 1:ddalpha$numPatterns){
        if (j != i){
          anotherClass <- rbind(anotherClass, ddalpha$patterns[[j]]$depths)
        }
      }
      classifier <- learn(ddalpha, i, -i,
                          ddalpha$patterns[[i]]$depths, 
                          anotherClass)
      
      classifier$index          = counter
      classifier$index1         = i
      classifier$index2         = -i
      
      if(class(classifier)=="list")
        class(classifier) <- paste0("ddalpha.", ddalpha$methodSeparator)
      else
        class(classifier) <- c(class(classifier), paste0("ddalpha.", ddalpha$methodSeparator))
      
      # Adding the classifier to the list of classifiers
      ddalpha$classifiers[[i]] <- classifier
    }
    ddalpha$numClassifiers <- ddalpha$numPatterns
  }
  
  return (ddalpha)
}

.alpha_learn <- function(ddalpha, index1, index2, depths1, depths2){
  points <- as.vector(t(rbind(depths1, depths2)))
  numClass1 <- nrow(depths1)
  numClass2 <- nrow(depths2)
  numPoints <- numClass1 + numClass2
  dimension <- ncol(depths1)
  cardinalities <- c(numClass1, numClass2)
  upToPower <- ddalpha$maxDegree
  minFeatures <- 2
  maxExtDimension <- (factorial(dimension + upToPower) / (factorial(dimension)*factorial(upToPower))) - 1;
  
  p <- .C("AlphaLearnCV", as.double(points), as.integer(numPoints), as.integer(dimension), as.integer(cardinalities),  as.integer(upToPower), as.integer(ddalpha$numChunks), as.integer(minFeatures), as.integer(ddalpha$debug), portrait=double(maxExtDimension + 1))$portrait
  
  degree <- p[1];
  extDimension <- (factorial(dimension + degree) / (factorial(dimension)*factorial(degree))) - 1;
  d <- 0
  for (i in 2:(extDimension + 1)){
    if (p[i] != 0){
      d <- d + 1
    }
  }
  
  return(list(
    hyperplane     = p,
    degree         = degree,
    dimProperties  = length(p) - 1,
    dimFeatures    = d
  ))
}

.alpha_classify <- function(ddalpha, classifier, depths){
  toClassify <- as.double(as.vector(t(depths)))
  m = as.integer(nrow(depths))
  q = as.integer(ncol(depths))
  result <- .C("AlphaClassify", 
               toClassify, 
               m, q, 
               as.integer(classifier$degree), 
               as.double(classifier$hyperplane), 
               output=integer(m))$output
  return(result)
}

.polynomial_learn <- function(ddalpha, index1, index2, depths1, depths2){
  
  polynomial <- .polynomial_learn_C(ddalpha$maxDegree, 
                                    rbind(depths1[,c(index1, index2)], depths2[,c(index1, index2)]),
                                    nrow(depths1), nrow(depths2), 
                                    ddalpha$numChunks, ddalpha$seed)
  return(list(
    polynomial     = polynomial$coefficients,
    degree         = polynomial$degree,
    axis           = polynomial$axis))
}

.polynomial_classify <- function(ddalpha, classifier, depths){

  x = ifelse(classifier$axis == 0, classifier$index1, classifier$index2)
  y = ifelse(classifier$axis == 0, classifier$index2, classifier$index1)
  
  result <- (- depths[,y])
  for (obj in 1:nrow(depths)){        
    val <- depths[obj,x]
    for(j in 1:classifier$degree){
      result[obj] <- result[obj] + classifier$polynomial[j]*val^j
    }
  }
  return(result)
}


.ddalpha.learn.multiclass <- function(ddalpha){
  
  fname = paste0(".", ddalpha$methodSeparator, "_learn")
  learn <- try(match.fun(fname), silent = T)
  if (!is.function(learn))
    stop("Wrong or absent function: ", fname)
  
  classifier <- learn(ddalpha)
  if(class(classifier)=="list")
    class(classifier) <- paste0("ddalpha.", ddalpha$methodSeparator)
  else
    class(classifier) <- c(class(classifier), paste0("ddalpha.", ddalpha$methodSeparator))
      
  ddalpha$classifiers[[1]] <- classifier
  ddalpha$numClassifiers <- 1
  
  return (ddalpha)
}

#.ddalpha.learn.knnlm 
.knnlm_learn <- function(ddalpha){
  
  x <- NULL
  y <- NULL
  for (i in 1:ddalpha$numPatterns){
    x <- rbind(x, ddalpha$patterns[[i]]$depths)
    y <- c(y, rep(i - 1, ddalpha$patterns[[i]]$cardinality))
  }
  x <- as.vector(t(x))
  y <- as.vector(y)

  k <- .C("KnnLearnJK", 
          as.double(x), 
          as.integer(y), 
          as.integer(ddalpha$numPoints), 
          as.integer(ddalpha$numPatterns), 
          as.integer(ddalpha$knnrange), 
          as.integer(2), 
          k=integer(1))$k

  return (list(knnK = k,
               knnX = x,
               knnY = y))
}

.knnlm_classify <- function(ddalpha, classifier, depths){
  z <- as.vector(t(depths))
  output <- .C("KnnClassify", 
               as.double(z), 
               as.integer(nrow(depths)), 
               as.double(classifier$knnX), 
               as.integer(classifier$knnY), 
               as.integer(ddalpha$numPoints), 
               as.integer(ddalpha$numPatterns), 
               as.integer(classifier$knnK), 
               as.integer(2), 
               output=integer(nrow(depths)))$output
  
  return(output+1)
}

.maxD_learn <- function(ddalpha) return(list())

.maxD_classify <- function(ddalpha, classifier, depths) apply(depths, 1, which.max)

.ddalpha.learn.outsiders <- function(ddalpha, methodsOutsider = "LDA", settingsOutsider = NULL){
  # Refine treatments
  if (is.null(settingsOutsider)){
    ddalpha$methodsOutsider <- .parse.methods(methodsOutsider)
  }else{
    ddalpha$methodsOutsider <- .parse.settings(ddalpha, settingsOutsider)
  }
  # Train treatments
  treatments = list(LDA = .lda_learn, QDA = .qda_learn, 
                    KNN = .knn_learn, KNNAff = .knnAff_learn, depth.Mahalanobis = .mah_learn,
                    Ignore = NA, Mark = NA, RandEqual = NA, RandProp = NA)
  for (i in 1:length(ddalpha$methodsOutsider)){
    .treatment = treatments[[ddalpha$methodsOutsider[[i]]$method]]
    if(is.null(.treatment)) stop("Unknown outsiders treatment method ", ddalpha$methodsOutsider[[i]]$method)
    if(!is.function(.treatment)) next; # need no training
    ddalpha$methodsOutsider[[i]] <- .treatment(ddalpha, ddalpha$methodsOutsider[[i]])
  }
  
  return(ddalpha)
}

.ddalpha.count.depths <- function(ddalpha, objects, ...){
  fname = paste0(".", ddalpha$methodDepth, "_depths")
  f <- (match.fun(fname))
  if (!is.function(f))
    stop(paste0("Wrong function for ", ddalpha$methodDepth))
  
  # Count for all data
  if (is.null(objects))  
    for (i in 1:ddalpha$numPatterns){
      objects <- rbind(objects, ddalpha$patterns[[i]]$points)
    }
  else ## if needtransform == 1 the data is already scaled
  # Transform the data once
  if (ddalpha$needtransform == 1){
    objects <- ddalpha$patterns[[1]]$transformer(objects)
  } 
  
  # Calculate depths for all classes together
  if (ddalpha$needtransform != 2){
    res <- f(ddalpha, objects, ...)
    return(res) 
  }
  # Calculate depths w.r.t. each class
  else {
    d <- NULL
    # w.r.t. each class
    for (cls in 1:ddalpha$numPatterns){      
      depth <- f(ddalpha, objects, class = cls, ...)
      
      d = cbind(d, depth)
    }
    
    return(d)    
  }
}

.ddalpha.classify.outsiders<- function (objects, ddalpha, settings){
  if (settings$method == "Ignore"){
    return (.ignore_classify(nrow(objects)))
  }
  if (settings$method == "Mark"){
    return (.right_classify(nrow(objects)))
  }
  if (settings$method == "RandEqual"){
    return (.randequal_classify(nrow(objects), ddalpha))
  }
  if (settings$method == "RandProp"){
    return (.randprop_classify(nrow(objects), ddalpha, settings))
  }
  
  treatments = list(LDA = .lda_classify, QDA = .qda_classify, 
                    KNN = .knn_classify, KNNAff = .knnAff_classify, 
                    depth.Mahalanobis = .mah_classify)

  .treatment = treatments[[settings$method]]
  if(is.null(.treatment)) stop("Unknown outsiders treatment method ", settings$method)
  classified = .treatment(objects, ddalpha, settings)
  return(classified)
}

################################################################################
# Functions used for intermediate calculations and checks are presented below
################################################################################

.are_classifiable <- function(objects, points, cardinalities){
  convexes <- .count_convexes(objects, points, cardinalities)
  return (ifelse(rowSums(convexes)>0,1,0))
}

.count_convexes <- function(objects, points, cardinalities, seed = 0){
  if (is.na(seed)) seed = 0
  x <- as.vector(t(points))
  dimension <- ncol(points)
  numClasses <- length(cardinalities)
  o <- as.vector(t(objects))
  numObjects <- nrow(objects)
  result <- .C("IsInConvexes", 
               as.double(x), 
               as.integer(dimension), 
               as.integer(cardinalities), 
               as.integer(numClasses), 
               as.double(o), 
               as.integer(numObjects), 
               as.integer(seed), 
               isInConvexes=integer(numObjects*numClasses))$isInConvexes
  result <- matrix(result, byrow = T, ncol = numClasses)
  return (result)
}

.halfspace_space <- function(ddalpha){
  points <- NULL
  cardinalities <- NULL
  for (i in 1:ddalpha$numPatterns){
    points <- rbind(points, ddalpha$patterns[[i]]$points)
    cardinalities <- c(cardinalities, ddalpha$patterns[[i]]$cardinality)
  }
  x <- as.vector(t(points))
  c <- as.vector(cardinalities)
  k <- ddalpha$numDirections

  method = ddalpha$dmethod
  
  if (method == 0){
    if (ddalpha$sameDirections){
      rez <- .C("HDSpace", as.double(x), as.integer(ncol(points)), as.integer(c), as.integer(ddalpha$numPatterns), as.integer(k), as.integer(1), as.integer(ddalpha$seed), dspc=double(nrow(points)*ddalpha$numPatterns), dirs=double(k*ncol(points)), prjs=double(k*nrow(points)))
      return (list(dspace=matrix(rez$dspc, nrow=nrow(points), ncol=ddalpha$numPatterns, byrow=TRUE), directions=rez$dirs, projections=rez$prjs))
    }else{
      rez <- .C("HDSpace", as.double(x), as.integer(ncol(points)), as.integer(c), as.integer(ddalpha$numPatterns), as.integer(k), as.integer(0), as.integer(ddalpha$seed), dspc=double(nrow(points)*ddalpha$numPatterns), dirs=double(k*ncol(points)), prjs=double(k*nrow(points)))
      return (list(dspace=matrix(rez$dspc, nrow=nrow(points), ncol=ddalpha$numPatterns, byrow=TRUE), directions=0, projections=0))
    }
  } else 
    if (method %in% 1:3){
      
      ds <- .C("HDepthSpaceEx", 
               as.double(x), 
               as.double(x), 
               as.integer(c), 
               as.integer(length(cardinalities)), 
               as.integer(nrow(points)),  
               as.integer(ncol(points)), 
               as.integer(method), 
               depths=double(nrow(points)*length(cardinalities)))$depths  
      
      return (list(dspace=matrix(ds, nrow=nrow(points), ncol=ddalpha$numPatterns, byrow=F)))
    }
  else 
    stop("wrong choise of the algorithm, method = ", method)
}

.halfspace_depths <- function(ddalpha, objects){
  points <- NULL
  cardinalities <- NULL
  for (i in 1:ddalpha$numPatterns){
    points <- rbind(points, ddalpha$patterns[[i]]$points)
    cardinalities <- c(cardinalities, ddalpha$patterns[[i]]$cardinality)
  }
  x <- as.vector(t(points))
  y <- as.vector(t(objects))
  c <- as.vector(cardinalities)
  k <- ddalpha$numDirections
  method <- ddalpha$dmethod
  if (method == 0){
    if (ddalpha$sameDirections){
      result <- .C("HDepth", 
                   as.double(x), 
                   as.double(y), 
                   as.integer(nrow(objects)), 
                   as.integer(ncol(points)), 
                   as.integer(c), 
                   as.integer(ddalpha$numPatterns), 
                   as.double(ddalpha$directions), 
                   as.double(ddalpha$projections), 
                   as.integer(k), 
                   as.integer(1), 
                   as.integer(ddalpha$seed),
                   depths=double(ddalpha$numPatterns*nrow(objects)))
    }else{
      result <- .C("HDepth", 
                   as.double(x), 
                   as.double(y), 
                   as.integer(nrow(objects)), 
                   as.integer(ncol(points)), 
                   as.integer(c), 
                   as.integer(ddalpha$numPatterns), 
                   dirs=double(k*ncol(points)), 
                   prjs=double(k*nrow(points)),
                   as.integer(k), 
                   as.integer(0), 
                   as.integer(ddalpha$seed),
                   depths=double(ddalpha$numPatterns*nrow(objects)))
    }
  }
  else
    if (method %in% 1:3){
      
      ds <- .C("HDepthSpaceEx", 
               as.double(x), 
               as.double(y), 
               as.integer(c), 
               as.integer(length(cardinalities)), 
               as.integer(nrow(objects)),  
               as.integer(ncol(points)), 
               as.integer(method), 
               depths=double(nrow(objects)*length(cardinalities)))$depths  
      
      depths <- matrix(ds, nrow=nrow(objects), ncol=length(cardinalities), byrow=F)
      return (depths)
  }
  else 
    stop("wrong choise of the algorithm, method = ", method)
  return (matrix(result$depths, nrow=nrow(objects), ncol=ddalpha$numPatterns, byrow=TRUE))
}

.zonoid_depths <- function(ddalpha, objects, ownPattern = 0){
  depths <- NULL
  for (i in 1:ddalpha$numPatterns){
    pattern <- ddalpha$patterns[[i]]$points
    x <- as.vector(t(pattern))
    y <- as.vector(t(objects))
    ds <- .C("ZDepth", as.double(x), as.double(y), as.integer(nrow(pattern)), as.integer(nrow(objects)), 
             as.integer(ncol(pattern)), as.integer(ddalpha$seed), depths=double(nrow(objects)))$depths
    if (i == ownPattern){
      ds <- replace(ds, which(ds < 1/nrow(pattern) - sqrt(.Machine$double.eps)), 1/nrow(pattern))
    }else{
      ds <- replace(ds, which(ds < 1/nrow(pattern) - sqrt(.Machine$double.eps)), 0)
    }
    depths <- cbind(depths, ds)
  }
  return (depths)
}

.estimate_moments <- function(ddalpha){
  for (i in 1:ddalpha$numPatterns){
    if (ddalpha$mahEstimate == "moment"){
      
      ddalpha$patterns[[i]]$center <- colMeans(ddalpha$patterns[[i]]$points)
      ddalpha$patterns[[i]]$cov    <- cov(ddalpha$patterns[[i]]$points)
      try(
        ddalpha$patterns[[i]]$sigma  <- solve(ddalpha$patterns[[i]]$cov)
      )        
    }
    if (ddalpha$mahEstimate == "MCD"){
      try(
        estimate <- covMcd(ddalpha$patterns[[i]]$points, ddalpha$mahParMcd)
      )
      try(
        ddalpha$patterns[[i]]$center <- estimate$center
      )
      try(
        ddalpha$patterns[[i]]$cov    <- estimate$cov
      )
      try(
        ddalpha$patterns[[i]]$sigma  <- solve(estimate$cov)
      )
    }
  }
  return(ddalpha)
}

.Mahalanobis_learn <- function(ddalpha){

  if (is.null(ddalpha$mahPriors)){
    ddalpha$mahPriors <- c()
    for (i in 1:ddalpha$numPatterns){
      ddalpha$mahPriors[i] <- ddalpha$patterns[[i]]$cardinality/ddalpha$numPoints
    }
  }
  
  ddalpha <- .estimate_moments(ddalpha)
  
  for (i in 1:ddalpha$numPatterns){
    ddalpha$patterns[[i]]$depths <- .Mahalanobis_depths(ddalpha, ddalpha$patterns[[i]]$points)
  }
  
  return(ddalpha)
}

.Mahalanobis_depths <- function(ddalpha, objects){
  depths <- NULL
  for (j in 1:ddalpha$numPatterns){
    depths <- cbind(depths, .Mahalanobis_depth (objects, center = ddalpha$patterns[[j]]$center, sigma = ddalpha$patterns[[j]]$sigma))
  }
  return (depths)
}

.Mahalanobis_depth <- function(points, center = colMeans(points), sigma = solve(cov(points))){
  if (is.data.frame(points))
    points <- as.matrix(points, drop = F)
  if(is.vector(points))
    points <- t(as.matrix(points, drop = F))
  if (!is.matrix(points))
    stop("Wrong format of 'points'")
  
  i = 1; step = 200
  d <- NULL
  while (i<=nrow(points)){
    tmp1 <- t(t(points[i:min(i+step, nrow(points)),, drop = F]) - center)

    dd <- diag(tmp1 %*% sigma %*% t(tmp1))
    d <- c(d,1/(1 + dd))
    i = i+1+step
  }

   # d <- 1/(1 + (points - center) %*% sigma %*% t(points - center))
  return (d)
}

.projection_space <- function(ddalpha){
  points <- NULL
  cardinalities <- NULL
  for (i in 1:ddalpha$numPatterns){
    points <- rbind(points, ddalpha$patterns[[i]]$points)
    cardinalities <- c(cardinalities, ddalpha$patterns[[i]]$cardinality)
  }
  if (ddalpha$dmethod == "random"){
    x <- as.vector(t(points))
    y <- as.vector(t(points))
    c <- as.vector(cardinalities)
    k <- ddalpha$numDirections
    result <- .C("ProjectionDepth", 
                 as.double(x), 
                 as.double(y), 
                 as.integer(nrow(points)), 
                 as.integer(ncol(points)), 
                 as.integer(c), 
                 as.integer(ddalpha$numPatterns), 
                 dirs=double(k*ncol(points)), 
                 prjs=double(k*nrow(points)), 
                 as.integer(k), 
                 as.integer(1), 
                 as.integer(ddalpha$seed),
                 dspc=double(ddalpha$numPatterns*nrow(points)))
    return (list(dspace=matrix(result$dspc, nrow=nrow(points), 
                               ncol=ddalpha$numPatterns, byrow=TRUE), 
                 directions=result$dirs, projections=result$prjs))
  }
  if (ddalpha$dmethod == "linearize"){
    depths <- NULL
    for (i in 1:ddalpha$numPatterns){
      ds <- .zdepth(ddalpha$patterns[[i]]$points, points)
      depths <- cbind(depths, ds)
    }
    return (list(dspace=depths))
  }
}

.projection_depths <- function(ddalpha, objects){
  points <- NULL
  cardinalities <- NULL
  for (i in 1:ddalpha$numPatterns){
    points <- rbind(points, ddalpha$patterns[[i]]$points)
    cardinalities <- c(cardinalities, ddalpha$patterns[[i]]$cardinality)
  }
  if (ddalpha$dmethod == "random"){
    x <- as.vector(t(points))
    y <- as.vector(t(objects))
    c <- as.vector(cardinalities)
    k <- ddalpha$numDirections
    result <- .C("ProjectionDepth", 
                 as.double(x), 
                 as.double(y), 
                 as.integer(nrow(objects)), 
                 as.integer(ncol(points)), 
                 as.integer(c), 
                 as.integer(ddalpha$numPatterns), 
                 as.double(ddalpha$directions), 
                 as.double(ddalpha$projections), 
                 as.integer(k), 
                 as.integer(0), 
                 as.integer(ddalpha$seed),
                 depths=double(ddalpha$numPatterns*nrow(objects)))
    return (matrix(result$depths, nrow=nrow(objects), ncol=ddalpha$numPatterns, byrow=TRUE))
  }
  if (ddalpha$dmethod == "linearize"){
    depths <- NULL
    for (i in 1:ddalpha$numPatterns){
      ds <- .zdepth(ddalpha$patterns[[i]]$points, objects)
      depths <- cbind(depths, ds)
    }
    return (depths)
  }
}

.simplicialVolume_depths <- function(ddalpha, objects){
  if (is.data.frame(objects))
    objects = as.matrix(objects)
  depths <- NULL
  for (i in 1:ddalpha$numPatterns){
    pattern <- ddalpha$patterns[[i]]$points
    
    points <- as.vector(t(pattern))
    x <- as.vector(t(objects))
    if (ddalpha$d_useCov == 0){
      covEst <- diag(ncol(pattern))
    } else if (ddalpha$d_useCov == 1){
      covEst <- cov(pattern)
    } else if (ddalpha$d_useCov == 2){
      covEst <- covMcd(pattern, ddalpha$d_parMcd)$cov
    }
    ds <- .C("OjaDepth", 
             as.double(points), 
             as.double(x), 
             as.integer(nrow(pattern)), 
             as.integer(nrow(objects)), 
             as.integer(ncol(pattern)), 
             as.integer(ddalpha$seed),
             as.integer(ddalpha$d_exact),
             as.integer(.longtoint(ddalpha$d_k)),
             as.integer(ddalpha$d_useCov),
             as.double(as.vector(t(covEst))),
             depths=double(nrow(objects)))$depths
    depths <- cbind(depths, ds, deparse.level = 0)
  }  
  return (depths) 
}

.simplicial_depths <- function(ddalpha, objects){
  if (is.data.frame(objects))
    objects = as.matrix(objects)
  depths <- NULL
  for (i in 1:ddalpha$numPatterns){
    pattern <- ddalpha$patterns[[i]]$points
    
    points <- as.vector(t(pattern))
    x <- as.vector(t(objects))
    ds <- .C("SimplicialDepth", 
             as.double(points), 
             as.double(x), 
             as.integer(nrow(pattern)), 
             as.integer(nrow(objects)), 
             as.integer(ncol(pattern)), 
             as.integer(ddalpha$seed),
             as.integer(ddalpha$d_exact),
             as.integer(.longtoint(ddalpha$d_k)),           
             depths=double(nrow(objects)))$depths
    depths <- cbind(depths, ds, deparse.level = 0)
  }  
  return (depths) 
}



.spatial_learn <- function(ddalpha){
  if (ddalpha$mahEstimate == "none"){
    for (i in 1:ddalpha$numPatterns){
        ddalpha$patterns[[i]]$center <- colMeans(ddalpha$patterns[[i]]$points)
        ddalpha$patterns[[i]]$cov <- NA
    }
  }else{
    ddalpha <- .estimate_moments(ddalpha)
  }
  
  for (i in 1:ddalpha$numPatterns){
    ddalpha$patterns[[i]]$depths <- .spatial_depths(ddalpha, ddalpha$patterns[[i]]$points)
  }
  
  return(ddalpha)
}

.spatialLocal_learn <- function(ddalpha){
  
  ddalpha <- .estimate_moments(ddalpha)
  
  for (i in 1:ddalpha$numPatterns){
    ddalpha$patterns[[i]]$depths <- .spatialLocal_depths(ddalpha, ddalpha$patterns[[i]]$points)
  }
  
  return(ddalpha)
}

.spatial_depths <- function(ddalpha, objects){
  if (is.data.frame(objects))
    objects = as.matrix(objects)
  depths <- NULL
  for (i in 1:ddalpha$numPatterns){
    pattern <- ddalpha$patterns[[i]]$points
    mean <- ddalpha$patterns[[i]]$center
    cov <- ddalpha$patterns[[i]]$cov
    #suppressWarnings(
    if(sum(is.na(cov)) == 0){
      cov.eig <- eigen(cov)
      B <- cov.eig$vectors %*% diag(sqrt(cov.eig$values))
      lambda <- solve(B)
    } else{
      lambda = diag(ncol(pattern))
    }#)
    ds <- rep(-1, nrow(objects))
    for (i in 1:nrow(objects)){
      tmp1 <- t(lambda %*% (objects[i,] - t(pattern)))
      tmp1 <- tmp1[which(rowSums(tmp1) != 0),]
      tmp2 <- 1/sqrt(rowSums(tmp1^2))
      ds[i] <- 1 - sqrt(sum((colSums(tmp2*tmp1)/nrow(pattern))^2))
    }
    depths <- cbind(depths, ds)
  }
  return (depths)
}

.spatialLocal_depths <- function(ddalpha, objects){
  depths <- NULL
  for (i in 1:ddalpha$numPatterns){
    depths <- cbind(depths, depth.spatial.local(objects, ddalpha$patterns[[i]]$points, ddalpha$kernel.bandwidth[i]))
  }
  return (depths)
}

.NONE_depths <- function(ddalpha, objects){
  depths <- matrix(0, ncol = ddalpha$numPatterns, nrow = nrow(objects))
  
  return (depths)
}
#==========================================================

.alpha_learnOLD <- function(maxDegree, data, numClass1, numClass2, numChunks, debug = F){
  points <- as.vector(t(data))
  numPoints <- numClass1 + numClass2
  dimension <- ncol(data)
  cardinalities <- c(numClass1, numClass2)
  upToPower <- maxDegree
  minFeatures <- 2
  maxExtDimension <- (factorial(dimension + maxDegree) / (factorial(dimension)*factorial(maxDegree))) - 1;
  
  p <- .C("AlphaLearnCV", as.double(points), as.integer(numPoints), as.integer(dimension), as.integer(cardinalities),  as.integer(upToPower), as.integer(numChunks), as.integer(minFeatures), as.integer(debug), portrait=double(maxExtDimension + 1))$portrait
  
  degree <- p[1];
  extDimension <- (factorial(dimension + degree) / (factorial(dimension)*factorial(degree))) - 1;
  d <- 0
  for (i in 2:(extDimension + 1)){
    if (p[i] != 0){
      d <- d + 1
    }
  }
  return(list(por=p,dim=d,deg=degree))
}

.parse.methods <- function(methods){
  if (!is.vector(methods) || length(methods) <= 0){return(list())}
  methods.refined <- unique(toupper(methods))
  
  treatments.settings <- list()
  counter <- 1
  for (i in 1:length(methods.refined)){
    supported <- FALSE
    if (methods.refined[i] == "LDA"){
      supported <- TRUE
      treatment.settings <- structure(
        list(name = "LDA", 
             method = "LDA", 
             priors = NULL, 
             lda = NULL), 
        .Names = c("name", "method", "priors", "lda"))
    }
    if (methods.refined[i] == "QDA"){
      supported <- TRUE
      treatment.settings <- structure(
        list(name = "QDA", 
             method = "QDA", 
             priors = NULL, 
             qda = NULL))
    }
    if (methods.refined[i] == "KNNAFF"){
      supported <- TRUE
      treatment.settings <- structure(
        list(name = "KNNAff", 
             method = "KNNAff", 
             knnAff.methodAggregation = "majority", 
             knnAff.range = -1, 
             knnAff.k = -1, 
             knnAff.classifiers = NULL), 
        .Names = c("name", "method", "knnAff.methodAggregation", "knnAff.range", "knnAff.k", "knnAff.classifiers"))
    }
    if (methods.refined[i] == "KNN"){
      supported <- TRUE
      treatment.settings <- structure(
        list(name = "KNN", 
             method = "KNN", 
             knn.range = -1, 
             knn.k = -1, 
             knn.train = NULL, 
             knn.cl = NULL), 
        .Names = c("name", "method", "knn.range", "knn.k", "knn.train", "knn.cl"))
    }
    if (methods.refined[i] == "DEPTH.MAHALANOBIS"){
      supported <- TRUE
      treatment.settings <- structure(
        list(name = "depth.Mahalanobis", 
             method = "depth.Mahalanobis", 
             mah.estimate = "moment", 
             priors = NULL, 
             mah.classes = NULL, 
             mcd.alpha = 0.5), 
        .Names = c("name", "method", "mah.estimate", "priors", "mah.classes", "mcd.alpha"))
    }
    if (methods.refined[i] == "RANDEQUAL"){
      supported <- TRUE
      treatment.settings <- structure(
        list(name = "RandEqual", 
             method = "RandEqual"), 
        .Names = c("name", "method"))
    }
    if (methods.refined[i] == "RANDPROP"){
      supported <- TRUE
      treatment.settings <- structure(
        list(name = "RandProp", 
             method = "RandProp", 
             priors = NULL), 
        .Names = c("name", "method", "priors"))
    }
    if (methods.refined[i] == "IGNORE"){
      supported <- TRUE
      treatment.settings <- structure(
        list(name = "Ignore", 
             method = "Ignore"), 
        .Names = c("name", "method"))
    }
    if (supported){
      treatments.settings[[counter]] <- treatment.settings
      counter <- counter + 1
    }
  }
  return(treatments.settings)
}

.parse.settings <- function(ddalpha, settings){
  if (!is.list(settings) || length(settings) <= 0){return(list())}
  treatments.names <- c()
  
  treatments.settings <- list()
  counter <- 1
  for (i in 1:length(settings)){
    supported <- FALSE
    if (!is.list(settings[[i]]) 
        || is.null(settings[[i]]$name) 
        || is.null(settings[[i]]$method)){
      warning("In treatment number ", i, ": The treatment has unacceptable format. The treatment will be ignored")
      next
    }
    if (!is.character(settings[[i]]$method) 
        || length(settings[[i]]$method) != 1 
        || !(toupper(settings[[i]]$method) %in% ddalpha$treatments)){
      warning("In treatment number ", i, ": The method name of the treatment is not acceptable. The treatment will be ignored")
      next
    }
    if (!is.character(settings[[i]]$name) 
        || length(settings[[i]]$name) != 1){
      warning("In treatment number ", i, ": The name of the treatment is not acceptable. The treatment will be ignored")
      next
    }
    if (settings[[i]]$name %in% treatments.names){
      warning("In treatment number ", i, ": Treatment with the name ", settings[[i]]$name, " already exists. The treatment will be ignored")
      next
    }else{
      treatments.names <- c(treatments.names, settings[[i]]$name)
    }
    if (toupper(settings[[i]]$method) == "LDA"){
      tmp.name <- settings[[i]]$name
      if (!is.vector(settings[[i]]$priors, mode = "double") 
          || is.na(min(settings[[i]]$priors)) 
          || length(settings[[i]]$priors) != ddalpha$numPatterns 
          || min(settings[[i]]$priors) <= 0 
          || max(settings[[i]]$priors) <= 0){
        warning("In treatment number ", i, ": Argument \"priors\" not specified correctly. Defaults in the form of class portions are applied")
        tmp.priors <- NULL
      }else{
        tmp.priors <- settings[[i]]$priors/sum(settings[[i]]$priors)
      }
      treatment.settings <- structure(
        list(name = tmp.name, 
             method = "LDA", 
             priors = tmp.priors, 
             lda = NULL), 
        .Names = c("name", "method", "priors", "lda"))
      supported <- TRUE
    }
    if (toupper(settings[[i]]$method) == "KNNAFF"){
      tmp.name <- settings[[i]]$name
      if (!is.character(settings[[i]]$knnAff.methodAggregation) 
          || length(settings[[i]]$knnAff.methodAggregation) != 1 
          || !(settings[[i]]$knnAff.methodAggregation %in% c("majority", "sequent"))){
        warning("In treatment number ", i, ": Argument \"knnAff.methodAggregation\" not specified correctly. \"majority\" is used as a default value")
        tmp.knnAff.methodAggregation <- "majority"
      }else{
        tmp.knnAff.methodAggregation <- settings[[i]]$knnAff.methodAggregation
      }
      if (!is.numeric(settings[[i]]$knnAff.range) 
          || is.na(settings[[i]]$knnAff.range) 
          || length(settings[[i]]$knnAff.range) != 1 
          || !.is.wholenumber(settings[[i]]$knnAff.range) 
          || !(settings[[i]]$knnAff.range >= 2 
               && settings[[i]]$knnAff.range <= (ddalpha$patterns[[ddalpha$numPatterns]]$cardinality + ddalpha$patterns[[ddalpha$numPatterns - 1]]$cardinality - 1) 
               || settings[[i]]$knnAff.range == -1)){
        warning("In treatment number ", i, ": Argument \"knnAff.range\" not specified correctly. Defaults are applied")
        tmp.knnAff.range <- -1
      }else{
        tmp.knnAff.range <- settings[[i]]$knnAff.range
      }
      if (!is.numeric(settings[[i]]$knnAff.k) 
          || is.na(settings[[i]]$knnAff.k) 
          || length(settings[[i]]$knnAff.k) != 1 
          || !.is.wholenumber(settings[[i]]$knnAff.k) 
          || !(settings[[i]]$knnAff.k >= 1 
               && settings[[i]]$knnAff.k <= (ddalpha$patterns[[ddalpha$numPatterns]]$cardinality + ddalpha$patterns[[ddalpha$numPatterns - 1]]$cardinality) 
               || settings[[i]]$knnAff.k == -1)){
        warning("In treatment number ", i, ": Argument \"knnAff.k\" not specified correctly. Defaults are applied")
        tmp.knnAff.k <- -1
      }else{
        tmp.knnAff.k <- settings[[i]]$knnAff.k
      }
      treatment.settings <- structure(
        list(name = tmp.name, 
             method = "KNNAff", 
             knnAff.methodAggregation = tmp.knnAff.methodAggregation, 
             knnAff.range = tmp.knnAff.range, 
             knnAff.k = tmp.knnAff.k, 
             knnAff.classifiers = NULL), 
        .Names = c("name", "method", "knnAff.methodAggregation", "knnAff.range", "knnAff.k", "knnAff.classifiers"))
      supported <- TRUE
    }
    if (toupper(settings[[i]]$method) == "KNN"){
      tmp.name <- settings[[i]]$name
      if (!is.numeric(settings[[i]]$knn.range) 
          || is.na(settings[[i]]$knn.range) 
          || length(settings[[i]]$knn.range) != 1 
          || !.is.wholenumber(settings[[i]]$knn.range) 
          || !(settings[[i]]$knn.range >= 2 
               && settings[[i]]$knn.range <= (ddalpha$numPoints - 1) 
               || settings[[i]]$knn.range == -1)){
        warning("In treatment number ", i, ": Argument \"knn.range\" not specified correctly. Defaults are applied")
        tmp.knn.range <- -1
      }else{
        tmp.knn.range <- settings[[i]]$knn.range
      }
      if (!is.numeric(settings[[i]]$knn.k) 
          || is.na(settings[[i]]$knn.k) 
          || length(settings[[i]]$knn.k) != 1 
          || !.is.wholenumber(settings[[i]]$knn.k) 
          || !(settings[[i]]$knn.k >= 1 
               && settings[[i]]$knn.k <= (ddalpha$numPoints) 
               || settings[[i]]$knn.k == -1)){
        warning("In treatment number ", i, ": Argument \"knn.k\" not specified correctly. Defaults are applied")
        tmp.knn.k <- -1
      }else{
        tmp.knn.k <- settings[[i]]$knn.k
      }
      treatment.settings <- structure(
        list(name = tmp.name, 
             method = "KNN", 
             knn.range = tmp.knn.range, 
             knn.k = tmp.knn.k, 
             knn.train = NULL, 
             knn.cl = NULL), 
        .Names = c("name", "method", "knn.range", "knn.k", "knn.train", "knn.cl"))
      supported <- TRUE
    }
    if (toupper(settings[[i]]$method) == "DEPTH.MAHALANOBIS"){
      tmp.name <- settings[[i]]$name
      if (!is.character(settings[[i]]$mah.estimate) 
          || length(settings[[i]]$mah.estimate) != 1 
          || !(settings[[i]]$mah.estimate %in% c("moment", "MCD"))){
        warning("In treatment number ", i, ": Argument \"mah.estimate\" not specified correctly. \"moment\" is used as a default value")
        tmp.mah.estimate <- "moment"
      }else{
        tmp.mah.estimate <- settings[[i]]$mah.estimate
      }
      if (!is.vector(settings[[i]]$priors, mode = "double") 
          || is.na(min(settings[[i]]$priors)) 
          || length(settings[[i]]$priors) != ddalpha$numPatterns 
          || min(settings[[i]]$priors) <= 0 
          || max(settings[[i]]$priors) <= 0){
        warning("In treatment number ", i, ": Argument \"priors\" not specified correctly. Defaults in the form of class portions are applied")
        tmp.priors <- NULL
      }else{
        tmp.priors <- settings[[i]]$priors/sum(settings[[i]]$priors)
      }
      if (!is.vector(settings[[i]]$mcd.alpha, mode = "double") 
          || is.na(min(settings[[i]]$mcd.alpha)) 
          || length(settings[[i]]$mcd.alpha) != 1 
          || settings[[i]]$mcd.alpha < 0.5 
          || settings[[i]]$mcd.alpha > 1){
        if (tmp.mah.estimate == "MCD"){
          warning("In treatment number ", i, ": Argument \"mcd.alpha\" not specified correctly. 0.75 is used as a default value")
        }
        tmp.mcd.alpha <- 0.75
      }else{
        tmp.mcd.alpha <- settings[[i]]$mcd.alpha
      }
      treatment.settings <- structure(
        list(name = tmp.name, 
             method = "depth.Mahalanobis", 
             mah.estimate = tmp.mah.estimate, 
             priors = tmp.priors, 
             mah.classes = NULL, 
             mcd.alpha = tmp.mcd.alpha), 
        .Names = c("name", "method", "mah.estimate", "priors", "mah.classes", "mcd.alpha"))
      supported <- TRUE
    }
    if (toupper(settings[[i]]$method) == "RANDEQUAL"){
      tmp.name <- settings[[i]]$name
      treatment.settings <- structure(
        list(name = tmp.name, 
             method = "RandEqual"), 
        .Names = c("name", "method"))
      supported <- TRUE
    }
    if (toupper(settings[[i]]$method) == "RANDPROP"){
      tmp.name <- settings[[i]]$name
      if (!is.vector(settings[[i]]$priors, mode = "double") 
          || is.na(min(settings[[i]]$priors)) 
          || length(settings[[i]]$priors) != ddalpha$numPatterns 
          || min(settings[[i]]$priors) <= 0 
          || max(settings[[i]]$priors) <= 0){
        warning("In treatment number ", i, ": Argument \"priors\" not specified correctly. Defaults in the form of class portions are applied")
        tmp.priors <- NULL
      }else{
        tmp.priors <- settings[[i]]$priors/sum(settings[[i]]$priors)
      }
      treatment.settings <- structure(
        list(name = tmp.name, 
             method = "RandProp", 
             priors = NULL), 
        .Names = c("name", "method", "priors"))
      supported <- TRUE
    }
    if (toupper(settings[[i]]$method) == "IGNORE"){
      tmp.name <- settings[[i]]$name
      treatment.settings <- structure(
        list(name = tmp.name, 
             method = "Ignore"), 
        .Names = c("name", "method"))
      supported <- TRUE
    }
    if (supported){
      treatments.settings[[counter]] <- treatment.settings
      counter <- counter + 1
    }
  }

  return (treatments.settings)
}

.lda_learn <- function(ddalpha, settings){
  settings$lda <- MASS::lda(formula=as.formula("CLASS ~ ."), data=ddalpha$raw, priors=settings$priors)
  settings$priors <- settings$lda$prior
  return (settings)
}

.lda_classify <- function(objects, ddalpha, settings){
  if (!is.data.frame(objects)) objects = as.data.frame(objects)
  names(objects) <- names(ddalpha$raw)[1:ncol(objects)]
  return (as.list(predict(settings$lda, objects)$class))
}

.qda_learn <- function(ddalpha, settings){
  settings$qda <- MASS::qda(formula=as.formula("CLASS ~ ."), data=ddalpha$raw, priors=settings$priors)
  settings$priors <- settings$qda$prior
  return (settings)
}

.qda_classify <- function(objects, ddalpha, settings){
  if (!is.data.frame(objects)) objects = as.data.frame(objects)
  names(objects) <- names(ddalpha$raw)[1:ncol(objects)]
  return (as.list(predict(settings$qda, objects)$class))
}

.knnAff_learn <- function(ddalpha, settings){
  counter <- 1
  # Determining multi-class behaviour
  if (settings$knnAff.methodAggregation == "majority"){
    for (i in 1:(ddalpha$numPatterns - 1)){
      for (j in (i + 1):ddalpha$numPatterns){
        # Creating a classifier
        classifier.index          <- counter
        classifier.index1         <- i
        classifier.index2         <- j
        classifier.points         <- as.double(t(rbind(ddalpha$patterns[[i]]$points, ddalpha$patterns[[j]]$points)))
        classifier.cardinalities  <- as.integer(c(ddalpha$patterns[[i]]$cardinality, ddalpha$patterns[[j]]$cardinality))
        if (settings$knnAff.k < 1 || settings$knnAff.k > (ddalpha$patterns[[i]]$cardinality + ddalpha$patterns[[j]]$cardinality - 1))
        {
          if (settings$knnAff.range < 2 || settings$knnAff.range > (ddalpha$patterns[[i]]$cardinality + ddalpha$patterns[[j]]$cardinality - 1)){
            maxk <- 10*( (ddalpha$numPoints)^(1/ddalpha$dimension) ) + 1
          }else{
            maxk <- settings$knnAff.range
          }
          maxk <- min(maxk, ddalpha$patterns[[i]]$cardinality + ddalpha$patterns[[j]]$cardinality - 1)
          maxk <- max(maxk, 2)
          classifier.range <- maxk
          classifier.k <- as.integer(.C("KnnAffInvLearnJK", 
                             classifier.points, 
                             as.integer(ddalpha$dimension), 
                             classifier.cardinalities, 
                             as.integer(maxk), 
                             k=integer(1))$k)
        }else{
          classifier.range <- settings$knnAff.range
          classifier.k <- as.integer(settings$knnAff.k)
        }
        # Adding the classifier to the list of classifiers
        settings$knnAff.classifiers[[counter]] <- 
          list(index = classifier.index, 
               index1 = classifier.index1, 
               index2 = classifier.index2, 
               points = classifier.points, 
               cardinalities = classifier.cardinalities, 
               k = classifier.k, 
               range = classifier.range)
        counter <- counter + 1
      }
    }
  }
  if (settings$knnAff.methodAggregation == "sequent"){
    for (i in 1:ddalpha$numPatterns){
      anotherClass <- NULL
      for (j in 1:ddalpha$numPatterns){
        if (j != i){
          anotherClass <- rbind(anotherClass, ddalpha$patterns[[j]]$points)
        }
      }
      classifier.index          <- counter
      classifier.index1         <- i
      classifier.index2         <- -1
      classifier.points         <- as.double(t(rbind(ddalpha$patterns[[i]]$points, anotherClass)))
      classifier.cardinalities  <- as.integer(c(ddalpha$patterns[[i]]$cardinality, nrow(anotherClass)))
      if (settings$knnAff.k < 1 || settings$knnAff.k > ddalpha$numPoints)
      {
        if (settings$knnAff.range < 2 || settings$knnAff.range > (ddalpha$numPoints - 1)){
          maxk <- 10*( (ddalpha$numPoints)^(1/ddalpha$dimension) ) + 1
        }else{
          maxk <- settings$knnAff.range
        }
        maxk <- min(maxk, ddalpha$numPoints - 1)
        maxk <- max(maxk, 2)
        classifier.range <- maxk
        classifier.k <- as.integer(.C("KnnAffInvLearnJK", 
                           classifier.points, 
                           as.integer(ddalpha$dimension), 
                           classifier.cardinalities, 
                           as.integer(maxk), 
                           k=integer(1))$k)
      }else{
        classifier.range <- settings$knnAff.range
        classifier.k <- as.integer(settings$knnAff.k)
      }
      # Adding the classifier to the list of classifiers
      settings$knnAff.classifiers[[counter]] <- 
        list(index = classifier.index, 
             index1 = classifier.index1, 
             index2 = classifier.index2, 
             points = classifier.points, 
             cardinalities = classifier.cardinalities, 
             k = classifier.k, 
             range = classifier.range)
      counter <- counter + 1
    }
  }
  return (settings)
}

.knnAff_classify <- function(objects, ddalpha, settings){
  # Correct input data
  if (!is.matrix(objects)){
    objects <- matrix(objects, nrow=1)
  }
  # Initialization of the vote array
  votes <- matrix(rep(0, nrow(objects)*ddalpha$numPatterns), nrow=nrow(objects), ncol=ddalpha$numPatterns)
  for (i in 1:length(settings$knnAff.classifiers)){
    res <- .C("KnnAffInvClassify", 
              as.double(t(objects)), 
              as.integer(nrow(objects)), 
              settings$knnAff.classifiers[[i]]$points, 
              as.integer(ddalpha$dimension), 
              settings$knnAff.classifiers[[i]]$cardinalities, 
              settings$knnAff.classifiers[[i]]$k, 
              output=integer(nrow(objects)))$output
    for (j in 1:nrow(objects)){
      if (res[j] == 0){
        votes[j,settings$knnAff.classifiers[[i]]$index1] <- votes[j,settings$knnAff.classifiers[[i]]$index1] + 1
      }else{
        votes[j,settings$knnAff.classifiers[[i]]$index2] <- votes[j,settings$knnAff.classifiers[[i]]$index2] + 1
      }
    }
  }
  # Collect results
  results <- list()
  for (i in 1:nrow(objects)){
    results[[i]] <- ddalpha$patterns[[which.max(votes[i,])]]$name
  }
  return (results)
}

.knn_learn <- function(ddalpha, settings){
  settings$knn.train <- ddalpha$raw[,1:ddalpha$dimension]
  settings$knn.cl <- ddalpha$raw[,ddalpha$dimension + 1]
  if (settings$knn.k < 1 || settings$knn.k > ddalpha$numPoints){
    if (settings$knn.range < 1 || settings$knn.range > ddalpha$numPoints - 1){
      settings$knn.range <- 10*( (ddalpha$numPoints)^(1/ddalpha$dimension) ) + 1
      settings$knn.range <- min(settings$knn.range, ddalpha$numPoints - 1)
      settings$knn.range <- max(settings$knn.range, 2)
    }
    cv.err <- c()
    ks <- 1:settings$knn.range
    for (i in ks){
      newpre <- as.vector(class::knn.cv(settings$knn.train, settings$knn.cl, k = i))
      cv.err <- c(cv.err, sum(settings$knn.cl != newpre))
    }
    settings$knn.k <- ks[which.min(cv.err)] 
  }
  return (settings)
}

.knn_classify <- function(objects, ddalpha, settings){
  return (class::knn(settings$knn.train, objects, settings$knn.cl, settings$knn.k))
}

.mah_learn <- function(ddalpha, settings){
  settings$mah.classes <- list()
  if (is.null(settings$priors)){
    settings$priors <- c()
    for (i in 1:ddalpha$numPatterns){
      settings$priors[i] <- ddalpha$patterns[[i]]$cardinality/ddalpha$numPoints
    }
  }
  for (i in 1:ddalpha$numPatterns){
    class.prior <- NULL
    class.mean <- NULL
    class.cov <- NULL
    class.sigma <- NULL
    class.prior <- settings$priors[i]
    if (settings$mah.estimate == "moment"){
      class.mean <- colMeans(ddalpha$patterns[[i]]$points)
      class.cov <- cov(ddalpha$patterns[[i]]$points)
      class.sigma <- solve(class.cov)
    }
    if (settings$mah.estimate == "MCD"){
      estimate <- robustbase::covMcd(ddalpha$patterns[[i]]$points, alpha=settings$mcd.alpha)
      class.mean <- estimate$center
      class.cov <- estimate$cov
      class.sigma <- solve(class.cov)
    }
    settings$mah.classes[[i]] <- structure(
      list(index = i, 
           prior = class.prior, 
           mean = class.mean, 
           cov = class.cov, 
           sigma = class.sigma), 
      .Names = c("index", "prior", "mean", "cov", "sigma"))
  }
  return (settings)
}

.mah_classify <- function(objects, ddalpha, settings){
  # Correct input data
  if (!is.matrix(objects)){
    objects <- matrix(objects, nrow=1)
  }
  # Initialization of the vote array
  votes <- matrix(rep(0, nrow(objects)*ddalpha$numPatterns), nrow=nrow(objects), ncol=ddalpha$numPatterns)
  for (i in 1:nrow(objects)){
    for (j in 1:length(settings$mah.classes)){
      votes[i,j] <- settings$mah.classes[[j]]$prior*.Mahalanobis_depth(objects[i,], settings$mah.classes[[j]]$mean, settings$mah.classes[[j]]$sigma)
    }
  }
  # Collect results
  results <- list()
  for (i in 1:nrow(objects)){
    results[[i]] <- ddalpha$patterns[[which.max(votes[i,])]]$name
  }
  return (results)
}

.ignore_classify <- function(nobjects){
  return (as.list(rep("Ignored", nobjects)))
}

.right_classify <- function(nobjects){
  return (as.list(rep("Outsider", nobjects)))
}

.randequal_classify <- function(nobjects, ddalpha){
  results <- list()
  for (i in 1:nobjects){
    results[[i]] <- ddalpha$patterns[[sample(1:ddalpha$numPatterns,1)]]$name
  }
  return (results)
}

.randprop_classify <- function(nobjects, ddalpha, settings){
  priors <- settings$priors
  if (is.null(priors)){
    priors <- c()
    for (i in 1:ddalpha$numPatterns){
      priors[i] <- ddalpha$patterns[[i]]$cardinality/ddalpha$numPoints
    }
  }
  results <- list()
  for (i in 1:nobjects){
    results[[i]] <- ddalpha$patterns[[sample(1:ddalpha$numPatterns,1,prob = priors)]]$name
  }
  return (results)
}

# Function is taken from the R-documentation, "Examples" to the function "is.integer"
.is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

summary.ddalpha <- function(object, printseparators = T, ...){
  x = object
  if(!is.null(x$call)){
      cat("Call: ", format(x$call), "\nResulting table: ", paste(x$colnames, collapse = ",   "),"\n")
  }
  
  cat("ddalpha:\n")
  cat("\t num.points = ", x$numPoints, 
      ", dimension = ", x$dimension, 
      ", num.patterns = ", x$numPatterns, "\n", sep="")
  cat("\t depth = \"", x$methodDepth, "\"\n", sep="")
  cat("\t separator = \"", x$methodSeparator, "\"\n", sep="")
  cat("\t aggregation.method = \"", x$methodAggregation, "\"\n", sep="")
  if (is.numeric(x$numChunks)) cat("\t num.chunks =", x$numChunks, "\n")
  if (is.numeric(x$numDirections)) cat("\t num.directions =", x$numDirections, "\n")
  cat("\t use.convex =", x$useConvex, "\n")
  if (is.numeric(x$maxDegree)) cat("\t max.degree =", x$maxDegree, "\n")
  cat("patterns:\n")
  for (i in 1:length(x$patterns)){
    cat("\t ");print(x$patterns[[i]])
  }
  cat("num.classifiers = ", x$numClassifiers, "\n")
  if(x$numClassifiers == 1 || printseparators)
    for (i in 1:x$numClassifiers){
      print(x$classifiers[[i]], prefix = "\t ")
    }
  cat("outsider.methods:\n")
  if (is.null(x$methodsOutsider)){
    cat ("\t Absent\n", sep="")
  }else{
    for (i in 1:length(x$methodsOutsider)){
      cat ("\t \"",x$methodsOutsider[[i]]$name, "\":\n", sep="")
      cat("\t\t method = \"", x$methodsOutsider[[i]]$method, "\"\n", sep="")
      if (   x$methodsOutsider[[i]]$method == "LDA" 
          || x$methodsOutsider[[i]]$method == "RandProp" 
          || x$methodsOutsider[[i]]$method == "depth.Mahalanobis"){
        cat("\t\t priors =", x$methodsOutsider[[i]]$priors, "\n")
      }
      if (x$methodsOutsider[[i]]$method == "KNNAff"){
        cat("\t\t aggregation.method = \"", 
            x$methodsOutsider[[i]]$knnAff.methodAggregation, "\"\n", sep="")
        for (j in 1:length(x$methodsOutsider[[i]]$knnAff.classifiers)){
          cat("\t\t k.range = ", format(paste("1:", 
                                              x$methodsOutsider[[i]]$knnAff.classifiers[[j]]$range, sep=""), 
                                        justify="right", width=10), sep="")
        }
        cat("\n")
        for (j in 1:length(x$methodsOutsider[[i]]$knnAff.classifiers)){
          cat("\t\t k       = ", format( 
            x$methodsOutsider[[i]]$knnAff.classifiers[[j]]$k, 
            justify="right", width=10), sep="")
        }
        cat("\n")
      }
      if (x$methodsOutsider[[i]]$method == "KNN"){
        cat("\t\t k.range = 1:", x$methodsOutsider[[i]]$knn.range, "\n", sep="")
        cat("\t\t k =", x$methodsOutsider[[i]]$knn.k, "\n")
      }
      if (x$methodsOutsider[[i]]$method == "depth.Mahalanobis"){
        cat("\t\t estimate = \"", x$methodsOutsider[[i]]$mah.estimate, "\"\n", sep="")
        if (x$methodsOutsider[[i]]$mah.estimate == "MCD"){
          cat("\t\t mcd.alpha = ", x$methodsOutsider[[i]]$mcd.alpha, "\n")
        }
      }
    }
  }
  invisible(x)
}

print.ddalpha <- function(x, printseparators = F, ...){
  if(!is.null(x$call)){
    cat("Call: ", format(x$call), "\nResulting table: ", paste(x$colnames, collapse = ",   "),"\n")
  }

  cat("depth = \"", x$methodDepth, "\"\n", sep="")
  cat("separator = \"", x$methodSeparator, "\"\n", sep="")
  cat("aggregation.method = \"", x$methodAggregation, "\"\n", sep="")
  cat("patterns:\n")
  for (i in 1:length(x$patterns)){
    cat("\t ");print(x$patterns[[i]])
  }
  cat("num.classifiers = ", x$numClassifiers, "\n")
  if(x$numClassifiers == 1 || printseparators)
    for (i in 1:x$numClassifiers){
      print(x$classifiers[[i]], prefix = "\t ", full = F)
    }
  cat("outsider.methods:\n")
  if (is.null(x$methodsOutsider)){
    cat ("\t Absent\n", sep="")
  }else{
    for (i in 1:length(x$methodsOutsider)){
      cat("\t method = \"", x$methodsOutsider[[i]]$method, "\"\n", sep="")
      if (   x$methodsOutsider[[i]]$method == "LDA" 
             || x$methodsOutsider[[i]]$method == "RandProp" 
             || x$methodsOutsider[[i]]$method == "depth.Mahalanobis"){
        cat("\t priors =", x$methodsOutsider[[i]]$priors, "\n")
      }
      if (x$methodsOutsider[[i]]$method == "KNNAff"){
        cat("\t aggregation.method = \"", 
            x$methodsOutsider[[i]]$knnAff.methodAggregation, "\"\n", sep="")
        for (j in 1:length(x$methodsOutsider[[i]]$knnAff.classifiers)){
          cat("\t k.range = ", format(paste("1:", 
                                              x$methodsOutsider[[i]]$knnAff.classifiers[[j]]$range, sep=""), 
                                        justify="right", width=10), sep="")
        }
        cat("\n")
        for (j in 1:length(x$methodsOutsider[[i]]$knnAff.classifiers)){
          cat("\t k       = ", format( 
            x$methodsOutsider[[i]]$knnAff.classifiers[[j]]$k, 
            justify="right", width=10), sep="")
        }
        cat("\n")
      }
      if (x$methodsOutsider[[i]]$method == "KNN"){
        cat("\t k.range = 1:", x$methodsOutsider[[i]]$knn.range, "\n", sep="")
        cat("\t k =", x$methodsOutsider[[i]]$knn.k, "\n")
      }
      if (x$methodsOutsider[[i]]$method == "depth.Mahalanobis"){
        cat("\t estimate = \"", x$methodsOutsider[[i]]$mah.estimate, "\"\n", sep="")
        if (x$methodsOutsider[[i]]$mah.estimate == "MCD"){
          cat("\t mcd.alpha = ", x$methodsOutsider[[i]]$mcd.alpha, "\n")
        }
      }
    }
  }
  invisible(x)
}

print.ddalpha.pattern <- function(x, ...){
  cat("pattern[", x$index, "]:", sep="")
  cat("\t ", x$cardinality, " points, label = \"", x$name, "\"\n", sep="")
  invisible(x)
}

print.ddalpha.alpha <- function(x, prefix = "", full = T, ...){
  if(full){
    cat(prefix,x$index, ". alpha:\n", sep="")
    cat(prefix,"   degree: ", x$degree,"; axes: ", x$index1, ", ", x$index2, "\n", sep="")
    cat(prefix,"   hyperplane: ", paste(x$hyperplane, collapse = ", "), "\n", sep="")
  } else{
    cat(prefix, x$index, ". degree: ", x$degree,"; axes: ", x$index1, ", ", x$index2, "\n", sep="")
    cat(prefix,"   hyperplane: ", paste(x$hyperplane, collapse = ", "), "\n", sep="")
  }
  invisible(x)
}
print.ddalpha.polynomial <- function(x, prefix = "", full = T,...){
  if(full){
    cat(prefix,x$index, ". polynomial:\n", sep="")
    cat(prefix,"   degree: ", x$degree,"; axes: ", x$index1, ", ", x$index2, ifelse(x$axis!=0, ", invert", ""), "\n", sep="")
    cat(prefix,"   polynomial: ", paste(x$polynomial, collapse = ", "), "\n", sep="")
  } else {
    cat(prefix,x$index, ". degree: ", x$degree,"; axes: ", x$index1, ", ", x$index2, ifelse(x$axis!=0, ", invert", ""), "\n", sep="")
    cat(prefix,"   polynomial: ", paste(x$polynomial, collapse = ", "), "\n", sep="")
  }
  invisible(x)
}
print.ddalpha.knnlm <- function(x, prefix = "", full = T,...){
  if(full)
    cat(prefix,"knnlm:  k = ", x$knnK, "\n", sep="")
  else
    cat(prefix,"k = ", x$knnK, "\n", sep="")
  invisible(x)
}
print.ddalpha.maxD <- function(x, prefix = "", full = T,...){
  if(full)
    cat(prefix,"maximum depth separator\n")
  invisible(x)
}

# .ddalpha.learn.knnlm <- function(ddalpha){
#   
#   # Prepare outputs and distance matrix
#   y <- NULL
#   allPoints <- NULL
#   for (i in 1:ddalpha$numPatterns){
#     allPoints<- rbind(allPoints, ddalpha$patterns[[i]]$depths)
#     y <- c(y, rep(ddalpha$patterns[[i]]$name, ddalpha$patterns[[i]]$cardinality))
#   }
#   #  print(y)
#   dists <- knn.dist(allPoints, dist.meth="maximum")
#   #  plot(ddalpha$patterns[[1]]$depths, col = "red", xlim=c(0,1), ylim=c(0,1))
#   #  points(ddalpha$patterns[[2]]$depths, col = "blue", xlim=c(0,1), ylim=c(0,1))
#   # Cross-validate knn
#   cvErr  <- c()
#   allIndices <- 1:ddalpha$numPoints
#   krange <- 10*( (ddalpha$numPoints)^(1/ddalpha$numPatterns) ) + 1
#   krange <- min(krange, ceiling(ddalpha$numPoints/2))
#   if (ddalpha$numPatterns == 2){krange <- min(krange, 50)}
#   krange <- max(krange, 2)
#   #  cat("Range: 1:", krange, ".\n", sep="")
#   for (i in 1:krange){
#     curPreErr <- 0
#     for (j in 0:(ddalpha$numChunks - 1)){
#       testSel <- allIndices[allIndices%%ddalpha$numChunks == j]
#       test <- allIndices[testSel]
#       train <- allIndices[-test]
#       #      cat("1. Train: ", train, ", test: ", test, ".\n")
#       curPreErr <- curPreErr + sum(knn.predict(train, test, y, dists, k=i, agg.meth="majority", ties.meth="first") != y[test])
#     }
#     cvErr <- c(cvErr, curPreErr)
#     #    cat("i = ", i, "done.\n")
#     #    print(cvErr)
#   }
#   # Collect results
#   ddalpha$knnK <- (1:krange)[which.min(cvErr)]
#   ddalpha$knnX <- allPoints
#   ddalpha$knnY <- y
#   ddalpha$knnD <- dists
#   
#   #  print(ddalpha$knnK)
#   
#   return (ddalpha)
# }
