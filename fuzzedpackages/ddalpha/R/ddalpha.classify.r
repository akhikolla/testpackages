################################################################################
# File:             ddalpha.classify.r
# Created by:       Pavlo Mozharovskyi
# First published:  28.02.2013
# Last revised:     15.05.2013
# 
# Contains the classification function of the DDalpha-classifier.
# 
# For a description of the algorithm, see:
#   Lange, T., Mosler, K. and Mozharovskyi, P. (2012). Fast nonparametric 
#     classification based on data depth. Statistical Papers.
#   Mozharovskyi, P., Mosler, K. and Lange, T. (2013). Classifying real-world 
#     data with the DDalpha-procedure. Mimeo.
################################################################################

ddalpha.classify <- function(ddalpha, 
                             objects,
                             subset,
                             outsider.method = NULL, 
                             use.convex = NULL){
  # Checks
  if (!is.matrix(objects) && !is.data.frame(objects)){
    objects <- matrix(objects, nrow=1)
  }
  if (!(is.matrix(objects) && is.numeric(objects)
        || is.data.frame(objects) && prod(sapply(objects, is.numeric)))){
    warning("Argument \"objects\" has unacceptable format. Classification can not be performed!!!")
    return (NULL)
  } 
  
  # convert using formula
  if(!is.null(ddalpha$classif.formula)){
    objects = model.frame(ddalpha$classif.formula, data = objects)
  }
  
  if(!missing(subset))
    objects = objects[subset,]
  
  if (ncol(objects) != ddalpha$dimension){
    warning("Dimension of the objects to be classified does not correspond to the dimension of the trained classifier. Classification can not be performed!!!")
    return (NULL)
  }
  
  if (ddalpha$methodSeparator == "Dknn")
    return(dknn.classify.trained(objects, ddalpha))
  
#  if (!is.character(outsider.method) 
#      || length(outsider.method) != 1){
#    warning("Argument \"outsidet.method\" not specified correctly. Outsiders will be ignored!!!")
#    outsider.method <- NULL
#  }
  if (is.null(use.convex)){
    use.convex <- ddalpha$useConvex
  }
  depths <- matrix(nrow=0, ncol=ddalpha$numPatterns) #?
  freePoints <- matrix(nrow=0, ncol=ncol(objects)) #?
  
  if (is.null(ddalpha$methodDepth)){ #use only outsiders treatment
    classifiableIndices <- c()
    resultsDepths <- list()
    freePoints <- objects
  } 
  else {
  # Define points that can be classified by the DD-Alpha and the outsiders
  if (use.convex){
    points <- ddalpha$patterns[[1]]$points
    cardinalities <- c(ddalpha$patterns[[1]]$cardinality)
    for (i in 2:ddalpha$numPatterns){
      points <- rbind(points, ddalpha$patterns[[i]]$points)
      cardinalities <- c(cardinalities, ddalpha$patterns[[i]]$cardinality)
    }
    classifiable <- .are_classifiable(objects, points, cardinalities)
    classifiableIndices <- which(classifiable == 1)
    if (length(classifiableIndices) == 0){
      depths <- matrix(nrow=0, ncol=ddalpha$numPatterns)
      freePoints <- objects
    }else{
      depths <- .ddalpha.count.depths(ddalpha, objects[classifiableIndices,])
      
      freePoints <- matrix(objects[-classifiableIndices,,drop=F], 
                           nrow=nrow(objects)-length(classifiableIndices))
    }
  }else{
    ifelse(ddalpha$methodDepth == "ddplot",
      depths <- objects,
      depths <- .ddalpha.count.depths(ddalpha, objects) )

    classifiableIndices <- c()
    for (i in 1:nrow(depths)){
      if (sum(depths[i,]) > 0){
        classifiableIndices <- c(classifiableIndices, i)
      }
    }
    if (length(classifiableIndices) == 0){
      depths <- matrix(nrow=0, ncol=ddalpha$numPatterns)
      freePoints <- objects
    }else{
      depths <- suppressWarnings( as.matrix(depths[classifiableIndices,,drop=F], 
                       nrow=length(classifiableIndices), ncol=ddalpha$numPatterns))
      freePoints <- # objects[-classifiableIndices,,drop=F]# 
                suppressWarnings( as.matrix(objects[-classifiableIndices,,drop=F], 
                           nrow=nrow(objects)-length(classifiableIndices), 
                           ncol=ncol(objects)))
    }
  }
  
  # Classify with the pure DD classifiers
  resultsDepths <- list()
  if (nrow(depths) > 0){

    fname = paste0(".", ddalpha$methodSeparator, "_classify")
    classify <- .getFunction(fname)

    resultsDepths1 <- list()
    if (ddalpha$methodSeparatorBinary){
      #### Binary classifiers
      
      votes <- matrix(rep(0, nrow(depths)*ddalpha$numPatterns), nrow=nrow(depths), ncol=ddalpha$numPatterns)
      for (i in 1:ddalpha$numClassifiers){
        xAxis <- ddalpha$classifiers[[i]]$index1
        yAxis <- ddalpha$classifiers[[i]]$index2
        
        result <- classify(ddalpha, ddalpha$classifiers[[i]], depths)
        
        for (obj in 1:nrow(depths)){
          if (result[obj] > 0){
            votes[obj,xAxis] <- votes[obj,xAxis] + 1
          }else{
            votes[obj,yAxis] <- votes[obj,yAxis] + 1
          }
        }
      }
      for (i in 1:nrow(depths)){
        resultsDepths[[i]] <- ddalpha$patterns[[which.max(votes[i,])]]$name
      }
    } else {
      #### Multiclass classifiers
      
      indexes <- classify(ddalpha, ddalpha$classifiers[[1]], depths)
      
      for (i in 1:nrow(depths)){
        resultsDepths[[i]] <- ddalpha$patterns[[indexes[i]]]$name
      }
    }

  }
  } # end if(!is.null(ddalpha$methodDepth))
  
  # Classify Outsiders
  resultsOutsiders <- as.list(rep("Ignored", nrow(freePoints)))
  freePoints
  if (is.null(outsider.method) && length(ddalpha$methodsOutsider) == 1)
    outsider.method = ddalpha$methodsOutsider[[1]]$name
  if (length(resultsOutsiders) > 0 && !is.null(outsider.method)){
    for (i in 1:length(ddalpha$methodsOutsider)){
      if (toupper(ddalpha$methodsOutsider[[i]]$name) == toupper(outsider.method)){
        resultsOutsiders <- .ddalpha.classify.outsiders(freePoints, ddalpha, ddalpha$methodsOutsider[[i]])
        break
      }
    }
  }
  
  # Merge classifiable and outsiders
  if (length(resultsOutsiders) == 0)
    results <- resultsDepths
  else if(length(resultsDepths) == 0)
    results <- resultsOutsiders
  else{
    
    if(is.factor(resultsOutsiders[[1]]) && !is.factor(resultsDepths[[1]]))
      resultsOutsiders = lapply(resultsOutsiders, as.character)
    
    if(is.numeric(resultsDepths[[1]]) && !is.numeric(resultsOutsiders[[1]]))
      resultsOutsiders = as.numeric(resultsOutsiders)
    
    results <- list()
    counterDepths <- 1
    counterOutsiders <- 1
    for (i in 1:nrow(objects)){
      if (i %in% classifiableIndices){
        results[[i]] <- resultsDepths[[counterDepths]]
        counterDepths <- counterDepths + 1
      }else{
        results[[i]] <- resultsOutsiders[[counterOutsiders]]
        counterOutsiders <- counterOutsiders + 1
      }
    }
  }

  if (length(results) == 1) return(results[[1]])
                       else return (results)
}

predict.ddalpha <- function(object, 
                            objects, 
                            subset,
                            outsider.method = NULL, 
                            use.convex = NULL, ...){
  return(ddalpha.classify(object, objects, subset, outsider.method, use.convex))
}
