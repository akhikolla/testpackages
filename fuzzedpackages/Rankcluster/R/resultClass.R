#' @title Constructor of Output class
#'
#' @description This class contains a result of a run. Let K be the total number of cluster,
#'  p the number of dimension m the p-vector containing the size of each dimension.
#'
#' @slot proportion a K-vector of proportions.
#' @slot pi a K*p-matrix composed of the scale parameters.
#' @slot mu a matrix with K lines and sum(m) columns in which line k is composed of the location
#'   parameters of cluster k.
#' @slot ll the estimated log-likelihood.
#' @slot bic the estimated BIC criterion.
#' @slot icl the estimated ICL criterion.
#' @slot tik a n*K-matrix containing the estimation of the conditional probabilities for the observed
#' ranks to belong to each cluster.
#' @slot partition a n-vector containing the partition estimation resulting from the clustering.
#' @slot entropy a n*2-matrix containing for each observation its estimated cluster and its entropy. The entropy
#'   output illustrates the confidence in the clustering of each observation (a high entropy means a
#'      low confidence in the clustering)..
#' @slot probability a n*2-matrix similar to the entropy output, containing for each observation
#' its estimated cluster and its probability p(xi; mk, pk) given its
#' cluster. This probability is estimated using the last simulation of the presentation orders used
#' for the likelihood approximation. The probability output exhibits the best representative of
#' each cluster.
#' @slot convergence a boolean indicating if none problem of empty class has been encountered.
#' @slot partial a boolean indicating the presence of partial rankings or ties.
#' @slot partialRank a matrix containing the full rankings, estimated using the within cluster ISR parameters
#' when the ranking is partial. When ranking is full, partialRank simply contains the
#' observed ranking. Available only in presence of at least one partial ranking.
#' @slot partialRankScore confidence score in estimated partial rank
#' @slot distanceProp Distances (MSE) between the final estimation and the current
#' value at each iteration of the SEM-Gibbs algorithm (except the burn-in phase) for proportions. A list of Qsem-Bsem elements,
#' each element being a K*p-matrix. 
#' @slot distancePi Distances (MSE) between the final estimation and the current
#' value at each iteration of the SEM-Gibbs algorithm (except the burn-in phase) for scale parameters. A list of Qsem-Bsem elements,
#' each element being a K*p-matrix.
#' @slot distanceMu Distances (Kendall distance) between the final estimation and the current
#' value at each iteration of the SEM-Gibbs algorithm (except the burn-in phase) for proportions. A list of Qsem-Bsem elements,
#' each element being a K*p-matrix.
#' @slot distanceZ a vector of size Qsem-Bsem containing the rand index between the final
#' estimated partition and the current value at each iteration of the SEM-Gibbs algorithm (except
#' the burn-in phase). Let precise that the rand index is not affected by label switching.
#' @slot distancePartialRank Kendall distance between the final estimation of the partial rankings
#' (missing positions in such rankings are estimated) and the current value at each iteration of the
#' SEM-Gibbs algorithm (except the burn-in phase). distancePartialRank is a list of Qsem-Bsem
#' elements, each element being a matrix of size n*p. Available only in presence of at least one
#' partial ranking.
#' @slot proportionInitial a vector containing the initialization of proportions in the algorithm.
#' @slot piInitial a matrix containing the initialization of the probabilities of good paired comparison in the algorithm.
#' @slot muInitial a matrix containing the initialization of modal rankings in the algorithm.
#' @slot partialRankInitial a matrix containing the initialization of the partial rankings in the algorithm.
#' 
#'
#'
#' @name Output-class
#' @rdname Output-class
#' @exportClass Output
#'
setClass(
  Class = "Output",
  representation = representation(
    proportion = "numeric",
    pi = "matrix",
    mu = "matrix",
    ll = "numeric",
    bic = "numeric",
    icl = "numeric",
    tik = "matrix",
    partition = "numeric",
    entropy = "matrix",
    probability = "matrix",
    convergence = "logical",
    partial = "logical",
    partialRank = "matrix",
    distanceProp = "list",
    distancePi = "list",
    distanceMu = "list",
    distanceZ = "numeric",
    distancePartialRank = "list",
    proportionInitial = "numeric",
    piInitial = "matrix",
    muInitial = "matrix",
    partialRankInitial = "matrix",
    partialRankScore = "matrix"
  ),
  prototype = prototype(
    proportion = numeric(0),
    pi = matrix(nrow = 0, ncol = 0),
    mu = matrix(nrow = 0, ncol = 0),
    ll = numeric(0),
    bic = numeric(0),
    icl = numeric(0),
    tik = matrix(nrow = 0, ncol = 0),
    partition = numeric(0),
    entropy = matrix(nrow = 0, ncol = 0),
    probability = matrix(nrow = 0, ncol = 0),
    convergence = logical(0),
    partial = logical(0),
    partialRank = matrix(nrow = 0, ncol = 0),
    distanceProp = list(),
    distancePi = list(),
    distanceMu = list(),
    distanceZ = numeric(0),
    distancePartialRank = list(),
    proportionInitial = numeric(0),
    piInitial = matrix(nrow = 0, ncol = 0),
    muInitial = matrix(nrow = 0, ncol = 0),
    partialRankInitial = matrix(nrow = 0, ncol = 0),
    partialRankScore = matrix(nrow = 0, ncol = 0)
  )
)



##################################################################################
#' @title Constructor of Rankclust class
#'
#' @description This class contains results of \link{rankclust} function.
#'
#' @slot K a vector of the number of clusters.
#' @slot data the data used for clustering.
#' @slot criterion the model selection criterion used.
#' @slot convergence a boolean indicating if none problem of empty class has been encountered (for
#' any number of clusters).
#' @slot results a list of \link{Output-class}, containing the results for each number of clusters (one
#' element of the list is associated to one number of clusters).
#' 
#'
#'
#' @details
#' If \code{res} is the result of \link{rankclust} function, each slot of results can be reached by \code{res[k]@@slotname}, where
#' \code{k} is the number of clusters and \code{slotname} is the name of the slot we want to reach (see \link{Output-class}).
#' For the slots, \code{ll}, \code{bic}, \code{icl}, \code{res["slotname"]} returns a vector of size \code{k} containing the values of the
#' slot for each number of clusters.
#'
#' @name Rankclust-class
#' @rdname Rankclust-class
#' @exportClass Rankclust
#'
setClass(
  Class = "Rankclust",
  representation = representation(
    K = "numeric",
    results = "list",
    data = "matrix",
    criterion = "character",
    convergence = "logical"
  ),
  prototype = prototype(
    results = list(),
    data = matrix(ncol = 0, nrow = 0),
    K = numeric(0),
    criterion = "bic",
    convergence = logical(0)
  )

)


#' @title Getter method for rankclust output
#'
#' @description Extract values of various
#' slots of the output from the function \code{\link{rankclust}}.
#'
#' @param x object from which to extract element(s) or in which to replace element(s).
#' @param i the number of cluster of the element we want to extract or "bic", "icl", "ll"
#' @param j,drop not used 
#' 
#' @export
setMethod(
  f = "[",
  signature = signature(x = "Rankclust"),
  definition = function(x, i, j, drop) {
    if (x@convergence)
    {
      if (is.numeric(i))
      {
        if (i %in% x@K)
        {
          return(x@results[[which(x@K == i)]])
        }
        else
          stop("Invalid number of cluster.")
      }
      else
      {
        if (i == "bic")
        {
          bic = rep(NA, length(x@K))
          for (iter in 1:length(x@K))
          {
            if (x@results[[iter]]@convergence)
              bic[iter] = x@results[[iter]]@bic
          }
          return(bic)
        }
        else
        {
          if (i == "icl")
          {
            icl = rep(NA, length(x@K))
            for (iter in 1:length(x@K))
            {
              if (x@results[[iter]]@convergence)
                icl[iter] = x@results[[iter]]@icl
            }
            return(icl)
          }
          else
          {
            if (i == "ll")
            {
              ll = rep(NA, length(x@K))
              for (iter in 1:length(x@K))
              {
                if (x@results[[iter]]@convergence)
                  ll[iter] = x@results[[iter]]@ll
              }
              return(ll)
            }
            else
            {
              stop("Invalid Name.")
            }
          }
        }
      }
    }
  }
)

#' @title Summary function.
#'
#' @description This function gives the summary of an output from \code{rankclust}.
#'
#' @param object output object from \code{\link{rankclust}}.
#' @param ... Not used.
#' @export
setMethod(
  f = "summary",
  signature = "Rankclust",
  definition = function(object, ...) {
    if (object@convergence)
    {
      if (object@criterion == "bic")
      {
        BIC = c()
        for (i in object@K)
        {
          BIC = c(BIC, object@results[[which(object@K == i)]]@bic)
        }
        index = which(BIC == min(BIC))
      }
      else
      {
        ICL = c()
        for (i in object@K)
        {
          ICL = c(ICL, object@results[[which(object@K == i)]]@icl)
        }
        index = which(ICL == min(ICL))

      }

      cat("******************************************************************\n")
      cat("NUMBER OF CLUSTERS: ", object@K[index], "\n")
      if (object@criterion == "bic")
        cat(object@criterion, "=", object[object@K[index]]@bic)
      else
        cat(object@criterion, "=", object[object@K[index]]@icl)
      cat("\nLoglikelihood =", object[object@K[index]]@ll)
      cat("\n\n************************ PARAMETERS ******************************\n")
      cat("Proportion:", object[object@K[index]]@proportion)
      cat("\nReference ranks mu:\n")
      print(object[object@K[index]]@mu)
      cat("\nProbabilities pi:\n")
      print(object[object@K[index]]@pi)
      cat("\n************************ CLUSTERING ******************************\n")
      cat("Ranks with the highest entropy for each cluster:\n")
      for (i in 1:object@K[index])
      {
        # classe=object[object@K[index]]@entropy[object[object@K[index]]@entropy[,2]==i,]
        classe = which(object[object@K[index]]@entropy[, 2] == i)
        if (length(classe) != 0)
        {
          classe = classe[order(object[object@K[index]]@entropy[classe, 1], decreasing = TRUE)][1:min(5, length(classe))]
          # if(object@algorithm=="SEM")
          print(cbind(object@data[classe, ], object[object@K[index]]@entropy[classe, ]))
          # else
          # 	print(cbind(object@data[classe,-ncol(object@data)],object[object@K[index]]@entropy[classe,]))
          # if(length(classe)==2)
          # {
          # 	best5=classe
          # 	print(cbind(object@data[classe,-ncol(object@data)],best5[2:3]))
          # }
          # else
          # {
          # 	best5=classe[order(classe[,1],decreasing=TRUE),][1:min(5,nrow(classe)),]
          # 	print(cbind(object@data[best5[,1],-ncol(object@data)],best5[,2:3]))
          # }
        }

      }
      # rm(best5)
      cat("Ranks with the highest probability for each cluster:\n")
      for (i in 1:object@K[index])
      {
        # classe=object[object@K[index]]@probability[object[object@K[index]]@probability[,2]==i,]
        classe = which(object[object@K[index]]@probability[, 2] == i)
        if (length(classe) != 0)
        {
          classe = classe[order(object[object@K[index]]@probability[classe, 1], decreasing = TRUE)][1:min(5, length(classe))]
          # if(object@algorithm=="SEM")
          print(cbind(object@data[classe, ], object[object@K[index]]@probability[classe, ]))
          # else
          # 	print(cbind(object@data[classe,-ncol(object@data)],object[object@K[index]]@probability[classe,]))
          # if(length(classe)==2)
          # {
          # 	best5=classe
          # 	print(cbind(object@data[best5[1],-ncol(object@data)],best5[2:3]))
          # }
          # else
          # {
          # 	best5=classe[order(classe[,2],decreasing=TRUE),][1:min(5,nrow(classe)),]
          # 	print(cbind(object@data[best5[,1],-ncol(object@data)],best5[,2:3]))
          # }

        }
      }
      rm(classe)
      # rm(best5)
      if (object[object@K[index]]@partial)
      {
        cat("\n************************ PARTIAL RANK ****************************\n")
        if (min(50, nrow(object[object@K[index]]@partialRank)) == 50)
          cat("\nOnly the first 50 are printed, total length:", nrow(object[object@K[index]]@partialRank), "\n")
        print(object[object@K[index]]@partialRank[1:min(50, nrow(object[object@K[index]]@partialRank)), ])
      }

      cat("\n******************************************************************\n")
    }
    else
      cat("No convergence. Please retry\n")

  }
)



#' @title Show function.
#'
#' @description This function shows the elements of a given object.
#'
#' @param object an object of class \link{Output-class} or \link{Rankclust-class}.
#'
#' @rdname show-methods
#' @aliases show,Output-method
setMethod(
  f = "show",
  signature = "Output",
  definition = function(object) {
    cat("ll=", object@ll)
    cat("\nbic =", object@bic)
    cat("\nicl =", object@icl)
    cat("\nproportion:", object@proportion)
    cat("\nmu:\n")
    print(object@mu)
    cat("\npi:\n")
    print(object@pi)
    cat("\npartition:\n")
    print(object@partition[1:min(50, length(object@partition))])
    if (min(50, length(object@partition)) == 50)
      cat("\nOnly the first 50 are printed, total length:", length(object@partition))
    cat("\ntik:\n")
    print(object@tik[1:min(50, nrow(object@tik)), ])
    if (min(50, nrow(object@tik)) == 50)
      cat("\nOnly the first 50 rows are printed, total rows:", nrow(object@tik))
  }
)



#' @rdname show-methods
#' @aliases show,Rankclust-method
setMethod(
  f = "show",
  signature = "Rankclust",
  definition = function(object) {
    if(object@convergence)
    {
      for (i in object@K)
      {
        cat("\n******************************************************************\n")
        cat("Number of clusters:", i)
        cat("\n******************************************************************\n")
        show(object@results[[which(object@K == i)]])
        cat("\n******************************************************************\n")
        
      }
    }else{
      cat("Algorithm did not converge.\n")
    }

  }
)
