




#' Creates a model
#' @param formula the model formula
#' @param cloneNet create a deep copy of the network within the model object
#' @param theta the model parameters.
#' @details
#' Creates a C++ Model object. In general this isn't needed by most users of the
#' package.
#' @examples
#' data(ukFaculty)
#' model <- createCppModel(ukFaculty ~ edges)
#' model$calculate()
#' model$statistics()
createCppModel <- function(formula,
                           cloneNet = TRUE,
                           theta = NULL) {
  modelClass <- "Model"
  form <- formula
  env <- environment(form)
  net <- as.BinaryNet(eval(form[[2]], envir = env))
  if (cloneNet)
    net <- net$clone()
  terms <- .prepModelTerms(formula)
  model <- .makeCppModelFromTerms(terms, net, theta, modelClass)
  model
}


#
# constructs a model from terms output by .prepModelTerms
#
.makeCppModelFromTerms <- function(terms,
                                   net,
                                   theta = NULL,
                                   modelClass = "Model") {
  net <- as.BinaryNet(net)
  
  clss <- class(net)
  networkEngine <- substring(clss, 6, nchar(clss) - 3)
  ModelType <-
    eval(parse(text = paste(
      "lolog::", networkEngine, modelClass, sep = ""
    )))
  
  model <- new(ModelType)
  model$setNetwork(net)
  
  stats <- rev(terms$stats)
  offsets <- rev(terms$offsets)
  
  if (length(stats) > 0)
    for (i in 1:length(stats)) {
      t <-
        try(model$addStatistic(names(stats)[i], stats[[i]]), silent = TRUE)
      if (inherits(t, "try-error")) {
        to <-
          try(model$addOffset(names(offsets)[i], offsets[[i]]), silent = TRUE)
        if (inherits(to, "try-error"))
          stop(t)
      }
    }
  if (length(offsets) > 0)
    for (i in 1:length(offsets))
      model$addOffset(names(offsets)[i], offsets[[i]])
  if (!is.null(theta))
    model$setThetas(theta)
  model$setVertexOrder(as.integer(rank(terms$vertexOrder, ties.method = "min")))
  model
  
}

