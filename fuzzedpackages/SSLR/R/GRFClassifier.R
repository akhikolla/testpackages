#' @title General Interface for GRFClassifier (Transductive SVM classifier using the convex concave procedure) model
#' @description model from RSSL package
#' Implements the approach proposed in Zhu et al. (2003) to label propagation over
#' an affinity graph. Note, as in the original paper, we consider the transductive
#' scenario, so the implementation does not generalize to out of sample predictions.
#' The approach minimizes the squared difference in labels assigned to different objects,
#' where the contribution of each difference to the loss is weighted by the affinity between
#' the objects. The default in this implementation is to use a knn adjacency matrix based on euclidean
#' distance to determine this weight. Setting adjacency="heat" will use an RBF kernel over
#' euclidean distances between objects to determine the weights.
#' @param adjacency character; "nn" for nearest neighbour graph or "heat" for radial basis adjacency matrix
#' @param adjacency_distance character; distance metric for nearest neighbour adjacency matrix
#' @param adjacency_k integer; number of neighbours for the nearest neighbour adjacency matrix
#' @param scale logical; Should the features be normalized? (default: FALSE)
#' @param x_center logical; Should the features be centered?
#' @param adjacency_sigma double; width of the rbf adjacency matrix
#' @param class_mass_normalization logical; Should the Class Mass Normalization heuristic be applied? (default: TRUE)
#' @references Collobert, R. et al., 2006. Large scale transductive SVMs.
#' Journal of Machine Learning Research, 7, pp.1687-1712.
#' @example demo/GRFClassifier.R
#' @importFrom RSSL GRFClassifier
#' @importFrom RSSL responsibilities
#' @export
GRFClassifierSSLR <- function(adjacency = "nn",
                     adjacency_distance = "euclidean", adjacency_k = 6,
                     adjacency_sigma = 0.1, class_mass_normalization = TRUE,
                     scale = FALSE, x_center = FALSE) {

  train_function <- function(x, y) {

    load_RSSL()

    number_classes <- length(levels(y))

    #Check binary problem
    'if (number_classes > 2) {
      stop("TSVMSSLR is for binary problems")
    }'

    list_values <- get_x_y_And_unlabeled(x, y)

    model <- RSSL::GRFClassifier(X = list_values$x, y = list_values$y, X_u = list_values$X_u,
                                 adjacency = adjacency, adjacency_distance = adjacency_distance,
                                 adjacency_k = adjacency_k, adjacency_sigma = adjacency_sigma,
                                 class_mass_normalization = class_mass_normalization,
                                 x_center = x_center, scale = scale)

    result <- list(
      model = model
    )

    assignment <- factor(apply(responsibilities(model),1,which.max))
    result$classes = levels(y)
    #result$pred.params = c("class","raw")
    result$mode = "classification"
    result$labels_unlabeled = assignment
    class(result) <- "GRFClassifierSSLR"

    return(result)
  }

  args <- list(
    adjacency = adjacency, adjacency_distance = adjacency_distance,
    adjacency_k = adjacency_k, adjacency_sigma = adjacency_sigma,
    class_mass_normalization = class_mass_normalization,
    x_center = x_center, scale = scale
  )

  new_model_sslr(train_function, "GRFClassifierSSLR", args)

}

#' Predictions
#' @title predictions unlabeled data
#' @param object object
#' @param ... other parameters to be passed
#' @export
predictions <- function(object, ...){
  UseMethod("predictions")
}

#' Predictions
#' @title predictions unlabeled data
#' @param object object
#' @param ... other parameters to be passed
#' @export predictions.GRFClassifierSSLR
#' @export
predictions.GRFClassifierSSLR <- function(object,...) {

  result <- object$labels_unlabeled
  result

}



