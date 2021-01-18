#' @title General Interface for LaplacianSVM model
#' @description model from RSSL package
#' Manifold regularization applied to the support vector machine as proposed in Belkin et al. (2006). As an adjacency matrix, we use the k nearest neighbour graph based on a chosen distance (default: euclidean).
#' @param adjacency_k integer; Number of of neighbours used to construct adjacency graph.
#' @param adjacency_distance character; distance metric used to construct adjacency graph from the dist function. Default: "euclidean"
#' @param normalized_laplacian logical; If TRUE use the normalized Laplacian, otherwise, the Laplacian is used
#' @param gamma numeric; Weight of the unlabeled data
#' @param eps numeric; Small value to ensure positive definiteness of the matrix in the QP formulation
#' @inheritParams RSSL::BaseClassifier
#' @references Belkin, M., Niyogi, P. & Sindhwani, V., 2006. Manifold regularization:
#' A geometric framework for learning from labeled and unlabeled examples. Journal of
#' Machine Learning Research, 7, pp.2399-2434.
#' @example demo/LaplacianSVM.R
#' @importFrom RSSL LaplacianSVM
#' @export
LaplacianSVMSSLR <- function(lambda = 1, gamma = 1, scale = TRUE,
                              kernel = kernlab::vanilladot(), adjacency_distance = "euclidean",
                              adjacency_k = 6, normalized_laplacian = FALSE, eps = 1e-09) {

  train_function <- function(x, y) {

    load_RSSL()

    number_classes <- length(levels(y))

    #Check binary problem
    if (number_classes > 2) {
      stop("LaplacianSVMSSLR is for binary problems")
    }

    list_values <- get_x_y_And_unlabeled(x, y)

    model <- RSSL::LaplacianSVM(X = list_values$x, y = list_values$y, X_u = list_values$X_u,
                                 lambda = lambda, gamma = gamma, scale = scale,
                                 kernel = kernel, adjacency_distance = adjacency_distance,
                                 adjacency_k = adjacency_k, normalized_laplacian = normalized_laplacian,
                                 eps = eps)

    result <- list(
      model = model
    )

    result$classes = levels(y)
    result$pred.params = c("class","raw")
    result$mode = "classification"
    class(result) <- "LaplacianSVMSSLR"

    return(result)
  }

  args <- list(
    lambda = lambda, gamma = gamma, scale = scale,
    kernel = kernel, adjacency_distance = adjacency_distance,
    adjacency_k = adjacency_k, normalized_laplacian = normalized_laplacian,
    eps = eps
  )

  new_model_sslr(train_function, "LaplacianSVMSSLR", args)

}


#' @title Predict LaplacianSVMSSLR
#' @param object is the object
#' @param x is the dataset
#' @param ... This parameter is included for compatibility reasons.
#' @method predict LaplacianSVMSSLR
#' @importFrom magrittr %>%
predict.LaplacianSVMSSLR <- function(object, x, ...) {

  result <- object$model %>% predict(x)

  result
}



