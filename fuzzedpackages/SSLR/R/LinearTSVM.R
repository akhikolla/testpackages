#' @title General Interface for LinearTSVM model
#' @description model from RSSL package
#' Implementation of the Linear Support Vector Classifier. Can be solved in the Dual formulation, which is equivalent to \code{\link{SVM}} or the Primal formulation.
#' @param C Cost variable
#' @param Cstar numeric; Cost parameter of the unlabeled objects
#' @param scale Whether a z-transform should be applied (default: TRUE)
#' @param eps Small value to ensure positive definiteness of the matrix in QP formulation
#' @param s numeric; parameter controlling the loss function of the unlabeled objects
#' @param init numeric; Initial classifier parameters to start the convex concave procedure
#' @inheritParams RSSL::BaseClassifier
#' @example demo/LinearTSVM.R
#' @importFrom RSSL LinearTSVM
#' @export
LinearTSVMSSLR <- function(C = 1, Cstar = 0.1, s = 0, x_center = FALSE,
                           scale = FALSE, eps = 1e-06,
                           verbose = FALSE, init = NULL) {

  train_function <- function(x, y) {

    load_RSSL()

    number_classes <- length(levels(y))

    #Check binary problem
    if (number_classes > 2) {
      stop("LinearTSVMSSLR is for binary problems")
    }

    list_values <- get_x_y_And_unlabeled(x, y)

    model <- RSSL::LinearTSVM(X = list_values$x, y = list_values$y, X_u = list_values$X_u,
                               C = C, Cstar = Cstar, s = s, x_center = x_center,
                              scale = scale, eps = eps,
                              verbose = verbose, init = init)

    result <- list(
      model = model
    )

    result$classes = levels(y)
    result$pred.params = c("class","raw")
    result$mode = "classification"
    class(result) <- "LinearTSVMSSLR"

    return(result)
  }

  args <- list(
    C = C, Cstar = Cstar, s = s, x_center = x_center,
    scale = scale, eps = eps,
    verbose = verbose, init = init
  )

  new_model_sslr(train_function, "LinearTSVMSSLR", args)

}

#' @title Predict LinearTSVMSSLR
#' @param object is the object
#' @param x is the dataset
#' @param ... This parameter is included for compatibility reasons.
#' @method predict LinearTSVMSSLR
#' @importFrom stats predict
#' @importFrom magrittr %>%
predict.LinearTSVMSSLR <- function(object, x, ...) {

  result <- object$model %>% predict(x)

  result

}



