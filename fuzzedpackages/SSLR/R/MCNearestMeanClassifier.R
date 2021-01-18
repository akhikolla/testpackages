#' @title General Interface for MCNearestMeanClassifier (Moment Constrained Semi-supervised Nearest Mean Classifier) model
#' @description model from RSSL package
#' Update the means based on the moment constraints as defined in Loog (2010).
#' The means estimated using the labeled data are updated by making sure their
#' weighted mean corresponds to the overall mean on all (labeled and unlabeled) data.
#' Optionally, the estimated variance of the classes can be re-estimated after this
#' update is applied by setting update_sigma to \code{TRUE}. To get the true nearest mean
#' classifier, rather than estimate the class priors, set them to equal priors using, for
#' instance \code{prior=matrix(0.5,2)}.
#' @param update_sigma logical; Whether the estimate of the variance should be updated
#' after the means have been updated using the unlabeled data
#' @param prior matrix; Class priors for the classes
#' @inheritParams RSSL::BaseClassifier
#' @references Loog, M., 2010. Constrained Parameter Estimation for
#' Semi-Supervised Learning: The Case of the Nearest Mean Classifier.
#' In Proceedings of the 2010 European Conference on Machine
#' learning and Knowledge Discovery in Databases. pp. 291-304.
#' @example demo/MCNearestMeanClassifier.R
#' @importFrom RSSL MCNearestMeanClassifier
#' @export
MCNearestMeanClassifierSSLR <- function(update_sigma = FALSE, prior = NULL,
                                         x_center = FALSE, scale = FALSE) {

  train_function <- function(x, y) {

    load_RSSL()

    number_classes <- length(levels(y))

    #Check binary problem
    if (number_classes > 2) {
      stop("MCNearestMeanClassifierSSLR is for binary problems")
    }

    list_values <- get_x_y_And_unlabeled(x, y)

    model <- RSSL::MCNearestMeanClassifier(X = list_values$x, y = list_values$y, X_u = list_values$X_u,
                                            update_sigma = update_sigma, prior = prior,
                                            x_center = x_center, scale = scale)

    result <- list(
      model = model
    )

    result$classes = levels(y)
    result$pred.params = c("class","raw")
    result$mode = "classification"
    class(result) <- "MCNearestMeanClassifierSSLR"

    return(result)
  }

  args <- list(
    update_sigma = update_sigma, prior = prior,
    x_center = x_center, scale = scale
  )

  new_model_sslr(train_function, "MCNearestMeanClassifierSSLR", args)

}

#' @title Predict MCNearestMeanClassifierSSLR
#' @param object is the object
#' @param x is the dataset
#' @param ... This parameter is included for compatibility reasons.
#' @method predict MCNearestMeanClassifierSSLR
#' @importFrom magrittr %>%
predict.MCNearestMeanClassifierSSLR <- function(object, x, ...) {

  result <- object$model %>% predict(x)

  result

}



