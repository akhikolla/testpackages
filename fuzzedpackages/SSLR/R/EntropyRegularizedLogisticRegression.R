#' @title General Interface for EntropyRegularizedLogisticRegression model
#' @description model from RSSL package
#' R Implementation of entropy regularized logistic regression implementation
#' as proposed by Grandvalet & Bengio (2005). An extra term is added to the objective
#' function of logistic regression that penalizes the entropy of the posterior measured
#' on the unlabeled examples.
#' @param lambda l2 Regularization
#' @param lambda_entropy Weight of the labeled observations compared to the unlabeled observations
#' @param init Initial parameters for the gradient descent
#' @inheritParams RSSL::BaseClassifier
#' @references Grandvalet, Y. & Bengio, Y., 2005. Semi-supervised learning by entropy
#' minimization. In L. K. Saul, Y. Weiss, & L. Bottou, eds. Advances in Neural Information
#' Processing Systems 17. Cambridge, MA: MIT Press, pp. 529-536.
#' @example demo/EntropyRegularizedLogisticRegression.R
#' @importFrom RSSL EntropyRegularizedLogisticRegression
#' @export
EntropyRegularizedLogisticRegressionSSLR <- function(lambda = 0,
                                                       lambda_entropy = 1, intercept = TRUE,
                                                       init = NA, scale = FALSE,
                                                       x_center = FALSE) {

  train_function <- function(x, y) {

    load_RSSL()

    list_values <- get_x_y_And_unlabeled(x, y)

    model <- RSSL::EntropyRegularizedLogisticRegression(X = list_values$x, y = list_values$y, X_u = list_values$X_u,
                                                         lambda = lambda,
                                                         lambda_entropy = lambda_entropy, intercept = intercept,
                                                         init = init, scale = scale,
                                                         x_center = x_center)

    result <- list(
      model = model
    )

    result$classes = levels(y)
    result$pred.params = c("class","raw")
    result$mode = "classification"
    class(result) <- "EntropyRegularizedLogisticRegressionSSLR"

    return(result)
  }

  args <- list(
    lambda = lambda,
    lambda_entropy = lambda_entropy, intercept = intercept,
    init = init, scale = scale,
    x_center = x_center
  )

  new_model_sslr(train_function, "EntropyRegularizedLogisticRegressionSSLR", args)

}


#' @title Predict EntropyRegularizedLogisticRegressionSSLR
#' @param object is the object
#' @param x is the dataset
#' @param ... This parameter is included for compatibility reasons.
#' @method predict EntropyRegularizedLogisticRegressionSSLR
#' @importFrom stats predict
#' @importFrom magrittr %>%
predict.EntropyRegularizedLogisticRegressionSSLR <- function(object, x, ...) {

  result <- object$model %>% predict(x)

  result
}



