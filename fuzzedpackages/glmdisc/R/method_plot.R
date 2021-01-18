#' @include allClasses.R
NULL

#' Plots for the discretized / grouped data.
#'
#'
#' @exportMethod plot
#' @name plot
#' @aliases plot
#' @docType methods
#' @concept test discretization plot
#' @importFrom gam s
#'
if (!isGeneric("plot")) {
  methods::setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
}

#' Plots for the discretized data.
#'
#'
#' @rdname plot
#' @aliases plot,glmdisc,ANY-method
#' @param x The S4 \code{\link{glmdisc}} object to plot.
#' @examples
#' # Simulation of a discretized logit model
#' set.seed(1)
#' x <- matrix(runif(300), nrow = 100, ncol = 3)
#' cuts <- seq(0, 1, length.out = 4)
#' xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
#' theta <- t(matrix(c(0, 0, 0, 2, 2, 2, -2, -2, -2), ncol = 3, nrow = 3))
#' log_odd <- rowSums(t(sapply(seq_along(xd[, 1]), function(row_id) {
#'   sapply(
#'     seq_along(xd[row_id, ]),
#'     function(element) theta[xd[row_id, element], element]
#'   )
#' })))
#' y <- rbinom(100, 1, 1 / (1 + exp(-log_odd)))
#'
#' sem_disc <- glmdisc(x, y,
#'   iter = 50, m_start = 4, test = FALSE,
#'   validation = FALSE, criterion = "aic"
#' )
#' plot(sem_disc)
plot.glmdisc <- function(x) {

  # Graph 1 : ROC CURVE
  glm_simple_roc <- simple_roc(x@disc.data$labels, predict(x, x@cont.data[, -ncol(x@cont.data), drop = FALSE]))

  with(glm_simple_roc, graphics::plot(1 - FPR, TPR, col = 1 + labels, xlim = c(1, 0)))
  graphics::title(main = "ROC Curve on total provided data")
  graphics::lines(c(1, 0), c(0, 1))
  graphics::legend(c(0.4, 0.4), c(paste("Gini:", round(normalizedGini(x@cont.data$labels, predict(x, x@cont.data[, -ncol(x@cont.data), drop = FALSE])), 2))))

  # installed_graphical_package(object@best.disc$bestLogisticRegression)
  graphics::par(ask = TRUE)

  # Group of Graphs 2 : risk levels for each discretized attribute
  for (j in 1:(ncol(x@disc.data) - 1)) {
    graphics::spineplot(factor(x@disc.data[, j]), factor(x@disc.data$labels), main = paste("Risk levels of attribute", colnames(x@disc.data)[j]), xlab = "Discretized / grouped feature", ylab = "Class proportions")
  }

  # Group of Graphs 3 : continuous features + glmnet
  # par(ask=FALSE)
  for (j in which(x@parameters$types_data == "numeric")) {
    modele_gam <- gam::gam(stats::as.formula(paste("labels ~", paste("s(", gsub("+", ",5) + s(", as.character(x@best.disc$formulaOfBestLogisticRegression)[2], fixed = TRUE), ",5)"))), data = x@cont.data, family = "binomial")
    gam::plot.Gam(modele_gam, ask = FALSE, terms = paste0("s(", colnames(x@disc.data)[j], ", 5)"), se = TRUE, main = paste("Plot of a GAM fit to attribute", colnames(x@disc.data)[j]))
    graphics::par(ask = TRUE)
  }
  # par(ask=TRUE)

  # Group of Graphs 4 : categorical features significance
  for (j in which(x@parameters$types_data == "factor")) {
    print(table(factor(x@disc.data[, j]), factor(x@cont.data[, j])))
  }

  on.exit(graphics::par(ask = FALSE))
}

#' @rdname plot
#' @name plot,glmdisc,missing-method
#' @aliases plot,glmdisc,missing-method
#' @description This defines the \code{\link{plot}} method which will plot some useful graphs for the discretization scheme of S4 class \code{\link{glmdisc}}
# #' @param glmdisc The S4 glmdisc object to plot.
methods::setMethod(f = "plot", signature = c(x = "glmdisc", y = "missing"), definition = plot.glmdisc)
