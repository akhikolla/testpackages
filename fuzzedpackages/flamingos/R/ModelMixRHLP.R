#' A Reference Class which represents a fitted mixture of RHLP model.
#'
#' ModelMixRHLP represents an estimated mixture of RHLP model.
#'
#' @field param A [ParamMixRHLP][ParamMixRHLP] object. It contains the estimated
#'   values of the parameters.
#' @field stat A [StatMixRHLP][StatMixRHLP] object. It contains all the
#'   statistics associated to the MixRHLP model.
#' @seealso [ParamMixRHLP], [StatMixRHLP]
#' @export
#'
#' @examples
#' data(toydataset)
#'
#' # Let's fit a mixRHLP model on a dataset containing 2 clusters:
#' data <- toydataset[1:190,1:21]
#' x <- data$x
#' Y <- t(data[,2:ncol(data)])
#'
#' mixrhlp <- cemMixRHLP(X = x, Y = Y, K = 2, R = 2, p = 1, verbose = TRUE)
#'
#' # mixrhlp is a ModelMixRHLP object. It contains some methods such as 'summary' and 'plot'
#' mixrhlp$summary()
#' mixrhlp$plot()
#'
#' # mixrhlp has also two fields, stat and param which are reference classes as well
#'
#' # Log-likelihood:
#' mixrhlp$stat$loglik
#'
#' # Parameters of the polynomial regressions:
#' mixrhlp$param$beta
ModelMixRHLP <- setRefClass(
  "ModelMixRHLP",
  fields = list(
    param = "ParamMixRHLP",
    stat = "StatMixRHLP"
  ),
  methods = list(
    plot = function(what = c("estimatedsignal", "regressors", "loglikelihood"), ...) {
      "Plot method.
      \\describe{
        \\item{\\code{what}}{The type of graph requested:
          \\itemize{
            \\item \\code{\"estimatedsignal\" = } Estimated signal (field
              \\code{Ey} of class \\link{StatMixRHLP}).
            \\item \\code{\"regressors\" = } Polynomial regression components
              (fields \\code{polynomials} and \\code{pi_jkr} of class
              \\link{StatMixRHLP}).
            \\item \\code{\"loglikelihood\" = } Value of the log-likelihood for
              each iteration (field \\code{stored_loglik} of class
              \\link{StatMixRHLP}).
          }
        }
        \\item{\\code{\\dots}}{Other graphics parameters.}
      }
      By default, all the above graphs are produced."

      what <- match.arg(what, several.ok = TRUE)

      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar), add = TRUE)

      # yaxislim <- c(min(modelMixRHLP$Y) - 2 * mean(sqrt(apply(modelMixRHLP$Y, 1, var))), max(modelMixRHLP$Y) + 2 * mean(sqrt(apply(modelMixRHLP$Y, 1, var))))

      colorsvector = rainbow(param$K)

      if (any(what == "estimatedsignal")) {
        # Cluster and means
        nonemptyclusters = length(unique(stat$klas))
        par(mfrow = c(ceiling(sqrt(nonemptyclusters + 1)), round(sqrt(nonemptyclusters + 1))), mai = c(0.6, 0.6, 0.5, 0.25), mgp = c(2, 1, 0))

        matplot(param$fData$X, t(param$fData$Y), type = "l", lty = "solid", col = "black", xlab = "x", ylab = "y", ...)
        title(main = "Dataset")

        for (k in 1:param$K) {
          cluster_k = param$fData$Y[stat$klas == k, , drop = FALSE]

          if (length(cluster_k) != 0) {
            matplot(param$fData$X, t(cluster_k), type = "l", lty = "dotted", col = colorsvector[k], xlab = "x", ylab = "y", ...)
            lines(param$fData$X, stat$Ey[, k, drop = FALSE], col = "black", lty = "solid", lwd = 1.5, ...)
            title(main = sprintf("Cluster %1.1i", k))
          }
        }
      }

      if (any(what == "regressors")) {
        par(mfrow = c(2, 1), mai = c(0.6, 0.8, 0.5, 0.5))
        for (k in 1:param$K) {
          cluster_k = param$fData$Y[stat$klas == k, , drop = FALSE]

          if (length(cluster_k) != 0) {
            matplot(param$fData$X, t(cluster_k), type = "l", lty = "dotted", col = colorsvector[k], xlab = "x", ylab = "y", ...)

            # Polynomial regressors
            for (r in 1:param$R) {
              lines(param$fData$X, stat$polynomials[, r, k], col = "black", lty = "dotted", lwd = 1.5, ...)
            }

            lines(param$fData$X, stat$Ey[, k, drop = FALSE], col = "black", lty = "solid", lwd = 1.5, ...)
            title(main = sprintf("Cluster %1.1i", k))

            matplot(param$fData$X, stat$pi_jkr[1:param$fData$m, , k], type = "l", lty = "solid", xlab = "x", ylab = "Logistic proportions", ylim = c(0, 1), ...)
          }
        }
      }

      if (any(what == "loglikelihood")) {
        par(mfrow = c(1, 1))
        plot.default(1:length(stat$stored_loglik), stat$stored_loglik, type = "l", col = "blue", xlab = "EM iteration number", ylab = "Log-likelihood", ...)
        title(main = "Log-likelihood")
      }

    },

    summary = function(digits = getOption("digits")) {
      "Summary method.
      \\describe{
        \\item{\\code{digits}}{The number of significant digits to use when
          printing.}
      }"

      title <- paste("Fitted mixRHLP model")
      txt <- paste(rep("-", min(nchar(title) + 4, getOption("width"))), collapse = "")

      # Title
      cat(txt)
      cat("\n")
      cat(title)
      cat("\n")
      cat(txt)

      cat("\n")
      cat("\n")
      cat(paste0("MixRHLP model with K = ", param$K, ifelse(param$K > 1, " clusters", " cluster"), " and R = ", param$R, ifelse(param$R > 1, " regimes", " regime"), ":"))
      cat("\n")
      cat("\n")

      tab <- data.frame("log-likelihood" = stat$loglik, "nu" = param$nu,
                        "AIC" = stat$AIC, "BIC" = stat$BIC, "ICL" = stat$ICL,
                        row.names = "", check.names = FALSE)
      print(tab, digits = digits)

      cat("\nClustering table (Number of curves in each clusters):\n")
      print(table(stat$klas))

      cat("\nMixing probabilities (cluster weights):\n")
      pro <- data.frame(param$alpha)
      colnames(pro) <- 1:param$K
      print(pro, digits = digits, row.names = FALSE)

      cat("\n\n")

      txt <- paste(rep("-", min(nchar(title), getOption("width"))), collapse = "")

      for (k in 1:param$K) {
        cat(txt)
        cat("\nCluster ", k, " (k = ", k, "):\n", sep = "")

        cat("\nRegression coefficients for each regime/segment r (r=1...R):\n\n")
        if (param$p > 0) {
          row.names = c("1", sapply(1:param$p, function(x) paste0("X^", x)))
        } else {
          row.names = "1"
        }

        betas <- data.frame(matrix(param$beta[, , k], nrow = param$p + 1), row.names = row.names)
        colnames(betas) <- sapply(1:param$R, function(x) paste0("Beta(r = ", x, ")"))
        print(betas, digits = digits)

        if (param$variance_type == "homoskedastic") {
          sigma2 <- data.frame(param$sigma2[k])
          colnames(sigma2) <- "Sigma2"
          cat(paste0("\nVariance:\n\n"))
          print(sigma2, digits = digits, row.names = FALSE)
        } else {
          sigma2 <- data.frame(t(param$sigma2[, k]))
          colnames(sigma2) <- sapply(1:param$R, function(x) paste0("Sigma2(r = ", x, ")"))
          cat(paste0("\nVariances:\n\n"))
          print(sigma2, digits = digits, row.names = FALSE)
        }
        cat("\n")
      }

    }

  )
)
