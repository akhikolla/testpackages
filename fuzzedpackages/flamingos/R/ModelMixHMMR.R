#' A Reference Class which represents a fitted mixture of HMMR model.
#'
#' ModelMixHMMR represents an estimated mixture of HMMR model.
#'
#' @field param A [ParamMixHMMR][ParamMixHMMR] object. It contains the estimated
#'   values of the parameters.
#' @field stat A [StatMixHMMR][StatMixHMMR] object. It contains all the
#'   statistics associated to the MixHMMR model.
#' @seealso [ParamMixHMMR], [StatMixHMMR]
#' @export
#'
#' @examples
#' data(toydataset)
#' x <- toydataset$x
#' Y <- t(toydataset[,2:ncol(toydataset)])
#'
#' mixhmmr <- emMixHMMR(X = x, Y = Y, K = 3, R = 3, p = 1, verbose = TRUE)
#'
#' # mixhmmr is a ModelMixHMMR object. It contains some methods such as 'summary' and 'plot'
#' mixhmmr$summary()
#' mixhmmr$plot()
#'
#' # mixhmmr has also two fields, stat and param which are reference classes as well
#'
#' # Log-likelihood:
#' mixhmmr$stat$loglik
#'
#' # Parameters of the polynomial regressions:
#' mixhmmr$param$beta
ModelMixHMMR <- setRefClass(
  "ModelMixHMMR",
  fields = list(
    param = "ParamMixHMMR",
    stat = "StatMixHMMR"
  ),
  methods = list(
    plot = function(what = c("clustered", "smoothed", "loglikelihood"), ...) {
      "Plot method
      \\describe{
        \\item{\\code{what}}{The type of graph requested:
          \\itemize{
            \\item \\code{\"clustered\" = } Clustered curves (field
              \\code{klas} of class \\link{StatMixHMMR}).
            \\item \\code{\"smoothed\" = } Smoothed signal (field
              \\code{smoothed} of class {StatMixHMMR}).
            \\item \\code{\"loglikelihood\" = } Value of the log-likelihood for
              each iteration (field \\code{stored_loglik} of class
              \\link{StatMixHMMR}).
          }
        }
        \\item{\\code{\\dots}}{Other graphics parameters.}
      }"

      what <- match.arg(what, several.ok = TRUE)

      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar), add = TRUE)

      # yaxislim <- c(min(param$fData$Y) - 2 * mean(sqrt(apply(param$fData$Y, 1, var))), max(param$fData$Y) + 2 * mean(sqrt(apply(param$fData$Y, 1, var))))

      colorsvec <- rainbow(param$K)

      if (any(what == "clustered")) {
        par(mfrow = c(1, 1))
        matplot(param$fData$X, t(param$fData$Y), type = "l", lty = "dotted", col = colorsvec[stat$klas], xlab = "x", ylab = "y", ...)
        legend("bottomright", legend = sapply(1:param$K, function(x) paste0("Cluster ", x)), col = colorsvec, lty = "dotted", cex = 0.8)
        title(main = "Clustered curves")
      }

      if (any(what == "smoothed")) {

        nonemptyclusters = length(unique(stat$klas))
        par(mfrow = c(ceiling(sqrt(nonemptyclusters + 1)), round(sqrt(nonemptyclusters + 1))), mai = c(0.6, 0.6, 0.5, 0.25), mgp = c(2, 1, 0))

        matplot(param$fData$X, t(param$fData$Y), type = "l", lty = "solid", col = "black", xlab = "x", ylab = "y", ...)
        title(main = "Original dataset")

        for (k in 1:param$K) {
          if (sum(stat$klas == k) >= 1) {# At least one curve belongs to cluster k

            if (sum(stat$klas == k) == 1) {# Only one curve in cluster k
              matplot(param$fData$X, param$fData$Y[stat$klas == k,], type = "l", lty = "dotted", col = colorsvec[k], xlab = "x", ylab = "y", ...)
            } else {
              matplot(param$fData$X, t(param$fData$Y[stat$klas == k,]), type = "l", lty = "dotted", col = colorsvec[k], xlab = "x", ylab = "y", ...)
            }
            title(main = sprintf("Cluster %1.1i", k))
            lines(param$fData$X, stat$smoothed[, k], lwd = 1.5, ...)
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

      title <- paste("Fitted mixHMMR model")
      txt <- paste(rep("-", min(nchar(title) + 4, getOption("width"))), collapse = "")

      # Title
      cat(txt)
      cat("\n")
      cat(title)
      cat("\n")
      cat(txt)

      cat("\n")
      cat("\n")
      cat(
        paste0("MixHMMR model with K = ", param$K, ifelse(param$K > 1, " clusters", " cluster"), " and R = ", param$R, ifelse(param$R > 1, " regimes", " regime"), ":"))
      cat("\n")
      cat("\n")

      tab <- data.frame("log-likelihood" = stat$loglik, "nu" = param$nu,
                        "AIC" = stat$AIC, "BIC" = stat$BIC, "ICL" = stat$ICL1,
                        row.names = "", check.names = FALSE)
      print(tab, digits = digits)

      cat("\nClustering table (Number of curves in each clusters):\n")
      print(table(stat$klas))

      cat("\nMixing probabilities (cluster weights):\n")
      pro <- data.frame(t(param$alpha))
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

        if (param$p > 0) {
          row.names = c("1", sapply(1:param$p, function(x) paste0("X^", x)))
          betas <- data.frame(param$beta[, , k], row.names = row.names)
        } else {
          row.names = "1"
          betas <- data.frame(t(param$beta[, , k]), row.names = row.names)
        }
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
