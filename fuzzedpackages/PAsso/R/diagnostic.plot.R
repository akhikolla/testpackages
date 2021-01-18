#' @title Residual-based diagnostic plots
#'
#' @description A set of visualization tools for the diagnostic of the fitted model in
#' the partial association analysis. It can provides a plot matrix including Q-Q plots,
#' residual-vs-fitted plots, residual-vs-covariate plots of all the fitted models.
#' This function also support the direct diagnostic of the cumulative link regression model
#' in the class of \code{\link[ordinal]{clm}}, \code{\link[stats]{glm}}, \code{\link[rms]{lrm}},
#' \code{\link[rms]{orm}}, \code{\link[MASS]{polr}}. Currently, \code{\link[VGAM]{vglm}}
#' is not supported.
#'
#' @param object The object in the support classes (This function is mainly designed
#' for \code{PAsso}).
#'
#' @param output A character string specifying what type of output to plot. Default is
#' \code{"qq"} which produces a plot matrix with quantile-quantile plots of the residuals.
#' \code{"fitted"} produces a plot matrix between residuals and all corresponding fitted responses.
#' \code{"covariates"} produces a plot matrix between residuals and corresponding covariate.
#'
#' @param ... This function is based on a modified version of \code{"autoplot"} function in
#' \code{"sure"}. Additional optional arguments can be passed onto for drawing various plots.
#'
#' @return A \code{"ggplot"} object for supported models. For class "PAsso", it returns a plot in
#' \code{"gtable"} object that combines diagnostic plots of all responses.
#'
#' @export diagnostic.plot
#'
#' @examples
#' # Import data for partial association analysis
#' data("ANES2016")
#' ANES2016$PreVote.num <- as.factor(ANES2016$PreVote.num)
#'
#' PAsso_3v <- PAsso(responses = c("PreVote.num", "PID", "selfLR"),
#'                   adjustments = c("income.num", "age", "edu.year"),
#'                   data = ANES2016, uni.model = "probit",
#'                   method = c("kendall"),
#'                   resids.type = "surrogate", jitter = "latent")
#'
#' diag_p1 <- diagnostic.plot(object = PAsso_3v, output = "qq")
#' diag_p2 <- diagnostic.plot(object = PAsso_3v, output = "fitted")
#' diag_p3 <- diagnostic.plot(object = PAsso_3v, output = "covariate")
#'
#' # Simply diagnose a model
#' # Fit cumulative link models
#'
#' fit1 <- ordinal::clm(PreVote.num ~ income.num + age + edu.year, data = ANES2016, link = "logit")
#'
#' # diagnostic.plot
#' plot_qq_1 <- diagnostic.plot(object = fit1, output = "qq")
#' plot_fit_1 <- diagnostic.plot(object = fit1, output = "fitted")
#' plot_cov_1 <- diagnostic.plot(object = fit1, output = "covariate")
#'
diagnostic.plot <- function(object, ...) {
  UseMethod("diagnostic.plot")
}

#' @return A "ggplot" object based on the input residuals.
#'
#' @rdname diagnostic.plot
#' @export
diagnostic.plot.default <- function(
  object, ...
){
  warn_str <- paste("diagnostic.plot does not know how to handle object of class ",
                    class(object),
                    "and can only be used on classes PAsso, PAsso.test, resid, clm, glm, lrm, orm, polr.")
  warning(paste(strwrap(warn_str), collapse = "\n"))
}

#' @return A "ggplot" object based on the input residuals.
#'
#' @inheritParams diagnostic.plot
#'
#' @rdname diagnostic.plot
#' @export
diagnostic.plot.resid <- function(
  object,
  output = c("qq", "fitted", "covariate"),
  ...
) {
  autoplot.resid(object=object, ...)
}


#' @param object The object in the support classes (This function is mainly designed
#' for \code{PAsso}).
#'
#' @inheritParams autoplot
#'
#' @param model_id A number refers to the index of the model that needs to be diagnosed. If NULL, all
#' models will be diagnosed.
#' @param x_name A string refers to the covariate name that needs to be diagnosed. If NULL, all adjustments
#' will be diagnosed.
#' @param ... Additional optional arguments can be passed onto \code{\link[ggplot2]{ggplot}} for drawing
#' various plots.
#'
#' @rdname diagnostic.plot
#'
#' @method diagnostic.plot PAsso
#'
#' @return A plot in "gtable" object that combines diagnostic plots of all responses.
#'
#' @export
diagnostic.plot.PAsso <- function(
  object,
  output = c("qq", "fitted", "covariate"),
  model_id = NULL,
  x_name = NULL,
  ...
) {
  # object = PAsso_2; output = "covariate"

  # What type of output plot to produce
  output <- match.arg(output, several.ok = FALSE)
  rep_SRs <- object$rep_SRs
  resp_name <- attr(object, "responses")

  x_name <- if (is.null(x_name)) {
      attr(object, "adjustments")[1]
    }

  # do.call("<-", list(x_name, object$data[,x_name]))
  assign(x = x_name, object$data[,x_name]) # Save covariate in a vector as its name.

  n_resp <- length(resp_name)
  nCol <- floor(sqrt(n_resp))
  plot_list <- list()

  if (is.null(model_id)) { # If the diagnostic model is not specified, return a combined plot.
    # return a matrix-plot including diagnostics of models.
    if (output == "qq") {
      for (i in 1:n_resp) {
        plot_list[[i]] <-
          # autoplot(rep_SRs[,1,i], output = output,
          autoplot(object$fitted.models[[i]], output = output,
                   resp_name = resp_name[i], ...)
      }
      # Save the combined plot
      return(do.call("grid.arrange", c(plot_list, ncol=nCol)))

    } else if (output == "fitted") {

      for (i in 1:n_resp) {
        plot_list[[i]] <-
          autoplot(object$fitted.models[[i]], output = output, resp_name = resp_name[i],
                   alpha = 0.5, ...)
      }
      # Save the combined plot
      return(do.call("grid.arrange", c(plot_list, ncol=nCol)))

    } else {
      adjust_name <- attr(object, "adjustments")
      n_adjust <- length(adjust_name)
      t_lenght <- n_resp*n_adjust
      adjust_id <- rep(1:n_adjust, times=n_resp) # make index for covariate name in the for loop
      resp_id <- rep(1:n_resp, each=n_adjust) # make index for response name in the for loop

      for (i in 1:(t_lenght)) {
        plot_list[[i]] <-
          autoplot(object$fitted.models[[resp_id[i]]], output = "covariate",
                   x = object$data[,adjust_name[adjust_id[i]]],
                   xlab = adjust_name[adjust_id[i]],
                   # resp_name = resp_name[resp_id[i]], ...)
                   resp_name = resp_name[resp_id[i]])

        if (i %% n_adjust != 1) { # First plot of each response, draw ylab, otherwise, no ylab
          plot_list[[i]] <- plot_list[[i]] + ylab("")
        }
        if ((i-1) %/% n_adjust != (n_resp-1)) { # last row of plot(last response), draw xlab, otherwise, no xlab
          plot_list[[i]] <- plot_list[[i]] + xlab("")
        }
      }
      # Save the combined plot
      return(do.call("grid.arrange", c(plot_list, ncol=n_adjust)))
    }
  } else if ((model_id > n_resp) | (model_id <= 0)) {
    stop("'model_id' needs to be between ", 1, " and number of responses ", n_resp, "!")
  } else { # If the diagnostic model IS specified, return corresponding plot.
    return(autoplot(object$fitted.models[[model_id]], output = output,
             resp_name = resp_name[model_id], x = get(x_name), ...))
  }

}


#' @inheritParams autoplot
#'
#' @return A "ggplot" object based on the residuals generated from glm object.
#'
#' @rdname diagnostic.plot
#' @method diagnostic.plot glm
#' @export
diagnostic.plot.glm <- function(
  object,
  output = c("qq", "fitted", "covariate"),
  x = NULL,
  fit = NULL,
  distribution = qnorm,
  ncol = NULL,
  alpha = 1,
  xlab = NULL,
  color = "#444444",
  shape = 19,
  size = 2,
  qqpoint.color = "#444444",
  qqpoint.shape = 19,
  qqpoint.size = 2,
  qqline.color = "#888888",
  qqline.linetype = "dashed",
  qqline.size = 1,
  smooth = TRUE,
  smooth.color = "red",
  smooth.linetype = 1,
  smooth.size = 1,
  fill = NULL,
  resp_name = NULL,
  ...
) {
  autoplot.glm(object=object, output= output,
               x = x,
               fit = fit,
               distribution = distribution,
               ncol = ncol,
               alpha = alpha,
               xlab = xlab,
               color = color,
               shape = shape,
               size = size,
               qqpoint.color = qqpoint.color,
               qqpoint.shape = qqpoint.shape,
               qqpoint.size = qqpoint.size,
               qqline.color = qqline.color,
               qqline.linetype = qqline.linetype,
               qqline.size = qqline.size,
               smooth = smooth,
               smooth.color = smooth.color,
               smooth.linetype = smooth.linetype,
               smooth.size = smooth.size,
               fill = fill,
               resp_name = resp_name, ...)
}

#' @return A "ggplot" object based on the residuals generated from clm object.
#'
#' @inheritParams autoplot
#'
#' @rdname diagnostic.plot
#' @method diagnostic.plot clm
#' @export
diagnostic.plot.clm <- function(
  object,
  output = c("qq", "fitted", "covariate"),
  ...
) {
  autoplot.clm(object=object, output= output, ...)
}


#' @return A "ggplot" object based on the residuals generated from lrm object.
#'
#' @inheritParams autoplot
#'
#' @rdname diagnostic.plot
#' @method diagnostic.plot lrm
#' @export
diagnostic.plot.lrm <- function(
  object,
  output = c("qq", "fitted", "covariate"),
  ...
) {
  autoplot.lrm(object=object, output= output, ...)
}

#' @return A "ggplot" object based on the residuals generated from orm object.
#'
#' @inheritParams autoplot
#'
#' @rdname diagnostic.plot
#' @method diagnostic.plot orm
#' @export
diagnostic.plot.orm <- function(
  object,
  output = c("qq", "fitted", "covariate"),
  ...
) {
  autoplot.orm(object=object, output= output, ...)
}

#' @return A "ggplot" object based on the residuals generated from polr object.
#'
#' @inheritParams autoplot
#'
#' @rdname diagnostic.plot
#' @method diagnostic.plot polr
#' @export
diagnostic.plot.polr <- function(
  object,
  output = c("qq", "fitted", "covariate"),
  ...
) {
  autoplot.polr(object=object, output= output, ...)
}
