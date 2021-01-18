#' Simplified population fit of gastric emptying data
#'
#' Compute coefficients v0, tempt and kappa of a mixed model fit to a linexp
#' function with one grouping variable
#' @param d A data frame with columns
#' \itemize{
#'   \item \code{record} Record descriptor as grouping variable, e.g. patient ID
#'   \item \code{minute} Time after meal or start of recording.
#'   \item \code{vol} Volume of meal or stomach
#'  }
#' @param pnlsTol The value of pnlsTol at the initial iteration.
#' See \code{\link[nlme]{nlmeControl}}  When the model does not converge, \code{pnlsTol} is multiplied by 5 and the iteration repeated until convergence or \code{pnlsTol >= 0.5}. The effective value of  \code{pnlsTol} is returned in a separate list item. When it is known that a data set converges badly, it is recommended to set the initial \code{pnlsTol} to a higher value, but below 0.5, for faster convergence.
#'
#' @param model \code{linexp} (default) or \code{powexp}
#' @param variant For both models, there are 3 variants
#' \itemize{
#'   \item \code{variant = 1} The most generic version with independent estimates
#'   of all three parameters per record
#'   (\code{random = v0 + tempt + kappa ~ 1 | record}). The most likely to fail
#'   for degenerate cases. If this variant converges, use it.
#'   \item \code{variant = 2 } Diagonal random effects (\code{random = pdDiag(v0 + tempt + kappa) ~ 1; groups =  ~record }). Better convergence in critical cases. Note:  I never found out why I have to use the \code{groups} parameter instead of the \code{|}; see also p. 380 of Pinheiro/Bates.
#'   \item \code{variant = 3} Since parameters \code{kappa} and \code{beta} respectively are the most difficult to estimate, these are fixed in this variant (\code{random = v0 + tempt ~ 1}). This variant converges in all reasonable cases, but the estimates of \code{kappa} and \code{beta} cannot be use for secondary between-group analysis. If you are only interested in \code{t50}, you can use this safe version.
#'   }
#'
#' @return A list of class nlme_gastempt with elements
#' \code{coef, summary, plot, pnlsTol, message}
#' \itemize{
#'   \item \code{coef} is a data frame with columns:
#'     \itemize{
#'       \item \code{record} Record descriptor, e.g. patient ID
#'       \item \code{v0} Initial volume at t=0
#'       \item \code{tempt} Emptying time constant
#'       \item \code{kappa} Parameter \code{kappa} for
#'             \code{model = linexp}
#'       \item \code{beta} Parameter \code{beta} for \code{model = powexp}
#'       \item \code{t50} Half-time of emptying
#'       \item \code{slope_t50} Slope in t50; typically in units of ml/minute
#'  }
#'  On error, coef is NULL
#'   \item \code{nlme_result} Result of the nlme fit; can be used for addition
#'      processing, e.g. to plot residuals or via \code{summary} to extract AIC.
#'      On error, nlme_result is NULL.
#'   \item \code{plot} A ggplot graph of data and prediction. Plot of raw data is
#'      returned even when convergence was not achieved.
#'   \item \code{pnlsTol} Effective value of pnlsTo after convergence or failure.
#'   \item \code{message} String "Ok" on success, and the error message of
#'      \code{nlme} on failure.
#'  }
#'
#' @export
#' @import assertthat
#' @importFrom nlme nlme nlmeControl pdDiag
#' @importFrom 	stats coef median na.omit predict rnorm	rt uniroot
#' @importFrom utils head	tail
#' @importFrom methods new
#' @importFrom tibble rownames_to_column
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' suppressWarnings(RNGversion("3.5.0"))
#' set.seed(4711)
#' d = simulate_gastempt(n_record = 10, kappa_mean = 0.9, kappa_std = 0.3,
#'                       model = linexp)$data
#' fit_d = nlme_gastempt(d)
#'# fit_d$coef # direct access
#' coef(fit_d) # better use accessor function
#' coef(fit_d, signif = 3) # Can also set number of digits
#'# Avoid ugly ggplot shading (not really needed...)
#' library(ggplot2)
#' theme_set(theme_bw() + theme(panel.spacing = grid::unit(0,"lines")))
#'# fit_d$plot  # direct access is possible
#' plot(fit_d) # better use accessor function

nlme_gastempt = function(d, pnlsTol = 1.E-3, model = linexp, variant = 1){
  assert_that(all(c("record", "minute","vol") %in% names(d)))
  # Only linexp and powexp are supported
  assert_that(identical(model, linexp) || identical(model, powexp))
  variant = as.integer(variant)
  # alternative models
  linexp_models = c(nlme_linexp_1, nlme_linexp_2, nlme_linexp_3)
  powexp_models = c(nlme_powexp_1, nlme_powexp_2, nlme_powexp_3)
  assert_that(variant > 0 && variant <= length(linexp_models) )

  minute = vol = 0 # Keep NOTEs off
  # Make sure that records are factors
  d$record = as.factor(d$record)
  if (nlevels(d$record) < 3)
    stop(paste("At least 3 records are required, but there were only",
                nlevels(d$record) ))
  # Estimated initial values
  v0 = median(tail(sort(d$vol)))
  # raw data plot
  plot = ggplot(d, aes(x = minute, y = vol)) + geom_point() +
    facet_wrap(~ record) +
    expand_limits(x = 0, y = 0) # force zeroes to be visible
  success = FALSE
  # Repeat with increasing pnlsTol
  if (identical(model, linexp)) {
    title = paste0("Fitted linexp function, variant ", variant)
    nlme_func = linexp_models[[variant]]
    start = list(fixed = c(v0 = v0,
         tempt = mean(head(d$minute[d$vol < v0/2]))/2.5, kappa = 0.7))
  } else {
    title = paste0("Fitted powexp function, variant ", variant)
    nlme_func = powexp_models[[variant]]
    start = list(fixed = c(v0 = v0,
       tempt = mean(head(d$minute[d$vol < v0/2]))/1, beta = 1))
  }
  while (!success && pnlsTol < 0.5) {
    d_nlme = nlme_func(d, start, pnlsTol)
    success = !inherits(d_nlme, "try-error")
    #if (!success) print(d_nlme)
    if (!success) pnlsTol = pnlsTol * 5
  }

  # Return only plot for checking if not successful
  if (!success) {
    # Plot start values
    newdata = expand.grid(record = levels(d$record),
              minute = seq(0, max(d$minute), length.out = 50))
    pl = c(list(t = newdata$minute), start$fixed)
    newdata$vol = do.call(model, pl)
    start_r = lapply(start$fixed, signif, 3)
    # TODO: Enable subtitle when version of ggplot > 2.1.0
    subtitle = paste0("Start values are shown in red: ",
           paste0(names(start$fixed), " = ", start_r , collapse = ","))
    plot = plot + geom_line(data = newdata, col = "red") +
      ggtitle(paste("No convergence:", title, ", pnlsTol = ", pnlsTol ),
              subtitle = subtitle)
    ret = list(coef = NULL, summary = NULL, pnlsTol = pnlsTol, plot = plot,
               message = paste("pnlsTol = ", signif(pnlsTol, 2),"\n",
                               as.character(d_nlme)))
    class(ret) = "nlme_gastempt"
    return(ret)
  }
  cf =  coef(d_nlme) %>%
    tibble::rownames_to_column("record")  %>%
    t50() # Adds columns t50 and slope_t50
  # Plot prediction
  newdata = expand.grid(record = levels(d$record),
                        minute = seq(0, max(d$minute), length.out = 50))
  newdata$vol = predict(d_nlme, newdata)
  plot = plot + geom_line(data = newdata, col = "#006400") +
    ggtitle(paste0(title, ", pnlsTol = ", pnlsTol), subtitle = comment(d))
  ret = list(coef = cf, nlme_result = d_nlme, plot = plot,
             pnlsTol = pnlsTol, message = "Ok")
  class(ret) = "nlme_gastempt"
  ret
}

# Core functions - only random part varies
nlme_linexp_c = function(d, start, pnlsTol, random ){
    suppressWarnings(try(
      nlme(vol~linexp(minute, v0, tempt, kappa),
           data = d,
           fixed = v0 + tempt + kappa ~1,
           random = random,
           groups = ~record,
           start = start,
           # We reduce maxIter: normally, only 4 iterations are required
           # rarely about 11, and everything else does not work
           control = nlmeControl(msMaxIter = 20,
             pnlsTol = pnlsTol, maxIter = 15),
           na.action = na.omit),
      silent = TRUE))
}

nlme_powexp_c = function(d, start, pnlsTol, random){
  suppressWarnings(try(
    nlme(vol~powexp(minute, v0, tempt, beta),
         data = d,
         fixed = v0 + tempt + beta ~1,
         random = random,
         groups = ~record,
         start = start,
         control = nlmeControl(pnlsTol = pnlsTol),
         na.action = na.omit),
    silent = TRUE))
}


nlme_linexp_1 = function(d, start, pnlsTol){
  nlme_linexp_c(d, start, pnlsTol, v0 + tempt + kappa ~ 1)
}

nlme_powexp_1 = function(d, start, pnlsTol){
  nlme_powexp_c(d, start, pnlsTol, v0 + tempt + beta ~ 1)
}

# Version with pdDiag
nlme_linexp_2 = function(d, start, pnlsTol){
  nlme_linexp_c(d, start, pnlsTol, pdDiag(v0 + tempt + kappa ~ 1))
}

nlme_powexp_2 = function(d, start, pnlsTol){
  nlme_powexp_c(d, start, pnlsTol, pdDiag(v0 + tempt + beta ~ 1))
}

# Version with fixed kappa and beta
nlme_linexp_3 = function(d, start, pnlsTol){
  nlme_linexp_c(d, start, pnlsTol, pdDiag(v0 + tempt ~ 1))
}

nlme_powexp_3 = function(d, start, pnlsTol){
  nlme_powexp_c(d, start, pnlsTol, pdDiag(v0 + tempt ~ 1))
}


#' Extract coefficients from nlme_gastempt result
#'
#' @param object Result of a call to nlme_gastempt
#' @param ... other arguments
#'
#' @return a data frame with coefficients. See \code{\link{nlme_gastempt}} for an example.
#' @export
coef.nlme_gastempt = function(object, ...){
  Call = match.call(expand.dots = TRUE)
  sigdig = as.integer(Call[["signif"]])
  cf = object$coef
  if (!is.null(Call[["signif"]])) {
    cf[,-1] = lapply(cf[,-1], signif, sigdig)
  }
  cf
}

#' Plot data points and fit curve of an nlme_gastempt fit
#'
#' @param x Result of a call to nlme_gastempt
#' @param ... other arguments
#'
#' @return a ggplot object. Use \code{print()} if used non-interactively to show the curve
#' @method plot nlme_gastempt
#' @export
plot.nlme_gastempt = function(x, ...){
  x$plot
}
