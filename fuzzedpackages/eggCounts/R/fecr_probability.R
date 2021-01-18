fecr_probs <- function(stanFit, threshold = 0.95, lessthan = TRUE, plot = TRUE, xlab, ylab, main, verbose = TRUE, ...){
  # check if it is a valid model for computing probability
  modelName <- stanFit@model_name
  if (modelName %in% c("Bayesian model without zero-inflation","nb","Zero-inflated Bayesian model","zinb")) stop(
    "There is no reduction parameter for this model.")
  
  # extract fecr
  if (modelName == "indefficacy"){
    deltaMeansSample <- extract(stanFit,"delta_mu")[[1]]
    deltaShapeSample <- extract(stanFit,"delta_shape")[[1]]
    fecr <- 1 - qgamma(0.5, shape = deltaShapeSample, scale = deltaMeansSample/deltaShapeSample)
  } else {
    fecr <- 1 - extract(stanFit,"delta")[[1]]
  }

  prob95 <- mean(fecr<threshold)
  if (verbose) cat(paste0(c("The probability that the reduction is less than",threshold,"is",if(lessthan){prob95*100}else{100-prob95*100},"%.\n")))
  if (plot){
    d <- density(fecr, to = 1)
    if (missing(main)){main <- ""}
    if (missing(xlab)){xlab <- expression(1-delta)}
    if (missing(ylab)){ylab <- "Density"}
    plot(d, xlab = xlab, ylab = ylab, main = main, ...)
    if (lessthan){
      polygon(c(d$x[d$x<=threshold],rev(d$x[d$x<=threshold])),c(d$y[d$x<=threshold],rep(0,sum(d$x<=threshold))), col="gray")
    } else {
      polygon(c(d$x[d$x>=threshold],rev(d$x[d$x>=threshold])),c(d$y[d$x>=threshold],rep(0,sum(d$x>=threshold))), col="gray")
    }
    abline(v = threshold, col = 8)
  }
  return(invisible(if(lessthan){prob95*100}else{100-prob95*100}))
}