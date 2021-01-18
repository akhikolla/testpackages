#' @title Methods for flexreg Objects
#'
#' @description Methods for extracting information from fitted  regression model objects of class \code{`flexreg`}.
#'
#' @param object an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}}.
#' @param digits an integer indicating the number of decimal places. Default equal to 4.
#' @param ... additional arguments.
#'
#' @details  The \code{summary.flexreg} method summarizes the results of \code{\link{flexreg}}, adding also information from the functions
#' \code{\link{residuals.flexreg}} and \code{\link{WAIC}}. The \code{summary.flexreg} method returns an object of class \code{`summary.flexreg`} containing the relevant summary statistics which can subsequently be
#' printed using the associated \code{print} method.
#'
#' @examples
#' data("Reading")
#' dataset <- Reading
#' FB <- flexreg(accuracy ~ iq, dataset, n.iter = 1000)
#' summary(FB)
#'
#' @import rstan
#' @method summary flexreg
#' @export
#'

summary.flexreg <- function(object, ..., digits=4){
  x <- object
  call <- x$call
  model.name <- x$model@model_name
  link.mu <- x$link.mu
  link.phi <- x$link.phi
  formula <- x$formula

  covariate.names.mu <- colnames(x$design.X)
  covariate.names.phi <- colnames(x$design.Z)

  posterior <- x$model
  pars <- extract.pars(posterior = posterior)
  n.pars <- length(pars)
  pp <- rstan::extract(posterior, pars)
  summa <- lapply(pp, function(x) c(mean(x), sd(x),quantile(x, probs = c(.025,.5,.975))))
  summa.mat <- round(matrix(unlist(summa), ncol=5, byrow = T ),digits)
  colnames(summa.mat) <- c("Post. Mean", "Post. SD", "2.5%", "Post. Median",  "97.5%")
  rownames(summa.mat) <- pars

  summ.mu <- summa.mat[grep("beta",rownames(summa.mat)),]
  rownames(summ.mu) <- covariate.names.mu

  if(is.null(covariate.names.phi)){
    summ.phi <- summa.mat[which(rownames(summa.mat)=="phi"),]
    dim(summ.phi) <- c(1,5)
    colnames(summ.phi) <- c("Post. Mean", "Post. SD", "2.5%", "Post. Median",  "97.5%")
    rownames(summ.phi) <- "phi"
    } else {
    summ.phi <- summa.mat[grep("psi",rownames(summa.mat)),]
    rownames(summ.phi) <- covariate.names.phi
    }

  n.parz <- nrow(summ.mu)+nrow(summ.phi)
  if( n.pars > n.parz) {
    summ.add <- summa.mat[(n.parz+1):n.pars,]
    } else summ.add <- NULL

  residuals <- residuals.flexreg(x, type = "raw", cluster=FALSE, estimate="mean")
  summ.res <- round(quantile(residuals), digits)
  names(summ.res) <- c("Min", "1Q", "Median", "3Q", "Max")

  waic_out <- suppressWarnings(WAIC(x))

output <- list(call=call, Model=model.name, formula=formula, link.mu=link.mu, link.phi=link.phi,
                 Summary.res=summ.res,
                 Summary.mu=summ.mu, Summary.phi=summ.phi, Summary.add=summ.add,
                 waic_out = waic_out)
  class(output) <- "summary.flexreg"
  return(output)
}


#' Print Methods for summary.flexreg Objects
#'
#' @param x an object of class \code{`flexreg`}, usually the result of \code{\link{flexreg}}.

#' @rdname summary.flexreg
#' @export
#'

print.summary.flexreg <- function(x, ...){
  cat("Call: ")
  print(x$call)
  cat("\nModel name: ", x$Model, "\n \n")

  cat("Residuals:\n")
  print(x$Summary.res)

  cat("\nCoefficients (mean model with", x$link.mu, "link):\n")
  print(x$Summary.mu)

  cat("\nCoefficients (precision model with", x$link.phi, "link):\n")
  print(x$Summary.phi)

  if(!is.null(x$Summary.add)){
    cat("\nAdditional Parameters:\n")
    print(x$Summary.add)
  }

  cat("\nWaic method:")
  suppressWarnings(print(x$waic_out$waic_out))

}


#' @title Plot method for flexreg Objects
#'
#' @description Method for plotting regression curves for the mean from fitted regression model objects of class \code{`flexreg`}.
#'
#' @param x an object of class flexreg, usually the result of \code{\link{flexreg}}.
#' @param name.x a character containing the name of the covariate to be plotted on the x-axis of the scatterplot.
#' @param additional.cov.default a list of additional covariates to be set as default.
#' @param ... additional arguments. Currently not used.
#'
#' @details The function produces a scatterplot of the covariate specified in \code{name.x} and the response \code{y}. Any other variable involved in the formula must be set to a default through the \code{additional.cov.default} argument.
#' If the regression model is of \code{FB} type the function returns a scatterplot with three curves, one corresponding to the overall mean and two corresponding to the component means of the FB distribution, i.e., \eqn{\lambda_1} and \eqn{\lambda_1}.
#'
#'
#' @examples
#' data("Reading")
#' FB <- flexreg(accuracy ~ iq+dyslexia, Reading, n.iter=800)
#' plot(FB, name.x="iq", additional.cov.default = list("dyslexia"=1))
#'
#'
#' @import ggplot2
#'
#' @method plot flexreg
#' @export
#'

plot.flexreg <- function(x, name.x, additional.cov.default = NA, ...)
  {
  group <- Response <- NULL
  #assign("group", "Response", envir = .GlobalEnv)
  object <- x
  y <- object$response
  x <- object$design.X[,which(colnames(object$design.X) == name.x)]

 # additional.cov.default <- list("x2" = 0)#list("dyslexia"= -1, "iq:dyslexia" = 0)

  #prova con modello senza intercetta
  intercept <- ifelse("(Intercept)"  %in% colnames(object$design.X), T, F)
  if( intercept == T){
  newdata <- data.frame(1, x, additional.cov.default)
  names(newdata) <- c("(Intercept)", name.x, names(additional.cov.default))
  } else {
  newdata <- data.frame(x, additional.cov.default)
  names(newdata) <- c(name.x, names(additional.cov.default))
  }
  cluster <- "FB" %in% object$call$type | is.null(object$call$type)
  if (cluster == T) {
  mu.hat <- predict(object, newdata = newdata, type = "response", cluster = T)
  data.plot <- data.frame(y = y, x = x, mu.hat)
  data.plot <- reshape(data.plot, direction = "long",
                       varying = 3:5, v.names=c('Response'))
  data.plot$group <- as.factor(data.plot$time)
  plot1 <- ggplot(data.plot, aes(x=x, y=y, group = group))+ geom_point()+
    geom_line(aes(x=x, y=Response, colour = group))
  } else {
    mu.hat <- predict(object, newdata = newdata, type = "response", cluster = F)
    data.plot <- data.frame(y = y, x = x,  mu.hat)
    names(data.plot)[3] <- "Response"
    plot1 <- ggplot(data.plot, aes(x=x, y=y))+
      geom_line(aes(x=x, y=Response))
  }
  return(plot1+ geom_point()+theme_bw()+
           scale_color_manual(labels = c(expression(mu), expression(lambda[1]),expression(lambda[2])), values = c("black", "red", "blue"))+
        theme(legend.title = element_blank()))
}





