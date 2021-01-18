#' Runs t-tests if model parameters differ from user defined values
#' @useDynLib stratEst,.registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @param model a fitted model of class \code{stratEst.model}.
#' @param par a character vector. The class of model parameters to be tested. Default is \code{c("shares","probs","trembles", "coefficients")}.
#' @param values a numeric vector. The values the parameter estimates are compared to. Default is zero.
#' @param alternative  a character string. The alternative hypothesis. Options are \code{"two.sided"}, \code{"greater"} or \code{"less"}. Default is \code{"two.sided"}.
#' @param digits an integer. The number of digits of the result.
#' @export
#' @return A \code{data.frame} with one row for each tested parameter and 6 variables:
#' \item{estimate}{the parameter estimate.}
#' \item{diff}{the difference between the estimated parameter and the numeric value.}
#' \item{std.error}{the standard error of the estimated parameter.}
#' \item{t.value}{the t statistic.}
#' \item{res.degrees}{the residual degrees of freedom of the model.}
#' \item{p.value}{the p value of the t statistic.}
#' @details The test function of the package.
#' @references
#' Wang Z, Xu B, Zhou HJ (2014). "Social Cycling and Conditional Responses in the Rock-Paper-Scissors Game." \emph{Scientific Reports}, 4(1), 2045-2322.
#' @examples
#' ## Test if the choice probabilities of a mixed strategy for rock-paper-scissors.
#' ## The rock-paper-scissors data is from Wang, Xu, and Zhou (2014).
#' model.mixed <- stratEst.model(data = data.WXZ2014, strategies = strategies.RPS["mixed"])
#' t.probs <- stratEst.test(model = model.mixed, par = "probs", values = 1/3)
#' print(t.probs)
#' @export
stratEst.test <- function( model, par = c("shares","probs","trembles","coefficients"), values = 0, alternative = "two.sided", digits = 4 ){

  # checks
  if( "stratEst.model" %in% class(model) == F ){
    stop("stratEst.test error: The object passed to the argument 'model' must be of class 'stratEst.model'.")
  }
  if( "character" %in% class(par) == F ){
    stop("stratEst.test error: The object passed to the argument 'par' must be of a character string or vector.")
  }else{
    for( i in 1:length(par) ){
      if( par[i] %in% c("shares","probs","trembles","coefficients") == F ){
        stop("stratEst.test error: The object passed to the argument 'par' can only contain the following chracters: 'shares','probs','trembles','coefficients'.")
      }
    }
  }
  if( "numeric" %in% class(values) == F ){
    stop("stratEst.test error: The object passed to the argument 'values' must be numeric.")
  }
  if( "numeric" %in% class(values) == F | length(digits) != 1 ){
    stop("stratEst.test error: The object passed to the argument 'digits' must be a positive integer.")
  }

  par_matrix <- NULL
  est <- NULL
  se <- NULL
  row_names <- NULL

  if( length(model$shares.par > 0) & is.null(model$coefficients) & "shares" %in% par ){
    est <- c(est,model$shares.par)
    ses <- model$shares.se
    ses[ ses == 0 ] = NA
    se <- c(se,ses)
    row_names <- c(row_names,paste("shares.par.",as.character(seq(1,length(model$shares.par),by = 1)),sep=""))
  }
  if( is.null(model$coefficients.par) == F  & "coefficients" %in% par ){
    est <- c(est,model$coefficients.par)
    ses <- model$coefficients.se
    ses[ ses == 0 ] = NA
    se <- c(se,ses)
    row_names <- c(row_names,paste("coefficients.par.",as.character(seq(1,length(model$coefficients.par),by = 1)),sep=""))
  }
  if( length(model$probs.par > 0) & length(model$probs.se > 0)  & "probs" %in% par ){
    est <- c(est,model$probs.par)
    ses <- model$probs.se
    ses[ ses == 0 ] = NA
    se <- c(se,ses)
    row_names <- c(row_names,paste("probs.par.",as.character(seq(1,length(model$probs.par),by = 1)),sep=""))
  }
  if( length(model$trembles.par > 0)  & "trembles" %in% par ){
    est <- c(est,model$trembles.par)
    ses <- model$trembles.se
    ses[ ses == 0 ] = NA
    se <- c(se,ses)
    row_names <- c(row_names,paste("trembles.par.",as.character(seq(1,length(model$trembles.par),by = 1)),sep=""))
  }
    diff <- est-values
    z <- diff/se
    # p-value
    if( alternative == "two.sided" ){
      p <- 2*stats::pt( abs(z) , model$res.degrees , lower = F )
    }else if( alternative == "greater" ){
      p <- stats::pt( z , model$res.degrees , lower = F )
    }else if( alternative == "smaller" ){
      p <- stats::pt( z , model$res.degrees , lower = T )
    }else{
      stop("stratEst error: The argument 'alternative' must be one of the following chracter strings: 'two.sided','greater' or 'smaller'.")
    }
    test_matrix = cbind( est , diff , se , z , rep(model$res.degrees, length(est)) , p )
    colnames(test_matrix) <- c("estimate","diff","std.error","t.value","df","p.value")
    rownames(test_matrix) <- row_names
    par_data <- data.frame(round(test_matrix,digits))
    colnames(par_data) <- c("estimate","diff","std.error","t.value","df","p.value")

  return(par_data)

}

.onUnload <- function (libpath) { library.dynam.unload("stratEst", libpath)}
