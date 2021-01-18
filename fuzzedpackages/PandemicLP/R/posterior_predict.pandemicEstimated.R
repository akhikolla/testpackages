# Environment created to exchange objects between package functions
pandemic_environment = new.env()
#pandemic_environment$fullPred = list()

#' Draw from the posterior predictive distribution for pandemic data
#'
#' The posterior predictive distribution is the distribution of the outcome
#' implied by the model after using the observed data to update our beliefs
#' about the unknown parameters in the model. Simulating data from the posterior
#' predictive distribution using the observed predictors is useful for checking
#' the fit of the model. Drawing from the posterior predictive distribution at
#' interesting values of the predictors also lets us visualize how a manipulation
#' of a predictor affects (a function of) the outcome(s). With new observations of
#' predictor variables we can use the posterior predictive distribution to generate
#' predicted outcomes.
#' @method posterior_predict pandemicEstimated
#' @param object An object of class \code{pandemicEstimated} created by function \code{\link{pandemic_model}}.
#' @param horizonLong How far into the future the long-term prediction is desired.
#' @param horizonShort How far into the future the short-term prediction is desired.
#' @param ... Currently unused.
#' @references
#' CovidLP Team, 2020. CovidLP: Short and Long-term Prediction for COVID-19. Departamento de Estatistica. UFMG,
#' Brazil. URL: \url{http://est.ufmg.br/covidlp/home/en/}
#' @return An object of class \code{pandemicPredicted}. It includes the sampled predictive distribution
#' the model used to predict, which is the same as the one used to estimate the data. This object can be used
#' directly into the plot function and contains the following elements:
#' \item{\code{predictive_Long}}{
#'   A \code{M x horizonLong} matrix with the full sample of the predictive distribution
#'   for the long-term prediction, where M is the sample size.
#'   The prediction is for daily new cases.
#'   }
#'   \item{\code{predictive_Short}}{
#'   A \code{M x horizonShort} matrix with the full sample of the predictive distribution
#'   for the short-term prediction, where M is the sample size.
#'   The prediction is for daily cumulative cases.
#'   }
#'   \item{\code{data}}{
#'   The data passed on from the \code{\link{pandemicEstimated-objects}} under the element \code{Y$data}.
#'   }
#'   \item{\code{location}}{
#'   A string with the name of the location.
#'   }
#'   \item{\code{cases_type}}{
#'   A string with either "confirmed" or "deaths" to represent the type of data that has been fitted and predicted.
#'   }
#'   \item{\code{pastMu}}{
#'   The fitted means of the data for the observed data points.
#'   }
#'   \item{\code{futMu}}{
#'   The predicted means of the data for the predicted data points.
#'   }
#' Function \code{\link{pandemic_stats}} provides a few useful statistics based on the predictions.
#'
#' @seealso \code{\link{pandemic_model}}, \code{\link{pandemic_stats}} and \code{\link{plot.pandemicPredicted}}. Details
#' about the models behind the calculations can be seen in \code{\link{models}}.
#' @examples
#' \dontrun{
#' dataMG = load_covid("Brazil","MG")
#' estimMG = pandemic_model(dataMG)
#' predMG = posterior_predict(estimMG)
#' predMG}
#' @importFrom rstantools posterior_predict
#' @importFrom methods slot
#' @export
posterior_predict.pandemicEstimated = function(object,horizonLong = 500, horizonShort = 14,...){

  if (class(object) != "pandemicEstimated") stop("Please use the output of the pandemic_model() function.")
  if (horizonShort <= 0) stop("Horizons must be positive.")
  if (horizonLong < horizonShort) stop("Long-term horizon may not be lesser than short term horizon.")
  #if (horizonLong > 1000) stop("Long-term horizon larger than 1000 is not currently supported.")

  chains = as.data.frame(object$fit)
  pop = object$Y$population
  longPred = max(1000,horizonLong)

  M = nrow(chains) ## Total iterations
  NA_replacement = 2*object$Y$population ## Set a NA replacement

  finalTime = sum(grepl("mu",names(chains))) ## How many mu's

  # generate points from the marginal predictive distribution
  if (grepl("seasonal",object$model_name)) s_code=c(seasonal_code(object$Y$data$date,object$seasonal_effect),0,0) else s_code = NULL
  pred = generatePredictedPoints_pandemic(M,chains,longPred, NA_replacement, object$model_name, finalTime, s_code)
  methods::slot(object$fit,"sim")$fullPred = list() # For internal use
  methods::slot(object$fit,"sim")$fullPred$thousandLongPred = pred$yL
  if (object$cases.type == "confirmed")
    methods::slot(object$fit,"sim")$fullPred$thousandShortPred = pred$yS + object$Y$data$cases[nrow(object$Y$data)]
  else
    methods::slot(object$fit,"sim")$fullPred$thousandShortPred = pred$yS + object$Y$data$deaths[nrow(object$Y$data)]
  methods::slot(object$fit,"sim")$fullPred$thousandMus = pred$mu
  y.futL = as.matrix(methods::slot(object$fit,"sim")$fullPred$thousandLongPred[,1:horizonLong])
  y.futS = as.matrix(methods::slot(object$fit,"sim")$fullPred$thousandShortPred[,1:horizonShort])
  errorCheck = which(methods::slot(object$fit,"sim")$fullPred$thousandShortPred[,ncol(methods::slot(object$fit,"sim")$fullPred$thousandShortPred)] > pop)
  if (length(errorCheck)){
    message(paste0(length(errorCheck)," samples were removed from the prediction due to unrealistic results."))
    methods::slot(object$fit,"sim")$fullPred$thousandShortPred = methods::slot(object$fit,"sim")$fullPred$thousandShortPred[-errorCheck,]
    methods::slot(object$fit,"sim")$fullPred$thousandLongPred = methods::slot(object$fit,"sim")$fullPred$thousandLongPred[-errorCheck,]
    y.futL = y.futL[-errorCheck,]
    y.futS = y.futS[-errorCheck,]
  }

  output <- list(predictive_Long = y.futL, predictive_Short = y.futS,
                 data = object$Y$data, location = object$Y$name, cases_type = object$cases.type,
                 pastMu = as.data.frame(object$fit)[grep("mu",names(chains))],
                 futMu = as.matrix(pred$mu[,1:horizonLong]),fit = object$fit,seasonal_effect = object$seasonal_effect,
                 errors = errorCheck)

  class(output) = "pandemicPredicted"
  return(output)
}

# Auxiliary function for posterior_predict.pandemicEstimated
# Provides the actual predictions
# c for chains, h for horizon, n for NA value, m for model, ft for final time and s for seasonal effects days
#' @importFrom stats rpois rgamma pnorm
generatePredictedPoints_pandemic = function(M,c,h,n,m,ft,s){

  y = mu = matrix(-Inf,ncol = h,nrow = M)
  if (m == "poisson: static generalized logistic")
    for (i in 1:h){
      mu[,i] = exp(log(c$f)+log(c$a)+log(c$c)-(c$c*(ft+i))-(c$f+1)*log(c$b+exp(-c$c*(ft+i)) ))
      y[,i] = stats::rpois(M,mu[,i])
    } else if (m == "negbin: static generalized logistic")
    for (i in 1:h){
      mu[,i] = exp(log(c$f)+log(c$a)+log(c$c)-(c$c*(ft+i))-(c$f+1)*log(c$b+exp(-c$c*(ft+i)) ))
      y[,i] = stats::rpois(M,stats::rgamma(M,mu[,i]*c$phi,c$phi))
    } else if (m == "poisson: static seasonal generalized logistic"){
      d_1 = c$d_1; d_2 = ifelse(is.null(c$d_2),1,c$d_2); d_3 = ifelse(is.null(c$d_3),1,c$d_3)
    for (i in 1:h){
      mu[,i] = exp(log(c$f)+log(c$a)+log(c$c)-(c$c*(ft+i))-(c$f+1)*log(c$b+exp(-c$c*(ft+i)) ))*
        d_1^((s[1]>0) * !((ft+i - s[1]) %% 7)) * d_2^((s[2]>0) * !((ft+i - s[2]) %% 7)) * d_3^((s[3]>0) * !((ft+i - s[3]) %% 7))
      y[,i] = stats::rpois(M,mu[,i])
    }} else if (m == "poisson: multi_waves(2)")
    for (i in 1:h){
      mu[,i] = exp(log(c$a1)+log(c$c1)-(c$c1*(ft+i))-2*log(c$b1+exp(-c$c1*(ft+i)) ) + stats::pnorm(c$alpha1*(ft+i-c$delta1),log.p = TRUE))+
          exp(log(c$a2)+log(c$c2)-(c$c2*(ft+i))-2*log(c$b2+exp(-c$c2*(ft+i)) ) + stats::pnorm(c$alpha2*(ft+i-c$delta2),log.p = TRUE))
      y[,i] = stats::rpois(M,mu[,i])
    } else stop(paste("Unknown model",m))

  if (any(is.na(y))){
    message(paste("Prediction had",sum(is.na(y)),"NA values. Replaced with large value for identification."))
    y.fut[is.na(y.fut)] = n
  }

  list(yL=y, yS = t(apply(y,1,cumsum)), mu = mu)
}
