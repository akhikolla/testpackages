#' Relevant Statistics of the Pandemic Model
#'
#' This function provides short and long-term predictions for the pandemic. 95\% credible intervals are
#' assigned for the number of cases for every future date predicted, as well as for the total number of cases,
#' and dates for the peak and end of the pandemic.\cr
#' \cr
#' Short-term predictions are made on the cumulative counts and long-term predictions
#' are based on the new case counts.
#'
#' @param object an object of S3 class \code{pandemicPredicted} created by function
#'  \code{\link{posterior_predict.pandemicEstimated}}.
#'
#' @return An object of S3 class \code{pandemicStats}. This object is a list containing the following elements:
#' \itemize{
#'    \item{\code{data}:}{
#'    A list with a data frame containing the observed pandemic data, a character string with the location name
#'    and a character string indicating the type of cases predicted.
#'    }
#'    \item{\code{ST_predict}:}{
#'    A data frame with the short-term predictions for the number of cumulative cases. For each future date predicted,
#'    the mean, median, 2.5 and 97.5 percentiles are provided.\cr
#'    The short-term horizon is determined by the \code{horizonShort} argument in the \code{posterior_predict} function.
#'    }
#'    \item{\code{LT_predict}:}{
#'    A data frame with the long-term predictions for the number of new cases. For each future date predicted,
#'    the mean, median, 2.5 and 97.5 percentiles are provided.\cr
#'    The long-term horizon is determined by the \code{horizonLong} argument in the \code{posterior_predict} function.
#'    }
#'    \item{\code{LT_summary}:}{
#'    A list with the estimated total number of cases and the dates for the peak and end of the pandemic.
#'    In each metric, the median, 2.5 and 97.5 percentiles are provided. For more information, see the
#'    \strong{Details} section.\cr
#'    }
#'    \item{\code{mu}:}{
#'    A data frame with the median values of the mean number of new cases for each date (starting from the first
#'    observed data point until the last date in the long-term horizon).
#'    }
#' }
#'
#' @details
#'  \subsection{Total Number of Cases}{
#'  The total number of cases is obtained by adding the cumulative total cases observed in the data to the
#'  predicted new cases for the next 1,000 days.
#'  }
#'  \subsection{Estimated Peak Dates}{
#'  The median, 2.5 and 97.5 percentiles are calculated on the mean number of new cases for the pandemic curve
#'  (starting from the first observed data point until 1,000 days after the last date observed in the data).\cr
#'  \cr
#'  The 95\% credible interval for the peak of cases is selected such that the two limiting dates of the 97.5
#'  percentile curve coincides with the highest value of the 2.5 percentile curve. This guarantees that all
#'  possible curves belonging to the confidence band will peak within the defined interval.
#'  }
#'  \subsection{End of the Pandemic Dates}{
#'  Represents the 99 percentile of the total number of cases in the pandemic.
#'  }
#'
#' @references
#' CovidLP Team, 2020. CovidLP: Short and Long-term Prediction for COVID-19. Departamento de Estatistica. UFMG,
#' Brazil. URL: \url{http://est.ufmg.br/covidlp/home/en/}
#'
#' @examples
#' \dontrun{
#' italy = load_covid("italy")
#' estim = pandemic_model(italy, case_type = "confirmed", covidLPconfig = TRUE)
#' pred = posterior_predict(estim)
#' stats = pandemic_stats(pred)
#' stats
#' }
#'
#' @seealso
#' \code{\link{load_covid}}, \code{\link{pandemic_model}}, \code{\link{posterior_predict.pandemicEstimated}}
#' and \code{\link{plot.pandemicPredicted}}.
#'
#' @importFrom stats quantile median
#' @importFrom methods slot
#' @export

pandemic_stats <- function(object){

  if (class(object) != "pandemicPredicted") stop("Please use the output of the posterior_predict() function.")
  if(missing(object)) stop("object is a required argument. See help(pandemic_stats) for more information.")

  t = length(object$data[[1]])
  ST_horizon = ncol(object$predictive_Short)
  LT_horizon = ncol(object$predictive_Long)
  longHorizon = ncol(methods::slot(object$fit,"sim")$fullPred$thousandLongPred)
  date_full <- as.Date(object$data$date[1]:(max(object$data$date) + longHorizon), origin = "1970-01-01")

  ### list output ST predicition:
  ST_predict <- data.frame( date  = date_full[(t+1):(t+ST_horizon)],
                            q2.5  = apply(object$predictive_Short,2,stats::quantile,.025),
                            med   = apply(object$predictive_Short,2,stats::median),
                            q97.5 = apply(object$predictive_Short,2,stats::quantile,.975),
                            mean  = colMeans(object$predictive_Short))
  row.names(ST_predict) <- NULL

  ### total number of cases
  if(object$cases_type == "confirmed"){
    cumulative_y =  object$data$cases[t]
  } else{
    cumulative_y =  object$data$deaths[t]
  }

  if(cumulative_y > 1000){
    lowquant <- apply(methods::slot(object$fit,"sim")$fullPred$thousandLongPred,2,stats::quantile,.025)
    medquant <- apply(methods::slot(object$fit,"sim")$fullPred$thousandLongPred,2,stats::median)
    highquant <- apply(methods::slot(object$fit,"sim")$fullPred$thousandLongPred,2,stats::quantile,.975)
  } else{
    lowquant <- c(cumulative_y , apply(methods::slot(object$fit,"sim")$fullPred$thousandShortPred,2,stats::quantile,.025))
    lowquant <- (lowquant - lag(lowquant, default = 0))[-1]
    medquant <- c(cumulative_y, apply(methods::slot(object$fit,"sim")$fullPred$thousandShortPred,2,stats::median))
    medquant <- (medquant - lag(medquant,default = 0))[-1]
    highquant <- c(cumulative_y, apply(methods::slot(object$fit,"sim")$fullPred$thousandShortPred,2,stats::quantile,.975))
    highquant <- (highquant - lag(highquant, default = 0))[-1]
  }

  TNC2.5 = sum(lowquant) + cumulative_y
  TNC50  = sum(medquant) + cumulative_y
  TNC97.5 = sum(highquant) + cumulative_y

  ### Calculates the peak and end dates
  peak2.5 <- peak50 <- peak97.5 <- NULL
  end2.5 <- end50 <- end97.5 <- NULL

  index_season <-NULL
  if(!is.null(object$seasonal_effect)){
    s_code <- seasonal_code(date_full, object$seasonal_effect)
    for (i in 1:length(s_code)){
      index_aux <- which((seq(1,(t+longHorizon),1) - s_code[i]) %% 7 == 0)
      index_season<-c(index_aux, index_season)
    }
    index_season<-sort(index_season)
  }


  chain_mu <- cbind(object$pastMu, methods::slot(object$fit,"sim")$fullPred$thousandMus)

  ### median dates
  mu50 <- apply(chain_mu, 2, stats::quantile, probs = 0.5)
  peak50 <- date_full[which.max(mu50)]

  q <- .99
  med_cumulative <- apply(as.matrix(mu50),2,cumsum)
  med_percent<- med_cumulative / med_cumulative[t + longHorizon]
  med_end <- which(med_percent - q > 0)[1]
  end50 <- date_full[med_end]

  ### calculates the upper and lower bound dates for the peak:
  mu25 <- apply(chain_mu, 2, stats::quantile, probs = 0.025)
  mu975 <- apply(chain_mu, 2, stats::quantile, probs = .975)

  mu25_aux <- if(is.null(object$seasonal_effect)) mu25 else mu25[-index_season]
  posMaxq2.5 <- which.max(mu25_aux)
  aux <- if (is.null(object$seasonal_effect)) mu975 - mu25_aux[posMaxq2.5] else  mu975[-index_season] - mu25_aux[posMaxq2.5]
  aux2 <- aux[ posMaxq2.5 : length(aux)]
  val <- ifelse( length(aux2[aux2 < 0]) > 0, min(aux2[aux2 > 0]), aux[length(aux)])
  date_max <- which(aux == val)

  aux <- if(is.null(object$seasonal_effect)) mu975 - mu25_aux[posMaxq2.5] else mu975[-index_season] - mu25_aux[posMaxq2.5]
  aux2 <- aux[1:posMaxq2.5]
  val <- min(aux2[aux2>0])
  date_min <- which(aux == val)

  date_full_aux <- if(is.null(object$seasonal_effect)) date_full else date_full[-index_season]
  peak2.5  <- date_full_aux[date_min]
  peak97.5 <- date_full_aux[date_max]

  ##calculates the upper and lower bound dates for the end of the pandemic:
  low_cumulative <- apply(as.matrix(mu25),2,cumsum)
  low_percent <- low_cumulative / low_cumulative[t + longHorizon]
  low_end <- which(low_percent - q > 0)[1]
  end2.5 <- date_full[low_end]

  high_cumulative <- apply(as.matrix(mu975),2,cumsum)
  high_percent <- high_cumulative / high_cumulative[t + longHorizon]
  high_end <- which( high_percent - q > 0)[1]
  end97.5 <- date_full[high_end]

  #### LT predictions:
  LT_predict <- data.frame( date  = date_full[(t+1):(t+LT_horizon)],
                            q2.5  = lowquant[1:LT_horizon],
                            med   = medquant[1:LT_horizon],
                            q97.5 = highquant[1:LT_horizon],
                            mean  = colMeans(object$predictive_Long))
  row.names(LT_predict) <- NULL

  ## Long-term summary
  LT_summary <- list(total_cases_LB = TNC2.5,
                     total_cases_med = TNC50,
                     total_cases_UB = TNC97.5,
                     peak_date_LB = peak2.5 ,
                     peak_date_med = peak50 ,
                     peak_date_UB = peak97.5,
                     end_date_LB = end2.5,
                     end_date_med = end50,
                     end_date_UB = end97.5)

  muplot <- data.frame(date = date_full[1:(t+LT_horizon)], mu = mu50[1:(t+LT_horizon)])
  row.names(muplot) <-NULL

  dataplot <- list(data = object$data, location = object$location, case_type = object$cases_type)

  output <- list( data = dataplot, ST_predict = ST_predict, LT_predict = LT_predict,
                  LT_summary = LT_summary, mu = muplot )

  class(output) = "pandemicStats"
  return(output)

}
