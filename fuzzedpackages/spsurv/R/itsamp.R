#' Inverse Transform Sampling To Generate Time-to-event Data From Parametric Models
#'
#' @export
#' @description Random survival times generation for the weibull or
#'  log-logistic distributions with parameters `scale` and `shape`.
#' @param n integer; sample size
#' @param beta vector of regression coefficients
#' @param event_scale,censor_scale  event and censoring scale parameters
#' @param features matrix of features (columns)
#' @param shape event and censoring distribution shape
#' @param model either "ph" (default) or "aft" for weibull and
#' "po" or "aft" for log-logistic distribution
#' @param censor logical; if `TRUE`, censoring is required, that is
#' mean(status) > 0
#' @param dist "weibull" or "llogis"
#'
#' @details sim_surv returns weibull (log-logistic) randomly
#' generated survival times. According to Collett (2003), the
#' accelerated failure time model encompasses a wide variety of parametric
#' models, including weibull and log-logistic models.
#'
#' @return data.frame of `ncol(x) +2` columns in which the
#'  survival times are the response variable denoted by `y`,
#'   `status` indicates failure (0 = failure) and the features
#'   are appended to the next columns.
#'
#'@examples rows <- 200
#'
#' categorical <- rbinom(rows, size = 3, prob = .5)
#' x <- data.frame(numerical = rnorm(rows),
#'            cat0 = as.numeric(categorical == 0),
#'            cat1 = as.numeric(categorical == 1),
#'            cat2 = as.numeric(categorical == 2),
#'            cat3 = as.numeric(categorical == 3))
#'
#' newdata <- itsamp(n = rows, beta = c(1, -2, .5, .1, 1),
#'   features = x, model = 'ph', dist = 'weibull')
#'
#' @seealso \code{\link[spsurv]{spbp}}

itsamp <- function(n,
  beta = c(2, -1),
  event_scale = 10, # baseline hazard
  censor_scale = 4, # hazard of censoring
  features = data.frame(x1 = rnorm(n, 0), x2 = rnorm(n, 0)),
  shape = 2,
  model = c("ph", "po", "aft"),
  dist = c("weibull", "llogis"),
  censor = TRUE){

  if(n %% 1 != 0) stop("n must be integer")
  if(dim(features)[1] != n) stop("Lengths differ: dim(x)[1] = ", dim(features)[1],
         " but n = ", n)
  if(class(features) != "data.frame") stop("x must be data.frame")

  x <-  model.matrix(~., data = features)[, -1]

  eta <- x %*% beta
  model <- match.arg(model)
  dist <- match.arg(dist)

  if(dist == "weibull"){
    if(model == "po"){
      stop("no method available, try another dist or model.")
    }
    else{
      if(model == "ph"){     #proportional hazards
        scale <- (event_scale * exp(eta) ^ (-1/shape))
      }
      else if(model == "aft"){    #accelerated failure time
        scale <- event_scale * exp(-shape * eta) ^ (-1/shape)
      }
      t <- rweibull(n, scale = scale, shape = shape) #event time
    }
    c <- rweibull(n, scale = censor_scale, shape = shape) #censoring time
  }
  else{
    if(model == "ph"){
      stop("no method available, try another dist or model.")
    }
    else{
      if(model == "po"){      #proportional odds
        location = log(event_scale) + eta ## logis location parameter
        scale <- exp(location)
      }
      else if(model == "aft"){    #accelerated failure time
        location = log(event_scale) - shape * eta ## logis location parameter
        scale <- exp(location)
      }
      t <- exp(rlogis(n, location = -location/shape, scale = 1/shape)) #event time
    }
    c <- exp(rlogis(n, location = -log(censor_scale)/shape, scale = 1/shape)) #censoring time
  }

  y <- pmin(t,c)   #response
  status <- as.numeric(y == t)   # event indicator

  if(mean(status) == 0 && censor == TRUE){
    stop("No censoring was generated with the chosen arguments")
  }
  else{
    message("The censoring percentage is ",
            mean(status) * 100, "%, generated from ",
            model, " ", dist)
  }
  db <- data.frame(y, status, x)
  attr(db, "censoring") <- mean(status)
  attr(db, "model") <- model
  attr(db, "dist") <- dist
  attr(db, "scale") <- scale
  attr(db, "shape") <- shape
  return(db)
}
