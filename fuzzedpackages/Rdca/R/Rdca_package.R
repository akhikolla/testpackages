

# *******************************Decline Model *********************************


#' Arps decline object
#'
#' Create an object of class 'decline'
#'
#' @param input_unit a unit system for parameters, a character string either 'SI' or 'Field'
#' @param output_unit a unit system for properties, a character string either 'SI' or 'Field'
#' @param fluid fluid type, a character string either 'oil' or 'gas'
#' @param model decline model, a character string. 'exponential', 'harmonic', 'hyperbolic', and 'modified_hyperbolic' models are currently available
#' @param qi Arp's decline parameter, a numeric value. Depending on the 'input_unit', 'fluid' type, and also the decline_time() 'unit' parameter, it has different units. 'm3/day' for gas production in 'SI' unit with daily data, 'm3/month' for gas production in 'SI' unit with monthly data, 'MSCF/day' for gas production in 'Field' unit with daily data, 'MSCF/month' for gas production in 'Field' unit with monthly data, 'm3/day' for oil production in 'SI' unit with daily data, 'm3/month' for oil production in 'SI' unit with monthly data, 'bbl/day' for oil production in 'Field' unit with daily data, and 'bbl/month' for oil production in 'Field' unit with monthly data
#' @param Di Arp's nominal decline parameter, a numeric value in '1/day', '1/month', or '1/year' depending on the decline_time() 'unit' parameter
#' @param b Arp's decline parameter. It is zero for the 'exponential' model, one for the 'harmonic' model, and a value between zero and one for the 'hyperbolic' model. For unconventional reservoirs, b values more than one are also reported
#' @param Dt Arp's "modified_hyperbolic" nominal terminal decline parameter, a numeric value in '1/day', '1/month', or '1/year' depending on the decline_time() 'unit' input
#' @param q_abnd abandonment rate, a numeric value defaulted to NULL. If present, the model predicts the time to reach to the abandonment rate and also the estimated ultimate recovery (EUR) till the abandonment time. It has the same unit as 'qi'.
#'
#' @return a list of class 'decline' with all the required parameters for the decline_predict() S3 methods
#'
#' @export
#'
#' @examples
#' decline_param_1 <- decline_param(input_unit = "Field", output_unit = "Field", fluid = "oil",
#' model = "exponential", qi = 1000, Di = 0.15, b = 0, q_abnd = NULL)
#'
#' decline_param_1
#'
#' decline_param_2 <- decline_param(input_unit = "Field", output_unit = "SI", fluid = "oil",
#' model = "hyperbolic", qi = 15000, Di = 0.1, b = 0.95, q_abnd = 200)
#'
#' decline_param_2
#'
#' decline_param_3 <- decline_param(input_unit = "Field", output_unit = "Field", fluid = "gas",
#' model = "modified_hyperbolic", qi = 100000, Di = 0.15, b = 0.85, Dt = 0.005, q_abnd = NULL)
#'
#' decline_param_3


decline_param <- function(input_unit = "Field", output_unit = "Field", fluid = "oil", model = "exponential", qi = NULL, Di = NULL, b = NULL, Dt = NULL, q_abnd = NULL) {

   if (!is.character(input_unit)) stop("'input_unit' must be a character string.")
   if (!is.character(output_unit)) stop("'output_unit' must be a character string.")
   if (!is.character(fluid)) stop("'fluid' must be a character string.")
   if (!is.character(model)) stop("'model' must be a character string.")
   if (!(input_unit %in% c("SI","Field"))) stop("'input_unit' must be either 'SI' or 'Field'.")
   if (!(output_unit %in% c("SI","Field"))) stop("'output_unit' must be either 'SI' or 'Field'.")
   if (!(fluid %in% c("oil","gas"))) stop("'fluid' must be either 'oil' or 'gas'.")
   if (!(model %in% c("exponential","harmonic",
                       "hyperbolic", "modified_hyperbolic"))) stop("DCA model is not in the list of available models.")
   if (length(model) > 1) stop("Only one model is acceptable for the decline analysis.")

   if (model == "exponential") {

      if (!is.numeric(qi)) stop("'qi' must be a numeric.")
      if (!is.numeric(Di)) stop("'Di' must be a numeric.")
      if (!is.numeric(b)) stop("'b' must be a numeric.")
      if (b != 0) stop("The value of 'b' is zero in the 'exponential' model.")
      if (is.null(q_abnd)) {
         decline_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, qi = qi, Di = Di, b = b)
      } else {
         if (!is.numeric(q_abnd)) stop("'q_abnd' must be a numeric or NULL.")
         decline_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, qi = qi, Di = Di, b = b, q_abnd = q_abnd)
      }
      class(decline_lst) <- c("exponential", "decline")
   }

   if (model == "harmonic") {

      if (!is.numeric(qi)) stop("'qi' must be a numeric.")
      if (!is.numeric(Di)) stop("'Di' must be a numeric.")
      if (!is.numeric(b)) stop("'b' must be a numeric.")
      if (b != 1) stop("The value of 'b' is one in the 'harmonic' model.")
      if (is.null(q_abnd)) {
         decline_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, qi = qi, Di = Di, b = b)
      } else {
         if (!is.numeric(q_abnd)) stop("'q_abnd' must be a numeric or NULL.")
         decline_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, qi = qi, Di = Di, b = b, q_abnd = q_abnd)
      }
      class(decline_lst) <- c("harmonic", "decline")
   }

   if (model == "hyperbolic") {

      if (!is.numeric(qi)) stop("'qi' must be a numeric.")
      if (!is.numeric(Di)) stop("'Di' must be a numeric.")
      if (!is.numeric(b)) stop("'b' must be a numeric.")
      if (b == 0) stop("'b' must be greater than zero for the 'hyperbolic' model.")
      if (b == 1) stop("The value of 'b' must be less than one.")
      if (is.null(q_abnd)) {
         decline_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, qi = qi, Di = Di, b = b)
      } else {
         if (!is.numeric(q_abnd)) stop("'q_abnd' must be a numeric or NULL.")
         decline_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, qi = qi, Di = Di, b = b, q_abnd = q_abnd)
      }
      class(decline_lst) <- c("hyperbolic", "decline")
   }

   if (model == "modified_hyperbolic") {

      if (!is.numeric(qi)) stop("'qi' must be a numeric.")
      if (!is.numeric(Di)) stop("'Di' must be a numeric.")
      if (!is.numeric(b)) stop("'b' must be a numeric.")
      if (!is.numeric(Dt)) stop("'Dt' must be a numeric.")
      if (b == 0) stop("'b' must be greater than zero for the 'modified-hyperbolic' model.")
      if (is.null(q_abnd)) {
         decline_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, qi = qi, Di = Di, b = b, Dt = Dt)
      } else {
         if (!is.numeric(q_abnd)) stop("'q_abnd' must be a numeric or NULL.")
         decline_lst <- list(input_unit = input_unit, output_unit = output_unit, fluid = fluid, qi = qi, Di = Di, b = b, Dt = Dt, q_abnd = q_abnd)
      }
      class(decline_lst) <- c("modified_hyperbolic", "decline")
   }
   return(decline_lst)
}



# ******************************** Time-Date Model *****************************

#' Arps time object
#'
#' Create an object of class 'time'
#'
#' @param x a vector of times or a daily sequence of dates
#' @param unit time/date unit of vector x
#'
#' @return a list of class 'time' with all the required parameters for the decline_predict() S3 methods
#'
#' @export
#'
#' @examples
#' decline_time_1 <- decline_time(c(1:730), unit = "day")
#'
#' decline_time_1
#'
#' decline_time_2 <- decline_time(c(1:240), unit = "month")
#'
#' decline_time_2
#'
#' decline_time_3 <- decline_time(c(1:15), unit = "year")
#'
#' decline_time_3
#'
#' decline_time_4 <- decline_time(seq(as.Date("2020/1/1"), by = "days",
#' length.out = 360), unit = "date")
#'
#' decline_time_4


decline_time <- function(x, unit = "day") {

   if (!is.character(unit)) stop("'unit' must be a character string of type 'day', 'month', 'year', or 'date'.")
   if (!(unit %in% c("day", "month", "year", "date"))) stop("'unit' must be a character string of type 'day', 'month', 'year', or 'date'.")
   if (unit == "date") {
      if (class(x) != "Date") stop("'x' must be a sequence of type 'Date'.")
   }

   if (unit %in% c("day", "month", "year")) {
      if (!is.vector(x)) stop("'x' must be a vector of days, months, or years.")
      if (!is.numeric(x)) stop("'x' must be a numeric vector of days, months, or years.")
      if (any(duplicated(x))) stop("'x' must be a non-duplicated numeric vector of days, months, or years.")
      if (is.unsorted(x)) stop("'x' must be sorted in an ascending order.")
      time_lst <- list(t = x, unit = unit, reference_date = Sys.Date())
      class(time_lst) <- c(unit, "time")
   }
   if (unit == "date") {
      x <- as.Date(x)
      if (any(duplicated(x))) stop("'x' must be a non-duplicated sequence of dates.")
      if (is.unsorted(x)) stop("'x' must be sorted in an ascending order.")
      time_day <- as.numeric(x - x[1] + 1)
      time_lst <- list(t = time_day, unit = unit, reference_date = x[1])
      class(time_lst) <- c("day", "time")
   }
   return(time_lst)
}



# *****************************Decline-Predict Model ***************************

#' Arps decline prediction
#'
#' Create a data frame of decline predictions according to the class of 'decline_lst' and 'time_lst' objects
#'
#' @param decline_lst a list object of class 'decline'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame of decline estimates according to the class of 'decline_lst' and 'time_lst' objects
#'
#' @export
#'
#' @examples
#' decline_param_1 <- decline_param(input_unit = "Field", output_unit = "Field",
#' fluid = "oil",
#' model = "exponential", qi = 1000, Di = 0.15, b = 0, q_abnd = NULL)
#' decline_time_1 <- decline_time(c(1:7300), unit = "day")
#' decline_predict_1 <- decline_predict(decline_param_1, decline_time_1)
#'
#' head(decline_predict_1)
#'
#' decline_param_2 <- decline_param(input_unit = "Field", output_unit = "SI",
#' fluid = "oil",
#' model = "hyperbolic", qi = 15000, Di = 0.1, b = 0.95, q_abnd = 200)
#' decline_time_2 <- decline_time(seq(as.Date("2016/04/15"), by = "days",
#' length.out = 3600), unit = "date")
#' decline_predict_2 <- decline_predict(decline_param_2, decline_time_2)
#'
#' head(decline_predict_2)
#'
#' @references
#' \insertRef{Arps1945}{Rdca}
#'
#' \insertRef{Robertson1988}{Rdca}
#'

decline_predict <- function(decline_lst, time_lst) {
   if ((inherits(decline_lst, "decline") == TRUE) & (inherits(time_lst, "time"))) {
      UseMethod("decline_predict")
   } else {
      if (!inherits(decline_lst, "decline")) {
         stop("A class of 'decline' must be assigned to the 'decline_lst' parameter of the decline_predict() function.")
      }
      if (!inherits(time_lst, "time")) {
         stop("A class of 'time' must be assigned to the 'time_lst' parameter of the decline_predict() function.")
      }
   }
}


#' S3 method for class 'decline_predict'
#'
#' Create a data frame of exponential decline predictions
#'
#' @param decline_lst a list object of class 'decline'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame of decline estimates using the Arps exponential model
#'
#' @export
#'
#' @references
#' \insertRef{Arps1945}{Rdca}

decline_predict.exponential <- function(decline_lst, time_lst) {
   month_to_day <- 30.41667
   year_to_day <- 365
   Date <- NULL
   `Time_(day)` <- NULL
   `Time_(month)` <- NULL
   `Time_(year)` <- NULL
   results <- as.data.frame(decline_predict_cpp(decline_lst, time_lst))
   if (time_lst$unit == "date") {
      results$Date <- (results$`Time_(day)` - 1) + time_lst$`reference_date`
      results <- results %>% dplyr::select(Date, dplyr::everything())
   }
   if (time_lst$unit == "day") {
      results$Date <- (results$`Time_(day)` - 1) + time_lst$`reference_date`
      results <- results %>% dplyr::select(Date, dplyr::everything())
   }
   if (time_lst$unit == "month") {
      results$Date <- (results$`Time_(month)` - 1) * month_to_day + time_lst$`reference_date`
      results <- results %>% dplyr::select(Date, `Time_(month)`, dplyr::everything())
   }
   if (time_lst$unit == "year") {
      results$Date <- (results$`Time_(year)` - 1) * year_to_day + time_lst$`reference_date`
      results <- results %>% dplyr::select(Date, `Time_(year)`, dplyr::everything())
   }
   return(results)
}



#' S3 method for class 'decline_predict'
#'
#' Create a data frame of harmonic decline predictions
#'
#' @param decline_lst a list object of class 'decline'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame of decline estimates using the Arps harmonic model
#'
#' @export
#'
#' @references
#' \insertRef{Arps1945}{Rdca}

decline_predict.harmonic <- function(decline_lst, time_lst) {
   month_to_day <- 30.41667
   year_to_day <- 365
   Date <- NULL
   `Time_(day)` <- NULL
   `Time_(month)` <- NULL
   `Time_(year)` <- NULL
   results <- as.data.frame(decline_predict_cpp(decline_lst, time_lst))
   if (time_lst$unit == "date") {
      results$Date <- (results$`Time_(day)` - 1) + time_lst$`reference_date`
      results <- results %>% dplyr::select(Date, dplyr::everything())
   }
   if (time_lst$unit == "day") {
      results$Date <- (results$`Time_(day)` - 1) + time_lst$`reference_date`
      results <- results %>% dplyr::select(Date, dplyr::everything())
   }
   if (time_lst$unit == "month") {
      results$Date <- (results$`Time_(month)` - 1) * month_to_day + time_lst$`reference_date`
      results <- results %>% dplyr::select(Date, `Time_(month)`, dplyr::everything())
   }
   if (time_lst$unit == "year") {
      results$Date <- (results$`Time_(year)` - 1) * year_to_day + time_lst$`reference_date`
      results <- results %>% dplyr::select(Date, `Time_(year)`, dplyr::everything())
   }
   return(results)
}


#' S3 method for class 'decline_predict'
#'
#' Create a data frame of hyperbolic decline predictions
#'
#' @param decline_lst a list object of class 'decline'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame of decline estimates using the Arps hyperbolic model
#'
#' @export
#'
#' @references
#' \insertRef{Arps1945}{Rdca}

decline_predict.hyperbolic <- function(decline_lst, time_lst) {
   month_to_day <- 30.41667
   year_to_day <- 365
   Date <- NULL
   `Time_(day)` <- NULL
   `Time_(month)` <- NULL
   `Time_(year)` <- NULL
   results <- as.data.frame(decline_predict_cpp(decline_lst, time_lst))
   if (time_lst$unit == "date") {
      results$Date <- (results$`Time_(day)` - 1) + time_lst$`reference_date`
      results <- results %>% dplyr::select(Date, dplyr::everything())
   }
   if (time_lst$unit == "day") {
      results$Date <- (results$`Time_(day)` - 1) + time_lst$`reference_date`
      results <- results %>% dplyr::select(Date, dplyr::everything())
   }
   if (time_lst$unit == "month") {
      results$Date <- (results$`Time_(month)` - 1) * month_to_day + time_lst$`reference_date`
      results <- results %>% dplyr::select(Date, `Time_(month)`, dplyr::everything())
   }
   if (time_lst$unit == "year") {
      results$Date <- (results$`Time_(year)` - 1) * year_to_day + time_lst$`reference_date`
      results <- results %>% dplyr::select(Date, `Time_(year)`, dplyr::everything())
   }
   return(results)
}



#' S3 method for class 'decline_predict'
#'
#' Create a data frame of modified_hyperbolic decline predictions
#'
#' @param decline_lst a list object of class 'decline'
#' @param time_lst a list object of class 'time'
#'
#' @return a data frame of decline estimates using the Arps modified_hyperbolic model
#'
#' @export
#'
#' @references
#' \insertRef{Arps1945}{Rdca}

decline_predict.modified_hyperbolic <- function(decline_lst, time_lst) {
   month_to_day <- 30.41667
   year_to_day <- 365
   Date <- NULL
   `Time_(day)` <- NULL
   `Time_(month)` <- NULL
   `Time_(year)` <- NULL
   results <- as.data.frame(decline_predict_cpp(decline_lst, time_lst))
   if (time_lst$unit == "date") {
      results$Date <- (results$`Time_(day)` - 1) + time_lst$`reference_date`
      results <- results %>% dplyr::select(Date, dplyr::everything())
   }
   if (time_lst$unit == "day") {
      results$Date <- (results$`Time_(day)` - 1) + time_lst$`reference_date`
      results <- results %>% dplyr::select(Date, dplyr::everything())
   }
   if (time_lst$unit == "month") {
      results$Date <- (results$`Time_(month)` - 1) * month_to_day + time_lst$`reference_date`
      results <- results %>% dplyr::select(Date, `Time_(month)`, dplyr::everything())
   }
   if (time_lst$unit == "year") {
      results$Date <- (results$`Time_(year)` - 1) * year_to_day + time_lst$`reference_date`
      results <- results %>% dplyr::select(Date, `Time_(year)`, dplyr::everything())
   }
   return(results)
}




# ****************************** CURVE FITTING *********************************


#' Arps decline_fit object
#'
#' Create an object of class 'decline_fit'
#'
#' @param input_unit a unit system for parameters, a character string either 'SI' or 'Field'
#' @param output_unit a unit system for properties, a character string either 'SI' or 'Field'
#' @param fluid fluid type, a character string either 'oil' or 'gas'
#' @param model decline model, a character string. 'exponential', 'harmonic', 'hyperbolic', and 'modified_hyperbolic' models are currently available
#' @param fit_data a character string, either 'rate', or 'cum'
#' @param prod_data a numeric vector of rates or cumulative according to 'fit_data'
#' @param initial_param a numeric vector of initial estimates for the Arps decline model
#' @param lower an optional numeric vector of lower bounds for the Arps decline model parameters. See 'minpack.lm' package for details
#' @param upper an optional numeric vector of upper bounds for the Arps decline model parameters. See 'minpack.lm' package for details
#' @param control an optional list of control settings. See 'minpack.lm' package for details
#'
#' @return a list of class 'decline_fit' with all the required parameters for the decline_fit() S3 methods
#'
#' @export
#'
#' @examples
#'
#' prod_data <- 3000 * exp(-0.00234 * c(1:300))
#' dcl_fit_param_exp <- decline_fit_param(input_unit = "Field", output_unit = "Field",
#' fluid = "oil", model = "exponential", fit_data = "rate", prod_data = prod_data,
#' initial_param = c(1000, 0.1, 0), lower = NULL, upper = NULL, control = NULL)
#'
#' dcl_fit_param_exp
#'
#' prod_data <- 4500 / (1 + 0.002 * 0.834 * c(1:400)) ^ (1 / 0.834)
#' dcl_fit_param_mod_hyp <- decline_fit_param(input_unit = "Field", output_unit = "Field",
#' fluid = "oil", model = "modified_hyperbolic", fit_data = "rate", prod_data = prod_data,
#' initial_param = c(10000, 0.1, 0.8, 0.01), lower = NULL,upper = NULL,
#' control = list(maxiter = 100))
#'
#' dcl_fit_param_mod_hyp

decline_fit_param <- function(input_unit = "Field", output_unit = "Field", fluid = "oil", model = "exponential", fit_data = "rate", prod_data, initial_param, lower = NULL, upper = NULL, control = NULL) {

   if (!is.character(input_unit)) stop("'input_unit' must be a character string either 'SI' or 'Field'.")
   if (!is.character(output_unit)) stop("'output_unit' must be a character string either 'SI' or 'Field'.")
   if (!is.character(fluid)) stop("'fluid' must be a character string.")
   if (!is.character(model)) stop("'model' must be a character string.")
   if (!is.character(fit_data)) stop("'fit_data' must be a character string.")
   if (!is.numeric(prod_data)) stop("'prod_data' must be a numeric vector.")
   if (!is.numeric(initial_param)) stop("'initial_param' must be a numeric vector.")
   if (!(is.numeric(lower) | is.null(lower))) stop("'lower' must be a numeric vector or 'NULL'.")
   if (!(is.numeric(upper) | is.null(upper))) stop("'upper' must be a numeric vector or 'NULL'.")
   if (!(is.list(control) | is.null(control))) stop("'control' is an optional 'list' of control settings.")
   if (!(input_unit %in% c("SI","Field"))) stop("'input_unit' must be either 'SI' or 'Field'.")
   if (!(output_unit %in% c("SI","Field"))) stop("'output_unit' must be either 'SI' or 'Field'.")
   if (!(fluid %in% c("oil","gas"))) stop("'fluid' must be either 'oil' or 'gas'.")
   if (!(model %in% c("exponential","harmonic", "hyperbolic", "modified_hyperbolic"))) stop("DCA model is not in the list of available models.")
   if (length(model) > 1) stop("Only one model is acceptable for the decline optimization.")
   if (fit_data %in% c("rate", "cum") == FALSE) stop("'fit_data' must be either 'rate' or 'cum'.")
   if (length(fit_data) > 1) stop("'fit_data' is either 'rate' or 'cum'.")
   if (model == "exponential") {
      if (length(initial_param) != 3) {
         stop("'initial_param' must be a numeric vector of length 3.")
      }
   }
   if (model == "harmonic") {
      if (length(initial_param) != 3) {
         stop("'initial_param' must be a numeric vector of length 3.")
      }
   }
   if (model == "hyperbolic") {
      if (length(initial_param) != 3) {
         stop("'initial_param' must be a numeric vector of length 3.")
      }
   }
   if (model == "modified_hyperbolic") {
      if (length(initial_param) != 4) {
         stop("'initial_param' must be a numeric vector of length 4.")
      }
   }
   lst <- list(input_unit, output_unit, fluid, model, fit_data, prod_data, initial_param, lower, upper, control)
   names(lst) <- c("input_unit", "output_unit", "fluid", "model", "fit_data", "prod_data", "initial_parameters", "lower", "upper", "control")
   if (model == "exponential") {
      class(lst) <- c("exponential_fit", "decline_fit")
   }
   if (model == "harmonic") {
      class(lst) <- c("harmonic_fit", "decline_fit")
   }
   if (model == "hyperbolic") {
      class(lst) <- c("hyperbolic_fit", "decline_fit")
   }
   if (model == "modified_hyperbolic") {
      class(lst) <- c("modified_hyperbolic_fit", "decline_fit")
   }
   return(lst)
}



#' Arps decline_fit prediction
#'
#' Generate a list of estimates for the Arps decline model according to the class of 'decline_fit_lst' and 'time_lst' objects
#'
#' @param decline_fit_lst a list object of class 'decline_fit'
#' @param time_lst a list object of class 'time'
#'
#' @return a list of estimates for the parameters of the Arps model according to the class of 'decline_fit_lst' and 'time_lst' objects
#'
#' @export
#'
#' @examples
#' dcl_time_hyp <- decline_time(1:10000, unit = "day")
#' prod_data <- 4500 / (1 + 0.002 * 0.834 * dcl_time_hyp$t) ^ (1 / 0.834)
#' dcl_fit_param_hyp <- decline_fit_param(input_unit = "Field", output_unit = "Field",
#' fluid = "gas", model = "hyperbolic", fit_data = "rate", prod_data = prod_data,
#' initial_param = c(1000, 0.01, 1.0), lower = c(0, 1e-6, 1e-6), upper = NULL,
#' control = list(maxiter = 100))
#' dcl_fit_hyp <- decline_fit(dcl_fit_param_hyp, dcl_time_hyp)
#'
#' dcl_fit_hyp
#'
decline_fit <- function(decline_fit_lst, time_lst) {
   if ((inherits(decline_fit_lst, "decline_fit") == TRUE) & (inherits(time_lst, "time"))) {
      UseMethod("decline_fit")
   } else {
      if (!inherits(decline_fit_lst, "decline_fit")) {
         stop("A class of 'decline_fit' must be assigned to the 'decline_fit_lst' parameter of the decline_fit() function.")
      }
      if (!inherits(time_lst, "time")) {
         stop("A class of 'time' must be assigned to the 'time_lst' parameter of the decline_fit() function.")
      }
   }
}



#' S3 method for class 'decline_fit'
#'
#' Return a list of estimated parameters for the Arps exponential decline model
#'
#' @param decline_fit_lst a list object of class 'decline_fit'
#' @param time_lst a list object of class 'time'
#'
#' @return a list of estimates for the parameters of the Arps exponential model
#'
#' @export
#'
decline_fit.exponential_fit <- function(decline_fit_lst, time_lst) {
   exponential_fit(decline_fit_lst, time_lst)
}




#' S3 method for class 'decline_fit'
#'
#' Return a list of estimated parameters for the Arps harmonic decline model
#'
#' @param decline_fit_lst a list object of class 'decline_fit'
#' @param time_lst a list object of class 'time'
#'
#' @return a list of estimates for the parameters of the Arps harmonic model
#'
#' @export
#'
decline_fit.harmonic_fit <- function(decline_fit_lst, time_lst) {
   harmonic_fit(decline_fit_lst, time_lst)
}




#' S3 method for class 'decline_fit'
#'
#' Return a list of estimated parameters for the Arps hyperbolic decline model
#'
#' @param decline_fit_lst a list object of class 'decline_fit'
#' @param time_lst a list object of class 'time'
#'
#' @return a list of estimates for the parameters of the Arps hyperbolic model
#'
#' @export
#'
decline_fit.hyperbolic_fit <- function(decline_fit_lst, time_lst) {
   hyperbolic_fit(decline_fit_lst, time_lst)
}



#' S3 method for class 'decline_fit'
#'
#' Return a list of estimated parameters for the Arps modified_hyperbolic decline model
#'
#' @param decline_fit_lst a list object of class 'decline_fit'
#' @param time_lst a list object of class 'time'
#'
#' @return a list of estimates for the parameters of the Arps modified_hyperbolic model
#'
#' @export
#'
decline_fit.modified_hyperbolic_fit <- function(decline_fit_lst, time_lst) {
   modified_hyperbolic_fit(decline_fit_lst, time_lst)
}





