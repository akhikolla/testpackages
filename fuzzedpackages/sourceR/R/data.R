# Data classes for sourceR models

Data_ = R6::R6Class(
  'Data',
  public = list(
    x = NULL,
    initialize = function(data, ...)
    {
      private$pack(data, ...)
      private$check()
      private$sort()
    }
  ),
  private = list(
    sort = function()
    {
      dn = dimnames(self$x)
      self$x = do.call('[', c(list(self$x), lapply(dn, function(a)
        gtools::mixedorder(a)), list(drop=FALSE)))
    }
  )
)


Y_ = R6::R6Class('Y',
                inherit = Data_,
                private = list(
                  pack = function(data,
                                  y,
                                  type,
                                  time,
                                  location) {
                    if (is.null(time)) {
                      data$time = rep('1', nrow(data))
                      time = 'time'
                    }
                    if (is.null(location)) {
                      data$location = rep('1', nrow(data))
                      location = 'location'
                    }
                    data = data[, c(y, type, time, location)]
                    names(data) = c('y', 'type', 'time', 'location')
                    data$type = as.character(data$type)
                    data$time = as.character(data$time)
                    data$location = as.character(data$location)
                    self$x = tryCatch(
                      reshape2::acast(data, type ~ time ~ location, value.var = 'y',drop=F),
                      condition = function(c)
                        stop(
                          "Case data must have single observation for each type/time/location combination"
                        )
                    )
                    names(dimnames(self$x)) = names(data)[-1]
                  },
                  check = function() {
                    if (!all(isFiniteInteger(self$x)))
                      stop('Missing values in case data.  Each type/time/location must have a value.')
                    if (any(self$x < 0))
                      stop('Negative value found in data')
                  }
                ))
#' Constructs disease count data
#'
#' The Y constructor function returns an R6 Y class
#' which feeds disease count data into sourceR models.
#'
#' @param data long-format data.frame containing source data
#' @param y character string giving name of disease counts column in data
#' @param type character string giving name of type column in data
#' @param time optional column denoting times of disease count observations
#' @param location optional column denoting location of disease count observations
#'
#' @return A Y disease count data structure for use in sourceR models
#' @export
Y = function(data, y, type, time = NULL, location = NULL)
  Y_$new(data, y, type, time, location)


X_ = R6::R6Class('X',
                inherit = Data_,
                private = list(
                  pack = function(data, x, type, source, time)
                  {
                    if (is.null(time)) {
                      data$time = rep('1', nrow(data))
                      time = 'time'
                    }
                    data = data[, c(x, type, source, time)]
                    names(data) = c('x', 'type', 'source', 'time')
                    data$type = as.character(data$type)
                    data$source = as.character(data$source)
                    data$time = as.character(data$time)
                    self$x = tryCatch(
                      reshape2::acast(data, type ~ source ~ time, value.var = 'x', drop=F),
                      condition = function(c)
                        stop(
                          "Source data must have a single observation for each source/type/time combination."
                        )
                    )
                    names(dimnames(self$x)) = names(data)[-1]
                  },
                  check = function()
                  {
                    if (!all(is.finite(self$x)))
                      stop("Missing value in source data.  Each source/type/time combination must have a value.")
                    if (any(self$x < 0))
                      stop("Negative value in source data.")
                  }
                ))

#' Constructs source data
#'
#' The X constructor function returns an R6 X class
#' which feeds source data into sourceR models.
#'
#' @param data long-format data.frame containing source data
#' @param x character string giving name of source counts column in data
#' @param type character string giving name of type column in data
#' @param source character string giving name of source column in data
#' @param time optional column denoting times of source observation
#'
#' @return A X source data structure for use in sourceR models
#' @export
X = function(data, x, type, source, time = NULL)
  X_$new(data, x, type, source, time)



Prev_ = R6::R6Class('Prev',
                   inherit = Data_,
                   private = list(
                     pack = function(data, prev, source, time)
                     {
                       if (is.null(time)) {
                         data$time = rep('1', nrow(data))
                         time = 'time'
                       }
                       data = data[, c(prev, source, time)]
                       names(data) = c('prev', 'source', 'time')
                       data$source = as.character(data$source)
                       data$time = as.character(data$time)
                       self$x = tryCatch(
                         reshape2::acast(data, source ~ time, value.var = 'prev', drop=F),
                         condition = function(c)
                           stop('Prevalence must have one value per source/time')
                       )
                       names(dimnames(self$x)) = c('source', 'time')
                     },
                     check = function()
                     {
                       if (any(!is.finite(self$x)))
                         stop('Non-finite value in prevalence')
                       if (any(self$x <= 0 | self$x > 1))
                         stop('Prevalence is not a probability')
                     }
                   ))

#' Constructs prevalence data
#'
#' The Prev constructor function returns an R6 Prevalence class
#' which feeds data into sourceR models.
#'
#' @param data long-format data.frame containing prevalence data by
#'             source and time.
#' @param prev character string giving name of prevalence column in data
#' @param source character string giving name of source column in data
#' @param time optional column denoting times of prevalence observation
#'
#' @return A Prev data structure for use in sourceR models
#' @export
Prev = function(data, prev, source, time = NULL)
  Prev_$new(data, prev, source, time)

#' Alpha prior hyperparameter class
Alpha_ = R6::R6Class('Alpha',
                    inherit = Data_,
                    private = list(
                      pack = function(data, alpha, source, time, location)
                      {
                        if (is.null(time)) {
                          data$time = rep('1', nrow(data))
                        }
                        if(is.null(data$location)) {
                          data$time = rep('1', nrow(data))
                        }
                        data = data[, c(alpha, source, time, location)]
                        names(data) = c('alpha', 'source', 'time', 'location')

                        self$x = tryCatch(
                          reshape2::acast(data, source ~ time ~ location, value.var = 'alpha', drop=F),
                          condition = function(c)
                            stop('Alpha must have one value per source/time/location')
                        )
                        names(dimnames(self$x)) = c('source', 'time','location')
                      },
                      check = function()
                      {
                        if (any(!is.finite(self$x)))
                          stop('Non-finite value in alpha')
                        if (any(self$x <= 0))
                          stop('Alpha must be strictly positive')
                      }
                    ))

#' Constructs alpha prior
#'
#' The Alpha constructure function returns an R6 Alpha_
#' class which feeds sanitised prior or initialisation values for alpha into the model.
#'
#' @param data long-format data.frame containing Dirichlet prior hyperparameter, source, time, and location columns
#' @param alpha name of hyperparameter column
#' @param source name of source column
#' @param time name of optional time column
#' @param location name of optional location column
#'
#' @return An Alpha_ data structure for use in sourceR models
#' @export
Alpha = function(data, alpha, source, time = NULL, location = NULL)
  Alpha_$new(data, alpha, source, time, location)


Q_ = R6::R6Class('Q',
                 inherit = Data_,
                 public = list(
                   s = NULL,
                   theta = NULL
                 ),
                 private = list(
                   pack = function(data, q, type)
                   {
                     self$theta = unique(data[,q])
                     self$s = vapply(data[,q], function(qi) which(self$theta == qi), integer(1))
                     names(self$s) = data[,type]
                   },
                   check = function()
                   {
                     if(any(is.na(self$theta)))
                       stop("NA in q value")
                     if(any(self$theta <= 0))
                       stop("q value <= 0. q must be positive.")
                   }
                 ))

#' Constructs initial values for q
#'
#' The Q constructor returns a R6 Q_ class which feeds sanitised
#' initial values for q into the model.
#'
#' @param data long-format data.frame containing type-name and q-value for each observation.
#' @param q name of 'q' column
#' @param type name of type column
#'
#' @return A Q_ data structure for use in sourceR models
#' @export
Q = function(data, q, type)
  Q_$new(data, q, type)
