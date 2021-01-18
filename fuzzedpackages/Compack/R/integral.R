
# @title
# Integration by trapezoidal
#
# @description
# sum up the areas under trapezoids.
#
# @param X a data frame or vector
#        \itemize{
#        \item If X is a data.frame, \code{dim(X)=nobs*(p+1)},
#              nobs is the total number of observations, p is size of
#              covariates
#        \item If X is a vector, X is a sequence of observed time points
#        }
# @param T.name a character specifying name for time varaible in X
# @param basis a matrix of basis
# @param sseq a sequence used as predictor variable to generate basis
# @param from,to two scalers defining domain for integration - defalt values
#                from = 1, to = 1.
#
# @return a scaler approximating integral. If X is a vector, taking integration
#   for basis on user provided time points sequence.


ITG_trap <- function(X,
                     basis, sseq,
                     T.name = "TIME",
                     from = 0, to = 1) {

  if(is.vector(X)) {
    X <- data.frame(TIME = sort(X))
    names(X) <- T.name
  } else {
    X <- X[order(X[, T.name]), ]
  }
  p <- dim(X)[2] - 1


  date_obs <- as.vector(t(X[, T.name]))
  date_ord <- match(date_obs, sseq)
  date_len <- length(date_obs)

  if("TRUE" %in% is.na(date_ord)) stop("Incomplete basis")
  if(date_obs[1] < from || date_obs[date_len] > to) stop(paste0("Integral out of range [", from, ",", to, "]") )

  sdiff <- date_obs - c(from, date_obs[-date_len])
  step_sum <- 0
  if(p > 0) {

    X.comp <- as.matrix(X[, ! colnames(X) %in% T.name])
    I <- diag(p)

    e1 <- X.comp[1, ] %*% kronecker(I, t(basis[date_ord[1], ]))

    if(date_obs[1] != from) step_sum <- e1 * sdiff[1]

    for (l in 2:date_len) {
      e2 <- X.comp[l, ] %*% kronecker(I, t(basis[date_ord[l], ]))
      step_sum <- step_sum + (e1 + e2) * sdiff[l] / 2
      e1 <- e2
    }   # Integration of X*phi over time

    if(date_obs[date_len] != to) {
      step_sum <- step_sum + e1 * (to - date_obs[date_len])
    }

  } else {

    e1 <- basis[date_ord[1], ]

    if(date_obs[1] != from) step_sum <- e1 * sdiff[1]

    for (l in 2:date_len) {
      e2 <- basis[date_ord[l], ]
      step_sum <- step_sum + (e1 + e2) * sdiff[l] / 2
      e1 <- e2
    }

    if(date_obs[date_len] != to) {
      step_sum <- step_sum + e1 * (to - date_obs[date_len])
    }

  }


  return(step_sum)

}


# Itegration by step function
#
# Sum up the areas under rectangles.
#
# @param X a data.frame or vector
# \itemize{
# \item If X is a data.frame, dim(X)=nobs*(p+1), nobs is the total number of observations, p is size of covariates
# \item If X is a vector, X is a sequence of observed time points}
# @param T.name a character specifying name for time varaible in X
# @param basis a matrix of basis
# @param sseq a sequence used as predictor variable to generate basis
# @param from,to two scalers defining domain for integration - defalt values from = 0, to = 1.
#
# @return a scaler approximating integral. If X is a vector, taking integration for basis on user provided time points sequence.


ITG_step <- function(X,
                     basis, sseq,
                     T.name = "TIME",
                     from = 0, to = 1) {


  if(is.vector(X)) {
    X <- data.frame(TIME = sort(X))
  } else {
    X <- X[order(X[, T.name]), ]
  }
  p <- dim(X)[2] - 1


  date_obs <- as.vector(t(X[, T.name]))
  date_ord <- match(date_obs, sseq)
  #date_ord <- round(sseq, digits = 4) %in% round(date_obs, digits = 4)
  #date_ord <- sseq[date_ord]
  date_len <- length(date_obs)

  if("TRUE" %in% is.na(date_ord)) stop("Incomplete basis")
  if(date_obs[1] < from || date_obs[date_len] > to) stop(paste0("Integration out of range [", from, ",", to, "]") )


  I <- diag(p)
  sdiff <- date_obs - c(from, date_obs[-date_len])

  step_sum <- 0
  if(p > 0) {
    X.comp <- as.matrix(X[, ! colnames(X) %in% T.name])

    for(l in 1:date_len) {
      e1 <- X.comp[l, ] %*% kronecker(I, t(basis[date_ord[l], ]))
      step_sum <- step_sum + e1 * sdiff[l]
    }

    if(date_obs[date_len] != to) {
      step_sum <- step_sum + e1 * (to - date_obs[date_len])
    }
  } else {

    for (l in 1:date_len) {
      e1 <- basis[date_ord[l], ]
      step_sum <- step_sum + e1 * sdiff[l]
    }

    if(date_obs[date_len] != to) {
      step_sum <- step_sum + e1 * (to - date_obs[date_len])
    }
  }

  return(step_sum)

}




ITG <- function(X, basis, sseq, T.name = "TIME", interval = c(0,1), insert = c("FALSE", "X", "basis"), method = c("step", "trapezoidal")) {
  digits = 10
  insert <- match.arg(insert)
  method <- match.arg(method)
  col.index <- colnames(X) %in% T.name
  a <- point.interp(X[,T.name], sseq)
  X <- as.matrix(X)
  X <- switch(insert,
              "FALSE" = X,
              "X" = cbind(sseq, diag(1 - a$frac) %*% X[a$left, !col.index, drop=FALSE] +
                            diag(a$frac) %*% X[a$right, !col.index , drop=FALSE]),
              "basis" = cbind(sseq, X[a$right, !col.index, drop=FALSE])
  )

  colnames(X)[col.index] <- T.name
  X[, col.index] <- round(X[, col.index], digits = digits)
  area <- switch(method,
                 "step" = ITG_step(X, basis, sseq, T.name, from = interval[1], to = interval[2]),
                 "trapezoidal" = ITG_trap(X, basis, sseq, T.name, from = interval[1], to = interval[2])
  )
  output <- list()
  output$integral <- area
  output$X <- X
  return(output)

}



