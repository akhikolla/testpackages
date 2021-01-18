# Summary method for objects of class "tssdr"
`summary.tssdr` <- function(object, type = c("rectangle", "alllag", "alldir", "big"), thres = 0.8, ...) {
  type <- match.arg(type)
  L <- object$L
  # Lag-direction combinations (type = "big" only)
  if (type == "big") {
    pk <- lagdir_big(L, thres = thres)
    p <- max(pk[, 2])
  } else {
    ld <- lagdir_rect(L, thres = thres)
    switch(type,
                rectangle = {
                  p <- ld$dir_rect  # Number of directions
                  k <- ld$lag_rect  # Number of lags
                },
                alllag = {
                  p <- ld$dir_alllag
                  k <- ld$lag_alllag
                  },
                alldir = {
                  p <- ld$dir_alldir
                  k <- ld$lag_alldir
                  }
                )
  }
  W <- object$W[1:p, , drop = F] #The p'*p matrix: first p' rows of a matrix W
  S <- object$S[, 1:(p + 1)] #The response and the first p' directions
  RES <- list(W = W, L = L, S = S, type = type, algorithm = object$algorithm,
              yname = object$yname, Xname = object$Xname)
  
  if (type == "big") {
    RES$pk = pk
    } else {
    RES$p <- p
    RES$k <- k
  }
  class(RES) <- "summary.tssdr"
  RES
}

# Printing method for objects of class "summary.tssdr"
`print.summary.tssdr` <- function(x, digits = 3, ...) {
  cat("\n")
  cat(strwrap(paste("Summary of", x$algorithm, "for response",
                    x$yname, "and predictors", x$Xname), prefix = "\t"), "\n")
  cat("\n")
  cat(paste("The signal separation matrix W is:"), "\n")
  cat("\n")
  print(x$W, digits = digits)
  cat("\n")
  cat("The L matrix is: \n")
  cat("\n")
  print(x$L, digits = digits)
  cat("\n")
  if (x$type == "rectangle") {
    cat("Using the rectangle method:\n")
    cat("\n")
    if (x$p != 1) {
      cat(paste("The first", x$p), "directions and ")
    } else { #xp!=1 if
      cat("The first direction and ")
    } #xp=1 else
    if (x$k == 1) {
      cat("the first lag")
    } else { #xk=1 if
      cat(paste("the first", x$k), "lags")
    } #xk=1 else
    cat(" are relevant.\n")
  } else { #rectangle if
    if (x$type == "alldir") {
      cat(paste("Choosing first all the", x$p), "directions:\n")
      cat("\n")
      if (x$k == 1) {
        cat("The first lag is relevant. \n")
      } else {  #xk=1 if
        cat(paste("The first", x$k), "lags are relevant. \n")
      } #xk=1 else
    } else { #alldir if
      if (x$type == "alllag") {
        cat(paste("Choosing first all the", x$k), "lags:\n")
        cat("\n")
        if (x$p == 1) {
          cat("The first direction is relevant. \n")
        } else { #xp=1 if
          cat(paste("The first", x$p), "directions are relevant. \n")
        } #xp=1 else
      } else { #alllag if
    cat("Using the biggest values in L to choose the directions and lags:\n")
    cat("\n")
    cat("The relevant combinations of lags and directions are \n")
    cat("\n")
    print(x$pk, digits = digits)
      } #alllag else
    } #alldir else
  } #rectangle else
  cat("\n")
} #print.summary.tssdr

# Coef method for objects of class "summary.tssdr"
# Extracts the signal separation matrix k'*p
`coef.summary.tssdr` <- function(object, ...) object$W

# Extracting the important estimated directions:
# Components method for objects of class "summary.tssdr"
`components.summary.tssdr` <- function(x, ...) x$S[, -1]

# Plotting method for objects of class "summary.tssdr" (R's basic time series plot)
`plot.summary.tssdr` <- function(x, main = "The response and the chosen directions", ...) {
  plot.tssdr(x, main = main, ...)
}
