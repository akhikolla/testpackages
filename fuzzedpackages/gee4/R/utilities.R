#' @title Extract or Get Generalized Components from a Fitted GEE-MCD/WGEE-MCD
#'   Model
#'
#' @description Extract (or get) "components" - in a generalized sense - from a
#'   fitted joint mean covariance model from an object of class "geerMod".
#'
#' @param object a fitted joint mean covariance model of class "geerMod", i.e.,
#' typically the result of geer().
#' @param name a character vector specifying the name(s) of the "component".
#'
#' When sub.num is not specified or equal to 0, possible values are:
#' \describe{
#'   \item{\code{"m"}}{a vector of number of measurement for each subject}
#'   \item{\code{"Y"}}{response vector}
#'   \item{\code{"X"}}{model matrix for mean structure}
#'   \item{\code{"Z"}}{model matrix for covariance structure (the diagonal
#'   matrix)}
#'   \item{\code{"W"}}{model matrix for covariance structure (the lower triangular
#'   matrix)}
#'   \item{\code{"H"}}{a vector of weights used in WGEE-MCD}
#'   \item{\code{"FIM"}}{Fisher Information matrix}
#'   \item{\code{"theta"}}{parameter estimates of GEE-MCD/WGEE-MCD model}
#'   \item{\code{"sd"}}{standard deviations for parameter estimates of
#'   GEE-MCD/WGEE-MCD model}
#'   \item{\code{"beta"}}{parameter estimates for mean structure model}
#'   \item{\code{"lambda"}}{parameter estimates for covariace structure (the
#'   diagonal matrix)}
#'   \item{\code{"gamma"}}{parameter estimates for covariance structure (the lower
#'   triangular matrix)}
#'   \item{\code{"quasilik"}}{quasi-likelihood}
#'   \item{\code{"QIC"}}{Quasi information criterion}
#'   \item{\code{"iter"}}{number of iterations until convergence}
#'   \item{\code{"triple"}}{(p, d, q)}
#'   \item{\code{"pij"}}{a vector of remaining probability}
#'   \item{\code{"cpij"}}{a vector of cumulative remaining probability}
#' }
#'
#' When sub.num is specified, possible values are:
#' \describe{
#'   \item{\code{"m"}}{number of measurements for subject i}
#'   \item{\code{"Y"}}{response vector for subject i}
#'   \item{\code{"X"}}{model matrix of subject i for mean structure }
#'   \item{\code{"Z"}}{model matrix of subject i for covariance structure (the
#'   diagonal matrix)}
#'   \item{\code{"W"}}{model matrix of subject i for covariance structure (the
#'   lower triangular matrix)}
#'   \item{\code{"D"}}{the estimated diagonal matrix for subject i}
#'   \item{\code{"T"}}{the estimated lower triangular matrix for subject i}
#'   \item{\code{"Sigma"}}{the estimated covariance matrix for subject i}
#'   \item{\code{"mu"}}{the estimated mean for subject i}
#'   \item{\code{"pij"}}{a vector of remaining probability for subject i}
#'   \item{\code{"cpij"}}{a vector of cumulative remaining probability for subject i}
#' }
#'
#' @param sub.num refer to i's subject
#'
#' @examples fitgee.ar1 <- geer(cd4|id|time ~ 1|1, data = aids, triple =
#'   c(6,3,3), method = 'gee-mcd', corr.struct = 'ar1', rho = 0.5, control =
#'   geerControl(trace=TRUE))
#'
#' sd  <- getGEER(fitgee.ar1, "sd")
#' QIC <- getGEER(fitgee.ar1, "QIC")
#' Di  <- getGEER(fitgee.ar1, "D", 10)
#'
#' @export
getGEER <- function(object, name, sub.num) UseMethod("getGEER")

#' @describeIn getGEER Extract or Get Generalized Components from a Fitted
#'   GEE-MCD/WGEE-MCD Model
#' @export
getGEER.geerMod <- function(object,
  name = c("m", "Y", "X", "Z", "W", "H", "D", "T", "Sigma", "mu",
           "theta", "beta", "lambda", "gamma", "alpha", "sd", "FIM",
           "quasilik", "BIC", "iter", "triple", "pij", "cpij"),
  sub.num = 0)
{
  if(missing(name)) stop("'name' must not be missing")

  stopifnot(is(object,"geerMod"))

  opt     <- object@opt
  args    <- object@args
  devcomp <- object@devcomp
  rho     <- object@rho

  if(sub.num < 0 || sub.num > length(args$m))
    stop("incorrect value for 'sub.num'")

  m = args$m
  Y = args$Y
  X = args$X
  Z = args$Z
  W = args$W
  H = args$H
  theta  = drop(opt$par)

  if (devcomp$dims["ID"] == 1) corrStruct <- "id"
  else if (devcomp$dims["CS"] == 1) corrStruct <- "cs"
  else if (devcomp$dims["AR1"] == 1) corrStruct <- "ar1"

  # obj <- new("gee_jmcm", m, Y, X, Z, W, corrStruct, rho)
  obj <- .Call("gee_jmcm__new", m, Y, X, Z, W, corrStruct, rho)

  if(sub.num == 0) {
    switch(name,
      "m" = args$m,
      "Y" = args$Y,
      "X" = args$X,
      "Z" = args$Z,
      "W" = args$W,
      "H" = args$H,
      "theta"  = drop(opt$par),
      "beta"   = drop(opt$beta),
      "lambda" = drop(opt$lambda),
      "gamma"  = drop(opt$gamma),
      "alpha"  = drop(opt$alpha),
      "sd"     = .Call("gee_jmcm__get_sd", obj, theta), # "sd" = obj$get_sd(theta),
      "FIM"    = .Call("gee_jmcm__get_fim", obj, theta), # "FIM" = obj$get_fim(theta),
      "quasilik" = opt$quasilik,
      "QIC"      = opt$QIC,
      "iter"   = opt$iter,
      "triple" = object$triple,
      "pij"    = drop(opt$pij),
      "cpij"    = drop(opt$cpij))
  } else {
      if (sub.num == 1) vindex = 1
      else vindex = sum(m[1:(sub.num-1)]) + 1
      switch(name,
             "m"     = .Call("gee_jmcm__get_m",     obj, sub.num),
             "Y"     = .Call("gee_jmcm__get_Y",     obj, sub.num),
             "X"     = .Call("gee_jmcm__get_X",     obj, sub.num),
             "Z"     = .Call("gee_jmcm__get_Z",     obj, sub.num),
             "W"     = .Call("gee_jmcm__get_W",     obj, sub.num),
             "D"     = .Call("gee_jmcm__get_D",     obj, theta, sub.num),
             "T"     = .Call("gee_jmcm__get_T",     obj, theta, sub.num),
             "Sigma" = .Call("gee_jmcm__get_Sigma", obj, theta, sub.num),
             "mu"    = .Call("gee_jmcm__get_mu",    obj, theta, sub.num),
      ## "m" = obj$get_m(sub.num),
      ## "Y" = obj$get_Y(sub.num),
      ## "X" = obj$get_X(sub.num),
      ## "Z" = obj$get_Z(sub.num),
      ## "W" = obj$get_W(sub.num),
      ## "D" = obj$get_D(theta, sub.num),
      ## "T" = obj$get_T(theta, sub.num),
      ## "Sigma" = obj$get_Sigma(theta, sub.num),
      ## "mu"    = obj$get_mu(theta, sub.num),
      "pij"   = drop(opt$pij)[vindex:(vindex+m[sub.num]-1)],
      "cpij"  = drop(opt$cpij)[vindex:(vindex+m[sub.num]-1)])
  }
}

lagseq <- function(time)
{
  res <- NULL
  if(length(time) != 1) {
    for(i in 2:length(time)) {
      for(j in 1:(i-1))
        res <- c(res, (time[i] - time[j]))
    }
  }
  res
}


#' @title Plot Fitted Curves for One or More geerMod Objects
#'
#' @description Plot fitted curves and corresponding 95\% confidence interval
#' for one or more geerMod objects
#'
#' @param object a fitted joint GEE-MCD/WGEE-MCD model of class "geerMod", i.e.,
#'   typically the result of geer().
#' @param text some corresponding descriptions for the objects.
#' @param ... additional pairs of 'object' and 'text'
#' @param include.CI whether or not 95% confidence interval of the fitted values
#'   for the first object should be plotted.
#'
#' @examples fitgee.ar1 <- geer(cd4|id|time ~ 1|1, data = aids, triple =
#'   c(6,3,3), method = 'gee-mcd', corr.struct = 'ar1', rho = 0.5)
#' fittedcurve(fitgee.ar1, text = "GEE-MCD fitted curve", include.CI = TRUE)
#'
#' @export
fittedcurve <- function(object, text = "fitted curve", ..., include.CI = FALSE)
{
  dots <- list(...)
  ndots <- length(dots)

  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

  # create a gee_jmcm object
  opt     <- object@opt
  args    <- object@args
  devcomp <- object@devcomp
  rho     <- object@rho

  m = args$m
  Y = args$Y
  X = args$X
  Z = args$Z
  W = args$W
  theta  = drop(opt$par)

  if (devcomp$dims["ID"] == 1) corrStruct <- "id"
  else if (devcomp$dims["CS"] == 1) corrStruct <- "cs"
  else if (devcomp$dims["AR1"] == 1) corrStruct <- "ar1"

  ## obj <- new("gee_jmcm", m, Y, X, Z, W, corrStruct, rho)
  ## sd  <- obj$get_sd(theta)
  ## fim <- obj$get_fim(theta)
  obj <- .Call("gee_jmcm__new", m, Y, X, Z, W, corrStruct, rho)
  sd  <- .Call("gee_jmcm__get_sd", obj, theta)
  fim <- .Call("gee_jmcm__get_fim", obj, theta)

  # initialization
  opt <- object@opt

  beta <- opt$beta
  lambda <- opt$lambda
  gamma  <- opt$gamma

  lbta <- length(beta)
  llmd <- length(lambda)
  lgma <- length(gamma)

  lmd_sd <- sd[(lbta+1):(lbta+llmd)]
  gma_sd <- sd[(lbta+llmd+1):(lbta+llmd+lgma)]

  args   <- object@args
  Y <- args[["Y"]]
  time <- args[["time"]]

  ts   <- seq(min(time), max(time), length.out = 100)
  tslag <- seq(0, max(time) - min(time), length.out = 100)

  ###############################
  # plot the fitted mean curves #
  ###############################
  nline <- 1

  X.ts    <- NULL
  for(i in 0:(lbta-1)) X.ts    <- cbind(X.ts, ts^i)

  Yest <- drop(X.ts %*% beta)

  # ylim.max <- max(Yest)
  # ylim.min <- min(Yest)
  # ylim.diff <- ylim.max - ylim.min
  # plot(ts, Yest, type = "l",
  #      ylim = c(ylim.min-0.5*ylim.diff, ylim.max+0.5*ylim.diff),
  #      xlab = "Time", ylab = "Response")
  plot(ts, Yest, type = "l", ylim = c(min(Y), max(Y)), xlab = "Time", ylab = "Response")
  # plot(time, Y, xlab = "Time", ylab = "Response")
  # lines(ts, Yest)

  if (include.CI) {
    nline <- nline + 1
    fim_bta <- fim[1:lbta,1:lbta]
    Y_var <- diag(X.ts %*% tcrossprod(solve(fim_bta),X.ts))
    Yest.u <- Yest + 1.96 * sqrt(Y_var)
    Yest.l <- Yest - 1.96 * sqrt(Y_var)
    lines(ts, Yest.u, lty = 2)
    lines(ts, Yest.l, lty = 2)
    text <- c(text, "95% confidence intervals")
  }

  if (ndots)
  for(i in 1:(ndots/2)) {
    nline <- nline + 1

    object2 <- dots[[2*i-1]]
    text2   <- dots[[2*i]]
    text    <- c(text, text2)

    opt2 <- object2@opt
    beta2 <- opt2$beta
    lbta2 <- length(beta2)

    Yest2 <- drop(X.ts %*% beta2)

    col <- 'black'
    if(nline > 6) col <- 'grey'
    lines(ts, Yest2, lty = nline, col = col)
  }

  col <- rep('black', nline)
  if (nline > 6) {
    col <- c(rep('black', 6), rep('grey', nline - 6))
  }
  ncol <- 1
  if (nline > 5) ncol <- 2
  legend(min(ts), max(Y), text, lty=1:nline, col=col, ncol=ncol)

  Z.ts    <- NULL
  W.tslag <- NULL

  for(i in 0:(llmd-1)) Z.ts    <- cbind(Z.ts, ts^i)
  for(i in 0:(lgma-1)) W.tslag <- cbind(W.tslag, tslag^i)

  Zlmd <- Z.ts %*% lambda
  Wgma <- W.tslag %*% gamma

  ####################################
  # plot the fitted curves for logIV #
  ####################################
  nline <- 1

  ylim.max <- max(Zlmd)
  ylim.min <- min(Zlmd)
  ylim.diff <- ylim.max - ylim.min
  plot(ts, Zlmd, type = "l",
       ylim = c(ylim.min-0.5*ylim.diff, ylim.max+0.5*ylim.diff),
       xlab="Time", ylab="Log-innovat. var.")

  if (include.CI) {
    nline <- nline + 1
    fim_lmd <- fim[(lbta+1):(lbta+llmd), (lbta+1):(lbta+llmd)]
    Zlmd_var <- diag(Z.ts %*% tcrossprod(solve(fim_lmd),Z.ts))
    Zlmd.u <- Zlmd + 1.96 * sqrt(Zlmd_var)
    Zlmd.l <- Zlmd - 1.96 * sqrt(Zlmd_var)
    lines(ts, Zlmd.u, lty = 2)
    lines(ts, Zlmd.l, lty = 2)
  }

  if(ndots)
  for(i in 1:(ndots/2)) {
    nline <- nline + 1
    object2 <- dots[[2*i-1]]
    text2   <- dots[[2*i]]

    opt2 <- object2@opt
    lambda2 <- opt2$lambda
    llmd2 <- length(lambda2)

    Zlmd2 <- Z.ts %*% lambda2

    col <- 'black'
    if(nline > 6) col <- 'grey'
    lines(ts, Zlmd2, lty = nline, col = col)
  }

  ###################################
  # plot the fitted curves for GARP #
  ###################################
  nline <- 1

  ylim.max <- max(Wgma)
  ylim.min <- min(Wgma)
  ylim.diff <- ylim.max - ylim.min
  plot(tslag, Wgma, type = "l",
       ylim = c(ylim.min-0.5*ylim.diff, ylim.max+0.5*ylim.diff),
       xlab="Lag", ylab="Autoregres. coeffic.")
  # phi <- -Tt[upper.tri(Tt, diag=FALSE)]
  # plot(tlag, phi, xlab="Lag", ylab="Autoregres. coeffic.")
  # lines(tslag, Wgma)

  if (include.CI) {
    nline <- nline + 1
    fim_gma <- fim[(lbta+llmd+1):(lbta+llmd+lgma), (lbta+llmd+1):(lbta+llmd+lgma)]
    Wgma_var <- diag(W.tslag %*% tcrossprod(solve(fim_gma),W.tslag))
    Wgma.u <- Wgma + 1.96 * sqrt(Wgma_var)
    Wgma.l <- Wgma - 1.96 * sqrt(Wgma_var)
    lines(tslag, Wgma.u, lty = 2)
    lines(tslag, Wgma.l, lty = 2)
  }

  if(ndots)
  for(i in 1:(ndots/2)) {
    nline <- nline + 1
    object2 <- dots[[2*i-1]]
    text2   <- dots[[2*i]]

    opt2 <- object2@opt
    gamma2 <- opt2$gamma
    lgma2 <- length(gamma2)

    Wgma2 <- W.tslag %*% gamma2

    col <- 'black'
    if(nline > 6) col <- 'grey'
    lines(tslag, Wgma2, lty = nline, col = col)
  }
}

## #' @title Plot Fitted Mean Curves
## #'
## #' @description plot fitted mean curves
## #'
## #' @param object a fitted joint mean covariance model of class "geerMod", i.e.,
## #' typically the result of geer().
## #'
## #' @examples
## #' cattleA <- cattle[cattle$group=='A', ]
## #' fit.mcd <- geer(weight | id | I(ceiling(day/14 + 1)) ~ 1 | 1, data=cattleA,
## #'   triple = c(8, 3, 4), cov.method = 'mcd')
## #' meanplot(fit.mcd)
## #'
## #' @export
## meanplot <- function(object)
## {
##   op <- par(mfrow = c(1, 1))

##   opt <- object@opt
##   beta <- opt$beta
##   lbta <- length(beta)

##   args   <- object@args
##   Y <- args[["Y"]]
##   time <- args[["time"]]

##   ts   <- seq(min(time), max(time), length.out = 100)

##   X.ts    <- NULL
##   for(i in 0:(lbta-1)) X.ts    <- cbind(X.ts, ts^i)

##   Yest <- drop(X.ts %*% beta)
##   plot(time, Y, xlab = "Time", ylab = "Response")
##   lines(ts, Yest)
## }

## #' @title Plot Sample Regressograms and Fitted Curves
## #'
## #' @description Plot the sample regressograms based on the sample covariance
## #' matrix and superimpose the corresponding fitted curves to check the model
## #' fitting when the longitudinal dataset is balanced.
## #'
## #' @param object a fitted joint mean covariance model of class "geerMod", i.e.,
## #' typically the result of geer().
## #' @param time a vector of obeservation time points
## #'
## #' @examples
## #' cattleA <- cattle[cattle$group=='A', ]
## #' fit.mcd <- geer(weight | id | I(ceiling(day/14 + 1)) ~ 1 | 1, data=cattleA,
## #'   triple = c(8, 3, 4), cov.method = 'mcd')
## #' regressogram(fit.mcd, time = 1:11)
## #'
## #' @export
## regressogram <- function(object, time)
## {
##   debug <- 0

##   op <- par(mfrow = c(1, 2))

##   opt <- object@opt

##   lambda <- opt$lambda
##   gamma  <- opt$gamma

##   llmd <- length(lambda)
##   lgma <- length(gamma)

##   args   <- object@args
##   dims   <- object@devcomp$dims

##   m <- args[["m"]]
##   Y <- args[["Y"]]
##   X <- args[["X"]]
##   Z <- args[["Z"]]
##   W <- args[["W"]]

##   if (length(unique(m)) != 1)
##     stop("No regressograms. Unbalanced longitudinal dataset.")

##   # create a data matrix
##   DataMat <- t(Y[1:m[1]])
##   for(i in 2:length(m))
##   {
##     DataMat <- rbind(DataMat, t(Y[(sum(m[1:(i-1)])+1):sum(m[1:i])]))
##   }
##   dimnames(DataMat) <- NULL

##   S <- cov(DataMat)  # sample covariance matrix
##   R <- cor(DataMat)  # sample correlation matrix

##   # FIXME: singularity check
##   C <- t(chol(S))    # Cholesky factor of S
##   D <- diag(diag(C))

##   # transpose of matrix T in MCD
##   Tt <- t(forwardsolve(C %*% diag(diag(C)^(-1)), diag(dim(D)[1])))

##   # transpose of matrix L in ACD
##   Lt <- t(diag(diag(C)^(-1)) %*% C)

##   ts    <- seq(min(time), max(time), length.out = 100)
##   tlag  <- lagseq(time)
##   tslag <- seq(min(tlag), max(tlag), length.out = 100)

##   Z.ts    <- NULL
##   W.tslag <- NULL
##   for(i in 0:(llmd-1)) Z.ts <- cbind(Z.ts, ts^i)
##   for(i in 0:(lgma-1)) W.tslag <- cbind(W.tslag, tslag^i)

##   Zlmd <- Z.ts %*% lambda
##   Wgma <- W.tslag %*% gamma

##   # dims["MCD"] = 1
##   # dims[]

##   # plot regressogram for MCD, ACD or HPC
##   if (dims["MCD"] == 1) {
##     # the first plot
##     plot(time, log(diag(D)^2), xlab="Time", ylab="Log-innovat. var.")
##     lines(ts, Zlmd)

##     # the second plot
##     phi <- -Tt[upper.tri(Tt, diag=FALSE)]
##     plot(tlag, phi, xlab="Lag", ylab="Autoregres. coeffic.")
##     lines(tslag, Wgma)

##   } else if (dims["ACD"] == 1) {
##     # the first plot
##     plot(time, log(diag(D)^2), xlab="Time", ylab="Log-innovat. var.")
##     lines(ts, Zlmd)

##     # the second plot
##     phi <- Lt[upper.tri(Lt, diag=FALSE)]
##     plot(tlag, phi, xlab="Lag", ylab="MA. coeffic.")
##     lines(tslag, Wgma)

##   } else if (dims["HPC"] == 1) {
##     # the first plot
##     H <- diag(sqrt(diag(S)))
##     plot(time, log(diag(H)^2), xlab="Time", ylab="Log-variance")
##     lines(ts, Zlmd)

##     # the second plot
##     B <- t(chol(R))
##     PhiMat <- matrix(0, dim(B)[1], dim(B)[2])
##     for(j in 2:dim(B)[1]) {
##       for(k in 1:(j-1)) {
##         tmp <- 1
##         if (k != 1) {
##           tmp <- prod(sin(PhiMat[j, 1:(k-1)]))
##         } # if
##         PhiMat[j,k] <- acos(B[j, k]/tmp)
##       } # for k
##     } # for j
##     PhiMatt <- t(PhiMat)

##     phi <- PhiMatt[upper.tri(PhiMatt, diag=FALSE)]
##     plot(tlag, phi, xlab="Lag", ylab="Angles")
##     lines(tslag, Wgma)
##   } # HPC
## }

#' @title Plot Fitted Curves and Corresponding Confidence Interval using
#' bootstrapping method
#'
#' @description Plot fitted curves and corresponding 95\% confidence interval
#' using bootstrapping method.
#'
#' @param object a fitted joint mean covariance model of class "geerMod", i.e.,
#' typically the result of geer().
#' @param nboot number of the bootstrap replications.
#'
#' @examples \dontrun{
#' # It may take hours for large bootstrap replications
#' fitgee.ar1 <- geer(cd4|id|time ~ 1|1, data = aids, triple = c(6,3,3),
#'   method = 'gee-mcd', corr.struct = 'ar1', rho = 0.5,
#'   control = geerControl(trace=TRUE))
#' bootcurve(fitgee.ar1, nboot = 1000)
#' }
#'
#' @export
bootcurve <- function(object, nboot)
{
  debug <- 0

  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

  opt <- object@opt

  theta  <- opt$par
  beta   <- opt$beta
  lambda <- opt$lambda
  gamma  <- opt$gamma

  ltht   <- length(theta)
  lbta   <- length(beta)
  llmd   <- length(lambda)
  lgma   <- length(gamma)

  args   <- object@args
  rho    <- object@rho
  dims   <- object@devcomp$dims

  m <- args[["m"]]
  Y <- args[["Y"]]
  X <- args[["X"]]
  Z <- args[["Z"]]
  W <- args[["W"]]
  time <- args[["time"]]

  ts    <- seq(min(time), max(time), length.out = 100)
  tslag <- seq(0, max(time) - min(time), length.out = 100)

  X.ts    <- NULL
  Z.ts    <- NULL
  W.tslag <- NULL

  for(i in 0:(lbta-1)) X.ts    <- cbind(X.ts, ts^i)
  for(i in 0:(llmd-1)) Z.ts    <- cbind(Z.ts, ts^i)
  for(i in 0:(lgma-1)) W.tslag <- cbind(W.tslag, tslag^i)

  Yest <- drop(X.ts %*% beta)
  Zlmd <- drop(Z.ts %*% lambda)
  Wgma <- drop(W.tslag %*% gamma)

  Yest.boot <- NULL
  Zlmd.boot <- NULL
  Wgma.boot <- NULL

  result <- NULL
  for(iter in 1:nboot) {
    # generate a bootstrap sample
    index <- sample(length(m), replace=T)

    # construct corresponding arguments
    m.boot <- m[index]
    Y.boot <- NULL
    X.boot <- NULL
    Z.boot <- NULL
    W.boot <- NULL
    for(i in 1:length(m)) {
      if (index[i] == 1) {
        Y.boot <- c(Y.boot, Y[1:m[1]])
        X.boot <- rbind(X.boot, X[1:m[1], ])
        Z.boot <- rbind(Z.boot, Z[1:m[1], ])

        if (m[1] != 1) {
          first <- 1
          last  <- m[1] * (m[1] - 1) / 2
          W.boot <- rbind(W.boot, W[first:last, ])
        }
      } else {
        first <- sum(m[1:(index[i]-1)]) + 1
        last  <- sum(m[1:index[i]])
        Y.boot <- c(Y.boot, Y[first:last])
        X.boot <- rbind(X.boot, X[first:last, ])
        Z.boot <- rbind(Z.boot, Z[first:last, ])

        if (m[index[i]] != 1) {
          first <- 0
          for(j in 1:(index[i]-1)) {
            first <- first + m[j] * (m[j] - 1) / 2
          }
          last  <- first  + m[index[i]] * (m[index[i]] - 1) / 2
          first <- first + 1

          W.boot <- rbind(W.boot, W[first:last, ])
        }
      }
    } # for

    dims["MCD"] = 1
    dims["ACD"] = 0
    dims["HPC"] = 0

    if (dims["ID"] == 1) corr.struct <- "id"
    if (dims["CS"] == 1) corr.struct <- "cs"
    if (dims["AR1"] == 1) corr.struct <- "ar1"

    control <- geerControl()

    opt <- optimizeGeer(m.boot, Y.boot, X.boot, Z.boot, W.boot, time, method = 'gee-mcd',
      corr.struct = corr.struct, rho = rho, ipw.order=1, control = control, start = theta)

    result <- rbind(result, drop(opt$par))
    cat("iter ", iter, ": ", format(round(result[iter, ], 4), nsmall=4), "\n")

    beta.boot   <- drop(opt$par)[1:lbta]
    lambda.boot <- drop(opt$par)[(lbta + 1):(lbta + llmd)]
    gamma.boot  <- drop(opt$par)[(lbta + llmd + 1):(lbta + llmd + lgma)]

    Yest.boot <- rbind(Yest.boot, drop(X.ts %*% beta.boot))
    Zlmd.boot <- rbind(Zlmd.boot, drop(Z.ts %*% lambda.boot))
    Wgma.boot <- rbind(Wgma.boot, drop(W.tslag %*% gamma.boot))
  }

  Yest.boot <- apply(Yest.boot, 2, function(x) sort(x))
  Zlmd.boot <- apply(Zlmd.boot, 2, function(x) sort(x))
  Wgma.boot <- apply(Wgma.boot, 2, function(x) sort(x))

  Yest.u <- drop(Yest.boot[floor(0.975 * nboot), ])
  Yest.l <- drop(Yest.boot[ceiling(0.025 * nboot), ])

  Zlmd.u <- drop(Zlmd.boot[floor(0.975 * nboot), ])
  Zlmd.l <- drop(Zlmd.boot[ceiling(0.025 * nboot), ])
  Wgma.u <- drop(Wgma.boot[floor(0.975 * nboot), ])
  Wgma.l <- drop(Wgma.boot[ceiling(0.025 * nboot), ])

  plot(time, Y, xlab = "Time", ylab = "Response")
  lines(ts, Yest)
  lines(ts, Yest.u, lty = 2, lwd = 2)
  lines(ts, Yest.l, lty = 2, lwd = 2)

  if (dims["MCD"] == 1 || dims["ACD"] == 1) {
    xlab="Time"
    ylab="Log-innovat. var."
  }
  if (dims["HPC"] == 1) {
    xlab="Time"
    ylab="Log-variance"
  }
  plot(ts, Zlmd, type = 'l', xlab = xlab, ylab = ylab)
  lines(ts, Zlmd.u, lty = 2, lwd = 2)
  lines(ts, Zlmd.l, lty = 2, lwd = 2)

  if (dims["MCD"] == 1) {
    xlab="Lag"
    ylab="Autoregres. coeffic."
  }
  if (dims["ACD"] == 1) {
    xlab="Lag"
    ylab="MA. coeffic."
  }
  if (dims["HPC"] == 1) {
    xlab="Lag"
    ylab="Angles"
  }
  plot(tslag, Wgma, type = 'l', xlab = xlab, ylab = ylab)
  lines(tslag, Wgma.u, lty = 2, lwd = 2)
  lines(tslag, Wgma.l, lty = 2, lwd = 2)

}

# #' @export
# globalSearch <- function(formula, data = NULL,
#   cov.method = c('mcd', 'acd', 'hpc'),
#   control = geerControl())
# {
#   args <- ldFormula(formula, data, triple = c(1, 1, 1),
#                     cov.method, control)
#
#   m <- max(args$m)
#
#   triple <- c(m-1, m-1, m-1)
#   full <- ans <- geer(formula, data, triple=triple, cov.method, control)
#
#   bta0 <- ans@opt$beta
#   lmd0 <- ans@opt$lambda
#   gma0 <- ans@opt$gamma
#
#   cat("-------------------------------------------------------\n")
#
#   for(i in (m-2):1) {
#     triple <- c(i, m-1, m-1)
#     fit <- geer(formula, data, triple=triple, cov.method, control)
#     if(ans@opt$BIC > fit@opt$BIC) {
#       p   <- i
#       ans <- fit
#       cat("triple: ")    ; print(triple)
#       cat("    logLik: "); print(fit@opt$loglik)
#       cat("    BIC   : "); print(fit@opt$BIC)
#     }
#   }
#
#   cat("-------------------------------------------------------\n")
#
#   {
#     ans <- full
#     for(i in (m-2):1) {
#       triple <- c(m-1, i, m-1)
#       fit <- geer(formula, data, triple=triple, cov.method, control)
#       if(ans@opt$BIC > fit@opt$BIC) {
#         d <- i
#         ans <- fit
#         cat("triple: "); print(triple)
#         cat("    logLik: "); print(fit@opt$loglik)
#         cat("    BIC   : "); print(fit@opt$BIC)
#       }
#     }
#
#     cat("-------------------------------------------------------\n")
#
#     ans <- full
#     for(i in (m-2):1) {
#       triple <- c(m-1, m-1, i)
#       fit <- geer(formula, data, triple=triple, cov.method, control)
#       if(ans@opt$BIC > fit@opt$BIC) {
#         q <- i
#         ans <- fit
#         cat("triple: "); print(triple)
#         cat("    logLik: "); print(fit@opt$loglik)
#         cat("    BIC   : "); print(fit@opt$BIC)
#       }
#     }
#   }
#
#   cat("-------------------------------------------------------\n")
#
#   cat("p = "); print(p)
#   cat("d = "); print(d)
#   cat("q = "); print(q)
#
#   c(p,d,q)
# }

#' @title Generate Cattle B Data with Pre-specified Dropout Rate
#'
#' @description We assume that all cattle have the first four observations and
#'   there is a certain chance that some cows will quit the study (dropout) at
#'   the fifth measurement if their weights are below a certain threshold. The
#'   threshold value for the weight is chosen so that some fixed rates of MAR
#'   dropout at the fifth measurement time.
#'
#' @param dropout.rate the dropout rate at the fifth measurement.
#'
#' @export
GenerateCattleMAR <- function(dropout.rate)
{
  # data("cattle")
  cattleB <- cattle[cattle$group == 'B', ]
  # utils::data(sysdata, envir = environment())
  # cattleB <- sysdata[sysdata$group == 'B', ]
  # data matrix for cattleB (without dropout)
  dm <- NULL
  for(i in 31:60) dm <- rbind(dm, cattleB[cattleB$id == i, ]$weight)
  # data matrix for cattleB (ordered by the 5th measurement )
  datamatrix <- dm[order(dm[,5]), ]

  time1 <- c(0, 14, 28, 42)
  time2 <- c(0, 14, 28, 42, 56, 70, 84, 98, 112, 126, 133)

  ndrop  <- 30 * dropout.rate
  id     <- c(rep(1:ndrop, each=4), rep((ndrop+1):30, each=11))

  day1 <- rep(time1, ndrop)
  day2 <- rep(time2, 30-ndrop)
  day <- c(day1, day2)

  group  <- rep('B', ndrop*4 + (30-ndrop)*11)

  weight1 <- t(datamatrix[1:ndrop, 1:4])
  weight2 <- t(datamatrix[(ndrop+1):30, ])
  weight <- c(c(weight1), c(weight2))

  data.frame(id, day, group, weight)
}


## # function to generate cattleB data with specified parameters for dropout model
## #' @export
# GenerateCattleIPW <- function(alpha)
# {
#   debug = 0
#
#   cattleB <- cattle[cattle$group == 'B', ]
#   # data matrix for cattleB (without dropout)
#   datamatrix <- NULL
#   for(i in 31:60)
#     datamatrix <- rbind(datamatrix, cattleB[cattleB$id == i, ]$weight)
#
#   time <- c(0, 14, 28, 42, 56, 70, 84, 98, 112, 126, 133)
#
#   m <- rep(11, 30)
#   Y <- cattleB$weight
#   ipwobj <- new("ipw", m, Y, 1)
#   p <- ipwobj$get_p(alpha)
#
#   id <- NULL
#   day <- NULL
#   group <- NULL
#   weight <- NULL
#   for(sub.num in 1:30) {
#
#     if (sub.num == 1) vindex = 1
#     else vindex = sum(m[1:(sub.num-1)]) + 1
#     pij <- p[vindex:(vindex+m[sub.num]-1)]
#
#     if(debug) cat("\npij =", pij)
#
#     obsindex <- 1
#     remain   <- 1
#     while(obsindex <= 11) {
#       remain <- rbinom(1,1,pij[obsindex])
#       if (remain == 0) break
#
#       id <- c(id, sub.num)
#       day <- c(day, time[obsindex])
#       group <- c(group, 'B')
#       weight <- c(weight, datamatrix[sub.num, obsindex])
#
#       obsindex <- obsindex + 1
#
#       if(debug) {
#         if(is.na(remain)) cat("i = ", sub.num, " j = ", obsindex)
#       }
#     }
#   }
#
#   data.frame(id, day, group, weight)
# }


#' @title Calculate the Remaining Probabilities in Inverse Probability Weights
#'
#' @description Calculate a vector of remaining probabilities in inverse
#'   probability weights for one specific or all subjects
#'
#' @param m an integer vector of number of measurements for each subject.
#' @param Y a vector of responses for all subjects.
#' @param order the order for MAR remaining model.
#' @param alpha parameters in MAR remaining model.
#' @param sub.num subject number (0 = all subjects).
#'
#' @export
CalculateIPWprob <- function(m, Y, order, alpha, sub.num = 0) {
  ## obj <- new("ipw", m, Y, order)
  ## pij <- obj$get_p(alpha)
  obj <- .Call("ipw__new", m, Y, order)
  pij <- .Call("ipw__get_p", obj, alpha)

  if (sub.num != 0) {
    if (sub.num == 1) vindex = 1
    else vindex = sum(m[1:(sub.num-1)]) + 1
    pij <- pij[vindex:(vindex+m[sub.num]-1)]
  }

  pij
}


#' @title Calculate the Cumulative Remaining Probabilities in Inverse
#'   Probability Weights
#'
#' @description Calculate a vector of cumulative remaining probabilities in
#'   inverse probability weights for one specific or all subjects.
#'
#' @param m an integer vector of number of measurements for each subject.
#' @param Y a vector of responses for all subjects.
#' @param order the order for MAR remaining model.
#' @param alpha parameters in MAR remaining model.
#' @param sub.num subject number (0 = all subjects).
#'
#' @export
CalculateIPWcumprob <- function(m, Y, order, alpha, sub.num = 0) {
  ## obj <- new("ipw", m, Y, order)
  ## Pi <- obj$get_Pi(alpha)
  obj <- .Call("ipw__new", m, Y, order)
  Pi <- .Call("ipw__get_Pi", obj, alpha)

  if (sub.num != 0) {
    if (sub.num == 1) vindex = 1
    else vindex = sum(m[1:(sub.num-1)]) + 1
    Pi <- Pi[vindex:(vindex+m[sub.num]-1)]
  }

  Pi
}
