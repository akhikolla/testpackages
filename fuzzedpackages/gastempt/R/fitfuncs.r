#' Functions for gastric emptying analysis
#'
#' The \code{linexp} and the power exponential (\code{powexp}) functions can
#' be used to fit gastric emptying curves.
#'
#' The \code{linexp} function can have an initial overshoot
#' to model secretion.
#'
#' \code{vol(t) = v0 * (1 + kappa * t / tempt) * exp(-t / tempt)}
#'
#' The \code{powexp} function introduced by  Elashof et al. is
#' montonously decreasing but has more freedom to model details in the
#' function tail.
#'
#' \code{vol(t) = v0 * exp(-(t / tempt) ^ beta)}
#'
#' The \code{_slope} functions return the first derivatives of \code{linexp}
#' and \code{powexp}.
#' Use the \code{_log} functions to enforce positive parameters
#' \code{tempt} and \code{beta}. Rarely required for gastric emptying curves.
#'
#' @param v0 Initial volume at t=0.
#' @param t Time after meal or start of scan, in minutes; can be a vector.
#' @param tempt Emptying time constant in minutes (scalar).
#' @param logtempt Logarithm of emptying time constant in minutes (scalar).
#' @param beta Power term for power exponential function (scalar).
#' @param logbeta Logarithm of power term for power exponential function (scalar).
#' @param kappa Overshoot term for linexp function (scalar).
#' @param logkappa Logarithm of overshoot term for linexp function (scalar).
#' @param pars Default NULL. If not NULL, the other parameters with exception
#' of \code{t} are not used and are retrieved as named parameters
#' from the numeric vector pars instead.
#' @return Vector of \code{length(t)} for computed volume.
#' @examples
#' t = seq(0,100, by=5)
#' kappa = 1.3
#' tempt = 60
#' v0 = 400
#' beta = 3
#' pars = c(v0 = v0, tempt = tempt, kappa = kappa)
#' par(mfrow=c(1,3))
#' plot(t, linexp(t, v0, tempt, kappa), type = "l", ylab = "volume",
#'    main = "linexp\nkappa = 1.3 and 1.0")
#' lines(t, linexp(t, v0, tempt, 1), type = "l", col = "green")
#' # This should give the same plot as above
#' plot(t, linexp(t, pars = pars), type = "l", ylab = "volume",
#'    main = "linexp\nkappa = 1.3 and 1.0\nwith vectored parameters")
#' lines(t, linexp(t, v0, tempt, 1), type = "l", col = "green")
#' plot(t, powexp(t, v0, tempt, beta), type = "l", ylab = "volume",
#'   main = "powexp\nbeta = 2 and 1")
#' lines(t, powexp(t, v0, tempt, 1), type = "l", col = "green")
#' @name gastemptfunc
#' @rdname gastemptfunc
#' @export
linexp = function(t, v0 = 1, tempt = NULL, kappa = NULL, pars = NULL) {
  use_pars = is.null(tempt) ||  is.null(kappa)
  if (use_pars == is.null(pars))
      stop("Either (tempt, kappa) or pars must be given in linexp")
  if (!is.null(pars)) {
    v0 = pars[["v0"]]
    tempt = pars[["tempt"]]
    kappa = pars[["kappa"]]
  }
  as.numeric(v0 * (1 + kappa * t / tempt) * exp(-t / tempt))
}

#' @rdname gastemptfunc
#' @export
linexp_slope = function(t, v0 = 1, tempt = NULL, kappa = NULL, pars = NULL){
  use_pars = is.null(tempt) ||  is.null(kappa)
  if (use_pars == is.null(pars))
    stop("Either (tempt, kappa) or pars must be given in linexp_slope")
  if (!is.null(pars)) {
    v0 = pars[["v0"]]
    tempt = pars[["tempt"]]
    kappa = pars[["kappa"]]
  }
  # d/dt v (1+(k t)/p) exp(-t/p)  Wolframalpha
  as.numeric((v0 * exp(-t/tempt)*((kappa - 1)*tempt - kappa*t))/(tempt*tempt))
}

#' @rdname gastemptfunc
#' @export
linexp_auc = function(v0 = 1, tempt = NULL, kappa = NULL, pars = NULL){
  use_pars = is.null(tempt) &&  is.null(kappa)
  if (use_pars == is.null(pars))
    stop("Either (tempt, kappa) or pars must be given in linexp_auc")
  if (!is.null(pars)) {
    v0 = pars[["v0"]]
    tempt = pars[["tempt"]]
    kappa = pars[["kappa"]]
  }
  if (is.null(v0) || is.null(tempt) || is.null(kappa))
    stop("Parameters v0, tempt and kappa must be given")
  as.numeric(v0 * (1 + kappa) * tempt)
}


#' @rdname gastemptfunc
#' @export
powexp =  function(t, v0 = 1, tempt = NULL, beta = NULL, pars = NULL){
  use_pars = is.null(tempt) &&  is.null(beta)
  if (use_pars == is.null(pars))
    stop("Either (tempt, beta) or pars must be given in powexp")
  if (!is.null(pars)) {
    v0 = pars[["v0"]]
    tempt = pars[["tempt"]]
    beta = pars[["beta"]]
  }
  if (is.null(v0) || is.null(tempt) || is.null(beta))
    stop("Parameters v0, tempt and beta must be given")
  as.numeric(v0 * exp(-(t / tempt) ^ beta))
}

#' @rdname gastemptfunc
#' @export
powexp_slope = function(t, v0 = 1, tempt = NULL, beta = NULL, pars = NULL){
  use_pars = is.null(tempt) &&  is.null(beta)
  if (use_pars == is.null(pars))
    stop("Either (tempt, beta) or pars must be given in powexp_slope")
  if (!is.null(pars)) {
    v0 = pars[["v0"]]
    tempt = pars[["tempt"]]
    beta = pars[["beta"]]
  }
  .expr1 <- t/tempt
  .expr4 <- v0 * exp(-.expr1^beta)
  -(.expr4 * (.expr1^(beta - 1) * (beta * (1/tempt))))
}


#' @rdname gastemptfunc
#' @export
linexp_log = function(t, v0 = 1, logtempt = NULL, logkappa = NULL, pars = NULL){
  use_pars = is.null(logtempt) &&  is.null(logkappa)
  if (use_pars == is.null(pars))
    stop("Either (v0, logtempt, logkappa) or pars must be given in linexp_log")
  if (!is.null(pars)) {
    v0 = as.numeric(pars["v0"])
    tempt = exp(as.numeric(pars[["logtempt"]]))
    kappa = exp(as.numeric(pars[["logkappa"]]))
  } else {
    tempt = exp(logtempt)
    kappa = exp(logkappa)
  }
  as.numeric(v0 * (1 + kappa * t / tempt) * exp(-t / tempt))
}

#' @rdname gastemptfunc
#' @export
powexp_log = function(t, v0 = 1, logtempt = NULL, logbeta = NULL, pars = NULL){
  use_pars =  is.null(logtempt) &&  is.null(logbeta)
  if (use_pars == is.null(pars))
    stop("Either (v0, tempt, beta) or pars must be given in powexp_log")
  if (!is.null(pars)) {
    v0 = as.numeric(pars["v0"])
    tempt = exp(as.numeric(pars[["logtempt"]]))
    beta = exp(as.numeric(pars[["logbeta"]]))
  } else {
    tempt = exp(logtempt)
    beta = exp(logbeta)
  }
  as.numeric(v0 * exp(-(t / tempt) ^ beta))
}

