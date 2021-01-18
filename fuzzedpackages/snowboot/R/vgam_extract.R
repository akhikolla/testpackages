#These function were taken from VGAM (v. 1.1-2) because the package was going to be archived.
# Thomas W. Yee (2019). VGAM: Vector Generalized Linear and Additive Models. R package version 1.1-2. URL
# https://CRAN.R-project.org/package=VGAM

VGAM_zeta = function (x, deriv = 0, shift = 1)
{
  deriv.arg <- deriv
  rm(deriv)
  if (!is.Numeric(deriv.arg, length.arg = 1, integer.valued = TRUE))
    stop("'deriv' must be a single non-negative integer")
  if (deriv.arg < 0 || deriv.arg > 2)
    stop("'deriv' must be 0, 1, or 2")
  if (deriv.arg > 0)
    return(zeta.specials(Zeta.derivative(x, deriv.arg = deriv.arg,
                                         shift = shift), x, deriv.arg, shift))
  if (any(special <- Re(x) <= 1)) {
    ans <- x
    ans[special] <- Inf
    special3 <- Re(x) < 1
    ans[special3] <- NA
    special4 <- (0 < Re(x)) & (Re(x) < 1) & (Im(x) == 0)
    ans[special4] <- Zeta.derivative(x[special4], deriv.arg = deriv.arg,
                                     shift = shift)
    special2 <- Re(x) < 0
    if (any(special2)) {
      x2 <- x[special2]
      cx <- 1 - x2
      ans[special2] <- 2^(x2) * pi^(x2 - 1) * sin(pi *
                                                    x2/2) * gamma(cx) * Recall(cx)
    }
    if (any(!special)) {
      ans[!special] <- Recall(x[!special])
    }
    return(zeta.specials(ans, x, deriv.arg, shift))
  }
  aa <- 12
  ans <- 0
  for (ii in 0:(aa - 1)) ans <- ans + 1/(shift + ii)^x
  ans <- ans + Zeta.aux(shape = x, aa, shift = shift)
  ans[shift <= 0] <- NaN
  zeta.specials(ans, x, deriv.arg = deriv.arg, shift = shift)
}


zeta.specials = function (ans, x, deriv.arg = 0, shift = 1)
{
  stieltjes <- c(
    +0.5772156649015328606065120900824024310421593359,
    -0.0728158454836767248605863758749013191377363383,
    -0.0096903631928723184845303860352125293590658061,
    +0.0020538344203033458661600465427533842857158044,
    +0.0023253700654673000574681701775260680009044694,
    +0.0007933238173010627017533348774444448307315394,
    -0.0002387693454301996098724218419080042777837151,
    -0.0005272895670577510460740975054788582819962534,
    -0.0003521233538030395096020521650012087417291805,
    -0.0000343947744180880481779146237982273906207895,
    +0.0002053328149090647946837222892370653029598537)
  names(stieltjes) <- as.character(0:10)

  ans[Im(x) == 0 & abs(x) < .Machine$double.eps & deriv.arg ==
        0 & shift == 1] <- -0.5
  ans[Im(x) == 0 & x > 1e+05 & deriv.arg == 0 & shift == 1] <- 1
  ans[Im(x) == 0 & x > 1e+05 & deriv.arg == 1 & shift == 1] <- 0
  ans[Im(x) == 0 & x > 1e+05 & deriv.arg == 2 & shift == 1] <- 0
  ans[Im(x) == 0 & abs(x) < .Machine$double.eps & deriv.arg ==
        1 & shift == 1] <- -0.5 * log(2 * pi)
  ans[Im(x) == 0 & abs(x) < .Machine$double.eps & deriv.arg ==
        2 & shift == 1] <- (-0.5 * (log(2 * pi))^2 - pi^2/24 +
                              0.5 * (stieltjes["0"])^2 + stieltjes["1"])
  ans
}


Zeta.derivative = function (x, deriv.arg = 0, shift = 1)
{
  if (!all(shift == 1))
    stop("currently 'shift' must all be 1")
  if (!is.Numeric(deriv.arg, length.arg = 1, integer.valued = TRUE))
    stop("'deriv.arg' must be a single non-negative integer")
  if (deriv.arg < 0 || deriv.arg > 2)
    stop("'deriv.arg' must be 0, 1, or 2")
  if (any(Im(x) != 0))
    stop("Sorry, currently can only handle x real, not complex")
  if (any(x < 0))
    stop("Sorry, currently cannot handle x < 0")
  ok <- is.finite(x) & x > 0 & x != 1
  ans <- rep_len(NA_real_, length(x))
  nn <- sum(ok)
  if (nn)
    ans[ok] <- .C("vzetawr", as.double(x[ok]), ans = double(nn),
                  as.integer(deriv.arg), as.integer(nn))$ans
  if (deriv.arg == 0)
    ans[is.finite(x) & abs(x) < 1e-12] <- -0.5
  ans
}


is.Numeric = function (x, length.arg = Inf, integer.valued = FALSE, positive = FALSE)
{
  if (all(is.numeric(x)) && all(is.finite(x)) && (if (is.finite(length.arg)) length(x) ==
                                                  length.arg else TRUE) && (if (integer.valued) all(x == round(x)) else TRUE) &&
      (if (positive) all(x > 0) else TRUE)) TRUE else FALSE
}


Zeta.aux = function (shape, qq, shift = 1)
{
  LLL <- max(length(shape), length(qq))
  if (length(shape) != LLL)
    shape <- rep_len(shape, LLL)
  if (length(qq) != LLL)
    qq <- rep_len(qq, LLL)
  if (any(qq < 12 - 1))
    warning("all values of argument 'q' should be 12 or more")
  aa <- qq
  B2 <- c(1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510)
  kk <- length(B2)
  ans <- 1/((shape - 1) * (shift + aa)^(shape - 1)) + 0.5/(shift +
                                                             aa)^shape
  term <- (shape/2)/(shift + aa)^(shape + 1)
  ans <- ans + term * B2[1]
  for (mm in 2:kk) {
    term <- term * (shape + 2 * mm - 3) * (shape + 2 * mm -
                                             2)/((2 * mm - 1) * 2 * mm * (shift + aa)^2)
    ans <- ans + term * B2[mm]
  }
  ifelse(aa - 1 <= qq, ans, rep(0, length(ans)))
}


