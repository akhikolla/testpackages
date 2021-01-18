

checkProportion <- function(proportion, paramName = "proportion", eps = 1e-10)
{
  if (missing(proportion))
    stop(paste0(paramName, " is missing"))
  if (!is.vector(proportion, mode = "numeric"))
    stop(paste0(paramName, " must be a vector of positive real whose sum equal 1"))
  if (min(proportion) < 0)
    stop(paste0(paramName, " must be a vector of positive real whose sum equal 1"))
  if (abs(1 - sum(proportion)) > eps)
    stop(paste0(paramName, " must be a vector of positive real whose sum equal 1"))
}

checkPi <- function(pi, paramName = "pi")
{
  if (missing(pi))
    stop(paste0(paramName, " is missing"))
  if (!is.numeric(pi) || !is.matrix(pi))
    stop(paste0(paramName, " must be a matrix of probabilities"))
  if ((min(pi) < 0) && (max(pi) > 1))
    stop(paste0(paramName, " must be a matrix of probabilities"))
}

checkM <- function(m)
{
  if(missing(m))
    stop("m is missing")
  if (!is.vector(m, mode = "numeric"))
    stop("m must be a (vector of) integer strictly greater than 1")
  if (length(m) != length(m[m > 1]))
    stop("m must be a (vector of) integer strictly greater than 1")
  if (!min(m == round(m)))
    stop("m must be a (vector of) integer strictly greater than 1")
}

checkM2 <- function(m, pi, mu, piName = "pi", muName = "mu")
{
  if (length(m) != ncol(pi))
    stop(paste0("The number of column of ", piName," and m do not match."))
  if (sum(m) != ncol(mu))
    stop(paste0("The number of column of ", muName," and sum(m) do not match."))
}

checkMu <- function(mu, proportion, pi, muName = "mu", proportionName = "proportion", piName = "pi")
{
  if (missing(mu))
    stop(paste0(muName, " is missing"))
  if (!is.numeric(mu) || !is.matrix(mu))
    stop(paste0(muName, " must be a matrix of positive integer"))
  if (min(mu) < 1)
    stop(paste0(muName, " must be a matrix of positive integer"))
  if (nrow(mu) != length(proportion))
    stop(paste0("The number of rows of ", muName, " and the length of ", proportionName , " do not match."))
  if (nrow(mu) != nrow(pi))
    stop(paste0("The number of rows of ", muName, " and ", piName, " do not match."))
}


checkData <- function(data)
{
  if (missing(data))
    stop("data is missing")
  if (!is.numeric(data) || !is.matrix(data))
    stop("X must be a matrix of positive integer")
  if (length(data[data >= 0]) != length(data))
    stop("data must be a matrix of positive integer")
}