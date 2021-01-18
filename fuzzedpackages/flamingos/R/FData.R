#' A Reference Class which represents functional data.
#'
#' FData is a reference class which represents general independent and
#' identically distributed (i.i.d.) functional objects. The data can be ordered
#' by time (functional time series). In the last case, the field `X` represents
#' the time.
#'
#' @field X Numeric vector of length \emph{m} representing the
#'   covariates/inputs.
#' @field Y Matrix of size \eqn{(n, m)} representing the observed
#'   responses/outputs. `Y` consists of \emph{n} functions of `X` observed at
#'   points \eqn{1,\dots,m}.
#' @export
FData <- setRefClass(
  "FData",
  fields = list(
    X = "numeric", # Covariates
    Y = "matrix", # Response
    m = "numeric",
    n = "numeric",
    vecY = "matrix"
  ),
  methods = list(

    initialize = function(X = numeric(1), Y = matrix(1)) {

      X <<- X
      Y <<- as.matrix(Y)

      n <<- nrow(Y)
      m <<- ncol(Y)

      vecY <<- matrix(t(Y), ncol = 1)

      if (n == 1) {
        Y <<- t(Y)
      }

    }
  )
)
