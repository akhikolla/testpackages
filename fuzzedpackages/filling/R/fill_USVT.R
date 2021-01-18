#' Matrix Completion by Universal Singular Value Thresholding
#'
#' \code{fill.USVT} is a matrix \emph{estimation} method suitable for low-rank structure. In the context of
#' our package, we provide this method under the matrix \emph{completion} problem category. It aims at
#' exploiting the idea of thresholding the singular values to minimize the mean-squared error, defined as
#' \deqn{\mathrm{MSE}(\hat{A}):={E} \left\{ \frac{1}{np} \sum_{i=1}^{n} \sum_{j=1}^{p} (\hat{a}_{ij} - a_{ij})^2  \right\}}
#' where \eqn{A} is an \eqn{(n\times p)} matrix with some missing values and \eqn{\hat{A}} is an estimate.
#'
#' @param A an \eqn{(n\times p)} partially observed matrix.
#' @param eta control for thresholding \eqn{\in (0,1)}.
#'
#' @return a named list containing \describe{
#' \item{X}{an \eqn{(n\times p)} estimated matrix after completion, which is \eqn{\hat{A}} in the equation above.}
#' }
#'
#' @examples
#' \dontrun{
#' ## load image data of 'lena128'
#' data(lena128)
#'
#' ## transform 5% of entries into missing
#' set.seed(5)
#' A <- aux.rndmissing(lena128, x=0.05)
#'
#' ## apply the method with 3 different control 'eta'
#' fill1 <- fill.USVT(A, eta=0.01)
#' fill2 <- fill.USVT(A, eta=0.5)
#' fill3 <- fill.USVT(A, eta=0.99)
#'
#' ## visualize only the last ones from each run
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' image(A, col=gray((0:100)/100), axes=FALSE, main="5% missing")
#' image(fill1$X, col=gray((0:100)/100), axes=FALSE, main="eta=0.01")
#' image(fill2$X, col=gray((0:100)/100), axes=FALSE, main="eta=0.5")
#' image(fill3$X, col=gray((0:100)/100), axes=FALSE, main="eta=0.99")
#' par(opar)
#' }
#'
#' @references
#' \insertRef{chatterjee2015}{filling}
#'
#' @export
fill.USVT <- function(A, eta=0.01){
  ##############################################################
  # Parameter
  #   1. data check, dimension, and missing-1 matrix
  X = check_data(A)
  if (check_bycol(X)==FALSE){   warning("* fill.USVT : there exists at least one column full of missing entries.")}
  if (check_bycol(t(X))==FALSE){warning("* fill.USVT : there exists at least one row full of missing entries.")}
  # 2. eta
  if ((length(eta)>1)||(!is.numeric(eta))||(eta<=0)||(eta>=1)){
    stop("* fill.USVT : a parameter 'eta' should be in (0,1).")
  }

  ##############################################################
  # Preparation Step
  # 1. scaling : later the output must be rescaled back : bounded by abs <= 1
  a = min(X[!is.na(X)])
  b = max(X[!is.na(X)])
  # 2. size must be (mxn) with m <= n
  if (nrow(X) > ncol(X)){
    tflag = TRUE   # must be transposed later
    X     = t(X)
  } else {
    tflag = FALSE
  }
  # 3. size parameters
  m = nrow(X)
  n = ncol(X)
  # 4. rescaling of the data
  M = ((X - array((a+b)/2, c(m,n)))/((b-a)/2))

  ##############################################################
  # Main Computation as in Paper, under $1.2
  # (1) fill in the proxy Y
  Y = M
  Y[is.na(M)] = 0.0
  # (2) singular value decomposition
  svdY = svd(Y)
  s    = svdY$d
  # (3) proportion of observed values of X
  phat = (sum(!is.na(A))/length(A))
  # (4) thresholding value
  thrval = ((2+eta)*sqrt(n*phat))
  S      = (svdY$d >= phat)
  # (5) Define W
  if (sum(S)==1){
    W = (outer(svdY$u[,], svdY$v[,S])*svdY$d[S])/phat
  } else {
    W = (svdY$u[,S]%*%diag(svdY$d[S])%*%t(svdY$v[,S]))/phat
  }
  # (6) Alter Mhat
  Mhat = W
  Mhat[(W>1)]  = 1.0
  Mhat[(W<-1)] = -1.0
  # (7) re-scale back
  result = ((Mhat*((b-a)/2)) + array((a+b)/2, c(m,n)))
  # (8) control the transpose
  if (tflag==TRUE){
    result = t(result)
  }


  ##############################################################
  ## RETURN
  output = list()
  output$X = result
  return(output)
}


# ## load image data of 'lena128'
# data(lena128)
#
# # some parameters
# ngray = 64
# missp = 0.1
#
# ## transform 5% of entries into missing
# set.seed(5)
# A = lena128
# B = aux.rndmissing(lena128, x=missp)
#
# B1 = fill.simple(B, method="mean")
# B2 = fill.HardImpute(B, lambdas=c(500,100,50), rk=100)
# B3 = fill.OptSpace(B, ropt=15)
# B4 = fill.nuclear(B)
# B5 = fill.KNNimpute(B, k=5)
# B6 = fill.SVDimpute(B, k=5)
#
# graphics.off()
# x11(height = 5.5, width=9)
# par(mfrow=c(2,4), mar=c(2, 2, 2, 2))
# image(A, col=gray((0:ngray)/ngray), axes=FALSE, main="(A)")
# image(B, col=gray((0:ngray)/ngray), axes=FALSE, main="(B)")
# image(B1$X, col=gray((0:ngray)/ngray), axes=FALSE, main="(C)")
# image(B2$X[,,3], col=gray((0:ngray)/ngray), axes=FALSE, main="(D)")
# image(B3$X, col=gray((0:ngray)/ngray), axes=FALSE, main="(E)")
# image(B4$X, col=gray((0:ngray)/ngray), axes=FALSE, main="(F)")
# image(B5$X, col=gray((0:ngray)/ngray), axes=FALSE, main="(G)")
# image(B6$X, col=gray((0:ngray)/ngray), axes=FALSE, main="(H)")
# savePlot(filename="/home/kisung/Desktop/filling/paper/figure-1.png")
