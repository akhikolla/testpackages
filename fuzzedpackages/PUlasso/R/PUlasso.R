#' PUlasso : An efficient algorithm to solve Positive and Unlabeled(PU) problem with lasso or group lasso penalty
#' @description The package efficiently solves PU problem in low or high dimensional setting using Maximization-Minorization and (block) coordinate descent. It allows simultaneous feature selection and parameter estimation for classification. Sparse calculation and parallel computing are supported for the further computational speed-up. See Hyebin Song, Garvesh Raskutti (2018) <\url{https://arxiv.org/abs/1711.08129}>.
#' @details
#' Main functions: grpPUlasso, cv.grpPUlasso, coef, predict
#' @author Hyebin Song, \email{hsong@@stat.wisc.edu}, Garvesh Raskutti, \email{raskutti@@stat.wisc.edu}.
#' @keywords PUlearning, Lasso, Group Lasso
#' @examples
#' data("simulPU")
#' fit<-grpPUlasso(X=simulPU$X,z=simulPU$z,py1=simulPU$truePY1)
#' \dontrun{
#' cvfit<-cv.grpPUlasso(X=simulPU$X,z=simulPU$z,py1=simulPU$truePY1)
#' }
#' coef(fit,lambda=fit$lambda[10])
#' predict(fit,newdata = head(simulPU$X), lambda= fit$lambda[10],type = "response")
"_PACKAGE"