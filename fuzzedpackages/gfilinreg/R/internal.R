#' @importFrom Rcpp evalCpp
#' @useDynLib gfilinreg
NULL

#' @importFrom arrangements icombinations
#' @importFrom EigenR Eigen_rank Eigen_inverse
#' @importFrom utils head
#' @importFrom stats dt qt dlogis qlogis model.matrix as.formula
#' @importFrom lazyeval f_eval_lhs f_rhs
#' @noRd
gfilinregR <- function(
  formula, data = NULL, distr = "student", df = Inf, L = 10L, lucky = FALSE
){
  distr <- match.arg(distr, c("student", "logistic"))
  if(distr == "student"){
    qdistr <- function(x) qt(x, df = df)
    ddistr <- function(x) dt(x, df = df, log = TRUE)
  }else{
    qdistr <- function(x) qlogis(x)
    ddistr <- function(x) dlogis(x, log = TRUE)
  }
  y <- f_eval_lhs(formula, data = data)
  X <- model.matrix(formula, data = data)
  betas <- colnames(X)
  X <- unname(X)
  n <- nrow(X)
  p <- ncol(X)
  if(Eigen_rank(X) < p){
    stop("Design is not of full rank.")
  }
  q <- p + 1L
  # centers of hypercubes (volume 1/L^p)
  centers <- as.matrix(
    do.call(
      expand.grid, rep(list(seq(0, 1, length.out = L+1L)[-1L] - 1/(2*L)), q)
    )
  )
  # remove centers having equal coordinates (H'H is not invertible)
  centers <-
    centers[apply(centers, 1L, function(row) length(unique(row)) > 1L),]
  # outputs
  M <- (L^q - L) / 2L # number of centers yielding sigma>0
  J <-  rep(NA_real_, M)
  Theta <- matrix(NA_real_, nrow = M, ncol = q)
  # algorithm
  Iiterator <- icombinations(n, q)
  I <- Iiterator$getnext()
  XI <- X[I, , drop = FALSE]
  while(Eigen_rank(XI) < p){
    I <- Iiterator$getnext()
    XI <- X[I, , drop = FALSE]
  }
  XmI <- X[-I, , drop = FALSE]
  yI <- y[I]
  ymI <- y[-I]
  counter <- 0L
  if(lucky){
    for(m in 1L:nrow(centers)){
      H <- cbind(XI, qdistr(centers[m, ]))
      theta <- Eigen_inverse(crossprod(H)) %*% t(H) %*% yI
      if(theta[q] > 0){ # sigma>0
        counter <- counter + 1L
        J[counter] <-
          sum(ddistr((ymI - XmI %*% head(theta, -1L))/theta[q])) -
          (n-q) * log(theta[q])
        Theta[counter,] <- theta
      }
    }
  }else{
    for(m in 1L:nrow(centers)){
      H <- cbind(XI, qdistr(centers[m, ]))
      if(Eigen_rank(H) < q){
        Theta <- head(Theta, -1L)
        J <- head(J, -1L)
        next
      }
      theta <- Eigen_inverse(crossprod(H)) %*% t(H) %*% yI
      if(theta[q] > 0){ # sigma>0
        counter <- counter + 1L
        J[counter] <-
          sum(ddistr((ymI - XmI %*% head(theta, -1L))/theta[q])) -
          (n-q) * log(theta[q])
        Theta[counter,] <- theta
      }
    }
  }
  J <- exp(J)
  out <- list(
    Theta = as.data.frame(`colnames<-`(Theta, c(betas, "sigma"))),
    weight = J/sum(J)
  )
  attr(out, "distr") <- distr
  attr(out, "df") <- df
  rhs <- as.character(f_rhs(formula))
  if(rhs[1L] == "+") rhs <- rhs[-1L]
  attr(out, "formula") <- as.formula(
    paste0("~ ", paste0(rhs, collapse = " + "))
  )
  class(out) <- "gfilinreg"
  out
}


