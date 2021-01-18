## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- message=FALSE-----------------------------------------------------------
library(Ryacas)
library(Matrix)

## -----------------------------------------------------------------------------
N <- 3
L1chr <- diag("1", 1 + N)
L1chr[cbind(1+(1:N), 1:N)] <- "-a"
L1s <- ysym(L1chr)
L1s

## -----------------------------------------------------------------------------
K1s <- L1s %*% t(L1s)
V1s <- solve(K1s)

## ---- results="asis"----------------------------------------------------------
cat(
  "\\begin{align} K_1 &= ", tex(K1s), " \\\\ 
                  V_1 &= ", tex(V1s), " \\end{align}", sep = "")

## -----------------------------------------------------------------------------
N <- 3
L2chr <- diag("1", 1 + 2*N)
L2chr[cbind(1+(1:N), 1:N)] <- "-a"
L2chr[cbind(1 + N + (1:N), 1 + 1:N)] <- "-b"
L2s <- ysym(L2chr)
L2s

## -----------------------------------------------------------------------------
K2s <- L2s %*% t(L2s)
V2s <- solve(K2s)
# Try simplify; causes timeout on CRAN Fedora, hence in try() call.
# Else, just use 
# V2s <- simplify(solve(K2s))
try(V2s <- simplify(V2s), silent = TRUE)

## ---- results="asis"----------------------------------------------------------
cat(
  "\\begin{align} K_2 &= ", tex(K2s), " \\\\ 
                  V_2 &= ", tex(V2s), " \\end{align}", sep = "")

## -----------------------------------------------------------------------------
sparsify <- function(x) {
  if (requireNamespace("Matrix", quietly = TRUE)) {
    library(Matrix)
    
    return(Matrix::Matrix(x, sparse = TRUE))
  }
  
  return(x)
}

alpha <- 0.5
beta <- -0.3

## AR(1)
N <- 3
L1 <- diag(1, 1 + N)
L1[cbind(1+(1:N), 1:N)] <- -alpha
K1 <- L1 %*% t(L1)
V1 <- solve(K1)
sparsify(K1)
sparsify(V1)

## Dynamic linear models
N <- 3
L2 <- diag(1, 1 + 2*N)
L2[cbind(1+(1:N), 1:N)] <- -alpha
L2[cbind(1 + N + (1:N), 1 + 1:N)] <- -beta
K2 <- L2 %*% t(L2)
V2 <- solve(K2)
sparsify(K2)
sparsify(V2)

## -----------------------------------------------------------------------------
V1s_eval <- eval(yac_expr(V1s), list(a = alpha))
V2s_eval <- eval(yac_expr(V2s), list(a = alpha, b = beta))
all.equal(V1, V1s_eval)
all.equal(V2, V2s_eval)

