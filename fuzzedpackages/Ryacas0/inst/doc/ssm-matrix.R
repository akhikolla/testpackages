## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- message=FALSE-----------------------------------------------------------
library(Ryacas0)
library(Matrix)

## -----------------------------------------------------------------------------
get_output_width()
set_output_width(120)
get_output_width()

## -----------------------------------------------------------------------------
N <- 3
L1chr <- diag("1", 1 + N)
L1chr[cbind(1+(1:N), 1:N)] <- "-a"
L1s <- as.Sym(L1chr)
L1s

## -----------------------------------------------------------------------------
# FIXME: * vs %*%
K1s <- Simplify(L1s * Transpose(L1s))
V1s <- Simplify(Inverse(K1s))

## ---- results="asis"----------------------------------------------------------
cat(
  "\\begin{align} K_1 &= ", TeXForm(K1s), " \\\\ 
                  V_1 &= ", TeXForm(V1s), " \\end{align}", sep = "")

## -----------------------------------------------------------------------------
N <- 3
L2chr <- diag("1", 1 + 2*N)
L2chr[cbind(1+(1:N), 1:N)] <- "-a"
L2chr[cbind(1 + N + (1:N), 1 + 1:N)] <- "-b"
L2s <- as.Sym(L2chr)
L2s

## -----------------------------------------------------------------------------
K2s <- Simplify(L2s * Transpose(L2s))
V2s <- Simplify(Inverse(K2s))

## ---- results="asis"----------------------------------------------------------
cat(
  "\\begin{align} K_2 &= ", TeXForm(K2s), " \\\\ 
                  V_2 &= ", TeXForm(V2s), " \\end{align}", sep = "")

## -----------------------------------------------------------------------------
sparsify <- function(x) {
  Matrix::Matrix(x, sparse = TRUE)
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
V1s_eval <- Eval(V1s, list(a = alpha))
V2s_eval <- Eval(V2s, list(a = alpha, b = beta))

all.equal(V1, V1s_eval)
all.equal(V2, V2s_eval)

