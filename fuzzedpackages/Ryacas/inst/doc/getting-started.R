## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- message=FALSE-----------------------------------------------------------
library(Ryacas)

## -----------------------------------------------------------------------------
eq <- "x^2 + 4 + 2*x + 2*x"
yac_str(eq) # No task was given to yacas, so we simply get the same returned
yac_str(paste0("Simplify(", eq, ")"))
yac_str(paste0("Factor(", eq, ")"))
yac_expr(paste0("Factor(", eq, ")"))
yac_str(paste0("TeXForm(Factor(", eq, "))"))

## -----------------------------------------------------------------------------
y_fn(eq, "Factor")
yac_str(y_fn(eq, "Factor"))
yac_str(y_fn(y_fn(eq, "Factor"), "TeXForm"))

## -----------------------------------------------------------------------------
eq %>% y_fn("Factor")
eq %>% y_fn("Factor") %>% yac_str()
eq %>% y_fn("Factor") %>% y_fn("TeXForm") %>% yac_str()

## -----------------------------------------------------------------------------
yac_str(paste0("Factor(", eq, ")"))
expr <- yac_expr(paste0("Factor(", eq, ")"))
expr
eval(expr, list(x = 2))

## -----------------------------------------------------------------------------
eqy <- ysym(eq)
eqy
as_r(eqy)
eqy %>% y_fn("Factor") # Notice how we do not need to call yac_str()/yac_expr()

## -----------------------------------------------------------------------------
A <- outer(0:3, 1:4, "-") + diag(2:5)
a <- 1:4
B <- ysym(A)
B
b <- ysym(a)
b

## -----------------------------------------------------------------------------
y_fn(B, "Transpose")
y_fn(B, "Inverse")
y_fn(B, "Trace")

## -----------------------------------------------------------------------------
A %*% a
B %*% b
t(A)
t(B)
A[, 2:3]
B[, 2:3]
A %*% solve(A)
B %*% solve(B)

## ---- eval = FALSE------------------------------------------------------------
#  hilbert <- function(n) {
#    i <- 1:n
#    H <- 1 / outer(i - 1, i, "+")
#    return(H)
#  }

## -----------------------------------------------------------------------------
hilbert_den <- function(n) { 
  i <- 1:n
  H <- outer(i - 1, i, "+")
  return(H)
}
Hden <- hilbert_den(4)
Hden
H <- 1/Hden
H

## -----------------------------------------------------------------------------
Hyden <- ysym(Hden)
Hyden
Hy <- 1/Hyden
Hy

## -----------------------------------------------------------------------------
as_r(Hy) # now floating-point and the related problems
diag(Hy)
Hy[upper.tri(Hy)]
Hy[1:2, ]
dim(Hy)
A <- Hy
A[lower.tri(A)] <- "x"
A
as_r(A)
eval(as_r(A), list(x = 999))

## -----------------------------------------------------------------------------
x <- ysym("x")
y <- ysym("y")
f <- (1 - x)^2 + 100*(y - x^2)^2
f
tex(f)

## -----------------------------------------------------------------------------
N <- 30
x <- seq(-1, 2, length=N)
y <- seq(-1, 2, length=N)
f_r <- as_r(f)
f_r
z <- outer(x, y, function(x, y) eval(f_r, list(x = x, y = y)))
levels <- c(0.001, .1, .3, 1:5, 10, 20, 30, 40, 50, 60, 80, 100, 500, 1000)
cols <- rainbow(length(levels))
contour(x, y, z, levels = levels, col = cols)

## -----------------------------------------------------------------------------
g <- deriv(f, c("x", "y"))
g

## -----------------------------------------------------------------------------
crit_sol_all <- solve(g, c("x", "y"))
crit_sol_all
crit_sol <- crit_sol_all[1, ] %>% y_rmvars()
crit_sol
crit <- crit_sol %>% as_r()
crit

## -----------------------------------------------------------------------------
H <- Hessian(f, c("x", "y"))
H
tex(H)

## -----------------------------------------------------------------------------
H_crit <- eval(as_r(H), list(x = crit[1], y = crit[2]))
H_crit
eigen(H_crit, only.values = TRUE)$values

