## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- message=FALSE-----------------------------------------------------------
library(Ryacas)

## -----------------------------------------------------------------------------
x <- ysym("x")
2*x^2 - 5
c(-2, 5)*x
c(2*x, -x^3)
as_r(c(-2, 5)*x) # or yac_expr(c(-2, 5)*x)

## -----------------------------------------------------------------------------
A <- outer(0:3, 1:4, "-") + diag(2:5)
a <- 1:4
A
a

## -----------------------------------------------------------------------------
B <- ysym(A)
B
as_r(B)
b <- ysym(a)
b
as_r(b)

## -----------------------------------------------------------------------------
y_fn(B, "Transpose")
y_fn(B, "Inverse")
y_fn(B, "Trace")

## -----------------------------------------------------------------------------
A %*% a
B %*% b
t(A)
t(B)
exp(B)
as_r(exp(B))
A[, 2:3]
B[, 2:3]
A[upper.tri(A)] <- 1
B[upper.tri(B)] <- 1
A
B
2*A - A
2*B - B
A %*% solve(A)
B %*% solve(B)
solve(A %*% t(A))
solve(B %*% t(B))
solve(A, a)
solve(B, b)

## -----------------------------------------------------------------------------
yac_str("W") # Get variable W if exists, or else just a symbol
yac_str("Variables()") # ...or list variables
B
yac_assign(B, "W") # assign B in R to W in yacas
yac_str("W") # Get variable W if exists, or else just a symbol
yac_str("Variables()") # ...or list variables
yac_silent("Clear(W)")
yac_str("Variables()") # List variables
yac_str("W") # Get variable W if exists, or else just a symbol

## -----------------------------------------------------------------------------
D <- diag(4) %>% ysym()
D
D <- D/2
D
D[2:3, 1] <- "d"
D[3, 4] <- "2*d + 2"
D
D %>% solve()
D %>% solve() %>% simplify()
D %>% solve() %>% simplify() %>% tex()

## -----------------------------------------------------------------------------
L <- ysym("x^2 * (y/4) - a*(3*x + 3*y/2 - 45)")
L

## -----------------------------------------------------------------------------
deriv(L, "x")
Hessian(L, "x")

## -----------------------------------------------------------------------------
deriv(L, c("x", "y", "a"))
H <- Hessian(L, c("x", "y", "a"))
H
as_r(H)
eval(as_r(H), list(x = 2, y = 2, a = 2))

## -----------------------------------------------------------------------------
L2 <- ysym(c("x^2 * (y/4) - a*(3*x + 3*y/2 - 45)", 
                   "x^3 + 4*a^2")) # just some function
L2
Jacobian(L2, "x")
Jacobian(L2, c("x", "y", "a"))

## -----------------------------------------------------------------------------
xs <- ysym("x")
poly <- xs^2 - xs - 6
poly
zeroes <- solve(poly, "x") # Solve(x^2 - x - 6 == 0, x)
zeroes
tex(zeroes)
zeroes %>% y_rmvars()

## -----------------------------------------------------------------------------
solve(poly, 3, "x") # Solve(x^2 - x - 6 == 3, x)
solve(poly, 3, "x") %>% tex()

## -----------------------------------------------------------------------------
x <- ysym("x")
y <- ysym("y")
lhs <- c(3*x*y - y, x)
rhs <- c(-5*x, y+4)

## -----------------------------------------------------------------------------
sol <- solve(lhs, rhs, c("x", "y"))
sol
sol_vals <- lapply(seq_len(nrow(sol)), function(sol_no) {
  y_rmvars(sol[sol_no, ])
})
sol_vals
sol_envir <- lapply(sol_vals, function(l) {
  list(x = as_r(l[1]), y = as_r(l[2]))
})
sol_envir
do.call(rbind, lapply(seq_along(sol_envir), function(sol_no) {
  sol_val <- sol_envir[[sol_no]]
  data.frame(sol_no = sol_no,
             eq_no = seq_along(sol_val),
             lhs = eval(as_r(lhs), sol_val),
             rhs = eval(as_r(rhs), sol_val))
}))

