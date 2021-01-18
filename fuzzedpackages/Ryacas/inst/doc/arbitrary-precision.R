## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- message=FALSE-----------------------------------------------------------
library(Ryacas)

## -----------------------------------------------------------------------------
(0.3 - 0.2) - 0.1
yac_str("(0.3 - 0.2) - 0.1") # decimal number does not 
                             # always work well; often 
                             # it is better to represent as 
                             # rational number if possible

# (1/3 - 1/5) = 5/15 - 3/15 = 2/15
(1/3 - 1/5) - 2/15
yac_str("(1/3 - 1/5) - 2/15")

## -----------------------------------------------------------------------------
yac_str("1/3")
yac_str("N(1/3)")
yac_str("N(1/3, 1)")
yac_str("N(1/3, 200)")

## -----------------------------------------------------------------------------
pol <- "(x-1)^7"
p1 <- pol %>% yac_expr()
p1
p2 <- pol %>% y_fn("Expand") %>% yac_expr()
p2

## -----------------------------------------------------------------------------
eval(p1, list(x = 1.001))
eval(p2, list(x = 1.001))

## -----------------------------------------------------------------------------
# First try with 1.001:
pol_val <- paste0("WithValue(x, 1.001, ", pol, ")")
pol_val
yac_str(pol_val)
#... to get result symbolically, use instead the number as a fraction
pol_val <- paste0("WithValue(x, 1001/1000, ", pol, ")")
pol_val
yac_str(pol_val)
pol_val %>% y_fn("Denom") %>% yac_str()
pol_val %>% y_fn("Denom") %>% y_fn("IntLog", "10") %>% yac_str()
pol_val <- paste0("WithValue(x, 1001/1000, Expand(", pol, "))")
pol_val
yac_str(pol_val)

## ---- fig.height=4, fig.width = 6, fig.pos="H", fig.cap="solid line: $(x-1)^7$; dashed line: factorization of $(x-1)^7$.", echo=FALSE----
xval <- seq(.99, 1.01, by = 0.001)
p1_val <- sapply(xval, function(x) eval(p1, list(x = x)))
p2_val <- sapply(xval, function(x) eval(p2, list(x = x)))
plot(xval, p1_val, type = "l", xlab = "x", ylab = "y")
lines(xval, p2_val, lty = 2)
legend("bottomright",
       legend = c(
         parse(text = paste0("p[1]: y == ", as.character(p1))),
         parse(text = paste0("p[2]: y == ", as.character(p2)))),
       lty = c(1, 2),
       cex = 0.7)

## -----------------------------------------------------------------------------
hilbert <- function(n) { 
  i <- 1:n
  H <- 1 / outer(i - 1, i, "+")
  return(H)
}

## -----------------------------------------------------------------------------
hilbert_sym <- function(n) { 
  mat <- matrix("", nrow = n, ncol = n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      mat[i, j] <- paste0("1 / (", (i-1), " + ", j, ")")
    }
  }
  
  return(mat)
}

## -----------------------------------------------------------------------------
A <- hilbert(8)
A
B1 <- hilbert_sym(8)
B1
B <- ysym(B1) # convert to yacas symbol
B

## -----------------------------------------------------------------------------
Ainv <- solve(A)
Ainv
Binv1 <- solve(B) # result is still yacas symbol
Binv1
Binv <- as_r(Binv1) # convert to R numeric matrix
Binv
Ainv - Binv
max(abs(Ainv - Binv))
max(abs((Ainv - Binv) / Binv))

