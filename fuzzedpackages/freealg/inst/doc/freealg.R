## ----set-options, echo = FALSE------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#", dev = "png", fig.width = 7, fig.height = 3.5, message = FALSE, warning = FALSE)
options(width = 80, tibble.width = Inf)

## -----------------------------------------------------------------------------
library("freealg")
X <- freealg(words = list(1, c(24,25), c(25,24), c(1,1,1,2)), coeffs = c(5, 43, 6, -17))
dput(X)
X

## -----------------------------------------------------------------------------
(X <- as.freealg("3aab -2abbax"))  # caret ("^") not yet implemented
(Y <- as.freealg("2 -3aab -aBBAA"))  # uppercase letters are inverses
(Z <- as.freealg(1:3))

## -----------------------------------------------------------------------------
X^2        # powers are implemented
X+Y        #  'aab' term cancels
1000+Y*Z   # algebra multiplication and addition works as expected

## -----------------------------------------------------------------------------
set.seed(0)
phi <- rfalg(n=5,inc=TRUE)
phi
options("usecaret" = TRUE)
phi
options("usecaret" = FALSE)  # reset to default

## -----------------------------------------------------------------------------
X <- as.freealg("x+y+X+Y")
X^2
constant(X^4)

## -----------------------------------------------------------------------------
f <- function(n){constant(as.freealg("x+y+X+Y")^n)}
sapply(c(0,2,4,6,8,10),f)

