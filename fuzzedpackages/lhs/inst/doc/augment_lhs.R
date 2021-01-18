## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
source("VignetteCommonCode.R")
require(lhs)

graph2DaugmentLHS1 <- function(sims, extras)
{
  A <- randomLHS(sims, 2)
  B <- augmentLHS(A, extras)

  plot.default(A[,1], A[,2], type = "n", ylim = c(0,1),
    xlim = c(0,1), xlab = "x1", ylab = "x2", xaxs = "i", yaxs = "i", main = "")
  for (i in 1:length(A[,1]))
  {
    rect(floor(A[i,1]*sims)/sims, floor(A[i,2]*sims)/sims,
      ceiling(A[i,1]*sims)/sims, ceiling(A[i,2]*sims)/sims, col = "grey")
  }
  points(A[,1], A[,2], pch = 19, col = "red")
  abline(v = (0:sims)/sims, h = (0:sims)/sims)
  
  return(list(A = A, B = B, sims = sims, extras = extras))
}

graph2DaugmentLHS2 <- function(X)
{
  A <- X$A
  B <- X$B
  sims <- X$sims
  extras <- X$extras

  plot.default(A[,1], A[,2], type = "n", ylim = c(0,1),
    xlim = c(0,1), xlab = "x1", ylab = "x2", xaxs = "i", yaxs = "i", main = "")
  N <- sims + extras
  for (i in 1:length(B[,1]))
  {
    rect(floor(B[i,1]*N)/N, floor(B[i,2]*N)/N,
      ceiling(B[i,1]*N)/N, ceiling(B[i,2]*N)/N, col = "grey")
  }
  points(A[,1], A[,2], pch = 19, col = "red")
  points(B[((sims + 1):(sims + extras)), 1], B[((sims + 1):(sims + extras)), 2],
    pch = 19, col = "blue")
  abline(v = (0:N)/N, h = (0:N)/N)
}

# X <- graph2DaugmentLHS1(5,5)
# graph2DaugmentLHS2(X)

## ----randomlhs----------------------------------------------------------------
A <- randomLHS(5,2)

## ----original5, echo=FALSE, fig.align='center', fig.height=5, fig.width=5-----
set.seed(10)
X <- graph2DaugmentLHS1(5, 5)

## ----augment5-----------------------------------------------------------------
B <- augmentLHS(A, 5)

## ----augmented10, fig.align='center', echo=FALSE, fig.height=5, fig.width=5----
graph2DaugmentLHS2(X)

## ----random_and_augment-------------------------------------------------------
A <- randomLHS(7, 2)
B <- augmentLHS(A, 3)

## ----Z, echo=FALSE, fig.align='center', fig.height=5, fig.width=5-------------
set.seed(12)
X <- graph2DaugmentLHS1(7, 3)

## ----W, echo=FALSE, fig.align='center', fig.height=5, fig.width=5-------------
graph2DaugmentLHS2(X)

