## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

require(lhs)
source("VignetteCommonCode.R")

graph2dLHS <- function(Alhs)
{
  stopifnot(ncol(Alhs) == 2)
  sims <- nrow(Alhs)
  par(mar = c(4,4,2,2))
  plot.default(Alhs[,1], Alhs[,2], type = "n", ylim = c(0,1),
    xlim = c(0,1), xlab = "Parameter 1", ylab = "Parameter 2", xaxs = "i", 
    yaxs = "i", main = "")
  for (i in 1:nrow(Alhs))
  {
    rect(floor(Alhs[i,1]*sims)/sims, floor(Alhs[i,2]*sims)/sims,
      ceiling(Alhs[i,1]*sims)/sims, ceiling(Alhs[i,2]*sims)/sims, col = "grey")
  }
  points(Alhs[,1], Alhs[,2], pch = 19, col = "red")
  abline(v = (0:sims)/sims, h = (0:sims)/sims)
}

# transform is a function of the kind that takes a number
# transform <- function(x){return(qnorm(x,mean=0, std=1))}
graph2dLHSTransform <- function(Alhs, transform1, transform2, min1, max1, min2, max2)
{
  stopifnot(ncol(Alhs) == 2)
  stopifnot(all(Alhs[,1] <= max1 && Alhs[,1] >= min1))
  stopifnot(all(Alhs[,2] <= max2 && Alhs[,2] >= min2))
  sims <- nrow(Alhs)
  breaks <- seq(0,1,length = sims + 1)[2:(sims)]
  breaksTransformed1 <- sapply(breaks, transform1)
  breaksTransformed2 <- sapply(breaks, transform2)

  par(mar = c(4,4,2,2))
  plot.default(Alhs[,1], Alhs[,2], type = "n", 
               ylim = c(min2, max2),
               xlim = c(min1, max1),
               xlab = "Parameter 1", ylab = "Parameter 2", 
               xaxs = "i", yaxs = "i", main = "")
  for (si in 1:sims)
  {
    temp <- Alhs[si,]
    for (i in 1:sims)
    {
      if ((i == 1 && min1 <= temp[1] && breaksTransformed1[i] >= temp[1]) ||
            (i == sims && max1 >= temp[1] && breaksTransformed1[i - 1] <= temp[1]) ||
            (breaksTransformed1[i - 1] <= temp[1] && breaksTransformed1[i] >= temp[1]))
      {
        for (j in 1:sims)
        {
          if ((j == 1 && min2 <= temp[2] && breaksTransformed2[j] >= temp[2]) ||
                (j == sims && max2 >= temp[2] && breaksTransformed2[j - 1] <= temp[2]) ||
                (breaksTransformed2[j - 1] <= temp[2] && breaksTransformed2[j] >= temp[2]))
          {
            if (i == 1)
            {
              xbot <- min1
              xtop <- breaksTransformed1[i]
            } else if (i == sims)
            {
              xbot <- breaksTransformed1[i - 1]
              xtop <- max1
            } else 
            {
              xbot <- breaksTransformed1[i - 1]
              xtop <- breaksTransformed1[i]
            }
            if (j == 1)
            {
              ybot <- min2
              ytop <- breaksTransformed2[j]
            } else if (j == sims)
            {
              ybot <- breaksTransformed2[j - 1]
              ytop <- max2
            } else 
            {
              ybot <- breaksTransformed2[j - 1]
              ytop <- breaksTransformed2[j]
            }
            rect(xbot, ybot, xtop, ytop, col = "grey")
          }
          
        }
      }
    }
  }
  points(Alhs[,1], Alhs[,2], pch = 19, col = "red")
  abline(v = breaksTransformed1, h = breaksTransformed2)
}

#set.seed(1111)
#A <- randomLHS(5,4)
#f <- function(x){qnorm(x)}
#g <- function(x){qlnorm(x, meanlog=0.5, sdlog=1)}
#B <- A
#B[,1] <- f(A[,1])
#B[,2] <- g(A[,2])
#graph2dLHSTransform(B[,1:2], f, g, -4, 4, 0, 8)
#f <- function(x){qunif(x, 3, 5)}
#B <- apply(A, 2, f)
#graph2dLHSTransform(B[,1:2], f)


## ----block1-------------------------------------------------------------------
# set the seed for reproducibility
set.seed(1111)
# a design with 5 samples from 4 parameters
A <- randomLHS(5, 4) 
A

## ----figureX, fig.align='center', fig.height=5, fig.width=5, echo=FALSE-------
graph2dLHS(A[,1:2])

## ----block 3------------------------------------------------------------------
B <- matrix(nrow = nrow(A), ncol = ncol(A))
B[,1] <- qnorm(A[,1], mean = 0, sd = 1)
B[,2] <- qlnorm(A[,2], meanlog = 0.5, sdlog = 1)
B[,3] <- A[,3]
B[,4] <- qunif(A[,4], min = 7, max = 10)
B

## ----figureY, fig.align='center', fig.height=5, fig.width=5, echo=FALSE-------
f <- function(x){qnorm(x)}
g <- function(x){qlnorm(x, meanlog = 0.5, sdlog = 1)}
graph2dLHSTransform(B[,1:2], f, g, -4, 4, 0, 8)

## ----block 4------------------------------------------------------------------
set.seed(101)
A <- randomLHS(30, 10)
A1 <- optimumLHS(30, 10, maxSweeps = 4, eps = 0.01)
A2 <- maximinLHS(30, 10, dup = 5)
A3 <- improvedLHS(30, 10, dup = 5)
A4 <- geneticLHS(30, 10, pop = 1000, gen = 8, pMut = 0.1, criterium = "S")
A5 <- geneticLHS(30, 10, pop = 1000, gen = 8, pMut = 0.1, criterium = "Maximin")

## ----Z, fig.align='center', fig.height=7, fig.width=7, echo=FALSE-------------
pairs(A, pch = 19, col = "blue", cex = 0.5)

## ----W, fig.align='center', fig.height=7, fig.width=7, echo=FALSE-------------
pairs(A1, pch = 19, col = "blue", cex = 0.5)

## ----G, fig.align='center', fig.height=7, fig.width=7, echo=FALSE-------------
pairs(A2, pch = 19, col = "blue", cex = 0.5)

