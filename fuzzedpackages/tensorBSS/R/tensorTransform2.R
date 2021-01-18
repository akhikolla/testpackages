tensorTransform2 <- function(x, A, mode, transpose=FALSE)
{
  r <- length(dim(x)) - 1
  ind <- (1:r)[mode]
  aB <- ifelse(transpose,1,2)
  
  for (i in ind)
  {
    x <- tensor(x, A[[i]], alongA = i, alongB = aB)
    x <- aperm(x, cyclicalPermute(i, r))
  }
  x
}


# # Test
# 
# n <- 5
# x <- array(rnorm(5*6*7), dim=c(7, 6, 5))
# A1 <- matrix(runif(14),ncol=7)
# A2 <- matrix(rexp(18),ncol=6)
# A  <- list(A1=A1, A2=A2)
# At <- list(tA1 = t(A1), tA2 = t(A2))
# 
# x1 <- tensorTransform2(x,A,1)
# x2 <- tensorTransform2(x,A,-2)
# x3 <- tensorTransform(x,A1,1)
# x1==x2
# x1==x3
# x4 <- tensorTransform2(x,At,-2, TRUE)
# x1 ==x4
# x5 <- tensorTransform2(x,A,1:2)