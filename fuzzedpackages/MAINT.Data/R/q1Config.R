q1Config <- function(Config)
{
  if (length(intersect(c(1,2,3),Config)>0)) C1 <- 1
  else  C1 <- NULL
  if (length(intersect(c(4,5),Config)>0)) C5 <- 5
  else C5 <- NULL
  c(C1,C5)
}

q1CovCase <- function(CovCase)
{
  if (length(intersect(c(1,2),CovCase)>0)) C1 <- 1
  else  C1 <- NULL
  if (length(intersect(c(3,4),CovCase)>0)) C4 <- 4
  else C4 <- NULL
  c(C1,C4)
}