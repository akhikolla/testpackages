# Helper function

######################
# Online covariance  #
# Helper functions   #
# ####################
makeOnlineCov <- function(d)
{
  e <- new.env()
  e$n <- 0
  e$means <- rep(0, d)
  e$M2 <- matrix(rep(0,d*d),nrow=d)
  e
}

updateOnlineCov <- function(e, x)
{
  x <- as.numeric(x)

  # n increases by 1
  e$n <- e$n + 1

  # Difference between old mean and new observation (see online mean)
  dx <- x - e$means

  # Update means
  e$means <- e$means + dx/e$n

  # Difference between new mean and observation
  dy <- x - e$means

  # Outer product of differences
  otr <- outer(dx,dy)

  # Add to the outer product matrix M2
  e$M2 <- e$M2 + otr
}

getOnlineCov <- function(e)
{
  e$M2 / (e$n - 1)
}

