# simulate non-linear clusters (half-circles)
simulHalfCircles <- function(n=200, prob=c(100,100), sd = 0.25)
{
  # simul class
  z <- sample(x=c(1,2), size=n, replace=TRUE, prob=prob)
  tau <- runif(n, -4, 4)
  eta <- rnorm(n, 0, sd)
  x<- matrix(ncol=2, nrow = n)
  for (i in 1:n)
  {
    if (z[i] == 1)
    { x [i,] <- t(c(-1 + tau[i]+eta[i],  7 - tau[i]*tau[i]/2 + eta[i])); }
    else
    { x [i,] <- t(c( 1 + tau[i]+eta[i], -7 + tau[i]*tau[i]/2 + eta[i]));}
  }
  list(x=x, z=z)
}
hc<-simulHalfCircles(200)
plot(hc$x,hc$y,xlab='x',ylab='y',col=hc$z+1)

# simul bullsEye
simulBullsEye <- function(n=320, prob = c(80,240), sd = c(0.25, 0.25), radius = c(1,5))
{
  z <- sample(x = c(1,2), size=n, replace=TRUE, prob=prob)
  x<- matrix(0, nrow =n, ncol=2)
  # small circle
  for (i in 1:n)
  {
    theta=runif(1)*2*pi
    r= rnorm(1, mean=radius[z[i]], sd=sd[z[i]])
    x[i,1]<-r*cos(theta)
    x[i,2]<-r*sin(theta)
  }
  list(x=x, z=z)
}
hc<-simulBullsEye()
plot(hc$x,hc$y,xlab='x',ylab='y',col=hc$z+1)

