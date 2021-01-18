
tvp <- function(y,x,V,lambda,W=NULL,kappa=NULL,c=NULL)
{

### estimates time-varying parameters regression model

### y - a numeric or a column matrix of a dependent variable,

### x - a matrix of independent variables (drivers), different columns correspond to different variables

### V - initial variance in the state space equation for the recursive moment estimator updating method,
###     as in the paper by Raftery et al. (2010),

### lambda - a forgetting factor between 0 and 1 used in variance approximations

### W - a method for setting the initial values of variance for the models equations,
###     by default (if W is not specified) the method based on the linear regression,
###     as in the paper by Raftery et al. (2010) is used,
###     alternatively an arbitrary positive number can be specified

### kappa - a parameter in the exponentially weighted moving average, between 0 and 1,
###         if not specified, the method from Raftery et al. (2010) is applied

### c - a parameter indicating whether constant is included,
###     by default c=TRUE (constant is included),
###     it is not possible to set c=FALSE if ncol(x)=0


if ( is.null(c) ) { c <- TRUE }

x <- as.matrix(x)
if ( ncol(x) == 0 ) { c <- TRUE }

xe <- cbind(rep.int(1,nrow(x)),x)


if ( c == TRUE )
  {
    theta <- matrix(0,ncol=1, nrow=ncol(xe))
  }
else
  {
    theta <- matrix(0,ncol=1, nrow=ncol(xe)-1)
  }

if (is.null(W))
{
  if (ncol(x)>0)
    {
      beta <- lm(y~x)$coefficients[1]
    }
  else
    {
      beta <- mean(y)
    }
  if (length(y)>1)
    {
      beta <- as.numeric(beta * beta + var(as.numeric(y)))
    }
  else
    {
      beta <- as.numeric(beta * beta)
    }

  E <- diag(beta,ncol=1+ncol(x), nrow=1+ncol(x))

  v <- vector()
  if (ncol(x)>0 && length(y)>1)
    {
      for (i in 1:ncol(x))
        {
          v[i] <- var(x[,i])
          if (v[i]==0) { v[i] <- 0.001 * (1/(2^ncol(x))) }
        }
    }
  if (length(y)>1)
    {
      vary <- as.numeric(var(as.numeric(y)))
    }

  if (ncol(x)>0 && length(y)>1)
    {
      for (j in 1:length(v))
        {
          E[j+1,j+1] <- vary/v[j]
        }
    }
  else
    {
      if (ncol(x)>0)
        {
          for (j in 1:ncol(x))
            {
              E[j+1,j+1] <- 0
            }
        }
    }
}

if (!is.null(W))
{
  E <- diag(W,ncol=1+ncol(x), nrow=1+ncol(x))
}

#########################
#########################

if (c == FALSE)
  {
    xe <- xe[,-1,drop=FALSE]
    E <- E[-1,-1,drop=FALSE]
  }

tvpcppout <- tvpcpp(x,y,xe,theta,E,lambda,V,kappa)

thetas <- tvpcppout[[1]]
y.tvp <- tvpcppout[[2]]
pdensi <- tvpcppout[[3]]

thetas <- t(thetas[,-ncol(thetas)])
if (!is.null(colnames(x)))
  {
    if ( c == TRUE )
      {
        colnames(thetas) <- c("const",colnames(x))
      }
    else
      {
        if (ncol(x)==1) { thetas <- t(thetas) }
        colnames(thetas) <- colnames(x)
      }
  }
if (ncol(x)==0)
  {
    thetas <- t(thetas)
    colnames(thetas) <- "const"
  }

y.tvp <- as.vector(y.tvp)
ret <- list(y.tvp,thetas,pdensi,as.matrix(y))
names(ret) <- c("y.hat","thetas","pred.dens.","y")
class(ret) <- "tvp"

return(ret)

}
