\name{covOrd}
\alias{covOrd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Warping-Based Covariance for an Ordinal Input}
\description{Creator function for the class \code{\link{covOrd-class}}}
\usage{
covOrd(ordered, 
       k1Fun1 = k1Fun1Matern5_2, 
       warpFun = c("norm", "unorm", "power", "spline1", "spline2"), 
       cov = c("corr", "homo"), 
       hasGrad = TRUE, inputs = "u", 
       par = NULL, parLower = NULL, parUpper = NULL, 
       label = "Ordinal kernel", 
       intAsChar = TRUE, 
       ...)
}
\arguments{
  \item{ordered}{
    An object coerced to \code{\link{ordered}} representing an ordinal input. Only the levels and their order will be used.
  }
  \item{k1Fun1}{
    A function representing a 1-dimensional stationary kernel function, with no or fixed parameters.
}
  \item{warpFun}{
    Name of the increasing warping function.
}
  \item{cov}{
    Character indicating whether a correlation or homoscedastic kernel is used.
}
  \item{hasGrad}{
    Object of class \code{"logical"}. If \code{TRUE}, both \code{k1Fun1} and \code{warpFun} must return the gradient as an attribute of the result.
}
  \item{inputs}{
    Character: name of the ordinal input.
}
  \item{par, parLower, parUpper}{
    Numeric vectors containing covariance parameter values/bounds in the following order: warping, range and variance if required (\code{cov == "homo"}).
}
  \item{label}{
    Character giving a brief description of the kernel.
}
  \item{intAsChar}{
  Logical. If \code{TRUE} (default), an integer-valued input will be coerced into a character. Otherwise, it will be coerced into a factor.
}
  \item{\dots}{
    Not used at this stage.
}

}
\details{
Covariance kernel for qualitative ordered inputs obtained by warping. 

Let \eqn{u} be an ordered factor with levels \eqn{u_1, \dots, u_L}{u[1], ..., u[L]}. 
Let \eqn{k_1}{k1} be a 1-dimensional stationary kernel (with no or fixed parameters), \eqn{F} a warping function i.e. an increasing function on the interval \eqn{[0,1]} and \eqn{\theta}{theta} a scale parameter. Then \eqn{k} is defined by:
\deqn{k(u_i, u_j) = k_1([F(z_i) - F(z_{j})]/\theta)}{k(u[i], u[j]) = k1([F(z[i]) - F(z[j])]/theta)}
where \eqn{z_1, \dots, z_L}{z[1], ..., z[L]} form a regular sequence from \eqn{0} to \eqn{1} (included). At this stage, the possible choices are:
\itemize{
\item A distribution function (cdf) truncated to \eqn{[0,1]}, among the Power and Normal cdfs. 
\item For the Normal distribution, an unnormalized version, corresponding to the restriction of the cdf on \eqn{[0,1]}, is also implemented (\code{warp = "unorm"}).
\item An increasing spline of degree 1 (piecewise linear function) or 2. In this case, \eqn{F} is unnormalized. For degree 2, the implementation depends on scaling functions from DiceKriging package, which must be loaded here.
}
Notice that for unnormalized \code{F}, we set \eqn{\theta}{theta} to 1, in order to avoid overparameterization.
}
\value{
An object of class \code{covOrd-class}, inheriting from \code{covQual-class}.
}
%\note{
%}

\seealso{
\code{\link{covOrd-class}}
}
\examples{
u <- ordered(1:6, labels = letters[1:6])

myCov <- covOrd(ordered = u, cov = "homo", intAsChar = FALSE)
myCov
coef(myCov) <- c(mean = 0.5, sd = 1, theta = 3, sigma2 = 2)
myCov

checkX(myCov, X = data.frame(u = c(1L, 3L)))
covMat(myCov, X = data.frame(u = c(1L, 3L)))

myCov2 <- covOrd(ordered = u, k1Fun1 = k1Fun1Cos, warpFun = "power")
coef(myCov2) <- c(pow = 1, theta = 1) 
myCov2

plot(myCov2, type = "cor", method = "ellipse")
plot(myCov2, type = "warp", col = "blue", lwd = 2)

myCov3 <- covOrd(ordered = u, k1Fun1 = k1Fun1Cos, warpFun = "spline1")
coef(myCov3) <- c(rep(0.5, 2), 2, rep(0.5, 2))
myCov3

plot(myCov3, type = "cor", method = "ellipse")
plot(myCov3, type = "warp", col = "blue", lwd = 2)


str(warpPower)  # details on the list describing the Power cdf
str(warpNorm)   # details on the list describing the Normal cdf
}

