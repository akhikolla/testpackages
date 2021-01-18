\name{kmc.bjtest}
\alias{kmc.bjtest}
\title{
Calculate the NPMLE with constriants for accelerated failure time model with given coefficients.
}
\description{
Use the empirical likelihood ratio and Wilks theorem to test if the regression coefficient equals
beta. 
\deqn{El(F)=\prod_{i=1}^{n}(\Delta F(T_i))^{\delta_i}(1-F(T_i))^{1-\delta_i}}
with constraints
\deqn{
	\sum_i g(T_i)\Delta F(T_i)=0,\quad,i=1,2,\ldots
}
Instead of EM algorithm, this function calculates the Kaplan-Meier estimator with mean constraints recursively to test \eqn{H_0:~\beta=\beta_0} in the accelerated failure time model:
\deqn{
	\log(T_i) = y_i = x_i\beta^\top + \epsilon_i,
}
where \eqn{\epsilon} is distribution free.
}
\usage{
kmc.bjtest(y, d, x, beta,init.st="naive")
}

\arguments{
  \item{y}{
      Response variable vector (length n).
}
  \item{d}{
      Status vector (length n), 0: right censored; 1 uncensored.
}
  \item{x}{
      n by p explanatory variable matrix.
}
  \item{beta}{
      The value of the regression coeffiecnt vector (length p) to be tested.
}
  
  \item{init.st}{
Type of methods to initialize the algorithm. By default, init.st is set to 1/n
}
}
\details{
The empirical likelihood is the likelihood of the error term when the coefficients are specified. Model assumptions are the same as requirements of a standard Buckley-James estimator.
}

\value{
 a list with the following components:
 \item{prob}{
 	the probabilities that max the empirical likelihood under estimating equation.
 }
 \item{logel1}{
 	the log empirical likelihood without constraints, i.e. under Kaplan-Merier of residuals'
 }
 \item{logel2}{
 	the log empirical likelihood with constraints, i.e. under null hypotheses or estimation equations.}
 \item{"-2LLR"}{
 	the -2 loglikelihood ratio; have approximate chisq distribution under null hypotheses
 }
 
}
\references{
Buckley, J. and James, I. (1979). Linear regression with censored data. \emph{Biometrika}, \bold{66} 429-36

Zhou, M., & Li, G. (2008).
Empirical likelihood analysis of the Buckley-James estimator. 
\emph{Journal of multivariate analysis}, \bold{99(4)}, 649-664.

Zhou, M. and Yang, Y. (2015).
A recursive formula for the Kaplan-Meier estimator with mean constraints and its application to empirical likelihood
\emph{Computational Statistics}. \bold{Online ISSN} 1613-9658.

}
\author{
    Mai Zhou(mai@ms.uky.edu), Yifan Yang(yifan.yang@uky.edu)

}


\seealso{
	  	\code{\link{plotkmc2D}}, \code{\link[emplik]{bjtest}}.
}
\examples{
x <- c( 1, 1.5, 2, 3, 4.2, 5.0, 6.1, 5.3, 4.5, 0.9, 2.1, 4.3) # positive time
d <- c( 1,   1, 0, 1, 0, 1, 1, 1, 1, 0, 0,   1)               # status censored/uncensored

}

