\name{cglm}

\alias{cglm}

\title{
Estimates the ratio of the regression coefficients and the dispersion parameter
in conditional generalized linear models for clustered data. 
}

\description{
\code{cglm} estimates the ratio of the regression coefficients and the 
dispersion parameter in conditional generalized linear models. This 
is of particular interest in the so-called case-time-control design.
}     

\usage{
cglm(method, formula, data, id, link, ...)
}

\arguments{
\item{method}{
a string specifying the desired estimation method; either \code{"ts"} for 
two-step estimation, or \code{"cml"} for conditional maximum likelihood estimation.  
}
\item{formula}{
a symbolic description of the model to be fitted. 
}
  \item{data}{
a data frame containing the variables in the model. 
}  
  \item{id}{
a string containing the name of the cluster identification variable. 
}
  \item{link}{
  a string specifying the desired link function. This argument is not used 
  when \code{method="cml"}.   
}
  \item{\dots}{
  optional arguments passed on to the \code{nleqslv} function or the \code{optim}
  function, which are used to solve the estimating equations when 
  \code{method="ts"} and \code{method="cml"}, respectively. See the help pages 
  for \code{nleqslv} and \code{optim}.      
}

}

\details{ 
Let \eqn{y_{ij}} and \eqn{x_{ij}} be the outcome and covariate(s) for subject 
\eqn{j} in cluster \eqn{i}, respectively. Consider the conditional generalized 
linear model
\deqn{p(y_{ij}|i,x_{ij})=\textrm{exp}\left[\frac{\theta_{ij}y_{ij}-A(\theta_{ij})}{\phi}+k(y_{ij},\phi)\right]}
where 
\deqn{\theta_{ij}=\eta\{E(y_{ij}|x_{ij})\}=b_i+\beta x_{ij}.}
\code{cglm} estimates the ratio 
\deqn{\beta / \phi.}
This ratio is of particular interest in so-called case-time-control designs;
see Sjolander (2016) and Sjolander and Ning (2018) for details. 
Two estimation methods are allowed; the two-step method proposed by Sjolander (2016) and the 
conditional maximum likelihood method proposed by Sjolander and Ning (2018).    

}

\value{
An object of class \code{"cglm"} is a list containing 
\item{call}{
  the matched call. 
  }
\item{coefficients}{
  the ratio of the estimated coefficients and the estimated dispersion parameter.     
  }
\item{var}{
  the variance-covariance matrix. 
  }
\item{convergence}{
  was a solution found to the estimating equations?  
  }

}

\note{

Missing data are not allowed.

}

\references{
Sjolander A. (2017). The case-time-control method for non-binary exposures. 
\emph{Sociological Methodology} \bold{47}(1), 182-211.

Sjolander A., Ning Y. (2018). A general and robust estimation method for 
the case-time-control design. 
\emph{Sociological Methodology} \bold{49}(1), 349-365.

}

\author{
Arvid Sjolander. 
}


\examples{

data(teenpov)

fit.ide <- cglm(method="ts", formula=hours~nonpov+inschool+spouse+age+mother, 
  data=teenpov, id="ID", link="identity")
summary(fit.ide)

fit.log <- cglm(method="ts", formula=hours~nonpov+inschool+spouse+age+mother, 
  data=teenpov, id="ID", link="log")
summary(fit.log)

fit.cglm <- cglm(method="cml", formula=hours~nonpov+inschool+spouse+age+mother,
  data=teenpov, id="ID")
summary(fit.cglm)

}
