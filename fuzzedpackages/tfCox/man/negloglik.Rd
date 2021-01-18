\name{negloglik}
\alias{negloglik}
\title{
Calculate the negative log likelihood from Cox model. 
}
\description{
Calculate the negative log likelihood from Cox model from the estimated coefficient matrix theta.  
}
\usage{
negloglik(dat, theta)
}
\arguments{
  \item{dat}{
A list that contains \code{time}, \code{status} and \code{X}. \code{time} is failure or censoring time, \code{status} is censoring indicator, and  \code{X} is n x p matrix and may have p > n.
}
  \item{theta}{
  An n x p matrix of coefficients corresponding to covariates \code{X}. 
}
}

\author{
Jiacheng Wu
}

\references{
Jiacheng Wu & Daniela Witten (2019) Flexible and Interpretable Models for Survival Data, Journal of Computational and Graphical Statistics, DOI: 10.1080/10618600.2019.1592758
}

\seealso{
\code{\link{predict_best_lambda}}, \code{\link{tfCox_choose_lambda}}
}

\examples{
#generate training and testing data
dat = sim_dat(n=100, zerof=0, scenario=1)
test_dat = sim_dat(n=100, zerof=0, scenario=1)

#choose the optimal tuning parameter
cv = tfCox_choose_lambda(dat, test_dat, ord=0, alpha=1)
plot(cv$lam_seq, cv$loss)

#optimal tuning parameter
cv$best_lambda

#predict the coefficients of testing covariates from the optimal tuning parameter
#from tfCox_choose_lambda object. 
theta_hat = predict_best_lambda(cv, test_dat$X)

#calculate the loss in the testing data based on the estimated coefficients theta
negloglik(test_dat, theta_hat)
}