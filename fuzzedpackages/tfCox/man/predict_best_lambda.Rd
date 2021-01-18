\name{predict_best_lambda}
\alias{predict_best_lambda}
\title{
Predict from the optimal lambda from tfCox_choose_lambda
}
\description{
Estimate the corresponding theta values from the optimal tuning parameter obtained by \code{\link{tfCox_choose_lambda}}.
}
\usage{
predict_best_lambda(cv, newX)
}
\arguments{
  \item{cv}{
  An object from tfCox_choose_lambda. 
}
  \item{newX}{
  The new covariate values. 
}
}

\value{
Estimated theta values. 
}

\author{
Jiacheng Wu
}

\references{
Jiacheng Wu & Daniela Witten (2019) Flexible and Interpretable Models for Survival Data, Journal of Computational and Graphical Statistics, DOI: 10.1080/10618600.2019.1592758
}

\seealso{
\code{\link{tfCox_choose_lambda}}, \code{\link{negloglik}}
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

#Estimate the theta values of testing covariates from the optimal tuning parameter
#from tfCox_choose_lambda object. 
theta_hat = predict_best_lambda(cv, test_dat$X)
}