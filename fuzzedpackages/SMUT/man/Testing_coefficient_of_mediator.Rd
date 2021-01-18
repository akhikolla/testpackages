\name{Testing_coefficient_of_mediator}
\alias{Testing_coefficient_of_mediator}
\docType{package}
\title{
Testing coefficient of mediator in the outcome model
}
\description{
Testing coefficient of mediator, namely theta, in the outcome model. The outcome model is the following. \cr
outcome ~ intercept + G*gamma + mediator*theta + error
}



\author{
Wujuan Zhong
}

\usage{
Testing_coefficient_of_mediator(G, mediator, outcome, 
outcome_type="continuous", method="score", approxi=TRUE, debug=FALSE)
}
\arguments{
    \item{G}{
  n by p matrix (n rows and p columns). Each row is one individual; each column is one SNP.
  }
    \item{mediator}{
  a vector length of n. It is the mediator variable.
  }
    \item{outcome}{
  a vector length of n. It is the outcome variable.
  }
    \item{outcome_type}{
  Type of the outcome variable. For now, this package only deals with continuous outcome. Default is "continuous".
  }
    \item{method}{
  The method of testing coefficient of mediator in the outcome model. The score test is used. Default is "score".
  }
    \item{approxi}{
  a boolean value. This is an indicator whether the approximation of the score statistic is applied to save computing time. Default is TRUE.
  }
     \item{debug}{
  a boolean value. If TRUE a lot of computing details is printed; otherwise the
function is completely silent. Default is FALSE.
  }
}

\value{
  P value for testing the coefficient of mediator in the outcome model.
}

\examples{
library(SMUT)
# load the Genotype data included in this R package
data("Genotype_data")

# generate one mediator and one outcome

set.seed(1)
beta=rnorm(ncol(Genotype_data),1,2)
e1 = rnorm(nrow(Genotype_data), 0, 1)
mediator = 1 + eigenMapMatMult(Genotype_data,beta) + e1

theta=0.8
gamma=rnorm(ncol(Genotype_data),0.5,2)
e2 = rnorm(nrow(Genotype_data), 0, 1)
outcome = 2 + eigenMapMatMult(Genotype_data,gamma) + theta*mediator + e2

p_value=Testing_coefficient_of_mediator(G=Genotype_data,mediator=mediator,outcome=outcome)
print(p_value)

# Thanks for using our R package SMUT 
}

