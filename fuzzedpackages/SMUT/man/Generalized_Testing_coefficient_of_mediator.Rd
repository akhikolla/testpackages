\name{Generalized_Testing_coefficient_of_mediator}
\alias{Generalized_Testing_coefficient_of_mediator}
\docType{package}
\title{
Testing coefficient of mediator in the outcome model for an outcome following an exponential family distribution or a survival outcome
}
\description{
Testing coefficient of mediator, namely theta, in the outcome model. The outcome model is the following. \cr
outcome ~ intercept + G*gamma + mediator*theta + error
}

\author{
Wujuan Zhong
}

\usage{
Generalized_Testing_coefficient_of_mediator(G,mediator,outcome,
covariates=NULL,outcome_type,
approxi=TRUE,verbose=FALSE)
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
    \item{covariates}{
   n by r matrix (n rows and r columns). Each row is one individual; each column is one covariate.
  }
    \item{outcome_type}{
  Type of the outcome variable. "continuous" for a continuous outcome; "binary" for a binary outcome; "count" for a count outcome; "survival" for a survival outcome.
  }
    \item{approxi}{
  a boolean value. This is an indicator whether the approximation of computing derivatives is applied to save computing time. Default is TRUE.
  }
     \item{verbose}{
  a boolean value. If TRUE a lot of computing details is printed. Default is FALSE.
  }
}

\value{
  \item{p_value}{
    P value for testing the coefficient of mediator in the outcome model.
  }
 
  \item{theta_hat}{
    The point estimate of theta (coefficient of mediator) in the outcome model. 
  }

}



