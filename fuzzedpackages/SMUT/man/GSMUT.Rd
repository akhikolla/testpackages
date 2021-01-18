\name{GSMUT}
\alias{GSMUT}
\docType{package}
\title{
Generalized Multi-SNP Mediation Intersection-Union Test
}
\description{
Testing the mediation effect of multiple SNPs on an outcome following an exponential family distribution or a survival outcome through a continuous mediator.
}

\author{
Wujuan Zhong
}

\usage{
GSMUT(G,mediator,outcome,covariates=NULL,outcome_type,
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
  \item{p_value_IUT}{
    The p value for testing the mediation effect (theta*beta) based on intersection-union test.
  }
  \item{p_value_theta}{
    The p value for testing theta in the outcome model.
    The outcome model is the following. \cr
    outcome ~ intercept + covariates*iota + G*gamma + mediator*theta 
  }
  \item{theta_hat}{
    The point estimate of theta (coefficient of mediator) in the outcome model. 
  }
  \item{p_value_beta}{
     The p value for testing beta in the mediator model.
     The mediator model is the following. \cr
     mediator ~ intercept + covariates*iota + G*beta + error
  }

}

\examples{
library(SMUT)
# load the Genotype data included in this R package
data("Genotype_data")


##### for a binary outcome #####
set.seed(1)

# generate two covariates
covariate_1=rnorm(nrow(Genotype_data),0,1)
covariate_2=sample(c(0,1),size=nrow(Genotype_data),replace = TRUE)
covariates=cbind(covariate_1,covariate_2)

# generate a mediator 
beta=rnorm(ncol(Genotype_data),0,0.5)
tau_M=c(-0.3,0.2)
e1 = rnorm(nrow(Genotype_data), 0, 1)
mediator = 1 + eigenMapMatMult(Genotype_data,beta) + 
           eigenMapMatMult(covariates, tau_M) + e1

#### generate a binary outcome ####
theta=1
gamma=rnorm(ncol(Genotype_data),0,0.5)
tau=c(-0.2,0.2)
eta=1 + eigenMapMatMult(Genotype_data, gamma) + 
    eigenMapMatMult(covariates, tau) + theta * mediator
pi=1/(1+exp( -(eta ) ))
outcome=rbinom(length(pi),size=1,prob=pi) 
result=GSMUT(G=Genotype_data,mediator=mediator,outcome=outcome,
covariates=covariates,outcome_type="binary")
print(result)
# p_value_IUT is the p value for the mediation effect.

\dontrun{

##### generate a count outcome #####
theta=1
gamma=rnorm(ncol(Genotype_data),0,0.5)
tau=c(-0.2,0.2)
eta=1 + eigenMapMatMult(Genotype_data, gamma) + 
    eigenMapMatMult(covariates, tau) + theta * mediator
mu_param=exp(eta) # the mean parameter
phi_param=10  # the shape parameter
outcome=rnbinom(length(mu_param),size=phi_param,mu=mu_param)
result=GSMUT(G=Genotype_data,mediator=mediator,outcome=outcome,
covariates=covariates,outcome_type="count")
print(result)
# p_value_IUT is the p value for the mediation effect.

##### generate a survival outcome #####
theta=2
gamma=rnorm(ncol(Genotype_data),0,0.5)
tau=c(-0.2,0.2)
eta=1 + eigenMapMatMult(Genotype_data, gamma) + 
    eigenMapMatMult(covariates, tau) + theta * mediator
v=runif(nrow(Genotype_data))
lambda=0.01; rho=1; rateC=0.001
Tlat=(- log(v) / (lambda * exp( eta  )))^(1 / rho)
# censoring times
C= rexp(nrow(Genotype_data), rate=rateC)
# follow-up times and event indicators
time= pmin(Tlat, C)
status= as.numeric(Tlat <= C)
outcome=cbind(time,status)      
colnames(outcome)=c("time","status")
result=GSMUT(G=Genotype_data,mediator=mediator,outcome=outcome,
covariates=covariates,outcome_type="survival")
print(result)
# p_value_IUT is the p value for the mediation effect.

} 


}

