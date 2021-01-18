\name{SMUT}
\alias{SMUT}
\docType{package}
\title{
Multi-SNP Mediation Intersection-Union Test
}
\description{
Testing the mediation effect of multiple SNPs on an outcome through a mediator.
}

\author{
Wujuan Zhong
}

\references{
Zhong, W., Spracklen, C. N., Mohlke, K. L., Zheng, X., Fine, J., & Li, Y. (2019). Multi-SNP mediation intersection-union test. \emph{Bioinformatics}.
}

\usage{
SMUT(G, mediator, outcome, 
outcome_type="continuous", method="score",approxi=TRUE, debug=FALSE)
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
  \item{p_value_IUT}{
    The p value for testing the mediation effect (theta*beta) based on intersection-union test.
  }
  \item{p_value_theta}{
    The p value for testing theta in the outcome model.
    The outcome model is the following. \cr
    outcome ~ intercept + G*gamma + mediator*theta + error.
  }
  \item{p_value_beta}{
     The p value for testing beta in the mediator model.
     The mediator model is the following. \cr
     mediator ~ intercept + G*beta + error
  }

}

\examples{
library(SMUT)
# load the Genotype data included in this R package
data("Genotype_data")

# generate one mediator and one outcome

# first example, the mediation effect is significant
set.seed(1)
beta=rnorm(ncol(Genotype_data),1,2)
e1 = rnorm(nrow(Genotype_data), 0, 1)
mediator = 1 + eigenMapMatMult(Genotype_data,beta) + e1

theta=0.8
gamma=rnorm(ncol(Genotype_data),0.5,2)
e2 = rnorm(nrow(Genotype_data), 0, 1)
outcome = 2 + eigenMapMatMult(Genotype_data,gamma) + theta*mediator + e2

p_value=SMUT(G=Genotype_data,mediator=mediator,outcome=outcome)
print(p_value)

# p_value_IUT is the p value for the mediation effect.
# we have significant(at alpha level 0.05) mediation effects (p_value_IUT = 0.001655787).

# second example, the mediation effect is non-significant
set.seed(1)
beta=rnorm(ncol(Genotype_data),1,2)
e1 = rnorm(nrow(Genotype_data), 0, 1)
mediator = 1 + eigenMapMatMult(Genotype_data,beta) + e1

theta=0
gamma=rnorm(ncol(Genotype_data),0.5,2)
e2 = rnorm(nrow(Genotype_data), 0, 1)
outcome = 2 + eigenMapMatMult(Genotype_data,gamma) + theta*mediator + e2

p_value=SMUT(G=Genotype_data,mediator=mediator,outcome=outcome)
print(p_value)

# p_value_IUT is the p value for the mediation effect.
# we have non-significant(at alpha level 0.05) mediation effects (p_value_IUT = 0.3281677).

# third example, the mediation effect is non-significant
set.seed(1)
beta=rep(0,ncol(Genotype_data))
e1 = rnorm(nrow(Genotype_data), 0, 1)
mediator = 1 + eigenMapMatMult(Genotype_data,beta) + e1

theta=0.8
gamma=rnorm(ncol(Genotype_data),0.5,2)
e2 = rnorm(nrow(Genotype_data), 0, 1)
outcome = 2 + eigenMapMatMult(Genotype_data,gamma) + theta*mediator + e2

p_value=SMUT(G=Genotype_data,mediator=mediator,outcome=outcome)
print(p_value)

# p_value_IUT is the p value for the mediation effect.
# we have non-significant(at alpha level 0.05) mediation effects (p_value_IUT = 0.5596977).

# Thanks for using our R package SMUT 
}

