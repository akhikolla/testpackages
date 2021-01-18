\name{init_mutmodel}
\alias{init_mutmodel}
\alias{print.mutmodel}
\title{
init_mutmodel
}
\description{
Method to initialise a mutation model.
}
\usage{
init_mutmodel(modeltype = 1, mutpars = NULL, ...)
\method{print}{mutmodel}(x, ...)
}
\arguments{
  \item{modeltype}{\preformatted{1: SMM (traditional single-step mutation model).
2: LMM (Logistic mutation model introduced in Jochens (2011) 
   'Empirical Evaluation Reveals Best Fit of a Logistic Mutation 
   Model for Human Y-Chromosomal Microsatellites').
3: Exponential mutation model (unpublished).}}
  \item{mutpars}{A matrix specifying the mutation parameters for each locus. Rows are parameters and columns are loci. If a vector, the same values are used for all loci.}
  \item{x}{A \code{mutmodel} object.}
  \item{\dots}{ not used.}
}
\details{
  Mutation parameters for each locus. \preformatted{
    Mutmodel 1 (SMM): 2 parameters per locus
      P(i -> i-1) = mu\_d
      P(i -> i+1) = mu\_u
      P(i -> i)   = 1 - P(i -> i-1) - P(i -> i+1)
                  = 1 - mu\_d - mu\_u

      mutpars[1, locus]: mu\_d
      mutpars[2, locus]: mu\_u

    Mutmodel 2 (LMM): 6 parameters per locus
      P(i -> i-1) = gamma\_d / (1 + exp(alpha\_d*(beta\_d - i)))
      P(i -> i+1) = gamma\_u / (1 + exp(alpha\_u*(beta\_u - i)))
      P(i -> i)   = 1 - P(i -> i-1) - P(i -> i+1)

      mutpars[1, locus]: gamma\_d
      mutpars[2, locus]: alpha\_d
      mutpars[3, locus]: beta\_d
      mutpars[4, locus]: gamma\_u
      mutpars[5, locus]: alpha\_u
      mutpars[6, locus]: beta\_u

   Mutmodel 3 (EMM): 4 parameters per locus
      P(i -> i-1) = 1/((1 + exp(a + b*i))*(1 + exp(alpha + beta*i)))
      P(i -> i+1) = exp(alpha+beta*i)/((1+exp(a+b*i))*(1+exp(alpha+beta*i)))
      P(i -> i)   = 1 - P(i -> i-1) - P(i -> i+1)
                  = exp(a + b*i)/(1 + exp(a + b*i))

      mutpars[1, locus]: a
      mutpars[2, locus]: b
      mutpars[3, locus]: alpha
      mutpars[4, locus]: beta
  }
}
\value{
  A \code{mutmodel} object (a list with entires \code{modeltype} and \code{mutpars}).
}
\examples{
mutpars <- matrix(c(c(0.003, 0.001), rep(0.004, 2), rep(0.001, 2)), ncol = 3)
colnames(mutpars) <- c("DYS19", "DYS389I", "DYS391")
mutmodel <- init_mutmodel(modeltype = 1L, mutpars = mutpars)
mutmodel

mutpars.locus1 <- c(0.149,   2.08,    18.3,   0.149,   0.374,   27.4) # DYS19
mutpars.locus2 <- c(0.500,   1.18,    18.0,   0.500,   0.0183,  349)  # DYS389I
mutpars.locus3 <- c(0.0163,  17.7,    11.1,   0.0163,  0.592,   14.1) # DYS391
mutpars <- matrix(c(mutpars.locus1, mutpars.locus2, mutpars.locus3), ncol = 3)
colnames(mutpars) <- c("DYS19", "DYS389I", "DYS391")
mutmodel <- init_mutmodel(modeltype = 2L, mutpars = mutpars)
mutmodel
}

