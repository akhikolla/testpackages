\name{flan-package}
\alias{flan-package}
\alias{flan}
\docType{package}
\encoding{UTF-8}



\title{Fluctuation analysis of mutant cell counts}

\description{
 Statistical tools for fluctuation analysis of mutant cells counts.
}

\details{
  \tabular{ll}{
    Package: \tab flan\cr
    Type: \tab Package\cr
    Version: \tab 0.8\cr
  Date: \tab 2020-04-29\cr
    License: \tab GPL\cr
  }

  This package contains functions dedicated to :

  - Distribution (dflan, pflan, qflan, and rflan): built as the classic distribution functions, like dnorm, pnorm, qnorm, and rnorm.

  - Parametric estimation (mutestim): estimates of the parameters of interest, which are the mean number of mutations (or the mutation probability) and the fitness parameter.

  - Parametric testing (flan.test): built as the classic test functions, like t.test.
}

\author{
  Adrien Mazoyer, Rémy Drouilhet, Stéphane Despréaux, and Bernard Ycart

  Maintainer: Adrien Mazoyer <mazoyer.adrien@courrier.uqam.ca>
}

\references{

  A. Mazoyer: Fluctuation analysis on mutation models with birth-date dependence
  Math. Biosci. (2018)

  A. Mazoyer: Time inhomogeneous mutation models with birth-date dependence
  B. Math. Biol. (2017)

  A. Mazoyer, R. Drouilhet, S. Despréaux and B. Ycart: flan: An R Package for Inference on Mutation Models.
  R Journal 9(1) (2017)

  B. Ycart and N. Veziris: Unbiased estimates of mutation rates under fluctuating final counts.
PLoS one 9(7) e101434 (2014)

  B. Ycart: Fluctuation analysis with cell deaths.
 J. Applied Probab. Statist, 9(1):12-28 (2014)

  B. Ycart: Fluctuation analysis: can estimates be trusted?
One PLoS one 8(12) e80958 (2013)

  A. Hamon and B. Ycart: Statistics for the Luria-Delbrück distribution.
Elect. J. Statist., 6:1251-1272 (2012)
}


\keyword{package}
