\name{gMWT-package}
\alias{gMWT-package}
\alias{gMWT}
\docType{package}

\title{
  Implementations of Generalized Mann-Whitney Type Tests
}

\description{
The package provides nonparametric tools for the comparison of several groups/treatments when 
the number of variables is large. The tools are the following.

(i) Computation of the so-called probabilistic index (PI) estimates for the group comparisons.
The probabilistic indices here are (a) the probability \eqn{\hat P_t} that a random observation from group \eqn{t}
is smaller than a random observation from other groups, (b) the probability \eqn{\hat P_{tt'}} that 
a random observation from group \eqn{t} is smaller than a random observation from group \eqn{t'},
and (c) the probability \eqn{\hat P_{tt't''}} that random observations from groups \eqn{t,t',t''} 
appear in this same order. The tools are also given to produce the plots of PIs for set of variables.

(ii) Computation of the p-values of some classical  and some new nonparametric tests for the comparison of
several groups/treatments.  The tests are based on the use of the  probabilistic 
indices (PI) \eqn{\hat P_t}, \eqn{\hat P_{tt'}} and \eqn{\hat P_{tt't''}}. Classical Wilcoxon-Mann-Whitney test, 
Kruskal-Wallis and Jonckheere-Terpstra tests are included.

(iii) Tools for the simultaneous testing problem. As the package is meant for the analysis of gene expression data, 
for example, the plots of the numbers of observed and expected rejected null hypotheses and their ratio as
a function of alpha are provided. Different multiple testing adjustment methods are implemented, i.a. Bonferroni,
Bejamini & Hochberg and Westfall & Young.  
}



\details{
  \tabular{ll}{
    Package: \tab GMW\cr
    Type: \tab Package\cr
    Version: \tab 1.1.1\cr
    Date: \tab 2020-02-11\cr
    License: \tab GPL\cr
    LazyLoad: \tab yes\cr
  }
}

\author{
  Daniel Fischer, Hannu Oja

  Maintainer: Daniel Fischer <daniel.fischer@luke.fi>
}

\references{
  Daniel Fischer, Hannu Oja (2015). Mann-Whitney Type Tests for Microarray Experiments: The R Package gMWT. Journal of Statistical Software, 65(9), 1-19. URL http://www.jstatsoft.org/v65/i01/.

  Fischer, D., Oja, H., Schleutker, J., Sen, P.K., Wahlfors, T. (2013): Generalized Mann-Whitney Type Tests for Microarray Experiments, Scandinavian Journal of Statistic. doi: 10.1111/sjos.12055.
}

\keyword{package}