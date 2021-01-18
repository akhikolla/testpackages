\name{hwe}
\alias{hwe}
\title{
Relative genotypic frequencies under HWE
}
\description{
  Compute relative genotypic frequencies for a given
  relative allelic frequency
}
\usage{
hwe(raf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{raf}{
    Relative minor allele frequency for the B allele.
  }
}
\details{
  For a bi-allelic SNP with genotypes \emph{AA}, \emph{AB} and
  \emph{ BB},
  with a relative allele frequency \eqn{q}{q} for the \emph{B} allele,
  this function returns the corresponding relative genotypic
  frequencies.
  It is not intended to be called directly by the user.
  To conduct power calculations, use \link{sim.snp.expsurv.power}
  or the convenience wrapper function \link{survSNP.power.table}.
}
\value{
A vector of length three relative genotypic frequencies.
}
\references{
  Kouros Owzar, Zhiguo Li, Nancy Cox and Sin-Ho Jung.
  Power and Sample Size Calculations for SNP Association
  Studies with Censored Time-to-Event Outcomes.
  \url{http://onlinelibrary.wiley.com/doi/10.1002/gepi.21645/abstract}
}

\author{
Kouros Owzar, Zhiguo Li, Nancy Cox, Sin-Ho Jung and Chanhee Yi
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
censbnd(0.1,hwe(0.1),0.9)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
