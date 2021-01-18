\name{PseudoR2.Cure}
\alias{PseudoR2.Cure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pseudo R2 criterion
}
\description{
Pseudo R2 criterion for a mixture of population (susceptible and nonsusceptible populations)
}
\usage{
PseudoR2.Cure(ygene, ydelai, yetat, strate, ordered = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ygene}{
The main variable of interest
}
  \item{ydelai}{
The right censored delay until the event
}
  \item{yetat}{
The censoring indicator
}
  \item{strate}{
The varaiables acting on the nonsusceptible or cured population 
}
  \item{ordered}{
A value of \code{TRUE} or \code{FALSE} indicating weither or not the times to event are ordered
}
}

\value{
A pseudo R2 value lying between 0 and 1.
}
\references{
Mbogning, C. and Broet, P. (2016). Bagging survival tree procedure for
variable selection and prediction in the presence of nonsusceptible patients.
BMC bioinformatics, 17(1), 1.
}
\author{
Cyprien Mbogning and Philippe Broet
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Bagg_Surv}} \code{\link{Bagg_pred_Surv}} \code{\link{improper_tree}}
}
\examples{
data(burn)
PseudoR2.Cure(ygene = burn$Z3, ydelai = burn$T3, yetat = burn$D3, strate = burn$Z2)
PseudoR2.Cure(ygene = burn$Z2, ydelai = burn$T3, yetat = burn$D3, strate = burn$Z2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ documentation }% use one of  RShowDoc("KEYWORDS")
\keyword{ tree }% __ONLY ONE__ keyword per line
\keyword{survival}