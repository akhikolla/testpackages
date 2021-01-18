\name{Selvar}
\alias{Selvar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of heritability in high dimensional sparse linear mixed models using variable selection.
}
\description{
This function selects active components in sparse linear mixed models in order to estimate heritability. The selection allows us to reduce the size of the data sets which improves the accuracy of the estimations. Our package also provides a confidence interval for the estimated heritability.
}
\usage{
Selvar(Y,Z,X,thresh_vect,nb_boot=80,nb_repli=50,CI_level=0.95,nb_cores=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
Vector of observations of size n.
}
  \item{Z}{
Matrix with genetic information of size n x N.
}
  \item{X}{
Matrix of fixed effects of size n x d.
}
  \item{thresh_vect}{
Vector of thresholds in the stability selection step: the higher the threshold, the smallest the set of selected components. 
}
  \item{nb_boot}{
Number of subsamples of Y to apply our bootstrap technique. The value by default is 80.
}
  \item{nb_repli}{
Number of replications in the stability selection. The value by default is 50.
}
  \item{CI_level}{
Level of the confidence interval for the estimation of the heritability. The value by default is 0.95.
}
\item{nb_cores}{
Number of cores of the computer. It is used for parallelizing the computations. The value by default is 1.
}
}

\value{
%%  If it is a LIST, use
  \item{heritability}{Estimation of the heritability}
 \item{CI_up}{Upper bound of the confidence interval for the estimated heritability}
\item{CI_low}{Lower bound of the confidence interval for the estimated heritability}
\item{selec_ind}{Indexes of the columns of the selected components}}

\author{
Anna Bonnet and Celine Levy-Leduc
}

\examples{
library(EstHer)
data(Y)
data(W)
data(X)
Z=scale(W,center=TRUE,scale=TRUE)
res=Selvar(Y,Z,X,thresh_vect=c(0.7,0.8,0.9),nb_boot=80,nb_repli=50,CI_level=0.95,nb_cores=1) 
res$heritability
res$CI_low
res$CI_up
}


\keyword{variable selection}
\keyword{heritability}
\keyword{linear mixed models}
