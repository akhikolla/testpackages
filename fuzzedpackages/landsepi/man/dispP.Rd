\name{dispP}
\alias{dispP}
\alias{dispP_1}
\alias{dispP_2}
\alias{dispP_3}
\alias{dispP_4}
\alias{dispP_5}
\docType{data}
\title{
  Dispersal matrices for rust fungi of cereal crops.
}
\description{
  Five vectorised dispersal matrices of pathogens as typified by rust fungi of cereal crops (genus \emph{Puccinia}), 
  and associated with landscapes 1 to 5 (composed of 155, 154, 152, 153 and 156 fields, respectively).
}
\usage{
 dispP_1
 dispP_2
 dispP_3
 dispP_4
 dispP_5
}
\format{
  The format is:
 num [1:24025] 8.81e-01 9.53e-04 7.08e-10 1.59e-10 3.29e-06 ...
}
\details{
  The pathogen dispersal matrix gives the probability for a pathogen in a field i (row) to migrate to field i' (column) through dispersal. It is computed based on a dispersal kernel and the euclidian distance between each point in fields i and i', using the CaliFloPP algorithm (Bouvier et al. 2009).  
  The dispersal kernel is an isotropic power-law function of equation:  
  f(x)=((b-2)*(b-1)/(2*pi*a^2)) * (1 + x/a)^(-b)  
  with a=40 a scale parameter and b=7 related to the weight of the dispersal tail.  
  The expected mean dispersal distance is given by 2*a/(b-3)=20 m.
%%  ~~ If necessary, more details than the __description__ above ~~
}
\references{
  Bouvier A, Kiêu K, Adamczyk K, Monod H. Computation of the integrated flow of particles between polygons. Environ. Model Softw. 2009;24(7):843-9. doi: http://dx.doi.org/10.1016/j.envsoft.2008.11.006.
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
dispP_1
summary(dispP_1)
## maybe str(dispP_1) ; plot(dispP_1) ...
}
\keyword{datasets}
