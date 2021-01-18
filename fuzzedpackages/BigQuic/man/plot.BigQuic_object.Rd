\name{plot.BigQuic_object}
\alias{plot.BigQuic_object}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Plot
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Makes plot of the precision matrix showing non-zero values.  The diagonal is 
shown in only black because the agreement with itself is not highly 
interesting.  Negative relations are shown in green and positive in red.  The 
saturation indicates the normalized strength of the relation.  The matrix is 
symmetric and technically only the lower or upper triangle would suffice to 
provide identical information.  
}
\usage{
  \method{plot}{BigQuic_object}(x, ...)
}
\arguments{
  \item{x}{
    The BigQuic object, which will have its optimal precision matrix plotted.     
  }
  \item{...}{
    plot can take a variety of arguments depending on the type, that is 
    represented by ...
  }
}

\keyword{BigQuic}
