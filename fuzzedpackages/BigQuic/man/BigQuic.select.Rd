\name{BigQuic.select}
\alias{BigQuic.select}
\title{
%%   ~~ data name/kind ... ~~
BigQuic Select
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Selects the optimal lambda value from those in the BigQuic_object, i.e. BigQuic
Result.  
}
\usage{BigQuic.select(BigQuic_result = NULL, stars.thresh = 0.1, 
       stars.subsample.ratio = NULL, rep.num = 20, verbose = TRUE, 
       verbose2 = 0)}
\arguments{
  \item{BigQuic_result}{
    A BigQuic_object returned from running BigQuic.  
  }
  \item{stars.thresh}{
    The threshold used in the Stars selection method for choosing a lambda
  }
  \item{stars.subsample.ratio}{
    The ratio giving how large the subsamples will be for Stars, if null 
    there is a heuristic calculation.  
  }
  \item{rep.num}{
    Number of times to do the repetition in Stars.  
  }
  \item{verbose}{
    Controls the level of verbosity in a part of the code.  
  }
  \item{verbose2}{
    Controls the level of verbosity in another section of code.  
  }
}
\keyword{BigQuic}
