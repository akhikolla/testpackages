\name{generate_sample}
\alias{generate_sample}
\title{
%%   ~~ data name/kind ... ~~
Generate Sample
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Generates a sample data set for using with BigQuic, the default seed is 1 for 
reproducibility.  For high dimensional data, choose p much larger than n.  
}
\usage{generate_sample(n = 200, p = 150, seed = NULL)}
\arguments{
  \item{n}{
    The number of rows in the resulting data set.  
  }
    \item{p}{
    The number of columns in the resulting data set.  
  }
  \item{seed}{
    A seed for the random number generator in R.    
  }
}
\keyword{BigQuic}
