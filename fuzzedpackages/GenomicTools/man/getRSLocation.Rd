\name{getRSLocation}
\alias{getRSLocation}

\title{
  Get the location of an SNP.
}

\description{
  This function gets the exact location of a SNPs, based on a specific database.
}
\usage{
  getRSLocation(rs, species=NULL)
}

\arguments{
  \item{rs}{SNP rs number}
  \item{species}{Species identifier}
}

\details{
  This function can download from any specifi ensembl database the location of an SNP (as given as RS number) for a specific species.
}
\value{
  A list with Chr and Location (BP) of an SNP.
}

\examples{
# Depending on the internet connection the following call
# can take up 10 seconds and hence does not pass the Cran
# checks. This is why it is wrapped into \dontrun 

  \dontrun{
    getRSLocation("rs16033432", species="Gallus_gallus")
  }
}

\author{
  Daniel Fischer
}

 
