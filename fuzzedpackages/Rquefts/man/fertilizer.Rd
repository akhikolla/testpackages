\name{Fertilizers}

\alias{fertilizers}
\alias{nutrientRates}

\title{
Helper functions to go from fertilizers to nutrients
}


\description{
Computes the amount of nutrients given a rate of fertilizer. 
}

\usage{
fertilizers()
nutrientRates(supply, treatment) 
}

\arguments{
  \item{supply}{data.frame with columns "N", "P", "K" expressed as percentage of the product (row)}
  \item{treatment}{amounts applied}
}



\examples{
# fertilizer product list
fert <- fertilizers()
# shortening some of the names for display
fert[,2] = substr(fert[,2], 1, 20)
# contents are expressed as a percentage.
fert

myferts <- fert[c(8,15), ] 
nutrientRates(myferts, c(100,50))
}

