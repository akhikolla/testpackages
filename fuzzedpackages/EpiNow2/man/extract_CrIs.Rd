% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarise.R
\name{extract_CrIs}
\alias{extract_CrIs}
\title{Extract Credible Intervals Present}
\usage{
extract_CrIs(summarised)
}
\arguments{
\item{summarised}{A data frame as processed by \code{calc_CrIs}}
}
\value{
A numeric vector of credible intervals detected in the data frame.
}
\description{
\ifelse{html}{\out{<a href='https://www.tidyverse.org/lifecycle/#stable'><img src='figures/lifecycle-stable.svg' alt='Stable lifecycle'></a>}}{\strong{Stable}}
Helper function to extract the credible intervals present in a data frame.
}
\examples{
samples <- data.frame(value = 1:10, type = "car")
summarised <- calc_CrIs(samples, summarise_by = "type",
                        CrIs = c(seq(0.05, 0.95, 0.05))) 
extract_CrIs(summarised)
}
