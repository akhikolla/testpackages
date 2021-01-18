\name{rcpp_forcelayout}
\alias{rcpp_forcelayout}
\title{
  ~ C++ called function: rcpp_forcelayout ~
}
\description{
  This function takes the inputs from \link{forcelayout}, and writes the data in Json array, then this function calls a windows cmd function to open an index.html in the extdata.
}
\usage{
rcpp_forcelayout(schedule, path)	
}
\arguments{
  \item{schedule}{A number column}

  \item{path}{The path for the library}
}
\examples{
 if(interactive()){
  rcpp_throwchart(weekschedule,path.package("DataViz"))
}
}
