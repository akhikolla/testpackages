%File datasailr-package.Rd

\name{datasailr-package}
\alias{datasailr-package}
\docType{package}
\title{
  Row by Row Data Processing Tool, Using 'DataSailr' Script
}
\description{
  A row by row data processing tool. You can write data manipulation code in 'DataSailr' script which is specially intended for the data manipulation. The package uses 'libsailr' (C/C++ library) for its 'DataSailr' code parsing and its execution.
}
\details{
  The main function, datasailr::sail(), takes dataframe as its 1st argument and 'DataSailr' script as its 2nd argument, and returns the processing results as dataframe. This package works with 'libsailr' library, which conducts arithmetic calculations and string/character manipulations following 'Sailr' script (which is a base version of 'DataSailr' script). 'DataSailr' package passes both 'DataSailr' script and values of each row of dataframe to the 'libsailr', and 'libsailr' processes those values of each row. The processing results are stored and are returned as dataframe.
}
\author{
Toshi Umehara, toshi@niceume.com.
Maintainer: Toshi Umehara <toshi@niceume.com>
}
\keyword{ package }
\seealso{
  \code{\link{sail}}.
}
\examples{
library(datasailr)
data(mtcars)

code_example = 'powerful = .
  if( 145 >= hp && hp > 0 ){
    powerful = 0
  }else if( hp > 145 ){
    powerful = 1
  }else{
    powerful = 2
  }

  germany = re/(^Merc|^Porsche|^Volvo)/
  japan = re/(^Mazda|^Honda|^Toyota)/

  if ( _rowname_ =~  germany ) { 
      country = "Germany" 
  }else if( _rowname_ =~ japan ){
      country = "Japan"
  }else {
      country = "Other"
  }
  company = rexp_matched(1)
'

sail(mtcars, code_example)
}

