%File sail.Rd

\name{sail}
\alias{sail}
\title{Run 'DataSailr' script and manipulate dataframe }
\description{
  datasailr::sail() is the main function of this 'DataSailr' package. This function takes data.frame as the 1st argument and 'DataSailr' script as the 2nd argument. 
}
\usage{
  sail( df, code, fullData = TRUE, rowname = "_rowname_", stringsAsFactors = FALSE)
}
\arguments{
  \item{df}{ data.frame to be processed.}
  \item{code}{ 'DataSailr' script that instructs how to manipulate values and strings.}
  \item{fullData}{ When fullData is set FALSE, only the columns corresponding to left-hand-side(LHS) variables in 'DataSailr' scripts are returned. (e.g. In 'DataSailr' script, bmi = weight / height / height , the variable bmi is the LHS.) If set TRUE, all the columns are returned. }
  \item{rowname}{ Dataframe's row names are accessed via the variable name specified by this option. The default variable name is _rowname_, meaning that _rowname_ represents dataframe's row names in 'DataSailr' script, though this does not allow assignment. Even if you assign some value to this variable, they are not reflected in the row names of the resultant dataframe. }
  \item{stringsAsFactors}{ Columns of character vectors (= string vectors) that appear as left-hand-side(LHS) variables in 'DataSailr' scripts can be returned as either character vectors or factor vectors. If stringsAsFactors is set to FALSE(default), character columns are returned as character vectors. If stringsAsFactors is set to TRUE, character columns are returned as factor vectors. }
}
\value{
  If fullData argument is set FALSE, return value is a data.frame. If fullData is TRUE, the result is created by cbind(), which concatenate the original input and newly created dataframe. For example, if the original input is 'data.table', the return value is 'data.table'.
}
\examples{
library(datasailr)
data(iris)

iris_code = '
iris_type = .
if ( Species == "setosa" ) { 
    type = 1
    type_char = "A"
}else if( Species == "versicolor" ) { 
    type = 2
    type_char = "B"
}else if(Species == "virginica" ) { 
    type = 3
    type_char = "C"
}
'
sail(iris, iris_code)
}

