\name{buildPointls}
\alias{buildPointls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Read a data file
}
\description{
  This reads a data file into a data frame, performs various checks, assign namesto columns, and can select rows.
}
\usage{
buildPointls(dataFile = blackbox.getOption("dataFile"), respCols = NULL,
             subsetRows = NULL, ycolname, cleanResu = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataFile}{
  Name of data file
}
  \item{respCols}{
  A way to select response columns in later analyses (see Details). NULL or a numeric vector.
}
  \item{subsetRows}{
  A set of rows to select. All rows are retained if this is NULL
  }
  \item{ycolname}{
  A name to be given to the response variable; willbe used in many further outputs.
}
  \item{cleanResu}{
  A connection, or a character string naming a file for some nicely formated output. If \code{""} (the default), print to the standard output connection.
}
}
\details{
  The input file is a an ASCII numeric data table with the following columns. The first columns contain values of all canonical parameters of the model in canonical order, as given by \code{blackbox.getOption{"ParameterNames"}}. Pairs of lines may have identical parameter vectors, but not more than pairs. The next columns may all be used as response variables.

  \code{respCols} identifies columns that will be used to construct the smoothed response (but all columns are retained in this function's return value). If it is NULL, then the last column will be used. If a numeric vector, it identifies response columns (where column 1 is the first column after the parameters columns) which values will be summed to construct the response variable.
}
\value{
A data frame with as many columns as the input table.
As a side effect, the function sets the \code{blackbox.options ycolname} and \code{respCols} to respectively the input \code{ycolname} and to the column names deduced from the input \code{respCols} indices.
}
