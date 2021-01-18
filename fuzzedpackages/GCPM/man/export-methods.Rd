\name{export-methods}
\docType{methods}
\alias{export-methods}
\alias{export}
\alias{export,GCPM-method}
\alias{export,GCPM,character,character,numeric-method}
\alias{export,GCPM,character,character,missing-method}
\alias{export,GCPM,character,missing,numeric-method}
\alias{export,GCPM,character,missing,missing-method}
\alias{export,GCPM,missing,character,numeric-method}
\alias{export,GCPM,missing,missing,numeric-method}
\alias{export,GCPM,missing,missing,missing-method}
\alias{export,GCPM,missing,character,missing-method}

\title{Export Main Results}
\description{
This method provides an easy way to export the main results of the portfolio
(i.e. after running \code{analyze}). A summary file and the portfolio loss
distribution (PDF and CDF) are exported to \code{path.out}. With the help of
\code{file.format} one can specify the csv format (\dQuote{csv1} or
\dQuote{csv2}). If a vector \code{alpha} of loss levels is specified,
risk contributions to EC, VaR and ES are also exported according to level(s)
\code{alpha}.}
\usage{
export(this,path.out,file.format,alpha)
}
\arguments{
  \item{this}{Object of class \code{GCPM}}
  \item{path.out}{string specifying the output path}
  \item{file.format}{string specifying the file format (i.e \dQuote{csv1} or
  \dQuote{csv2})}
  \item{alpha}{numeric vector with loss levels between 0 and 1}
}
\keyword{methods}
