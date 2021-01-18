\name{getFit.regmed.grid}
\alias{getFit.regmed.grid}
\title{Retrieve regmed object from regmed.grid object}
\description{Retrieve regmed object from regmed.grid object using a selected lambda value.  If no lambda value is provided, the best fit (based on minimum BIC) is returned.}
\usage{
  getFit.regmed.grid(obj, lambda)  
}
\arguments{
  \item{obj}{regmed.grid object}
  \item{lambda}{lambda value used to select regmed object from
    regmed.grid object (lambda must be a value used in obj).
    If left NULL, returned value will be based on minimum BIC.}
}
\value{
  regmed object
}
