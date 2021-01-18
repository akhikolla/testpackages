\name{forcelayout}
\alias{forcelayout}
\title{
~ Main function: forcelayout ~
}
\description{
 forcelayout method is a dynamic method showing longitudinal data set evolution.
}
\usage{
forcelayout(schedule, webinteract, ttime)
}
\arguments{
 \item{schedule}{[numeric] or [integer]: A (non-empty) vector of data values.}
 
 \item{webinteract}{[bool]: Is the function used in interactive mode?}

 \item{ttime}{[string]: A (non-empty) time-unit value to fetch data from data.frame.}
}
\examples{
  if(interactive()){
    forcelayout(weekschedule)
  }
  if(!interactive()){
    forcelayout(weekschedule, webinteract = FALSE, ttime = "Monday")
  }
}
