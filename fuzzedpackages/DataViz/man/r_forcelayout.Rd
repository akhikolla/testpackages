\name{r_forcelayout}
\alias{r_forcelayout}
\title{
~ R graphics function: r_forcelayout ~
}
\description{
 Used when interactive is false and creates a plot through R of this data visualisation method. 
}
\usage{
r_forcelayout(schedule, ttime)
}
\arguments{
 \item{schedule}{[string]: A (non-empty) data.frame of data values.}
 
 \item{ttime}{[string]: A (non-empty) time-unit value to fetch data from data.frame.}

}
\examples{
 if(interactive()){
    r_forcelayout(weekschedule,ttime = "Tuesday") 
 }
}