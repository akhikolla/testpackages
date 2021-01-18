  \name{Simp}
\alias{Simp}
\docType{data}
\title{Simulated Chimp (Simp) Data}
\description{This is a simulated movement track from a continuous auto-correlated process of total duration 60 time units with four behavioral phases that switch at times T=10, 20 and 50 from, variously, higher or lower velocity and longer or shorter characteristic time scale of autocorrelation. The original data was simulated with time intervals of 0.01 (leading to a complete track with 6000 data points). 200 points were randomly sampled from the "true" movement track, such that the intervals between locations are random with mean interval 0.3 units. 
}
\usage{data(Simp)}
\format{ Data frame (of class "track") containing 200 rows and three columns:
\describe{
\item{Time}{times of observation} 
\item{X, Y}{coordinates}   
}}
\source{
For more details, see: \url{http://wiki.cbr.washington.edu/qerm/index.php/Behavioral_Change_Point_Analysis}
}
\examples{
data(Simp)
plot(Simp)
}
\keyword{datasets}
